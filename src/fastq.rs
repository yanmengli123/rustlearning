//! ============================================================
//! 模块2: FASTQ质量控制 (fastq)
//! ============================================================
//! 本模块提供FASTQ格式测序数据的质量控制和分析功能。
//! 包括：质量统计、过滤、GC含量分布、长度分布、
//! 滑动窗口过滤等。
//!
//! FASTQ格式说明：
//! 第1行: @header（以@开头的序列标识符）
//! 第2行: 序列碱基（A/T/G/C/N）
//! 第3行: +（分隔行，可选重复header）
//! 第4行: 质量分数（ASCII编码，Phred+33格式）
//!
//! 质量分数Q与错误率P的关系:
//!   Q = -10 * log10(P)
//!   Q20 = 1%错误率，Q30 = 0.1%错误率
//!
//! 设计原则：
//! - 支持gzip压缩文件（.gz后缀自动检测）
//! - 使用Phred+33编码（Illumina 1.8+标准）
//! - 返回Python字典便于数据分析
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use std::collections::HashMap;  // Rust 标准哈希表
use std::fs::File;              // 文件操作
use std::io::{self, BufRead, BufReader, Write, BufWriter};  // IO操作

/// -----------------------------------------------------------
/// 单条FASTQ记录结构体
/// -----------------------------------------------------------
/// 存储FASTQ格式的单条测序读段信息
/// 包含三个字段：
///   header - 序列标识符（不含@前缀）
///   seq    - 碱基序列（A/T/G/C/N）
///   qual   - 质量分数字符串（ASCII编码）
struct FastqRecord {
    header: String,  // 序列标识符
    seq: String,     // 碱基序列
    qual: String,    // 质量分数字符串
}

/// -----------------------------------------------------------
/// 解析FASTQ文件，返回所有记录
/// -----------------------------------------------------------
/// 参数: path - FASTQ文件路径（支持.gz压缩格式）
/// 返回: Vec<FastqRecord> - 所有FASTQ记录列表
/// 算法:
///   1. 根据文件后缀判断是否为gzip压缩格式
///   2. 逐行读取，每次读取4行（header/seq/+/qual）
///   3. 以@开头的行为header行
///   4. 跳过不符合格式的行
/// 用途: 批量读取FASTQ文件用于后续统计分析
fn parse_fastq_records(path: &str) -> io::Result<Vec<FastqRecord>> {
    let file = File::open(path)?;  // 打开文件
    // 根据后缀名判断是否为gzip压缩文件
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        let gz = flate2::read::GzDecoder::new(file);  // gzip解码器
        Box::new(BufReader::new(gz))  // 使用缓冲读取器包装
    } else {
        Box::new(BufReader::new(file))  // 普通文件缓冲读取
    };

    let mut records = Vec::new();  // 存储所有记录
    let mut lines = reader.lines();  // 行迭代器

    // 循环读取，每次处理一条FASTQ记录（4行）
    while let Some(Ok(header)) = lines.next() {
        if !header.starts_with('@') {
            continue;  // 跳过非header行
        }
        // 读取后续3行：序列行、+行、质量行
        if let (Some(Ok(seq)), Some(Ok(_plus)), Some(Ok(qual))) =
            (lines.next(), lines.next(), lines.next())
        {
            // 构建FastqRecord，header去掉前导@符号
            records.push(FastqRecord {
                header: header[1..].to_string(),  // 去掉@前缀
                seq,
                qual,
            });
        }
    }
    Ok(records)
}

/// -----------------------------------------------------------
/// 计算Phred质量分数
/// -----------------------------------------------------------
/// 参数: qual - 质量分数ASCII字符串
/// 返回: Vec<u8> - 质量分数数组
/// 算法:
///   Phred+33编码：ASCII值减去33得到质量分数
///   例如: '!' (ASCII 33) = Q0, 'A' (ASCII 65) = Q32
///   Q值越高表示碱基识别越可靠
fn qual_scores(qual: &str) -> Vec<u8> {
    qual.bytes().map(|b| b - 33).collect()  // ASCII减33得到Phred分数
}

/// -----------------------------------------------------------
/// FASTQ质控统计
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用（PyO3自动传入）
///   path - FASTQ文件路径
/// 返回: Python字典，包含全面的质控指标
/// 统计指标:
///   - total_reads: 总读段数
///   - total_bases: 总碱基数
///   - avg_length: 平均读段长度
///   - min_length / max_length: 最短/最长读段
///   - avg_qual: 平均质量分数
///   - q20_bases / q30_bases: Q20/Q30碱基数
///   - q20_rate / q30_rate: Q20/Q30比率
///   - gc_content: GC含量
///   - n_bases: N碱基数量
/// 用途: 评估测序数据质量，决定是否需要过滤或重新测序
#[pyfunction]
pub fn fastq_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total_reads = records.len();  // 总读段数
    let mut total_bases: u64 = 0;    // 总碱基数
    let mut q20_bases: u64 = 0;      // Q≥20的碱基数
    let mut q30_bases: u64 = 0;      // Q≥30的碱基数
    let mut gc_count: u64 = 0;       // GC碱基总数
    let mut n_count: u64 = 0;        // N碱基数
    let mut length_sum: u64 = 0;     // 长度总和（用于计算平均值）
    let mut min_len = usize::MAX;    // 最短长度
    let mut max_len = 0usize;        // 最长长度
    let mut total_qual_sum: f64 = 0.0;  // 质量分数总和
    let mut per_base_qual: HashMap<usize, Vec<f64>> = HashMap::new();  // 每个位点的质量

    // 遍历每条记录进行统计
    for rec in &records {
        let slen = rec.seq.len();
        total_bases += slen as u64;     // 累加碱基数
        length_sum += slen as u64;      // 累加长度
        min_len = min_len.min(slen);    // 更新最短长度
        max_len = max_len.max(slen);    // 更新最长长度

        // 获取质量分数数组
        let scores = qual_scores(&rec.qual);
        for (i, &q) in scores.iter().enumerate() {
            total_qual_sum += q as f64;     // 累加质量分数
            per_base_qual.entry(i).or_default().push(q as f64);  // 记录每个位点
            if q >= 20 {
                q20_bases += 1;  // 统计Q20碱基
            }
            if q >= 30 {
                q30_bases += 1;  // 统计Q30碱基
            }
        }
        // 统计GC和N碱基
        for c in rec.seq.chars() {
            match c {
                'G' | 'g' | 'C' | 'c' => gc_count += 1,  // GC碱基
                'N' | 'n' => n_count += 1,  // 未知碱基
                _ => {}
            }
        }
    }

    // 计算各项平均值和比率
    let avg_len = if total_reads > 0 {
        length_sum as f64 / total_reads as f64
    } else {
        0.0
    };
    let avg_qual = if total_bases > 0 {
        total_qual_sum / total_bases as f64
    } else {
        0.0
    };
    let q20_rate = if total_bases > 0 {
        q20_bases as f64 / total_bases as f64
    } else {
        0.0
    };
    let q30_rate = if total_bases > 0 {
        q30_bases as f64 / total_bases as f64
    } else {
        0.0
    };
    let gc_rate = if total_bases > 0 {
        gc_count as f64 / total_bases as f64
    } else {
        0.0
    };

    // 计算每个位点的平均质量（用于质量分布图）
    let mut per_base_mean_qual: Vec<f64> = Vec::new();
    let max_pos = per_base_qual.keys().max().copied().unwrap_or(0);
    for i in 0..=max_pos {
        if let Some(scores) = per_base_qual.get(&i) {
            let mean = scores.iter().sum::<f64>() / scores.len() as f64;
            per_base_mean_qual.push(mean);
        }
    }

    // 构建Python字典返回所有统计指标
    let dict = PyDict::new_bound(py);
    dict.set_item("total_reads", total_reads)?;
    dict.set_item("total_bases", total_bases)?;
    dict.set_item("avg_length", avg_len)?;
    dict.set_item("min_length", if min_len == usize::MAX { 0 } else { min_len })?;
    dict.set_item("max_length", max_len)?;
    dict.set_item("avg_qual", avg_qual)?;
    dict.set_item("q20_bases", q20_bases)?;
    dict.set_item("q30_bases", q30_bases)?;
    dict.set_item("q20_rate", q20_rate)?;
    dict.set_item("q30_rate", q30_rate)?;
    dict.set_item("gc_content", gc_rate)?;
    dict.set_item("n_bases", n_count)?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// FASTQ质量过滤
/// -----------------------------------------------------------
/// 参数:
///   input    - 输入FASTQ文件路径
///   output   - 输出FASTQ文件路径
///   min_len  - 最小读段长度阈值
///   min_qual - 最小平均质量分数阈值
/// 返回: 通过过滤的读段数量
/// 算法:
///   1. 遍历所有读段
///   2. 检查长度是否≥min_len
///   3. 计算平均质量分数是否≥min_qual
///   4. 通过过滤的读段写入输出文件
/// 用途: 去除低质量读段，提高后续分析准确性
#[pyfunction]
pub fn fastq_filter(
    input: &str,
    output: &str,
    min_len: usize,
    min_qual: f64,
) -> PyResult<usize> {
    let records = parse_fastq_records(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(file);  // 使用缓冲写入提高性能

    let mut kept = 0usize;  // 通过过滤的读段计数
    for rec in &records {
        // 检查长度过滤
        if rec.seq.len() < min_len {
            continue;  // 长度不足，跳过
        }
        // 计算平均质量分数
        let scores = qual_scores(&rec.qual);
        let mean_q = scores.iter().map(|&q| q as f64).sum::<f64>() / scores.len() as f64;
        if mean_q < min_qual {
            continue;  // 质量不足，跳过
        }
        // 写入FASTQ格式（4行）
        writeln!(writer, "@{}", rec.header).unwrap();
        writeln!(writer, "{}", rec.seq).unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "{}", rec.qual).unwrap();
        kept += 1;
    }
    Ok(kept)
}

/// -----------------------------------------------------------
/// 每个位点的质量分布
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   path - FASTQ文件路径
/// 返回: Python字典，包含每个位点的平均和中位质量
/// 算法:
///   1. 收集所有读段在每个位点的质量分数
///   2. 计算每个位点的平均值和中位数
///   3. 用于生成类似FastQC的per-base quality plot
/// 用途: 识别质量下降的位点（如测序末尾）
#[pyfunction]
pub fn per_base_quality(py: Python, path: &str) -> PyResult<PyObject> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // 收集每个位点的所有质量分数
    let mut per_pos: HashMap<usize, Vec<f64>> = HashMap::new();
    for rec in &records {
        let scores = qual_scores(&rec.qual);
        for (i, &q) in scores.iter().enumerate() {
            per_pos.entry(i).or_default().push(q as f64);
        }
    }

    // 计算每个位点的统计指标
    let dict = PyDict::new_bound(py);
    let mut positions: Vec<usize> = per_pos.keys().copied().collect();
    positions.sort();  // 按位点排序
    let mut means = Vec::new();    // 平均质量
    let mut medians = Vec::new();  // 中位质量
    for &pos in &positions {
        let mut scores = per_pos[&pos].clone();
        scores.sort_by(|a, b| a.partial_cmp(b).unwrap());  // 排序用于求中位数
        let mean = scores.iter().sum::<f64>() / scores.len() as f64;
        let median = scores[scores.len() / 2];  // 中位数
        means.push(mean);
        medians.push(median);
    }
    dict.set_item("positions", positions)?;
    dict.set_item("means", means)?;
    dict.set_item("medians", medians)?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// Read长度分布统计
/// -----------------------------------------------------------
/// 参数: path - FASTQ文件路径
/// 返回: HashMap<usize, usize> - 长度 -> 读段数量
/// 用途: 检测测序读段长度是否一致，识别接头污染
#[pyfunction]
pub fn length_distribution(path: &str) -> PyResult<HashMap<usize, usize>> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut dist: HashMap<usize, usize> = HashMap::new();
    for rec in &records {
        *dist.entry(rec.seq.len()).or_insert(0) += 1;  // 统计每种长度的读段数
    }
    Ok(dist)
}

/// -----------------------------------------------------------
/// GC含量分布统计
/// -----------------------------------------------------------
/// 参数:
///   path     - FASTQ文件路径
///   bin_size - GC含量分箱大小（如0.1表示10%为一个箱）
/// 返回: HashMap<String, usize> - GC区间 -> 读段数量
/// 算法:
///   1. 计算每条读段的GC含量
///   2. 按bin_size分箱统计
///   3. 用于检测GC偏好、污染等问题
/// 用途: 识别异常GC分布（如接头污染、物种混合）
#[pyfunction]
pub fn gc_distribution(path: &str, bin_size: f64) -> PyResult<HashMap<String, usize>> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut dist: HashMap<String, usize> = HashMap::new();
    for rec in &records {
        let gc = super::sequence::gc_content(&rec.seq);  // 使用sequence模块计算GC
        // 将GC含量分箱
        let bin = (gc / bin_size).floor() * bin_size;
        let key = format!("{:.1}-{:.1}", bin, bin + bin_size);  // 区间标签
        *dist.entry(key).or_insert(0) += 1;
    }
    Ok(dist)
}

/// -----------------------------------------------------------
/// 滑动窗口质量过滤
/// -----------------------------------------------------------
/// 参数:
///   path          - FASTQ文件路径
///   window_size   - 滑动窗口大小
///   min_avg_qual  - 窗口内最小平均质量阈值
/// 返回: (总读段数, 通过过滤读段数)
/// 算法:
///   1. 对每条读段，使用滑动窗口扫描质量分数
///   2. 如果任意窗口的平均质量低于阈值，则过滤该读段
///   3. 类似于Trimmomatic的滑动窗口修剪策略
/// 用途: 去除3'端质量下降的读段，保留高质量区域
#[pyfunction]
pub fn sliding_window_filter(
    path: &str,
    window_size: usize,
    min_avg_qual: f64,
) -> PyResult<(usize, usize)> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total = records.len();  // 总读段数
    let mut passed = 0usize;   // 通过过滤的读段数

    for rec in &records {
        let scores = qual_scores(&rec.qual);
        let mut ok = true;  // 标记是否通过过滤
        // 处理短于窗口的读段
        if scores.len() < window_size {
            let mean = scores.iter().map(|&q| q as f64).sum::<f64>() / scores.len() as f64;
            if mean < min_avg_qual {
                ok = false;
            }
        } else {
            // 滑动窗口扫描
            for window in scores.windows(window_size) {
                let mean = window.iter().map(|&q| q as f64).sum::<f64>() / window_size as f64;
                if mean < min_avg_qual {
                    ok = false;  // 发现低质量窗口
                    break;       // 提前终止
                }
            }
        }
        if ok {
            passed += 1;  // 通过过滤
        }
    }
    Ok((total, passed))
}
