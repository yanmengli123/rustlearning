//! ============================================================
//! 模块5: SAM/BAM文件处理 (sam_bam)
//! ============================================================
//! 本模块提供SAM（Sequence Alignment/Map）格式文件的
//! 解析和统计功能。
//!
//! SAM格式说明：
//!   SAM是存储高通量测序reads比对结果的标准文本格式
//!   由头部（@开头）和比对记录组成
//!   每条记录包含11个必选字段和可选标签
//!
//! 字段说明（11个必选字段）：
//!   QNAME - 查询序列名称（read ID）
//!   FLAG  - 比对标志（位运算编码各种信息）
//!   RNAME - 参考序列名称（染色体）
//!   POS   - 比对起始位置（1-based）
//!   MAPQ  - 比对质量（0-60，越高越可靠）
//!   CIGAR - 比对描述字符串
//!   RNEXT - 配对read的参考序列
//!   PNEXT - 配对read的位置
//!   TLEN  - 插入片段长度
//!   SEQ   - 序列碱基
//!   QUAL  - 质量分数
//!
//! FLAG位说明（二进制位）：
//!   0x1   - read是配对的
//!   0x2   - read正确配对
//!   0x4   - read未比对
//!   0x8   - mate未比对
//!   0x10  - read反向互补比对
//!   0x20  - mate反向互补比对
//!   0x40  - read是第一个配对
//!   0x80  - read是第二个配对
//!   0x100 - 非主要比对
//!   0x200 - 未通过质量控制
//!   0x400 - PCR或光学重复
//!   0x800 - 补充比对
//!
//! 设计原则：
//! - 支持SAM文本格式解析
//! - 提供丰富的统计和过滤功能
//! - 使用位运算解析FLAG字段
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use std::collections::HashMap;  // Rust 标准哈希表
use std::fs::File;              // 文件操作
use std::io::{self, BufRead, BufReader};  // IO操作

/// -----------------------------------------------------------
/// SAM记录结构体
/// -----------------------------------------------------------
/// 存储SAM格式单条比对记录的所有字段
/// 包含11个必选字段和可选标签（TAG）
#[allow(dead_code)]
#[derive(Debug)]
pub struct SamRecord {
    pub qname: String,              // 查询序列名称（read ID）
    pub flag: u16,                  // 比对标志（位运算编码）
    pub rname: String,              // 参考序列名称（染色体）
    pub pos: i64,                   // 比对起始位置（1-based）
    pub mapq: u8,                   // 比对质量（0-60）
    pub cigar: String,              // CIGAR字符串
    pub rnext: String,              // 配对read参考序列
    pub pnext: i64,                 // 配对read位置
    pub tlen: i64,                  // 插入片段长度
    pub seq: String,                // 序列碱基
    pub qual: String,               // 质量分数
    pub tags: HashMap<String, String>,  // 可选标签（如NM:i:0）
}

/// -----------------------------------------------------------
/// 解析SAM行
/// -----------------------------------------------------------
/// 参数: line - SAM格式的一行文本
/// 返回: Option<SamRecord> - 解析成功返回Some，否则None
/// 算法:
///   1. 跳过以@开头的头部行
///   2. 按制表符分割，提取11个必选字段
///   3. 剩余字段解析为标签（TAG:TYPE:VALUE格式）
/// 用途: 逐行解析SAM文件
pub fn parse_sam_line(line: &str) -> Option<SamRecord> {
    if line.starts_with('@') {
        return None;  // 跳过头部行
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return None;  // 字段不足，无效行
    }

    // 解析可选标签（第12个字段开始）
    let mut tags = HashMap::new();
    for field in &fields[11..] {
        if let Some(idx) = field.find(':') {
            let tag = &field[..idx];       // 标签名（如"NM"）
            let val = &field[idx + 2..];   // 标签值（跳过类型标识符）
            tags.insert(tag.to_string(), val.to_string());
        }
    }

    // 构建SamRecord结构体
    Some(SamRecord {
        qname: fields[0].to_string(),    // QNAME
        flag: fields[1].parse().unwrap_or(0),  // FLAG
        rname: fields[2].to_string(),    // RNAME
        pos: fields[3].parse().unwrap_or(0),   // POS
        mapq: fields[4].parse().unwrap_or(0),  // MAPQ
        cigar: fields[5].to_string(),    // CIGAR
        rnext: fields[6].to_string(),    // RNEXT
        pnext: fields[7].parse().unwrap_or(0), // PNEXT
        tlen: fields[8].parse().unwrap_or(0),  // TLEN
        seq: fields[9].to_string(),      // SEQ
        qual: fields[10].to_string(),    // QUAL
        tags,                             // 可选标签
    })
}

/// -----------------------------------------------------------
/// 读取SAM文件所有记录
/// -----------------------------------------------------------
/// 参数: path - SAM文件路径
/// 返回: io::Result<Vec<SamRecord>> - 所有记录列表
/// 用途: 批量读取SAM文件用于统计分析
fn read_sam(path: &str) -> io::Result<Vec<SamRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = parse_sam_line(&line) {
            records.push(rec);
        }
    }
    Ok(records)
}

/// -----------------------------------------------------------
/// 解析SAM flag（位运算）
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   flag - SAM FLAG值（u16整数）
/// 返回: Python字典，每个标志位的布尔值
/// 算法: 使用位与运算（&）检查各个标志位
/// 用途: 理解read的比对状态（是否配对、是否mapped等）
#[pyfunction]
pub fn parse_flag(py: Python, flag: u16) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    dict.set_item("read_paired", flag & 0x1 != 0)?;          // 0x1: read是配对的
    dict.set_item("read_mapped_in_pair", flag & 0x2 != 0)?;  // 0x2: 正确配对
    dict.set_item("read_unmapped", flag & 0x4 != 0)?;        // 0x4: 未比对
    dict.set_item("mate_unmapped", flag & 0x8 != 0)?;        // 0x8: mate未比对
    dict.set_item("read_reverse", flag & 0x10 != 0)?;        // 0x10: 反向互补比对
    dict.set_item("mate_reverse", flag & 0x20 != 0)?;        // 0x20: mate反向互补
    dict.set_item("first_in_pair", flag & 0x40 != 0)?;       // 0x40: 第一个配对read
    dict.set_item("second_in_pair", flag & 0x80 != 0)?;      // 0x80: 第二个配对read
    dict.set_item("not_primary", flag & 0x100 != 0)?;        // 0x100: 非主要比对
    dict.set_item("read_fails_qc", flag & 0x200 != 0)?;      // 0x200: 未通过QC
    dict.set_item("read_is_duplicate", flag & 0x400 != 0)?;  // 0x400: 重复read
    dict.set_item("supplementary", flag & 0x800 != 0)?;      // 0x800: 补充比对
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// SAM文件统计
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   path - SAM文件路径
/// 返回: Python字典，包含全面的比对统计
/// 统计指标:
///   - total_reads: 总read数
///   - mapped: 已比对read数
///   - unmapped: 未比对read数
///   - paired: 配对read数
///   - duplicates: 重复read数
///   - primary: 主要比对数
///   - supplementary: 补充比对数
///   - avg_mapq: 平均比对质量
///   - chrom_counts: 各染色体read分布
/// 用途: 评估比对质量，检测偏倚
#[pyfunction]
pub fn sam_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total = records.len();
    // 使用位运算统计各种状态
    let mapped = records.iter().filter(|r| r.flag & 0x4 == 0).count();      // 未设置0x4位=已比对
    let unmapped = total - mapped;
    let paired = records.iter().filter(|r| r.flag & 0x1 != 0).count();      // 设置0x1位=配对
    let duplicates = records.iter().filter(|r| r.flag & 0x400 != 0).count(); // 设置0x400位=重复
    let primary = records.iter().filter(|r| r.flag & 0x100 == 0).count();    // 未设置0x100位=主要
    let supplementary = records.iter().filter(|r| r.flag & 0x800 != 0).count(); // 设置0x800位=补充

    // 统计染色体分布和平均MAPQ
    let mut chrom_counts: HashMap<String, usize> = HashMap::new();
    let mut mapq_sum: u64 = 0;
    let mut mapped_with_mapq = 0usize;
    for rec in &records {
        if rec.flag & 0x4 == 0 && rec.rname != "*" {  // 已比对且有参考序列
            *chrom_counts.entry(rec.rname.clone()).or_insert(0) += 1;
            mapq_sum += rec.mapq as u64;
            mapped_with_mapq += 1;
        }
    }

    let avg_mapq = if mapped_with_mapq > 0 {
        mapq_sum as f64 / mapped_with_mapq as f64
    } else {
        0.0
    };

    // 构建返回字典
    let dict = PyDict::new_bound(py);
    dict.set_item("total_reads", total)?;
    dict.set_item("mapped", mapped)?;
    dict.set_item("unmapped", unmapped)?;
    dict.set_item("paired", paired)?;
    dict.set_item("duplicates", duplicates)?;
    dict.set_item("primary", primary)?;
    dict.set_item("supplementary", supplementary)?;
    dict.set_item("avg_mapq", avg_mapq)?;
    dict.set_item("chrom_counts", chrom_counts)?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// 按mapping quality过滤read
/// -----------------------------------------------------------
/// 参数:
///   path     - SAM文件路径
///   min_mapq - 最小比对质量阈值
/// 返回: Vec<String> - 通过过滤的read名称列表
/// 用途: 去除低质量比对，提高下游分析准确性
#[pyfunction]
pub fn filter_by_mapq(path: &str, min_mapq: u8) -> PyResult<Vec<String>> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let filtered: Vec<String> = records
        .iter()
        .filter(|r| r.flag & 0x4 == 0 && r.mapq >= min_mapq)  // 已比对且MAPQ达标
        .map(|r| r.qname.clone())
        .collect();
    Ok(filtered)
}

/// -----------------------------------------------------------
/// 按染色体区间提取read
/// -----------------------------------------------------------
/// 参数:
///   path  - SAM文件路径
///   chrom - 染色体名称
///   start - 起始位置
///   end   - 结束位置
/// 返回: Vec<String> - 区间内的read信息（QNAME, POS, CIGAR）
/// 用途: 查看特定区域的比对情况
#[pyfunction]
pub fn fetch_region(path: &str, chrom: &str, start: i64, end: i64) -> PyResult<Vec<String>> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let hits: Vec<String> = records
        .iter()
        .filter(|r| {
            r.rname == chrom        // 匹配染色体
                && r.pos >= start   // 位置在区间内
                && r.pos <= end
                && r.flag & 0x4 == 0  // 已比对
        })
        .map(|r| format!("{}\t{}\t{}", r.qname, r.pos, r.cigar))
        .collect();
    Ok(hits)
}

/// -----------------------------------------------------------
/// 单个位置的coverage统计
/// -----------------------------------------------------------
/// 参数:
///   path     - SAM文件路径
///   chrom    - 染色体名称
///   position - 目标位置
/// 返回: 覆盖该位置的read数量
/// 算法:
///   检查每个read的起始位置和长度，判断是否覆盖目标位置
/// 用途: 特定位点的深度分析
#[pyfunction]
pub fn coverage_at_position(path: &str, chrom: &str, position: i64) -> PyResult<usize> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut count = 0usize;
    for rec in &records {
        if rec.rname != chrom || rec.flag & 0x4 != 0 {
            continue;  // 跳过不同染色体或未比对的read
        }
        let read_len = rec.seq.len() as i64;
        // 检查read是否覆盖目标位置
        if rec.pos <= position && rec.pos + read_len > position {
            count += 1;
        }
    }
    Ok(count)
}

/// -----------------------------------------------------------
/// 区间coverage（返回每个位置的coverage）
/// -----------------------------------------------------------
/// 参数:
///   path  - SAM文件路径
///   chrom - 染色体名称
///   start - 起始位置
///   end   - 结束位置
/// 返回: Vec<i64> - 每个位置的coverage值
/// 算法:
///   1. 初始化长度为(end-start)的数组
///   2. 对于每个read，计算与区间的重叠部分
///   3. 重叠部分的coverage+1
/// 用途: 生成coverage图、检测拷贝数变异
#[pyfunction]
pub fn region_coverage(
    path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> PyResult<Vec<i64>> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let size = (end - start) as usize;
    let mut cov = vec![0i64; size];  // 初始化coverage数组

    for rec in &records {
        if rec.rname != chrom || rec.flag & 0x4 != 0 {
            continue;  // 跳过不匹配的read
        }
        let read_start = rec.pos;
        let read_end = rec.pos + rec.seq.len() as i64;
        // 计算重叠区间
        let overlap_start = read_start.max(start);
        let overlap_end = read_end.min(end);
        if overlap_start < overlap_end {
            // 重叠部分coverage+1
            for pos in overlap_start..overlap_end {
                cov[(pos - start) as usize] += 1;
            }
        }
    }
    Ok(cov)
}

/// -----------------------------------------------------------
/// Insert size统计（paired-end测序）
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   path - SAM文件路径
/// 返回: Python字典，包含insert size统计信息
/// 统计指标:
///   - count: 有效paired-end read数
///   - mean: 平均insert size
///   - median: 中位insert size
///   - min / max: 最小/最大insert size
/// 用途: 评估文库质量、检测嵌合体读段
#[pyfunction]
pub fn insert_size_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut sizes: Vec<i64> = Vec::new();
    for rec in &records {
        // 筛选条件：配对、已比对、正向tlen
        if rec.flag & 0x1 != 0 && rec.flag & 0x4 == 0 && rec.tlen > 0 {
            sizes.push(rec.tlen.abs());  // 使用绝对值
        }
    }
    sizes.sort_unstable();  // 排序用于计算中位数

    // 计算统计指标
    let mean = if sizes.is_empty() {
        0.0
    } else {
        sizes.iter().sum::<i64>() as f64 / sizes.len() as f64
    };
    let median = if sizes.is_empty() {
        0.0
    } else {
        sizes[sizes.len() / 2] as f64
    };
    let min = sizes.first().copied().unwrap_or(0);
    let max = sizes.last().copied().unwrap_or(0);

    let dict = PyDict::new_bound(py);
    dict.set_item("count", sizes.len())?;
    dict.set_item("mean", mean)?;
    dict.set_item("median", median)?;
    dict.set_item("min", min)?;
    dict.set_item("max", max)?;
    Ok(dict.into())
}
