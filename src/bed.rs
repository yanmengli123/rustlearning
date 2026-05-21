//! ============================================================
//! 模块6: 基因组区间操作 (bed)
//! ============================================================
//! 本模块提供基因组区间文件的解析和操作功能。
//! 包括：BED格式解析和统计、区间交集/并集/差集、
//! VCF格式解析和统计、GTF格式解析和统计等。
//!
//! BED格式说明（Browser Extensible Data）：
//!   用于定义基因组区间的标准格式
//!   必选字段（最少3列）：
//!     chrom  - 染色体名称
//!     start  - 起始位置（0-based）
//!     end    - 结束位置（不包含）
//!   可选字段：
//!     name   - 区间名称
//!     score  - 分数
//!     strand - 链方向（+/-/.）
//!
//! VCF格式说明（Variant Call Format）：
//!   存储基因组变异信息的标准格式
//!   包含变异位置、参考/替代等位基因、质量分数等
//!
//! GTF格式说明（Gene Transfer Format）：
//!   描述基因结构（基因、转录本、外显子等）
//!   广泛用于RNA-seq分析
//!
//! 设计原则：
//! - 支持BED/VCF/GTF三种格式
//! - 区间操作使用0-based坐标系
//! - 自动跳过注释行（#开头）和空行
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use std::collections::HashMap;  // Rust 标准哈希表
use std::fs::File;              // 文件操作
use std::io::{self, BufRead, BufReader};  // IO操作

/// -----------------------------------------------------------
/// BED区间结构体
/// -----------------------------------------------------------
/// 存储BED格式的区间信息
/// 使用0-based坐标系（start包含，end不包含）
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BedInterval {
    pub chrom: String,   // 染色体名称
    pub start: i64,      // 起始位置（0-based）
    pub end: i64,        // 结束位置（不包含）
    pub name: String,    // 区间名称
    pub score: f64,      // 分数
    pub strand: char,    // 链方向（+/-/.）
}

/// -----------------------------------------------------------
/// 解析BED行
/// -----------------------------------------------------------
/// 参数: line - BED格式的一行文本
/// 返回: Option<BedInterval> - 解析成功返回Some
/// 算法:
///   1. 跳过注释行（#）和空行
///   2. 按制表符分割字段
///   3. 至少需要3个字段（chrom, start, end）
///   4. 可选字段使用默认值填充
pub fn parse_bed_line(line: &str) -> Option<BedInterval> {
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 3 {
        return None;  // 字段不足
    }
    Some(BedInterval {
        chrom: fields[0].to_string(),
        start: fields[1].parse().unwrap_or(0),   // 起始位置
        end: fields[2].parse().unwrap_or(0),     // 结束位置
        name: fields.get(3).unwrap_or(&".").to_string(),   // 默认"."
        score: fields.get(4).unwrap_or(&"0").parse().unwrap_or(0.0),  // 默认0
        strand: fields.get(5).unwrap_or(&".").chars().next().unwrap_or('.'),  // 默认"."
    })
}

/// -----------------------------------------------------------
/// 读取BED文件
/// -----------------------------------------------------------
/// 参数: path - BED文件路径
/// 返回: io::Result<Vec<BedInterval>> - 区间列表
pub fn read_bed(path: &str) -> io::Result<Vec<BedInterval>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(iv) = parse_bed_line(&line) {
            intervals.push(iv);
        }
    }
    Ok(intervals)
}

/// -----------------------------------------------------------
/// BED文件统计
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   path - BED文件路径
/// 返回: Python字典，包含区间统计信息
/// 统计指标:
///   - total_intervals: 总区间数
///   - total_bases: 总碱基数（区间长度之和）
///   - avg_length: 平均区间长度
///   - min_length / max_length: 最短/最长区间
///   - chrom_counts: 各染色体区间分布
/// 用途: 评估BED文件内容、检查区间分布
#[pyfunction]
pub fn bed_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let intervals = read_bed(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total = intervals.len();
    let mut chrom_counts: HashMap<String, usize> = HashMap::new();
    let mut total_bases: i64 = 0;
    let mut min_len = i64::MAX;
    let mut max_len = 0i64;

    for iv in &intervals {
        *chrom_counts.entry(iv.chrom.clone()).or_insert(0) += 1;
        let len = iv.end - iv.start;  // 区间长度
        total_bases += len;
        min_len = min_len.min(len);
        max_len = max_len.max(len);
    }

    let avg_len = if total > 0 {
        total_bases as f64 / total as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("total_intervals", total)?;
    dict.set_item("total_bases", total_bases)?;
    dict.set_item("avg_length", avg_len)?;
    dict.set_item("min_length", if min_len == i64::MAX { 0 } else { min_len })?;
    dict.set_item("max_length", max_len)?;
    dict.set_item("chrom_counts", chrom_counts)?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// 区间交集（Intersection）
/// -----------------------------------------------------------
/// 参数:
///   path_a - BED文件A路径
///   path_b - BED文件B路径
/// 返回: Vec<String> - 交集结果（制表符分隔）
/// 算法:
///   1. 读取两个BED文件
///   2. 对于A中的每个区间，检查与B中所有区间的重叠
///   3. 重叠条件：同染色体且 overlap_start < overlap_end
///   4. 输出重叠区间信息
/// 用途: 寻找两个基因组特征的共同区域
#[pyfunction]
pub fn bed_intersect(path_a: &str, path_b: &str) -> PyResult<Vec<String>> {
    let a = read_bed(path_a)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let b = read_bed(path_b)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut results = Vec::new();
    for iv_a in &a {
        for iv_b in &b {
            if iv_a.chrom == iv_b.chrom {
                // 计算重叠区间
                let overlap_start = iv_a.start.max(iv_b.start);
                let overlap_end = iv_a.end.min(iv_b.end);
                if overlap_start < overlap_end {
                    // 有重叠
                    results.push(format!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        iv_a.chrom, overlap_start, overlap_end, iv_a.name, iv_b.name,
                        overlap_end - overlap_start  // 重叠长度
                    ));
                }
            }
        }
    }
    Ok(results)
}

/// -----------------------------------------------------------
/// 区间并集（Merge）
/// -----------------------------------------------------------
/// 参数: path - BED文件路径
/// 返回: Vec<String> - 合并后的区间
/// 算法:
///   1. 按染色体和起始位置排序
///   2. 遍历区间，如果当前区间与上一个重叠，则合并
///   3. 合并条件：同染色体且 last.end >= current.start
/// 用途: 合并重叠区间，生成连续区域
#[pyfunction]
pub fn bed_merge(path: &str) -> PyResult<Vec<String>> {
    let mut intervals = read_bed(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // 排序：先按染色体，再按起始位置，最后按结束位置
    intervals.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    // 合并重叠区间
    let mut merged: Vec<BedInterval> = Vec::new();
    for iv in intervals {
        if let Some(last) = merged.last_mut() {
            if last.chrom == iv.chrom && last.end >= iv.start {
                // 重叠，合并（扩展结束位置）
                last.end = last.end.max(iv.end);
                continue;
            }
        }
        merged.push(iv);  // 不重叠，添加新区间
    }

    // 格式化输出
    let result: Vec<String> = merged
        .iter()
        .map(|iv| format!("{}\t{}\t{}", iv.chrom, iv.start, iv.end))
        .collect();
    Ok(result)
}

/// -----------------------------------------------------------
/// 区间差集（Subtract）
/// -----------------------------------------------------------
/// 参数:
///   path_a - BED文件A路径（被减区间）
///   path_b - BED文件B路径（减去区间）
/// 返回: Vec<String> - A中不与B重叠的部分
/// 算法:
///   1. 对于A中的每个区间，遍历B中所有区间
///   2. 如果B区间与A区间重叠，从A中移除重叠部分
///   3. 处理可能的多次切割（B区间可能在A内部）
/// 用途: 提取特有区间、去除已知区域
#[pyfunction]
pub fn bed_subtract(path_a: &str, path_b: &str) -> PyResult<Vec<String>> {
    let a = read_bed(path_a)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let b = read_bed(path_b)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut results = Vec::new();
    for iv_a in &a {
        // remaining存储A区间被切割后的剩余部分
        let mut remaining: Vec<(i64, i64)> = vec![(iv_a.start, iv_a.end)];
        for iv_b in &b {
            if iv_a.chrom != iv_b.chrom {
                continue;  // 不同染色体，跳过
            }
            let mut new_remaining = Vec::new();
            for (s, e) in remaining {
                if iv_b.end <= s || iv_b.start >= e {
                    // B区间在A外部，A区间完整保留
                    new_remaining.push((s, e));
                } else {
                    // B区间与A重叠，切割A区间
                    if iv_b.start > s {
                        new_remaining.push((s, iv_b.start));  // 左侧剩余
                    }
                    if iv_b.end < e {
                        new_remaining.push((iv_b.end, e));    // 右侧剩余
                    }
                }
            }
            remaining = new_remaining;
        }
        // 输出所有剩余区间
        for (s, e) in remaining {
            results.push(format!("{}\t{}\t{}", iv_a.chrom, s, e));
        }
    }
    Ok(results)
}

/// -----------------------------------------------------------
/// 最近特征查找（Closest）
/// -----------------------------------------------------------
/// 参数:
///   path_a - BED文件A路径（查询区间）
///   path_b - BED文件B路径（参考区间）
/// 返回: Vec<String> - 每个A区间最近的B区间信息
/// 算法:
///   1. 对于A中的每个区间，遍历B中所有区间
///   2. 计算距离（不重叠时为间隔距离，重叠时为0）
///   3. 选择距离最小的B区间
/// 用途: 寻找最近的基因、调控元件等
#[pyfunction]
pub fn bed_closest(path_a: &str, path_b: &str) -> PyResult<Vec<String>> {
    let a = read_bed(path_a)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let b = read_bed(path_b)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut results = Vec::new();
    for iv_a in &a {
        let mut min_dist = i64::MAX;
        let mut closest_name = String::new();
        for iv_b in &b {
            if iv_a.chrom != iv_b.chrom {
                continue;  // 不同染色体，跳过
            }
            // 计算距离
            let dist = if iv_a.end <= iv_b.start {
                iv_b.start - iv_a.end  // A在B左侧
            } else if iv_b.end <= iv_a.start {
                iv_a.start - iv_b.end  // A在B右侧
            } else {
                0  // 重叠
            };
            if dist < min_dist {
                min_dist = dist;
                closest_name = iv_b.name.clone();
            }
        }
        results.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            iv_a.chrom, iv_a.start, iv_a.end, iv_a.name, closest_name, min_dist
        ));
    }
    Ok(results)
}

/// -----------------------------------------------------------
/// VCF记录结构体
/// -----------------------------------------------------------
/// 存储VCF格式的变异信息
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,             // 染色体
    pub pos: i64,                  // 位置（1-based）
    pub id: String,                // 变异ID
    pub ref_allele: String,        // 参考等位基因
    pub alt_allele: String,        // 替代等位基因
    pub qual: f64,                 // 质量分数
    pub filter: String,            // 过滤状态
    pub info: HashMap<String, String>,  // INFO字段
}

/// -----------------------------------------------------------
/// 解析VCF行
/// -----------------------------------------------------------
/// 参数: line - VCF格式的一行
/// 返回: Option<VcfRecord>
pub fn parse_vcf_line(line: &str) -> Option<VcfRecord> {
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return None;
    }
    // 解析INFO字段（key=value格式，分号分隔）
    let mut info = HashMap::new();
    if fields[7] != "." {
        for item in fields[7].split(';') {
            if let Some(idx) = item.find('=') {
                info.insert(item[..idx].to_string(), item[idx + 1..].to_string());
            } else {
                info.insert(item.to_string(), "true".to_string());  // 标志位
            }
        }
    }
    Some(VcfRecord {
        chrom: fields[0].to_string(),
        pos: fields[1].parse().unwrap_or(0),
        id: fields[2].to_string(),
        ref_allele: fields[3].to_string(),
        alt_allele: fields[4].to_string(),
        qual: fields[5].parse().unwrap_or(0.0),
        filter: fields[6].to_string(),
        info,
    })
}

/// -----------------------------------------------------------
/// 读取VCF文件
/// -----------------------------------------------------------
/// 参数: path - VCF文件路径
/// 返回: io::Result<Vec<VcfRecord>>
pub fn read_vcf(path: &str) -> io::Result<Vec<VcfRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = parse_vcf_line(&line) {
            records.push(rec);
        }
    }
    Ok(records)
}

/// -----------------------------------------------------------
/// VCF统计
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   path - VCF文件路径
/// 返回: Python字典，包含变异统计
/// 统计指标:
///   - total_variants: 总变异数
///   - snps: SNP数量（单核苷酸多态性）
///   - indels: InDel数量（插入/缺失）
///   - chrom_counts: 各染色体变异分布
/// 用途: 变异数据概览、质量评估
#[pyfunction]
pub fn vcf_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let records = read_vcf(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total = records.len();
    let mut snps = 0usize;
    let mut indels = 0usize;
    let mut chrom_counts: HashMap<String, usize> = HashMap::new();

    for rec in &records {
        *chrom_counts.entry(rec.chrom.clone()).or_insert(0) += 1;
        // SNP: 参考和替代等位基因都是单碱基
        if rec.ref_allele.len() == 1 && rec.alt_allele.len() == 1 {
            snps += 1;
        } else {
            indels += 1;  // 插入或缺失
        }
    }

    let dict = PyDict::new_bound(py);
    dict.set_item("total_variants", total)?;
    dict.set_item("snps", snps)?;
    dict.set_item("indels", indels)?;
    dict.set_item("chrom_counts", chrom_counts)?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// VCF质量过滤
/// -----------------------------------------------------------
/// 参数:
///   path     - VCF文件路径
///   min_qual - 最小质量分数阈值
/// 返回: Vec<String> - 通过过滤的变异记录
/// 用途: 去除低质量变异调用
#[pyfunction]
pub fn vcf_filter(path: &str, min_qual: f64) -> PyResult<Vec<String>> {
    let records = read_vcf(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let filtered: Vec<String> = records
        .iter()
        .filter(|r| r.qual >= min_qual)  // 质量过滤
        .map(|r| {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{:.1}",
                r.chrom, r.pos, r.id, r.ref_allele, r.alt_allele, r.qual
            )
        })
        .collect();
    Ok(filtered)
}

/// -----------------------------------------------------------
/// GTF记录结构体
/// -----------------------------------------------------------
/// 存储GTF格式的基因结构信息
/// 包含基因、转录本、外显子等特征
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GtfRecord {
    pub chrom: String,                    // 染色体
    pub source: String,                   // 来源（如Ensembl）
    pub feature: String,                  // 特征类型（gene/transcript/exon等）
    pub start: i64,                       // 起始位置（1-based）
    pub end: i64,                         // 结束位置
    pub score: String,                    // 分数
    pub strand: char,                     // 链方向
    pub frame: String,                    // 阅读框（0/1/2）
    pub attributes: HashMap<String, String>,  // 属性（gene_id, transcript_id等）
}

/// -----------------------------------------------------------
/// 解析GTF行
/// -----------------------------------------------------------
/// 参数: line - GTF格式的一行
/// 返回: Option<GtfRecord>
/// 算法:
///   1. 解析9个制表符分隔的字段
///   2. 第9个字段是属性，格式为：key "value"; key "value";
pub fn parse_gtf_line(line: &str) -> Option<GtfRecord> {
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }
    // 解析属性字段
    let mut attrs = HashMap::new();
    for item in fields[8].split(';') {
        let item = item.trim();
        if item.is_empty() {
            continue;
        }
        let parts: Vec<&str> = item.splitn(2, ' ').collect();
        if parts.len() == 2 {
            let val = parts[1].trim_matches('"');  // 去除引号
            attrs.insert(parts[0].to_string(), val.to_string());
        }
    }
    Some(GtfRecord {
        chrom: fields[0].to_string(),
        source: fields[1].to_string(),
        feature: fields[2].to_string(),
        start: fields[3].parse().unwrap_or(0),
        end: fields[4].parse().unwrap_or(0),
        score: fields[5].to_string(),
        strand: fields[6].chars().next().unwrap_or('.'),
        frame: fields[7].to_string(),
        attributes: attrs,
    })
}

/// -----------------------------------------------------------
/// GTF统计
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   path - GTF文件路径
/// 返回: Python字典，包含基因结构统计
/// 统计指标:
///   - gene_count: 基因数
///   - transcript_count: 转录本数
///   - feature_counts: 各特征类型计数
/// 用途: 了解注释文件结构、验证完整性
#[pyfunction]
pub fn gtf_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let file = File::open(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let mut feature_counts: HashMap<String, usize> = HashMap::new();
    let mut gene_count = 0usize;
    let mut transcript_count = 0usize;

    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = parse_gtf_line(&line) {
            *feature_counts.entry(rec.feature.clone()).or_insert(0) += 1;
            match rec.feature.as_str() {
                "gene" => gene_count += 1,
                "transcript" => transcript_count += 1,
                _ => {}
            }
        }
    }

    let dict = PyDict::new_bound(py);
    dict.set_item("gene_count", gene_count)?;
    dict.set_item("transcript_count", transcript_count)?;
    dict.set_item("feature_counts", feature_counts)?;
    Ok(dict.into())
}
