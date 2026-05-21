//! ============================================================
//! 模块10: 表观基因组学 (epigenomics)
//! ============================================================
//! 本模块提供表观基因组学数据的分析功能。
//! 包括：ATAC-seq Tn5插入位点、Peak read计数、
//! Fragment长度分布、核小体定位、ChIP-seq覆盖度、
//! CpG甲基化、CUT&Tag片段统计等。
//!
//! 表观基因组学研究DNA和组蛋白的化学修饰：
//! - ATAC-seq: 检测开放染色质区域
//! - ChIP-seq: 检测蛋白质-DNA结合位点
//! - Bisulfite-seq: 检测DNA甲基化
//! - CUT&Tag/CUT&RUN: 新一代染色质分析技术
//!
//! 关键概念：
//! - Tn5转座酶：ATAC-seq中用于标记开放染色质
//! - NFR (Nucleosome-Free Region)：无核小体区域
//! - Peak: 信号富集区域（如转录因子结合位点）
//! - FRiP (Fraction of Reads in Peaks)：Peak中reads比例
//! - CpG甲基化：CpG位点的胞嘧啶甲基化
//!
//! 设计原则：
//! - 支持多种表观基因组学技术
//! - 返回Python字典或列表
//! - 计算质量控制指标（如FRiP）
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use std::collections::HashMap;  // Rust 标准哈希表

/// -----------------------------------------------------------
/// Tn5插入位点计算（ATAC-seq）
/// -----------------------------------------------------------
/// 参数:
///   chrom      - 染色体名称
///   pos        - read比对位置
///   is_reverse - 是否为反向链
/// 返回: (染色体, 插入位点位置)
/// 算法:
///   Tn5转座酶在DNA开放区域插入接头
///   正向链：插入位点 = read起始位置
///   反向链：插入位点 = read起始位置 - 1
///   这是因为Tn5以二聚体形式插入，切割位点有4bp偏移
/// 用途: ATAC-seq开放染色质分析
#[pyfunction]
pub fn tn5_insertion_sites(
    chrom: &str,
    pos: i64,
    is_reverse: bool,
) -> (String, i64) {
    let insertion = if is_reverse {
        pos - 1  // 反向链：位置-1
    } else {
        pos      // 正向链：位置不变
    };
    (chrom.to_string(), insertion)
}

/// -----------------------------------------------------------
/// Peak reads统计
/// -----------------------------------------------------------
/// 参数:
///   read_positions - read位置列表 (chrom, pos)
///   peaks          - peak区间列表 (chrom, start, end)
/// 返回: HashMap<usize, i64> - peak索引 -> read计数
/// 算法:
///   对于每个read，检查它落在哪些peak中
///   一个read可能落在多个peak中（重叠peak）
/// 用途: 计算每个peak的read覆盖度
#[pyfunction]
pub fn count_reads_in_peaks(
    read_positions: Vec<(String, i64)>,
    peaks: Vec<(String, i64, i64)>,
) -> HashMap<usize, i64> {
    let mut counts: HashMap<usize, i64> = HashMap::new();
    for (chrom, pos) in &read_positions {
        for (i, (pc, ps, pe)) in peaks.iter().enumerate() {
            if chrom == pc && *pos >= *ps && *pos < *pe {
                *counts.entry(i).or_insert(0) += 1;
            }
        }
    }
    counts
}

/// -----------------------------------------------------------
/// Fragment长度分布统计
/// -----------------------------------------------------------
/// 参数: fragment_lengths - fragment长度列表
/// 返回: HashMap<i64, usize> - 长度 -> 数量
/// 用途:
///   fragment长度分布可以揭示：
///   - 核小体定位（约200bp周期性）
///   - 开放染色质（短fragment）
///   - 核小体包裹区域（长fragment）
#[pyfunction]
pub fn fragment_length_distribution(
    fragment_lengths: Vec<i64>,
) -> HashMap<i64, usize> {
    let mut dist: HashMap<i64, usize> = HashMap::new();
    for len in fragment_lengths {
        *dist.entry(len).or_insert(0) += 1;
    }
    dist
}

/// -----------------------------------------------------------
/// 核小体-free reads分类
/// -----------------------------------------------------------
/// 参数:
///   fragment_lengths - fragment长度列表
///   max_nfr_len      - NFR最大长度阈值（通常150bp）
/// 返回: (NFR数量, 核小体数量)
/// 算法:
///   - 长度≤max_nfr_len的fragment来自无核小体区域（NFR）
///   - 长度>max_nfr_len的fragment来自核小体包裹区域
/// 用途: 评估ATAC-seq文库质量
#[pyfunction]
pub fn classify_nucleosome_free(
    fragment_lengths: Vec<i64>,
    max_nfr_len: i64,
) -> (usize, usize) {
    let nfr = fragment_lengths.iter().filter(|&&l| l <= max_nfr_len).count();
    let nucleosomal = fragment_lengths.len() - nfr;
    (nfr, nucleosomal)
}

/// -----------------------------------------------------------
/// ChIP-seq peak覆盖度计算
/// -----------------------------------------------------------
/// 参数:
///   py            - Python解释器引用
///   read_positions - read位置列表
///   peak_chrom    - peak染色体
///   peak_start    - peak起始位置
///   peak_end      - peak结束位置
/// 返回: Python字典，包含覆盖度统计
/// 统计指标:
///   - total_reads: peak内总read数
///   - max_coverage: 最大覆盖度
///   - mean_coverage: 平均覆盖度
///   - peak_length: peak长度
/// 用途: 评估ChIP-seq peak质量
#[pyfunction]
pub fn chipseq_peak_coverage(
    py: Python,
    read_positions: Vec<(String, i64)>,
    peak_chrom: &str,
    peak_start: i64,
    peak_end: i64,
) -> PyResult<PyObject> {
    let size = (peak_end - peak_start) as usize;
    let mut coverage = vec![0i64; size];  // 初始化覆盖度数组

    // 统计每个位置的read数
    for (chrom, pos) in &read_positions {
        if chrom == peak_chrom && *pos >= peak_start && *pos < peak_end {
            let idx = (pos - peak_start) as usize;
            coverage[idx] += 1;
        }
    }

    // 计算统计指标
    let total: i64 = coverage.iter().sum();
    let max_cov = *coverage.iter().max().unwrap_or(&0);
    let mean_cov = if size > 0 {
        total as f64 / size as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("total_reads", total)?;
    dict.set_item("max_coverage", max_cov)?;
    dict.set_item("mean_coverage", mean_cov)?;
    dict.set_item("peak_length", size)?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// CpG甲基化水平聚合
/// -----------------------------------------------------------
/// 参数:
///   methylation_calls - 甲基化检测结果 (chrom, pos, is_methylated)
///   cpg_sites         - CpG位点列表 (chrom, pos)
///   max_distance      - 最大匹配距离
/// 返回: HashMap<(chrom, pos), (甲基化数, 总数)>
/// 算法:
///   1. 对于每个甲基化检测结果，找到最近的CpG位点
///   2. 聚合同一位点的多次检测
///   3. 计算甲基化比例 = 甲基化数 / 总数
/// 用途: 计算CpG位点的甲基化水平
#[pyfunction]
pub fn cpg_methylation_levels(
    methylation_calls: Vec<(String, i64, bool)>,
    cpg_sites: Vec<(String, i64)>,
    max_distance: i64,
) -> HashMap<(String, i64), (usize, usize)> {
    let mut site_counts: HashMap<(String, i64), (usize, usize)> = HashMap::new();

    for (call_chrom, call_pos, is_methylated) in &methylation_calls {
        for (site_chrom, site_pos) in &cpg_sites {
            if call_chrom == site_chrom && (call_pos - site_pos).abs() <= max_distance {
                let entry = site_counts
                    .entry((site_chrom.clone(), *site_pos))
                    .or_insert((0, 0));
                if *is_methylated {
                    entry.0 += 1;  // 甲基化计数
                }
                entry.1 += 1;  // 总计数
            }
        }
    }
    site_counts
}

/// -----------------------------------------------------------
/// Fragment与peak重叠检测
/// -----------------------------------------------------------
/// 参数:
///   fragments - fragment列表 (chrom, start, end)
///   peaks     - peak列表 (chrom, start, end)
/// 返回: Vec<Vec<usize>> - 每个fragment重叠的peak索引列表
/// 用途: 检查哪些fragment落在哪些peak中
#[pyfunction]
pub fn fragment_overlap_peaks(
    fragments: Vec<(String, i64, i64)>,
    peaks: Vec<(String, i64, i64)>,
) -> Vec<Vec<usize>> {
    fragments
        .iter()
        .map(|(fc, fs, fe)| {
            peaks
                .iter()
                .enumerate()
                .filter(|(_, (pc, ps, pe))| fc == pc && *fs < *pe && *ps < *fe)
                .map(|(i, _)| i)
                .collect()
        })
        .collect()
}

/// -----------------------------------------------------------
/// 基因组分箱覆盖度
/// -----------------------------------------------------------
/// 参数:
///   read_positions - read位置列表
///   bins           - 分箱区间列表
/// 返回: Vec<i64> - 每个分箱的read数
/// 用途: 计算基因组分bin的覆盖度，用于CNV分析
#[pyfunction]
pub fn coverage_over_bins(
    read_positions: Vec<(String, i64)>,
    bins: Vec<(String, i64, i64)>,
) -> Vec<i64> {
    bins.iter()
        .map(|(chrom, start, end)| {
            read_positions
                .iter()
                .filter(|(rc, rp)| rc == chrom && rp >= start && rp < end)
                .count() as i64
        })
        .collect()
}

/// -----------------------------------------------------------
/// CUT&Tag/CUT&RUN片段统计
/// -----------------------------------------------------------
/// 参数:
///   py        - Python解释器引用
///   fragments - fragment列表 (chrom, start, end)
///   peaks     - peak列表 (chrom, start, end)
/// 返回: Python字典，包含片段统计
/// 统计指标:
///   - total_fragments: 总片段数
///   - fragments_in_peaks: peak内片段数
///   - frip: Fraction of Reads in Peaks
///   - median_fragment_length: 中位片段长度
/// 用途: 评估CUT&Tag实验质量
#[pyfunction]
pub fn cuttag_fragment_stats(
    py: Python,
    fragments: Vec<(String, i64, i64)>,
    peaks: Vec<(String, i64, i64)>,
) -> PyResult<PyObject> {
    let total_fragments = fragments.len();
    let mut in_peaks = 0usize;
    let mut fragment_lengths: Vec<i64> = Vec::new();

    // 统计片段长度和peak内片段数
    for (fc, fs, fe) in &fragments {
        let len = fe - fs;
        fragment_lengths.push(len);
        for (pc, ps, pe) in &peaks {
            if fc == pc && *fs < *pe && *ps < *fe {
                in_peaks += 1;
                break;  // 一个片段只计数一次
            }
        }
    }

    // 计算中位数
    fragment_lengths.sort_unstable();
    let median_len = if fragment_lengths.is_empty() {
        0
    } else {
        fragment_lengths[fragment_lengths.len() / 2]
    };
    // FRiP = Fraction of Reads in Peaks
    let frip = if total_fragments > 0 {
        in_peaks as f64 / total_fragments as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("total_fragments", total_fragments)?;
    dict.set_item("fragments_in_peaks", in_peaks)?;
    dict.set_item("frip", frip)?;
    dict.set_item("median_fragment_length", median_len)?;
    Ok(dict.into())
}
