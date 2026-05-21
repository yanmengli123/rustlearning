//! ============================================================
//! 模块9: 单细胞RNA-seq分析 (single_cell)
//! ============================================================
//! 本模块提供单细胞RNA-seq（scRNA-seq）数据处理功能。
//! 包括：barcode提取、UMI提取、barcode纠错、barcode rescue、
//! UMI去重、稀疏矩阵构建、Matrix Market格式输出、
//! 细胞QC指标、抗体标签去重、CRISPR guide计数等。
//!
//! 单细胞RNA-seq基本流程：
//!   1. 测序得到包含cell barcode和UMI的reads
//!   2. 提取barcode和UMI
//!   3. 纠正barcode错误（白名单匹配）
//!   4. 比对到参考基因组
//!   5. UMI去重
//!   6. 构建细胞×基因计数矩阵
//!   7. 下游分析（聚类、降维等）
//!
//! 关键概念：
//! - Cell barcode: 标记单个细胞的短序列（如10x Genomics的16bp）
//! - UMI (Unique Molecular Identifier): 标记原始分子的短序列
//! - 白名单: 已知的有效barcode列表
//! - 稀疏矩阵: 大部分元素为0的矩阵，高效存储
//! - Matrix Market (.mtx): 稀疏矩阵的标准存储格式
//!
//! 设计原则：
//! - 使用Levenshtein距离进行barcode纠错
//! - 使用Hamming距离进行barcode rescue
//! - 支持稀疏矩阵和MTX格式输出
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use std::collections::HashMap;  // Rust 标准哈希表

/// -----------------------------------------------------------
/// Cell barcode提取
/// -----------------------------------------------------------
/// 参数:
///   read         - 测序read序列
///   barcode_len  - barcode长度（如10x Genomics为16bp）
///   barcode_start - barcode起始位置（通常为0）
/// 返回: Option<String> - 提取的barcode，越界返回None
/// 用途: 从read中提取cell barcode
#[pyfunction]
pub fn extract_barcode(read: &str, barcode_len: usize, barcode_start: usize) -> Option<String> {
    if barcode_start + barcode_len > read.len() {
        return None;  // 越界
    }
    Some(read[barcode_start..barcode_start + barcode_len].to_string())
}

/// -----------------------------------------------------------
/// UMI提取
/// -----------------------------------------------------------
/// 参数:
///   read      - 测序read序列
///   umi_len   - UMI长度（如10x Genomics为12bp）
///   umi_start - UMI起始位置
/// 返回: Option<String> - 提取的UMI
/// 用途: 从read中提取UMI用于去重
#[pyfunction]
pub fn extract_umi(read: &str, umi_len: usize, umi_start: usize) -> Option<String> {
    if umi_start + umi_len > read.len() {
        return None;  // 越界
    }
    Some(read[umi_start..umi_start + umi_len].to_string())
}

/// -----------------------------------------------------------
/// Barcode纠错（白名单匹配）
/// -----------------------------------------------------------
/// 参数:
///   barcode      - 测序得到的barcode（可能有错误）
///   whitelist    - 白名单barcode列表
///   max_distance - 最大编辑距离阈值
/// 返回: Option<String> - 纠正后的barcode，无法纠正返回None
/// 算法:
///   1. 计算barcode与所有白名单barcode的Levenshtein距离
///   2. 选择距离最小且≤max_distance的白名单barcode
///   3. 如果精确匹配（距离=0），立即返回
/// 用途: 纠正测序错误导致的barcode错误
#[pyfunction]
pub fn correct_barcode(
    barcode: &str,
    whitelist: Vec<String>,
    max_distance: usize,
) -> Option<String> {
    let mut best_match = None;
    let mut best_dist = max_distance + 1;

    for wb in &whitelist {
        let dist = super::alignment::levenshtein_distance(barcode, wb);
        if dist < best_dist {
            best_dist = dist;
            best_match = Some(wb.clone());
        }
        if dist == 0 {
            return best_match;  // 精确匹配，立即返回
        }
    }
    // 检查是否在阈值内
    if best_dist <= max_distance {
        best_match
    } else {
        None
    }
}

/// -----------------------------------------------------------
/// Barcode rescue（Hamming距离）
/// -----------------------------------------------------------
/// 参数:
///   barcode     - 测序barcode
///   whitelist   - 白名单barcode列表
///   max_hamming - 最大Hamming距离阈值
/// 返回: Vec<(String, usize)> - 匹配的barcode及距离
/// 算法:
///   使用Hamming距离（只允许替换，不允许插入/删除）
///   要求barcode和白名单barcode长度相同
/// 用途: 当Levenshtein纠错失败时的备选方案
#[pyfunction]
pub fn barcode_rescue(
    barcode: &str,
    whitelist: Vec<String>,
    max_hamming: usize,
) -> Vec<(String, usize)> {
    let mut matches = Vec::new();
    for wb in &whitelist {
        if barcode.len() != wb.len() {
            continue;  // 长度不同，跳过
        }
        let dist = super::alignment::hamming_distance(barcode, wb).unwrap_or(999);
        if dist <= max_hamming {
            matches.push((wb.clone(), dist));
        }
    }
    matches
}

/// -----------------------------------------------------------
/// UMI去重
/// -----------------------------------------------------------
/// 参数:
///   umis              - UMI列表
///   max_edit_distance - 最大编辑距离阈值
/// 返回: Vec<(String, Vec<usize>)> - (代表UMI, 索引列表)
/// 说明: 委托给rnaseq模块的umi_collapse函数
/// 用途: 合并测序错误导致的UMI重复
#[pyfunction]
pub fn umi_dedup(
    umis: Vec<String>,
    max_edit_distance: usize,
) -> Vec<(String, Vec<usize>)> {
    super::rnaseq::umi_collapse(umis, max_edit_distance)
}

/// -----------------------------------------------------------
/// 稀疏矩阵条目结构体
/// -----------------------------------------------------------
/// 存储稀疏矩阵的非零元素
/// 使用(row, col, value)三元组表示
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SparseEntry {
    pub row: usize,   // 行索引（特征/基因）
    pub col: usize,   // 列索引（细胞/barcode）
    pub value: f64,   // 计数值
}

/// -----------------------------------------------------------
/// Feature-barcode矩阵构建
/// -----------------------------------------------------------
/// 参数:
///   py       - Python解释器引用
///   features - 特征名称列表（基因名）
///   barcodes - 细胞barcode列表
///   entries  - 稀疏矩阵条目 (row, col, value)
/// 返回: Python字典，包含矩阵信息
/// 输出格式:
///   - features: 特征名列表
///   - barcodes: barcode列表
///   - entries: 稀疏条目列表
///   - n_features: 特征数
///   - n_barcodes: barcode数
/// 用途: 构建scRNA-seq计数矩阵
#[pyfunction]
pub fn build_feature_barcode_matrix(
    py: Python,
    features: Vec<String>,
    barcodes: Vec<String>,
    entries: Vec<(usize, usize, f64)>,
) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    dict.set_item("features", &features)?;
    dict.set_item("barcodes", &barcodes)?;

    // 将条目转为Python字典列表
    let sparse_data: Vec<PyObject> = entries
        .iter()
        .map(|(r, c, v)| {
            let d = PyDict::new_bound(py);
            d.set_item("row", *r).unwrap();
            d.set_item("col", *c).unwrap();
            d.set_item("value", *v).unwrap();
            d.into()
        })
        .collect();

    dict.set_item("entries", sparse_data)?;
    dict.set_item("n_features", features.len())?;
    dict.set_item("n_barcodes", barcodes.len())?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// Matrix Market格式输出
/// -----------------------------------------------------------
/// 参数:
///   n_rows   - 行数（特征数）
///   n_cols   - 列数（细胞数）
///   entries  - 稀疏矩阵条目 (row, col, value)
/// 返回: String - Matrix Market格式文本
/// 格式说明:
///   第1行: %%MatrixMarket matrix coordinate real general
///   第2行: 行数 列数 非零元素数
///   后续行: 行索引 列索引 值（1-based索引）
/// 用途: 生成标准MTX格式文件，用于Scanr/Seurat等工具
#[pyfunction]
pub fn to_matrix_market(
    n_rows: usize,
    n_cols: usize,
    entries: Vec<(usize, usize, f64)>,
) -> String {
    let mut lines = vec![
        "%%MatrixMarket matrix coordinate real general".to_string(),
        format!("{} {} {}", n_rows, n_cols, entries.len()),
    ];
    // 注意：Matrix Market使用1-based索引
    for (row, col, val) in entries {
        lines.push(format!("{} {} {}", row + 1, col + 1, val));
    }
    lines.join("\n")
}

/// -----------------------------------------------------------
/// 细胞QC指标计算
/// -----------------------------------------------------------
/// 参数:
///   py              - Python解释器引用
///   counts_per_cell - 每个细胞的总UMI计数
///   genes_per_cell  - 每个细胞检测到的基因数
/// 返回: Python字典，包含QC指标
/// QC指标:
///   - n_cells: 细胞总数
///   - total_counts: 总UMI数
///   - total_genes: 总基因数
///   - avg_counts_per_cell: 平均UMI数/细胞
///   - avg_genes_per_cell: 平均基因数/细胞
///   - max/min_counts: 最大/最小UMI数
/// 用途: 评估单细胞数据质量，过滤低质量细胞
#[pyfunction]
pub fn cell_qc(
    py: Python,
    counts_per_cell: Vec<usize>,
    genes_per_cell: Vec<usize>,
) -> PyResult<PyObject> {
    let total_counts: usize = counts_per_cell.iter().sum();
    let total_genes: usize = genes_per_cell.iter().sum();
    let n_cells = counts_per_cell.len();
    let avg_counts = if n_cells > 0 {
        total_counts as f64 / n_cells as f64
    } else {
        0.0
    };
    let avg_genes = if n_cells > 0 {
        total_genes as f64 / n_cells as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("n_cells", n_cells)?;
    dict.set_item("total_counts", total_counts)?;
    dict.set_item("total_genes", total_genes)?;
    dict.set_item("avg_counts_per_cell", avg_counts)?;
    dict.set_item("avg_genes_per_cell", avg_genes)?;
    dict.set_item("max_counts", counts_per_cell.iter().max().unwrap_or(&0))?;
    dict.set_item("min_counts", counts_per_cell.iter().min().unwrap_or(&0))?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// Antibody/hashtag标签去重
/// -----------------------------------------------------------
/// 参数:
///   tag_counts - 每个细胞的标签计数矩阵
///   tag_names  - 标签名称列表
///   threshold  - 阈值（低于此值视为阴性）
/// 返回: Vec<String> - 每个细胞的标签分配结果
/// 算法:
///   1. 对于每个细胞，找出计数最高的标签
///   2. 如果最高计数>threshold，则分配该标签
///   3. 否则标记为"unknown"
/// 用途: CITE-seq/Cell Hashing样本去重
#[pyfunction]
pub fn demux_tags(
    tag_counts: Vec<Vec<f64>>,
    tag_names: Vec<String>,
    threshold: f64,
) -> Vec<String> {
    let mut assignments = Vec::new();
    for cell_tags in &tag_counts {
        let mut max_tag = "unknown".to_string();
        let mut max_val = 0.0f64;
        for (i, &val) in cell_tags.iter().enumerate() {
            if val > max_val && val > threshold {
                max_val = val;
                if i < tag_names.len() {
                    max_tag = tag_names[i].clone();
                }
            }
        }
        assignments.push(max_tag);
    }
    assignments
}

/// -----------------------------------------------------------
/// CRISPR guide计数
/// -----------------------------------------------------------
/// 参数: guide_assignments - guide分配列表
/// 返回: HashMap<String, usize> - guide -> 细胞数
/// 用途: Perturb-seq实验中统计每个guide的细胞数
#[pyfunction]
pub fn count_guides(
    guide_assignments: Vec<String>,
) -> HashMap<String, usize> {
    let mut counts: HashMap<String, usize> = HashMap::new();
    for guide in guide_assignments {
        *counts.entry(guide).or_insert(0) += 1;
    }
    counts
}
