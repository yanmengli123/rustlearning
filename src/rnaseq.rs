//! ============================================================
//! 模块8: RNA-seq分析 (rnaseq)
//! ============================================================
//! 本模块提供RNA-seq数据处理和分析功能。
//! 包括：基因计数矩阵构建、外显子计数、剪接位点统计、
//! 内含子保留检测、UMI去重、基因类型统计等。
//!
//! RNA-seq基本流程：
//!   1. 测序得到FASTQ文件
//!   2. 比对到参考基因组（产生SAM/BAM文件）
//!   3. 基因/转录本定量（产生计数矩阵）
//!   4. 差异表达分析
//!
//! 关键概念：
//! - GTF文件：描述基因结构（外显子、内含子、转录本等）
//! - 基因计数矩阵：每个基因在每个样本中的read计数
//! - UMI（Unique Molecular Identifier）：消除PCR扩增偏差
//! - 剪接位点：外显子-内含子边界
//! - 内含子保留：read覆盖了整个内含子区域
//!
//! 设计原则：
//! - 使用GTF注释文件定义基因区间
//! - 支持UMI去重（Levenshtein距离）
//! - 返回Python字典或列表
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use std::collections::HashMap;  // Rust 标准哈希表
use std::fs::File;              // 文件操作
use std::io::{self, BufRead, BufReader};  // IO操作

/// -----------------------------------------------------------
/// 从GTF构建基因区间索引
/// -----------------------------------------------------------
/// 参数: path - GTF注释文件路径
/// 返回: HashMap<染色体, Vec<(chrom, start, end, gene_id)>>
/// 算法:
///   1. 解析GTF文件，提取feature=="gene"的行
///   2. 获取gene_id属性
///   3. 按染色体分组存储基因区间
/// 用途: 快速查找特定染色体上的基因
pub fn build_gene_index(path: &str) -> io::Result<HashMap<String, Vec<(String, i64, i64, String)>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut genes: HashMap<String, Vec<(String, i64, i64, String)>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == "gene" {
                // 获取gene_id属性
                let gene_id = rec
                    .attributes
                    .get("gene_id")
                    .cloned()
                    .unwrap_or_default();
                // 按染色体分组
                genes
                    .entry(rec.chrom.clone())
                    .or_default()
                    .push((rec.chrom, rec.start, rec.end, gene_id));
            }
        }
    }
    Ok(genes)
}

/// -----------------------------------------------------------
/// 基因计数矩阵构建
/// -----------------------------------------------------------
/// 参数:
///   gtf_path       - GTF注释文件路径
///   read_positions - read位置列表 (chrom, pos)
/// 返回: HashMap<String, i64> - gene_id -> read计数
/// 算法:
///   1. 构建基因区间索引
///   2. 对于每个read，查找它落在哪个基因区间内
///   3. 累加基因计数
/// 用途: RNA-seq基因定量（生成计数矩阵）
#[pyfunction]
pub fn gene_count_matrix(
    gtf_path: &str,
    read_positions: Vec<(String, i64)>,
) -> PyResult<HashMap<String, i64>> {
    let gene_index = build_gene_index(gtf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut counts: HashMap<String, i64> = HashMap::new();

    // 对每个read查找所属基因
    for (chrom, pos) in &read_positions {
        if let Some(genes) = gene_index.get(chrom) {
            for (_, start, end, gene_id) in genes {
                if pos >= start && pos < end {
                    *counts.entry(gene_id.clone()).or_insert(0) += 1;
                }
            }
        }
    }
    Ok(counts)
}

/// -----------------------------------------------------------
/// 外显子计数
/// -----------------------------------------------------------
/// 参数:
///   gtf_path       - GTF注释文件路径
///   read_positions - read位置列表
/// 返回: HashMap<String, i64> - exon_id -> read计数
/// 算法:
///   1. 解析GTF，提取feature=="exon"的行
///   2. 获取exon_id属性
///   3. 对于每个read，查找它落在哪个外显子内
/// 用途: 外显子水平的定量分析
#[pyfunction]
pub fn exon_count(
    gtf_path: &str,
    read_positions: Vec<(String, i64)>,
) -> PyResult<HashMap<String, i64>> {
    let file = File::open(gtf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    // 提取所有外显子区间
    let mut exons: Vec<(String, i64, i64, String)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == "exon" {
                let exon_id = rec
                    .attributes
                    .get("exon_id")
                    .cloned()
                    .unwrap_or_default();
                exons.push((rec.chrom, rec.start, rec.end, exon_id));
            }
        }
    }

    // 计数
    let mut counts: HashMap<String, i64> = HashMap::new();
    for (chrom, pos) in &read_positions {
        for (ec, es, ee, eid) in &exons {
            if chrom == ec && pos >= es && pos < ee {
                *counts.entry(eid.clone()).or_insert(0) += 1;
            }
        }
    }
    Ok(counts)
}

/// -----------------------------------------------------------
/// 剪接位点统计
/// -----------------------------------------------------------
/// 参数:
///   cigar - CIGAR字符串
///   pos   - 比对起始位置
/// 返回: Vec<(String, i64, i64)> - 剪接位点列表
/// 算法:
///   CIGAR中的N操作表示跳过参考序列（即内含子）
///   从CIGAR中提取所有N操作的位置和长度
///   格式：(junction描述, 起始位置, 结束位置)
/// 用途: 检测剪接事件、识别新剪接位点
#[pyfunction]
pub fn splice_junction_stats(cigar: &str, pos: i64) -> Vec<(String, i64, i64)> {
    let ops = super::alignment::parse_cigar(cigar);
    let mut junctions = Vec::new();
    let mut current_pos = pos;

    for (op, len) in &ops {
        match op {
            'M' | '=' | 'X' => current_pos += *len as i64,  // 匹配/错配：移动位置
            'N' => {
                // 内含子（剪接位点）
                let junction_start = current_pos;
                let junction_end = current_pos + *len as i64;
                junctions.push((
                    format!("{}-{}", junction_start, junction_end),
                    junction_start,
                    junction_end,
                ));
                current_pos += *len as i64;
            }
            'D' => current_pos += *len as i64,  // 删除：移动位置
            _ => {}
        }
    }
    junctions
}

/// -----------------------------------------------------------
/// 内含子保留检测
/// -----------------------------------------------------------
/// 参数:
///   cigar        - CIGAR字符串
///   pos          - 比对起始位置
///   intron_start - 内含子起始位置
///   intron_end   - 内含子结束位置
/// 返回: bool - 是否发生内含子保留
/// 算法:
///   检查read的匹配区域是否完全覆盖内含子
///   如果M操作覆盖了[intron_start, intron_end]区间，
///   则认为发生了内含子保留
/// 用途: 研究选择性剪接、内含子保留事件
#[pyfunction]
pub fn detect_intron_retention(
    cigar: &str,
    pos: i64,
    intron_start: i64,
    intron_end: i64,
) -> bool {
    let ops = super::alignment::parse_cigar(cigar);
    let mut current_pos = pos;

    for (op, len) in &ops {
        match op {
            'M' | '=' | 'X' => {
                let read_end = current_pos + *len as i64;
                // 检查read是否覆盖整个内含子
                if current_pos <= intron_start && read_end >= intron_end {
                    return true;  // 内含子保留
                }
                current_pos += *len as i64;
            }
            'N' => current_pos += *len as i64,  // 跳过内含子
            'D' => current_pos += *len as i64,
            _ => {}
        }
    }
    false
}

/// -----------------------------------------------------------
/// UMI去重（Collapsing）
/// -----------------------------------------------------------
/// 参数:
///   umis          - UMI列表
///   max_distance  - 最大编辑距离（视为同一UMI的阈值）
/// 返回: Vec<(String, Vec<usize>)> - (代表UMI, 索引列表)
/// 算法:
///   1. 遍历所有UMI
///   2. 对于每个UMI，检查是否可以合并到已有组
///   3. 合并条件：与组代表UMI的编辑距离≤max_distance
///   4. 无法合并则创建新组
/// 用途:
///   PCR扩增会产生重复reads
///   UMI可以标记原始分子
///   编辑距离1可以容忍测序错误
#[pyfunction]
pub fn umi_collapse(
    umis: Vec<String>,
    max_distance: usize,
) -> Vec<(String, Vec<usize>)> {
    let mut groups: Vec<(String, Vec<usize>)> = Vec::new();

    for (i, umi) in umis.iter().enumerate() {
        let mut merged = false;
        // 尝试合并到已有组
        for group in groups.iter_mut() {
            let dist = super::alignment::levenshtein_distance(&group.0, umi);
            if dist <= max_distance {
                group.1.push(i);  // 添加到组
                merged = true;
                break;
            }
        }
        if !merged {
            // 创建新组
            groups.push((umi.clone(), vec![i]));
        }
    }
    groups
}

/// -----------------------------------------------------------
/// 基因类型（biotype）统计
/// -----------------------------------------------------------
/// 参数: gtf_path - GTF文件路径
/// 返回: HashMap<String, usize> - 基因类型 -> 数量
/// 算法:
///   从GTF中提取gene类型的行，获取gene_biotype属性
///   统计每种biotype的数量
/// 常见biotype:
///   protein_coding: 蛋白编码基因
///   lncRNA: 长非编码RNA
///   miRNA: 微小RNA
///   rRNA: 核糖体RNA
/// 用途: 了解注释文件中的基因组成
#[pyfunction]
pub fn gene_biotype_stats(gtf_path: &str) -> PyResult<HashMap<String, usize>> {
    let file = File::open(gtf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let mut biotypes: HashMap<String, usize> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == "gene" {
                let biotype = rec
                    .attributes
                    .get("gene_biotype")
                    .cloned()
                    .unwrap_or_else(|| "unknown".to_string());
                *biotypes.entry(biotype).or_insert(0) += 1;
            }
        }
    }
    Ok(biotypes)
}

/// -----------------------------------------------------------
/// 转录本长度统计
/// -----------------------------------------------------------
/// 参数: gtf_path - GTF文件路径
/// 返回: HashMap<String, i64> - transcript_id -> 总外显子长度
/// 算法:
///   1. 解析GTF中feature=="exon"的行
///   2. 获取transcript_id属性
///   3. 累加每个转录本的所有外显子长度
/// 用途:
///   转录本长度用于：
///   - RPKM/FPKM标准化
///   - 长度偏倚校正
///   - 转录本完整性评估
#[pyfunction]
pub fn transcript_length_stats(gtf_path: &str) -> PyResult<HashMap<String, i64>> {
    let file = File::open(gtf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let mut exon_lengths: HashMap<String, i64> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == "exon" {
                let tx_id = rec
                    .attributes
                    .get("transcript_id")
                    .cloned()
                    .unwrap_or_default();
                // 累加外显子长度
                *exon_lengths.entry(tx_id).or_insert(0) += rec.end - rec.start;
            }
        }
    }
    Ok(exon_lengths)
}
