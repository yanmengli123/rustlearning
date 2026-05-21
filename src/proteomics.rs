//! ============================================================
//! 模块12: 蛋白质组学 (proteomics)
//! //! 本模块提供蛋白质组学数据分析功能。
//! 包括：胰蛋白酶酶切、肽段质量计算、missed cleavage枚举、
//! 肽段唯一性判断、decoy数据库生成、修饰位点枚举、
//! 肽段m/z计算、FASTA解析、分子量计算、等电点估计等。
//!
//! 蛋白质组学研究细胞/组织中的蛋白质组成和修饰：
//! - Bottom-up蛋白质组学：先酶切为肽段，再质谱分析
//! - Top-down蛋白质组学：直接分析完整蛋白质
//!
//! 关键概念：
//! - Trypsin：最常用的蛋白酶，切割K/R后（P前不切）
//! - Missed cleavage：酶切不完全产生的长肽段
//! - Decoy database：反向序列数据库，用于FDR控制
//! - m/z：质荷比，质谱检测的核心参数
//! - pI (等电点)：蛋白质净电荷为0时的pH值
//!
//! 设计原则：
//! - 使用标准氨基酸质量和分子量
//! - 支持多种蛋白质组学计算
//! - 返回Python兼容类型
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use std::collections::HashMap;  // Rust 标准哈希表

/// -----------------------------------------------------------
/// Trypsin酶切模拟
/// -----------------------------------------------------------
/// 参数:
///   protein               - 蛋白质序列
///   max_missed_cleavages  - 最大允许missed cleavage数
/// 返回: Vec<String> - 酶切产生的肽段列表
/// 算法:
///   Trypsin切割规则：
///   1. 在K或R后切割
///   2. 如果下一个氨基酸是P，则不切割
///   例如：K/R↓X（X≠P）
///
///   生成所有可能的missed cleavage肽段：
///   - 0 missed cleavage: 连续切割
///   - 1 missed cleavage: 跳过一个切割位点
///   - n missed cleavage: 跳过n个切割位点
/// 用途: 质谱数据库搜索前的肽段生成
#[pyfunction]
pub fn trypsin_digest(protein: &str, max_missed_cleavages: usize) -> Vec<String> {
    let mut cleavage_sites: Vec<usize> = vec![0];  // 起始位置
    let chars: Vec<char> = protein.chars().collect();

    // 找到所有切割位点
    for i in 0..chars.len().saturating_sub(1) {
        // K或R后切割，但P前不切
        if (chars[i] == 'K' || chars[i] == 'R') && chars[i + 1] != 'P' {
            cleavage_sites.push(i + 1);
        }
    }
    cleavage_sites.push(chars.len());  // 结束位置

    // 生成所有可能的肽段（包括missed cleavage）
    let mut peptides = Vec::new();
    for mc in 0..=max_missed_cleavages {
        for i in 0..cleavage_sites.len().saturating_sub(1 + mc) {
            let start = cleavage_sites[i];
            let end = cleavage_sites[i + 1 + mc];
            if end > start && end <= chars.len() {
                let peptide: String = chars[start..end].iter().collect();
                if !peptide.is_empty() {
                    peptides.push(peptide);
                }
            }
        }
    }
    peptides
}

/// -----------------------------------------------------------
/// 肽段质量计算（单同位素质量）
/// -----------------------------------------------------------
/// 参数: peptide - 肽段序列
/// 返回: f64 - 单同位素质量（Da）
/// 算法:
///   质量 = H2O + Σ(氨基酸残基质量)
///   H2O质量 = 18.01056 Da
///   使用单同位素质量（最轻同位素）
/// 用途: 质谱峰匹配、肽段鉴定
#[pyfunction]
pub fn peptide_mass(peptide: &str) -> f64 {
    // 氨基酸残基单同位素质量表（Da）
    let aa_masses: HashMap<char, f64> = HashMap::from([
        ('G', 57.02146), ('A', 71.03711), ('V', 99.06841),
        ('L', 113.08406), ('I', 113.08406), ('P', 97.05276),
        ('F', 147.06841), ('W', 186.07931), ('M', 131.04049),
        ('S', 87.03203), ('T', 101.04768), ('C', 103.00919),
        ('Y', 163.06333), ('H', 137.05891), ('D', 115.02694),
        ('E', 129.04259), ('N', 114.04293), ('Q', 128.05858),
        ('K', 128.09496), ('R', 156.10111),
    ]);

    let mut mass = 18.01056; // H2O
    for c in peptide.chars() {
        mass += aa_masses.get(&c).unwrap_or(&0.0);
    }
    mass
}

/// -----------------------------------------------------------
/// Missed cleavage枚举
/// -----------------------------------------------------------
/// 参数: peptide - 肽段序列
/// 返回: Vec<String> - 所有可能的子肽段
/// 算法: 枚举所有可能的连续子串
/// 用途: 生成理论碎片离子谱
#[pyfunction]
pub fn enumerate_missed_cleavages(peptide: &str) -> Vec<String> {
    let mut result = Vec::new();
    let chars: Vec<char> = peptide.chars().collect();
    // 枚举所有长度的子串
    for len in 1..=chars.len() {
        for start in 0..=chars.len() - len {
            let sub: String = chars[start..start + len].iter().collect();
            result.push(sub);
        }
    }
    result
}

/// -----------------------------------------------------------
/// Peptide uniqueness判断
/// -----------------------------------------------------------
/// 参数:
///   peptide          - 肽段序列
///   protein_database - 蛋白质数据库
/// 返回: Vec<usize> - 包含该肽段的蛋白质索引列表
/// 用途:
///   唯一肽段可以用于蛋白质定量
///   非唯一肽段可能来自同源蛋白质
#[pyfunction]
pub fn is_unique_peptide(peptide: &str, protein_database: Vec<String>) -> Vec<usize> {
    let mut protein_indices = Vec::new();
    for (i, protein) in protein_database.iter().enumerate() {
        if protein.contains(peptide) {
            protein_indices.push(i);
        }
    }
    protein_indices
}

/// -----------------------------------------------------------
/// Decoy database生成（反转序列）
/// -----------------------------------------------------------
/// 参数: sequences - 原始蛋白质序列列表
/// 返回: Vec<String> - 反转后的序列列表
/// 算法:
///   将蛋白质序列反转（如 "ACDE" -> "EDCA"）
///   保持氨基酸组成不变，但序列顺序随机化
/// 用途:
///   Decoy数据库用于控制假发现率（FDR）
///   在目标-诱饵策略中，同时搜索target和decoy数据库
#[pyfunction]
pub fn generate_decoy(sequences: Vec<String>) -> Vec<String> {
    sequences
        .iter()
        .map(|seq| seq.chars().rev().collect())
        .collect()
}

/// -----------------------------------------------------------
/// Modification site枚举
/// -----------------------------------------------------------
/// 参数:
///   peptide       - 肽段序列
///   mod_residues  - 可修饰的氨基酸列表（如["M", "C"]）
/// 返回: Vec<usize> - 可修饰位点的位置索引
/// 用途: 翻译后修饰（PTM）位点预测
#[pyfunction]
pub fn enumerate_modification_sites(
    peptide: &str,
    mod_residues: Vec<String>,
) -> Vec<usize> {
    let mut sites = Vec::new();
    for (i, c) in peptide.chars().enumerate() {
        // 检查是否为可修饰氨基酸
        if mod_residues.iter().any(|s| s.chars().next() == Some(c)) {
            sites.push(i);
        }
    }
    sites
}

/// -----------------------------------------------------------
/// 肽段m/z计算
/// -----------------------------------------------------------
/// 参数:
///   peptide - 肽段序列
///   charge  - 电荷状态
/// 返回: f64 - 质荷比 (m/z)
/// 算法:
///   m/z = (M + z * H) / z
///   M = 肽段质量
///   z = 电荷数
///   H = 质子质量 = 1.007276 Da
/// 用途: 质谱峰预测、母离子质量计算
#[pyfunction]
pub fn peptide_mz(peptide: &str, charge: i32) -> f64 {
    let mass = peptide_mass(peptide);
    let proton = 1.007276;  // 质子质量
    (mass + charge as f64 * proton) / charge as f64
}

/// -----------------------------------------------------------
/// FASTA蛋白库解析
/// -----------------------------------------------------------
/// 参数: path - FASTA文件路径
/// 返回: Vec<(String, String)> - (header, 序列)列表
/// 算法:
///   1. 逐行读取FASTA文件
///   2. 以>开头的行是header
///   3. 后续行是序列（可跨多行）
/// 用途: 加载蛋白质数据库用于搜索
#[pyfunction]
pub fn parse_fasta(path: &str) -> PyResult<Vec<(String, String)>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut sequences = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();

    for line in content.lines() {
        if line.starts_with('>') {
            // 保存前一条序列
            if !current_header.is_empty() {
                sequences.push((current_header.clone(), current_seq.clone()));
            }
            current_header = line[1..].to_string();  // 去掉>前缀
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());  // 拼接序列行
        }
    }
    // 保存最后一条序列
    if !current_header.is_empty() {
        sequences.push((current_header, current_seq));
    }
    Ok(sequences)
}

/// -----------------------------------------------------------
/// 蛋白质分子量计算
/// -----------------------------------------------------------
/// 参数: protein - 蛋白质序列
/// 返回: f64 - 分子量（Da）
/// 算法:
///   分子量 = H2O + Σ(氨基酸平均分子量)
///   使用平均分子量（而非单同位素质量）
/// 用途: 蛋白质鉴定、SDS-PAGE预测
#[pyfunction]
pub fn protein_molecular_weight(protein: &str) -> f64 {
    // 氨基酸平均分子量表（Da）
    let aa_weights: HashMap<char, f64> = HashMap::from([
        ('G', 57.052), ('A', 71.079), ('V', 99.133),
        ('L', 113.160), ('I', 113.160), ('P', 97.117),
        ('F', 147.177), ('W', 186.213), ('M', 131.199),
        ('S', 87.078), ('T', 101.105), ('C', 103.144),
        ('Y', 163.176), ('H', 137.142), ('D', 115.089),
        ('E', 129.116), ('N', 114.104), ('Q', 128.131),
        ('K', 128.174), ('R', 156.188),
    ]);

    let mut weight = 18.015; // H2O
    for c in protein.chars() {
        weight += aa_weights.get(&c).unwrap_or(&0.0);
    }
    weight
}

/// -----------------------------------------------------------
/// 肽段等电点（pI）估计
/// -----------------------------------------------------------
/// 参数: peptide - 肽段序列
/// 返回: f64 - 估计的等电点
/// 算法（简化版）：
///   1. 统计带正电荷的氨基酸（K, R, H）
///   2. 统计带负电荷的氨基酸（D, E）
///   3. pI ≈ 7.0 + (正电荷数 - 负电荷数) * 0.5
/// 说明: 这是简化计算，精确计算需要Henderson-Hasselbalch方程
/// 用途: 离子交换色谱条件优化
#[pyfunction]
pub fn peptide_pi(peptide: &str) -> f64 {
    let mut positive = 0.0f64;  // 正电荷氨基酸数
    let mut negative = 0.0f64;  // 负电荷氨基酸数

    for c in peptide.chars() {
        match c {
            'K' | 'R' | 'H' => positive += 1.0,  // 碱性氨基酸
            'D' | 'E' => negative += 1.0,         // 酸性氨基酸
            _ => {}
        }
    }
    // 简化的pI计算
    7.0 + (positive - negative) * 0.5
}
