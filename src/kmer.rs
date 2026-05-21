//! ============================================================
//! 模块3: k-mer分析与序列草图 (kmer)
//! ============================================================
//! 本模块提供k-mer相关的序列分析功能。
//! 包括：k-mer计数、canonical k-mer、minimizer提取、
//! MinHash草图、Jaccard相似度、Mash距离、k-mer频谱等。
//!
//! k-mer是基因组学中的基本概念：
//! - k-mer: 长度为k的连续子序列
//! - 对于长度为L的序列，有 L-k+1 个k-mer
//! - k-mer频率可用于物种鉴定、序列比较、组装等
//!
//! Canonical k-mer:
//!   DNA是双链的，一个k-mer和它的反向互补应该被视为同一个
//!   canonical k-mer取两者中字典序较小的那个
//!
//! MinHash/Minimizer:
//!   用于快速估算序列相似度的数据结构
//!   通过哈希函数将k-mer映射到整数，取最小值作为签名
//!
//! 设计原则：
//! - 统一转大写处理，保证大小写一致性
//! - 参数验证：k必须>0且不超过序列长度
//! - 返回Python兼容类型（HashMap/Vec）
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use std::collections::HashMap;  // Rust 标准哈希表
use std::hash::{Hash, Hasher};  // 哈希相关trait
use std::collections::hash_map::DefaultHasher;  // 默认哈希器

/// -----------------------------------------------------------
/// 计算k-mer计数
/// -----------------------------------------------------------
/// 参数:
///   seq - DNA序列字符串
///   k   - k-mer长度
/// 返回: HashMap<String, usize> - 每个k-mer的出现次数
/// 算法:
///   使用滑动窗口，每次提取长度为k的子序列，统计出现次数
///   时间复杂度: O(L)，其中L为序列长度
/// 用途: 序列组成分析、物种鉴定、污染检测
#[pyfunction]
pub fn count_kmers(seq: &str, k: usize) -> PyResult<HashMap<String, usize>> {
    // 参数验证
    if k == 0 || k > seq.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "k must be > 0 and <= seq length",
        ));
    }
    let seq = seq.to_uppercase();  // 统一转大写
    let mut kmers: HashMap<String, usize> = HashMap::new();
    let bytes = seq.as_bytes();
    // 滑动窗口提取所有k-mer
    for i in 0..=seq.len() - k {
        let kmer = String::from_utf8_lossy(&bytes[i..i + k]).to_string();
        *kmers.entry(kmer).or_insert(0) += 1;  // 计数累加
    }
    Ok(kmers)
}

/// -----------------------------------------------------------
/// 计算canonical k-mer
/// -----------------------------------------------------------
/// 参数:
///   seq - DNA序列字符串
///   k   - k-mer长度
/// 返回: HashMap<String, usize> - canonical k-mer计数
/// 算法:
///   对每个k-mer，计算其反向互补序列
///   取k-mer和反向互补中字典序较小的作为canonical k-mer
///   这样正链和负链的同一个位置会映射到同一个canonical k-mer
/// 生物学意义:
///   DNA是双链的，测序可能来自任一链
///   使用canonical k-mer可以避免链方向的影响
/// 用途: 基因组比较、序列相似度计算
#[pyfunction]
pub fn canonical_kmer(seq: &str, k: usize) -> PyResult<HashMap<String, usize>> {
    if k == 0 || k > seq.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "k must be > 0 and <= seq length",
        ));
    }
    let seq = seq.to_uppercase();
    let mut kmers: HashMap<String, usize> = HashMap::new();
    let bytes = seq.as_bytes();
    for i in 0..=seq.len() - k {
        let kmer = String::from_utf8_lossy(&bytes[i..i + k]).to_string();
        let rc = super::sequence::reverse_complement(&kmer);  // 计算反向互补
        // 取字典序较小的作为canonical k-mer
        let canonical = if kmer <= rc { kmer } else { rc };
        *kmers.entry(canonical).or_insert(0) += 1;
    }
    Ok(kmers)
}

/// -----------------------------------------------------------
/// Minimizer提取
/// -----------------------------------------------------------
/// 参数:
///   seq - DNA序列字符串
///   k   - k-mer长度
///   w   - 窗口大小（w个连续k-mer中取最小的）
/// 返回: Vec<(usize, String)> - (位置, minimizer k-mer)
/// 算法:
///   1. 在每个长度为w的窗口内，找出字典序最小的k-mer
///   2. 如果当前窗口的minimizer与上一个不同，则记录
///   3. 这样可以大幅减少需要存储的k-mer数量
/// 用途: 序列索引、快速比对（minimap2的核心思想）
#[pyfunction]
pub fn minimizers(seq: &str, k: usize, w: usize) -> PyResult<Vec<(usize, String)>> {
    if k == 0 || w == 0 || k > seq.len() || w > seq.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Invalid k or w",
        ));
    }
    let seq = seq.to_uppercase();
    let mut result = Vec::new();
    let mut last_min = String::new();  // 记录上一个minimizer

    // 滑动窗口
    for i in 0..=seq.len().saturating_sub(w + k - 1) {
        let mut min_kmer = String::new();  // 当前窗口的最小k-mer
        let mut min_pos = 0;               // 最小k-mer的位置
        // 在窗口内找最小k-mer
        for j in 0..w {
            let pos = i + j;
            if pos + k > seq.len() {
                break;
            }
            let kmer = &seq[pos..pos + k];
            if min_kmer.is_empty() || kmer < min_kmer.as_str() {
                min_kmer = kmer.to_string();
                min_pos = pos;
            }
        }
        // 去重：只记录变化的minimizer
        if min_kmer != last_min {
            result.push((min_pos, min_kmer.clone()));
            last_min = min_kmer;
        }
    }
    Ok(result)
}

/// -----------------------------------------------------------
/// MinHash sketch生成
/// -----------------------------------------------------------
/// 参数:
///   seq        - DNA序列字符串
///   k          - k-mer长度
///   num_hashes - 保留的哈希值数量（sketch大小）
/// 返回: Vec<u64> - 哈希值签名（排序去重后的前num_hashes个）
/// 算法:
///   1. 计算所有k-mer的哈希值
///   2. 对哈希值排序去重
///   3. 保留最小的num_hashes个哈希值
/// 用途: 快速估算两个序列的相似度（MinHash算法）
#[pyfunction]
pub fn minhash_sketch(seq: &str, k: usize, num_hashes: usize) -> PyResult<Vec<u64>> {
    if k == 0 || k > seq.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "k must be > 0 and <= seq length",
        ));
    }
    let seq = seq.to_uppercase();
    let bytes = seq.as_bytes();
    let mut hashes: Vec<u64> = Vec::new();

    // 计算所有k-mer的哈希值
    for i in 0..=seq.len() - k {
        let kmer = &bytes[i..i + k];
        let hash = hash_kmer(kmer);
        hashes.push(hash);
    }
    // 排序去重，保留最小的num_hashes个
    hashes.sort_unstable();
    hashes.dedup();
    hashes.truncate(num_hashes);
    Ok(hashes)
}

/// -----------------------------------------------------------
/// k-mer哈希函数（内部辅助函数）
/// -----------------------------------------------------------
/// 参数: kmer - k-mer字节序列
/// 返回: u64 - 哈希值
/// 算法: 使用Rust标准库的DefaultHasher进行哈希
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    kmer.hash(&mut hasher);
    hasher.finish()
}

/// -----------------------------------------------------------
/// Jaccard相似度计算
/// -----------------------------------------------------------
/// 参数:
///   seq1 - 序列1
///   seq2 - 序列2
///   k    - k-mer长度
/// 返回: f64 - Jaccard相似度 (0.0 ~ 1.0)
/// 算法:
///   J(A,B) = |A ∩ B| / |A ∪ B|
///   其中A和B分别是两个序列的k-mer集合
/// 用途: 序列相似度比较、聚类分析
#[pyfunction]
pub fn jaccard_similarity(seq1: &str, seq2: &str, k: usize) -> PyResult<f64> {
    let kmers1 = count_kmers(seq1, k)?;
    let kmers2 = count_kmers(seq2, k)?;

    // 提取k-mer集合（去重）
    let set1: std::collections::HashSet<&String> = kmers1.keys().collect();
    let set2: std::collections::HashSet<&String> = kmers2.keys().collect();

    // 计算交集和并集大小
    let intersection = set1.intersection(&set2).count();
    let union = set1.union(&set2).count();

    if union == 0 {
        return Ok(0.0);
    }
    Ok(intersection as f64 / union as f64)
}

/// -----------------------------------------------------------
/// Containment index（包含指数）
/// -----------------------------------------------------------
/// 参数:
///   seq1 - 序列1（被查询序列）
///   seq2 - 序列2（参考序列）
///   k    - k-mer长度
/// 返回: f64 - Containment index (0.0 ~ 1.0)
/// 算法:
///   C(A,B) = |A ∩ B| / |A|
///   表示A中有多少比例的k-mer出现在B中
///   与Jaccard不同，Containment对序列长度不敏感
/// 用途: 检测序列是否包含在另一个序列中（如污染检测）
#[pyfunction]
pub fn containment_index(seq1: &str, seq2: &str, k: usize) -> PyResult<f64> {
    let kmers1 = count_kmers(seq1, k)?;
    let kmers2 = count_kmers(seq2, k)?;

    let set1: std::collections::HashSet<&String> = kmers1.keys().collect();
    let set2: std::collections::HashSet<&String> = kmers2.keys().collect();

    let intersection = set1.intersection(&set2).count();

    if set1.is_empty() {
        return Ok(0.0);
    }
    Ok(intersection as f64 / set1.len() as f64)
}

/// -----------------------------------------------------------
/// Mash-like距离估算
/// -----------------------------------------------------------
/// 参数:
///   seq1        - 序列1
///   seq2        - 序列2
///   k           - k-mer长度
///   sketch_size - MinHash草图大小
/// 返回: f64 - Mash距离 (0.0 ~ 1.0)
/// 算法:
///   1. 使用MinHash为两个序列生成草图
///   2. 计算两个草图的Jaccard相似度
///   3. Mash距离 = -1/k * ln(J)
///   4. Jaccard为0时距离为1.0
/// 用途: 快速估算基因组距离，用于大规模比较
#[pyfunction]
pub fn mash_distance(seq1: &str, seq2: &str, k: usize, sketch_size: usize) -> PyResult<f64> {
    // 生成两个序列的MinHash草图
    let sketch1 = minhash_sketch(seq1, k, sketch_size)?;
    let sketch2 = minhash_sketch(seq2, k, sketch_size)?;

    // 转为HashSet计算Jaccard
    let set1: std::collections::HashSet<u64> = sketch1.into_iter().collect();
    let set2: std::collections::HashSet<u64> = sketch2.into_iter().collect();

    let intersection = set1.intersection(&set2).count();
    let union = set1.union(&set2).count();

    if union == 0 {
        return Ok(1.0);  // 无共同k-mer，距离最大
    }
    let jaccard = intersection as f64 / union as f64;
    if jaccard <= 0.0 {
        return Ok(1.0);
    }
    // Mash距离公式: D = -1/k * ln(J)
    let distance = -1.0 / k as f64 * jaccard.ln();
    Ok(distance.min(1.0))  // 限制在[0,1]范围内
}

/// -----------------------------------------------------------
/// k-mer频谱分析
/// -----------------------------------------------------------
/// 参数:
///   seq - DNA序列
///   k   - k-mer长度
/// 返回: HashMap<usize, usize> - 频次 -> 出现该频次的k-mer种类数
/// 算法:
///   1. 统计所有k-mer的出现次数
///   2. 统计每个出现次数对应的k-mer种类数
///   3. 例如：{1: 100, 2: 50} 表示有100种k-mer出现1次，50种出现2次
/// 用途: 基因组组装质量评估、重复序列检测
#[pyfunction]
pub fn kmer_spectrum(seq: &str, k: usize) -> PyResult<HashMap<usize, usize>> {
    let kmers = count_kmers(seq, k)?;
    let mut spectrum: HashMap<usize, usize> = HashMap::new();
    // 统计频次分布
    for (_, &count) in &kmers {
        *spectrum.entry(count).or_insert(0) += 1;
    }
    Ok(spectrum)
}

/// -----------------------------------------------------------
/// Syncmer提取
/// -----------------------------------------------------------
/// 参数:
///   seq - DNA序列
///   k   - k-mer长度
///   s   - s-mer长度（s必须<k）
/// 返回: Vec<(usize, String)> - (位置, syncmer k-mer)
/// 算法:
///   1. 对每个k-mer，计算其所有s-mer子串
///   2. 找出字典序最小的s-mer
///   3. 如果最小s-mer出现在k-mer的首尾位置，则该k-mer是syncmer
///   4. Syncmer比minimizer有更好的均匀性
/// 用途: 序列索引、去重、比对（比minimizer更均匀的采样）
#[pyfunction]
pub fn syncmers(seq: &str, k: usize, s: usize) -> PyResult<Vec<(usize, String)>> {
    if s >= k || k > seq.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Invalid k or s (s must be < k)",
        ));
    }
    let seq = seq.to_uppercase();
    let mut result = Vec::new();

    // 遍历所有k-mer
    for i in 0..=seq.len() - k {
        let kmer = &seq[i..i + k];
        // 找出k-mer中字典序最小的s-mer
        let mut min_smer = String::new();
        for j in 0..=k - s {
            let smer = &kmer[j..j + s];
            if min_smer.is_empty() || smer < min_smer.as_str() {
                min_smer = smer.to_string();
            }
        }
        // 检查首尾s-mer是否为最小值
        let first_smer = &kmer[0..s];
        let last_smer = &kmer[k - s..k];
        if first_smer == min_smer || last_smer == min_smer {
            result.push((i, kmer.to_string()));
        }
    }
    Ok(result)
}
