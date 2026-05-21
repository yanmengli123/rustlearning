//! ============================================================
//! 模块11: 宏基因组学 (metagenomics)
//! //! 本模块提供宏基因组学数据分析功能。
//! 包括：k-mer分类、ANI近似计算、基因组草图、
//! 16S序列预处理、marker gene匹配、污染筛查、
//! OTU/ASV聚类、丰度估计等。
//!
//! 宏基因组学研究环境样本中的微生物群落：
//! - 宏基因组测序（Shotgun metagenomics）：直接测序所有DNA
//! - 16S rRNA测序：针对标记基因的扩增子测序
//! - OTU (Operational Taxonomic Unit)：操作分类单元
//! - ASV (Amplicon Sequence Variant)：扩增子序列变体
//!
//! 关键概念：
//! - ANI (Average Nucleotide Identity)：平均核苷酸相似度
//! - Marker gene: 用于物种鉴定的保守基因（如16S rRNA）
//! - k-mer分类：基于k-mer频率的物种分类
//! - Greedy clustering: 贪心聚类算法
//!
//! 设计原则：
//! - 使用k-mer方法进行快速比较
//! - 支持MinHash草图用于大规模比较
//! - 返回Python兼容类型
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use std::collections::HashMap;  // Rust 标准哈希表

/// -----------------------------------------------------------
/// k-mer分类器（简化版）
/// -----------------------------------------------------------
/// 参数:
///   query_seq        - 查询序列
///   reference_kmers  - 参考k-mer数据库 (kmer -> taxon列表)
///   k                - k-mer长度
/// 返回: HashMap<String, f64> - taxon -> 相对丰度
/// 算法:
///   1. 计算查询序列的k-mer集合
///   2. 与参考数据库比对，统计每个taxon的命中数
///   3. 计算相对丰度 = taxon命中数 / 总命中数
/// 用途: 物种组成分析
#[pyfunction]
pub fn kmer_classify(
    query_seq: &str,
    reference_kmers: HashMap<String, Vec<String>>,
    k: usize,
) -> HashMap<String, f64> {
    let query_kmers = super::kmer::count_kmers(query_seq, k).unwrap_or_default();
    let query_set: std::collections::HashSet<&String> = query_kmers.keys().collect();

    let mut taxon_hits: HashMap<String, usize> = HashMap::new();
    let mut total_hits = 0usize;

    // 统计每个taxon的命中数
    for (kmer, taxa) in &reference_kmers {
        if query_set.contains(kmer) {
            for taxon in taxa {
                *taxon_hits.entry(taxon.clone()).or_insert(0) += 1;
                total_hits += 1;
            }
        }
    }

    // 计算相对丰度
    let mut abundances: HashMap<String, f64> = HashMap::new();
    if total_hits > 0 {
        for (taxon, hits) in &taxon_hits {
            abundances.insert(taxon.clone(), *hits as f64 / total_hits as f64);
        }
    }
    abundances
}

/// -----------------------------------------------------------
/// ANI近似计算（Average Nucleotide Identity）
/// -----------------------------------------------------------
/// 参数:
///   seq1         - 序列1
///   seq2         - 序列2
///   k            - k-mer长度
///   sketch_size  - MinHash草图大小
/// 返回: f64 - ANI近似值 (0.0 ~ 1.0)
/// 算法:
///   使用MinHash草图估算Jaccard相似度
///   ANI ≈ 1 + (1/k) * ln(2*J/(1+J))
///   这是Mash distance的逆运算
/// 用途:
///   ANI是物种鉴定的金标准
///   ANI>95%通常认为是同一物种
#[pyfunction]
pub fn ani_approximate(seq1: &str, seq2: &str, k: usize, sketch_size: usize) -> PyResult<f64> {
    // 生成两个序列的MinHash草图
    let sketch1 = super::kmer::minhash_sketch(seq1, k, sketch_size)?;
    let sketch2 = super::kmer::minhash_sketch(seq2, k, sketch_size)?;

    // 计算Jaccard相似度
    let set1: std::collections::HashSet<u64> = sketch1.into_iter().collect();
    let set2: std::collections::HashSet<u64> = sketch2.into_iter().collect();

    let intersection = set1.intersection(&set2).count();
    let union = set1.union(&set2).count();

    if union == 0 {
        return Ok(0.0);
    }
    let jaccard = intersection as f64 / union as f64;
    // ANI ≈ 1 + (1/k) * ln(2*J/(1+J))
    if jaccard <= 0.0 {
        return Ok(0.0);
    }
    let ani = 1.0 + (1.0 / k as f64) * (2.0 * jaccard / (1.0 + jaccard)).ln();
    Ok(ani.max(0.0).min(1.0))  // 限制在[0,1]范围内
}

/// -----------------------------------------------------------
/// 基因组草图生成
/// -----------------------------------------------------------
/// 参数:
///   seq         - 基因组序列
///   k           - k-mer长度
///   sketch_size - 草图大小
/// 返回: Vec<u64> - MinHash草图（哈希值列表）
/// 用途: 用于大规模基因组比较（如Mash工具）
#[pyfunction]
pub fn genome_sketch(
    seq: &str,
    k: usize,
    sketch_size: usize,
) -> PyResult<Vec<u64>> {
    super::kmer::minhash_sketch(seq, k, sketch_size)
}

/// -----------------------------------------------------------
/// 16S序列预处理（去引物、清洗）
/// -----------------------------------------------------------
/// 参数:
///   seq             - 原始16S序列
///   forward_primer  - 正向引物序列
///   reverse_primer  - 反向引物序列（5'->3'方向）
///   max_mismatch    - 最大允许错配数
/// 返回: Option<String> - 处理后的序列（引物之间）
/// 算法:
///   1. 在序列中查找正向引物
///   2. 查找反向引物的反向互补序列
///   3. 提取两个引物之间的序列
///   4. 如果未找到引物，返回None
/// 用途: 16S扩增子测序数据预处理
#[pyfunction]
pub fn preprocess_16s(
    seq: &str,
    forward_primer: &str,
    reverse_primer: &str,
    max_mismatch: usize,
) -> Option<String> {
    // 查找正向引物
    let fwd_hits = super::alignment::primer_match(seq, forward_primer, max_mismatch);
    // 查找反向引物的反向互补
    let rev_rc = super::sequence::reverse_complement(reverse_primer);
    let rev_hits = super::alignment::primer_match(seq, &rev_rc, max_mismatch);

    if let (Some((fwd_pos, _)), Some((rev_pos, _))) = (fwd_hits.first(), rev_hits.first()) {
        let start = fwd_pos + forward_primer.len();  // 正向引物之后
        let end = *rev_pos;                          // 反向引物之前
        if start < end && end <= seq.len() {
            return Some(seq[start..end].to_string());
        }
    }
    None
}

/// -----------------------------------------------------------
/// Marker gene匹配
/// -----------------------------------------------------------
/// 参数:
///   query    - 查询序列
///   markers  - marker基因列表 (名称, 序列)
///   k        - k-mer长度
///   threshold - 相似度阈值
/// 返回: Vec<(String, f64)> - 匹配的marker及相似度
/// 算法:
///   使用Jaccard相似度比较查询序列与每个marker
///   返回相似度≥threshold的所有匹配
/// 用途: 基于marker基因的物种鉴定
#[pyfunction]
pub fn marker_gene_match(
    query: &str,
    markers: Vec<(String, String)>,
    k: usize,
    threshold: f64,
) -> Vec<(String, f64)> {
    let mut matches = Vec::new();
    for (name, marker_seq) in &markers {
        let sim = super::kmer::jaccard_similarity(query, marker_seq, k).unwrap_or(0.0);
        if sim >= threshold {
            matches.push((name.clone(), sim));
        }
    }
    matches
}

/// -----------------------------------------------------------
/// 污染筛查
/// -----------------------------------------------------------
/// 参数:
///   seq         - 查询序列
///   host_kmers  - 宿主基因组k-mer哈希值列表
///   k           - k-mer长度
///   sketch_size - 草图大小
/// 返回: f64 - 宿主污染比例 (0.0 ~ 1.0)
/// 算法:
///   1. 生成查询序列的MinHash草图
///   2. 计算与宿主k-mer的重叠比例
///   3. 高重叠比例表示可能存在宿主污染
/// 用途: 检测测序数据中的宿主DNA污染
#[pyfunction]
pub fn contamination_screen(
    seq: &str,
    host_kmers: Vec<u64>,
    k: usize,
    sketch_size: usize,
) -> PyResult<f64> {
    let read_sketch = super::kmer::minhash_sketch(seq, k, sketch_size)?;
    let host_set: std::collections::HashSet<u64> = host_kmers.into_iter().collect();
    let read_set: std::collections::HashSet<u64> = read_sketch.into_iter().collect();

    let overlap = read_set.intersection(&host_set).count();
    if read_set.is_empty() {
        return Ok(0.0);
    }
    Ok(overlap as f64 / read_set.len() as f64)
}

/// -----------------------------------------------------------
/// OTU/ASV聚类辅助（Greedy Clustering）
/// -----------------------------------------------------------
/// 参数:
///   sequences         - 序列列表
///   identity_threshold - 相似度阈值（如0.97表示97%相似度）
///   k                 - k-mer长度
/// 返回: Vec<(String, Vec<usize>)> - (代表序列, 索引列表)
/// 算法:
///   1. 遍历每个序列
///   2. 计算与已有聚类代表序列的Jaccard相似度
///   3. 如果相似度≥阈值，加入该聚类
///   4. 否则创建新聚类
/// 用途: 16S序列聚类为OTU
#[pyfunction]
pub fn greedy_cluster(
    sequences: Vec<String>,
    identity_threshold: f64,
    k: usize,
) -> Vec<(String, Vec<usize>)> {
    let mut clusters: Vec<(String, Vec<usize>)> = Vec::new();

    for (i, seq) in sequences.iter().enumerate() {
        let mut assigned = false;
        // 尝试分配到已有聚类
        for cluster in clusters.iter_mut() {
            let sim = super::kmer::jaccard_similarity(&cluster.0, seq, k).unwrap_or(0.0);
            if sim >= identity_threshold {
                cluster.1.push(i);
                assigned = true;
                break;
            }
        }
        if !assigned {
            // 创建新聚类
            clusters.push((seq.clone(), vec![i]));
        }
    }
    clusters
}

/// -----------------------------------------------------------
/// 丰度估计辅助
/// -----------------------------------------------------------
/// 参数:
///   query_kmers     - 查询序列k-mer计数
///   reference_kmers - 参考基因组k-mer计数
/// 返回: f64 - 重叠比例（0.0 ~ 1.0）
/// 算法:
///   计算查询序列与参考基因组共有k-mer的比例
///   重叠比例越高，表示该基因组丰度越高
/// 用途: 估计宏基因组中各物种的相对丰度
#[pyfunction]
pub fn estimate_abundance(
    query_kmers: HashMap<String, usize>,
    reference_kmers: HashMap<String, usize>,
) -> f64 {
    let query_set: std::collections::HashSet<&String> = query_kmers.keys().collect();
    let ref_set: std::collections::HashSet<&String> = reference_kmers.keys().collect();
    let common = query_set.intersection(&ref_set).count();

    if ref_set.is_empty() {
        return 0.0;
    }
    common as f64 / ref_set.len() as f64
}
