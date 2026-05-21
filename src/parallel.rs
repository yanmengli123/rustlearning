//! ============================================================
//! 模块15: 批量/并行处理 (parallel)
//! ============================================================
//! 本模块提供批量处理和并行计算的辅助功能。
//! 包括：并行GC含量计算、并行k-mer计数、并行反向互补、
//! 距离矩阵计算、并行统计、批量FASTQ过滤等。
//!
//! 批量处理的意义：
//! - 减少Python-Rust调用开销
//! - 统一处理多个序列/文件
//! - 简化用户代码
//!
//! 并行计算：
//! - 当前实现为顺序执行（未来可改为真正的并行）
//! - 利用Rust的高性能单线程处理
//! - Python端可使用multiprocessing实现真正的并行
//!
//! 设计原则：
//! - 批量接口：一次处理多个序列/文件
//! - 返回Python兼容类型
//! - 错误处理：使用默认值或跳过失败项
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use std::collections::HashMap;  // Rust 标准哈希表

/// -----------------------------------------------------------
/// 并行GC含量计算
/// -----------------------------------------------------------
/// 参数: seqs - 序列列表
/// 返回: Vec<f64> - 每个序列的GC含量
/// 用途: 批量计算多个序列的GC含量
#[pyfunction]
pub fn parallel_gc_content(seqs: Vec<String>) -> Vec<f64> {
    seqs.iter().map(|s| super::sequence::gc_content(s)).collect()
}

/// -----------------------------------------------------------
/// 并行k-mer计数
/// -----------------------------------------------------------
/// 参数:
///   seqs - 序列列表
///   k    - k-mer长度
/// 返回: HashMap<String, usize> - 所有序列的k-mer总计数
/// 算法:
///   1. 对每个序列计算k-mer计数
///   2. 合并所有序列的计数结果
/// 用途: 批量k-mer分析、宏基因组物种组成
#[pyfunction]
pub fn parallel_count_kmers(seqs: Vec<String>, k: usize) -> PyResult<HashMap<String, usize>> {
    let mut total: HashMap<String, usize> = HashMap::new();
    for seq in &seqs {
        let kmers = super::kmer::count_kmers(seq, k)?;
        // 合并计数
        for (kmer, count) in kmers {
            *total.entry(kmer).or_insert(0) += count;
        }
    }
    Ok(total)
}

/// -----------------------------------------------------------
/// 并行反向互补
/// -----------------------------------------------------------
/// 参数: seqs - 序列列表
/// 返回: Vec<String> - 反向互补序列列表
/// 用途: 批量生成反向互补序列
#[pyfunction]
pub fn parallel_reverse_complement(seqs: Vec<String>) -> Vec<String> {
    seqs.iter().map(|s| super::sequence::reverse_complement(s)).collect()
}

/// -----------------------------------------------------------
/// 并行Hamming距离矩阵计算
/// -----------------------------------------------------------
/// 参数: seqs - 序列列表（必须等长）
/// 返回: Vec<Vec<usize>> - 距离矩阵
/// 算法:
///   计算所有序列对之间的Hamming距离
///   结果矩阵是对称的，对角线为0
/// 用途: 序列聚类、系统发育分析
#[pyfunction]
pub fn parallel_hamming_matrix(seqs: Vec<String>) -> PyResult<Vec<Vec<usize>>> {
    let n = seqs.len();
    let mut matrix = vec![vec![0usize; n]; n];
    // 只计算上三角矩阵，对称填充
    for i in 0..n {
        for j in i + 1..n {
            let dist = super::alignment::hamming_distance(&seqs[i], &seqs[j]).unwrap_or(999);
            matrix[i][j] = dist;
            matrix[j][i] = dist;  // 对称
        }
    }
    Ok(matrix)
}

/// -----------------------------------------------------------
/// 并行Levenshtein距离矩阵计算
/// -----------------------------------------------------------
/// 参数: seqs - 序列列表
/// 返回: Vec<Vec<usize>> - 距离矩阵
/// 用途: 不等长序列的距离计算、序列聚类
#[pyfunction]
pub fn parallel_levenshtein_matrix(seqs: Vec<String>) -> Vec<Vec<usize>> {
    let n = seqs.len();
    let mut matrix = vec![vec![0usize; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let dist = super::alignment::levenshtein_distance(&seqs[i], &seqs[j]);
            matrix[i][j] = dist;
            matrix[j][i] = dist;  // 对称
        }
    }
    matrix
}

/// -----------------------------------------------------------
/// 并行序列统计
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   seqs - 序列列表
/// 返回: Python字典，包含统计列表
/// 统计指标（每个序列）:
///   - length: 长度
///   - gc_content: GC含量
///   - n_count: N碱基数
///   - is_valid_dna: 是否合法DNA
/// 用途: 批量序列质控
#[pyfunction]
pub fn parallel_seq_stats(py: Python, seqs: Vec<String>) -> PyResult<PyObject> {
    let results: Vec<PyObject> = seqs
        .iter()
        .map(|seq| {
            Python::with_gil(|py| {
                super::sequence::seq_stats(py, seq).unwrap()
            })
        })
        .collect();

    let dict = PyDict::new_bound(py);
    dict.set_item("n_sequences", seqs.len())?;
    dict.set_item("stats", results)?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// 并行ORF查找
/// -----------------------------------------------------------
/// 参数: seqs - DNA序列列表
/// 返回: Vec<Vec<(usize, usize, String)>> - 每个序列的ORF列表
/// 用途: 批量基因预测
#[pyfunction]
pub fn parallel_find_orfs(seqs: Vec<String>) -> Vec<Vec<(usize, usize, String)>> {
    seqs.iter().map(|s| super::sequence::find_orfs(s)).collect()
}

/// -----------------------------------------------------------
/// 并行FASTA文件解析
/// -----------------------------------------------------------
/// 参数: paths - FASTA文件路径列表
/// 返回: Vec<Vec<(String, String)>> - 每个文件的序列列表
/// 用途: 批量加载多个FASTA文件
#[pyfunction]
pub fn parallel_parse_fasta(paths: Vec<String>) -> PyResult<Vec<Vec<(String, String)>>> {
    let results: Vec<Vec<(String, String)>> = paths
        .iter()
        .map(|p| super::proteomics::parse_fasta(p).unwrap_or_default())
        .collect();
    Ok(results)
}

/// -----------------------------------------------------------
/// 并行BED交集计算
/// -----------------------------------------------------------
/// 参数:
///   path_a  - BED文件A路径
///   paths_b - BED文件B路径列表
/// 返回: Vec<Vec<String>> - 每个B文件与A的交集结果
/// 用途: 批量区间交集分析
#[pyfunction]
pub fn parallel_bed_intersect(
    path_a: &str,
    paths_b: Vec<String>,
) -> PyResult<Vec<Vec<String>>> {
    let results: Vec<Vec<String>> = paths_b
        .iter()
        .map(|p| super::bed::bed_intersect(path_a, p).unwrap_or_default())
        .collect();
    Ok(results)
}

/// -----------------------------------------------------------
/// 并行Jaccard相似度矩阵
/// -----------------------------------------------------------
/// 参数:
///   seqs - 序列列表
///   k    - k-mer长度
/// 返回: Vec<Vec<f64>> - 相似度矩阵
/// 用途: 序列聚类、物种相似度分析
#[pyfunction]
pub fn parallel_jaccard_matrix(seqs: Vec<String>, k: usize) -> PyResult<Vec<Vec<f64>>> {
    let n = seqs.len();
    let mut matrix = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        matrix[i][i] = 1.0;  // 自身相似度为1
        for j in i + 1..n {
            let sim = super::kmer::jaccard_similarity(&seqs[i], &seqs[j], k)?;
            matrix[i][j] = sim;
            matrix[j][i] = sim;  // 对称
        }
    }
    Ok(matrix)
}

/// -----------------------------------------------------------
/// 批量FASTQ质量过滤
/// -----------------------------------------------------------
/// 参数:
///   inputs   - 输入文件路径列表
///   outputs  - 输出文件路径列表
///   min_len  - 最小长度阈值
///   min_qual - 最小质量阈值
/// 返回: Vec<usize> - 每个文件通过过滤的read数
/// 用途: 批量处理多个FASTQ文件
#[pyfunction]
pub fn batch_fastq_filter(
    inputs: Vec<String>,
    outputs: Vec<String>,
    min_len: usize,
    min_qual: f64,
) -> PyResult<Vec<usize>> {
    let mut results = Vec::new();
    for (inp, out) in inputs.iter().zip(outputs.iter()) {
        let kept = super::fastq::fastq_filter(inp, out, min_len, min_qual)?;
        results.push(kept);
    }
    Ok(results)
}

/// -----------------------------------------------------------
/// 并行描述统计
/// -----------------------------------------------------------
/// 参数:
///   py       - Python解释器引用
///   datasets - 数据集列表
/// 返回: Python字典，包含每个数据集的统计
/// 用途: 批量统计分析
#[pyfunction]
pub fn parallel_descriptive_stats(
    py: Python,
    datasets: Vec<Vec<f64>>,
) -> PyResult<PyObject> {
    let results: Vec<PyObject> = datasets
        .iter()
        .map(|data| {
            Python::with_gil(|py| {
                super::statistics::descriptive_stats(py, data.clone()).unwrap()
            })
        })
        .collect();

    let dict = PyDict::new_bound(py);
    dict.set_item("n_datasets", datasets.len())?;
    dict.set_item("stats", results)?;
    Ok(dict.into())
}
