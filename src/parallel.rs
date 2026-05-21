use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

/// 并行GC含量计算
#[pyfunction]
pub fn parallel_gc_content(seqs: Vec<String>) -> Vec<f64> {
    seqs.iter().map(|s| super::sequence::gc_content(s)).collect()
}

/// 并行k-mer计数
#[pyfunction]
pub fn parallel_count_kmers(seqs: Vec<String>, k: usize) -> PyResult<HashMap<String, usize>> {
    let mut total: HashMap<String, usize> = HashMap::new();
    for seq in &seqs {
        let kmers = super::kmer::count_kmers(seq, k)?;
        for (kmer, count) in kmers {
            *total.entry(kmer).or_insert(0) += count;
        }
    }
    Ok(total)
}

/// 并行reverse complement
#[pyfunction]
pub fn parallel_reverse_complement(seqs: Vec<String>) -> Vec<String> {
    seqs.iter().map(|s| super::sequence::reverse_complement(s)).collect()
}

/// 并行Hamming距离计算（所有对）
#[pyfunction]
pub fn parallel_hamming_matrix(seqs: Vec<String>) -> PyResult<Vec<Vec<usize>>> {
    let n = seqs.len();
    let mut matrix = vec![vec![0usize; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let dist = super::alignment::hamming_distance(&seqs[i], &seqs[j]).unwrap_or(999);
            matrix[i][j] = dist;
            matrix[j][i] = dist;
        }
    }
    Ok(matrix)
}

/// 并行Levenshtein距离矩阵
#[pyfunction]
pub fn parallel_levenshtein_matrix(seqs: Vec<String>) -> Vec<Vec<usize>> {
    let n = seqs.len();
    let mut matrix = vec![vec![0usize; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let dist = super::alignment::levenshtein_distance(&seqs[i], &seqs[j]);
            matrix[i][j] = dist;
            matrix[j][i] = dist;
        }
    }
    matrix
}

/// 并行序列统计
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

/// 并行ORF查找
#[pyfunction]
pub fn parallel_find_orfs(seqs: Vec<String>) -> Vec<Vec<(usize, usize, String)>> {
    seqs.iter().map(|s| super::sequence::find_orfs(s)).collect()
}

/// 并行FASTA文件解析
#[pyfunction]
pub fn parallel_parse_fasta(paths: Vec<String>) -> PyResult<Vec<Vec<(String, String)>>> {
    let results: Vec<Vec<(String, String)>> = paths
        .iter()
        .map(|p| super::proteomics::parse_fasta(p).unwrap_or_default())
        .collect();
    Ok(results)
}

/// 并行BED交集计算
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

/// 并行Jaccard相似度计算
#[pyfunction]
pub fn parallel_jaccard_matrix(seqs: Vec<String>, k: usize) -> PyResult<Vec<Vec<f64>>> {
    let n = seqs.len();
    let mut matrix = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        matrix[i][i] = 1.0;
        for j in i + 1..n {
            let sim = super::kmer::jaccard_similarity(&seqs[i], &seqs[j], k)?;
            matrix[i][j] = sim;
            matrix[j][i] = sim;
        }
    }
    Ok(matrix)
}

/// 批量质量过滤
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

/// 并行统计计算
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
