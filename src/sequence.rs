use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

/// 计算GC含量
#[pyfunction]
pub fn gc_content(seq: &str) -> f64 {
    let seq = seq.to_uppercase();
    let total = seq.len() as f64;
    if total == 0.0 {
        return 0.0;
    }
    let gc = seq.chars().filter(|&c| c == 'G' || c == 'C').count() as f64;
    gc / total
}

/// 反向互补序列
#[pyfunction]
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            'N' | 'n' => 'N',
            'U' | 'u' => 'A',
            _ => c,
        })
        .collect()
}

/// 统计各碱基数量
#[pyfunction]
pub fn count_bases(seq: &str) -> HashMap<char, usize> {
    let mut counts = HashMap::new();
    for c in seq.to_uppercase().chars() {
        *counts.entry(c).or_insert(0) += 1;
    }
    counts
}

/// 序列长度
#[pyfunction]
pub fn seq_length(seq: &str) -> usize {
    seq.len()
}

/// 检查是否为合法DNA序列
#[pyfunction]
pub fn is_valid_dna(seq: &str) -> bool {
    !seq.is_empty()
        && seq
            .chars()
            .all(|c| matches!(c, 'A' | 'a' | 'T' | 't' | 'G' | 'g' | 'C' | 'c' | 'N' | 'n'))
}

/// 检查是否为合法RNA序列
#[pyfunction]
pub fn is_valid_rna(seq: &str) -> bool {
    !seq.is_empty()
        && seq
            .chars()
            .all(|c| matches!(c, 'A' | 'a' | 'U' | 'u' | 'G' | 'g' | 'C' | 'c' | 'N' | 'n'))
}

/// DNA转RNA
#[pyfunction]
pub fn transcribe(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'T' => 'U',
            't' => 'u',
            _ => c,
        })
        .collect()
}

/// RNA转DNA
#[pyfunction]
pub fn reverse_transcribe(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'U' => 'T',
            'u' => 't',
            _ => c,
        })
        .collect()
}

/// 计算N碱基数
#[pyfunction]
pub fn count_n(seq: &str) -> usize {
    seq.chars().filter(|&c| c == 'N' || c == 'n').count()
}

/// 碱基大小写标准化（转大写）
#[pyfunction]
pub fn normalize_seq(seq: &str) -> String {
    seq.to_uppercase()
}

/// 序列切片
#[pyfunction]
pub fn seq_slice(seq: &str, start: usize, end: usize) -> PyResult<String> {
    if start > seq.len() || end > seq.len() || start > end {
        return Err(pyo3::exceptions::PyIndexError::new_err(format!(
            "Invalid range: {}..{} for seq of length {}",
            start,
            end,
            seq.len()
        )));
    }
    Ok(seq[start..end].to_string())
}

/// 拼接多条序列
#[pyfunction]
pub fn concat_seqs(seqs: Vec<String>, separator: &str) -> String {
    seqs.join(separator)
}

/// DNA序列压缩比（无损压缩信息量估计）
#[pyfunction]
pub fn compression_ratio(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let bases = count_bases(seq);
    let len = seq.len() as f64;
    let mut entropy = 0.0;
    for (_, &count) in &bases {
        let p = count as f64 / len;
        if p > 0.0 {
            entropy -= p * p.log2();
        }
    }
    entropy / 2.0 // DNA有4种碱基，最大熵为2
}

/// 低复杂度序列检测
#[pyfunction]
pub fn is_low_complexity(seq: &str, threshold: f64) -> bool {
    if seq.is_empty() {
        return true;
    }
    let bases = count_bases(seq);
    let max_count = bases.values().max().unwrap_or(&0);
    (*max_count as f64 / seq.len() as f64) > threshold
}

/// 寻找ORF（开放阅读框）
#[pyfunction]
pub fn find_orfs(seq: &str) -> Vec<(usize, usize, String)> {
    let seq = seq.to_uppercase();
    let bytes = seq.as_bytes();
    let mut orfs = Vec::new();
    let stop_codons: &[&[u8]] = &[b"TAA", b"TAG", b"TGA"];

    for frame in 0..3 {
        let mut i = frame;
        while i + 3 <= bytes.len() {
            if &bytes[i..i + 3] == b"ATG" {
                let start = i;
                let mut j = i + 3;
                while j + 3 <= bytes.len() {
                    if stop_codons.contains(&&bytes[j..j + 3]) {
                        let orf_seq =
                            String::from_utf8_lossy(&bytes[start..j + 3]).to_string();
                        orfs.push((start, j + 3, orf_seq));
                        break;
                    }
                    j += 3;
                }
            }
            i += 3;
        }
    }
    orfs
}

/// 密码子翻译为氨基酸
#[pyfunction]
pub fn translate(seq: &str) -> String {
    let codon_table: HashMap<&str, char> = HashMap::from([
        ("TTT", 'F'), ("TTC", 'F'), ("TTA", 'L'), ("TTG", 'L'),
        ("CTT", 'L'), ("CTC", 'L'), ("CTA", 'L'), ("CTG", 'L'),
        ("ATT", 'I'), ("ATC", 'I'), ("ATA", 'I'), ("ATG", 'M'),
        ("GTT", 'V'), ("GTC", 'V'), ("GTA", 'V'), ("GTG", 'V'),
        ("TCT", 'S'), ("TCC", 'S'), ("TCA", 'S'), ("TCG", 'S'),
        ("CCT", 'P'), ("CCC", 'P'), ("CCA", 'P'), ("CCG", 'P'),
        ("ACT", 'T'), ("ACC", 'T'), ("ACA", 'T'), ("ACG", 'T'),
        ("GCT", 'A'), ("GCC", 'A'), ("GCA", 'A'), ("GCG", 'A'),
        ("TAT", 'Y'), ("TAC", 'Y'), ("TAA", '*'), ("TAG", '*'),
        ("CAT", 'H'), ("CAC", 'H'), ("CAA", 'Q'), ("CAG", 'Q'),
        ("AAT", 'N'), ("AAC", 'N'), ("AAA", 'K'), ("AAG", 'K'),
        ("GAT", 'D'), ("GAC", 'D'), ("GAA", 'E'), ("GAG", 'E'),
        ("TGT", 'C'), ("TGC", 'C'), ("TGA", '*'), ("TGG", 'W'),
        ("CGT", 'R'), ("CGC", 'R'), ("CGA", 'R'), ("CGG", 'R'),
        ("AGT", 'S'), ("AGC", 'S'), ("AGA", 'R'), ("AGG", 'R'),
        ("GGT", 'G'), ("GGC", 'G'), ("GGA", 'G'), ("GGG", 'G'),
    ]);
    let seq = seq.to_uppercase();
    let mut protein = String::new();
    let bytes = seq.as_bytes();
    let mut i = 0;
    while i + 3 <= bytes.len() {
        let codon = &seq[i..i + 3];
        if let Some(&aa) = codon_table.get(codon) {
            protein.push(aa);
            if aa == '*' {
                break;
            }
        } else {
            protein.push('X');
        }
        i += 3;
    }
    protein
}

/// 统计序列综合信息
#[pyfunction]
pub fn seq_stats(py: Python, seq: &str) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    let bases = count_bases(seq);
    dict.set_item("length", seq.len())?;
    dict.set_item("gc_content", gc_content(seq))?;
    dict.set_item("n_count", count_n(seq))?;
    dict.set_item("is_valid_dna", is_valid_dna(seq))?;
    dict.set_item("a_count", bases.get(&'A').unwrap_or(&0))?;
    dict.set_item("t_count", bases.get(&'T').unwrap_or(&0))?;
    dict.set_item("g_count", bases.get(&'G').unwrap_or(&0))?;
    dict.set_item("c_count", bases.get(&'C').unwrap_or(&0))?;
    Ok(dict.into())
}
