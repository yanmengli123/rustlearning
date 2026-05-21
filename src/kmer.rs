use pyo3::prelude::*;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

/// 计算k-mer计数
#[pyfunction]
pub fn count_kmers(seq: &str, k: usize) -> PyResult<HashMap<String, usize>> {
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
        *kmers.entry(kmer).or_insert(0) += 1;
    }
    Ok(kmers)
}

/// 计算canonical k-mer（取lexicographically较小的）
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
        let rc = super::sequence::reverse_complement(&kmer);
        let canonical = if kmer <= rc { kmer } else { rc };
        *kmers.entry(canonical).or_insert(0) += 1;
    }
    Ok(kmers)
}

/// Minimizer提取
#[pyfunction]
pub fn minimizers(seq: &str, k: usize, w: usize) -> PyResult<Vec<(usize, String)>> {
    if k == 0 || w == 0 || k > seq.len() || w > seq.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Invalid k or w",
        ));
    }
    let seq = seq.to_uppercase();
    let mut result = Vec::new();
    let mut last_min = String::new();

    for i in 0..=seq.len().saturating_sub(w + k - 1) {
        let mut min_kmer = String::new();
        let mut min_pos = 0;
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
        if min_kmer != last_min {
            result.push((min_pos, min_kmer.clone()));
            last_min = min_kmer;
        }
    }
    Ok(result)
}

/// MinHash sketch
#[pyfunction]
pub fn minhash_sketch(seq: &str, k: usize, num_hashes: usize) -> PyResult<Vec<u64> > {
    if k == 0 || k > seq.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "k must be > 0 and <= seq length",
        ));
    }
    let seq = seq.to_uppercase();
    let bytes = seq.as_bytes();
    let mut hashes: Vec<u64> = Vec::new();

    for i in 0..=seq.len() - k {
        let kmer = &bytes[i..i + k];
        let hash = hash_kmer(kmer);
        hashes.push(hash);
    }
    hashes.sort_unstable();
    hashes.dedup();
    hashes.truncate(num_hashes);
    Ok(hashes)
}

fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    kmer.hash(&mut hasher);
    hasher.finish()
}

/// Jaccard相似度
#[pyfunction]
pub fn jaccard_similarity(seq1: &str, seq2: &str, k: usize) -> PyResult<f64> {
    let kmers1 = count_kmers(seq1, k)?;
    let kmers2 = count_kmers(seq2, k)?;

    let set1: std::collections::HashSet<&String> = kmers1.keys().collect();
    let set2: std::collections::HashSet<&String> = kmers2.keys().collect();

    let intersection = set1.intersection(&set2).count();
    let union = set1.union(&set2).count();

    if union == 0 {
        return Ok(0.0);
    }
    Ok(intersection as f64 / union as f64)
}

/// Containment index
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

/// Mash-like distance估算
#[pyfunction]
pub fn mash_distance(seq1: &str, seq2: &str, k: usize, sketch_size: usize) -> PyResult<f64> {
    let sketch1 = minhash_sketch(seq1, k, sketch_size)?;
    let sketch2 = minhash_sketch(seq2, k, sketch_size)?;

    let set1: std::collections::HashSet<u64> = sketch1.into_iter().collect();
    let set2: std::collections::HashSet<u64> = sketch2.into_iter().collect();

    let intersection = set1.intersection(&set2).count();
    let union = set1.union(&set2).count();

    if union == 0 {
        return Ok(1.0);
    }
    let jaccard = intersection as f64 / union as f64;
    if jaccard <= 0.0 {
        return Ok(1.0);
    }
    let distance = -1.0 / k as f64 * jaccard.ln();
    Ok(distance.min(1.0))
}

/// k-mer频谱
#[pyfunction]
pub fn kmer_spectrum(seq: &str, k: usize) -> PyResult<HashMap<usize, usize>> {
    let kmers = count_kmers(seq, k)?;
    let mut spectrum: HashMap<usize, usize> = HashMap::new();
    for (_, &count) in &kmers {
        *spectrum.entry(count).or_insert(0) += 1;
    }
    Ok(spectrum)
}

/// Syncmer提取
#[pyfunction]
pub fn syncmers(seq: &str, k: usize, s: usize) -> PyResult<Vec<(usize, String)>> {
    if s >= k || k > seq.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Invalid k or s (s must be < k)",
        ));
    }
    let seq = seq.to_uppercase();
    let mut result = Vec::new();

    for i in 0..=seq.len() - k {
        let kmer = &seq[i..i + k];
        let mut min_smer = String::new();
        for j in 0..=k - s {
            let smer = &kmer[j..j + s];
            if min_smer.is_empty() || smer < min_smer.as_str() {
                min_smer = smer.to_string();
            }
        }
        // check if first or last s-mer is the minimum
        let first_smer = &kmer[0..s];
        let last_smer = &kmer[k - s..k];
        if first_smer == min_smer || last_smer == min_smer {
            result.push((i, kmer.to_string()));
        }
    }
    Ok(result)
}
