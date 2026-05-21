use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

/// k-mer分类器（简化版）
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

    for (kmer, taxa) in &reference_kmers {
        if query_set.contains(kmer) {
            for taxon in taxa {
                *taxon_hits.entry(taxon.clone()).or_insert(0) += 1;
                total_hits += 1;
            }
        }
    }

    let mut abundances: HashMap<String, f64> = HashMap::new();
    if total_hits > 0 {
        for (taxon, hits) in &taxon_hits {
            abundances.insert(taxon.clone(), *hits as f64 / total_hits as f64);
        }
    }
    abundances
}

/// ANI近似计算（Average Nucleotide Identity）
#[pyfunction]
pub fn ani_approximate(seq1: &str, seq2: &str, k: usize, sketch_size: usize) -> PyResult<f64> {
    let sketch1 = super::kmer::minhash_sketch(seq1, k, sketch_size)?;
    let sketch2 = super::kmer::minhash_sketch(seq2, k, sketch_size)?;

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
    Ok(ani.max(0.0).min(1.0))
}

/// Genome sketching
#[pyfunction]
pub fn genome_sketch(
    seq: &str,
    k: usize,
    sketch_size: usize,
) -> PyResult<Vec<u64>> {
    super::kmer::minhash_sketch(seq, k, sketch_size)
}

/// 16S序列预处理（去引物、清洗）
#[pyfunction]
pub fn preprocess_16s(
    seq: &str,
    forward_primer: &str,
    reverse_primer: &str,
    max_mismatch: usize,
) -> Option<String> {
    let fwd_hits = super::alignment::primer_match(seq, forward_primer, max_mismatch);
    let rev_rc = super::sequence::reverse_complement(reverse_primer);
    let rev_hits = super::alignment::primer_match(seq, &rev_rc, max_mismatch);

    if let (Some((fwd_pos, _)), Some((rev_pos, _))) = (fwd_hits.first(), rev_hits.first()) {
        let start = fwd_pos + forward_primer.len();
        let end = *rev_pos;
        if start < end && end <= seq.len() {
            return Some(seq[start..end].to_string());
        }
    }
    None
}

/// Marker gene匹配
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

/// Contamination screening
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

/// OTU/ASV聚类辅助（简单 greedy clustering）
#[pyfunction]
pub fn greedy_cluster(
    sequences: Vec<String>,
    identity_threshold: f64,
    k: usize,
) -> Vec<(String, Vec<usize>)> {
    let mut clusters: Vec<(String, Vec<usize>)> = Vec::new();

    for (i, seq) in sequences.iter().enumerate() {
        let mut assigned = false;
        for cluster in clusters.iter_mut() {
            let sim = super::kmer::jaccard_similarity(&cluster.0, seq, k).unwrap_or(0.0);
            if sim >= identity_threshold {
                cluster.1.push(i);
                assigned = true;
                break;
            }
        }
        if !assigned {
            clusters.push((seq.clone(), vec![i]));
        }
    }
    clusters
}

/// Abundance estimation辅助（基于k-mer覆盖度）
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
