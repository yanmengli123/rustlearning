use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

/// Tn5 insertion site计算（ATAC-seq）
#[pyfunction]
pub fn tn5_insertion_sites(
    chrom: &str,
    pos: i64,
    is_reverse: bool,
) -> (String, i64) {
    let insertion = if is_reverse {
        pos - 1  // reverse strand: pos - 1
    } else {
        pos      // forward strand: pos
    };
    (chrom.to_string(), insertion)
}

/// Peak reads统计
#[pyfunction]
pub fn count_reads_in_peaks(
    read_positions: Vec<(String, i64)>,
    peaks: Vec<(String, i64, i64)>,
) -> HashMap<usize, i64> {
    let mut counts: HashMap<usize, i64> = HashMap::new();
    for (chrom, pos) in &read_positions {
        for (i, (pc, ps, pe)) in peaks.iter().enumerate() {
            if chrom == pc && *pos >= *ps && *pos < *pe {
                *counts.entry(i).or_insert(0) += 1;
            }
        }
    }
    counts
}

/// Fragment长度统计
#[pyfunction]
pub fn fragment_length_distribution(
    fragment_lengths: Vec<i64>,
) -> HashMap<i64, usize> {
    let mut dist: HashMap<i64, usize> = HashMap::new();
    for len in fragment_lengths {
        *dist.entry(len).or_insert(0) += 1;
    }
    dist
}

/// Nucleosome-free reads分类
#[pyfunction]
pub fn classify_nucleosome_free(
    fragment_lengths: Vec<i64>,
    max_nfr_len: i64,
) -> (usize, usize) {
    let nfr = fragment_lengths.iter().filter(|&&l| l <= max_nfr_len).count();
    let nucleosomal = fragment_lengths.len() - nfr;
    (nfr, nucleosomal)
}

/// ChIP-seq peak coverage计算
#[pyfunction]
pub fn chipseq_peak_coverage(
    py: Python,
    read_positions: Vec<(String, i64)>,
    peak_chrom: &str,
    peak_start: i64,
    peak_end: i64,
) -> PyResult<PyObject> {
    let size = (peak_end - peak_start) as usize;
    let mut coverage = vec![0i64; size];

    for (chrom, pos) in &read_positions {
        if chrom == peak_chrom && *pos >= peak_start && *pos < peak_end {
            let idx = (pos - peak_start) as usize;
            coverage[idx] += 1;
        }
    }

    let total: i64 = coverage.iter().sum();
    let max_cov = *coverage.iter().max().unwrap_or(&0);
    let mean_cov = if size > 0 {
        total as f64 / size as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("total_reads", total)?;
    dict.set_item("max_coverage", max_cov)?;
    dict.set_item("mean_coverage", mean_cov)?;
    dict.set_item("peak_length", size)?;
    Ok(dict.into())
}

/// CpG methylation level聚合
#[pyfunction]
pub fn cpg_methylation_levels(
    methylation_calls: Vec<(String, i64, bool)>,
    cpg_sites: Vec<(String, i64)>,
    max_distance: i64,
) -> HashMap<(String, i64), (usize, usize)> {
    let mut site_counts: HashMap<(String, i64), (usize, usize)> = HashMap::new();

    for (call_chrom, call_pos, is_methylated) in &methylation_calls {
        for (site_chrom, site_pos) in &cpg_sites {
            if call_chrom == site_chrom && (call_pos - site_pos).abs() <= max_distance {
                let entry = site_counts
                    .entry((site_chrom.clone(), *site_pos))
                    .or_insert((0, 0));
                if *is_methylated {
                    entry.0 += 1; // methylated
                }
                entry.1 += 1; // total
            }
        }
    }
    site_counts
}

/// Fragment overlap peaks
#[pyfunction]
pub fn fragment_overlap_peaks(
    fragments: Vec<(String, i64, i64)>,
    peaks: Vec<(String, i64, i64)>,
) -> Vec<Vec<usize>> {
    fragments
        .iter()
        .map(|(fc, fs, fe)| {
            peaks
                .iter()
                .enumerate()
                .filter(|(_, (pc, ps, pe))| fc == pc && *fs < *pe && *ps < *fe)
                .map(|(i, _)| i)
                .collect()
        })
        .collect()
}

/// Coverage over genomic bins
#[pyfunction]
pub fn coverage_over_bins(
    read_positions: Vec<(String, i64)>,
    bins: Vec<(String, i64, i64)>,
) -> Vec<i64> {
    bins.iter()
        .map(|(chrom, start, end)| {
            read_positions
                .iter()
                .filter(|(rc, rp)| rc == chrom && rp >= start && rp < end)
                .count() as i64
        })
        .collect()
}

/// CUT&Tag / CUT&RUN 片段统计
#[pyfunction]
pub fn cuttag_fragment_stats(
    py: Python,
    fragments: Vec<(String, i64, i64)>,
    peaks: Vec<(String, i64, i64)>,
) -> PyResult<PyObject> {
    let total_fragments = fragments.len();
    let mut in_peaks = 0usize;
    let mut fragment_lengths: Vec<i64> = Vec::new();

    for (fc, fs, fe) in &fragments {
        let len = fe - fs;
        fragment_lengths.push(len);
        for (pc, ps, pe) in &peaks {
            if fc == pc && *fs < *pe && *ps < *fe {
                in_peaks += 1;
                break;
            }
        }
    }

    fragment_lengths.sort_unstable();
    let median_len = if fragment_lengths.is_empty() {
        0
    } else {
        fragment_lengths[fragment_lengths.len() / 2]
    };
    let frip = if total_fragments > 0 {
        in_peaks as f64 / total_fragments as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("total_fragments", total_fragments)?;
    dict.set_item("fragments_in_peaks", in_peaks)?;
    dict.set_item("frip", frip)?;
    dict.set_item("median_fragment_length", median_len)?;
    Ok(dict.into())
}
