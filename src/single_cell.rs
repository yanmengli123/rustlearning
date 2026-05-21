use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

/// Cell barcode提取
#[pyfunction]
pub fn extract_barcode(read: &str, barcode_len: usize, barcode_start: usize) -> Option<String> {
    if barcode_start + barcode_len > read.len() {
        return None;
    }
    Some(read[barcode_start..barcode_start + barcode_len].to_string())
}

/// UMI提取
#[pyfunction]
pub fn extract_umi(read: &str, umi_len: usize, umi_start: usize) -> Option<String> {
    if umi_start + umi_len > read.len() {
        return None;
    }
    Some(read[umi_start..umi_start + umi_len].to_string())
}

/// Barcode correction（白名单匹配）
#[pyfunction]
pub fn correct_barcode(
    barcode: &str,
    whitelist: Vec<String>,
    max_distance: usize,
) -> Option<String> {
    let mut best_match = None;
    let mut best_dist = max_distance + 1;

    for wb in &whitelist {
        let dist = super::alignment::levenshtein_distance(barcode, wb);
        if dist < best_dist {
            best_dist = dist;
            best_match = Some(wb.clone());
        }
        if dist == 0 {
            return best_match;
        }
    }
    if best_dist <= max_distance {
        best_match
    } else {
        None
    }
}

/// Barcode rescue（Hamming距离）
#[pyfunction]
pub fn barcode_rescue(
    barcode: &str,
    whitelist: Vec<String>,
    max_hamming: usize,
) -> Vec<(String, usize)> {
    let mut matches = Vec::new();
    for wb in &whitelist {
        if barcode.len() != wb.len() {
            continue;
        }
        let dist = super::alignment::hamming_distance(barcode, wb).unwrap_or(999);
        if dist <= max_hamming {
            matches.push((wb.clone(), dist));
        }
    }
    matches
}

/// UMI deduplication
#[pyfunction]
pub fn umi_dedup(
    umis: Vec<String>,
    max_edit_distance: usize,
) -> Vec<(String, Vec<usize>)> {
    super::rnaseq::umi_collapse(umis, max_edit_distance)
}

/// Sparse matrix条目
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SparseEntry {
    pub row: usize,
    pub col: usize,
    pub value: f64,
}

/// Feature-barcode matrix构建
#[pyfunction]
pub fn build_feature_barcode_matrix(
    py: Python,
    features: Vec<String>,
    barcodes: Vec<String>,
    entries: Vec<(usize, usize, f64)>,
) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    dict.set_item("features", &features)?;
    dict.set_item("barcodes", &barcodes)?;

    let sparse_data: Vec<PyObject> = entries
        .iter()
        .map(|(r, c, v)| {
            let d = PyDict::new_bound(py);
            d.set_item("row", *r).unwrap();
            d.set_item("col", *c).unwrap();
            d.set_item("value", *v).unwrap();
            d.into()
        })
        .collect();

    dict.set_item("entries", sparse_data)?;
    dict.set_item("n_features", features.len())?;
    dict.set_item("n_barcodes", barcodes.len())?;
    Ok(dict.into())
}

/// Matrix Market格式输出
#[pyfunction]
pub fn to_matrix_market(
    n_rows: usize,
    n_cols: usize,
    entries: Vec<(usize, usize, f64)>,
) -> String {
    let mut lines = vec![
        "%%MatrixMarket matrix coordinate real general".to_string(),
        format!("{} {} {}", n_rows, n_cols, entries.len()),
    ];
    for (row, col, val) in entries {
        lines.push(format!("{} {} {}", row + 1, col + 1, val));
    }
    lines.join("\n")
}

/// Cell QC指标计算
#[pyfunction]
pub fn cell_qc(
    py: Python,
    counts_per_cell: Vec<usize>,
    genes_per_cell: Vec<usize>,
) -> PyResult<PyObject> {
    let total_counts: usize = counts_per_cell.iter().sum();
    let total_genes: usize = genes_per_cell.iter().sum();
    let n_cells = counts_per_cell.len();
    let avg_counts = if n_cells > 0 {
        total_counts as f64 / n_cells as f64
    } else {
        0.0
    };
    let avg_genes = if n_cells > 0 {
        total_genes as f64 / n_cells as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("n_cells", n_cells)?;
    dict.set_item("total_counts", total_counts)?;
    dict.set_item("total_genes", total_genes)?;
    dict.set_item("avg_counts_per_cell", avg_counts)?;
    dict.set_item("avg_genes_per_cell", avg_genes)?;
    dict.set_item("max_counts", counts_per_cell.iter().max().unwrap_or(&0))?;
    dict.set_item("min_counts", counts_per_cell.iter().min().unwrap_or(&0))?;
    Ok(dict.into())
}

/// Antibody/hashtag demultiplexing
#[pyfunction]
pub fn demux_tags(
    tag_counts: Vec<Vec<f64>>,
    tag_names: Vec<String>,
    threshold: f64,
) -> Vec<String> {
    let mut assignments = Vec::new();
    for cell_tags in &tag_counts {
        let mut max_tag = "unknown".to_string();
        let mut max_val = 0.0f64;
        for (i, &val) in cell_tags.iter().enumerate() {
            if val > max_val && val > threshold {
                max_val = val;
                if i < tag_names.len() {
                    max_tag = tag_names[i].clone();
                }
            }
        }
        assignments.push(max_tag);
    }
    assignments
}

/// CRISPR guide计数
#[pyfunction]
pub fn count_guides(
    guide_assignments: Vec<String>,
) -> HashMap<String, usize> {
    let mut counts: HashMap<String, usize> = HashMap::new();
    for guide in guide_assignments {
        *counts.entry(guide).or_insert(0) += 1;
    }
    counts
}
