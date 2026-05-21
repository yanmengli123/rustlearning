use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// 从GTF构建gene区间索引
pub fn build_gene_index(path: &str) -> io::Result<HashMap<String, Vec<(String, i64, i64, String)>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut genes: HashMap<String, Vec<(String, i64, i64, String)>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == "gene" {
                let gene_id = rec
                    .attributes
                    .get("gene_id")
                    .cloned()
                    .unwrap_or_default();
                genes
                    .entry(rec.chrom.clone())
                    .or_default()
                    .push((rec.chrom, rec.start, rec.end, gene_id));
            }
        }
    }
    Ok(genes)
}

/// 基因计数矩阵构建
#[pyfunction]
pub fn gene_count_matrix(
    gtf_path: &str,
    read_positions: Vec<(String, i64)>,
) -> PyResult<HashMap<String, i64>> {
    let gene_index = build_gene_index(gtf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut counts: HashMap<String, i64> = HashMap::new();

    for (chrom, pos) in &read_positions {
        if let Some(genes) = gene_index.get(chrom) {
            for (_, start, end, gene_id) in genes {
                if pos >= start && pos < end {
                    *counts.entry(gene_id.clone()).or_insert(0) += 1;
                }
            }
        }
    }
    Ok(counts)
}

/// Exon计数
#[pyfunction]
pub fn exon_count(
    gtf_path: &str,
    read_positions: Vec<(String, i64)>,
) -> PyResult<HashMap<String, i64>> {
    let file = File::open(gtf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let mut exons: Vec<(String, i64, i64, String)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == "exon" {
                let exon_id = rec
                    .attributes
                    .get("exon_id")
                    .cloned()
                    .unwrap_or_default();
                exons.push((rec.chrom, rec.start, rec.end, exon_id));
            }
        }
    }

    let mut counts: HashMap<String, i64> = HashMap::new();
    for (chrom, pos) in &read_positions {
        for (ec, es, ee, eid) in &exons {
            if chrom == ec && pos >= es && pos < ee {
                *counts.entry(eid.clone()).or_insert(0) += 1;
            }
        }
    }
    Ok(counts)
}

/// Splice junction统计
#[pyfunction]
pub fn splice_junction_stats(cigar: &str, pos: i64) -> Vec<(String, i64, i64)> {
    let ops = super::alignment::parse_cigar(cigar);
    let mut junctions = Vec::new();
    let mut current_pos = pos;

    for (op, len) in &ops {
        match op {
            'M' | '=' | 'X' => current_pos += *len as i64,
            'N' => {
                let junction_start = current_pos;
                let junction_end = current_pos + *len as i64;
                junctions.push((
                    format!("{}-{}", junction_start, junction_end),
                    junction_start,
                    junction_end,
                ));
                current_pos += *len as i64;
            }
            'D' => current_pos += *len as i64,
            _ => {}
        }
    }
    junctions
}

/// Intron retention检测
#[pyfunction]
pub fn detect_intron_retention(
    cigar: &str,
    pos: i64,
    intron_start: i64,
    intron_end: i64,
) -> bool {
    let ops = super::alignment::parse_cigar(cigar);
    let mut current_pos = pos;

    for (op, len) in &ops {
        match op {
            'M' | '=' | 'X' => {
                let read_end = current_pos + *len as i64;
                // 如果read覆盖了整个intron区域，认为是intron retention
                if current_pos <= intron_start && read_end >= intron_end {
                    return true;
                }
                current_pos += *len as i64;
            }
            'N' => current_pos += *len as i64,
            'D' => current_pos += *len as i64,
            _ => {}
        }
    }
    false
}

/// UMI collapsing
#[pyfunction]
pub fn umi_collapse(
    umis: Vec<String>,
    max_distance: usize,
) -> Vec<(String, Vec<usize>)> {
    let mut groups: Vec<(String, Vec<usize>)> = Vec::new();

    for (i, umi) in umis.iter().enumerate() {
        let mut merged = false;
        for group in groups.iter_mut() {
            let dist = super::alignment::levenshtein_distance(&group.0, umi);
            if dist <= max_distance {
                group.1.push(i);
                merged = true;
                break;
            }
        }
        if !merged {
            groups.push((umi.clone(), vec![i]));
        }
    }
    groups
}

/// Gene biotype统计
#[pyfunction]
pub fn gene_biotype_stats(gtf_path: &str) -> PyResult<HashMap<String, usize>> {
    let file = File::open(gtf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let mut biotypes: HashMap<String, usize> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == "gene" {
                let biotype = rec
                    .attributes
                    .get("gene_biotype")
                    .cloned()
                    .unwrap_or_else(|| "unknown".to_string());
                *biotypes.entry(biotype).or_insert(0) += 1;
            }
        }
    }
    Ok(biotypes)
}

/// Transcript length统计
#[pyfunction]
pub fn transcript_length_stats(gtf_path: &str) -> PyResult<HashMap<String, i64>> {
    let file = File::open(gtf_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let mut exon_lengths: HashMap<String, i64> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == "exon" {
                let tx_id = rec
                    .attributes
                    .get("transcript_id")
                    .cloned()
                    .unwrap_or_default();
                *exon_lengths.entry(tx_id).or_insert(0) += rec.end - rec.start;
            }
        }
    }
    Ok(exon_lengths)
}
