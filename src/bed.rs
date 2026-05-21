use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// BED区间
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BedInterval {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub name: String,
    pub score: f64,
    pub strand: char,
}

/// 解析BED行
pub fn parse_bed_line(line: &str) -> Option<BedInterval> {
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 3 {
        return None;
    }
    Some(BedInterval {
        chrom: fields[0].to_string(),
        start: fields[1].parse().unwrap_or(0),
        end: fields[2].parse().unwrap_or(0),
        name: fields.get(3).unwrap_or(&".").to_string(),
        score: fields.get(4).unwrap_or(&"0").parse().unwrap_or(0.0),
        strand: fields.get(5).unwrap_or(&".").chars().next().unwrap_or('.'),
    })
}

/// 读取BED文件
pub fn read_bed(path: &str) -> io::Result<Vec<BedInterval>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(iv) = parse_bed_line(&line) {
            intervals.push(iv);
        }
    }
    Ok(intervals)
}

/// BED文件统计
#[pyfunction]
pub fn bed_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let intervals = read_bed(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total = intervals.len();
    let mut chrom_counts: HashMap<String, usize> = HashMap::new();
    let mut total_bases: i64 = 0;
    let mut min_len = i64::MAX;
    let mut max_len = 0i64;

    for iv in &intervals {
        *chrom_counts.entry(iv.chrom.clone()).or_insert(0) += 1;
        let len = iv.end - iv.start;
        total_bases += len;
        min_len = min_len.min(len);
        max_len = max_len.max(len);
    }

    let avg_len = if total > 0 {
        total_bases as f64 / total as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("total_intervals", total)?;
    dict.set_item("total_bases", total_bases)?;
    dict.set_item("avg_length", avg_len)?;
    dict.set_item("min_length", if min_len == i64::MAX { 0 } else { min_len })?;
    dict.set_item("max_length", max_len)?;
    dict.set_item("chrom_counts", chrom_counts)?;
    Ok(dict.into())
}

/// 区间交集
#[pyfunction]
pub fn bed_intersect(path_a: &str, path_b: &str) -> PyResult<Vec<String>> {
    let a = read_bed(path_a)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let b = read_bed(path_b)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut results = Vec::new();
    for iv_a in &a {
        for iv_b in &b {
            if iv_a.chrom == iv_b.chrom {
                let overlap_start = iv_a.start.max(iv_b.start);
                let overlap_end = iv_a.end.min(iv_b.end);
                if overlap_start < overlap_end {
                    results.push(format!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        iv_a.chrom, overlap_start, overlap_end, iv_a.name, iv_b.name,
                        overlap_end - overlap_start
                    ));
                }
            }
        }
    }
    Ok(results)
}

/// 区间并集（合并重叠区间）
#[pyfunction]
pub fn bed_merge(path: &str) -> PyResult<Vec<String>> {
    let mut intervals = read_bed(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    intervals.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    let mut merged: Vec<BedInterval> = Vec::new();
    for iv in intervals {
        if let Some(last) = merged.last_mut() {
            if last.chrom == iv.chrom && last.end >= iv.start {
                last.end = last.end.max(iv.end);
                continue;
            }
        }
        merged.push(iv);
    }

    let result: Vec<String> = merged
        .iter()
        .map(|iv| format!("{}\t{}\t{}", iv.chrom, iv.start, iv.end))
        .collect();
    Ok(result)
}

/// 区间差集（A中不与B重叠的部分）
#[pyfunction]
pub fn bed_subtract(path_a: &str, path_b: &str) -> PyResult<Vec<String>> {
    let a = read_bed(path_a)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let b = read_bed(path_b)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut results = Vec::new();
    for iv_a in &a {
        let mut remaining: Vec<(i64, i64)> = vec![(iv_a.start, iv_a.end)];
        for iv_b in &b {
            if iv_a.chrom != iv_b.chrom {
                continue;
            }
            let mut new_remaining = Vec::new();
            for (s, e) in remaining {
                if iv_b.end <= s || iv_b.start >= e {
                    new_remaining.push((s, e));
                } else {
                    if iv_b.start > s {
                        new_remaining.push((s, iv_b.start));
                    }
                    if iv_b.end < e {
                        new_remaining.push((iv_b.end, e));
                    }
                }
            }
            remaining = new_remaining;
        }
        for (s, e) in remaining {
            results.push(format!("{}\t{}\t{}", iv_a.chrom, s, e));
        }
    }
    Ok(results)
}

/// 最近特征查找
#[pyfunction]
pub fn bed_closest(path_a: &str, path_b: &str) -> PyResult<Vec<String>> {
    let a = read_bed(path_a)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let b = read_bed(path_b)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut results = Vec::new();
    for iv_a in &a {
        let mut min_dist = i64::MAX;
        let mut closest_name = String::new();
        for iv_b in &b {
            if iv_a.chrom != iv_b.chrom {
                continue;
            }
            let dist = if iv_a.end <= iv_b.start {
                iv_b.start - iv_a.end
            } else if iv_b.end <= iv_a.start {
                iv_a.start - iv_b.end
            } else {
                0
            };
            if dist < min_dist {
                min_dist = dist;
                closest_name = iv_b.name.clone();
            }
        }
        results.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            iv_a.chrom, iv_a.start, iv_a.end, iv_a.name, closest_name, min_dist
        ));
    }
    Ok(results)
}

/// VCF记录
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: i64,
    pub id: String,
    pub ref_allele: String,
    pub alt_allele: String,
    pub qual: f64,
    pub filter: String,
    pub info: HashMap<String, String>,
}

/// 解析VCF行
pub fn parse_vcf_line(line: &str) -> Option<VcfRecord> {
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return None;
    }
    let mut info = HashMap::new();
    if fields[7] != "." {
        for item in fields[7].split(';') {
            if let Some(idx) = item.find('=') {
                info.insert(item[..idx].to_string(), item[idx + 1..].to_string());
            } else {
                info.insert(item.to_string(), "true".to_string());
            }
        }
    }
    Some(VcfRecord {
        chrom: fields[0].to_string(),
        pos: fields[1].parse().unwrap_or(0),
        id: fields[2].to_string(),
        ref_allele: fields[3].to_string(),
        alt_allele: fields[4].to_string(),
        qual: fields[5].parse().unwrap_or(0.0),
        filter: fields[6].to_string(),
        info,
    })
}

/// 读取VCF
pub fn read_vcf(path: &str) -> io::Result<Vec<VcfRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = parse_vcf_line(&line) {
            records.push(rec);
        }
    }
    Ok(records)
}

/// VCF统计
#[pyfunction]
pub fn vcf_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let records = read_vcf(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total = records.len();
    let mut snps = 0usize;
    let mut indels = 0usize;
    let mut chrom_counts: HashMap<String, usize> = HashMap::new();

    for rec in &records {
        *chrom_counts.entry(rec.chrom.clone()).or_insert(0) += 1;
        if rec.ref_allele.len() == 1 && rec.alt_allele.len() == 1 {
            snps += 1;
        } else {
            indels += 1;
        }
    }

    let dict = PyDict::new_bound(py);
    dict.set_item("total_variants", total)?;
    dict.set_item("snps", snps)?;
    dict.set_item("indels", indels)?;
    dict.set_item("chrom_counts", chrom_counts)?;
    Ok(dict.into())
}

/// VCF过滤
#[pyfunction]
pub fn vcf_filter(path: &str, min_qual: f64) -> PyResult<Vec<String>> {
    let records = read_vcf(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let filtered: Vec<String> = records
        .iter()
        .filter(|r| r.qual >= min_qual)
        .map(|r| {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{:.1}",
                r.chrom, r.pos, r.id, r.ref_allele, r.alt_allele, r.qual
            )
        })
        .collect();
    Ok(filtered)
}

/// GTF记录
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GtfRecord {
    pub chrom: String,
    pub source: String,
    pub feature: String,
    pub start: i64,
    pub end: i64,
    pub score: String,
    pub strand: char,
    pub frame: String,
    pub attributes: HashMap<String, String>,
}

/// 解析GTF行
pub fn parse_gtf_line(line: &str) -> Option<GtfRecord> {
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }
    let mut attrs = HashMap::new();
    for item in fields[8].split(';') {
        let item = item.trim();
        if item.is_empty() {
            continue;
        }
        let parts: Vec<&str> = item.splitn(2, ' ').collect();
        if parts.len() == 2 {
            let val = parts[1].trim_matches('"');
            attrs.insert(parts[0].to_string(), val.to_string());
        }
    }
    Some(GtfRecord {
        chrom: fields[0].to_string(),
        source: fields[1].to_string(),
        feature: fields[2].to_string(),
        start: fields[3].parse().unwrap_or(0),
        end: fields[4].parse().unwrap_or(0),
        score: fields[5].to_string(),
        strand: fields[6].chars().next().unwrap_or('.'),
        frame: fields[7].to_string(),
        attributes: attrs,
    })
}

/// GTF统计
#[pyfunction]
pub fn gtf_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let file = File::open(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let mut feature_counts: HashMap<String, usize> = HashMap::new();
    let mut gene_count = 0usize;
    let mut transcript_count = 0usize;

    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = parse_gtf_line(&line) {
            *feature_counts.entry(rec.feature.clone()).or_insert(0) += 1;
            match rec.feature.as_str() {
                "gene" => gene_count += 1,
                "transcript" => transcript_count += 1,
                _ => {}
            }
        }
    }

    let dict = PyDict::new_bound(py);
    dict.set_item("gene_count", gene_count)?;
    dict.set_item("transcript_count", transcript_count)?;
    dict.set_item("feature_counts", feature_counts)?;
    Ok(dict.into())
}
