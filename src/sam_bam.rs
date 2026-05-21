use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// SAM记录结构
#[allow(dead_code)]
#[derive(Debug)]
pub struct SamRecord {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub pos: i64,
    pub mapq: u8,
    pub cigar: String,
    pub rnext: String,
    pub pnext: i64,
    pub tlen: i64,
    pub seq: String,
    pub qual: String,
    pub tags: HashMap<String, String>,
}

/// 解析SAM行
pub fn parse_sam_line(line: &str) -> Option<SamRecord> {
    if line.starts_with('@') {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return None;
    }

    let mut tags = HashMap::new();
    for field in &fields[11..] {
        if let Some(idx) = field.find(':') {
            let tag = &field[..idx];
            let val = &field[idx + 2..];
            tags.insert(tag.to_string(), val.to_string());
        }
    }

    Some(SamRecord {
        qname: fields[0].to_string(),
        flag: fields[1].parse().unwrap_or(0),
        rname: fields[2].to_string(),
        pos: fields[3].parse().unwrap_or(0),
        mapq: fields[4].parse().unwrap_or(0),
        cigar: fields[5].to_string(),
        rnext: fields[6].to_string(),
        pnext: fields[7].parse().unwrap_or(0),
        tlen: fields[8].parse().unwrap_or(0),
        seq: fields[9].to_string(),
        qual: fields[10].to_string(),
        tags,
    })
}

/// 读取SAM文件所有记录
fn read_sam(path: &str) -> io::Result<Vec<SamRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = parse_sam_line(&line) {
            records.push(rec);
        }
    }
    Ok(records)
}

/// 解析SAM flag
#[pyfunction]
pub fn parse_flag(py: Python, flag: u16) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    dict.set_item("read_paired", flag & 0x1 != 0)?;
    dict.set_item("read_mapped_in_pair", flag & 0x2 != 0)?;
    dict.set_item("read_unmapped", flag & 0x4 != 0)?;
    dict.set_item("mate_unmapped", flag & 0x8 != 0)?;
    dict.set_item("read_reverse", flag & 0x10 != 0)?;
    dict.set_item("mate_reverse", flag & 0x20 != 0)?;
    dict.set_item("first_in_pair", flag & 0x40 != 0)?;
    dict.set_item("second_in_pair", flag & 0x80 != 0)?;
    dict.set_item("not_primary", flag & 0x100 != 0)?;
    dict.set_item("read_fails_qc", flag & 0x200 != 0)?;
    dict.set_item("read_is_duplicate", flag & 0x400 != 0)?;
    dict.set_item("supplementary", flag & 0x800 != 0)?;
    Ok(dict.into())
}

/// SAM文件统计
#[pyfunction]
pub fn sam_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total = records.len();
    let mapped = records.iter().filter(|r| r.flag & 0x4 == 0).count();
    let unmapped = total - mapped;
    let paired = records.iter().filter(|r| r.flag & 0x1 != 0).count();
    let duplicates = records.iter().filter(|r| r.flag & 0x400 != 0).count();
    let primary = records.iter().filter(|r| r.flag & 0x100 == 0).count();
    let supplementary = records.iter().filter(|r| r.flag & 0x800 != 0).count();

    let mut chrom_counts: HashMap<String, usize> = HashMap::new();
    let mut mapq_sum: u64 = 0;
    let mut mapped_with_mapq = 0usize;
    for rec in &records {
        if rec.flag & 0x4 == 0 && rec.rname != "*" {
            *chrom_counts.entry(rec.rname.clone()).or_insert(0) += 1;
            mapq_sum += rec.mapq as u64;
            mapped_with_mapq += 1;
        }
    }

    let avg_mapq = if mapped_with_mapq > 0 {
        mapq_sum as f64 / mapped_with_mapq as f64
    } else {
        0.0
    };

    let dict = PyDict::new_bound(py);
    dict.set_item("total_reads", total)?;
    dict.set_item("mapped", mapped)?;
    dict.set_item("unmapped", unmapped)?;
    dict.set_item("paired", paired)?;
    dict.set_item("duplicates", duplicates)?;
    dict.set_item("primary", primary)?;
    dict.set_item("supplementary", supplementary)?;
    dict.set_item("avg_mapq", avg_mapq)?;
    dict.set_item("chrom_counts", chrom_counts)?;
    Ok(dict.into())
}

/// 按mapping quality过滤
#[pyfunction]
pub fn filter_by_mapq(path: &str, min_mapq: u8) -> PyResult<Vec<String>> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let filtered: Vec<String> = records
        .iter()
        .filter(|r| r.flag & 0x4 == 0 && r.mapq >= min_mapq)
        .map(|r| r.qname.clone())
        .collect();
    Ok(filtered)
}

/// 按染色体区间提取reads
#[pyfunction]
pub fn fetch_region(path: &str, chrom: &str, start: i64, end: i64) -> PyResult<Vec<String>> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let hits: Vec<String> = records
        .iter()
        .filter(|r| {
            r.rname == chrom
                && r.pos >= start
                && r.pos <= end
                && r.flag & 0x4 == 0
        })
        .map(|r| format!("{}\t{}\t{}", r.qname, r.pos, r.cigar))
        .collect();
    Ok(hits)
}

/// Coverage统计（按染色体位置）
#[pyfunction]
pub fn coverage_at_position(path: &str, chrom: &str, position: i64) -> PyResult<usize> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut count = 0usize;
    for rec in &records {
        if rec.rname != chrom || rec.flag & 0x4 != 0 {
            continue;
        }
        let read_len = rec.seq.len() as i64;
        if rec.pos <= position && rec.pos + read_len > position {
            count += 1;
        }
    }
    Ok(count)
}

/// 区间coverage（返回每个位置的coverage）
#[pyfunction]
pub fn region_coverage(
    path: &str,
    chrom: &str,
    start: i64,
    end: i64,
) -> PyResult<Vec<i64>> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let size = (end - start) as usize;
    let mut cov = vec![0i64; size];

    for rec in &records {
        if rec.rname != chrom || rec.flag & 0x4 != 0 {
            continue;
        }
        let read_start = rec.pos;
        let read_end = rec.pos + rec.seq.len() as i64;
        let overlap_start = read_start.max(start);
        let overlap_end = read_end.min(end);
        if overlap_start < overlap_end {
            for pos in overlap_start..overlap_end {
                cov[(pos - start) as usize] += 1;
            }
        }
    }
    Ok(cov)
}

/// Insert size统计（paired-end）
#[pyfunction]
pub fn insert_size_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let records = read_sam(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut sizes: Vec<i64> = Vec::new();
    for rec in &records {
        if rec.flag & 0x1 != 0 && rec.flag & 0x4 == 0 && rec.tlen > 0 {
            sizes.push(rec.tlen.abs());
        }
    }
    sizes.sort_unstable();

    let mean = if sizes.is_empty() {
        0.0
    } else {
        sizes.iter().sum::<i64>() as f64 / sizes.len() as f64
    };
    let median = if sizes.is_empty() {
        0.0
    } else {
        sizes[sizes.len() / 2] as f64
    };
    let min = sizes.first().copied().unwrap_or(0);
    let max = sizes.last().copied().unwrap_or(0);

    let dict = PyDict::new_bound(py);
    dict.set_item("count", sizes.len())?;
    dict.set_item("mean", mean)?;
    dict.set_item("median", median)?;
    dict.set_item("min", min)?;
    dict.set_item("max", max)?;
    Ok(dict.into())
}
