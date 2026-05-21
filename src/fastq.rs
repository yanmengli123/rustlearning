use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};

/// 单条FASTQ记录
struct FastqRecord {
    header: String,
    seq: String,
    qual: String,
}

/// 解析FASTQ文件，返回记录迭代器
fn parse_fastq_records(path: &str) -> io::Result<Vec<FastqRecord>> {
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        let gz = flate2::read::GzDecoder::new(file);
        Box::new(BufReader::new(gz))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut records = Vec::new();
    let mut lines = reader.lines();

    while let Some(Ok(header)) = lines.next() {
        if !header.starts_with('@') {
            continue;
        }
        if let (Some(Ok(seq)), Some(Ok(_plus)), Some(Ok(qual))) =
            (lines.next(), lines.next(), lines.next())
        {
            records.push(FastqRecord {
                header: header[1..].to_string(),
                seq,
                qual,
            });
        }
    }
    Ok(records)
}

/// 计算Phred质量分数
fn qual_scores(qual: &str) -> Vec<u8> {
    qual.bytes().map(|b| b - 33).collect()
}

/// FASTQ质控统计
#[pyfunction]
pub fn fastq_stats(py: Python, path: &str) -> PyResult<PyObject> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total_reads = records.len();
    let mut total_bases: u64 = 0;
    let mut q20_bases: u64 = 0;
    let mut q30_bases: u64 = 0;
    let mut gc_count: u64 = 0;
    let mut n_count: u64 = 0;
    let mut length_sum: u64 = 0;
    let mut min_len = usize::MAX;
    let mut max_len = 0usize;
    let mut total_qual_sum: f64 = 0.0;
    let mut per_base_qual: HashMap<usize, Vec<f64>> = HashMap::new();

    for rec in &records {
        let slen = rec.seq.len();
        total_bases += slen as u64;
        length_sum += slen as u64;
        min_len = min_len.min(slen);
        max_len = max_len.max(slen);

        let scores = qual_scores(&rec.qual);
        for (i, &q) in scores.iter().enumerate() {
            total_qual_sum += q as f64;
            per_base_qual.entry(i).or_default().push(q as f64);
            if q >= 20 {
                q20_bases += 1;
            }
            if q >= 30 {
                q30_bases += 1;
            }
        }
        for c in rec.seq.chars() {
            match c {
                'G' | 'g' | 'C' | 'c' => gc_count += 1,
                'N' | 'n' => n_count += 1,
                _ => {}
            }
        }
    }

    let avg_len = if total_reads > 0 {
        length_sum as f64 / total_reads as f64
    } else {
        0.0
    };
    let avg_qual = if total_bases > 0 {
        total_qual_sum / total_bases as f64
    } else {
        0.0
    };
    let q20_rate = if total_bases > 0 {
        q20_bases as f64 / total_bases as f64
    } else {
        0.0
    };
    let q30_rate = if total_bases > 0 {
        q30_bases as f64 / total_bases as f64
    } else {
        0.0
    };
    let gc_rate = if total_bases > 0 {
        gc_count as f64 / total_bases as f64
    } else {
        0.0
    };

    // 计算每个位点的平均质量
    let mut per_base_mean_qual: Vec<f64> = Vec::new();
    let max_pos = per_base_qual.keys().max().copied().unwrap_or(0);
    for i in 0..=max_pos {
        if let Some(scores) = per_base_qual.get(&i) {
            let mean = scores.iter().sum::<f64>() / scores.len() as f64;
            per_base_mean_qual.push(mean);
        }
    }

    let dict = PyDict::new_bound(py);
    dict.set_item("total_reads", total_reads)?;
    dict.set_item("total_bases", total_bases)?;
    dict.set_item("avg_length", avg_len)?;
    dict.set_item("min_length", if min_len == usize::MAX { 0 } else { min_len })?;
    dict.set_item("max_length", max_len)?;
    dict.set_item("avg_qual", avg_qual)?;
    dict.set_item("q20_bases", q20_bases)?;
    dict.set_item("q30_bases", q30_bases)?;
    dict.set_item("q20_rate", q20_rate)?;
    dict.set_item("q30_rate", q30_rate)?;
    dict.set_item("gc_content", gc_rate)?;
    dict.set_item("n_bases", n_count)?;
    Ok(dict.into())
}

/// FASTQ过滤
#[pyfunction]
pub fn fastq_filter(
    input: &str,
    output: &str,
    min_len: usize,
    min_qual: f64,
) -> PyResult<usize> {
    let records = parse_fastq_records(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(file);

    let mut kept = 0usize;
    for rec in &records {
        if rec.seq.len() < min_len {
            continue;
        }
        let scores = qual_scores(&rec.qual);
        let mean_q = scores.iter().map(|&q| q as f64).sum::<f64>() / scores.len() as f64;
        if mean_q < min_qual {
            continue;
        }
        writeln!(writer, "@{}", rec.header).unwrap();
        writeln!(writer, "{}", rec.seq).unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "{}", rec.qual).unwrap();
        kept += 1;
    }
    Ok(kept)
}

/// 每个位点的质量分布
#[pyfunction]
pub fn per_base_quality(py: Python, path: &str) -> PyResult<PyObject> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut per_pos: HashMap<usize, Vec<f64>> = HashMap::new();
    for rec in &records {
        let scores = qual_scores(&rec.qual);
        for (i, &q) in scores.iter().enumerate() {
            per_pos.entry(i).or_default().push(q as f64);
        }
    }

    let dict = PyDict::new_bound(py);
    let mut positions: Vec<usize> = per_pos.keys().copied().collect();
    positions.sort();
    let mut means = Vec::new();
    let mut medians = Vec::new();
    for &pos in &positions {
        let mut scores = per_pos[&pos].clone();
        scores.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mean = scores.iter().sum::<f64>() / scores.len() as f64;
        let median = scores[scores.len() / 2];
        means.push(mean);
        medians.push(median);
    }
    dict.set_item("positions", positions)?;
    dict.set_item("means", means)?;
    dict.set_item("medians", medians)?;
    Ok(dict.into())
}

/// Read长度分布
#[pyfunction]
pub fn length_distribution(path: &str) -> PyResult<HashMap<usize, usize>> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut dist: HashMap<usize, usize> = HashMap::new();
    for rec in &records {
        *dist.entry(rec.seq.len()).or_insert(0) += 1;
    }
    Ok(dist)
}

/// GC含量分布
#[pyfunction]
pub fn gc_distribution(path: &str, bin_size: f64) -> PyResult<HashMap<String, usize>> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut dist: HashMap<String, usize> = HashMap::new();
    for rec in &records {
        let gc = super::sequence::gc_content(&rec.seq);
        let bin = (gc / bin_size).floor() * bin_size;
        let key = format!("{:.1}-{:.1}", bin, bin + bin_size);
        *dist.entry(key).or_insert(0) += 1;
    }
    Ok(dist)
}

/// Sliding window质量过滤
#[pyfunction]
pub fn sliding_window_filter(
    path: &str,
    window_size: usize,
    min_avg_qual: f64,
) -> PyResult<(usize, usize)> {
    let records = parse_fastq_records(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let total = records.len();
    let mut passed = 0usize;

    for rec in &records {
        let scores = qual_scores(&rec.qual);
        let mut ok = true;
        if scores.len() < window_size {
            let mean = scores.iter().map(|&q| q as f64).sum::<f64>() / scores.len() as f64;
            if mean < min_avg_qual {
                ok = false;
            }
        } else {
            for window in scores.windows(window_size) {
                let mean = window.iter().map(|&q| q as f64).sum::<f64>() / window_size as f64;
                if mean < min_avg_qual {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            passed += 1;
        }
    }
    Ok((total, passed))
}
