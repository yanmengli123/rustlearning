use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

/// 基因组区间
#[derive(Debug, Clone)]
pub struct GenomicInterval {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
}

/// 区间重叠检测
#[pyfunction]
pub fn intervals_overlap(
    chrom1: &str, start1: i64, end1: i64,
    chrom2: &str, start2: i64, end2: i64,
) -> bool {
    chrom1 == chrom2 && start1 < end2 && start2 < end1
}

/// 区间重叠长度
#[pyfunction]
pub fn overlap_length(
    chrom1: &str, start1: i64, end1: i64,
    chrom2: &str, start2: i64, end2: i64,
) -> i64 {
    if chrom1 != chrom2 {
        return 0;
    }
    let overlap_start = start1.max(start2);
    let overlap_end = end1.min(end2);
    if overlap_start < overlap_end {
        overlap_end - overlap_start
    } else {
        0
    }
}

/// 批量区间overlap查询
#[pyfunction]
pub fn batch_overlap(
    intervals_a: Vec<(String, i64, i64)>,
    intervals_b: Vec<(String, i64, i64)>,
) -> Vec<(usize, usize)> {
    let mut results = Vec::new();
    for (i, (c1, s1, e1)) in intervals_a.iter().enumerate() {
        for (j, (c2, s2, e2)) in intervals_b.iter().enumerate() {
            if c1 == c2 && s1 < e2 && s2 < e1 {
                results.push((i, j));
            }
        }
    }
    results
}

/// Sliding genomic bins
#[pyfunction]
pub fn sliding_bins(
    chrom: &str,
    start: i64,
    end: i64,
    bin_size: i64,
    step: i64,
) -> Vec<(String, i64, i64)> {
    let mut bins = Vec::new();
    let mut pos = start;
    while pos + bin_size <= end {
        bins.push((chrom.to_string(), pos, pos + bin_size));
        pos += step;
    }
    bins
}

/// Window-based counting
#[pyfunction]
pub fn window_count(
    intervals: Vec<(String, i64, i64)>,
    reads: Vec<(String, i64)>,
    window_size: i64,
) -> Vec<i64> {
    let mut counts = Vec::new();
    for (chrom, start, end) in &intervals {
        let n_windows = ((end - start) / window_size) as usize;
        let mut window_counts = vec![0i64; n_windows + 1];
        for (r_chrom, r_pos) in &reads {
            if r_chrom == chrom && r_pos >= start && r_pos < end {
                let idx = ((r_pos - start) / window_size) as usize;
                if idx < window_counts.len() {
                    window_counts[idx] += 1;
                }
            }
        }
        counts.extend(window_counts);
    }
    counts
}

/// TSS距离计算
#[pyfunction]
pub fn tss_distance(
    read_chrom: &str,
    read_pos: i64,
    gene_chrom: &str,
    gene_start: i64,
    gene_end: i64,
    gene_strand: char,
) -> Option<i64> {
    if read_chrom != gene_chrom {
        return None;
    }
    let tss = if gene_strand == '+' { gene_start } else { gene_end };
    Some((read_pos - tss).abs())
}

/// 区间长度
#[pyfunction]
pub fn interval_length(start: i64, end: i64) -> i64 {
    (end - start).abs()
}

/// Interval Tree查询辅助（简化版：线性扫描）
#[pyfunction]
pub fn interval_tree_query(
    query_chrom: &str,
    query_start: i64,
    query_end: i64,
    intervals: Vec<(String, i64, i64)>,
) -> Vec<(usize, String, i64, i64)> {
    let mut results = Vec::new();
    for (i, (chrom, start, end)) in intervals.iter().enumerate() {
        if chrom == query_chrom && *start < query_end && query_start < *end {
            results.push((i, chrom.clone(), *start, *end));
        }
    }
    results
}

/// Blacklist区域过滤
#[pyfunction]
pub fn filter_blacklist(
    intervals: Vec<(String, i64, i64)>,
    blacklist: Vec<(String, i64, i64)>,
) -> Vec<(String, i64, i64)> {
    intervals
        .into_iter()
        .filter(|(chrom, start, end)| {
            !blacklist
                .iter()
                .any(|(bc, bs, be)| chrom == bc && *start < *be && *bs < *end)
        })
        .collect()
}

/// Nearest feature查找
#[pyfunction]
pub fn nearest_feature(
    query_chrom: &str,
    query_pos: i64,
    features: Vec<(String, i64, i64, String)>,
) -> Option<(String, i64, String)> {
    let mut min_dist = i64::MAX;
    let mut nearest = None;

    for (chrom, start, end, name) in &features {
        if chrom != query_chrom {
            continue;
        }
        let dist = if query_pos < *start {
            *start - query_pos
        } else if query_pos > *end {
            query_pos - *end
        } else {
            0
        };
        if dist < min_dist {
            min_dist = dist;
            nearest = Some((name.clone(), min_dist, if query_pos < *start { "upstream".to_string() } else if query_pos > *end { "downstream".to_string() } else { "overlapping".to_string() }));
        }
    }
    nearest
}

/// Enhancer-Promoter匹配
#[pyfunction]
pub fn match_enhancer_promoter(
    enhancers: Vec<(String, i64, i64, String)>,
    promoters: Vec<(String, i64, i64, String)>,
    max_distance: i64,
) -> Vec<(String, String, i64)> {
    let mut matches = Vec::new();
    for (ec, es, ee, eid) in &enhancers {
        for (pc, ps, pe, pid) in &promoters {
            if ec != pc {
                continue;
            }
            let dist = if *ee < *ps {
                *ps - *ee
            } else if *pe < *es {
                *es - *pe
            } else {
                0
            };
            if dist <= max_distance {
                matches.push((eid.clone(), pid.clone(), dist));
            }
        }
    }
    matches
}

/// 覆盖度统计（区间内的reads数）
#[pyfunction]
pub fn coverage_over_intervals(
    intervals: Vec<(String, i64, i64)>,
    read_positions: Vec<(String, i64)>,
) -> Vec<i64> {
    intervals
        .iter()
        .map(|(chrom, start, end)| {
            read_positions
                .iter()
                .filter(|(rc, rp)| rc == chrom && rp >= start && rp < end)
                .count() as i64
        })
        .collect()
}
