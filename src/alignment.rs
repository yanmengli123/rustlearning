use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::cmp;

/// Hamming距离（等长序列）
#[pyfunction]
pub fn hamming_distance(a: &str, b: &str) -> PyResult<usize> {
    if a.len() != b.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Sequences must have the same length",
        ));
    }
    Ok(a.chars().zip(b.chars()).filter(|(x, y)| x != y).count())
}

/// Levenshtein编辑距离
#[pyfunction]
pub fn levenshtein_distance(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let n = a.len();
    let m = b.len();

    let mut dp = vec![vec![0usize; m + 1]; n + 1];
    for i in 0..=n {
        dp[i][0] = i;
    }
    for j in 0..=m {
        dp[0][j] = j;
    }
    for i in 1..=n {
        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            dp[i][j] = cmp::min(
                cmp::min(dp[i - 1][j] + 1, dp[i][j - 1] + 1),
                dp[i - 1][j - 1] + cost,
            );
        }
    }
    dp[n][m]
}

/// Needleman-Wunsch全局比对
#[pyfunction]
pub fn needleman_wunsch(
    a: &str,
    b: &str,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> PyResult<(String, String, i32)> {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let n = a.len();
    let m = b.len();

    let mut dp = vec![vec![0i32; m + 1]; n + 1];
    for i in 0..=n {
        dp[i][0] = i as i32 * gap_penalty;
    }
    for j in 0..=m {
        dp[0][j] = j as i32 * gap_penalty;
    }

    for i in 1..=n {
        for j in 1..=m {
            let score = if a[i - 1] == b[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };
            dp[i][j] = cmp::max(
                cmp::max(dp[i - 1][j - 1] + score, dp[i - 1][j] + gap_penalty),
                dp[i][j - 1] + gap_penalty,
            );
        }
    }

    // traceback
    let mut aligned_a = String::new();
    let mut aligned_b = String::new();
    let mut i = n;
    let mut j = m;
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let score = if a[i - 1] == b[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };
            if dp[i][j] == dp[i - 1][j - 1] + score {
                aligned_a.push(a[i - 1]);
                aligned_b.push(b[j - 1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i][j] == dp[i - 1][j] + gap_penalty {
            aligned_a.push(a[i - 1]);
            aligned_b.push('-');
            i -= 1;
        } else {
            aligned_a.push('-');
            aligned_b.push(b[j - 1]);
            j -= 1;
        }
    }

    aligned_a = aligned_a.chars().rev().collect();
    aligned_b = aligned_b.chars().rev().collect();

    Ok((aligned_a, aligned_b, dp[n][m]))
}

/// Smith-Waterman局部比对
#[pyfunction]
pub fn smith_waterman(
    a: &str,
    b: &str,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> PyResult<(String, String, i32)> {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let n = a.len();
    let m = b.len();

    let mut dp = vec![vec![0i32; m + 1]; n + 1];
    let mut max_score = 0i32;
    let mut max_i = 0;
    let mut max_j = 0;

    for i in 1..=n {
        for j in 1..=m {
            let score = if a[i - 1] == b[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };
            dp[i][j] = cmp::max(
                0,
                cmp::max(
                    cmp::max(dp[i - 1][j - 1] + score, dp[i - 1][j] + gap_penalty),
                    dp[i][j - 1] + gap_penalty,
                ),
            );
            if dp[i][j] > max_score {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // traceback
    let mut aligned_a = String::new();
    let mut aligned_b = String::new();
    let mut i = max_i;
    let mut j = max_j;
    while i > 0 && j > 0 && dp[i][j] > 0 {
        let score = if a[i - 1] == b[j - 1] {
            match_score
        } else {
            mismatch_penalty
        };
        if dp[i][j] == dp[i - 1][j - 1] + score {
            aligned_a.push(a[i - 1]);
            aligned_b.push(b[j - 1]);
            i -= 1;
            j -= 1;
        } else if dp[i][j] == dp[i - 1][j] + gap_penalty {
            aligned_a.push(a[i - 1]);
            aligned_b.push('-');
            i -= 1;
        } else {
            aligned_a.push('-');
            aligned_b.push(b[j - 1]);
            j -= 1;
        }
    }

    aligned_a = aligned_a.chars().rev().collect();
    aligned_b = aligned_b.chars().rev().collect();

    Ok((aligned_a, aligned_b, max_score))
}

/// Banded Needleman-Wunsch（带状全局比对）
#[pyfunction]
pub fn banded_alignment(
    a: &str,
    b: &str,
    bandwidth: usize,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> PyResult<i32> {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let n = a.len();
    let m = b.len();

    if n.abs_diff(m) > bandwidth {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Length difference exceeds bandwidth",
        ));
    }

    let mut dp = vec![vec![i32::MIN / 2; m + 1]; n + 1];
    dp[0][0] = 0;

    for i in 0..=n {
        let j_start = if i > bandwidth { i - bandwidth } else { 0 };
        let j_end = cmp::min(m, i + bandwidth);
        for j in j_start..=j_end {
            if i > 0 && j > 0 {
                let score = if a[i - 1] == b[j - 1] {
                    match_score
                } else {
                    mismatch_penalty
                };
                dp[i][j] = dp[i - 1][j - 1] + score;
            }
            if i > 0 {
                dp[i][j] = cmp::max(dp[i][j], dp[i - 1][j] + gap_penalty);
            }
            if j > 0 {
                dp[i][j] = cmp::max(dp[i][j], dp[i][j - 1] + gap_penalty);
            }
        }
    }

    Ok(dp[n][m])
}

/// CIGAR字符串解析
#[pyfunction]
pub fn parse_cigar(cigar: &str) -> Vec<(char, usize)> {
    let mut ops = Vec::new();
    let mut num = 0usize;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + c.to_digit(10).unwrap() as usize;
        } else {
            ops.push((c, num.max(1)));
            num = 0;
        }
    }
    ops
}

/// CIGAR操作统计
#[pyfunction]
pub fn cigar_stats(py: Python, cigar: &str) -> PyResult<PyObject> {
    let ops = parse_cigar(cigar);
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut insertions = 0usize;
    let mut deletions = 0usize;
    let mut soft_clips = 0usize;
    let mut hard_clips = 0usize;

    for (op, len) in &ops {
        match op {
            'M' | '=' => matches += len,
            'X' => mismatches += len,
            'I' => insertions += len,
            'D' => deletions += len,
            'S' => soft_clips += len,
            'H' => hard_clips += len,
            _ => {}
        }
    }

    let dict = PyDict::new_bound(py);
    dict.set_item("matches", matches)?;
    dict.set_item("mismatches", mismatches)?;
    dict.set_item("insertions", insertions)?;
    dict.set_item("deletions", deletions)?;
    dict.set_item("soft_clips", soft_clips)?;
    dict.set_item("hard_clips", hard_clips)?;
    Ok(dict.into())
}

/// 引物序列模糊匹配（允许一定mismatch）
#[pyfunction]
pub fn primer_match(seq: &str, primer: &str, max_mismatch: usize) -> Vec<(usize, usize)> {
    let primer_len = primer.len();
    let mut hits = Vec::new();

    if primer_len > seq.len() {
        return hits;
    }

    for i in 0..=seq.len() - primer_len {
        let subseq = &seq[i..i + primer_len];
        let dist = levenshtein_distance(subseq, primer);
        if dist <= max_mismatch {
            hits.push((i, dist));
        }
    }
    hits
}
