use pyo3::prelude::*;
use pyo3::types::PyDict;
use rand::Rng;

/// Permutation test
#[pyfunction]
pub fn permutation_test(
    py: Python,
    group_a: Vec<f64>,
    group_b: Vec<f64>,
    n_permutations: usize,
) -> PyResult<PyObject> {
    let observed_diff = mean(&group_a) - mean(&group_b);
    let mut rng = rand::thread_rng();
    let mut combined = group_a.clone();
    combined.extend(group_b.iter());
    let n_a = group_a.len();
    let mut extreme_count = 0usize;

    for _ in 0..n_permutations {
        // Fisher-Yates shuffle
        for i in (1..combined.len()).rev() {
            let j = rng.gen_range(0..=i);
            combined.swap(i, j);
        }
        let perm_diff = mean(&combined[..n_a]) - mean(&combined[n_a..]);
        if perm_diff.abs() >= observed_diff.abs() {
            extreme_count += 1;
        }
    }

    let p_value = extreme_count as f64 / n_permutations as f64;

    let dict = PyDict::new_bound(py);
    dict.set_item("observed_diff", observed_diff)?;
    dict.set_item("p_value", p_value)?;
    dict.set_item("n_permutations", n_permutations)?;
    Ok(dict.into())
}

fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    data.iter().sum::<f64>() / data.len() as f64
}

/// Bootstrap置信区间
#[pyfunction]
pub fn bootstrap_ci(
    py: Python,
    data: Vec<f64>,
    n_bootstrap: usize,
    confidence: f64,
) -> PyResult<PyObject> {
    let mut rng = rand::thread_rng();
    let n = data.len();
    let mut means = Vec::with_capacity(n_bootstrap);

    for _ in 0..n_bootstrap {
        let mut sum = 0.0f64;
        for _ in 0..n {
            let idx = rng.gen_range(0..n);
            sum += data[idx];
        }
        means.push(sum / n as f64);
    }

    means.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let alpha = (1.0 - confidence) / 2.0;
    let lower_idx = (alpha * n_bootstrap as f64) as usize;
    let upper_idx = ((1.0 - alpha) * n_bootstrap as f64) as usize;

    let dict = PyDict::new_bound(py);
    dict.set_item("mean", mean(&data))?;
    dict.set_item("ci_lower", means[lower_idx.min(means.len() - 1)])?;
    dict.set_item("ci_upper", means[upper_idx.min(means.len() - 1)])?;
    dict.set_item("confidence", confidence)?;
    dict.set_item("n_bootstrap", n_bootstrap)?;
    Ok(dict.into())
}

/// Hypergeometric test（富集分析用）
#[pyfunction]
pub fn hypergeometric_test(
    k: u64,     // 样本中属于该类别的数量
    n: u64,     // 样本大小
    k_total: u64, // 总体中属于该类别的数量
    n_total: u64, // 总体大小
) -> f64 {
    // P(X >= k) 的生存函数近似
    let p = k_total as f64 / n_total as f64;
    let expected = n as f64 * p;
    if k as f64 <= expected {
        return 1.0;
    }
    // 使用正态近似
    let variance = n as f64 * p * (1.0 - p) * (n_total as f64 - n as f64) / (n_total as f64 - 1.0);
    let z = (k as f64 - expected) / variance.sqrt();
    // 简化的p值近似
    (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// Fisher exact test（2x2列联表）
#[pyfunction]
pub fn fisher_exact(a: u64, b: u64, c: u64, d: u64) -> f64 {
    // 使用log-odds的正态近似
    let a = a as f64;
    let b = b as f64;
    let c = c as f64;
    let d = d as f64;
    let n = a + b + c + d;

    let odds_ratio = (a * d) / (b * c).max(1e-300);
    let se = (1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d).sqrt();
    let z = odds_ratio.ln() / se;
    // 简化的p值
    (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// 多重检验校正（Benjamini-Hochberg FDR）
#[pyfunction]
pub fn bh_correction(p_values: Vec<f64>) -> Vec<f64> {
    let n = p_values.len();
    if n == 0 {
        return Vec::new();
    }

    let mut indexed: Vec<(usize, f64)> = p_values.iter().enumerate().map(|(i, &p)| (i, p)).collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut adjusted = vec![0.0f64; n];
    let mut min_adj = 1.0f64;

    for (rank, (idx, p)) in indexed.iter().enumerate() {
        let adj = p * n as f64 / (rank + 1) as f64;
        let adj = adj.min(min_adj).min(1.0);
        adjusted[*idx] = adj;
        min_adj = min_adj.min(adj);
    }
    adjusted
}

/// Rolling correlation
#[pyfunction]
pub fn rolling_correlation(x: Vec<f64>, y: Vec<f64>, window: usize) -> Vec<f64> {
    let mut correlations = Vec::new();
    for i in 0..=x.len().saturating_sub(window) {
        let x_win = &x[i..i + window];
        let y_win = &y[i..i + window];
        correlations.push(pearson_corr(x_win, y_win));
    }
    correlations
}

/// Pearson相关系数
#[pyfunction]
pub fn pearson_correlation(x: Vec<f64>, y: Vec<f64>) -> f64 {
    pearson_corr(&x, &y)
}

pub fn pearson_corr(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return 0.0;
    }
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;

    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    let denom = (var_x * var_y).sqrt();
    if denom == 0.0 {
        return 0.0;
    }
    cov / denom
}

/// Mann-Whitney U检验（Wilcoxon rank-sum）
#[pyfunction]
pub fn mann_whitney_u(py: Python, group_a: Vec<f64>, group_b: Vec<f64>) -> PyResult<PyObject> {
    let n_a = group_a.len() as f64;
    let n_b = group_b.len() as f64;
    let n = n_a + n_b;

    let mut all: Vec<(f64, bool)> = Vec::new();
    for v in &group_a {
        all.push((*v, true));
    }
    for v in &group_b {
        all.push((*v, false));
    }
    all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut rank_sum_a = 0.0f64;
    for (i, (_, is_a)) in all.iter().enumerate() {
        if *is_a {
            rank_sum_a += (i + 1) as f64;
        }
    }

    let u_a = rank_sum_a - n_a * (n_a + 1.0) / 2.0;
    let u_b = n_a * n_b - u_a;
    let u = u_a.min(u_b);

    // 正态近似
    let mean_u = n_a * n_b / 2.0;
    let sd_u = (n_a * n_b * (n + 1.0) / 12.0).sqrt();
    let z = if sd_u > 0.0 { (u - mean_u) / sd_u } else { 0.0 };
    let p_value = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();

    let dict = PyDict::new_bound(py);
    dict.set_item("u_statistic", u)?;
    dict.set_item("p_value", p_value)?;
    Ok(dict.into())
}

/// 基本描述统计
#[pyfunction]
pub fn descriptive_stats(py: Python, data: Vec<f64>) -> PyResult<PyObject> {
    if data.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err("Empty data"));
    }
    let n = data.len() as f64;
    let sum: f64 = data.iter().sum();
    let mean = sum / n;

    let mut sorted = data.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let median = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    let variance = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    let std_dev = variance.sqrt();
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];
    let q1 = sorted[sorted.len() / 4];
    let q3 = sorted[sorted.len() * 3 / 4];

    let dict = PyDict::new_bound(py);
    dict.set_item("count", data.len())?;
    dict.set_item("mean", mean)?;
    dict.set_item("median", median)?;
    dict.set_item("std_dev", std_dev)?;
    dict.set_item("variance", variance)?;
    dict.set_item("min", min)?;
    dict.set_item("max", max)?;
    dict.set_item("q1", q1)?;
    dict.set_item("q3", q3)?;
    Ok(dict.into())
}
