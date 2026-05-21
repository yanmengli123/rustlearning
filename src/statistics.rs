//! ============================================================
//! 模块13: 统计分析 (statistics)
//! ============================================================
//! 本模块提供生物信息学常用的统计分析功能。
//! 包括：置换检验、Bootstrap置信区间、超几何检验、
//! Fisher精确检验、BH多重检验校正、Pearson相关系数、
//! 滚动相关、Mann-Whitney U检验、描述统计等。
//!
//! 统计方法在生物信息学中的应用：
//! - 置换检验：差异表达分析、富集分析
//! - Bootstrap：估计统计量的不确定性
//! - 多重检验校正：控制假发现率（FDR）
//! - 相关分析：基因共表达、表型关联
//! - 非参数检验：不满足正态假设的数据
//!
//! 设计原则：
//! - 使用正态近似（避免复杂的精确计算）
//! - 返回Python字典便于结果展示
//! - 处理边界情况（空数据、零方差等）
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use rand::Rng;                  // 随机数生成

/// -----------------------------------------------------------
/// Permutation test（置换检验）
/// -----------------------------------------------------------
/// 参数:
///   py            - Python解释器引用
///   group_a       - 组A数据
///   group_b       - 组B数据
///   n_permutations - 置换次数
/// 返回: Python字典，包含观察差值和p值
/// 算法:
///   1. 计算观察到的组间差值 (mean_A - mean_B)
///   2. 合并两组数据
///   3. 随机打乱n_permutations次
///   4. 每次计算随机分组的差值
///   5. p值 = 极端情况数 / 总置换次数
/// 用途: 差异表达分析、假设检验
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
            extreme_count += 1;  // 计数极端情况
        }
    }

    let p_value = extreme_count as f64 / n_permutations as f64;

    let dict = PyDict::new_bound(py);
    dict.set_item("observed_diff", observed_diff)?;
    dict.set_item("p_value", p_value)?;
    dict.set_item("n_permutations", n_permutations)?;
    Ok(dict.into())
}

/// 辅助函数：计算均值
fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    data.iter().sum::<f64>() / data.len() as f64
}

/// -----------------------------------------------------------
/// Bootstrap置信区间
/// -----------------------------------------------------------
/// 参数:
///   py          - Python解释器引用
///   data        - 原始数据
///   n_bootstrap - Bootstrap重采样次数
///   confidence  - 置信水平（如0.95表示95%）
/// 返回: Python字典，包含均值和置信区间
/// 算法:
///   1. 从原始数据中有放回地抽样n_bootstrap次
///   2. 计算每次抽样的均值
///   3. 对均值排序
///   4. 取alpha/2和1-alpha/2分位数作为置信区间
/// 用途: 估计统计量的不确定性
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
            sum += data[idx];  // 有放回抽样
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

/// -----------------------------------------------------------
/// Hypergeometric test（超几何检验）
/// -----------------------------------------------------------
/// 参数:
///   k       - 样本中属于该类别的数量（观测值）
///   n       - 样本大小
///   k_total - 总体中属于该类别的数量
///   n_total - 总体大小
/// 返回: f64 - p值（单尾）
/// 算法: 使用正态近似
///   期望值 E = n * (k_total / n_total)
///   如果观测值k小于期望值，返回1.0
///   否则计算Z分数并转换为p值
/// 用途: 富集分析、GO分析
#[pyfunction]
pub fn hypergeometric_test(
    k: u64,         // 样本中属于该类别的数量
    n: u64,         // 样本大小
    k_total: u64,   // 总体中属于该类别的数量
    n_total: u64,   // 总体大小
) -> f64 {
    // P(X >= k) 的生存函数近似
    let p = k_total as f64 / n_total as f64;
    let expected = n as f64 * p;
    if k as f64 <= expected {
        return 1.0;  // 观测值不显著
    }
    // 使用正态近似
    let variance = n as f64 * p * (1.0 - p) * (n_total as f64 - n as f64) / (n_total as f64 - 1.0);
    let z = (k as f64 - expected) / variance.sqrt();
    // 简化的p值近似
    (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// -----------------------------------------------------------
/// Fisher exact test（Fisher精确检验）
/// -----------------------------------------------------------
/// 参数:
///   a, b, c, d - 2x2列联表的四个单元格
/// 返回: f64 - p值
/// 算法: 使用log-odds的正态近似
///   OR = (a*d) / (b*c)
///   SE = sqrt(1/a + 1/b + 1/c + 1/d)
///   Z = ln(OR) / SE
/// 用途: 小样本的2x2列联表检验
#[pyfunction]
pub fn fisher_exact(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let a = a as f64;
    let b = b as f64;
    let c = c as f64;
    let d = d as f64;
    let odds_ratio = (a * d) / (b * c).max(1e-300);  // 避免除零
    let se = (1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d).sqrt();
    let z = odds_ratio.ln() / se;
    // 简化的p值
    (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// -----------------------------------------------------------
/// 多重检验校正（Benjamini-Hochberg FDR）
/// -----------------------------------------------------------
/// 参数: p_values - 原始p值列表
/// 返回: Vec<f64> - 调整后的p值（q值）
/// 算法:
///   1. 将p值排序（保持原始索引）
///   2. 计算调整后的p值: adj_p = p * n / rank
///   3. 确保单调递减：adj_p[i] = min(adj_p[i], adj_p[i+1])
///   4. 限制在[0,1]范围内
/// 用途: 控制假发现率（FDR），差异表达分析
#[pyfunction]
pub fn bh_correction(p_values: Vec<f64>) -> Vec<f64> {
    let n = p_values.len();
    if n == 0 {
        return Vec::new();
    }

    // 创建索引-p值对并排序
    let mut indexed: Vec<(usize, f64)> = p_values.iter().enumerate().map(|(i, &p)| (i, p)).collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut adjusted = vec![0.0f64; n];
    let mut min_adj = 1.0f64;

    // BH校正
    for (rank, (idx, p)) in indexed.iter().enumerate() {
        let adj = p * n as f64 / (rank + 1) as f64;
        let adj = adj.min(min_adj).min(1.0);
        adjusted[*idx] = adj;
        min_adj = min_adj.min(adj);  // 确保单调性
    }
    adjusted
}

/// -----------------------------------------------------------
/// Rolling correlation（滚动相关）
/// -----------------------------------------------------------
/// 参数:
///   x      - 序列x
///   y      - 序列y
///   window - 窗口大小
/// 返回: Vec<f64> - 每个窗口位置的相关系数
/// 算法:
///   在每个窗口位置计算Pearson相关系数
///   返回 (n - window + 1) 个相关系数
/// 用途: 时序数据分析、基因共表达动态
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

/// -----------------------------------------------------------
/// Pearson相关系数（公共接口）
/// -----------------------------------------------------------
/// 参数:
///   x - 序列x
///   y - 序列y
/// 返回: f64 - 相关系数 (-1.0 ~ 1.0)
/// 用途: 衡量两个变量的线性关系强度
#[pyfunction]
pub fn pearson_correlation(x: Vec<f64>, y: Vec<f64>) -> f64 {
    pearson_corr(&x, &y)
}

/// -----------------------------------------------------------
/// Pearson相关系数（内部实现）
/// -----------------------------------------------------------
/// 参数:
///   x - 序列x
///   y - 序列y
/// 返回: f64 - 相关系数
/// 算法:
///   r = Σ((x_i - x_mean)(y_i - y_mean)) / sqrt(Σ(x_i - x_mean)² * Σ(y_i - y_mean)²)
///   r ∈ [-1, 1]
///   r=1: 完全正相关
///   r=-1: 完全负相关
///   r=0: 无线性相关
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
        cov += dx * dy;      // 协方差
        var_x += dx * dx;    // x方差
        var_y += dy * dy;    // y方差
    }

    let denom = (var_x * var_y).sqrt();
    if denom == 0.0 {
        return 0.0;  // 避免除零
    }
    cov / denom
}

/// -----------------------------------------------------------
/// Mann-Whitney U检验（Wilcoxon rank-sum）
/// -----------------------------------------------------------
/// 参数:
///   py      - Python解释器引用
///   group_a - 组A数据
///   group_b - 组B数据
/// 返回: Python字典，包含U统计量和p值
/// 算法:
///   1. 合并两组数据并排序
///   2. 计算组A的秩和
///   3. 计算U统计量: U_A = R_A - n_A(n_A+1)/2
///   4. 使用正态近似计算p值
/// 用途: 非参数两样本检验（不假设正态分布）
#[pyfunction]
pub fn mann_whitney_u(py: Python, group_a: Vec<f64>, group_b: Vec<f64>) -> PyResult<PyObject> {
    let n_a = group_a.len() as f64;
    let n_b = group_b.len() as f64;
    let n = n_a + n_b;

    // 合并并标记来源
    let mut all: Vec<(f64, bool)> = Vec::new();
    for v in &group_a {
        all.push((*v, true));
    }
    for v in &group_b {
        all.push((*v, false));
    }
    all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // 计算秩和
    let mut rank_sum_a = 0.0f64;
    for (i, (_, is_a)) in all.iter().enumerate() {
        if *is_a {
            rank_sum_a += (i + 1) as f64;
        }
    }

    // 计算U统计量
    let u_a = rank_sum_a - n_a * (n_a + 1.0) / 2.0;
    let u_b = n_a * n_b - u_a;
    let u = u_a.min(u_b);

    // 正态近似计算p值
    let mean_u = n_a * n_b / 2.0;
    let sd_u = (n_a * n_b * (n + 1.0) / 12.0).sqrt();
    let z = if sd_u > 0.0 { (u - mean_u) / sd_u } else { 0.0 };
    let p_value = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();

    let dict = PyDict::new_bound(py);
    dict.set_item("u_statistic", u)?;
    dict.set_item("p_value", p_value)?;
    Ok(dict.into())
}

/// -----------------------------------------------------------
/// 描述统计
/// -----------------------------------------------------------
/// 参数:
///   py   - Python解释器引用
///   data - 数据列表
/// 返回: Python字典，包含各种描述统计量
/// 统计指标:
///   - count: 样本数
///   - mean: 均值
///   - median: 中位数
///   - std_dev: 标准差
///   - variance: 方差
///   - min / max: 最小/最大值
///   - q1 / q3: 第一/第三四分位数
/// 用途: 数据探索、质量评估
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

    // 中位数
    let median = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    // 方差和标准差
    let variance = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    let std_dev = variance.sqrt();

    // 最小/最大值
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];

    // 四分位数
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
