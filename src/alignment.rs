//! ============================================================
//! 模块4: 序列比对 (alignment)
//! ============================================================
//! 本模块提供序列比对相关的算法实现。
//! 包括：Hamming距离、Levenshtein编辑距离、
//! Needleman-Wunsch全局比对、Smith-Waterman局部比对、
//! 带状比对、CIGAR字符串解析、引物匹配等。
//!
//! 序列比对是生物信息学的核心算法：
//! - 全局比对（Needleman-Wunsch）：比较两条序列的整体相似性
//! - 局部比对（Smith-Waterman）：寻找最相似的局部区域
//! - 编辑距离：衡量两个字符串的差异程度
//!
//! CIGAR字符串（Compact Idiosyncratic Gapped Alignment Report）：
//!   SAM格式中描述比对结果的标准表示法
//!   M=匹配/错配, I=插入, D=缺失, S=软裁剪, H=硬裁剪
//!
//! 设计原则：
//! - 使用动态规划算法，保证最优解
//! - 返回Python兼容类型（元组/字典）
//! - 参数验证，错误输入返回PyErr
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use std::cmp;                   // 比较函数

/// -----------------------------------------------------------
/// Hamming距离（等长序列比较）
/// -----------------------------------------------------------
/// 参数:
///   a - 序列1
///   b - 序列2（必须与a等长）
/// 返回: 不同字符的数量
/// 错误: 如果序列长度不等，抛出ValueError
/// 算法:
///   逐位比较两个等长序列，统计不同位置的数量
///   Hamming距离只能用于等长序列，不支持插入/缺失
/// 用途: 快速比较相似序列（如同一基因的不同变异）
#[pyfunction]
pub fn hamming_distance(a: &str, b: &str) -> PyResult<usize> {
    if a.len() != b.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Sequences must have the same length",
        ));
    }
    // 逐字符比较，统计不同的数量
    Ok(a.chars().zip(b.chars()).filter(|(x, y)| x != y).count())
}

/// -----------------------------------------------------------
/// Levenshtein编辑距离
/// -----------------------------------------------------------
/// 参数:
///   a - 序列1
///   b - 序列2
/// 返回: 最小编辑距离（插入、删除、替换操作数）
/// 算法: 动态规划
///   dp[i][j] = 将a[0..i]转换为b[0..j]所需的最小操作数
///   转移方程:
///     dp[i][j] = min(
///       dp[i-1][j] + 1,    // 删除a[i]
///       dp[i][j-1] + 1,    // 插入b[j]
///       dp[i-1][j-1] + cost // 替换或匹配
///     )
///   时间复杂度: O(n*m)
/// 用途: 模糊字符串匹配、拼写纠错、DNA序列比较
#[pyfunction]
pub fn levenshtein_distance(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let n = a.len();
    let m = b.len();

    // 初始化DP矩阵
    let mut dp = vec![vec![0usize; m + 1]; n + 1];
    // 边界条件：空字符串转换为另一个字符串
    for i in 0..=n {
        dp[i][0] = i;  // 删除i个字符
    }
    for j in 0..=m {
        dp[0][j] = j;  // 插入j个字符
    }
    // 填充DP矩阵
    for i in 1..=n {
        for j in 1..=m {
            // 匹配代价：相同为0，不同为1
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            dp[i][j] = cmp::min(
                cmp::min(dp[i - 1][j] + 1, dp[i][j - 1] + 1),  // 删除或插入
                dp[i - 1][j - 1] + cost,  // 替换或匹配
            );
        }
    }
    dp[n][m]
}

/// -----------------------------------------------------------
/// Needleman-Wunsch全局比对
/// -----------------------------------------------------------
/// 参数:
///   a               - 序列1
///   b               - 序列2
///   match_score     - 匹配得分（正数）
///   mismatch_penalty - 错配罚分（负数）
///   gap_penalty     - 空位罚分（负数）
/// 返回: (比对后序列1, 比对后序列2, 总得分)
/// 算法: 动态规划 + 回溯
///   1. 初始化DP矩阵，第一行/列为gap罚分累加
///   2. 填充矩阵：dp[i][j] = max(
///        dp[i-1][j-1] + match/mismatch_score,
///        dp[i-1][j] + gap_penalty,
///        dp[i][j-1] + gap_penalty
///      )
///   3. 从右下角回溯到左上角，构建比对结果
/// 用途: 同源基因比对、蛋白质序列比较
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

    // 初始化DP矩阵
    let mut dp = vec![vec![0i32; m + 1]; n + 1];
    for i in 0..=n {
        dp[i][0] = i as i32 * gap_penalty;  // 第一列
    }
    for j in 0..=m {
        dp[0][j] = j as i32 * gap_penalty;  // 第一行
    }

    // 填充DP矩阵
    for i in 1..=n {
        for j in 1..=m {
            // 匹配或错配得分
            let score = if a[i - 1] == b[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };
            // 三种转移：对角线（匹配/错配）、上方（删除）、左方（插入）
            dp[i][j] = cmp::max(
                cmp::max(dp[i - 1][j - 1] + score, dp[i - 1][j] + gap_penalty),
                dp[i][j - 1] + gap_penalty,
            );
        }
    }

    // 回溯构建比对结果
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
            // 对角线转移：匹配或错配
            if dp[i][j] == dp[i - 1][j - 1] + score {
                aligned_a.push(a[i - 1]);
                aligned_b.push(b[j - 1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        // 上方转移：删除（a中有gap）
        if i > 0 && dp[i][j] == dp[i - 1][j] + gap_penalty {
            aligned_a.push(a[i - 1]);
            aligned_b.push('-');  // gap
            i -= 1;
        } else {
            // 左方转移：插入（b中有gap）
            aligned_a.push('-');  // gap
            aligned_b.push(b[j - 1]);
            j -= 1;
        }
    }

    // 反转得到正向比对
    aligned_a = aligned_a.chars().rev().collect();
    aligned_b = aligned_b.chars().rev().collect();

    Ok((aligned_a, aligned_b, dp[n][m]))
}

/// -----------------------------------------------------------
/// Smith-Waterman局部比对
/// -----------------------------------------------------------
/// 参数:
///   a               - 序列1
///   b               - 序列2
///   match_score     - 匹配得分（正数）
///   mismatch_penalty - 错配罚分（负数）
///   gap_penalty     - 空位罚分（负数）
/// 返回: (比对后序列1, 比对后序列2, 最佳局部得分)
/// 算法: 动态规划 + 回溯
///   与NW全局比对的主要区别：
///   1. 矩阵初始化为0（而非gap累加）
///   2. 转移时取max(0, ...)，允许从任意位置开始
///   3. 回溯从最大得分位置开始，到得分为0结束
/// 用途: 寻找序列中的保守结构域、局部相似性搜索
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

    // 初始化DP矩阵（全0）
    let mut dp = vec![vec![0i32; m + 1]; n + 1];
    let mut max_score = 0i32;  // 记录最大得分
    let mut max_i = 0;         // 最大得分位置
    let mut max_j = 0;

    // 填充DP矩阵
    for i in 1..=n {
        for j in 1..=m {
            let score = if a[i - 1] == b[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };
            // 局部比对：允许从0开始（丢弃负分区域）
            dp[i][j] = cmp::max(
                0,
                cmp::max(
                    cmp::max(dp[i - 1][j - 1] + score, dp[i - 1][j] + gap_penalty),
                    dp[i][j - 1] + gap_penalty,
                ),
            );
            // 记录最大得分位置
            if dp[i][j] > max_score {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // 从最大得分位置回溯
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

/// -----------------------------------------------------------
/// Banded Needleman-Wunsch（带状全局比对）
/// -----------------------------------------------------------
/// 参数:
///   a               - 序列1
///   b               - 序列2
///   bandwidth       - 带宽（限制对角线附近的搜索范围）
///   match_score     - 匹配得分
///   mismatch_penalty - 错配罚分
///   gap_penalty     - 空位罚分
/// 返回: 最佳比对得分
/// 算法:
///   只在主对角线附近bandwidth范围内计算DP
///   大幅减少计算量：O(n*m) -> O(n*bandwidth)
///   要求：|len(a) - len(b)| <= bandwidth
/// 用途: 长序列快速比对（如基因组mapping）
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

    // 检查长度差异是否在带宽范围内
    if n.abs_diff(m) > bandwidth {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Length difference exceeds bandwidth",
        ));
    }

    // 初始化DP矩阵，带外区域用极小值填充
    let mut dp = vec![vec![i32::MIN / 2; m + 1]; n + 1];
    dp[0][0] = 0;

    // 只在带内区域计算
    for i in 0..=n {
        // 计算当前行的带宽范围
        let j_start = if i > bandwidth { i - bandwidth } else { 0 };
        let j_end = cmp::min(m, i + bandwidth);
        for j in j_start..=j_end {
            if i > 0 && j > 0 {
                let score = if a[i - 1] == b[j - 1] {
                    match_score
                } else {
                    mismatch_penalty
                };
                dp[i][j] = dp[i - 1][j - 1] + score;  // 对角线
            }
            if i > 0 {
                dp[i][j] = cmp::max(dp[i][j], dp[i - 1][j] + gap_penalty);  // 上方
            }
            if j > 0 {
                dp[i][j] = cmp::max(dp[i][j], dp[i][j - 1] + gap_penalty);  // 左方
            }
        }
    }

    Ok(dp[n][m])
}

/// -----------------------------------------------------------
/// CIGAR字符串解析
/// -----------------------------------------------------------
/// 参数: cigar - CIGAR字符串（如"10M2I5M3D"）
/// 返回: Vec<(char, usize)> - (操作符, 长度)列表
/// CIGAR操作符说明:
///   M = 匹配或错配（0表示比对）
///   I = 插入（参考序列中有gap）
///   D = 删除（查询序列中有gap）
///   S = 软裁剪（序列被裁剪但保留在比对中）
///   H = 硬裁剪（序列被裁剪且不保留在比对中）
///   = = 精确匹配
///   X = 错配
/// 用途: SAM/BAM文件解析、比对结果分析
#[pyfunction]
pub fn parse_cigar(cigar: &str) -> Vec<(char, usize)> {
    let mut ops = Vec::new();
    let mut num = 0usize;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            // 数字字符：累积长度
            num = num * 10 + c.to_digit(10).unwrap() as usize;
        } else {
            // 操作符：记录操作和长度
            ops.push((c, num.max(1)));  // 默认长度为1
            num = 0;  // 重置数字
        }
    }
    ops
}

/// -----------------------------------------------------------
/// CIGAR操作统计
/// -----------------------------------------------------------
/// 参数:
///   py    - Python解释器引用
///   cigar - CIGAR字符串
/// 返回: Python字典，包含各种操作的统计
/// 统计指标:
///   - matches: 匹配/错配碱基数
///   - mismatches: 错配碱基数
///   - insertions: 插入碱基数
///   - deletions: 删除碱基数
///   - soft_clips: 软裁剪碱基数
///   - hard_clips: 硬裁剪碱基数
/// 用途: 评估比对质量、计算比对覆盖率
#[pyfunction]
pub fn cigar_stats(py: Python, cigar: &str) -> PyResult<PyObject> {
    let ops = parse_cigar(cigar);
    let mut matches = 0usize;      // 匹配/错配
    let mut mismatches = 0usize;   // 错配
    let mut insertions = 0usize;   // 插入
    let mut deletions = 0usize;    // 删除
    let mut soft_clips = 0usize;   // 软裁剪
    let mut hard_clips = 0usize;   // 硬裁剪

    for (op, len) in &ops {
        match op {
            'M' | '=' => matches += len,   // 匹配（M包含匹配和错配）
            'X' => mismatches += len,      // 错配
            'I' => insertions += len,      // 插入
            'D' => deletions += len,       // 删除
            'S' => soft_clips += len,      // 软裁剪
            'H' => hard_clips += len,      // 硬裁剪
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

/// -----------------------------------------------------------
/// 引物序列模糊匹配
/// -----------------------------------------------------------
/// 参数:
///   seq          - 目标序列
///   primer       - 引物序列
///   max_mismatch - 最大允许错配数
/// 返回: Vec<(usize, usize)> - (位置, 错配数)列表
/// 算法:
///   1. 在目标序列上滑动窗口（窗口大小=引物长度）
///   2. 计算每个位置与引物的编辑距离
///   3. 返回编辑距离≤max_mismatch的所有位置
/// 用途: PCR引物设计、序列定位
#[pyfunction]
pub fn primer_match(seq: &str, primer: &str, max_mismatch: usize) -> Vec<(usize, usize)> {
    let primer_len = primer.len();
    let mut hits = Vec::new();

    if primer_len > seq.len() {
        return hits;  // 引物比序列长，无法匹配
    }

    // 滑动窗口扫描
    for i in 0..=seq.len() - primer_len {
        let subseq = &seq[i..i + primer_len];
        let dist = levenshtein_distance(subseq, primer);  // 计算编辑距离
        if dist <= max_mismatch {
            hits.push((i, dist));  // 记录匹配位置和错配数
        }
    }
    hits
}
