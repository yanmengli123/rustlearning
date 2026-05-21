//! ============================================================
//! 模块7: 基因组区间高级操作 (genomic_intervals)
//! ============================================================
//! 本模块提供基因组区间的高级操作功能。
//! 包括：区间重叠检测、滑动窗口计数、TSS距离计算、
//! 区间树查询、Blacklist过滤、增强子-启动子匹配等。
//!
//! 核心概念：
//! - 基因组区间: 由染色体名、起始位置、结束位置定义
//! - 重叠: 两个区间在同一条染色体上有共同区域
//! - TSS: 转录起始位点（Transcription Start Site）
//! - Blacklist: 已知的不可靠基因组区域（如重复序列、着丝粒）
//!
//! 坐标系说明：
//! - 本模块使用0-based坐标系（start包含，end不包含）
//! - 与BED格式一致
//!
//! 设计原则：
//! - 简化版实现，适合中小规模数据
//! - 使用线性扫描代替区间树（简单但O(n)查询）
//! - 返回Python兼容类型
//! ============================================================

use pyo3::prelude::*;  // PyO3 核心宏和类型

/// -----------------------------------------------------------
/// 基因组区间结构体
/// -----------------------------------------------------------
/// 存储基因组区间的基本信息
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GenomicInterval {
    pub chrom: String,  // 染色体名称
    pub start: i64,     // 起始位置（0-based）
    pub end: i64,       // 结束位置（不包含）
}

/// -----------------------------------------------------------
/// 区间重叠检测
/// -----------------------------------------------------------
/// 参数:
///   chrom1, start1, end1 - 区间1
///   chrom2, start2, end2 - 区间2
/// 返回: bool - 是否重叠
/// 算法:
///   两个区间重叠的条件：
///   1. 在同一条染色体上
///   2. start1 < end2 且 start2 < end1
/// 用途: 快速判断两个区间是否有交集
#[pyfunction]
pub fn intervals_overlap(
    chrom1: &str, start1: i64, end1: i64,
    chrom2: &str, start2: i64, end2: i64,
) -> bool {
    chrom1 == chrom2 && start1 < end2 && start2 < end1
}

/// -----------------------------------------------------------
/// 区间重叠长度
/// -----------------------------------------------------------
/// 参数:
///   chrom1, start1, end1 - 区间1
///   chrom2, start2, end2 - 区间2
/// 返回: i64 - 重叠长度（不重叠返回0）
/// 算法:
///   overlap_start = max(start1, start2)
///   overlap_end = min(end1, end2)
///   重叠长度 = overlap_end - overlap_start（如果>0）
#[pyfunction]
pub fn overlap_length(
    chrom1: &str, start1: i64, end1: i64,
    chrom2: &str, start2: i64, end2: i64,
) -> i64 {
    if chrom1 != chrom2 {
        return 0;  // 不同染色体，无重叠
    }
    let overlap_start = start1.max(start2);
    let overlap_end = end1.min(end2);
    if overlap_start < overlap_end {
        overlap_end - overlap_start
    } else {
        0
    }
}

/// -----------------------------------------------------------
/// 批量区间重叠查询
/// -----------------------------------------------------------
/// 参数:
///   intervals_a - 区间集合A (chrom, start, end)
///   intervals_b - 区间集合B (chrom, start, end)
/// 返回: Vec<(usize, usize)> - 重叠的区间对索引
/// 算法: 暴力搜索，O(n*m)时间复杂度
/// 用途: 找出两个区间集合之间所有的重叠关系
#[pyfunction]
pub fn batch_overlap(
    intervals_a: Vec<(String, i64, i64)>,
    intervals_b: Vec<(String, i64, i64)>,
) -> Vec<(usize, usize)> {
    let mut results = Vec::new();
    for (i, (c1, s1, e1)) in intervals_a.iter().enumerate() {
        for (j, (c2, s2, e2)) in intervals_b.iter().enumerate() {
            if c1 == c2 && s1 < e2 && s2 < e1 {
                results.push((i, j));  // 记录重叠对
            }
        }
    }
    results
}

/// -----------------------------------------------------------
/// 滑动基因组分箱
/// -----------------------------------------------------------
/// 参数:
///   chrom    - 染色体名称
///   start    - 起始位置
///   end      - 结束位置
///   bin_size - 箱大小
///   step     - 步长（可以与bin_size不同，形成重叠窗口）
/// 返回: Vec<(String, i64, i64)> - 分箱区间列表
/// 算法:
///   从start开始，每次移动step，生成长度为bin_size的区间
///   当 pos + bin_size > end 时停止
/// 用途: 基因组分bin统计、GC含量分布、覆盖度分析
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
        pos += step;  // 移动步长
    }
    bins
}

/// -----------------------------------------------------------
/// 基于窗口的计数
/// -----------------------------------------------------------
/// 参数:
///   intervals   - 目标区间列表
///   reads       - reads位置列表 (chrom, pos)
///   window_size - 窗口大小
/// 返回: Vec<i64> - 每个窗口内的reads数
/// 算法:
///   1. 将每个区间划分为多个窗口
///   2. 对每个read，找到它所属的窗口
///   3. 累加窗口计数
/// 用途: 生成覆盖度分布图、peak calling
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
                    window_counts[idx] += 1;  // 对应窗口计数+1
                }
            }
        }
        counts.extend(window_counts);
    }
    counts
}

/// -----------------------------------------------------------
/// TSS距离计算
/// -----------------------------------------------------------
/// 参数:
///   read_chrom  - read染色体
///   read_pos    - read位置
///   gene_chrom  - 基因染色体
///   gene_start  - 基因起始位置
///   gene_end    - 基因结束位置
///   gene_strand - 基因链方向（+/ -）
/// 返回: Option<i64> - 到TSS的距离（不同染色体返回None）
/// 算法:
///   TSS（转录起始位点）的确定：
///   - 正链(+)基因：TSS = gene_start
///   - 负链(-)基因：TSS = gene_end
///   距离 = |read_pos - TSS|
/// 用途: 研究read相对于基因起始位点的分布
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
        return None;  // 不同染色体
    }
    // 根据链方向确定TSS位置
    let tss = if gene_strand == '+' { gene_start } else { gene_end };
    Some((read_pos - tss).abs())  // 返回绝对距离
}

/// -----------------------------------------------------------
/// 区间长度计算
/// -----------------------------------------------------------
/// 参数:
///   start - 起始位置
///   end   - 结束位置
/// 返回: i64 - 区间长度（绝对值）
#[pyfunction]
pub fn interval_length(start: i64, end: i64) -> i64 {
    (end - start).abs()
}

/// -----------------------------------------------------------
/// 区间树查询辅助（简化版：线性扫描）
/// -----------------------------------------------------------
/// 参数:
///   query_chrom - 查询染色体
///   query_start - 查询起始位置
///   query_end   - 查询结束位置
///   intervals   - 数据库区间列表
/// 返回: Vec<(usize, String, i64, i64)> - 匹配的区间索引和信息
/// 说明:
///   这是一个简化版实现，使用线性扫描而非真正的区间树
///   对于大规模数据，应使用区间树或BEDTools
/// 用途: 查找与查询区间重叠的所有区间
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

/// -----------------------------------------------------------
/// Blacklist区域过滤
/// -----------------------------------------------------------
/// 参数:
///   intervals - 待过滤区间列表
///   blacklist - 黑名单区间列表
/// 返回: Vec<(String, i64, i64)> - 不在黑名单中的区间
/// 算法:
///   保留与任何黑名单区间都不重叠的区间
/// 用途:
///   ENCODE项目定义了基因组中的不可靠区域（Blacklist）
///   如重复序列、着丝粒、端粒等
///   在ChIP-seq/ATAC-seq分析中需要过滤这些区域
#[pyfunction]
pub fn filter_blacklist(
    intervals: Vec<(String, i64, i64)>,
    blacklist: Vec<(String, i64, i64)>,
) -> Vec<(String, i64, i64)> {
    intervals
        .into_iter()
        .filter(|(chrom, start, end)| {
            // 保留与所有黑名单区间都不重叠的区间
            !blacklist
                .iter()
                .any(|(bc, bs, be)| chrom == bc && *start < *be && *bs < *end)
        })
        .collect()
}

/// -----------------------------------------------------------
/// 最近特征查找
/// -----------------------------------------------------------
/// 参数:
///   query_chrom - 查询染色体
///   query_pos   - 查询位置
///   features    - 特征列表 (chrom, start, end, name)
/// 返回: Option<(名称, 距离, 方向)> - 最近特征信息
/// 算法:
///   1. 计算查询位置到每个特征的距离
///   2. 区间内：距离=0
///   3. 区间左侧：距离=start-pos
///   4. 区间右侧：距离=pos-end
///   5. 返回距离最小的特征
/// 用途: 寻找最近的基因、调控元件
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
            continue;  // 不同染色体，跳过
        }
        // 计算距离
        let dist = if query_pos < *start {
            *start - query_pos  // 位于特征左侧
        } else if query_pos > *end {
            query_pos - *end    // 位于特征右侧
        } else {
            0  // 在特征内部
        };
        if dist < min_dist {
            min_dist = dist;
            nearest = Some((name.clone(), min_dist, if query_pos < *start { "upstream".to_string() } else if query_pos > *end { "downstream".to_string() } else { "overlapping".to_string() }));
        }
    }
    nearest
}

/// -----------------------------------------------------------
/// 增强子-启动子匹配
/// -----------------------------------------------------------
/// 参数:
///   enhancers    - 增强子列表 (chrom, start, end, name)
///   promoters    - 启动子列表 (chrom, start, end, name)
///   max_distance - 最大允许距离
/// 返回: Vec<(String, String, i64)> - (增强子名, 启动子名, 距离)
/// 算法:
///   对于每个增强子，找到距离≤max_distance的启动子
///   距离计算：不重叠时为间隔距离，重叠时为0
/// 用途: 研究基因调控关系、构建调控网络
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
                continue;  // 不同染色体，跳过
            }
            // 计算距离
            let dist = if *ee < *ps {
                *ps - *ee  // 增强子在启动子左侧
            } else if *pe < *es {
                *es - *pe  // 增强子在启动子右侧
            } else {
                0  // 重叠
            };
            if dist <= max_distance {
                matches.push((eid.clone(), pid.clone(), dist));
            }
        }
    }
    matches
}

/// -----------------------------------------------------------
/// 区间覆盖度统计
/// -----------------------------------------------------------
/// 参数:
///   intervals      - 目标区间列表
///   read_positions - read位置列表 (chrom, pos)
/// 返回: Vec<i64> - 每个区间内的read数量
/// 算法:
///   对于每个区间，统计落在其中的read数
///   read在区间内的条件：同染色体且 start <= pos < end
/// 用途: 计算peak/gene的read覆盖度
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
