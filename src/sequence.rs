//! ============================================================
//! 模块1: 序列处理 (sequence)
//! ============================================================
//! 本模块提供DNA/RNA/蛋白质序列的基础处理功能。
//! 包括：GC含量计算、反向互补、碱基统计、序列验证、
//! 转录/翻译、ORF查找、低复杂度检测等。
//!
//! 设计原则：
//! - 所有输入序列统一转大写处理，保证结果一致
//! - 空序列返回安全默认值（0.0、空字符串等），不报错
//! - 非法输入（如非ATGCN字符）由调用方自行检查
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use pyo3::types::PyDict;        // Python 字典类型
use std::collections::HashMap;  // Rust 标准哈希表

/// -----------------------------------------------------------
/// 计算GC含量
/// -----------------------------------------------------------
/// 参数:
///   seq - DNA/RNA序列字符串
/// 返回:
///   GC含量比例 (0.0 ~ 1.0)
/// 算法:
///   统计 G 和 C 的数量，除以总碱基数。
///   空序列返回 0.0 避免除零错误。
/// 用途:
///   质控评估、物种鉴定、测序偏差检测
#[pyfunction]
pub fn gc_content(seq: &str) -> f64 {
    let seq = seq.to_uppercase();  // 统一转大写，处理大小写混合输入
    let total = seq.len() as f64;  // 序列总长度
    if total == 0.0 {
        return 0.0;  // 空序列安全返回
    }
    // 过滤出所有 G 和 C 碱基，计数
    let gc = seq.chars().filter(|&c| c == 'G' || c == 'C').count() as f64;
    gc / total  // 返回比例
}

/// -----------------------------------------------------------
/// 反向互补序列
/// -----------------------------------------------------------
/// 参数:
///   seq - DNA序列
/// 返回:
///   反向互补序列
/// 算法:
///   1. 将序列反转
///   2. 每个碱基取互补碱基:
///      A <-> T, G <-> C, N <-> N, U -> A
/// 用途:
///   引物设计、序列比对、双链处理
#[pyfunction]
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()  // 反转序列
        .map(|c| match c {
            'A' | 'a' => 'T',   // A 互补为 T
            'T' | 't' => 'A',   // T 互补为 A
            'G' | 'g' => 'C',   // G 互补为 C
            'C' | 'c' => 'G',   // C 互补为 G
            'N' | 'n' => 'N',   // N 保持 N（未知碱基）
            'U' | 'u' => 'A',   // RNA的U互补为A
            _ => c,             // 其他字符保持不变
        })
        .collect()  // 收集为字符串
}

/// -----------------------------------------------------------
/// 统计各碱基数量
/// -----------------------------------------------------------
/// 参数:
///   seq - 序列字符串
/// 返回:
///   HashMap<char, usize> - 每种碱基的计数
/// 算法:
///   遍历序列每个字符，使用 HashMap 累加计数。
///   大小写统一转大写后统计。
/// 用途:
///   碱基组成分析、测序质量评估
#[pyfunction]
pub fn count_bases(seq: &str) -> HashMap<char, usize> {
    let mut counts = HashMap::new();  // 创建空哈希表
    for c in seq.to_uppercase().chars() {
        // entry API: 如果键不存在则插入0，然后+1
        *counts.entry(c).or_insert(0) += 1;
    }
    counts
}

/// -----------------------------------------------------------
/// 序列长度
/// -----------------------------------------------------------
/// 参数: seq - 序列字符串
/// 返回: 序列长度（字节数，对于ASCII字符等于字符数）
#[pyfunction]
pub fn seq_length(seq: &str) -> usize {
    seq.len()
}

/// -----------------------------------------------------------
/// 检查是否为合法DNA序列
/// -----------------------------------------------------------
/// 参数: seq - 序列字符串
/// 返回: true 如果只包含 A/T/G/C/N（大小写均可）
/// 用途: 输入验证、FASTA文件预检查
#[pyfunction]
pub fn is_valid_dna(seq: &str) -> bool {
    !seq.is_empty()  // 非空
        && seq.chars().all(|c| matches!(c,
            'A' | 'a' | 'T' | 't' |   // 标准DNA碱基
            'G' | 'g' | 'C' | 'c' |
            'N' | 'n'                   // N表示未知碱基
        ))
}

/// -----------------------------------------------------------
/// 检查是否为合法RNA序列
/// -----------------------------------------------------------
/// 参数: seq - 序列字符串
/// 返回: true 如果只包含 A/U/G/C/N（大小写均可）
/// 说明: RNA使用U（尿嘧啶）代替T（胸腺嘧啶）
#[pyfunction]
pub fn is_valid_rna(seq: &str) -> bool {
    !seq.is_empty()
        && seq.chars().all(|c| matches!(c,
            'A' | 'a' | 'U' | 'u' |   // RNA的U代替T
            'G' | 'g' | 'C' | 'c' |
            'N' | 'n'
        ))
}

/// -----------------------------------------------------------
/// DNA转RNA（转录）
/// -----------------------------------------------------------
/// 参数: seq - DNA序列
/// 返回: RNA序列（T -> U）
/// 生物学意义: DNA模板链转录为mRNA时，T变为U
#[pyfunction]
pub fn transcribe(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'T' => 'U',   // 胸腺嘧啶 -> 尿嘧啶
            't' => 'u',
            _ => c,
        })
        .collect()
}

/// -----------------------------------------------------------
/// RNA转DNA（逆转录）
/// -----------------------------------------------------------
/// 参数: seq - RNA序列
/// 返回: DNA序列（U -> T）
/// 生物学意义: 逆转录酶将RNA逆转录为cDNA
#[pyfunction]
pub fn reverse_transcribe(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'U' => 'T',   // 尿嘧啶 -> 胸腺嘧啶
            'u' => 't',
            _ => c,
        })
        .collect()
}

/// -----------------------------------------------------------
/// 统计N碱基数
/// -----------------------------------------------------------
/// 参数: seq - 序列字符串
/// 返回: N碱基的数量
/// 用途: 评估测序质量，N表示该位置无法确定碱基
#[pyfunction]
pub fn count_n(seq: &str) -> usize {
    seq.chars().filter(|&c| c == 'N' || c == 'n').count()
}

/// -----------------------------------------------------------
/// 碱基大小写标准化（转大写）
/// -----------------------------------------------------------
/// 参数: seq - 序列字符串
/// 返回: 全大写的序列
/// 用途: 统一处理不同来源的序列数据
#[pyfunction]
pub fn normalize_seq(seq: &str) -> String {
    seq.to_uppercase()
}

/// -----------------------------------------------------------
/// 序列切片
/// -----------------------------------------------------------
/// 参数:
///   seq   - 原始序列
///   start - 起始位置（包含）
///   end   - 结束位置（不包含）
/// 返回: 切片后的子序列
/// 错误: 如果范围越界，抛出 IndexError
/// 用途: 提取序列的特定区域（如引物区域、CDS区域）
#[pyfunction]
pub fn seq_slice(seq: &str, start: usize, end: usize) -> PyResult<String> {
    // 边界检查
    if start > seq.len() || end > seq.len() || start > end {
        return Err(pyo3::exceptions::PyIndexError::new_err(format!(
            "Invalid range: {}..{} for seq of length {}",
            start, end, seq.len()
        )));
    }
    Ok(seq[start..end].to_string())
}

/// -----------------------------------------------------------
/// 拼接多条序列
/// -----------------------------------------------------------
/// 参数:
///   seqs      - 序列列表
///   separator - 分隔符（可以为空字符串）
/// 返回: 拼接后的序列
/// 用途: 合并contig、构建consensus序列
#[pyfunction]
pub fn concat_seqs(seqs: Vec<String>, separator: &str) -> String {
    seqs.join(separator)
}

/// -----------------------------------------------------------
/// DNA序列信息压缩比
/// -----------------------------------------------------------
/// 参数: seq - 序列字符串
/// 返回: 压缩比 (0.0 ~ 1.0)
/// 算法:
///   计算序列的Shannon熵，除以最大熵（log2(4)=2）。
///   碱基分布越均匀，熵越高，压缩比越接近1。
///   低熵序列（如polyA）压缩比低。
/// 用途: 评估序列复杂度、检测低复杂度区域
#[pyfunction]
pub fn compression_ratio(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let bases = count_bases(seq);
    let len = seq.len() as f64;
    let mut entropy = 0.0;
    for (_, &count) in &bases {
        let p = count as f64 / len;  // 每种碱基的频率
        if p > 0.0 {
            entropy -= p * p.log2();  // Shannon熵公式: -Σ p*log2(p)
        }
    }
    entropy / 2.0  // DNA有4种碱基，最大熵为log2(4)=2
}

/// -----------------------------------------------------------
/// 低复杂度序列检测
/// -----------------------------------------------------------
/// 参数:
///   seq       - 序列字符串
///   threshold - 阈值 (0.0 ~ 1.0)
/// 返回: true 如果任一碱基占比超过阈值
/// 用途: 过滤低复杂度序列（如microsatellite、polyA尾）
#[pyfunction]
pub fn is_low_complexity(seq: &str, threshold: f64) -> bool {
    if seq.is_empty() {
        return true;  // 空序列视为低复杂度
    }
    let bases = count_bases(seq);
    // 找出出现次数最多的碱基
    let max_count = bases.values().max().unwrap_or(&0);
    // 如果最频繁碱基的占比超过阈值，认为是低复杂度
    (*max_count as f64 / seq.len() as f64) > threshold
}

/// -----------------------------------------------------------
/// 寻找ORF（开放阅读框）
/// -----------------------------------------------------------
/// 参数: seq - DNA序列
/// 返回: Vec<(start, end, orf_seq)> - ORF列表
/// 算法:
///   在3个reading frame中扫描：
///   1. 找到起始密码子 ATG
///   2. 每次移动3个碱基（一个密码子）
///   3. 遇到终止密码子 TAA/TAG/TGA 记录ORF
///   4. 继续扫描下一个ATG
/// 用途: 基因预测、蛋白质编码区域识别
#[pyfunction]
pub fn find_orfs(seq: &str) -> Vec<(usize, usize, String)> {
    let seq = seq.to_uppercase();
    let bytes = seq.as_bytes();
    let mut orfs = Vec::new();
    // 三个终止密码子
    let stop_codons: &[&[u8]] = &[b"TAA", b"TAG", b"TGA"];

    // 遍历3个reading frame (0, 1, 2)
    for frame in 0..3 {
        let mut i = frame;
        while i + 3 <= bytes.len() {
            // 找到起始密码子 ATG
            if &bytes[i..i + 3] == b"ATG" {
                let start = i;
                let mut j = i + 3;
                // 向后扫描，每次移动3个碱基
                while j + 3 <= bytes.len() {
                    // 检查是否为终止密码子
                    if stop_codons.contains(&&bytes[j..j + 3]) {
                        // 记录ORF（包含终止密码子）
                        let orf_seq = String::from_utf8_lossy(&bytes[start..j + 3]).to_string();
                        orfs.push((start, j + 3, orf_seq));
                        break;
                    }
                    j += 3;
                }
            }
            i += 3;
        }
    }
    orfs
}

/// -----------------------------------------------------------
/// 密码子翻译为氨基酸
/// -----------------------------------------------------------
/// 参数: seq - DNA序列（密码子序列）
/// 返回: 氨基酸序列（单字母缩写）
/// 算法:
///   使用标准遗传密码表，每次读取3个碱基，
///   查表转换为对应的氨基酸。
///   遇到终止密码子(*)停止翻译。
///   无法识别的密码子翻译为X。
/// 用途: 蛋白质序列预测、CDS翻译
#[pyfunction]
pub fn translate(seq: &str) -> String {
    // 标准遗传密码表
    let codon_table: HashMap<&str, char> = HashMap::from([
        ("TTT", 'F'), ("TTC", 'F'), ("TTA", 'L'), ("TTG", 'L'),  // Phe, Leu
        ("CTT", 'L'), ("CTC", 'L'), ("CTA", 'L'), ("CTG", 'L'),  // Leu
        ("ATT", 'I'), ("ATC", 'I'), ("ATA", 'I'), ("ATG", 'M'),  // Ile, Met(起始)
        ("GTT", 'V'), ("GTC", 'V'), ("GTA", 'V'), ("GTG", 'V'),  // Val
        ("TCT", 'S'), ("TCC", 'S'), ("TCA", 'S'), ("TCG", 'S'),  // Ser
        ("CCT", 'P'), ("CCC", 'P'), ("CCA", 'P'), ("CCG", 'P'),  // Pro
        ("ACT", 'T'), ("ACC", 'T'), ("ACA", 'T'), ("ACG", 'T'),  // Thr
        ("GCT", 'A'), ("GCC", 'A'), ("GCA", 'A'), ("GCG", 'A'),  // Ala
        ("TAT", 'Y'), ("TAC", 'Y'), ("TAA", '*'), ("TAG", '*'),  // Tyr, 终止
        ("CAT", 'H'), ("CAC", 'H'), ("CAA", 'Q'), ("CAG", 'Q'),  // His, Gln
        ("AAT", 'N'), ("AAC", 'N'), ("AAA", 'K'), ("AAG", 'K'),  // Asn, Lys
        ("GAT", 'D'), ("GAC", 'D'), ("GAA", 'E'), ("GAG", 'E'),  // Asp, Glu
        ("TGT", 'C'), ("TGC", 'C'), ("TGA", '*'), ("TGG", 'W'),  // Cys, 终止, Trp
        ("CGT", 'R'), ("CGC", 'R'), ("CGA", 'R'), ("CGG", 'R'),  // Arg
        ("AGT", 'S'), ("AGC", 'S'), ("AGA", 'R'), ("AGG", 'R'),  // Ser, Arg
        ("GGT", 'G'), ("GGC", 'G'), ("GGA", 'G'), ("GGG", 'G'),  // Gly
    ]);
    let seq = seq.to_uppercase();
    let mut protein = String::new();
    let bytes = seq.as_bytes();
    let mut i = 0;
    // 每次读取3个碱基（一个密码子）
    while i + 3 <= bytes.len() {
        let codon = &seq[i..i + 3];
        if let Some(&aa) = codon_table.get(codon) {
            protein.push(aa);
            if aa == '*' {
                break;  // 遇到终止密码子停止翻译
            }
        } else {
            protein.push('X');  // 无法识别的密码子用X表示
        }
        i += 3;
    }
    protein
}

/// -----------------------------------------------------------
/// 统计序列综合信息
/// -----------------------------------------------------------
/// 参数:
///   py  - Python解释器引用（PyO3自动传入）
///   seq - 序列字符串
/// 返回: Python字典，包含各项统计指标
/// 用途: 快速获取序列的全面统计信息
#[pyfunction]
pub fn seq_stats(py: Python, seq: &str) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);  // 创建Python字典
    let bases = count_bases(seq);

    // 设置各项统计指标
    dict.set_item("length", seq.len())?;           // 序列长度
    dict.set_item("gc_content", gc_content(seq))?; // GC含量
    dict.set_item("n_count", count_n(seq))?;       // N碱基数
    dict.set_item("is_valid_dna", is_valid_dna(seq))?; // 是否合法DNA

    // 各碱基计数（不存在则返回0）
    dict.set_item("a_count", bases.get(&'A').unwrap_or(&0))?;
    dict.set_item("t_count", bases.get(&'T').unwrap_or(&0))?;
    dict.set_item("g_count", bases.get(&'G').unwrap_or(&0))?;
    dict.set_item("c_count", bases.get(&'C').unwrap_or(&0))?;

    Ok(dict.into())  // 转为PyObject返回
}
