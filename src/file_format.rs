//! ============================================================
//! 模块14: 文件格式转换 (file_format)
//! ============================================================
//! 本模块提供生物信息学文件格式的转换功能。
//! 包括：FASTQ→FASTA、SAM→TSV、GTF→BED、VCF→表格、
//! 计数矩阵→MTX、矩阵合并、表达矩阵解析、gzip读写等。
//!
//! 常见生物信息学文件格式：
//! - FASTQ: 测序原始数据（含质量分数）
//! - FASTA: 序列数据（无质量分数）
//! - SAM/BAM: 比对结果
//! - BED: 基因组区间
//! - GTF/GFF: 基因注释
//! - VCF: 变异信息
//! - MTX: 稀疏矩阵（Matrix Market格式）
//!
//! 设计原则：
//! - 支持gzip压缩文件的自动检测
//! - 使用缓冲读写提高性能
//! - 返回转换的记录数
//! ============================================================

use pyo3::prelude::*;           // PyO3 核心宏和类型
use std::fs::File;              // 文件操作
use std::io::{BufRead, BufReader, Write, BufWriter};  // IO操作
use flate2::read::GzDecoder;    // gzip解码器
use flate2::write::GzEncoder;   // gzip编码器
use flate2::Compression;        // 压缩级别

/// -----------------------------------------------------------
/// FASTQ转FASTA
/// -----------------------------------------------------------
/// 参数:
///   input  - 输入FASTQ文件路径（支持.gz）
///   output - 输出FASTA文件路径
/// 返回: 转换的序列数量
/// 算法:
///   1. 读取FASTQ文件（每4行为一条记录）
///   2. 去掉@前缀，改为>前缀
///   3. 丢弃质量分数行
/// 用途: 当不需要质量信息时，转为更紧凑的FASTA格式
#[pyfunction]
pub fn fastq_to_fasta(input: &str, output: &str) -> PyResult<usize> {
    let file = File::open(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    // 根据后缀判断是否为gzip
    let reader: Box<dyn BufRead> = if input.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let out_file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(out_file);

    let mut count = 0usize;
    let mut lines = reader.lines();
    while let Some(Ok(header)) = lines.next() {
        if !header.starts_with('@') {
            continue;  // 跳过非header行
        }
        // 读取4行：header, seq, +, qual
        if let (Some(Ok(seq)), Some(Ok(_)), Some(Ok(_qual))) =
            (lines.next(), lines.next(), lines.next())
        {
            // FASTA格式：>header（去掉@改为>）
            writeln!(writer, ">{}", &header[1..]).unwrap();
            writeln!(writer, "{}", seq).unwrap();
            count += 1;
        }
    }
    Ok(count)
}

/// -----------------------------------------------------------
/// SAM转简化TSV
/// -----------------------------------------------------------
/// 参数:
///   input  - 输入SAM文件路径
///   output - 输出TSV文件路径
/// 返回: 转换的记录数
/// 算法:
///   提取SAM记录的7个关键字段：
///   QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, SEQ
///   输出为制表符分隔的文本文件
/// 用途: 简化SAM文件用于下游分析
#[pyfunction]
pub fn sam_to_tsv(input: &str, output: &str) -> PyResult<usize> {
    let file = File::open(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let out_file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(out_file);

    // 写入表头
    writeln!(writer, "QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tSEQ").unwrap();

    let mut count = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;  // 跳过头部
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 10 {
            // 提取7个字段
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[9]
            )
            .unwrap();
            count += 1;
        }
    }
    Ok(count)
}

/// -----------------------------------------------------------
/// GTF转BED
/// -----------------------------------------------------------
/// 参数:
///   input   - 输入GTF文件路径
///   output  - 输出BED文件路径
///   feature - 特征类型过滤（如"gene", "exon"）
/// 返回: 转换的记录数
/// 算法:
///   1. 解析GTF文件
///   2. 只保留指定feature类型的记录
///   3. 转换为BED格式（chrom, start, end, name, score, strand）
///   4. GTF的1-based坐标转为BED的0-based坐标（start-1）
/// 用途: 使用BEDTools等工具进行区间操作
#[pyfunction]
pub fn gtf_to_bed(input: &str, output: &str, feature: &str) -> PyResult<usize> {
    let file = File::open(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let out_file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(out_file);

    let mut count = 0usize;
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = super::bed::parse_gtf_line(&line) {
            if rec.feature == feature {
                // 提取gene_id作为名称
                let name = rec
                    .attributes
                    .get("gene_id")
                    .cloned()
                    .unwrap_or_else(|| ".".to_string());
                // GTF坐标转BED坐标（start-1）
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t0\t{}",
                    rec.chrom, rec.start, rec.end, name, rec.strand
                )
                .unwrap();
                count += 1;
            }
        }
    }
    Ok(count)
}

/// -----------------------------------------------------------
/// VCF转表格
/// -----------------------------------------------------------
/// 参数:
///   input  - 输入VCF文件路径
///   output - 输出TSV文件路径
/// 返回: 转换的记录数
/// 算法:
///   提取VCF的7个关键字段：
///   CHROM, POS, ID, REF, ALT, QUAL, FILTER
///   输出为制表符分隔的文本文件
/// 用途: 简化VCF文件用于数据分析
#[pyfunction]
pub fn vcf_to_table(input: &str, output: &str) -> PyResult<usize> {
    let file = File::open(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let out_file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(out_file);

    // 写入表头
    writeln!(writer, "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER").unwrap();

    let mut count = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;  // 跳过注释和空行
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 7 {
            // 提取7个字段
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6]
            )
            .unwrap();
            count += 1;
        }
    }
    Ok(count)
}

/// -----------------------------------------------------------
/// 计数矩阵转MTX格式
/// -----------------------------------------------------------
/// 参数:
///   matrix - 二维矩阵（行=基因，列=细胞）
///   output - 输出MTX文件路径
/// 返回: 非零元素数量
/// 算法:
///   1. 写入Matrix Market头部
///   2. 遍历矩阵，只输出非零元素
///   3. 使用1-based索引
/// 用途: 生成标准MTX格式文件（用于Scanpy/Seurat）
#[pyfunction]
pub fn count_matrix_to_mtx(
    matrix: Vec<Vec<f64>>,
    output: &str,
) -> PyResult<usize> {
    let out_file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(out_file);

    let n_rows = matrix.len();
    let n_cols = if n_rows > 0 { matrix[0].len() } else { 0 };

    // 写入Matrix Market头部
    writeln!(writer, "%%MatrixMarket matrix coordinate real general").unwrap();
    writeln!(writer, "{} {} {}", n_rows, n_cols, 0).unwrap(); // 先写占位符

    let mut nnz = 0usize;
    for (i, row) in matrix.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            if val != 0.0 {
                // 1-based索引
                writeln!(writer, "{} {} {}", i + 1, j + 1, val).unwrap();
                nnz += 1;
            }
        }
    }

    Ok(nnz)
}

/// -----------------------------------------------------------
/// 多样本计数矩阵合并
/// -----------------------------------------------------------
/// 参数: matrices - 多个计数矩阵列表
/// 返回: 合并后的大矩阵
/// 算法:
///   将多个矩阵按列（细胞）拼接
///   假设所有矩阵的行（基因）顺序相同
/// 用途: 合并多个样本的单细胞数据
#[pyfunction]
pub fn merge_count_matrices(
    matrices: Vec<Vec<Vec<f64>>>,
) -> Vec<Vec<f64>> {
    if matrices.is_empty() {
        return Vec::new();
    }

    let n_genes = matrices[0].len();
    // 计算总细胞数
    let total_cells: usize = matrices.iter().map(|m| if m.is_empty() { 0 } else { m[0].len() }).sum();

    // 创建合并矩阵
    let mut merged = vec![vec![0.0f64; total_cells]; n_genes];
    let mut col_offset = 0;

    // 逐个矩阵复制
    for matrix in &matrices {
        if matrix.is_empty() {
            continue;
        }
        let n_cols = matrix[0].len();
        for (i, row) in matrix.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                if i < n_genes && col_offset + j < total_cells {
                    merged[i][col_offset + j] = val;
                }
            }
        }
        col_offset += n_cols;
    }
    merged
}

/// -----------------------------------------------------------
/// TSV/CSV表达矩阵解析
/// -----------------------------------------------------------
/// 参数:
///   path      - 文件路径
///   delimiter - 分隔符（'\t'或','）
/// 返回: (基因名列表, 细胞名列表, 矩阵)
/// 算法:
///   1. 第一行为细胞名（跳过第一列）
///   2. 后续行：第一列为基因名，其余为表达值
/// 用途: 加载表达矩阵用于下游分析
#[pyfunction]
pub fn parse_expression_matrix(
    path: &str,
    delimiter: char,
) -> PyResult<(Vec<String>, Vec<String>, Vec<Vec<f64>>)> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut lines = content.lines();
    // 第一行是表头（细胞名）
    let header = lines.next().unwrap_or("").to_string();
    let cell_names: Vec<String> = header.split(delimiter).skip(1).map(|s| s.to_string()).collect();

    let mut gene_names = Vec::new();
    let mut matrix = Vec::new();

    // 后续行是基因数据
    for line in lines {
        let fields: Vec<&str> = line.split(delimiter).collect();
        if fields.is_empty() {
            continue;
        }
        gene_names.push(fields[0].to_string());  // 第一列是基因名
        // 后续列是表达值
        let row: Vec<f64> = fields[1..]
            .iter()
            .map(|s| s.parse().unwrap_or(0.0))
            .collect();
        matrix.push(row);
    }

    Ok((gene_names, cell_names, matrix))
}

/// -----------------------------------------------------------
/// 压缩写入（gzip格式）
/// -----------------------------------------------------------
/// 参数:
///   content - 要写入的文本内容
///   output  - 输出.gz文件路径
/// 用途: 压缩大文件节省存储空间
#[pyfunction]
pub fn write_gz(content: &str, output: &str) -> PyResult<()> {
    let file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder
        .write_all(content.as_bytes())
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    encoder
        .finish()
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    Ok(())
}

/// -----------------------------------------------------------
/// 读取压缩文件（gzip格式）
/// -----------------------------------------------------------
/// 参数: path - .gz文件路径
/// 返回: 解压后的文本内容
/// 用途: 读取压缩的测序数据文件
#[pyfunction]
pub fn read_gz(path: &str) -> PyResult<String> {
    let file = File::open(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut decoder = GzDecoder::new(file);
    let mut content = String::new();
    std::io::Read::read_to_string(&mut decoder, &mut content)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    Ok(content)
}
