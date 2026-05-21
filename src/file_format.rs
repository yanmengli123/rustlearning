use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

/// FASTQ转FASTA
#[pyfunction]
pub fn fastq_to_fasta(input: &str, output: &str) -> PyResult<usize> {
    let file = File::open(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
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
            continue;
        }
        if let (Some(Ok(seq)), Some(Ok(_)), Some(Ok(_qual))) =
            (lines.next(), lines.next(), lines.next())
        {
            writeln!(writer, ">{}", &header[1..]).unwrap();
            writeln!(writer, "{}", seq).unwrap();
            count += 1;
        }
    }
    Ok(count)
}

/// SAM转简化TSV
#[pyfunction]
pub fn sam_to_tsv(input: &str, output: &str) -> PyResult<usize> {
    let file = File::open(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let out_file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(out_file);

    writeln!(writer, "QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tSEQ").unwrap();

    let mut count = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 10 {
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

/// GTF转BED
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
                let name = rec
                    .attributes
                    .get("gene_id")
                    .cloned()
                    .unwrap_or_else(|| ".".to_string());
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

/// VCF转表格
#[pyfunction]
pub fn vcf_to_table(input: &str, output: &str) -> PyResult<usize> {
    let file = File::open(input)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);

    let out_file = File::create(output)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let mut writer = BufWriter::new(out_file);

    writeln!(writer, "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER").unwrap();

    let mut count = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 7 {
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

/// Count matrix转MTX格式
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

    writeln!(writer, "%%MatrixMarket matrix coordinate real general").unwrap();
    writeln!(writer, "{} {} {}", n_rows, n_cols, 0).unwrap(); // count non-zeros later

    let mut nnz = 0usize;
    for (i, row) in matrix.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            if val != 0.0 {
                writeln!(writer, "{} {} {}", i + 1, j + 1, val).unwrap();
                nnz += 1;
            }
        }
    }

    Ok(nnz)
}

/// 多样本计数矩阵合并
#[pyfunction]
pub fn merge_count_matrices(
    matrices: Vec<Vec<Vec<f64>>>,
) -> Vec<Vec<f64>> {
    if matrices.is_empty() {
        return Vec::new();
    }

    let n_genes = matrices[0].len();
    let total_cells: usize = matrices.iter().map(|m| if m.is_empty() { 0 } else { m[0].len() }).sum();

    let mut merged = vec![vec![0.0f64; total_cells]; n_genes];
    let mut col_offset = 0;

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

/// TSV/CSV表达矩阵解析
#[pyfunction]
pub fn parse_expression_matrix(
    path: &str,
    delimiter: char,
) -> PyResult<(Vec<String>, Vec<String>, Vec<Vec<f64>>)> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut lines = content.lines();
    let header = lines.next().unwrap_or("").to_string();
    let cell_names: Vec<String> = header.split(delimiter).skip(1).map(|s| s.to_string()).collect();

    let mut gene_names = Vec::new();
    let mut matrix = Vec::new();

    for line in lines {
        let fields: Vec<&str> = line.split(delimiter).collect();
        if fields.is_empty() {
            continue;
        }
        gene_names.push(fields[0].to_string());
        let row: Vec<f64> = fields[1..]
            .iter()
            .map(|s| s.parse().unwrap_or(0.0))
            .collect();
        matrix.push(row);
    }

    Ok((gene_names, cell_names, matrix))
}

/// 压缩写入
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

/// 读取压缩文件
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
