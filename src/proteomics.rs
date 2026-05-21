use pyo3::prelude::*;
use std::collections::HashMap;

/// Trypsin酶切模拟
#[pyfunction]
pub fn trypsin_digest(protein: &str, max_missed_cleavages: usize) -> Vec<String> {
    let mut cleavage_sites: Vec<usize> = vec![0];
    let chars: Vec<char> = protein.chars().collect();

    for i in 0..chars.len().saturating_sub(1) {
        if (chars[i] == 'K' || chars[i] == 'R') && chars[i + 1] != 'P' {
            cleavage_sites.push(i + 1);
        }
    }
    cleavage_sites.push(chars.len());

    let mut peptides = Vec::new();
    for mc in 0..=max_missed_cleavages {
        for i in 0..cleavage_sites.len().saturating_sub(1 + mc) {
            let start = cleavage_sites[i];
            let end = cleavage_sites[i + 1 + mc];
            if end > start && end <= chars.len() {
                let peptide: String = chars[start..end].iter().collect();
                if !peptide.is_empty() {
                    peptides.push(peptide);
                }
            }
        }
    }
    peptides
}

/// 肽段质量计算（单同位素质量）
#[pyfunction]
pub fn peptide_mass(peptide: &str) -> f64 {
    let aa_masses: HashMap<char, f64> = HashMap::from([
        ('G', 57.02146), ('A', 71.03711), ('V', 99.06841),
        ('L', 113.08406), ('I', 113.08406), ('P', 97.05276),
        ('F', 147.06841), ('W', 186.07931), ('M', 131.04049),
        ('S', 87.03203), ('T', 101.04768), ('C', 103.00919),
        ('Y', 163.06333), ('H', 137.05891), ('D', 115.02694),
        ('E', 129.04259), ('N', 114.04293), ('Q', 128.05858),
        ('K', 128.09496), ('R', 156.10111),
    ]);

    let mut mass = 18.01056; // H2O
    for c in peptide.chars() {
        mass += aa_masses.get(&c).unwrap_or(&0.0);
    }
    mass
}

/// Missed cleavage枚举
#[pyfunction]
pub fn enumerate_missed_cleavages(peptide: &str) -> Vec<String> {
    // 简单返回所有可能的子肽段
    let mut result = Vec::new();
    let chars: Vec<char> = peptide.chars().collect();
    for len in 1..=chars.len() {
        for start in 0..=chars.len() - len {
            let sub: String = chars[start..start + len].iter().collect();
            result.push(sub);
        }
    }
    result
}

/// Peptide uniqueness判断
#[pyfunction]
pub fn is_unique_peptide(peptide: &str, protein_database: Vec<String>) -> Vec<usize> {
    let mut protein_indices = Vec::new();
    for (i, protein) in protein_database.iter().enumerate() {
        if protein.contains(peptide) {
            protein_indices.push(i);
        }
    }
    protein_indices
}

/// Decoy database生成（反转序列）
#[pyfunction]
pub fn generate_decoy(sequences: Vec<String>) -> Vec<String> {
    sequences
        .iter()
        .map(|seq| seq.chars().rev().collect())
        .collect()
}

/// Modification site枚举
#[pyfunction]
pub fn enumerate_modification_sites(
    peptide: &str,
    mod_residues: Vec<String>,
) -> Vec<usize> {
    let mut sites = Vec::new();
    for (i, c) in peptide.chars().enumerate() {
        if mod_residues.iter().any(|s| s.chars().next() == Some(c)) {
            sites.push(i);
        }
    }
    sites
}

/// 肽段m/z计算
#[pyfunction]
pub fn peptide_mz(peptide: &str, charge: i32) -> f64 {
    let mass = peptide_mass(peptide);
    let proton = 1.007276;
    (mass + charge as f64 * proton) / charge as f64
}

/// FASTA蛋白库解析
#[pyfunction]
pub fn parse_fasta(path: &str) -> PyResult<Vec<(String, String)>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut sequences = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();

    for line in content.lines() {
        if line.starts_with('>') {
            if !current_header.is_empty() {
                sequences.push((current_header.clone(), current_seq.clone()));
            }
            current_header = line[1..].to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_header.is_empty() {
        sequences.push((current_header, current_seq));
    }
    Ok(sequences)
}

/// 蛋白质分子量
#[pyfunction]
pub fn protein_molecular_weight(protein: &str) -> f64 {
    let aa_weights: HashMap<char, f64> = HashMap::from([
        ('G', 57.052), ('A', 71.079), ('V', 99.133),
        ('L', 113.160), ('I', 113.160), ('P', 97.117),
        ('F', 147.177), ('W', 186.213), ('M', 131.199),
        ('S', 87.078), ('T', 101.105), ('C', 103.144),
        ('Y', 163.176), ('H', 137.142), ('D', 115.089),
        ('E', 129.116), ('N', 114.104), ('Q', 128.131),
        ('K', 128.174), ('R', 156.188),
    ]);

    let mut weight = 18.015; // H2O
    for c in protein.chars() {
        weight += aa_weights.get(&c).unwrap_or(&0.0);
    }
    weight
}

/// 肽段等电点（简化版）
#[pyfunction]
pub fn peptide_pi(peptide: &str) -> f64 {
    let mut positive = 0.0f64;
    let mut negative = 0.0f64;

    for c in peptide.chars() {
        match c {
            'K' | 'R' | 'H' => positive += 1.0,
            'D' | 'E' => negative += 1.0,
            _ => {}
        }
    }
    // 简化计算
    7.0 + (positive - negative) * 0.5
}
