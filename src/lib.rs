use pyo3::prelude::*;

mod sequence;
mod fastq;
mod kmer;
mod alignment;
mod sam_bam;
mod bed;
mod genomic_intervals;
mod rnaseq;
mod single_cell;
mod epigenomics;
mod metagenomics;
mod proteomics;
mod statistics;
mod file_format;
mod parallel;

#[pymodule]
fn my_python_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // === 1. 序列处理 (sequence) ===
    m.add_function(wrap_pyfunction!(sequence::gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::count_bases, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::seq_length, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::is_valid_dna, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::is_valid_rna, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::transcribe, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::reverse_transcribe, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::count_n, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::normalize_seq, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::seq_slice, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::concat_seqs, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::compression_ratio, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::is_low_complexity, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::find_orfs, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::translate, m)?)?;
    m.add_function(wrap_pyfunction!(sequence::seq_stats, m)?)?;

    // === 2. FASTQ质控 (fastq) ===
    m.add_function(wrap_pyfunction!(fastq::fastq_stats, m)?)?;
    m.add_function(wrap_pyfunction!(fastq::fastq_filter, m)?)?;
    m.add_function(wrap_pyfunction!(fastq::per_base_quality, m)?)?;
    m.add_function(wrap_pyfunction!(fastq::length_distribution, m)?)?;
    m.add_function(wrap_pyfunction!(fastq::gc_distribution, m)?)?;
    m.add_function(wrap_pyfunction!(fastq::sliding_window_filter, m)?)?;

    // === 3. k-mer与sketch算法 (kmer) ===
    m.add_function(wrap_pyfunction!(kmer::count_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::canonical_kmer, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::minimizers, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::minhash_sketch, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::jaccard_similarity, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::containment_index, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::mash_distance, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::kmer_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::syncmers, m)?)?;

    // === 4. 比对算法 (alignment) ===
    m.add_function(wrap_pyfunction!(alignment::hamming_distance, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::levenshtein_distance, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::needleman_wunsch, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::smith_waterman, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::banded_alignment, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::parse_cigar, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::cigar_stats, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::primer_match, m)?)?;

    // === 5. SAM/BAM处理 (sam_bam) ===
    m.add_function(wrap_pyfunction!(sam_bam::parse_flag, m)?)?;
    m.add_function(wrap_pyfunction!(sam_bam::sam_stats, m)?)?;
    m.add_function(wrap_pyfunction!(sam_bam::filter_by_mapq, m)?)?;
    m.add_function(wrap_pyfunction!(sam_bam::fetch_region, m)?)?;
    m.add_function(wrap_pyfunction!(sam_bam::coverage_at_position, m)?)?;
    m.add_function(wrap_pyfunction!(sam_bam::region_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(sam_bam::insert_size_stats, m)?)?;

    // === 6. BED/GTF/VCF处理 (bed) ===
    m.add_function(wrap_pyfunction!(bed::bed_stats, m)?)?;
    m.add_function(wrap_pyfunction!(bed::bed_intersect, m)?)?;
    m.add_function(wrap_pyfunction!(bed::bed_merge, m)?)?;
    m.add_function(wrap_pyfunction!(bed::bed_subtract, m)?)?;
    m.add_function(wrap_pyfunction!(bed::bed_closest, m)?)?;
    m.add_function(wrap_pyfunction!(bed::vcf_stats, m)?)?;
    m.add_function(wrap_pyfunction!(bed::vcf_filter, m)?)?;
    m.add_function(wrap_pyfunction!(bed::gtf_stats, m)?)?;

    // === 7. 基因组区间算法 (genomic_intervals) ===
    m.add_function(wrap_pyfunction!(genomic_intervals::intervals_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::overlap_length, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::batch_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::sliding_bins, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::window_count, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::tss_distance, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::interval_length, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::interval_tree_query, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::filter_blacklist, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::nearest_feature, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::match_enhancer_promoter, m)?)?;
    m.add_function(wrap_pyfunction!(genomic_intervals::coverage_over_intervals, m)?)?;

    // === 8. RNA-seq分析 (rnaseq) ===
    m.add_function(wrap_pyfunction!(rnaseq::gene_count_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(rnaseq::exon_count, m)?)?;
    m.add_function(wrap_pyfunction!(rnaseq::splice_junction_stats, m)?)?;
    m.add_function(wrap_pyfunction!(rnaseq::detect_intron_retention, m)?)?;
    m.add_function(wrap_pyfunction!(rnaseq::umi_collapse, m)?)?;
    m.add_function(wrap_pyfunction!(rnaseq::gene_biotype_stats, m)?)?;
    m.add_function(wrap_pyfunction!(rnaseq::transcript_length_stats, m)?)?;

    // === 9. 单细胞测序 (single_cell) ===
    m.add_function(wrap_pyfunction!(single_cell::extract_barcode, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::extract_umi, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::correct_barcode, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::barcode_rescue, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::umi_dedup, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::build_feature_barcode_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::to_matrix_market, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::cell_qc, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::demux_tags, m)?)?;
    m.add_function(wrap_pyfunction!(single_cell::count_guides, m)?)?;

    // === 10. 表观组学 (epigenomics) ===
    m.add_function(wrap_pyfunction!(epigenomics::tn5_insertion_sites, m)?)?;
    m.add_function(wrap_pyfunction!(epigenomics::count_reads_in_peaks, m)?)?;
    m.add_function(wrap_pyfunction!(epigenomics::fragment_length_distribution, m)?)?;
    m.add_function(wrap_pyfunction!(epigenomics::classify_nucleosome_free, m)?)?;
    m.add_function(wrap_pyfunction!(epigenomics::chipseq_peak_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(epigenomics::cpg_methylation_levels, m)?)?;
    m.add_function(wrap_pyfunction!(epigenomics::fragment_overlap_peaks, m)?)?;
    m.add_function(wrap_pyfunction!(epigenomics::coverage_over_bins, m)?)?;
    m.add_function(wrap_pyfunction!(epigenomics::cuttag_fragment_stats, m)?)?;

    // === 11. 宏基因组 (metagenomics) ===
    m.add_function(wrap_pyfunction!(metagenomics::kmer_classify, m)?)?;
    m.add_function(wrap_pyfunction!(metagenomics::ani_approximate, m)?)?;
    m.add_function(wrap_pyfunction!(metagenomics::genome_sketch, m)?)?;
    m.add_function(wrap_pyfunction!(metagenomics::preprocess_16s, m)?)?;
    m.add_function(wrap_pyfunction!(metagenomics::marker_gene_match, m)?)?;
    m.add_function(wrap_pyfunction!(metagenomics::contamination_screen, m)?)?;
    m.add_function(wrap_pyfunction!(metagenomics::greedy_cluster, m)?)?;
    m.add_function(wrap_pyfunction!(metagenomics::estimate_abundance, m)?)?;

    // === 12. 蛋白质组学 (proteomics) ===
    m.add_function(wrap_pyfunction!(proteomics::trypsin_digest, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::peptide_mass, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::enumerate_missed_cleavages, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::is_unique_peptide, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::generate_decoy, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::enumerate_modification_sites, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::peptide_mz, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::parse_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::protein_molecular_weight, m)?)?;
    m.add_function(wrap_pyfunction!(proteomics::peptide_pi, m)?)?;

    // === 13. 统计计算 (statistics) ===
    m.add_function(wrap_pyfunction!(statistics::permutation_test, m)?)?;
    m.add_function(wrap_pyfunction!(statistics::bootstrap_ci, m)?)?;
    m.add_function(wrap_pyfunction!(statistics::hypergeometric_test, m)?)?;
    m.add_function(wrap_pyfunction!(statistics::fisher_exact, m)?)?;
    m.add_function(wrap_pyfunction!(statistics::bh_correction, m)?)?;
    m.add_function(wrap_pyfunction!(statistics::rolling_correlation, m)?)?;
    m.add_function(wrap_pyfunction!(statistics::pearson_correlation, m)?)?;
    m.add_function(wrap_pyfunction!(statistics::mann_whitney_u, m)?)?;
    m.add_function(wrap_pyfunction!(statistics::descriptive_stats, m)?)?;

    // === 14. 文件格式转换 (file_format) ===
    m.add_function(wrap_pyfunction!(file_format::fastq_to_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(file_format::sam_to_tsv, m)?)?;
    m.add_function(wrap_pyfunction!(file_format::gtf_to_bed, m)?)?;
    m.add_function(wrap_pyfunction!(file_format::vcf_to_table, m)?)?;
    m.add_function(wrap_pyfunction!(file_format::count_matrix_to_mtx, m)?)?;
    m.add_function(wrap_pyfunction!(file_format::merge_count_matrices, m)?)?;
    m.add_function(wrap_pyfunction!(file_format::parse_expression_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(file_format::write_gz, m)?)?;
    m.add_function(wrap_pyfunction!(file_format::read_gz, m)?)?;

    // === 15. 并行批处理 (parallel) ===
    m.add_function(wrap_pyfunction!(parallel::parallel_gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_count_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_hamming_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_levenshtein_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_seq_stats, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_find_orfs, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_parse_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_bed_intersect, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_jaccard_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::batch_fastq_filter, m)?)?;
    m.add_function(wrap_pyfunction!(parallel::parallel_descriptive_stats, m)?)?;

    Ok(())
}
