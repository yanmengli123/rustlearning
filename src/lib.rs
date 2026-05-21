//! ============================================================
//! 生物信息学Rust扩展模块 - 主入口文件
//! ============================================================
//! 本文件是Rust + PyO3生物信息学扩展库的入口文件。
//! 负责：
//!   1. 声明所有子模块
//!   2. 注册所有Python可调用函数
//!   3. 创建Python模块
//!
//! 模块列表（15个模块，100+函数）：
//!   1.  sequence          - 序列处理（GC含量、反向互补、转录翻译等）
//!   2.  fastq             - FASTQ质量控制和过滤
//!   3.  kmer              - k-mer分析和MinHash草图
//!   4.  alignment         - 序列比对算法（NW、SW、编辑距离等）
//!   5.  sam_bam           - SAM/BAM文件解析和统计
//!   6.  bed               - BED/GTF/VCF区间操作
//!   7.  genomic_intervals - 基因组区间高级操作
//!   8.  rnaseq            - RNA-seq分析
//!   9.  single_cell       - 单细胞RNA-seq分析
//!   10. epigenomics       - 表观基因组学分析
//!   11. metagenomics      - 宏基因组学分析
//!   12. proteomics        - 蛋白质组学分析
//!   13. statistics        - 统计分析
//!   14. file_format       - 文件格式转换
//!   15. parallel          - 批量/并行处理
//!
//! 使用方法：
//!   import my_python_module as bio
//!   gc = bio.gc_content("ATCGATCG")
//!
//! 构建方法：
//!   maturin develop --release
//! ============================================================

use pyo3::prelude::*;  // PyO3 核心宏和类型

// ============================================================
// 模块声明
// ============================================================
// 每个mod声明对应src目录下的一个.rs文件
// 这些模块包含具体的函数实现

mod sequence;           // 模块1: 序列处理
mod fastq;              // 模块2: FASTQ质控
mod kmer;               // 模块3: k-mer分析
mod alignment;          // 模块4: 序列比对
mod sam_bam;            // 模块5: SAM/BAM处理
mod bed;                // 模块6: BED/GTF/VCF处理
mod genomic_intervals;  // 模块7: 基因组区间操作
mod rnaseq;             // 模块8: RNA-seq分析
mod single_cell;        // 模块9: 单细胞测序
mod epigenomics;        // 模块10: 表观组学
mod metagenomics;       // 模块11: 宏基因组
mod proteomics;         // 模块12: 蛋白质组学
mod statistics;         // 模块13: 统计计算
mod file_format;        // 模块14: 文件格式转换
mod parallel;           // 模块15: 并行批处理

// ============================================================
// Python模块定义
// ============================================================
// #[pymodule]宏将这个函数标记为Python模块的入口
// m: &Bound<'_, PyModule>是PyO3 0.22的新API
// 使用wrap_pyfunction!宏将Rust函数包装为Python可调用对象

#[pymodule]
fn my_python_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // ============================================================
    // 模块1: 序列处理 (sequence)
    // 提供DNA/RNA序列的基础处理功能
    // 包括：GC含量、反向互补、碱基统计、转录翻译、ORF查找等
    // ============================================================
    m.add_function(wrap_pyfunction!(sequence::gc_content, m)?)?;           // GC含量计算
    m.add_function(wrap_pyfunction!(sequence::reverse_complement, m)?)?;  // 反向互补
    m.add_function(wrap_pyfunction!(sequence::count_bases, m)?)?;         // 碱基计数
    m.add_function(wrap_pyfunction!(sequence::seq_length, m)?)?;          // 序列长度
    m.add_function(wrap_pyfunction!(sequence::is_valid_dna, m)?)?;        // DNA验证
    m.add_function(wrap_pyfunction!(sequence::is_valid_rna, m)?)?;        // RNA验证
    m.add_function(wrap_pyfunction!(sequence::transcribe, m)?)?;          // DNA转RNA
    m.add_function(wrap_pyfunction!(sequence::reverse_transcribe, m)?)?;  // RNA转DNA
    m.add_function(wrap_pyfunction!(sequence::count_n, m)?)?;             // N碱基计数
    m.add_function(wrap_pyfunction!(sequence::normalize_seq, m)?)?;       // 大小写标准化
    m.add_function(wrap_pyfunction!(sequence::seq_slice, m)?)?;           // 序列切片
    m.add_function(wrap_pyfunction!(sequence::concat_seqs, m)?)?;         // 序列拼接
    m.add_function(wrap_pyfunction!(sequence::compression_ratio, m)?)?;   // 压缩比
    m.add_function(wrap_pyfunction!(sequence::is_low_complexity, m)?)?;   // 低复杂度检测
    m.add_function(wrap_pyfunction!(sequence::find_orfs, m)?)?;           // ORF查找
    m.add_function(wrap_pyfunction!(sequence::translate, m)?)?;           // 密码子翻译
    m.add_function(wrap_pyfunction!(sequence::seq_stats, m)?)?;           // 序列统计

    // ============================================================
    // 模块2: FASTQ质控 (fastq)
    // 提供FASTQ格式测序数据的质量控制功能
    // 包括：质量统计、过滤、分布分析等
    // ============================================================
    m.add_function(wrap_pyfunction!(fastq::fastq_stats, m)?)?;            // FASTQ统计
    m.add_function(wrap_pyfunction!(fastq::fastq_filter, m)?)?;           // FASTQ过滤
    m.add_function(wrap_pyfunction!(fastq::per_base_quality, m)?)?;       // 位点质量
    m.add_function(wrap_pyfunction!(fastq::length_distribution, m)?)?;    // 长度分布
    m.add_function(wrap_pyfunction!(fastq::gc_distribution, m)?)?;        // GC分布
    m.add_function(wrap_pyfunction!(fastq::sliding_window_filter, m)?)?;  // 滑动窗口过滤

    // ============================================================
    // 模块3: k-mer与sketch算法 (kmer)
    // 提供k-mer计数、MinHash草图、相似度计算等功能
    // ============================================================
    m.add_function(wrap_pyfunction!(kmer::count_kmers, m)?)?;             // k-mer计数
    m.add_function(wrap_pyfunction!(kmer::canonical_kmer, m)?)?;          // canonical k-mer
    m.add_function(wrap_pyfunction!(kmer::minimizers, m)?)?;              // minimizer提取
    m.add_function(wrap_pyfunction!(kmer::minhash_sketch, m)?)?;          // MinHash草图
    m.add_function(wrap_pyfunction!(kmer::jaccard_similarity, m)?)?;      // Jaccard相似度
    m.add_function(wrap_pyfunction!(kmer::containment_index, m)?)?;       // 包含指数
    m.add_function(wrap_pyfunction!(kmer::mash_distance, m)?)?;           // Mash距离
    m.add_function(wrap_pyfunction!(kmer::kmer_spectrum, m)?)?;           // k-mer频谱
    m.add_function(wrap_pyfunction!(kmer::syncmers, m)?)?;                // syncmer提取

    // ============================================================
    // 模块4: 比对算法 (alignment)
    // 提供序列比对相关的算法实现
    // 包括：编辑距离、全局/局部比对、CIGAR解析等
    // ============================================================
    m.add_function(wrap_pyfunction!(alignment::hamming_distance, m)?)?;   // Hamming距离
    m.add_function(wrap_pyfunction!(alignment::levenshtein_distance, m)?)?; // Levenshtein距离
    m.add_function(wrap_pyfunction!(alignment::needleman_wunsch, m)?)?;  // NW全局比对
    m.add_function(wrap_pyfunction!(alignment::smith_waterman, m)?)?;    // SW局部比对
    m.add_function(wrap_pyfunction!(alignment::banded_alignment, m)?)?;  // 带状比对
    m.add_function(wrap_pyfunction!(alignment::parse_cigar, m)?)?;       // CIGAR解析
    m.add_function(wrap_pyfunction!(alignment::cigar_stats, m)?)?;       // CIGAR统计
    m.add_function(wrap_pyfunction!(alignment::primer_match, m)?)?;      // 引物匹配

    // ============================================================
    // 模块5: SAM/BAM处理 (sam_bam)
    // 提供SAM格式比对结果的解析和统计功能
    // ============================================================
    m.add_function(wrap_pyfunction!(sam_bam::parse_flag, m)?)?;           // FLAG解析
    m.add_function(wrap_pyfunction!(sam_bam::sam_stats, m)?)?;            // SAM统计
    m.add_function(wrap_pyfunction!(sam_bam::filter_by_mapq, m)?)?;      // MAPQ过滤
    m.add_function(wrap_pyfunction!(sam_bam::fetch_region, m)?)?;        // 区间提取
    m.add_function(wrap_pyfunction!(sam_bam::coverage_at_position, m)?)?; // 位点覆盖
    m.add_function(wrap_pyfunction!(sam_bam::region_coverage, m)?)?;      // 区间覆盖
    m.add_function(wrap_pyfunction!(sam_bam::insert_size_stats, m)?)?;   // 插入片段统计

    // ============================================================
    // 模块6: BED/GTF/VCF处理 (bed)
    // 提供基因组区间文件的解析和操作功能
    // ============================================================
    m.add_function(wrap_pyfunction!(bed::bed_stats, m)?)?;                // BED统计
    m.add_function(wrap_pyfunction!(bed::bed_intersect, m)?)?;            // 区间交集
    m.add_function(wrap_pyfunction!(bed::bed_merge, m)?)?;                // 区间合并
    m.add_function(wrap_pyfunction!(bed::bed_subtract, m)?)?;             // 区间差集
    m.add_function(wrap_pyfunction!(bed::bed_closest, m)?)?;              // 最近特征
    m.add_function(wrap_pyfunction!(bed::vcf_stats, m)?)?;                // VCF统计
    m.add_function(wrap_pyfunction!(bed::vcf_filter, m)?)?;               // VCF过滤
    m.add_function(wrap_pyfunction!(bed::gtf_stats, m)?)?;                // GTF统计

    // ============================================================
    // 模块7: 基因组区间算法 (genomic_intervals)
    // 提供基因组区间的高级操作功能
    // ============================================================
    m.add_function(wrap_pyfunction!(genomic_intervals::intervals_overlap, m)?)?; // 重叠检测
    m.add_function(wrap_pyfunction!(genomic_intervals::overlap_length, m)?)?;    // 重叠长度
    m.add_function(wrap_pyfunction!(genomic_intervals::batch_overlap, m)?)?;     // 批量重叠
    m.add_function(wrap_pyfunction!(genomic_intervals::sliding_bins, m)?)?;       // 滑动分箱
    m.add_function(wrap_pyfunction!(genomic_intervals::window_count, m)?)?;       // 窗口计数
    m.add_function(wrap_pyfunction!(genomic_intervals::tss_distance, m)?)?;       // TSS距离
    m.add_function(wrap_pyfunction!(genomic_intervals::interval_length, m)?)?;    // 区间长度
    m.add_function(wrap_pyfunction!(genomic_intervals::interval_tree_query, m)?)?; // 区间树查询
    m.add_function(wrap_pyfunction!(genomic_intervals::filter_blacklist, m)?)?;   // 黑名单过滤
    m.add_function(wrap_pyfunction!(genomic_intervals::nearest_feature, m)?)?;    // 最近特征
    m.add_function(wrap_pyfunction!(genomic_intervals::match_enhancer_promoter, m)?)?; // 增强子匹配
    m.add_function(wrap_pyfunction!(genomic_intervals::coverage_over_intervals, m)?)?; // 区间覆盖

    // ============================================================
    // 模块8: RNA-seq分析 (rnaseq)
    // 提供RNA-seq数据处理和分析功能
    // ============================================================
    m.add_function(wrap_pyfunction!(rnaseq::gene_count_matrix, m)?)?;     // 基因计数矩阵
    m.add_function(wrap_pyfunction!(rnaseq::exon_count, m)?)?;            // 外显子计数
    m.add_function(wrap_pyfunction!(rnaseq::splice_junction_stats, m)?)?; // 剪接位点
    m.add_function(wrap_pyfunction!(rnaseq::detect_intron_retention, m)?)?; // 内含子保留
    m.add_function(wrap_pyfunction!(rnaseq::umi_collapse, m)?)?;          // UMI去重
    m.add_function(wrap_pyfunction!(rnaseq::gene_biotype_stats, m)?)?;    // 基因类型统计
    m.add_function(wrap_pyfunction!(rnaseq::transcript_length_stats, m)?)?; // 转录本长度

    // ============================================================
    // 模块9: 单细胞测序 (single_cell)
    // 提供单细胞RNA-seq数据处理功能
    // ============================================================
    m.add_function(wrap_pyfunction!(single_cell::extract_barcode, m)?)?;  // barcode提取
    m.add_function(wrap_pyfunction!(single_cell::extract_umi, m)?)?;      // UMI提取
    m.add_function(wrap_pyfunction!(single_cell::correct_barcode, m)?)?;  // barcode纠错
    m.add_function(wrap_pyfunction!(single_cell::barcode_rescue, m)?)?;   // barcode rescue
    m.add_function(wrap_pyfunction!(single_cell::umi_dedup, m)?)?;        // UMI去重
    m.add_function(wrap_pyfunction!(single_cell::build_feature_barcode_matrix, m)?)?; // 矩阵构建
    m.add_function(wrap_pyfunction!(single_cell::to_matrix_market, m)?)?; // MTX格式
    m.add_function(wrap_pyfunction!(single_cell::cell_qc, m)?)?;          // 细胞QC
    m.add_function(wrap_pyfunction!(single_cell::demux_tags, m)?)?;       // 标签去重
    m.add_function(wrap_pyfunction!(single_cell::count_guides, m)?)?;     // guide计数

    // ============================================================
    // 模块10: 表观组学 (epigenomics)
    // 提供表观基因组学数据分析功能
    // ============================================================
    m.add_function(wrap_pyfunction!(epigenomics::tn5_insertion_sites, m)?)?; // Tn5插入位点
    m.add_function(wrap_pyfunction!(epigenomics::count_reads_in_peaks, m)?)?; // peak read计数
    m.add_function(wrap_pyfunction!(epigenomics::fragment_length_distribution, m)?)?; // 片段长度
    m.add_function(wrap_pyfunction!(epigenomics::classify_nucleosome_free, m)?)?; // NFR分类
    m.add_function(wrap_pyfunction!(epigenomics::chipseq_peak_coverage, m)?)?; // ChIP-seq覆盖
    m.add_function(wrap_pyfunction!(epigenomics::cpg_methylation_levels, m)?)?; // CpG甲基化
    m.add_function(wrap_pyfunction!(epigenomics::fragment_overlap_peaks, m)?)?; // 片段重叠
    m.add_function(wrap_pyfunction!(epigenomics::coverage_over_bins, m)?)?;     // 分箱覆盖
    m.add_function(wrap_pyfunction!(epigenomics::cuttag_fragment_stats, m)?)?;  // CUT&Tag统计

    // ============================================================
    // 模块11: 宏基因组 (metagenomics)
    // 提供宏基因组学数据分析功能
    // ============================================================
    m.add_function(wrap_pyfunction!(metagenomics::kmer_classify, m)?)?;    // k-mer分类
    m.add_function(wrap_pyfunction!(metagenomics::ani_approximate, m)?)?; // ANI近似
    m.add_function(wrap_pyfunction!(metagenomics::genome_sketch, m)?)?;   // 基因组草图
    m.add_function(wrap_pyfunction!(metagenomics::preprocess_16s, m)?)?;  // 16S预处理
    m.add_function(wrap_pyfunction!(metagenomics::marker_gene_match, m)?)?; // marker匹配
    m.add_function(wrap_pyfunction!(metagenomics::contamination_screen, m)?)?; // 污染筛查
    m.add_function(wrap_pyfunction!(metagenomics::greedy_cluster, m)?)?;  // 贪心聚类
    m.add_function(wrap_pyfunction!(metagenomics::estimate_abundance, m)?)?; // 丰度估计

    // ============================================================
    // 模块12: 蛋白质组学 (proteomics)
    // 提供蛋白质组学数据分析功能
    // ============================================================
    m.add_function(wrap_pyfunction!(proteomics::trypsin_digest, m)?)?;    // 胰蛋白酶酶切
    m.add_function(wrap_pyfunction!(proteomics::peptide_mass, m)?)?;      // 肽段质量
    m.add_function(wrap_pyfunction!(proteomics::enumerate_missed_cleavages, m)?)?; // missed cleavage
    m.add_function(wrap_pyfunction!(proteomics::is_unique_peptide, m)?)?; // 肽段唯一性
    m.add_function(wrap_pyfunction!(proteomics::generate_decoy, m)?)?;    // decoy数据库
    m.add_function(wrap_pyfunction!(proteomics::enumerate_modification_sites, m)?)?; // 修饰位点
    m.add_function(wrap_pyfunction!(proteomics::peptide_mz, m)?)?;        // 肽段m/z
    m.add_function(wrap_pyfunction!(proteomics::parse_fasta, m)?)?;       // FASTA解析
    m.add_function(wrap_pyfunction!(proteomics::protein_molecular_weight, m)?)?; // 分子量
    m.add_function(wrap_pyfunction!(proteomics::peptide_pi, m)?)?;        // 等电点

    // ============================================================
    // 模块13: 统计计算 (statistics)
    // 提供生物信息学常用统计分析功能
    // ============================================================
    m.add_function(wrap_pyfunction!(statistics::permutation_test, m)?)?;   // 置换检验
    m.add_function(wrap_pyfunction!(statistics::bootstrap_ci, m)?)?;       // Bootstrap CI
    m.add_function(wrap_pyfunction!(statistics::hypergeometric_test, m)?)?; // 超几何检验
    m.add_function(wrap_pyfunction!(statistics::fisher_exact, m)?)?;       // Fisher精确检验
    m.add_function(wrap_pyfunction!(statistics::bh_correction, m)?)?;      // BH FDR校正
    m.add_function(wrap_pyfunction!(statistics::rolling_correlation, m)?)?; // 滚动相关
    m.add_function(wrap_pyfunction!(statistics::pearson_correlation, m)?)?; // Pearson相关
    m.add_function(wrap_pyfunction!(statistics::mann_whitney_u, m)?)?;     // Mann-Whitney U
    m.add_function(wrap_pyfunction!(statistics::descriptive_stats, m)?)?;  // 描述统计

    // ============================================================
    // 模块14: 文件格式转换 (file_format)
    // 提供生物信息学文件格式的转换功能
    // ============================================================
    m.add_function(wrap_pyfunction!(file_format::fastq_to_fasta, m)?)?;   // FASTQ→FASTA
    m.add_function(wrap_pyfunction!(file_format::sam_to_tsv, m)?)?;       // SAM→TSV
    m.add_function(wrap_pyfunction!(file_format::gtf_to_bed, m)?)?;       // GTF→BED
    m.add_function(wrap_pyfunction!(file_format::vcf_to_table, m)?)?;     // VCF→表格
    m.add_function(wrap_pyfunction!(file_format::count_matrix_to_mtx, m)?)?; // 矩阵→MTX
    m.add_function(wrap_pyfunction!(file_format::merge_count_matrices, m)?)?; // 矩阵合并
    m.add_function(wrap_pyfunction!(file_format::parse_expression_matrix, m)?)?; // 表达矩阵
    m.add_function(wrap_pyfunction!(file_format::write_gz, m)?)?;         // 压缩写入
    m.add_function(wrap_pyfunction!(file_format::read_gz, m)?)?;          // 压缩读取

    // ============================================================
    // 模块15: 并行批处理 (parallel)
    // 提供批量处理和并行计算的辅助功能
    // ============================================================
    m.add_function(wrap_pyfunction!(parallel::parallel_gc_content, m)?)?;       // 并行GC
    m.add_function(wrap_pyfunction!(parallel::parallel_count_kmers, m)?)?;      // 并行k-mer
    m.add_function(wrap_pyfunction!(parallel::parallel_reverse_complement, m)?)?; // 并行反向互补
    m.add_function(wrap_pyfunction!(parallel::parallel_hamming_matrix, m)?)?;   // Hamming矩阵
    m.add_function(wrap_pyfunction!(parallel::parallel_levenshtein_matrix, m)?)?; // Levenshtein矩阵
    m.add_function(wrap_pyfunction!(parallel::parallel_seq_stats, m)?)?;        // 并行统计
    m.add_function(wrap_pyfunction!(parallel::parallel_find_orfs, m)?)?;        // 并行ORF
    m.add_function(wrap_pyfunction!(parallel::parallel_parse_fasta, m)?)?;      // 并行FASTA
    m.add_function(wrap_pyfunction!(parallel::parallel_bed_intersect, m)?)?;    // 并行BED
    m.add_function(wrap_pyfunction!(parallel::parallel_jaccard_matrix, m)?)?;   // Jaccard矩阵
    m.add_function(wrap_pyfunction!(parallel::batch_fastq_filter, m)?)?;        // 批量过滤
    m.add_function(wrap_pyfunction!(parallel::parallel_descriptive_stats, m)?)?; // 并行描述统计

    Ok(())
}
