import sys
import io
import os
import tempfile
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import my_python_module as bio

# ============================================================
# 生信功能测试 — 示意数据
# ============================================================

print("=" * 70)
print("  Rust + PyO3 生信功能库 全模块测试")
print("=" * 70)

# ── 1. 序列处理 ─────────────────────────────────────────────
print("\n" + "─" * 50)
print("1. 序列处理 (sequence)")
print("─" * 50)

seq = "ATGCGCATATGCGCGATATAGCTAGCTAGCTAGCTAGCN"
print(f"  序列: {seq}")
print(f"  GC含量: {bio.gc_content(seq):.4f}")
print(f"  反向互补: {bio.reverse_complement(seq)}")
print(f"  碱基统计: {bio.count_bases(seq)}")
print(f"  序列长度: {bio.seq_length(seq)}")
print(f"  是否合法DNA: {bio.is_valid_dna(seq)}")
print(f"  是否合法RNA: {bio.is_valid_rna('AUGCGCAU')}")
print(f"  DNA转RNA: {bio.transcribe('ATGCGC')}")
print(f"  RNA转DNA: {bio.reverse_transcribe('AUGCGC')}")
print(f"  N碱基数: {bio.count_n(seq)}")
print(f"  标准化序列: {bio.normalize_seq('atgc')}")
print(f"  序列切片[0:6]: {bio.seq_slice(seq, 0, 6)}")
print(f"  拼接序列: {bio.concat_seqs(['ATG', 'CGC', 'ATA'], '---')}")
print(f"  信息压缩比: {bio.compression_ratio(seq):.4f}")
print(f"  低复杂度检测: {bio.is_low_complexity('AAAAAAACCCCCC', 0.6)}")

orfs = bio.find_orfs("ATGAAATTTAAACCCGGGATGCCC")
print(f"  ORF查找: {orfs}")
print(f"  翻译蛋白质: {bio.translate('ATGAAATTTAAACCC')}")

stats = bio.seq_stats("ATGCGCGATATAGCTAGC")
print(f"  综合统计: {stats}")

# ── 2. FASTQ质控 ─────────────────────────────────────────────
print("\n" + "─" * 50)
print("2. FASTQ质控 (fastq)")
print("─" * 50)

# 创建示例FASTQ文件
tmpdir = tempfile.mkdtemp()
fastq_path = os.path.join(tmpdir, "sample.fq")
with open(fastq_path, "w") as f:
    for i in range(100):
        f.write(f"@read_{i}\n")
        f.write("ATGCGCATATGCGCGATATAGCTAGCTAGCTAGCTAGCN\n")
        f.write("+\n")
        f.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!\n")

fqc = bio.fastq_stats(fastq_path)
print(f"  总reads数: {fqc['total_reads']}")
print(f"  总碱基数: {fqc['total_bases']}")
print(f"  平均长度: {fqc['avg_length']:.1f}")
print(f"  Q20率: {fqc['q20_rate']:.4f}")
print(f"  Q30率: {fqc['q30_rate']:.4f}")
print(f"  GC含量: {fqc['gc_content']:.4f}")

# 过滤
filtered_path = os.path.join(tmpdir, "filtered.fq")
kept = bio.fastq_filter(fastq_path, filtered_path, min_len=10, min_qual=20.0)
print(f"  过滤后保留reads: {kept}")

# sliding window
total, passed = bio.sliding_window_filter(fastq_path, window_size=5, min_avg_qual=20.0)
print(f"  滑窗过滤: {passed}/{total} 通过")

# ── 3. k-mer与sketch ────────────────────────────────────────
print("\n" + "─" * 50)
print("3. k-mer与sketch算法 (kmer)")
print("─" * 50)

test_seq = "ATGCGCATGCGCATGCGC"
kmers = bio.count_kmers(test_seq, 3)
print(f"  k-mer计数(k=3): {dict(kmers)}")

canonical = bio.canonical_kmer(test_seq, 3)
print(f"  Canonical k-mer数量: {len(canonical)}")

mins = bio.minimizers(test_seq, 3, 4)
print(f"  Minimizers: {mins}")

sketch = bio.minhash_sketch(test_seq, 3, 5)
print(f"  MinHash sketch: {sketch}")

jaccard = bio.jaccard_similarity("ATGCGCATG", "ATGCGAATG", 3)
print(f"  Jaccard相似度: {jaccard:.4f}")

contain = bio.containment_index("ATGCGCATG", "ATGCGAATG", 3)
print(f"  Containment index: {contain:.4f}")

mash = bio.mash_distance("ATGCGCATGCGC", "ATGCGAATGCGC", 3, 10)
print(f"  Mash distance: {mash:.4f}")

spectrum = bio.kmer_spectrum("ATGCGCATGCGC", 3)
print(f"  k-mer频谱: {dict(spectrum)}")

# ── 4. 比对算法 ──────────────────────────────────────────────
print("\n" + "─" * 50)
print("4. 比对算法 (alignment)")
print("─" * 50)

ham = bio.hamming_distance("ATGCGC", "ATGCCC")
print(f"  Hamming距离: {ham}")

lev = bio.levenshtein_distance("kitten", "sitting")
print(f"  Levenshtein距离: {lev}")

nw_a, nw_b, nw_score = bio.needleman_wunsch("ATGC", "ATCC", 2, -1, -2)
print(f"  NW比对: {nw_a}")
print(f"          {nw_b}")
print(f"  NW得分: {nw_score}")

sw_a, sw_b, sw_score = bio.smith_waterman("AAATGC", "ATGC", 2, -1, -2)
print(f"  SW比对: {sw_a}")
print(f"          {sw_b}")
print(f"  SW得分: {sw_score}")

banded = bio.banded_alignment("ATGCGC", "ATGCGC", 3, 2, -1, -2)
print(f"  带状比对得分: {banded}")

cigar = "10M2I5M3D8M"
cigar_ops = bio.parse_cigar(cigar)
print(f"  CIGAR解析: {cigar_ops}")
cigar_s = bio.cigar_stats(cigar)
print(f"  CIGAR统计: {cigar_s}")

primers = bio.primer_match("AAATGCGCGATATAGCTAGC", "ATGCGC", 1)
print(f"  引物匹配: {primers}")

# ── 5. SAM/BAM处理 ──────────────────────────────────────────
print("\n" + "─" * 50)
print("5. SAM/BAM处理 (sam_bam)")
print("─" * 50)

sam_path = os.path.join(tmpdir, "sample.sam")
with open(sam_path, "w") as f:
    f.write("@HD\tVN:1.6\tSO:coordinate\n")
    f.write("read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tATGCGCATATGCGCGATATAGCTAGCTAGCTAGCTAGCNATGCGCATAT\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!\n")
    f.write("read2\t16\tchr1\t200\t30\t30M\t*\t0\t0\tATGCGCATATGCGCGATATAGCTAGCTAG\tIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    f.write("read3\t4\t*\t0\t0\t*\t*\t0\t0\tATGCGC\tIIIIII\n")
    f.write("read4\t0\tchr1\t150\t60\t40M\tchr2\t300\t100\tATGCGCATATGCGCGATATAGCTAGCTAGCTAGCTAGCN\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")

flag_info = bio.parse_flag(0)
print(f"  Flag=0 解析: read_paired={flag_info['read_paired']}, unmapped={flag_info['read_unmapped']}")

flag_info2 = bio.parse_flag(16)
print(f"  Flag=16 解析: reverse={flag_info2['read_reverse']}")

sam_stats = bio.sam_stats(sam_path)
print(f"  SAM统计: total={sam_stats['total_reads']}, mapped={sam_stats['mapped']}, unmapped={sam_stats['unmapped']}")

insert_stats = bio.insert_size_stats(sam_path)
print(f"  Insert size: {insert_stats}")

# ── 6. BED/GTF/VCF处理 ─────────────────────────────────────
print("\n" + "─" * 50)
print("6. BED/GTF/VCF处理 (bed)")
print("─" * 50)

bed_a = os.path.join(tmpdir, "a.bed")
bed_b = os.path.join(tmpdir, "b.bed")
with open(bed_a, "w") as f:
    f.write("chr1\t100\t500\tpeak1\t1000\t+\n")
    f.write("chr1\t800\t1200\tpeak2\t900\t+\n")
    f.write("chr2\t200\t600\tpeak3\t800\t-\n")
with open(bed_b, "w") as f:
    f.write("chr1\t300\t700\tgene1\t500\t+\n")
    f.write("chr1\t1000\t1500\tgene2\t600\t-\n")
    f.write("chr2\t100\t300\tgene3\t400\t+\n")

bed_info = bio.bed_stats(bed_a)
print(f"  BED A统计: total={bed_info['total_intervals']}, bases={bed_info['total_bases']}")

intersect = bio.bed_intersect(bed_a, bed_b)
print(f"  区间交集: {intersect}")

merged = bio.bed_merge(bed_a)
print(f"  区间合并: {merged}")

subtracted = bio.bed_subtract(bed_a, bed_b)
print(f"  区间差集: {subtracted}")

closest = bio.bed_closest(bed_a, bed_b)
print(f"  最近特征: {closest}")

# VCF
vcf_path = os.path.join(tmpdir, "sample.vcf")
with open(vcf_path, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t100\trs1\tA\tG\t30\tPASS\tDP=50\n")
    f.write("chr1\t200\trs2\tAT\tA\t25\tPASS\tDP=30\n")
    f.write("chr2\t300\trs3\tG\tC\t15\tLowQual\tDP=10\n")
    f.write("chr2\t500\trs4\tC\tT\t40\tPASS\tDP=100\n")

vcf_info = bio.vcf_stats(vcf_path)
print(f"  VCF统计: total={vcf_info['total_variants']}, SNPs={vcf_info['snps']}, indels={vcf_info['indels']}")

filtered_vcf = bio.vcf_filter(vcf_path, 25.0)
print(f"  VCF过滤(QUAL>=25): {filtered_vcf}")

# GTF
gtf_path = os.path.join(tmpdir, "sample.gtf")
with open(gtf_path, "w") as f:
    f.write('#!genome-build GRCh38\n')
    f.write('chr1\tensembl\tgene\t100\t500\t.\t+\t.\tgene_id "GENE1"; gene_biotype "protein_coding";\n')
    f.write('chr1\tensembl\ttranscript\t100\t500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";\n')
    f.write('chr1\tensembl\texon\t100\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_id "EX1";\n')
    f.write('chr1\tensembl\texon\t300\t500\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_id "EX2";\n')
    f.write('chr2\tensembl\tgene\t1000\t2000\t.\t-\t.\tgene_id "GENE2"; gene_biotype "lncRNA";\n')

gtf_info = bio.gtf_stats(gtf_path)
print(f"  GTF统计: genes={gtf_info['gene_count']}, transcripts={gtf_info['transcript_count']}")

# ── 7. 基因组区间算法 ────────────────────────────────────────
print("\n" + "─" * 50)
print("7. 基因组区间算法 (genomic_intervals)")
print("─" * 50)

print(f"  区间重叠: {bio.intervals_overlap('chr1', 100, 500, 'chr1', 300, 700)}")
print(f"  重叠长度: {bio.overlap_length('chr1', 100, 500, 'chr1', 300, 700)}")

batch = bio.batch_overlap(
    [("chr1", 100, 500), ("chr1", 800, 1200)],
    [("chr1", 300, 700), ("chr1", 900, 1100)]
)
print(f"  批量重叠: {batch}")

bins = bio.sliding_bins("chr1", 0, 1000, 100, 50)
print(f"  滑动窗口bins数: {len(bins)}")

tss = bio.tss_distance("chr1", 200, "chr1", 100, 500, '+')
print(f"  TSS距离: {tss}")

nearest = bio.nearest_feature("chr1", 600, [("chr1", 100, 200, "gene1"), ("chr1", 500, 800, "gene2")])
print(f"  最近特征: {nearest}")

ep_match = bio.match_enhancer_promoter(
    [("chr1", 100, 200, "E1"), ("chr1", 600, 700, "E2")],
    [("chr1", 400, 500, "P1"), ("chr1", 800, 900, "P2")],
    300
)
print(f"  Enhancer-Promoter匹配: {ep_match}")

cov = bio.coverage_over_intervals(
    [("chr1", 100, 500), ("chr1", 800, 1200)],
    [("chr1", 150), ("chr1", 200), ("chr1", 900)]
)
print(f"  区间覆盖: {cov}")

# ── 8. RNA-seq ───────────────────────────────────────────────
print("\n" + "─" * 50)
print("8. RNA-seq分析 (rnaseq)")
print("─" * 50)

junctions = bio.splice_junction_stats("50M100N30M", 1000)
print(f"  剪接junction: {junctions}")

intron_ret = bio.detect_intron_retention("100M", 1000, 1050, 1150)
print(f"  Intron retention: {intron_ret}")

umis = bio.umi_collapse(["AAGG", "AACG", "TTCC", "TTCG", "GGGG"], 1)
print(f"  UMI collapse: {umis}")

# ── 9. 单细胞 ───────────────────────────────────────────────
print("\n" + "─" * 50)
print("9. 单细胞测序 (single_cell)")
print("─" * 50)

barcode = bio.extract_barcode("AAGCTAGCTAGCTAGCTAGCTAGCTAGC", 16, 0)
print(f"  提取barcode: {barcode}")

umi = bio.extract_umi("AAGCTAGCTAGCTAGCTAGC", 10, 16)
print(f"  提取UMI: {umi}")

corrected = bio.correct_barcode("AAGG", ["AAGG", "TTCC", "GGCC"], 1)
print(f"  Barcode校正: {corrected}")

rescue = bio.barcode_rescue("AACC", ["AAGG", "TTCC", "AACC"], 2)
print(f"  Barcode rescue: {rescue}")

mm = bio.to_matrix_market(3, 3, [(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0)])
print(f"  Matrix Market格式:\n{mm}")

qc = bio.cell_qc([100, 200, 300, 400], [10, 20, 30, 40])
print(f"  Cell QC: {qc}")

tags = bio.demux_tags(
    [[100.0, 5.0, 3.0], [2.0, 200.0, 1.0], [1.0, 3.0, 150.0]],
    ["HTO1", "HTO2", "HTO3"],
    10.0
)
print(f"  Tag分配: {tags}")

guides = bio.count_guides(["g1", "g1", "g2", "g1", "g3", "g2"])
print(f"  Guide计数: {guides}")

# ── 10. 表观组学 ─────────────────────────────────────────────
print("\n" + "─" * 50)
print("10. 表观组学 (epigenomics)")
print("─" * 50)

tn5 = bio.tn5_insertion_sites("chr1", 1000, False)
print(f"  Tn5插入位点(+): {tn5}")
tn5_r = bio.tn5_insertion_sites("chr1", 1000, True)
print(f"  Tn5插入位点(-): {tn5_r}")

peak_counts = bio.count_reads_in_peaks(
    [("chr1", 150), ("chr1", 250), ("chr1", 500)],
    [("chr1", 100, 300), ("chr1", 200, 400)]
)
print(f"  Peak reads计数: {dict(peak_counts)}")

frag_dist = bio.fragment_length_distribution([100, 150, 200, 100, 150, 150, 300])
print(f"  片段长度分布: {dict(frag_dist)}")

nfr, nuc = bio.classify_nucleosome_free([100, 150, 200, 50, 80, 250, 300], 150)
print(f"  Nucleosome-free: {nfr}, Nucleosomal: {nuc}")

frag_overlap = bio.fragment_overlap_peaks(
    [("chr1", 100, 300), ("chr1", 500, 700)],
    [("chr1", 200, 400), ("chr1", 600, 800)]
)
print(f"  片段-Peak重叠: {frag_overlap}")

bins_cov = bio.coverage_over_bins(
    [("chr1", 150), ("chr1", 250), ("chr1", 600)],
    [("chr1", 0, 200), ("chr1", 200, 400), ("chr1", 400, 600)]
)
print(f"  Bin覆盖度: {bins_cov}")

ct_stats = bio.cuttag_fragment_stats(
    [("chr1", 100, 300), ("chr1", 400, 600), ("chr1", 700, 800)],
    [("chr1", 50, 350), ("chr1", 380, 620)]
)
print(f"  CUT&Tag统计: {ct_stats}")

# ── 11. 宏基因组 ─────────────────────────────────────────────
print("\n" + "─" * 50)
print("11. 宏基因组 (metagenomics)")
print("─" * 50)

ani = bio.ani_approximate("ATGCGCATGCGCATGCGC", "ATGCGCATGCGAATGCGC", 3, 20)
print(f"  ANI近似: {ani:.4f}")

sketch = bio.genome_sketch("ATGCGCATGCGCATGCGCATGCGC", 3, 5)
print(f"  Genome sketch: {sketch}")

preprocessed = bio.preprocess_16s(
    "AAATGCGCGATATAGCTAGCCCCCCCCCCCTAGCTAGCTAGC",
    "ATGCGC", "CTAGCT", 1
)
print(f"  16S预处理: {preprocessed}")

markers = bio.marker_gene_match(
    "ATGCGCATGCGC",
    [("gene1", "ATGCGCATGCGC"), ("gene2", "TTTTAAAACCC"), ("gene3", "ATGCGAATGCGC")],
    3, 0.3
)
print(f"  Marker gene匹配: {markers}")

clusters = bio.greedy_cluster(
    ["ATGCGC", "ATGCGA", "TTTAAA", "TTTAAAG", "ATGCGT"],
    0.5, 3
)
print(f"  OTU聚类结果: {[(c[0][:10], len(c[1])) for c in clusters]}")

# ── 12. 蛋白质组学 ───────────────────────────────────────────
print("\n" + "─" * 50)
print("12. 蛋白质组学 (proteomics)")
print("─" * 50)

peptides = bio.trypsin_digest("MPEPTIDEKFRANOTHER", 2)
print(f"  Trypsin酶切: {peptides}")

mass = bio.peptide_mass("MPEPTIDE")
print(f"  肽段质量: {mass:.4f}")

mz = bio.peptide_mz("MPEPTIDE", 2)
print(f"  肽段m/z(+2): {mz:.4f}")

decoys = bio.generate_decoy(["MPEPTIDE", "KFRANOTHER"])
print(f"  Decoy序列: {decoys}")

mod_sites = bio.enumerate_modification_sites("MPEPTIDECM", ['M', 'C'])
print(f"  修饰位点: {mod_sites}")
print(f"  分子量(Human serum albumin理论): {bio.protein_molecular_weight('MKWVTFISLLFLFSSAYS'):,.2f} Da")
print(f"  等电点估算: {bio.peptide_pi('MPEPTIDE'):.2f}")

# ── 13. 统计计算 ─────────────────────────────────────────────
print("\n" + "─" * 50)
print("13. 统计计算 (statistics)")
print("─" * 50)

desc = bio.descriptive_stats([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
print(f"  描述统计: mean={desc['mean']:.2f}, median={desc['median']:.2f}, std={desc['std_dev']:.2f}")

perm = bio.permutation_test([1.0, 2.0, 3.0, 4.0, 5.0], [3.0, 4.0, 5.0, 6.0, 7.0], 1000)
print(f"  Permutation test: diff={perm['observed_diff']:.2f}, p={perm['p_value']:.4f}")

boot = bio.bootstrap_ci([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], 1000, 0.95)
print(f"  Bootstrap CI: [{boot['ci_lower']:.2f}, {boot['ci_upper']:.2f}]")

fdr = bio.bh_correction([0.01, 0.04, 0.03, 0.50, 0.80, 0.10])
print(f"  BH FDR校正: {[f'{x:.4f}' for x in fdr]}")

corr = bio.pearson_correlation([1.0, 2.0, 3.0, 4.0, 5.0], [2.0, 4.0, 5.0, 4.0, 5.0])
print(f"  Pearson相关: {corr:.4f}")

roll = bio.rolling_correlation([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [2.0, 3.0, 5.0, 4.0, 6.0, 7.0], 3)
print(f"  Rolling相关: {[f'{x:.4f}' for x in roll]}")

fisher = bio.fisher_exact(10, 5, 3, 12)
print(f"  Fisher exact p: {fisher:.6f}")

hgt = bio.hypergeometric_test(5, 20, 50, 1000)
print(f"  Hypergeometric p: {hgt:.6f}")

mw = bio.mann_whitney_u([1.0, 2.0, 3.0, 4.0, 5.0], [6.0, 7.0, 8.0, 9.0, 10.0])
print(f"  Mann-Whitney U: U={mw['u_statistic']:.2f}, p={mw['p_value']:.6f}")

# ── 14. 文件格式转换 ─────────────────────────────────────────
print("\n" + "─" * 50)
print("14. 文件格式转换 (file_format)")
print("─" * 50)

fasta_path = os.path.join(tmpdir, "output.fa")
n_converted = bio.fastq_to_fasta(fastq_path, fasta_path)
print(f"  FASTQ->FASTA: {n_converted} reads")

sam_tsv = os.path.join(tmpdir, "output.tsv")
n_sam = bio.sam_to_tsv(sam_path, sam_tsv)
print(f"  SAM->TSV: {n_sam} records")

gtf_bed = os.path.join(tmpdir, "genes.bed")
n_gtf = bio.gtf_to_bed(gtf_path, gtf_bed, "gene")
print(f"  GTF->BED(gene): {n_gtf} entries")

vcf_tsv = os.path.join(tmpdir, "variants.tsv")
n_vcf = bio.vcf_to_table(vcf_path, vcf_tsv)
print(f"  VCF->TSV: {n_vcf} variants")

mtx = bio.count_matrix_to_mtx([[1.0, 0.0, 3.0], [0.0, 2.0, 0.0]], os.path.join(tmpdir, "matrix.mtx"))
print(f"  Matrix->MTX: {mtx} non-zeros")

merged_m = bio.merge_count_matrices([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
print(f"  合并矩阵: {merged_m}")

# ── 15. 并行批处理 ───────────────────────────────────────────
print("\n" + "─" * 50)
print("15. 并行批处理 (parallel)")
print("─" * 50)

gc_list = bio.parallel_gc_content(["ATGCGC", "AAAAAA", "GCGCGC", "ATATAT"])
print(f"  并行GC含量: {[f'{x:.4f}' for x in gc_list]}")

total_kmers = bio.parallel_count_kmers(["ATGCGC", "ATGCGA"], 3)
print(f"  并行k-mer总计: {len(total_kmers)} 种")

rc_list = bio.parallel_reverse_complement(["ATGCGC", "AAAAA"])
print(f"  并行反向互补: {rc_list}")

ham_matrix = bio.parallel_hamming_matrix(["ATGC", "ATCC", "TTCC"])
print(f"  Hamming距离矩阵: {ham_matrix}")

lev_matrix = bio.parallel_levenshtein_matrix(["ATGC", "ATCC", "TTCC"])
print(f"  Levenshtein距离矩阵: {lev_matrix}")

jaccard_m = bio.parallel_jaccard_matrix(["ATGCGC", "ATGCGA", "TTTAAA"], 3)
print(f"  Jaccard相似度矩阵:")
for row in jaccard_m:
    print(f"    {[f'{x:.4f}' for x in row]}")

orfs_list = bio.parallel_find_orfs(["ATGAAATTTAAACCC", "ATGCCCTTTGGGAAATGA"])
print(f"  并行ORF查找: {orfs_list}")

# ── 清理临时文件 ─────────────────────────────────────────────
import shutil
shutil.rmtree(tmpdir, ignore_errors=True)

print("\n" + "=" * 70)
print("  全部测试完成!")
print("=" * 70)
