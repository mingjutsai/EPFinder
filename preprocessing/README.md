# EPFinder Preprocessing Pipeline

This directory contains the preprocessing pipeline that constructs the 29-feature input matrix required by EPFinder for enhancer-promoter interaction prediction. The pipeline integrates osteoblast multi-omics data from Hi-C, ChIP-seq, DNase-seq, and RNA-seq at candidate GWAS SNP loci.

## Overview

For each candidate GWAS SNP, the pipeline:

1. Identifies all genomic loci in 3D contact with the SNP using osteoblast Hi-C data
2. Intersects contacting loci with annotated transcription start sites (TSS)
3. Defines enhancer regions (centered on the SNP) and promoter regions (centered on each TSS)
4. Quantifies chromatin mark signals and expression levels at each enhancer-promoter pair
5. Assembles a unified feature matrix for EPFinder inference

The output is one row per candidate enhancer-promoter pair, with 29 features spanning Hi-C contact, histone modifications, chromatin accessibility, and gene expression.

## Contents

| File | Description |
|------|-------------|
| `EPFinder_preprocessing.py` | Main preprocessing script (9-step pipeline) |
| `config.yaml` | Configuration file — all input paths and parameters |
| `feature_list` | Ordered list of chromatin feature names and signal file paths |

## Requirements

- Python 3.8+
- bedtools (≥2.30): must be installed and accessible at the path specified in `config.yaml`
- Python packages: `pyyaml`, `pandas`, `concurrent.futures` (standard library)

Install Python dependencies:

```bash
pip install pyyaml pandas
```

## Configuration

Edit `config.yaml` before running. Key parameters:

```yaml
# Input GWAS SNP file (TSV: chr, pos, ...)
input_gwas: "/path/to/gwas_snps.tsv"

# Path to bedtools executable
bedtools_path: "/usr/bin/bedtools"

# Hi-C data: folder containing per-chromosome KR-normalized contact files
hic_folder: "/path/to/hic/3k"
hic_prefix: "OB13.bam.pre.filter.q30.3k.hic.KR."

# Genome annotations
tss_file: "/path/to/gencode.v29.annotation.gtf.pcg.tss_sorted_lexicographical"

# Expression data (ENCODE osteoblast RNA-seq ENCSR000CUF / ENCFF220VSM)
tx_expression: "/path/to/osteoblast_RNAseq_tx_ENCFF220VSM.tsv.simple"
gene_list:     "/path/to/gencode.v29.annotation.gtf_gene.txt"
gene_expression: "/path/to/osteoblast_RNAseq_gene_ENCFF220VSM.tsv.simple"

# Ordered list of chromatin feature signal files (Roadmap E129 imputed bigwigs)
feature_list: "/path/to/preprocessing/feature_list"

# Window sizes
enhancer_window: 250   # bp around SNP center
promoter_window: 250   # bp around TSS center
hic_bin_size: 3000     # Hi-C contact bin size (bp)

# Output
output_dir:  "/path/to/output"
output_file: "eBMD_SNPs_29features_ML.tsv"

# Parallelism for Step 1 (Hi-C mapping)
step1_nproc: 24
```

## Usage

```bash
python EPFinder_preprocessing.py config.yaml
```

The script creates a `tmp/` subdirectory inside `output_dir` for intermediate files. The final output file is written to `output_dir/output_file`.

## Pipeline Steps

| Step | Function | Description |
|------|----------|-------------|
| 1 | `step1_hic_matrix_mapping` | Map each GWAS SNP to its Hi-C contact bin; extract all interacting loci and contact frequencies. Chromosomes are processed in parallel (controlled by `step1_nproc`). |
| 2 | `step2_hic_prom` | Filter Hi-C contacts; remove entries with missing values. |
| 3 | `step3_hicbin2_tss_overlap` | Intersect Hi-C interacting regions with annotated TSS positions using bedtools. |
| 4 | `step4_format` | Format SNP-TSS pairs into a unified intermediate file; define enhancer and promoter window coordinates. |
| 5 | `step5_mapping_tx_quantifications` | Map transcript-level RNA-seq expression (FPKQ) to each candidate enhancer-promoter pair. |
| 6 | `step6_mapping_gene_quantification` | Map gene-level RNA-seq expression (FPKQ) to each candidate enhancer-promoter pair. |
| 7 | `step7_split_class_enh_prom` | Split enhancer and promoter region coordinates into separate BED files for chromatin feature scoring. |
| 8 | `step8_bedtool_overlap` | Compute cumulative chromatin signal intensity (signal × bp overlap) at each enhancer and promoter window using bedtools intersect. |
| 9 | `step9_add_allfeatures` | Merge all chromatin feature scores into the final feature matrix, following the order specified in `feature_list`. |

## Input Data Sources

### GWAS SNP input

The input TSV file should contain one SNP per row with at minimum: chromosome (col 1) and position (col 2). The pipeline was developed for eBMD GWAS candidate causal variants (COJO p < 5 × 10⁻⁸ or log₁₀BF > 3; n = 2,695 SNPs from Morris et al. 2019).

### Hi-C data

Per-chromosome KR-normalized contact files from Day 13 osteoblast Hi-C (in-house; OB13), processed with the Juicer pipeline at 3 kb resolution. Files are named: `{hic_prefix}{chr}` (e.g., `OB13.bam.pre.filter.q30.3k.hic.KR.chr1`).

### Chromatin features (`feature_list`)

Signal files are Roadmap Epigenomics E129 (osteoblast) imputed p-value bigwig tracks converted to bedGraph format at 25 bp non-sliding windows, except CTCF which uses ENCODE file ENCFF643JJS directly. The 13 features (each scored at enhancer and promoter = 26 values) are:

| Feature | Source | ENCODE / Roadmap ID |
|---------|--------|---------------------|
| H2A.Z | Roadmap E129 imputed | E129-H2A.Z.imputed.pval.signal.hg38.bigwig |
| H3K27ac | Roadmap E129 imputed | E129-H3K27ac.imputed.pval.signal.hg38.bigwig |
| H3K27me3 | Roadmap E129 imputed | E129-H3K27me3.imputed.pval.signal.hg38.bigwig |
| H3K36me3 | Roadmap E129 imputed | E129-H3K36me3.imputed.pval.signal.hg38.bigwig |
| H3K4me1 | Roadmap E129 imputed | E129-H3K4me1.imputed.pval.signal.hg38.bigwig |
| H3K4me2 | Roadmap E129 imputed | E129-H3K4me2.imputed.pval.signal.hg38.bigwig |
| H3K4me3 | Roadmap E129 imputed | E129-H3K4me3.imputed.pval.signal.hg38.bigwig |
| H3K79me2 | Roadmap E129 imputed | E129-H3K79me2.imputed.pval.signal.hg38.bigwig |
| H3K9ac | Roadmap E129 imputed | E129-H3K9ac.imputed.pval.signal.hg38.bigwig |
| H3K9me3 | Roadmap E129 imputed | E129-H3K9me3.imputed.pval.signal.hg38.bigwig |
| H4K20me1 | Roadmap E129 imputed | E129-H4K20me1.imputed.pval.signal.hg38.bigwig |
| DNase | Roadmap E129 imputed | E129-DNase.imputed.pval.signal.hg38.bigwig |
| CTCF | ENCODE direct | ENCFF643JJS |

### RNA-seq expression

ENCODE osteoblast RNA-seq (ENCSR000CUF, file ENCFF220VSM). Transcript- and gene-level FPKQ values processed by the ENCODE pipeline.

## Output

The final output is a TSV file with one candidate enhancer-promoter pair per row and 29 columns:

| Column(s) | Description |
|-----------|-------------|
| 1–12 | SNP and TSS metadata (chr, pos, gene ID, TSS position, etc.) |
| HiC_Contact | KR-normalized Hi-C contact frequency between SNP bin and TSS bin |
| Tx_expression | Transcript-level expression at the promoter gene (FPKQ) |
| Gene_expression | Gene-level expression at the promoter gene (FPKQ) |
| {Mark}_Enh | Chromatin signal at the enhancer window (13 marks) |
| {Mark}_Prom | Chromatin signal at the promoter window (13 marks) |

Total: 1 (Hi-C) + 2 (expression) + 13 × 2 (enhancer + promoter marks) = **29 features**

## Notes

- The pipeline assumes hg38 throughout. Hi-C files, bigwig/bedGraph signals, TSS annotations, and GWAS coordinates must all be in hg38.
- Intermediate files in `tmp/` are not automatically deleted and can be used for debugging individual steps.
- Step 1 is the most compute-intensive. Adjust `step1_nproc` to match available CPU cores.
- The `feature_list` file controls both the identity and order of chromatin features. The order must match the feature order used during EPFinder model training.
