# EPFinder input formats

EPFinder expects all genomic coordinates to use hg38.

## Candidate SNP input

The preprocessing workflow reads a tab-delimited file with one candidate
variant per row. The first two columns are required:

| Column | Description |
| --- | --- |
| 1 | Chromosome, with or without `chr` prefix |
| 2 | 1-based SNP position |

Additional columns are preserved during preprocessing and can include SNP ID,
alleles, GWAS statistics or fine-mapping statistics.

## Hi-C contact files

Hi-C contacts are stored as one file per chromosome. The file path is assembled
from `hic_folder`, `hic_prefix` and chromosome name:

```text
{hic_folder}/{hic_prefix}chr1
{hic_folder}/{hic_prefix}chr2
...
```

Each file must contain three tab-delimited columns:

| Column | Description |
| --- | --- |
| 1 | Contact bin start |
| 2 | Contact bin start for the paired bin |
| 3 | KR-normalized contact value |

The default EPFinder model was built with 3 kb osteoblast Hi-C bins.

## TSS annotation

The TSS annotation file should be a sorted BED-like table compatible with
`bedtools intersect`. It must include transcript ID and gene symbol columns in
the positions expected by `preprocessing/EPFinder_preprocessing.py`.

## Expression matrices

The transcript expression file maps transcript IDs to expression values. The
gene expression file maps Ensembl gene IDs to expression values. The gene-list
file maps gene symbols to Ensembl gene IDs.

## Chromatin feature list

`preprocessing/feature_list` defines the feature order used by the trained
model. Keep this order unchanged unless you also retrain EPFinder.

Each row has:

```text
feature_name<TAB>path_to_hg38_bedGraph
```

The preprocessing workflow scores each feature at both the SNP-centered
enhancer window and the TSS-centered promoter window.

## Prediction matrix

The final preprocessing output contains metadata plus 29 model features:

- `HiC_Contact`
- `Tx_expression`
- `Gene_expression`
- 13 enhancer chromatin features
- 13 promoter chromatin features

The prediction CLI drops standard metadata columns (`#Class`, `ID`, `Enh`,
`Prom`, `TX`) when present and applies the trained PyCaret model to the
remaining feature columns.
