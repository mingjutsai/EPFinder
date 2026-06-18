# EPFinder

EPFinder is a machine-learning workflow that uses enhancer-promoter regulatory
features to prioritize target genes for noncoding SNPs. This repository provides
the first public pipeline version, from SNP preprocessing to model prediction.

## Repository contents

| Path | Purpose |
| --- | --- |
| `preprocessing/EPFinder_preprocessing.py` | Builds the 29-feature EPFinder input matrix from candidate SNPs and regulatory genomics files. |
| `preprocessing/config.yaml` | Public configuration template. Copy and edit paths before running. |
| `preprocessing/feature_list` | Required chromatin feature order for the trained model. |
| `scripts/EPFinder_predict.py` | Command-line EPFinder prediction using the trained PyCaret model. |
| `scripts/run_EPFinder_pipeline.py` | Convenience wrapper to run preprocessing and prediction together. |
| `finalize_EPFinder_model.pkl` | Trained EPFinder PyCaret model. |
| `dataset/gm12878_29features_ML.tsv` | Small benchmark-format matrix useful for testing prediction. |
| `docs/input_formats.md` | Input file requirements and expected formats. |

## Installation

Create the conda environment:

```bash
conda env create -f conda/EPFinder_env.yml
conda activate EPFinder_env
```

Register the environment as a Jupyter kernel only if you plan to run the
notebooks:

```bash
python -m ipykernel install --user --name EPFinder_env --display-name "EPFinder_env"
```

The preprocessing workflow also requires `bedtools` in `PATH` or at the path
set in `preprocessing/config.yaml`.

## Quick prediction test

Run EPFinder prediction on the included benchmark-format matrix:

```bash
python scripts/EPFinder_predict.py \
  --input dataset/gm12878_29features_ML.tsv \
  --model finalize_EPFinder_model \
  --output examples/gm12878_EPFinder_predictions.tsv
```

If the input contains `#Class`, the script also reports AUROC and AUPRC.

## Full preprocessing-to-prediction workflow

1. Copy the config template and edit paths for your server:

```bash
cp preprocessing/config.yaml preprocessing/config.local.yaml
```

2. Edit `preprocessing/config.local.yaml`.

Required inputs include:

- candidate SNP file in hg38 coordinates
- in-house Day 13 osteoblast Hi-C contact files
- hg38 transcript TSS annotation
- osteoblast transcript- and gene-level expression files
- hg38 bedGraph files for the chromatin features listed in `preprocessing/feature_list`

3. Run preprocessing:

```bash
python preprocessing/EPFinder_preprocessing.py preprocessing/config.local.yaml
```

The final preprocessing output is written to:

```text
{output_dir}/{output_file}
```

4. Run prediction:

```bash
python scripts/EPFinder_predict.py \
  --input /path/to/output/EPFinder_29features_ML.tsv \
  --model finalize_EPFinder_model \
  --output /path/to/output/EPFinder_predictions.tsv
```

Alternatively, run both steps with the wrapper:

```bash
python scripts/run_EPFinder_pipeline.py \
  --config preprocessing/config.local.yaml \
  --model finalize_EPFinder_model
```

## Output

`scripts/EPFinder_predict.py` appends an `EPFinder_score` column to the input
matrix. Higher scores indicate stronger model support for the SNP-promoter pair.

For GWAS applications, downstream ranking is typically performed per SNP or per
GWAS locus, depending on the biological question.

## Notes

- Keep the feature order in `preprocessing/feature_list` unchanged unless you
  retrain the model.
- All coordinates and regulatory files must use hg38.
- The public config uses placeholder paths. Use `*.local.yaml` files for server
  paths; these are ignored by Git.
- Large local analysis outputs are intentionally ignored to keep the GitHub
  repository focused on the reproducible EPFinder workflow.
