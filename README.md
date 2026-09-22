# EPFinder

[![CI](https://github.com/mingjutsai/EPFinder/actions/workflows/ci.yml/badge.svg)](https://github.com/mingjutsai/EPFinder/actions/workflows/ci.yml)

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
| `tests/` | Test suite; runs on synthetic data with no external inputs. |
| `conda/EPFinder_env.yml` | Recommended runtime environment, including bedtools. |
| `conda/EPFinder_env.lock.yml` | Full frozen export for exact reproduction. |

## Installation

Create the conda environment. This installs bedtools as well as the Python
stack, so no separate bedtools install is needed:

```bash
conda env create -f conda/EPFinder_env.yml
conda activate EPFinder_env
```

Register the environment as a Jupyter kernel only if you plan to run the
notebooks:

```bash
python -m ipykernel install --user --name EPFinder_env --display-name "EPFinder_env"
```

### Which environment file to use

| File | Use it when |
| --- | --- |
| `conda/EPFinder_env.yml` | Normal use. Pins the model stack and installs bedtools, without PyCaret's optional extras. |
| `conda/EPFinder_env.lock.yml` | You need byte-level reproduction of the environment the release model was validated in. |
| `requirements.txt` | You are managing Python yourself. Covers the Python stack only; install bedtools separately. |

`conda/EPFinder_env.yml` is authoritative for normal use. It resolves to roughly
half the packages of the full export, because EPFinder imports none of PyCaret's
dashboard, model-explainer, cloud or NLP extras.

Whichever path you take, `xgboost`, `lightgbm` and `catboost` must be present:
they are the estimators inside the release model's voting ensemble, and
`load_model()` fails with `ModuleNotFoundError` without them. PyCaret's base
install does not pull in xgboost or catboost.

### Why Python 3.8 is pinned

`finalize_EPFinder_model.pkl` was serialized under Python 3.8.18 with
scikit-learn 1.2.2 and PyCaret 3.2.0, and PyCaret warns on any deviation at load
time. Python 3.8 reached end of life in October 2024, so this pin is a property
of the released model artifact, not of EPFinder itself — the source in this
repository runs on current Python. Moving to a supported Python requires
re-exporting the model and re-validating it against the GM12878 benchmark below.

If bedtools is already installed elsewhere on your system, point at it with
`bedtools_path` in `preprocessing/config.yaml`.

## Tests

```bash
pip install pytest
pytest tests/
```

The suite builds its own synthetic inputs, so it needs no reference genomics
data. Tests skip rather than fail when bedtools or the PyCaret stack is absent.
See `tests/README.md`.

## Quick prediction test

Run EPFinder prediction on the included benchmark-format matrix:

```bash
python scripts/EPFinder_predict.py \
  --input dataset/gm12878_29features_ML.tsv \
  --output examples/gm12878_EPFinder_predictions.tsv
```

If the input contains `#Class`, the script also reports AUROC and AUPRC.

To get a spreadsheet-friendly file instead, give the output a `.csv` extension:

```bash
python scripts/EPFinder_predict.py \
  --input dataset/gm12878_29features_ML.tsv \
  --output examples/gm12878_EPFinder_predictions.csv
```

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
  --output /path/to/output/EPFinder_predictions.tsv
```

Alternatively, run both steps with the wrapper:

```bash
python scripts/run_EPFinder_pipeline.py \
  --config preprocessing/config.local.yaml
```

## Output

`scripts/EPFinder_predict.py` appends an `EPFinder_score` column to the input
matrix. Higher scores indicate stronger model support for the SNP-promoter pair.

For GWAS applications, downstream ranking is typically performed per SNP or per
GWAS locus, depending on the biological question.

### Output format

The prediction CLI writes TSV by default and switches to CSV when the output
path ends in `.csv`. Use `--output-format` to set the delimiter explicitly when
the filename does not carry the extension:

```bash
python scripts/EPFinder_predict.py \
  --input  /path/to/output/EPFinder_29features_ML.tsv \
  --output /path/to/output/EPFinder_predictions.csv

# or, for an output path with a different extension
python scripts/EPFinder_predict.py \
  --input  /path/to/output/EPFinder_29features_ML.tsv \
  --output /path/to/output/EPFinder_predictions.txt \
  --output-format csv
```

The input delimiter is detected the same way, and `--input-format` overrides it.
Preprocessing itself always writes TSV, because the workflow uses commas
internally to build the enhancer/promoter merge keys.

### Opening results in Excel

CSV output opens directly in Excel, but import it with **Data > From Text/CSV**
and set the `Prom_gene` column type to **Text** rather than double-clicking the
file. Excel's default conversion silently rewrites gene symbols such as `SEPT2`,
`MARCH1` and `DEC1` into dates, and the original symbols cannot be recovered
once the file is saved.

## Notes

- Keep the feature order in `preprocessing/feature_list` unchanged unless you
  retrain the model.
- All coordinates and regulatory files must use hg38.
- The public config uses placeholder paths. Use `*.local.yaml` files for server
  paths; these are ignored by Git.
- Large local analysis outputs are intentionally ignored to keep the GitHub
  repository focused on the reproducible EPFinder workflow.
