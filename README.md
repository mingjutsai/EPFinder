# EPFinder

## Environment Setup

Create and activate the conda environment:

```bash
conda env create -f conda/EPFinder_env.yml
conda activate EPFinder_env
```

Ensure that all the following steps are run within this activated environment.

## Usage

### Step 1: Data Preprocessing
Before running the preprocessing script, configure the settings in `preprocessing/config.yaml`. This YAML file contains paths to input files, Hi-C data, genome annotations, and other parameters required for the analysis.

Once configured, run the preprocessing script:

```bash
python preprocessing/EPFinder_preprocessing.py preprocessing/config.yaml
```

This script will generate the necessary input file for the prediction step.

### Step 2: Prediction
The user needs to register the kernel first. To register the kernel, run:

```bash
python -m ipykernel install --user --name EPFinder_env --display-name "EPFinder_env"
```

Then, they can open the `EPFinder_prediction.ipynb` notebook, choose the EPFinder_env as the kernel, and run it using the input file generated from Step 1 to perform the prediction.
