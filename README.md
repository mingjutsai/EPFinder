# EPFinder

## Usage

### Step 1: Data Preprocessing
Before running the preprocessing script, configure the settings in `preprocessing/config.yaml`. This YAML file contains paths to input files, Hi-C data, genome annotations, and other parameters required for the analysis.

Once configured, run the preprocessing script:

```bash
python preprocessing/EPFinder_preprocessing.py
```

This script will generate the necessary input file for the prediction step.

### Step 2: Prediction
Open the `EPFinder.ipynb` notebook and run it, using the input file generated from Step 1 to perform the prediction.
