#!/usr/bin/env python3
"""Run EPFinder preprocessing followed by model prediction."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import yaml


PATH_CONFIG_KEYS = (
    "input_gwas",
    "hic_folder",
    "tss_file",
    "tx_expression",
    "gene_list",
    "gene_expression",
    "feature_list",
    "output_dir",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the EPFinder preprocessing-to-prediction workflow."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="YAML config consumed by preprocessing/EPFinder_preprocessing.py.",
    )
    parser.add_argument(
        "--model",
        default="finalize_EPFinder_model",
        help="PyCaret EPFinder model path, with or without .pkl suffix.",
    )
    parser.add_argument(
        "--prediction-output",
        default=None,
        help="Output TSV for EPFinder predictions. Defaults to <output_dir>/<output_file>.predictions.tsv.",
    )
    parser.add_argument(
        "--skip-preprocessing",
        action="store_true",
        help="Skip preprocessing and run prediction on the configured output matrix.",
    )
    return parser.parse_args()


def load_config(path: Path) -> dict:
    with path.open("r") as handle:
        config = yaml.safe_load(handle)
    for key in PATH_CONFIG_KEYS:
        value = config.get(key)
        if value and not Path(str(value)).is_absolute():
            config[key] = str((path.parent / str(value)).resolve())
    return config


def run(cmd: list[str]) -> None:
    print("+ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    config_path = Path(args.config).resolve()
    config = load_config(config_path)

    output_dir = Path(config["output_dir"]).expanduser().resolve()
    feature_matrix = output_dir / config["output_file"]
    prediction_output = (
        Path(args.prediction_output).expanduser().resolve()
        if args.prediction_output
        else feature_matrix.with_suffix(feature_matrix.suffix + ".predictions.tsv")
    )

    if not args.skip_preprocessing:
        run(
            [
                sys.executable,
                str(repo_root / "preprocessing" / "EPFinder_preprocessing.py"),
                str(config_path),
            ]
        )

    run(
        [
            sys.executable,
            str(repo_root / "scripts" / "EPFinder_predict.py"),
            "--input",
            str(feature_matrix),
            "--model",
            args.model,
            "--output",
            str(prediction_output),
        ]
    )


if __name__ == "__main__":
    main()
