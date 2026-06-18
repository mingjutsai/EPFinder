#!/usr/bin/env python3
"""Run EPFinder model prediction on a 29-feature matrix."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd
from pycaret.classification import load_model, predict_model
from sklearn.metrics import average_precision_score, roc_auc_score


DEFAULT_METADATA_COLUMNS = ("#Class", "ID", "Enh", "Prom", "TX")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply a trained EPFinder PyCaret model to a preprocessing output matrix."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input TSV containing EPFinder metadata columns plus the 29 model features.",
    )
    parser.add_argument(
        "--model",
        default="finalize_EPFinder_model",
        help=(
            "PyCaret model path. Provide the path without .pkl, or provide a .pkl path "
            "and the suffix will be removed for pycaret.load_model()."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV with EPFinder prediction scores appended.",
    )
    parser.add_argument(
        "--metadata-columns",
        default=",".join(DEFAULT_METADATA_COLUMNS),
        help=(
            "Comma-separated columns to exclude from model features and retain as metadata. "
            "Default: #Class,ID,Enh,Prom,TX."
        ),
    )
    parser.add_argument(
        "--score-column",
        default="EPFinder_score",
        help="Name of the output EPFinder score column.",
    )
    parser.add_argument(
        "--label-column",
        default="#Class",
        help="Optional truth label column used to report AUROC/AUPRC if present.",
    )
    return parser.parse_args()


def normalize_model_path(model_path: str) -> str:
    path = Path(model_path)
    if path.suffix == ".pkl":
        return str(path.with_suffix(""))
    return model_path


def select_score_column(columns: Iterable[str]) -> str:
    columns = list(columns)
    preferred = (
        "prediction_score_1",
        "Score_1",
        "prediction_score",
        "Label_score",
        "score",
    )
    for name in preferred:
        if name in columns:
            return name
    score_like = [c for c in columns if "score" in c.lower()]
    if len(score_like) == 1:
        return score_like[0]
    raise ValueError(
        "Could not identify prediction score column in PyCaret output. "
        f"Available columns: {columns}"
    )


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    metadata_columns = [c for c in args.metadata_columns.split(",") if c]

    data = pd.read_csv(input_path, sep="\t")
    missing = [c for c in metadata_columns if c not in data.columns]
    if missing:
        print(f"Metadata columns not present and ignored: {', '.join(missing)}")

    feature_columns = [c for c in data.columns if c not in metadata_columns]
    if not feature_columns:
        raise ValueError("No feature columns remain after metadata-column removal.")

    model = load_model(normalize_model_path(args.model))
    predictions = predict_model(model, raw_score=True, data=data[feature_columns])
    score_source = select_score_column(predictions.columns)

    output = data.copy()
    output[args.score_column] = predictions[score_source].values
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False)

    print(f"Wrote predictions: {output_path}")
    print(f"Rows: {len(output):,}")
    print(f"Score column: {args.score_column} <- {score_source}")

    if args.label_column in output.columns:
        y_true = output[args.label_column]
        y_score = output[args.score_column]
        print(f"AUROC: {roc_auc_score(y_true, y_score):.4f}")
        print(f"AUPRC: {average_precision_score(y_true, y_score):.4f}")


if __name__ == "__main__":
    main()
