"""Prediction CLI tests.

Beyond scoring, these cover the install path itself: the ensemble inside
finalize_EPFinder_model.pkl embeds xgboost, lightgbm and catboost estimators,
so an environment missing any of them cannot even load the model.
"""

import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCHMARK = os.path.join(REPO, 'dataset', 'gm12878_29features_ML.tsv')
PREDICT = os.path.join(REPO, 'scripts', 'EPFinder_predict.py')

# Published performance of the release model on the bundled benchmark.
EXPECTED_AUROC = 0.9389
EXPECTED_AUPRC = 0.9237

pytest.importorskip('pycaret', reason='pycaret not installed')


def run_predict(output, *extra):
    result = subprocess.run(
        [sys.executable, PREDICT, '--input', BENCHMARK, '--output', str(output)] + list(extra),
        cwd=REPO, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def read(path, sep):
    import pandas as pd
    return pd.read_csv(path, sep=sep)


@pytest.mark.parametrize('module', ['xgboost', 'lightgbm', 'catboost'])
def test_estimator_libraries_are_installed(module):
    """Missing any of these makes load_model() raise ModuleNotFoundError."""
    pytest.importorskip(module, reason='%s missing: the model cannot load' % module)


def test_benchmark_reproduces_published_performance(tmp_path):
    out = tmp_path / 'pred.tsv'
    stdout = run_predict(out)
    assert 'AUROC: %.4f' % EXPECTED_AUROC in stdout, stdout
    assert 'AUPRC: %.4f' % EXPECTED_AUPRC in stdout, stdout

    frame = read(out, '\t')
    assert len(frame) == 89
    assert 'EPFinder_score' in frame.columns
    assert frame['EPFinder_score'].between(0, 1).all()


def test_csv_and_tsv_outputs_carry_identical_scores(tmp_path):
    tsv, csv = tmp_path / 'pred.tsv', tmp_path / 'pred.csv'
    run_predict(tsv)
    run_predict(csv)

    a, b = read(tsv, '\t'), read(csv, ',')
    assert list(a.columns) == list(b.columns)
    assert len(a) == len(b)
    assert (a['EPFinder_score'] - b['EPFinder_score']).abs().max() == 0


def test_output_format_flag_overrides_the_extension(tmp_path):
    out = tmp_path / 'pred.txt'
    stdout = run_predict(out, '--output-format', 'csv')
    assert 'Format: csv' in stdout
    assert ',' in open(str(out)).readline()
