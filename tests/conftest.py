import os
import shutil
import subprocess
import sys

import pytest
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREPROCESSING = os.path.join(REPO, 'preprocessing')
sys.path.insert(0, PREPROCESSING)

needs_bedtools = pytest.mark.skipif(
    shutil.which('bedtools') is None, reason='bedtools not installed')


@pytest.fixture
def run_workflow(tmp_path):
    """Write a config and run the whole preprocessing workflow on it."""

    def run(config, name='config.yaml', expect_failure=False):
        path = str(tmp_path / name)
        with open(path, 'w') as handle:
            yaml.safe_dump(config, handle)
        result = subprocess.run(
            [sys.executable, os.path.join(PREPROCESSING, 'EPFinder_preprocessing.py'), path],
            cwd=str(tmp_path), capture_output=True, text=True)
        if not expect_failure and result.returncode != 0:
            raise AssertionError(
                'workflow failed (%s)\nstdout:\n%s\nstderr:\n%s'
                % (result.returncode, result.stdout, result.stderr))
        return result

    return run
