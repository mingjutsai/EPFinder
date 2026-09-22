"""Preflight tests.

The validator's job is to turn input problems into loud failures before the
nine steps start, so each test asserts on the message as well as the failure.
"""

import shutil

import pytest

import fixture
from input_validation import validate_inputs
from conftest import needs_bedtools


@pytest.fixture
def config(tmp_path):
    cfg = fixture.build(tmp_path / 'inputs')
    cfg['bedtools_path'] = shutil.which('bedtools') or 'bedtools'
    return cfg


@needs_bedtools
def test_accepts_the_fixture_and_records_the_detected_layout(config):
    validate_inputs(config)
    assert config['gwas_columns'] == 3
    assert config['tss_chr_prefix'] is False
    assert config['feature_chr_prefix'] is True


@needs_bedtools
def test_detects_chromosome_conventions_per_file(tmp_path):
    cfg = fixture.build(tmp_path / 'chr', snp_style='chr', tss_style='chr',
                        signal_style='bare')
    cfg['bedtools_path'] = shutil.which('bedtools')
    validate_inputs(cfg)
    assert cfg['tss_chr_prefix'] is True
    assert cfg['feature_chr_prefix'] is False


@needs_bedtools
def test_rejects_a_ragged_snp_file(config):
    with open(config['input_gwas'], 'a') as handle:
        handle.write('3\t70000\trs333\textra\n')
    with pytest.raises(ValueError, match='same number of columns'):
        validate_inputs(config)


@needs_bedtools
def test_rejects_feature_files_that_disagree_on_chromosome_naming(config, tmp_path):
    odd = str(tmp_path / 'odd.bedGraph')
    with open(odd, 'w') as handle:
        handle.write('1\t9000\t11000\t2.0\n')
    with open(config['feature_list'], 'a') as handle:
        handle.write('MarkZ\t%s\n' % odd)
    with pytest.raises(ValueError, match="disagree on the 'chr' prefix"):
        validate_inputs(config)


@needs_bedtools
def test_rejects_expression_files_with_no_matching_identifiers(config):
    with open(config['tx_expression'], 'w') as handle:
        handle.write('ENST99999999999.1\t1.0\n')
    with pytest.raises(ValueError, match='no matching stable Ensembl IDs'):
        validate_inputs(config)


@needs_bedtools
def test_rejects_a_missing_hic_chromosome(config, tmp_path):
    import os
    os.remove(os.path.join(config['hic_folder'], 'TEST.hic.KR.chr2'))
    with pytest.raises(ValueError, match='Hi-C chromosome 2'):
        validate_inputs(config)
