"""End-to-end workflow tests.

These exist because the bugs they cover were silent: the workflow exited 0 and
wrote a full-looking matrix containing wrong numbers, so only comparing the
whole output against hand-derived values catches them.
"""

import os

import pytest

import fixture
from conftest import needs_bedtools


def read_matrix(path):
    with open(path) as handle:
        rows = [line.rstrip('\n').split('\t') for line in handle if line.strip()]
    header, body = rows[0], rows[1:]
    return [dict(zip(header, row)) for row in body]


def assert_matches_expected(path):
    rows = read_matrix(path)
    assert len(rows) == len(fixture.EXPECTED), 'wrong number of E-P pairs'
    for actual, expected in zip(rows, fixture.EXPECTED):
        for column, want in expected.items():
            assert column in actual, 'missing column %s' % column
            got = actual[column]
            if isinstance(want, float):
                assert float(got) == pytest.approx(want), \
                    '%s: expected %s, got %s' % (column, want, got)
            elif isinstance(want, int):
                assert int(got) == want, \
                    '%s: expected %s, got %s' % (column, want, got)
            else:
                assert got == want, '%s: expected %s, got %s' % (column, want, got)


# Every combination must give the same matrix: chromosome naming and SNP-file
# width are input conventions, not part of the biology.
SHAPES = [
    ('production', dict(snp_style='bare', tss_style='bare', signal_style='chr')),
    ('all-bare', dict(snp_style='bare', tss_style='bare', signal_style='bare')),
    ('all-chr', dict(snp_style='chr', tss_style='chr', signal_style='chr')),
    ('chr-snp-bare-tss', dict(snp_style='chr', tss_style='bare', signal_style='chr')),
    ('four-column-snps', dict(snp_columns=4)),
]


@needs_bedtools
@pytest.mark.parametrize('name,kwargs', SHAPES, ids=[s[0] for s in SHAPES])
def test_workflow_output_is_independent_of_input_conventions(
        name, kwargs, tmp_path, run_workflow):
    config = fixture.build(tmp_path / name, **kwargs)
    run_workflow(config)
    assert_matches_expected(os.path.join(config['output_dir'], config['output_file']))


@needs_bedtools
def test_chr_prefixed_input_does_not_zero_the_signal_features(tmp_path, run_workflow):
    """Regression: step 7 used to emit 'chrchr1' for chr-prefixed input.

    bedtools reported no overlap instead of failing, so every chromatin
    feature silently became 0.0 while the run exited 0.
    """
    config = fixture.build(tmp_path / 'chr', snp_style='chr', tss_style='chr',
                           signal_style='chr')
    run_workflow(config)
    rows = read_matrix(os.path.join(config['output_dir'], config['output_file']))
    signals = [float(row[c]) for row in rows
               for c in ('MarkX_Enh', 'MarkX_Prom', 'MarkY_Enh', 'MarkY_Prom')]
    assert all(value > 0 for value in signals), \
        'chromatin features collapsed to zero: %s' % signals


@needs_bedtools
def test_extra_snp_columns_do_not_corrupt_the_contact_value(tmp_path, run_workflow):
    """Regression: steps 2 and 4 used fixed column indices.

    With a 4-column SNP file, HiC_Contact was read from the Hi-C bin end, so a
    contact frequency became a genomic coordinate and rows were dropped.
    """
    config = fixture.build(tmp_path / 'wide', snp_columns=4)
    run_workflow(config)
    rows = read_matrix(os.path.join(config['output_dir'], config['output_file']))
    contacts = [float(row['HiC_Contact']) for row in rows]
    assert contacts == [5.5, 7.25], 'contact values corrupted: %s' % contacts
