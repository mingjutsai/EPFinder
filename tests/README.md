# Tests

Run everything:

```bash
pytest tests/
```

The suite needs no reference genomics data. `fixture.py` writes a two-SNP
synthetic dataset whose every feature value can be derived by hand, and the
prediction tests use the benchmark matrix and model already in the repository.

| File | Covers | Requires |
| --- | --- | --- |
| `test_preprocessing.py` | All nine steps end to end, across five input shapes | bedtools, pandas |
| `test_validation.py` | Preflight acceptance and rejection | bedtools |
| `test_prediction.py` | Model loading, benchmark performance, TSV/CSV parity | pycaret stack |

Tests skip rather than fail when bedtools or pycaret is missing, so the
preprocessing half can be developed without installing the model stack.

## Why the end-to-end test exists

The workflow's worst failures have been silent: it exited 0 and wrote a
full-looking matrix containing wrong numbers. Checking a return code or a row
count would not have caught either of the two bugs these tests now pin down, so
`test_preprocessing.py` compares the whole matrix against hand-derived values.

The fixture deliberately varies things that must not affect the result --
chromosome naming (`chr1` vs `1`) in each input file independently, and the
width of the SNP file -- because both were previously assumed rather than
detected.
