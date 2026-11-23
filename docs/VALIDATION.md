# Release validation

Date: 2026-11-23.

## Completed checks

- **96 synthetic unit/integration tests passed**, with no skips (35.712 seconds
  in the local environment). The final run unset external PINT clock and
  observatory overrides and disabled Python bytecode output.
- All 37 Python modules parsed successfully.
- An independent copy in a temporary directory successfully ran
  `run_analysis.py --help`, `scripts/prepare.py --help`, and
  `scripts/replot.py --help`. Imports of all three custom model packages resolved
  to that copy, not to its original parent analysis directory.
- The distribution was scanned for observed-source names, original source
  parameter signatures, personal machine paths, symlinks, binary artifacts,
  data files and ephemerides. No observed data or source ephemerides are included.
- Source SHA256 values and packaging edits are recorded in
  `source_manifest.json`. Runtime calculations are unchanged; the readable
  recovery report was translated to English and test fixtures were generalized.

The original injection amplitudes were too small to test all Taylor coefficients
at the original relative tolerance with the replacement synthetic orbit and the
sub-nanosecond TOA simulation floor. The artificial injected amplitudes were
increased; the original recovery assertions were retained. The orthometric
acceleration fixture now uses phase-connected artificial TOAs rather than a
regular time grid with large phase residuals. The original acceleration
acceptance threshold was retained. These are changes to synthetic test design,
not to the scientific likelihood, model equations, or acceptance thresholds.

## Local environment

Python 3.11, macOS. Distribution metadata reported:

| Distribution | Version |
|---|---|
| pint-pulsar | 1.1.7 |
| enterprise-pulsar | 0.0.0 (local distribution metadata; not a usable release pin) |
| emcee | 3.1.6 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| astropy | 8.0.1 |
| matplotlib | 3.11.1 |
| h5py | 3.16.0 |
| loguru | 0.7.3 |
| corner | 2.3.0 |
| platformdirs | 4.11.5 |
| scikit-sparse | 0.4.16 |

`requirements.txt` retains the source toolkit's dependency constraints,
including enterprise-pulsar >=3.4,<4. The local enterprise installation reports
non-release metadata, so this record is **not** proof of a clean installation
against every version allowed by those constraints. A fresh dependency install
and a remote GitHub Actions run were not performed. Run the synthetic suite in
the environment used for your own analysis and archive that environment.

## Reproduce the suite

From the repository root, with dependencies installed and external PINT
observatory/clock override variables unset:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 \
python -B -m unittest discover -s tests -v
```

The tests use artificial barycentric TOAs with Astropy's built-in ephemeris and
create temporary inputs/results. No scientific source data, ephemeris download,
or observatory clock table is needed. Logs and generated test results are not
included in the distribution.

These checks establish the exercised software behavior. They do not establish
MCMC convergence, model adequacy, prior robustness, clock accuracy, or a
scientific measurement for a real observed system.
