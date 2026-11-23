# Binary Timing and Noise Toolkit

A source-based Python workflow for binary-pulsar timing and stochastic-noise inference with **PINT**, **enterprise**, and **emcee**. It provides independent **DDAB**, **DDHAB**, and **DDGRAB** timing-model extensions, conditional recovery of timing parameters, and checks of the local timing approximation against full nonlinear PINT residuals.

**The code was generated and organized with assistance from OpenAI Codex.**

This repository contains software and synthetic tests only. It includes **no observed TOAs, source ephemerides, fitted solutions, posterior chains, observatory clock tables, or private runtime configuration**. Supply your own phase-connected narrowband radio timing data. The in-memory models and TOAs in the tests are artificial software fixtures.

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Models and parameter conventions](#models-and-parameter-conventions)
- [Noise model and configuration](#noise-model-and-configuration)
- [Sampling, resuming, and output safety](#sampling-resuming-and-output-safety)
- [Reusing samples](#reusing-samples)
- [Results and statistical interpretation](#results-and-statistical-interpretation)
- [Testing](#testing)
- [Troubleshooting](#troubleshooting)
- [Repository scope and licensing](#repository-scope-and-licensing)

## What the workflow does

1. Load a user-supplied `.par` and `.tim`, register a requested custom binary model, and archive the input files, including recursive TIM `INCLUDE` dependencies.
2. Construct telescope-dependent white noise, an achromatic red Gaussian process (GP), and a chromatic dispersion-measure (DM) GP.
3. Sample noise hyperparameters with `emcee`. Optionally sample the Shapiro parameters nonlinearly alongside them. Other free timing parameters and GP coefficients are locally Gaussian-marginalized at the current reference model.
4. Recover paired timing/noise draws from the conditional Gaussian distributions.
5. Compare the linearized predictions with full PINT residuals. If needed, prepare a new reference solution and repeat the analysis.
6. Export the accepted timing solution, stochastic-model description, diagnostics, posterior summaries, and figures.

The standard DD/DDH/DDGR paths retain PINT's forward delay definitions. Instance-local derivative implementations handle the derivative issues addressed by this toolkit; the installed PINT source and its original global model registrations are not replaced. Custom models are registered under new names.

## Installation

Python **3.11** is the tested interpreter version. Run commands from the repository root; this release is a source checkout, not a wheel or a package with a console entry point.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python run_analysis.py --help
```

`pint-pulsar==1.1.7` is the pulsar timing library. The unrelated units package named `pint` is not a substitute. Dependencies include `enterprise-pulsar`, `emcee`, NumPy, SciPy, Astropy, Matplotlib, h5py, Loguru, corner, and platformdirs. The PINT version is pinned because the binary extensions use its internal model interfaces.

Some platforms need a compiler and SuiteSparse development libraries for enterprise's `scikit-sparse` dependency. If installation fails there, follow the upstream enterprise/scikit-sparse installation instructions or use a suitable conda environment; do not skip this dependency and assume the likelihood is equivalent. See [PINT](https://github.com/nanograv/PINT), [enterprise](https://github.com/nanograv/enterprise), and [emcee](https://emcee.readthedocs.io/).

The original dependency ranges are retained. They are not a fully locked environment. After installing and testing your environment, archive its package versions alongside your analysis. Local test-environment details are in [docs/VALIDATION.md](docs/VALIDATION.md).

### Clock, ephemeris, and observatory data

The workflow defaults to offline operation. PINT must already be able to resolve the ephemeris and clock corrections requested by your inputs. On an initial run, `--allow-downloads` permits Astropy/PINT downloads; IERS automatic downloads are still disabled by the loader. Use the same setting when resuming a chain.

Custom observatories and clocks can be supplied using PINT's `PINT_OBS_OVERRIDE` and `PINT_CLOCK_OVERRIDE` environment variables. None are shipped here. The toolkit does not patch or synthesize observatory clocks. A successful load is not proof that every historical epoch has measured clock coverage: inspect PINT warnings and the actual clock definitions. Matching a `CLK` label alone does not establish identical time conversions across timing programs.

## Quick start

All paths below refer to files **you provide**. There is deliberately no bundled observed `.par`/`.tim` pair.

```bash
python run_analysis.py /path/to/source.par /path/to/source.tim \
  --outdir results/standard \
  --config examples/default.json
```

For a fresh environment that needs runtime files, append `--allow-downloads`.

Prepare the data without launching MCMC:

```bash
python scripts/prepare.py /path/to/source.par /path/to/source.tim \
  --outdir results/prepared \
  --config examples/default.json
```

The input must be a phase-connected **binary**, with finite positive radio frequencies and more TOAs than free timing parameters. Timing fit flags in the par file determine which ordinary timing parameters are varied. During preparation, `START` and `FINISH` and old summary fields are removed from the working copy: **the workflow uses all TOAs actually loaded from the TIM**. Filter the TIM yourself if you intend to select a time interval.

Existing PINT noise components are removed from the in-memory model and replaced by the configured noise model. They are not automatically inherited. TIME offsets and JUMPs are handled by PINT; the toolkit does not pre-shift observations to emulate another program's JUMP convention.

### Custom models

```bash
python run_analysis.py /path/to/dd.par /path/to/source.tim \
  --outdir results/ddab --use-ddab

python run_analysis.py /path/to/ddh.par /path/to/source.tim \
  --outdir results/ddhab --use-ddhab

python run_analysis.py /path/to/ddgr.par /path/to/source.tim \
  --outdir results/ddgrab --use-ddgrab
```

Use exactly one custom-model flag. The flag is also required when the input already names the custom model. The DDHAB loader accepts DDH, DDHAB, and the TEMPO alias DDFWHE. This repository does not install a TEMPO or TEMPO2 binary; `.par`/`.tim` files are interpreted through PINT.

## Models and parameter conventions

Let `dt = t - DTHEPOCH`, measured in **seconds**, where `t` is the barycentric TDB time passed to the binary model.

| Model | Definition | Interpretation |
|---|---|---|
| DDAB | `DTH(t) = DTH + DTH1*dt + DTH2*dt**2/2` | Total deformation, with DD's independent post-Keplerian parameters |
| DDHAB | Same total-DTH expansion | DDH Shapiro parameterization using H3/STIGMA |
| DDGRAB | `DTH(t) = DTH_GR + DTHAB + DTHAB1*dt + DTHAB2*dt**2/2` | Additive deformation relative to the GR prediction at the current masses/orbit |

Each model uses `eTheta(t) = ecc(t) * (1 + DTH(t))`. The coefficients of order zero are dimensionless; first derivatives have units `s^-1`, and second derivatives `s^-2`. **There is no hidden `10^-12` scale** for these extension parameters. The factor `1/2` means the second coefficient is the actual second derivative.

`DTHEPOCH` is a **fixed TDB MJD**. The custom loaders preserve an explicit input epoch; if absent, they choose the midpoint of the actual barycentric TOA span. A par-only load cannot infer that span and requires an explicit epoch. Changing the epoch requires transforming all Taylor coefficients to represent the same curve.

Values, fit flags, and uncertainties supplied for the extension coefficients are preserved. Missing coefficients default to zero and free in the custom models; explicitly set the fit flags you want. For example, the following fragments fit a constant and a slope while fixing curvature:

```text
BINARY DDGRAB
DTHEPOCH 55000 0
DTHAB  0 1
DTHAB1 0 1
DTHAB2 0 0
```

For DDAB or DDHAB, use `DTH`, `DTH1`, and `DTH2` instead. These fragments are parameter syntax examples, not source ephemerides.

DDGRAB's GR contribution is recomputed as masses and orbital parameters change. DDAB/DDHAB parameterize total DTH and do not separate GR and aberration. The inherited DDGR `DTH` output is a GR-derived quantity, not the full time-dependent DDGRAB signal. The tested orbital parameterization uses `PB`; FB/ORBWAVES combinations are not established by this release's tests.

### Python API

From the repository root (or with this checkout on `PYTHONPATH`):

```python
from ddgrab import get_model_and_toas
from pint.fitter import DownhillWLSFitter

model, toas = get_model_and_toas(
    "/path/to/source.par", "/path/to/source.tim",
    planets=False, usepickle=False,
)
fitter = DownhillWLSFitter(toas, model)
fitter.fit_toas(maxiter=8)
```

Replace `ddgrab` with `ddab` or `ddhab` for those extensions. This API example is a conventional timing fit; it does not execute the joint noise workflow. For a custom par file with an explicit epoch, import the extension before calling PINT's native `get_model()`.

## Noise model and configuration

White noise is grouped **by telescope**, not by backend flags. The per-TOA standard deviation is

```text
sigma_white = EFAC * sqrt(sigma_TOA**2 + 10**(2*log10_T2EQUAD))
```

Both `sigma_TOA` and `10**log10_T2EQUAD` are in seconds. This is the enterprise `MeasurementNoise(..., log10_t2equad=...)` convention; it is not interchangeable with an additive term outside EFAC without converting the parameters.

Red noise is a global achromatic Fourier GP. The DM GP uses the enterprise dispersion-delay basis and power-law convention. No early-time nuisance GP or ECORR model is automatically added.

JSON configuration keys:

| Key | Default | Meaning |
|---|---:|---|
| `red_modes` | 20 | Number of red-noise Fourier modes; zero disables this GP |
| `dm_modes` | 20 | Number of DM Fourier modes; zero disables this GP |
| `dm_selection` | `"all"` | `"all"` or `"informative"` observing epochs |
| `dm_gamma` | `null` | Sample the DM spectral index; a number fixes it |
| `dm_span_factor` | 1.0 | Multiplier for the DM Fourier timespan |
| `dm_delay_threshold_us` | 20.0 | DM-information selection threshold in microseconds |
| `dd_shapiro_method` | `"pint"` | Full PINT evaluation or validated `"cached_delay"` |
| `priors` | See below | Partial overrides of uniform bounds |

Unknown keys are rejected. Configuration files override defaults; command-line `--red-components`/`--red-modes` and `--dm-components`/`--dm-modes` override the corresponding JSON values.

```bash
python run_analysis.py /path/to/source.par /path/to/source.tim \
  --outdir results/custom_modes --red-components 30 --dm-components 40
```

Default uniform prior bounds:

| Parameter | Bounds | Units/convention |
|---|---|---|
| `efac` | `[0.1, 5]` | Dimensionless |
| `log10_t2equad` | `[-8, -3]` | Log10 seconds |
| `red_log10_A` | `[-20, -10]` | Enterprise power-law amplitude |
| `red_gamma` | `[0, 7]` | Spectral index |
| `dm_log10_A` | `[-20, -8]` | Enterprise DM power-law amplitude |
| `dm_gamma` | `[0, 7]` | Spectral index when sampled |
| `M2` | `[0.01, 5]` | Solar masses; nonlinear DD/DDAB mode |
| `SINI` | `[0.01, 1.0]` | Nonlinear DD/DDAB mode |
| `H3` | `[0, 3e-6]` | Seconds; nonlinear DDH/DDHAB mode |
| `STIGMA` | `[0.1, 1.0]` | Nonlinear DDH/DDHAB mode |

These defaults are computational starting points, not scientifically justified priors for every source. `examples/white_only.json` disables both GPs; `examples/small_smoke.json` uses two modes each for inexpensive checks, not a validated science configuration.

With `dm_selection="informative"`, frequency leverage is evaluated after allowing an achromatic epoch offset. A single-frequency epoch cannot independently determine both an offset and DM. The epoch DM error bars and the GP reconstruction use the same observations and are not independent measurements.

### Nonlinear Shapiro sampling

```bash
python run_analysis.py /path/to/dd.par /path/to/source.tim \
  --outdir results/dd_joint --sample-dd-shapiro

python run_analysis.py /path/to/ddh.par /path/to/source.tim \
  --outdir results/ddhab_joint --use-ddhab --sample-ddh-shapiro
```

`--sample-dd-shapiro` samples M2/SINI for DD or DDAB. `--sample-ddh-shapiro` samples H3/STIGMA for DDH or DDHAB. These options are mutually exclusive and remove those parameters from the local linear timing block to avoid counting them twice. They sample the pair even if the input fit flags were fixed; initial values must lie strictly inside the specified priors. With neither option, ordinary free Shapiro parameters remain in the locally linearized timing block.

Uniform priors in M2/SINI and H3/STIGMA imply different probability measures. In DDH/DDHAB, M2/SINI are derived from H3/STIGMA and should not be added as extra independent fitting parameters. The cached-delay method retains nonlinear Shapiro dependence but uses local orbital geometry; full-residual checks still apply.

## Sampling, resuming, and output safety

Main defaults: 48 walkers, 4 workers, 5,000 minimum steps, 30,000 maximum steps, 4,000 recovered draws, and at most three relinearization passes. The sampler rejects a walker count below twice the sampled dimension. Random seeds are set in the scripts; the top-level workflow does not expose a seed option. Record the code revision and environment for reproducibility.

By default, convergence requires finite positive autocorrelation times, relative changes below 0.02, and more than 50 maximum autocorrelation times **after burn-in**, as well as the minimum-step and proposal-history conditions. Approximate walker-split R-hat is recorded as a diagnostic, not used as the stopping gate. Walkers in one ensemble are not independent chains. Reaching `--max-steps` alone is not convergence.

Run the same command and output directory again to resume a compatible saved `pass_XX/noise/chain.h5`. `--max-steps` is the **total** iteration cap, not additional steps. Increase it to extend an unfinished chain. Changes to inputs, TIM dependencies, configuration, model flags, runtime overrides, or walker count are rejected on resume. Worker count, recovery draws, convergence criteria, and total step cap can be changed where supported by the checks.

**Output directories can be cleared.** A fresh run with no resumable chain clears the selected output directory. `--restart` clears it even if chains exist. A `--reuse-noise` run also starts fresh and replaces output. Use a dedicated disposable results directory, not a folder containing other work. The entry point guards input, working, source-sample, and tool directories, but it cannot identify arbitrary unrelated files you put in an otherwise valid output directory.

```bash
python run_analysis.py /path/to/source.par /path/to/source.tim \
  --outdir results/standard --max-steps 60000
```

If all passes fail numerical verification, the command exits with an error and retains each pass for inspection. Numerical acceptance requires central closure WRMS at most 10 ns, checked-draw closure WRMS at most 100 ns, and no invalid physical draws among the quantities checked by the implementation. It does not establish adequacy of the noise model or historical calibration.

## Reusing samples

`--reuse-noise` expects a directory containing at least `posterior_samples.npy` and `pars.txt`, with provenance and convergence information for posterior reuse. Column names must match **all** sampled parameters; safe named permutations are supported. Nonlinear Shapiro modes also require the paired Shapiro columns.

| Mode | Meaning |
|---|---|
| `auto` | Identical-likelihood reuse if the prepared-data SHA256 matches; otherwise attempt reweighting |
| `identical` | Require a matching prepared-data SHA256 and documented convergence |
| `reweight` | Require the old likelihood, equal prior bounds/dimensions, and adequate posterior overlap |
| `conditional` | Recover timing conditional on the supplied noise/other sampled-parameter distribution; does not claim a newly sampled joint posterior |

```bash
python run_analysis.py /path/to/source.par /path/to/source.tim \
  --outdir results/reweighted \
  --reuse-noise /path/to/old/pass_01/noise \
  --reuse-mode reweight \
  --old-data /path/to/old/pass_01/prepared/likelihood_data.npz
```

Importance reweighting is rejected when ESS is below `max(50, 0.1 * number_of_source_samples)`. This ESS detects weight collapse, not the original chain's mixing. Conditional reuse is explicit and may use a nonconverged supplied distribution, which must not be presented as a converged joint posterior.

For a **new MCMC** initialized from old paired samples, use `--initial-samples /path/to/old/noise` instead. This changes initialization, not the target priors or likelihood.

## Results and statistical interpretation

Each `pass_XX` preserves preparation, sampling/reuse, and recovery products. After numerical acceptance, final products are also copied directly into `--outdir`.

| Product | Purpose |
|---|---|
| `result.par` | Timing posterior medians and half central 68% intervals, preserving fixed parameters |
| `stochastic_model.json` | Noise medians and configuration; retain with the par file |
| `analysis_summary.json` | Inference type, convergence, closure, physical checks, WRMS definitions, and correlated-noise goodness of fit |
| `analysis_report.md` | English readable summary |
| `final_result.json`, `passes.json` | Accepted-result and pass history indexes |
| `timing_parameter_ranges.csv`, `all_parameter_ranges.csv` | Timing/noise/derived intervals and parameter status |
| `residuals_per_toa.npz` | Chronologically ordered noise-subtracted residuals, frequencies, telescopes, and uncertainties |
| `residual_reference_solution.npz` | Coherent conditional timing and GP reference solution |
| `pass_XX/noise/chain.h5` | Resumable emcee chain |
| `pass_XX/noise/convergence_diagnostics.json` | Sampling or reuse diagnostics |
| `pass_XX/prepared/metadata.json` | Input provenance, versions, configuration, and removed noise components |

Figures are saved as PNG/PDF with plotting arrays. Products include timing/noise corner plots, residual diagnostics, DM evolution, and DTH/DTHAB evolution for applicable models. Taylor-curve uncertainty is propagated from paired samples. Replot existing output with:

```bash
python scripts/replot.py /path/to/results
```

Residuals subtract coherent conditional timing and GP means at median sampled parameters. That reference can differ from the component-wise median par file. The reported WRMS is weighted and mean-subtracted, with both original-TOA and adjusted-white-noise weights documented. Subtracting a GP conditional mean is not equivalent to a full covariance-whitening transformation.

The correlated-noise chi-squared uses `r_timing.T @ C^-1 @ r_timing`, with white noise and GP covariance in `C`. GP coefficients are marginalized rather than counted as separate deterministic fitted parameters. Degrees of freedom are nominal because noise hyperparameters were estimated from the same data; do not read the summary as an exact chi-squared-test probability.

Keep prepared data, paired samples, reference solutions, and diagnostics when archiving an analysis. `result.par` alone is not a complete stochastic posterior or a portable encoding of the runtime GP model. Wide/strongly nonlinear posteriors, physical boundaries, poor DM leverage, incorrect phase connection, and uncertain clock calibration require additional analysis beyond passing software checks.

## Testing

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 \
python -B -m unittest discover -s tests -v
```

Tests generate synthetic barycentric TOAs in temporary directories and use Astropy's `builtin` ephemeris. They need no observed data, downloaded DE ephemeris, or observatory clock files. Unset external `PINT_CLOCK_OVERRIDE` / `PINT_OBS_OVERRIDE` values for an isolated test run. Temporary test files and Python/plotting caches must be writable.

Coverage includes model registration isolation, zero-extension limits, full-delay finite differences, epoch changes, PINT injection/recovery, serialization, enterprise design matrices, likelihood comparisons with dense Gaussian algebra, sampler resume, reuse diagnostics, and input/output protection. See [docs/VALIDATION.md](docs/VALIDATION.md) for the performed checks and their limits.

## Troubleshooting

- **Unknown custom binary:** use the matching `--use-ddab`, `--use-ddhab`, or `--use-ddgrab` flag; in Python import the extension before loading a custom par file.
- **Missing DTHEPOCH:** use the model's loader with TOAs, or supply a fixed epoch. Do not fit it alongside the Taylor coefficients.
- **Clock/ephemeris download blocked:** prepopulate the appropriate cache or use `--allow-downloads`. Confirm observatory aliases, clock coverage and reference convention separately.
- **Unrecognized custom observatory:** supply your own PINT observatory definition; the original local observatory override is intentionally not distributed.
- **Rank-deficient timing matrix:** remove unconstrained or redundant free parameters and check phase connection and frequency/time coverage.
- **No convergence:** inspect traces, priors, autocorrelation estimates and boundaries; increase the total cap or improve parameterization. Short smoke runs are not scientific measurements.
- **Closure rejected:** improve the reference model, simplify poorly constrained parameters, use an appropriate nonlinear Shapiro option, or inspect posterior boundaries; do not merely suppress the checks.
- **Resume rejected:** restore the exact inputs/settings or choose a new directory. `--restart` destroys the old output in the chosen directory.
- **Outputs expose private information:** output products intentionally contain input snapshots, source identifiers, local provenance and fitted parameters. They belong outside this code-only repository and must not be committed accidentally.

## Repository scope and licensing

```text
run_analysis.py       orchestration CLI
requirements.txt      scientific Python dependencies
scripts/              preparation, likelihood, sampling, recovery, diagnostics
ddab/, ddhab/, ddgrab/  independent PINT extensions
examples/             generic JSON configurations only
tests/                synthetic unit/integration tests
docs/                 release scope and validation record
licenses/             retained upstream notice
```

The release omits all original observed data, par files, validation runs, chains, logs, plots, notebooks, archives, cached runtime data, and machine-specific launchers. `.gitignore` adds protection for common scientific data/output formats; it is not a substitute for checking staged files before publishing.

This project is licensed under the [MIT License](LICENSE). See [docs/RELEASE_NOTES.md](docs/RELEASE_NOTES.md) for changes made during packaging and [NOTICE.md](NOTICE.md) for upstream attribution and applicable third-party licenses.

Third-party code and adaptations remain subject to their respective licenses and attribution notices. See `LICENSE`, `NOTICE.md`, and the `licenses/` directory for details.

