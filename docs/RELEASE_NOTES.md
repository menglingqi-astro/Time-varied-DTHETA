# Code-only distribution notes

Prepared 2026-11-23 from the local binary timing/noise toolkit.

## Included

- Main orchestration entry point and the complete Python `scripts/` workflow.
- DDAB, DDHAB, and DDGRAB packages, including the reusable fixed-noise utility.
- Seven synthetic `test_*.py` regression modules.
- Three generic noise configurations; no source-specific prior/configuration examples.
- A new English README, English model summaries, output report and publication exclusions.

## Packaging changes

- Replaced all source-derived synthetic fixture identifiers, coordinates,
  orbital/mass/spin values and high-precision reference values with arbitrary
  artificial choices. These fixtures are not a timing solution for an observed
  source. Existing test assertions are retained except where synthetic numerical
  values must change consistently.
- Strengthened artificial injection signals above the simulation precision floor
  and phase-connected the orthometric acceleration fixture; kept the original
  numerical acceptance assertions. See `VALIDATION.md`.
- Replaced the tests' DE ephemeris dependency with Astropy `builtin`, using
  barycentric synthetic TOAs. Removed hard-coded interpreter paths.
- Translated the generated readable report to English (`analysis_report.md`).
  JSON/array schemas and science calculations are preserved.
- Excluded legacy comparison utilities tied to external parent-directory code
  and actual validation datasets. General execution paths remain included.
- Kept the original scientific algorithms, model definitions, configuration
  defaults and dependency ranges; this is not a reanalysis or model revision.

The original analysis directory is not modified by these packaging edits.
No measurements, noise samples, fitted values, timing residuals, private clocks,
observatory overrides or input par/tim files are included. Generated test output
is outside the distribution tree. `source_manifest.json` records relative source
paths and hashes of copied Python files, not observations or machine paths.

The maintainer selected the MIT License for original toolkit contributions.
Upstream notices remain applicable; see `LICENSE` and `NOTICE.md`.
