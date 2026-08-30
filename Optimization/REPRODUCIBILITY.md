# Computational claims and exact sources

The recorded environment was Linux x86-64 with Python 3.12.13, NumPy 2.3.5,
SciPy 1.17.0, Matplotlib 3.10.8, and python-flint 0.9.0. Dependency versions
are pinned in `certificate/requirements.txt` and
`experiments/requirements-quantitative.txt`.

## Program-to-claim map

| Program | Input/output | Scientific status |
|---|---|---|
| `certificate/verify_n3_certificate.py` | reads `n3_certificate.json`; writes `n3_certificate_summary.json` | **Proof-producing.** Outward-rounded Arb verification of four event boxes, the complete 23,973-cell path cover, strict signs, exact adjacency, and the stabilized tail. |
| `certificate/generate_n3_certificate.py` | proposes `n3_certificate.generated.json` | **Proposal only.** Floating-point generator; never proof. |
| `certificate/independent/rational_interval_check.py` | prints the report recorded in `verification/fixed_lattice_output.txt` | **Proof-producing independent cross-check.** Exact rational intervals and proved exponential remainder bounds; predictor values are accepted only after exact tests. |
| `certificate/independent/compute_decimal.py` | prints 90-digit centers and derivatives | **Diagnostic only.** Decimal arithmetic is not outward-rounded. |
| `experiments/search_switches.py` | prints numerical `n=3` event estimates and states | **Numerical baseline only.** Also supplies the branch solver imported by the path-data program. |
| `experiments/linear_family_search.py` | writes `explicit_family_m*_L16.json` | **Finite-scale numerical evidence only.** Exact rational family construction followed by dense sampling and Brent refinement. The all-size theorem is analytic. |
| `experiments/check_linear_family_decimal.py` | writes `explicit_family_m*_L16_decimal90.json` | **High-precision numerical cross-check only.** Not an interval certificate. |
| `experiments/grid_resolution_study.py` | writes `grid_resolution_study.json` | **Exact Decimal postprocessing of numerical inputs.** Supports the grid-resolution table but does not certify the input event radii. |
| `experiments/make_n3_path_figure.py` | writes `n3_path_figure_data.json` and, on request, an untracked PDF | **Mixed provenance.** Event markers and state itinerary come from the passing Arb summary; smooth branch-gap samples are numerical. |

Every program mentioned in the article's computational discussion is present.
No other executable is needed for a reported computation.

## Integrity

The frozen theorem-level certificate has SHA-256 digest

```text
8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2
```

Run:

```bash
sha256sum -c SHA256SUMS
printf '%s  %s\n' \
  '8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2' \
  'certificate/n3_certificate.json' | sha256sum -c -
```

Any failure means that the checkout does not match the archived release.

## Rerun modes

```bash
./reproduce.sh --certificate
./reproduce.sh --quick
./reproduce.sh --figure
./reproduce.sh --full
```

- `--certificate` runs only the primary Arb verifier.
- `--quick` runs the primary verifier and grid-resolution postprocessing.
- `--figure` regenerates `experiments/n3_path_figure_data.json` and an untracked
  rendered PDF.
- `--full` runs the two rigorous verifiers, the floating-point `n=3` baseline,
  all three finite-family constructions, their 90-digit checks, the grid
  study, and path-plot generation.

A successful primary verifier reports `status=PASS`, `precision_bits=256`, and
`validated_path_cells=23973`. A successful independent run reports
`status=PASS` and `validated_rational_sign_cells=62501`.

The six finite-family JSON files and `grid_resolution_study.json` are numerical
evidence and postprocessing, not theorem-level interval certificates. Their
recorded event counts for `m=4,5,6` are respectively 4, 6, and 8.

