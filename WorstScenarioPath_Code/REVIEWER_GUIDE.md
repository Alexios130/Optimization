# Reviewer guide

This archive contains code and machine-readable evidence only. Please read the
manuscript supplied by the journal separately; no manuscript PDF or LaTeX
source is included here.

## Suggested review order

1. Check the frozen certificate digest shown below.
2. Run `certificate/verify_n3_certificate.py`, the primary Arb verifier.
3. If an implementation-independent cross-check is desired, run
   `certificate/independent/rational_interval_check.py`.
4. Consult `experiments/REPRODUCIBILITY_MANIFEST.md` only for the explicitly
   numerical examples. Those examples are not interval certificates.

## Computational claim-to-evidence map

| Manuscript item | Program and input | Evidentiary status |
|---|---|---|
| Four certified event-radius and derivative enclosures | `certificate/verify_n3_certificate.py` with `certificate/n3_certificate.json` | Primary rigorous certificate |
| Complete 23,973-cell continuation cover and finite sign chart | Same primary verifier and certificate | Primary rigorous certificate |
| Certified intercept-only root, stabilization radius, tail gap, and bridge bound | Same primary verifier and certificate | Primary rigorous certificate |
| Independent event, sign-cell, and tail checks | `certificate/independent/rational_interval_check.py` | Rigorous independent cross-check |
| Floating-point reproduction of the certified path | `experiments/search_switches.py` | Numerical reproduction only |
| Finite members of the analytic lower-bound family | `experiments/linear_family_search.py` and `experiments/check_linear_family_decimal.py` | Numerical illustration only |
| Direct-grid resolution table | `experiments/grid_resolution_study.py` | Decimal postprocessing of numerical event radii; not a certificate |
| Certified-path figure | `experiments/make_n3_path_figure.py` | Certified markers and ribbon; numerical smooth curves |

The generator and 90-digit center calculator are proposal and diagnostic
programs. They are not accepted as proof. The full inventory of exactly nine
programs appears in `REPRODUCIBILITY.md`.

## Exact certified instance

```text
n=3, p=1, S=2, lambda=1/2
x=(7/2, 3, 4)
y=(-1, +1, -1)
pi_1=(12/19, 1/19, 6/19)
pi_2=(1/8, 1/16, 13/16)
```

The four radius enclosures reported by the primary verifier are

```text
rho_1: [0.5245653271531518, 0.5245653271531519]
rho_2: [1.0369491958728015, 1.0369491958728016]
rho_3: [2.6961208910188072, 2.6961208910188073]
rho_4: [2.8324298315608840, 2.8324298315608841]
```

## Primary verification

The frozen certificate has SHA-256 digest

```text
8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2
```

From the extracted archive root, run:

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r certificate/requirements.txt
printf '%s  %s\n' \
  '8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2' \
  'certificate/n3_certificate.json' | sha256sum -c -
python3 certificate/verify_n3_certificate.py
```

The required successful fields are

```json
"status": "PASS",
"precision_bits": 256,
"validated_path_cells": 23973
```

The verifier uses `python-flint`/Arb outward-rounded balls. It does not import
NumPy or SciPy, refuses execution under `python -O`, and terminates on a failed
inclusion, non-strict required sign, uncovered boundary, or wrong event order.

## Independent verification

```bash
python3 certificate/independent/rational_interval_check.py
```

This implementation uses exact `Fraction` endpoints on the lattice with
denominator `10^18`, together with 32 Taylor terms and proved exponential
remainders. NumPy and SciPy only propose affine predictors. A proposed cell is
accepted solely by the exact rational interval tests. Decimal values printed
during execution are labelled as display summaries and are not proof inputs.

The required final fields are

```text
status: PASS
lattice_denominator: 1000000000000000000
exp_terms: 32
validated_rational_sign_cells: 62501
```

The recorded run took approximately 612 seconds as a single CPU-bound process;
runtime is hardware-dependent.

## Evidence boundaries

- `certificate/generate_n3_certificate.py` proposes boxes using floating-point
  solvers. Its output is not accepted as proof.
- `certificate/independent/compute_decimal.py` is diagnostic. Ninety decimal
  digits do not replace outward-rounded interval arithmetic.
- The finite lower-family JSON files illustrate the analytic construction;
  they do not prove its all-dimension statement.
- Dense grids, Brent refinement, high precision, small residuals, and nonzero
  numerical derivatives do not exclude every unobserved event.

For package integrity, run `sha256sum -c SHA256SUMS`. A public tag identifying
the exact version reviewed is recommended; a permanent archive and DOI can be
added before final publication.
