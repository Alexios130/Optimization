# Reproducibility guide

**Article:** *Switch Complexity of Worst-Scenario Paths in Robust Logistic Regression*  
**Author:** Alexis Seferlis  
**Email:** sefalexis@gmail.com  
**Repository:** <https://github.com/Alexios130/mannalen>

This repository contains every program and machine-readable input used for a
computer-assisted proof, numerical table, or generated figure in the article.
No external dataset, SQL database, web service, proprietary solver, or
undisclosed program is required.

## Evidence boundary

The theorem-level computer-assisted object is
`certificate/n3_certificate.json`, with SHA-256 digest

```text
8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2
```

`certificate/verify_n3_certificate.py` is the primary verifier. It uses Arb
outward-rounded ball arithmetic through `python-flint` and checks the complete
finite sign cover, the four interval-Newton event boxes, and the bridge to the
exact tail. `certificate/independent/rational_interval_check.py` is a
separately implemented rigorous cross-check.

All remaining calculations are explicitly labelled proposal, diagnostic, or
numerical evidence. In particular, the finite `L=16` family members are not
interval-certified theorems.

## Claim-to-program map

| Article item | Program | Principal input/output | Status |
|---|---|---|---|
| Four certified event enclosures and exhaustive `n=3` path | `certificate/verify_n3_certificate.py` | `certificate/n3_certificate.json`, `certificate/n3_certificate_summary.json` | Primary rigorous certificate |
| Independent event, sign-cell, and tail check | `certificate/independent/rational_interval_check.py` | `verification/fixed_lattice_output.txt` | Rigorous independent cross-check |
| Proposal of certificate boxes and slabs | `certificate/generate_n3_certificate.py` | generated JSON ledger | Floating-point proposal only |
| 90-digit `n=3` event diagnostics | `certificate/independent/compute_decimal.py` | terminal report | Numerical diagnostic only |
| Basic floating-point reproduction of the `n=3` path | `experiments/search_switches.py` | terminal JSON | Numerical only |
| Finite members of the analytic lower family | `experiments/linear_family_search.py` | `experiments/explicit_family_m*_L16.json` | Numerical illustration only |
| 90-digit check of finite family members | `experiments/check_linear_family_decimal.py` | `experiments/explicit_family_m*_L16_decimal90.json` | Numerical illustration only |
| Direct-grid resolution table | `experiments/grid_resolution_study.py` | `experiments/grid_resolution_study.json` | Exact Decimal postprocessing of numerical event radii; not certification |
| Certified-path figure | `experiments/make_n3_path_figure.py` | `experiments/n3_path_figure_data.json`, `article/figures/n3_certified_path.pdf` | Certified markers/path ribbon; numerical smooth curves |

## Rigorous verification

Run from the repository root:

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r certificate/requirements.txt
printf '%s  %s\n' \
  '8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2' \
  'certificate/n3_certificate.json' | sha256sum -c -
python3 certificate/verify_n3_certificate.py
python3 certificate/independent/rational_interval_check.py
```

The primary verifier must report `status=PASS`, `precision_bits=256`, and
`validated_path_cells=23973`. The independent checker must report
`status: PASS`, `lattice_denominator: 1000000000000000000`, `exp_terms: 32`,
and `validated_rational_sign_cells: 62501`.

Alternatively, `./reproduce.sh --quick` verifies package checksums, runs the
primary Arb certificate, and regenerates the article's numerical table and
figure.  `./reproduce.sh --full` also runs the independent rigorous checker
and reconstructs all finite-family numerical outputs; this mode can take many
minutes, depending on the machine.

## Numerical table and figure

Install the pinned quantitative dependencies and run:

```bash
python3 -m pip install -r experiments/requirements-quantitative.txt
python3 experiments/grid_resolution_study.py
python3 experiments/make_n3_path_figure.py
```

The grid study reads the three frozen 90-digit finite-family files and writes
`experiments/grid_resolution_study.json`. The figure program checks the
certificate hash before drawing and writes both its source data and the vector
PDF used by the article.

The complete commands for regenerating the finite-family inputs, the expected
terminal output, software versions, and per-file checksums are in
`experiments/REPRODUCIBILITY_MANIFEST.md`. Recorded text output is supplied
for comparison; executing the verifier is authoritative.

## Repository layout

- `certificate/`: exact model data, frozen proof object, primary verifier,
  and independent checker;
- `experiments/`: numerical illustration scripts, inputs, outputs, and pinned
  quantitative dependencies;
- `verification/`: preserved terminal reports from reviewed rigorous runs;
- `article/figures/`: the reproducible vector figure included in the paper.

Every path in the article is relative to this repository root.
