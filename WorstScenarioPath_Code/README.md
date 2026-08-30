# Reproducibility code for worst-scenario switch complexity

This code-only archive accompanies:

> Alexis Seferlis, *Switch Complexity of Worst-Scenario Paths in Robust
> Logistic Regression* (2026).

Contact: [sefalexis@gmail.com](mailto:sefalexis@gmail.com)

Repository: <https://github.com/Alexios130/mannalen>

The upload archive is named `WorstScenarioPath_Code.zip`. It contains the
exact rational inputs, the proof-producing certificate, recorded verification
outputs, and reproducible exploratory computations. The manuscript PDF and
LaTeX source are submitted separately and are deliberately not part of this
archive. Literature audits, internal review notes, and legacy ZIP files are
also excluded.

Start with `REPRODUCIBILITY.md` for the claim-to-program map and
`REVIEWER_GUIDE.md` for the shortest verification route.

For a one-command primary-certificate/table/figure check, run
`./reproduce.sh --quick`.  The `--full` mode additionally runs the independent
rational-interval verifier and rebuilds the finite-family numerical files; it
is intentionally much slower.  The script verifies `SHA256SUMS` before it
runs any program.

## The nine disclosed programs

The archive contains exactly the following nine Python programs. No
undisclosed program, external dataset, SQL database, web service, or
proprietary solver is needed.

| Program | Purpose | Evidentiary status |
|---|---|---|
| `certificate/verify_n3_certificate.py` | Verifies four event boxes, the complete 23,973-cell cover, derivative signs, bridge, and tail | Primary rigorous proof; 256-bit outward-rounded Arb |
| `certificate/independent/rational_interval_check.py` | Independently checks event windows, 62,501 sign cells, and the tail | Rigorous cross-check with exact rational intervals and proved Taylor remainders |
| `certificate/generate_n3_certificate.py` | Proposes a new rational certificate ledger | Floating-point proposal only; not proof |
| `certificate/independent/compute_decimal.py` | Prints 90-digit exploratory event centers | Diagnostic only; not proof |
| `experiments/search_switches.py` | Reproduces and searches switch paths numerically | Floating-point evidence only |
| `experiments/linear_family_search.py` | Constructs and scans finite members of the lower-bound family | Numerical illustration only |
| `experiments/check_linear_family_decimal.py` | Independently checks those finite members at high precision | Numerical illustration only |
| `experiments/grid_resolution_study.py` | Quantifies which numerical tie intervals a direct uniform grid samples | Decimal postprocessing of numerical evidence; not proof |
| `experiments/make_n3_path_figure.py` | Reproduces the article's vector path figure after checking the certificate hash | Certified markers and path ribbon; numerical smooth curves |

The analytic lower-bound result in the manuscript does not depend on the
finite-family JSON files.

## Primary certificate

The frozen proof object is `certificate/n3_certificate.json`, with SHA-256
digest

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

A successful run reports `status=PASS`, `precision_bits=256`, and
`validated_path_cells=23973`.

## Independent cross-check

```bash
python3 certificate/independent/rational_interval_check.py
```

A successful run reports `status: PASS`, lattice denominator `10^18`,
`exp_terms: 32`, and `validated_rational_sign_cells: 62501`. NumPy and SciPy
only propose affine predictors; every accepted enclosure and every pass/fail
decision uses exact `Fraction` interval endpoints and proved exponential
remainders. Decimal numbers printed during the run are labelled as display
summaries and are not used as proof data.

## Package integrity and layout

```bash
sha256sum -c SHA256SUMS
```

`SHA256SUMS` does not hash itself. The authoritative proof objects are the
executable verifiers and the frozen certificate, not the recorded text logs.

- `certificate/` contains the frozen certificate, two verifiers, and two
  non-proof proposal/diagnostic programs.
- `experiments/` contains five explicitly non-certified numerical programs
  and their frozen JSON inputs and outputs.
- `verification/` contains recorded successful runs for comparison.
- The root documentation explains execution, evidentiary boundaries, and
  GitHub upload.

## Public release

No software licence has been selected. Ordinary copyright therefore applies
unless the author later adds a licence. For peer review, a public repository
and a tag identifying the reviewed version are recommended. Before final
publication, the reviewed tag can also be archived with a permanent identifier
such as a Zenodo DOI.
