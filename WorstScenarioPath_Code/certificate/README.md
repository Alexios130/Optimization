# Rigorous certificate for the three-observation path

The model inputs are exact rationals. The frozen proof object is
`n3_certificate.json`; its SHA-256 digest is
`8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2`.
Verify this file before running any proposal code.

## Clean verification

From the repository root:

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r certificate/requirements.txt
printf '%s  %s\n' \
  '8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2' \
  'certificate/n3_certificate.json' | sha256sum -c -
python3 certificate/verify_n3_certificate.py
```

The verifier writes a JSON summary whose exact successful fields include:

```json
"status": "PASS",
"arithmetic": "python-flint Arb outward-rounded balls",
"precision_bits": 256,
"validated_path_cells": 23973
```

The verifier uses `python-flint==0.9.0` and Python's standard library. It
checks four interval-Newton boxes, the complete 23,973-cell continuation
cover, exact rational adjacency, derivative and gap signs, and the nonsmooth
tail. A failed inclusion, non-strict sign, uncovered boundary, or wrong event
order terminates the program.

## Independent rational check

Run:

```bash
python3 certificate/independent/rational_interval_check.py
```

This separately written program validates 62,501 rational sign cells using
`fractions.Fraction`, a fixed decimal lattice, and proved Taylor/geometric
remainder bounds. NumPy and SciPy are used only to propose affine predictors.
Event brackets are printed exactly as fractions. Other decimal interval and
width values are explicitly labelled as display summaries; they are not used
by any pass/fail decision.
The latest tested run took 611.93 seconds as a single CPU-bound process and
printed its summary at completion; runtime is hardware-dependent. The primary
Arb verifier took 1.86 seconds in the corresponding check.

## Optional proposal and diagnostic programs

`certificate/generate_n3_certificate.py` uses NumPy and SciPy floating-point solves to
propose a new ledger. It is not proof-producing. Its safe default output is
`certificate/n3_certificate.generated.json`, leaving the frozen proof object
unchanged:

```bash
python3 certificate/generate_n3_certificate.py
```

`certificate/independent/compute_decimal.py` prints 90-digit exploratory centres and
derivatives. Decimal precision is not outward-rounded interval arithmetic, so
its output is diagnostic only.

Run it from the repository root with:

```bash
python3 certificate/independent/compute_decimal.py
```

The environment used for the recorded run was Linux x86-64, Python 3.12.13,
NumPy 2.3.5, SciPy 1.17.0, and python-flint 0.9.0.
