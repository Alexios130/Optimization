# Quantitative research supplement: reproducibility manifest

Date checked: 2026-08-30

This manifest covers the computational evidence accompanying the quantitative
switch-complexity work. The analytic all-(m) lower-bound proof appears in the
separately submitted manuscript and is not part of this code-only archive.
The finite (L=16) examples below are **numerical evidence, not
interval-certified results**.
They must not be described as certified unless a separate outward-rounded
verifier is supplied and passes.

## Minimal supplement contents

Include exactly the following research files; the random-search logs,
differential-evolution runs, and exploratory `alt_*`, `family_*`, `l16_*`,
and `target_*` files are not needed.

### Documentation and environment

- `README.md`
- `REPRODUCIBILITY_MANIFEST.md`
- `requirements-quantitative.txt`

### Executable scripts

- `search_switches.py` — two-branch solver and manuscript (n=3) baseline
  reproduction.
- `linear_family_search.py` — exact rational construction plus dense/Brent
  continuation for finite paired-family members.
- `check_linear_family_decimal.py` — independent 90-digit Decimal event and
  sign checks for the finite family members.
- `grid_resolution_study.py` — exact Decimal postprocessing of the recorded
  event radii, open tie widths, and uniform-grid coverage.
- `make_n3_path_figure.py` — numerical branch-gap rendering with certified
  event markers read from the passing (n=3) Arb summary.

### Inputs and generated evidence

- `explicit_family_m4_L16.json`
- `explicit_family_m4_L16_decimal90.json`
- `explicit_family_m5_L16.json`
- `explicit_family_m5_L16_decimal90.json`
- `explicit_family_m6_L16.json`
- `explicit_family_m6_L16_decimal90.json`
- `grid_resolution_study.json`
- `n3_path_figure_data.json`
- `../article/figures/n3_certified_path.pdf`

## Environment used for this rerun

```text
Python 3.12.13
NumPy 2.3.5
SciPy 1.17.0
Matplotlib 3.10.8
```

Install the three external packages with:

```bash
python -m pip install -r requirements-quantitative.txt
```

The remaining imports are from the Python standard library.

## Exact rerun commands and observed results

Run every command from the repository's `experiments/` directory.

### 1. Syntax/import check

```bash
python -m py_compile \
  search_switches.py \
  linear_family_search.py check_linear_family_decimal.py \
  grid_resolution_study.py make_n3_path_figure.py
```

Observed result: exit status 0, no output.

### 2. Numerical reproduction of the certified manuscript baseline

```bash
python search_switches.py --check-n3 --grid 301
```

Observed numerical result:

```text
K = 4
roots = 0.5245653271531585,
        1.0369491958727957,
        2.6961208910188073,
        2.8324298315608853
states = 2, 12, 1, 12, 2
```

This command is only a floating-point reproduction of the manuscript result; the
separate Arb certificate remains the rigorous source for the (n=3) claim.

### 3. Finite (L=16) members of the analytic family

```bash
python linear_family_search.py \
  --m 4 --L 16 --epsilon 1/30 \
  --roots 1/3,5/12 --grid 200001 \
  --output explicit_family_m4_L16.json

python linear_family_search.py \
  --m 5 --L 16 --epsilon 1/30 \
  --roots 1/3,5/12,9/20 --grid 200001 \
  --output explicit_family_m5_L16.json

python linear_family_search.py \
  --m 6 --L 16 --epsilon 1/30 \
  --roots 1/3,3/8,5/12,9/20 --grid 200001 \
  --output explicit_family_m6_L16.json
```

Observed results:

| (m) | (n=2m) | detected (K) | interval-state sequence | minimum scaled separation |
|---:|---:|---:|---|---:|
| 4 | 8 | 4 | `2,12,1,12,2` | (2.118266721740092\times10^{-5}) |
| 5 | 10 | 6 | `2,12,1,12,2,12,1` | (8.71451897754838\times10^{-7}) |
| 6 | 12 | 8 | `2,12,1,12,2,12,1,12,2` | (4.59159626986505\times10^{-8}) |

The rerun regenerated all three JSON files byte-for-byte.

### 4. Independent high-precision checks of those JSON files

```bash
python check_linear_family_decimal.py \
  explicit_family_m4_L16.json explicit_family_m4_L16_decimal90.json

python check_linear_family_decimal.py \
  explicit_family_m5_L16.json explicit_family_m5_L16_decimal90.json

python check_linear_family_decimal.py \
  explicit_family_m6_L16.json explicit_family_m6_L16_decimal90.json
```

Observed results:

| (m) | events | minimum Decimal event separation | minimum absolute event derivative |
|---:|---:|---:|---:|
| 4 | 4 | (2.1182667256420550\times10^{-5}) | (5.4839400558469065\times10^{-8}) |
| 5 | 6 | (8.714523565050413\times10^{-7}) | (2.0528435126150636\times10^{-9}) |
| 6 | 8 | (4.592404336953724\times10^{-8}) | (1.1523342476508126\times10^{-10}) |

All reconstructed event systems had residuals below
(2.3\times10^{-90}), and the Decimal interval probes reproduced the state
sequences in the preceding table.  The rerun regenerated all three Decimal
JSON files byte-for-byte.

### 5. Uniform-grid resolution of the finite-family tie intervals

```bash
python grid_resolution_study.py
```

The reference grid has 200,000 equal subintervals, hence 200,001 nodes and
spacing (2.5\times10^{-6}) on the scaled-radius interval ([0,0.5]).  A tie is
counted as hit only if a node lies strictly between its two recorded boundary
radii.  The observed postprocessing result is:

| (m) | (n) | (K) | open tie intervals | ties hit by 200,001 nodes | minimum tie width | grid points sufficient from (h<\text{minimum width}) |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 8 | 4 | 2 | 2 | (2.1182667256420550\times10^{-5}) | 23,606 |
| 5 | 10 | 6 | 3 | 1 | (8.714523565050413\times10^{-7}) | 573,756 |
| 6 | 12 | 8 | 4 | 0 | (4.592404336953724\times10^{-8}) | 10,887,545 |

The last column is an alignment-independent sufficient count, not the least
node count that happens to hit every interval for the displayed alignment.
This calculation concerns direct sampling of the active set.  It does not say
that `linear_family_search.py` misses the event boundaries: that program
tracks the two branch gaps and refines their sign changes.  All radii remain
90-digit numerical evidence rather than outward-rounded certificates.

### 6. Reproducible (n=3) path figure

Run this command from `experiments/`:

```bash
python make_n3_path_figure.py
```

Observed result:

```text
certificate_sha256=8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2
event_midpoints=['0.5245653271531519', '1.036949195872802',
                 '2.696120891018807', '2.832429831560884']
state_sequence=['2', '12', '1', '12', '2']
```

The command reproduces `n3_path_figure_data.json` and the vector PDF
`../article/figures/n3_certified_path.pdf` byte-for-byte in the stated
environment.  The program verifies the certificate SHA-256 before drawing.
The vertical markers and active-set ribbon come from the Arb certificate;
the two smooth branch-gap curves are explicitly labelled numerical samples.

## SHA-256 checksums

These hashes identify the rerun reviewed on 2026-08-30.

| file | SHA-256 |
|---|---|
| `README.md` | `e38728ccd7d55ae973d0db92c5ed8fa0c0b1c8ec5a7f6cfd7dd441f86c1f0f42` |
| `requirements-quantitative.txt` | `815360a00e0e1d7a128cade3de11b66d3fe997cbd6e2c69f8081eb9bc8b88197` |
| `search_switches.py` | `3c9137c20be945f20afbbd8feeb01173dbe14c5e2fe6284d8dea8b6a6aec4789` |
| `linear_family_search.py` | `cdc1ba2b9970a017a95f77b93d1aaee61adb99e0a559da7cb88739f3490716e1` |
| `check_linear_family_decimal.py` | `632ad011e2537058af99550b91f1273bcca6c1ab83f4e2afc8613595039dd014` |
| `grid_resolution_study.py` | `6087d791bda972a6e9031085acfcf43be18797098c93bf7d62fad82a19d03eae` |
| `make_n3_path_figure.py` | `fb7a55a0263a2f3e8ce498552534ba1c9387cae641aab9df83f65e8d022c2f1d` |
| `explicit_family_m4_L16.json` | `ee8305ad562e720d7bcbac51264208656982fabdd7fedda28b79c18effc7fa60` |
| `explicit_family_m4_L16_decimal90.json` | `49ffe81a3e93bf5c2c93e6ac2a1d678fbd7f1397977e17b9d4056e8c3bef6e46` |
| `explicit_family_m5_L16.json` | `c26b9a2c72abbc901a4e0b385f91d726f220da1a710acfd332c58a8614cc561f` |
| `explicit_family_m5_L16_decimal90.json` | `25bbccf5751494b642f2b74292a4a04d22dbdade1ba807d9ad68621cc6c6bf04` |
| `explicit_family_m6_L16.json` | `6e38808dd06c72477da217297e517fa5ba091694d5ef56483812abe8c77c23af` |
| `explicit_family_m6_L16_decimal90.json` | `3a6b22898c2b5ac5e9bac6036f6c014e1dda92170cbbeb344182105655936fbc` |
| `grid_resolution_study.json` | `9c100beb7b2965861d752a2acbb2f5918cd0eaf03ccb19e71e839b5acfefb561` |
| `n3_path_figure_data.json` | `f4443be7f958f70c251c614fe2084496b9e51a12b76ed63f958ab9db71f5ae11` |
| `../article/figures/n3_certified_path.pdf` | `e2a18dd0403202923adc37ab3063a166dac60564a069f7292576b0faa70011e4` |

## Claim boundary

- The analytic asymptotic family proof supplies the theorem-level lower
  bound, subject to mathematical review of that proof.
- The finite (L=16) paths illustrate that mechanism and provide
  reproducible numerical evidence.
- Dense grids, Brent refinement, high precision, small residuals, and
  nonzero numerical derivatives do not by themselves exclude every
  unobserved event and therefore are never labelled certification here.
- Direct grid sampling can miss a nonempty tie interval even when branch-gap
  continuation still detects and refines both of its boundary events.
