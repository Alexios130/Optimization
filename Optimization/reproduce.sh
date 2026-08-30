#!/usr/bin/env bash
set -euo pipefail

PACKAGE_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$PACKAGE_ROOT"

if [[ -f SHA256SUMS ]]; then
  sha256sum -c SHA256SUMS
fi

printf '%s  %s\n' \
  '8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2' \
  'certificate/n3_certificate.json' | sha256sum -c -

if [[ ! -x .venv/bin/python ]]; then
  python3 -m venv .venv
fi

.venv/bin/python -m pip install -r certificate/requirements.txt
.venv/bin/python -m pip install -r experiments/requirements-quantitative.txt

run_figure() {
  .venv/bin/python experiments/make_n3_path_figure.py \
    --pdf experiments/generated/n3_certified_path.pdf
}

run_quick() {
  .venv/bin/python certificate/verify_n3_certificate.py
  .venv/bin/python experiments/grid_resolution_study.py
}

run_full_numerics() {
  .venv/bin/python experiments/search_switches.py --check-n3 --grid 301
  .venv/bin/python experiments/linear_family_search.py \
    --m 4 --L 16 --epsilon 1/30 --roots 1/3,5/12 --grid 200001 \
    --output experiments/explicit_family_m4_L16.json
  .venv/bin/python experiments/linear_family_search.py \
    --m 5 --L 16 --epsilon 1/30 --roots 1/3,5/12,9/20 --grid 200001 \
    --output experiments/explicit_family_m5_L16.json
  .venv/bin/python experiments/linear_family_search.py \
    --m 6 --L 16 --epsilon 1/30 --roots 1/3,3/8,5/12,9/20 --grid 200001 \
    --output experiments/explicit_family_m6_L16.json
  .venv/bin/python experiments/check_linear_family_decimal.py \
    experiments/explicit_family_m4_L16.json \
    experiments/explicit_family_m4_L16_decimal90.json
  .venv/bin/python experiments/check_linear_family_decimal.py \
    experiments/explicit_family_m5_L16.json \
    experiments/explicit_family_m5_L16_decimal90.json
  .venv/bin/python experiments/check_linear_family_decimal.py \
    experiments/explicit_family_m6_L16.json \
    experiments/explicit_family_m6_L16_decimal90.json
  .venv/bin/python experiments/grid_resolution_study.py
  run_figure
}

mode=${1:---quick}

case "$mode" in
  --certificate)
    .venv/bin/python certificate/verify_n3_certificate.py
    ;;
  --quick)
    run_quick
    ;;
  --figure)
    run_figure
    ;;
  --full)
    .venv/bin/python certificate/verify_n3_certificate.py
    .venv/bin/python certificate/independent/rational_interval_check.py
    run_full_numerics
    ;;
  *)
    echo "usage: ./reproduce.sh [--certificate|--quick|--figure|--full]" >&2
    exit 2
    ;;
esac
