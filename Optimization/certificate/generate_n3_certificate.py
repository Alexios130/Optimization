#!/usr/bin/env python3
"""Generate rational boxes for the rigorous n=3 path certificate.

This script is a *generator*, not the verifier.  Floating-point SciPy solves
are used only to propose rational centres and boxes.  The companion script
``verify_n3_certificate.py`` recomputes every mathematical inclusion with
Arb ball arithmetic.  Consequently a poor floating-point proposal can make
generation or verification fail, but it cannot make an invalid certificate
pass.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import root


X = np.array([3.5, 3.0, 4.0])
Y = np.array([-1.0, 1.0, -1.0])
P = np.array(
    [
        [12.0 / 19.0, 1.0 / 19.0, 6.0 / 19.0],
        [1.0 / 8.0, 1.0 / 16.0, 13.0 / 16.0],
    ]
)
LAM = 0.5


def sigmoid(q: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-q))


def gradient(scenario: int, theta: np.ndarray, rho: float) -> np.ndarray:
    beta, intercept = theta
    q = -Y * (X * beta + intercept) - rho * beta
    a = sigmoid(q)
    q_beta = -Y * X - rho
    return np.array(
        [
            LAM * beta + np.sum(P[scenario] * a * q_beta),
            LAM * intercept + np.sum(P[scenario] * a * (-Y)),
        ]
    )


def solve_single(scenario: int, rho: float, guess: np.ndarray | None = None) -> np.ndarray:
    if guess is None:
        guess = np.array([-0.50, -0.20])
    ans = root(lambda z: gradient(scenario, z, rho), guess, tol=1.0e-12)
    if not ans.success and np.linalg.norm(gradient(scenario, ans.x, rho), ord=np.inf) > 1.0e-10:
        raise RuntimeError(f"single-scenario solve failed at s={scenario + 1}, rho={rho}")
    return ans.x


def decimal(value: float) -> str:
    """A short decimal that denotes an exact rational in the verifier."""
    return format(float(value), ".17g")


def proposed_path_box(scenario: int, lo: float, hi: float) -> dict[str, str | int]:
    mid = 0.5 * (lo + hi)
    theta_mid = solve_single(scenario, mid)
    theta_lo = solve_single(scenario, lo, theta_mid)
    theta_hi = solve_single(scenario, hi, theta_mid)

    beta_radius = 1.35 * max(abs(theta_lo[0] - theta_mid[0]), abs(theta_hi[0] - theta_mid[0]))
    b_radius = 1.35 * max(abs(theta_lo[1] - theta_mid[1]), abs(theta_hi[1] - theta_mid[1]))
    beta_radius = max(beta_radius + 2.0e-12, 2.0e-11)
    b_radius = max(b_radius + 2.0e-12, 2.0e-11)

    return {
        "scenario": scenario + 1,
        "rho_lo": decimal(lo),
        "rho_hi": decimal(hi),
        "beta_mid": decimal(theta_mid[0]),
        "b_mid": decimal(theta_mid[1]),
        "beta_radius": decimal(beta_radius),
        "b_radius": decimal(b_radius),
    }


def subdivide(lo: float, hi: float, step: float) -> list[tuple[float, float]]:
    count = max(1, math.ceil((hi - lo) / step))
    edges = np.linspace(lo, hi, count + 1)
    return [(float(edges[k]), float(edges[k + 1])) for k in range(count)]


def path_cells(scenario: int, lo: float, hi: float, step: float, purpose: str, sign: int = 0):
    cells = []
    for a, b in subdivide(lo, hi, step):
        cell = proposed_path_box(scenario, a, b)
        cell["purpose"] = purpose
        if sign:
            cell["sign"] = sign
        cells.append(cell)
    return cells


def main(output: Path) -> None:
    # Wide boxes are used for exhaustive root isolation.  Interval Newton
    # contracts them to much narrower certified enclosures in the verifier.
    roots = [
        {
            "name": "rho_1",
            "scenario": 2,
            "derivative_sign": 1,
            "report_rho_lo": "0.5245653271531518",
            "report_rho_hi": "0.5245653271531519",
            "report_derivative_lo": "0.0029355",
            "report_derivative_hi": "0.0029357",
            "beta_mid": "-0.49466473910219629",
            "b_mid": "-0.14980130042270346",
            "rho_mid": "0.52456532715315181",
            "beta_radius": "0.00001",
            "b_radius": "0.00001",
            "rho_radius": "0.00001",
        },
        {
            "name": "rho_2",
            "scenario": 1,
            "derivative_sign": 1,
            "report_rho_lo": "1.0369491958728015",
            "report_rho_hi": "1.0369491958728016",
            "report_derivative_lo": "0.0040839",
            "report_derivative_hi": "0.0040841",
            "beta_mid": "-0.49312133891920980",
            "b_mid": "-0.24118658957399815",
            "rho_mid": "1.03694919587280154",
            "beta_radius": "0.00001",
            "b_radius": "0.00001",
            "rho_radius": "0.00001",
        },
        {
            "name": "rho_3",
            "scenario": 1,
            "derivative_sign": -1,
            "report_rho_lo": "2.6961208910188072",
            "report_rho_hi": "2.6961208910188073",
            "report_derivative_lo": "-0.0133944",
            "report_derivative_hi": "-0.0133942",
            "beta_mid": "-0.13802252132535599",
            "b_mid": "-0.55151056569378579",
            "rho_mid": "2.69612089101880727",
            "beta_radius": "0.00001",
            "b_radius": "0.00001",
            "rho_radius": "0.00001",
        },
        {
            "name": "rho_4",
            "scenario": 2,
            "derivative_sign": -1,
            "report_rho_lo": "2.8324298315608840",
            "report_rho_hi": "2.8324298315608841",
            "report_derivative_lo": "-0.0139753",
            "report_derivative_hi": "-0.0139750",
            "beta_mid": "-0.12604252405635794",
            "b_mid": "-0.53537389564131285",
            "rho_mid": "2.83242983156088409",
            "beta_radius": "0.00001",
            "b_radius": "0.00001",
            "rho_radius": "0.00001",
        },
    ]

    cells: list[dict[str, str | int]] = []

    # The derivative windows deliberately leave a short neighbourhood of
    # each positive maximum.  Those neighbourhoods are certified directly
    # by the ``positive`` cells below.
    cells += path_cells(0, 0.0, 2.07, 0.00025, "derivative", +1)
    cells += path_cells(0, 2.07, 2.12, 0.00050, "positive")
    cells += path_cells(0, 2.12, 2.85, 0.00025, "derivative", -1)

    cells += path_cells(1, 0.0, 2.05, 0.00025, "derivative", +1)
    cells += path_cells(1, 2.05, 2.10, 0.00050, "positive")
    cells += path_cells(1, 2.10, 3.1890, 0.00025, "derivative", -1)

    # The final 8e-4 before the nonsmooth zero-slope threshold is handled by
    # short negative-gap cells.  Their verification uses the same validated
    # path boxes but does not require differentiability at beta=0.
    cells += path_cells(1, 3.1890, 3.18970, 0.00005, "negative")

    certificate = {
        "format": "n3-path-certificate-v1",
        "precision_bits": 256,
        "model": {
            "x": [[7, 2], [3, 1], [4, 1]],
            "y": [-1, 1, -1],
            "pi": [
                [[12, 19], [1, 19], [6, 19]],
                [[1, 8], [1, 16], [13, 16]],
            ],
            "lambda": [1, 2],
        },
        "root_boxes": roots,
        "path_cells": cells,
        "intercept_root": {
            "b_mid": "-0.58881403292955227",
            "b_radius": "0.00000001",
            "prethreshold_lower": "3.18970",
            "report_b_lo": "-0.5888140329295523",
            "report_b_hi": "-0.5888140329295522",
            "report_rho_infty_lo": "3.1897944978297369",
            "report_rho_infty_hi": "3.1897944978297370",
            "report_tail_gap_upper": "-0.0058106647986468",
            "report_prethreshold_gap_upper": "-0.0053081614236373",
        },
    }

    output.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output} with {len(cells)} path-cell proposals (not yet validated)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Propose a new n=3 certificate ledger (not a proof)."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("n3_certificate.generated.json"),
        help="output path (default: n3_certificate.generated.json)",
    )
    args = parser.parse_args()
    main(args.output)
