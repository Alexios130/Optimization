#!/usr/bin/env python3
"""Instantiate the analytic linear-switch family at finite rational scales.

Calculations are performed in the stable bounded-feature variables

    a_j = 1+j/L,  lambda_L = 2^{-L}/L,  r in [0,1/2].

For L=4q^2 this path is exactly equivalent, via v=C_L beta and R=C_L r,
to the requested fixed-lambda=1 rational instance with
C_L=2q*2^(L/2) and features A_j=C_L(1+j/L).

This script supplies dense floating-point/Brent evidence only.  It is not an
interval certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


def multiply_linear(coeffs: list[Fraction], root: Fraction) -> list[Fraction]:
    out = [Fraction(0)] * (len(coeffs) + 1)
    for j, value in enumerate(coeffs):
        out[j] -= root * value
        out[j + 1] += value
    return out


def polynomial_coefficients(roots: list[Fraction]) -> list[Fraction]:
    coeffs = [Fraction(1)]
    for z in [Fraction(1)] + roots:
        coeffs = multiply_linear(coeffs, z)
    return coeffs


def softplus_negative(z: np.ndarray) -> np.ndarray:
    # z is positive and often large; log1p(exp(-z)) is stable here.
    return np.log1p(np.exp(-z))


class ScaledFamily:
    def __init__(self, m: int, roots: list[Fraction], eps: Fraction, L: int):
        self.m, self.roots, self.eps, self.L = m, roots, eps, L
        self.c = polynomial_coefficients(roots)
        if len(self.c) != m:
            raise ValueError("need exactly m-2 prescribed roots")
        self.w_fraction = [
            [Fraction(1, m) + eps * cj for cj in self.c],
            [Fraction(1, m) - eps * cj for cj in self.c],
        ]
        if min(min(row) for row in self.w_fraction) <= 0:
            raise ValueError("epsilon makes a weight nonpositive")
        self.w = np.asarray([[float(v) for v in row] for row in self.w_fraction])
        self.D = self.w[0] - self.w[1]
        self.a = 1.0 + np.arange(m, dtype=float) / L
        self.lam = math.ldexp(1.0 / L, -L)

    def branch(self, scenario: int, r: float) -> float:
        d = self.a - r

        def grad(v: float) -> float:
            z = v * d
            # sigmoid(-z), stable for z>=0
            em = np.exp(-z)
            sigm = em / (1.0 + em)
            return self.lam * v - float(np.dot(self.w[scenario], d * sigm))

        hi = max(8.0, 4.0 * self.L * math.log(2.0))
        while grad(hi) <= 0.0:
            hi *= 2.0
        return float(brentq(grad, 0.0, hi, xtol=2e-12, rtol=2e-14, maxiter=200))

    def gaps(self, r: float) -> np.ndarray:
        out = np.empty(2)
        d = self.a - r
        for s in range(2):
            v = self.branch(s, r)
            out[s] = float(np.dot(self.D, softplus_negative(v * d)))
        return out

    def scan(self, grid_size: int = 40001) -> dict:
        grid = np.linspace(0.0, 0.5, grid_size)
        vals = np.empty((grid_size, 2))
        for k, r in enumerate(grid):
            vals[k] = self.gaps(float(r))
        raw: list[tuple[float, int]] = []
        for s in range(2):
            for k in range(grid_size - 1):
                if vals[k, s] * vals[k + 1, s] < 0.0:
                    root = brentq(lambda rr: self.gaps(rr)[s], grid[k], grid[k + 1], xtol=2e-13, rtol=2e-13)
                    raw.append((float(root), s + 1))
        raw.sort()
        events: list[tuple[float, int]] = []
        for event in raw:
            if not events or abs(event[0] - events[-1][0]) > 1e-10 or event[1] != events[-1][1]:
                events.append(event)

        edges = [0.0] + [r for r, _ in events] + [0.5]
        states = []
        probes = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            r = 0.5 * (lo + hi)
            g = self.gaps(r)
            state = "1" if g[0] > 0 else ("2" if g[1] < 0 else "12")
            states.append(state)
            probes.append({"r": r, "gaps": g.tolist(), "state": state})

        q = int(round(math.sqrt(self.L / 4)))
        if 4 * q * q != self.L:
            raise ValueError("fixed-lambda scaling requires L=4q^2")
        C = 2 * q * (1 << (self.L // 2))
        feature_fractions = [Fraction(C * (self.L + j), self.L) for j in range(self.m)]
        observation_features = []
        labels = []
        scenario_observation_weights = [[], []]
        for j, A in enumerate(feature_fractions):
            observation_features.extend([str(A), str(-A)])
            labels.extend([1, -1])
            for s in range(2):
                half = self.w_fraction[s][j] / 2
                scenario_observation_weights[s].extend([str(half), str(half)])

        return {
            "status": "DENSE_NUMERICAL_EVIDENCE_NOT_INTERVAL_CERTIFICATION",
            "m": self.m,
            "n": 2 * self.m,
            "L": self.L,
            "q": q,
            "epsilon": str(self.eps),
            "prescribed_z_roots": [str(z) for z in self.roots],
            "polynomial_coefficients_ascending": [str(v) for v in self.c],
            "pair_weights": [[str(v) for v in row] for row in self.w_fraction],
            "scaled_computation": {
                "features": [str(Fraction(self.L + j, self.L)) for j in range(self.m)],
                "lambda": str(Fraction(1, self.L * (1 << self.L))),
                "r_interval": [0.0, 0.5],
                "grid_points": grid_size,
            },
            "fixed_lambda_one_instance": {
                "C_L": str(C),
                "features": observation_features,
                "labels": labels,
                "scenario_weights": scenario_observation_weights,
                "lambda": "1",
                "R_interval": ["0", str(Fraction(C, 2))],
            },
            "K_detected": len(events),
            "predicted_lower_bound": 2 * (self.m - 2),
            "event_radii_scaled": [r for r, _ in events],
            "event_radii_fixed_lambda": [float(C * r) for r, _ in events],
            "event_branches": [s for _, s in events],
            "interval_states": states,
            "interval_probes": probes,
            "minimum_scaled_event_separation": min(np.diff([r for r, _ in events])) if len(events) > 1 else None,
        }


DEFAULT_ROOTS = {
    4: [Fraction(3, 10), Fraction(2, 5)],
    5: [Fraction(7, 25), Fraction(1, 3), Fraction(2, 5)],
    6: [Fraction(7, 25), Fraction(3, 10), Fraction(1, 3), Fraction(2, 5)],
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, choices=[4, 5, 6], required=True)
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--epsilon", default="1/100")
    ap.add_argument("--roots", help="comma-separated rational z roots; defaults to report choices")
    ap.add_argument("--grid", type=int, default=40001)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    roots = DEFAULT_ROOTS[args.m] if not args.roots else [Fraction(v) for v in args.roots.split(",")]
    result = ScaledFamily(args.m, roots, Fraction(args.epsilon), args.L).scan(args.grid)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(
        f"m={args.m} L={args.L} K={result['K_detected']} "
        f"states={result['interval_states']} sep={result['minimum_scaled_event_separation']}"
    )


if __name__ == "__main__":
    main()
