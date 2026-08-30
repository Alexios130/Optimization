#!/usr/bin/env python3
"""High-precision Decimal Newton check for finite members of the linear family.

The input JSON must come from ``linear_family_search.py``.  Exact rational
weights/features are reconstructed, each proposed event is solved as the
two-equation system (branch stationarity, branch gap=0), and transversality
is evaluated.  Decimal arithmetic is high-precision numerical evidence, not
outward-rounded interval certification.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path


def dec_fraction(text: str) -> Decimal:
    f = Fraction(text)
    return Decimal(f.numerator) / Decimal(f.denominator)


def solve2(A, rhs):
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    return [
        (rhs[0] * A[1][1] - A[0][1] * rhs[1]) / det,
        (A[0][0] * rhs[1] - rhs[0] * A[1][0]) / det,
    ]


class DecimalFamily:
    def __init__(self, raw: dict):
        getcontext().prec = 90
        self.raw = raw
        self.m = raw["m"]
        self.L = raw["L"]
        self.a = [Decimal(1) + Decimal(j) / Decimal(self.L) for j in range(self.m)]
        self.w = [[dec_fraction(v) for v in row] for row in raw["pair_weights"]]
        self.D = [self.w[0][j] - self.w[1][j] for j in range(self.m)]
        self.lam = Decimal(1) / (Decimal(self.L) * (Decimal(2) ** self.L))

    @staticmethod
    def sig_minus(z: Decimal) -> Decimal:
        return Decimal(1) / (Decimal(1) + z.exp())

    def equations(self, scenario: int, v: Decimal, r: Decimal):
        d = [a - r for a in self.a]
        sig = [self.sig_minus(v * dj) for dj in d]
        ap = [s * (Decimal(1) - s) for s in sig]
        grad = self.lam * v - sum(self.w[scenario][j] * d[j] * sig[j] for j in range(self.m))
        gap = sum(
            self.D[j] * (Decimal(1) + (-v * d[j]).exp()).ln() for j in range(self.m)
        )
        grad_v = self.lam + sum(self.w[scenario][j] * d[j] * d[j] * ap[j] for j in range(self.m))
        grad_r = sum(self.w[scenario][j] * (sig[j] - v * d[j] * ap[j]) for j in range(self.m))
        gap_v = -sum(self.D[j] * d[j] * sig[j] for j in range(self.m))
        gap_r = v * sum(self.D[j] * sig[j] for j in range(self.m))
        return [grad, gap], [[grad_v, grad_r], [gap_v, gap_r]]

    def solve_v(self, scenario: int, r: Decimal) -> Decimal:
        lo = Decimal(0)
        hi = Decimal(4 * self.L)
        for _ in range(320):
            mid = (lo + hi) / 2
            grad = self.equations(scenario, mid, r)[0][0]
            if grad < 0:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    def probe(self, r0: float):
        r = Decimal(str(r0))
        gaps = []
        for scenario in range(2):
            v = self.solve_v(scenario, r)
            gaps.append(self.equations(scenario, v, r)[0][1])
        state = "1" if gaps[0] > 0 else ("2" if gaps[1] < 0 else "12")
        return {
            "r": format(r, ".18g"),
            "gap_branch_1": format(gaps[0], ".35g"),
            "gap_branch_2": format(gaps[1], ".35g"),
            "state": state,
        }

    def solve_event(self, scenario: int, r0: float):
        # Double-asymptotic initialization is already close; a scalar Newton
        # solve at fixed r provides the branch v seed.
        r = Decimal(str(r0))
        v = self.solve_v(scenario, r)
        for _ in range(40):
            F, J = self.equations(scenario, v, r)
            dv, dr = solve2(J, [-F[0], -F[1]])
            v += dv
            r += dr
            if max(abs(F[0]), abs(F[1])) < Decimal("1e-75"):
                break
        F, J = self.equations(scenario, v, r)
        dv_dr = -J[0][1] / J[0][0]
        total_gap_derivative = J[1][1] + J[1][0] * dv_dr
        return {
            "branch": scenario + 1,
            "r": format(r, ".75g"),
            "v": format(v, ".75g"),
            "max_residual": format(max(abs(F[0]), abs(F[1])), ".8g"),
            "gap_derivative_dr": format(total_gap_derivative, ".30g"),
        }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    raw = json.loads(args.input.read_text(encoding="utf-8"))
    family = DecimalFamily(raw)
    events = [
        family.solve_event(branch - 1, r)
        for r, branch in zip(raw["event_radii_scaled"], raw["event_branches"])
    ]
    high_r = [Decimal(e["r"]) for e in events]
    separations = [high_r[j + 1] - high_r[j] for j in range(len(high_r) - 1)]
    result = {
        "status": "90_DIGIT_NUMERICAL_EVIDENCE_NOT_INTERVAL_CERTIFICATION",
        "source": args.input.name,
        "events": events,
        "interval_probes": [family.probe(item["r"]) for item in raw["interval_probes"]],
        "minimum_event_separation_r": format(min(separations), ".40g"),
        "minimum_absolute_event_derivative": format(
            min(abs(Decimal(e["gap_derivative_dr"])) for e in events), ".30g"
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        f"{args.input.name}: events={len(events)} minsep={result['minimum_event_separation_r']} "
        f"minderiv={result['minimum_absolute_event_derivative']}"
    )


if __name__ == "__main__":
    main()
