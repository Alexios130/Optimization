#!/usr/bin/env python3
"""Verify the complete n=3 path certificate with Arb ball arithmetic.

Only rational input data and rational decimal boxes are trusted.  Every
transcendental evaluation is carried out by Arb with outward rounding.
The verifier checks:

* interval-Newton uniqueness for the four endpoint systems;
* validated enclosures of each single-scenario solution branch;
* the signs of the two branch-gap derivatives outside small maximum windows;
* strict positivity of both branch gaps inside those windows;
* the exact zero-slope tail threshold and its strict active scenario.

A successful run writes a compact JSON summary.  Failure of any inclusion or
strict sign terminates the run with an assertion error.
"""

from __future__ import annotations

import json
import hashlib
from decimal import Decimal
from pathlib import Path

from flint import arb, arb_mat, ctx


HERE = Path(__file__).resolve().parent
INPUT = HERE / "n3_certificate.json"
OUTPUT = HERE / "n3_certificate_summary.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "8218a69bcce225540785109591f3b8db6dc71892a5f508b491731e85a03fcfb2"
)

EXPECTED_MODEL = {
    "x": [[7, 2], [3, 1], [4, 1]],
    "y": [-1, 1, -1],
    "pi": [
        [[12, 19], [1, 19], [6, 19]],
        [[1, 8], [1, 16], [13, 16]],
    ],
    "lambda": [1, 2],
}


def ball(mid: str, radius: str) -> arb:
    return arb(mid, radius)


def interval(lo: str, hi: str) -> arb:
    lower = arb(lo)
    upper = arb(hi)
    return arb((lower + upper) / 2, (upper - lower) / 2)


def sigmoid(q: arb) -> arb:
    return 1 / (1 + (-q).exp())


def softplus(q: arb) -> arb:
    return (1 + q.exp()).log()


def abs_upper_float(x: arb) -> float:
    return float(x.abs_upper())


class Model:
    def __init__(self, raw: dict):
        self.x = [arb(a) / b for a, b in raw["x"]]
        self.y = [arb(v) for v in raw["y"]]
        self.pi = [[arb(a) / b for a, b in row] for row in raw["pi"]]
        self.lam = arb(raw["lambda"][0]) / raw["lambda"][1]
        self.delta = [self.pi[0][i] - self.pi[1][i] for i in range(3)]

    def branch_quantities(self, scenario: int, beta: arb, b: arb, rho: arb):
        g_beta = self.lam * beta
        g_b = self.lam * b
        gap = arb(0)

        a_values = []
        w_values = []
        c_values = []
        d_values = []
        for i in range(3):
            q = -self.y[i] * (self.x[i] * beta + b) - rho * beta
            a = sigmoid(q)
            w = a * (1 - a)
            c = -self.y[i] * self.x[i] - rho
            d = -self.y[i]
            g_beta += self.pi[scenario][i] * a * c
            g_b += self.pi[scenario][i] * a * d
            gap += self.delta[i] * softplus(q)
            a_values.append(a)
            w_values.append(w)
            c_values.append(c)
            d_values.append(d)

        return g_beta, g_b, gap, a_values, w_values, c_values, d_values

    def event_function(self, scenario: int, z: list[arb]) -> arb_mat:
        g_beta, g_b, gap, *_ = self.branch_quantities(scenario, z[0], z[1], z[2])
        return arb_mat([[g_beta], [g_b], [gap]])

    def event_jacobian(self, scenario: int, z: list[arb]) -> arb_mat:
        beta, b, rho = z
        _, _, _, aa, ww, cc, dd = self.branch_quantities(scenario, beta, b, rho)

        rows = [[arb(0) for _ in range(3)] for _ in range(3)]
        rows[0][0] = self.lam
        rows[1][1] = self.lam
        for i in range(3):
            p = self.pi[scenario][i]
            rows[0][0] += p * ww[i] * cc[i] * cc[i]
            rows[0][1] += p * ww[i] * cc[i] * dd[i]
            rows[0][2] += p * (-beta * ww[i] * cc[i] - aa[i])
            rows[1][0] += p * ww[i] * cc[i] * dd[i]
            rows[1][1] += p * ww[i] * dd[i] * dd[i]
            rows[1][2] += p * (-beta * ww[i] * dd[i])
            rows[2][0] += self.delta[i] * aa[i] * cc[i]
            rows[2][1] += self.delta[i] * aa[i] * dd[i]
            rows[2][2] += self.delta[i] * aa[i] * (-beta)
        return arb_mat(rows)

    def path_derivative(self, scenario: int, beta: arb, b: arb, rho: arb):
        _, _, _, aa, ww, cc, dd = self.branch_quantities(scenario, beta, b, rho)

        h00 = self.lam
        h01 = arb(0)
        h11 = self.lam
        gr0 = arb(0)
        gr1 = arb(0)
        gap_beta = arb(0)
        gap_b = arb(0)
        gap_rho = arb(0)
        for i in range(3):
            p = self.pi[scenario][i]
            h00 += p * ww[i] * cc[i] * cc[i]
            h01 += p * ww[i] * cc[i] * dd[i]
            h11 += p * ww[i] * dd[i] * dd[i]
            gr0 += p * (ww[i] * (-beta) * cc[i] - aa[i])
            gr1 += p * ww[i] * (-beta) * dd[i]
            gap_beta += self.delta[i] * aa[i] * cc[i]
            gap_b += self.delta[i] * aa[i] * dd[i]
            gap_rho += self.delta[i] * aa[i] * (-beta)

        theta_prime = -arb_mat([[h00, h01], [h01, h11]]).inv() * arb_mat([[gr0], [gr1]])
        gap_prime = gap_beta * theta_prime[0, 0] + gap_b * theta_prime[1, 0] + gap_rho
        return theta_prime[0, 0], theta_prime[1, 0], gap_prime

    def gap_lipschitz(self, rho_upper: arb) -> arb:
        result = arb(0)
        for i in range(3):
            result += abs(self.delta[i]) * (((abs(self.x[i]) + rho_upper) ** 2 + 1).sqrt())
        return result


def verify_root_box(model: Model, item: dict) -> dict:
    scenario = int(item["scenario"]) - 1
    assert Decimal(item["beta_radius"]) > 0
    assert Decimal(item["b_radius"]) > 0
    assert Decimal(item["rho_radius"]) > 0
    assert (
        Decimal(item["rho_mid"]) - Decimal(item["rho_radius"])
        >= 0
    )
    assert Decimal(item["report_rho_lo"]) < Decimal(item["report_rho_hi"])
    assert (
        Decimal(item["report_derivative_lo"])
        < Decimal(item["report_derivative_hi"])
    )
    box = [
        ball(item["beta_mid"], item["beta_radius"]),
        ball(item["b_mid"], item["b_radius"]),
        ball(item["rho_mid"], item["rho_radius"]),
    ]
    assert box[0].upper() < 0
    midpoint = [v.mid() for v in box]
    newton = arb_mat([[v] for v in midpoint]) - model.event_jacobian(scenario, box).inv() * model.event_function(scenario, midpoint)
    assert all(box[i].contains_interior(newton[i, 0]) for i in range(3)), item["name"]

    # Every zero of the branch gap in the rho projection must use the unique
    # branch minimizer.  Strong convexity encloses that minimizer inside the
    # beta/beta box, making the Newton uniqueness exhaustive for this radius
    # interval rather than merely local.
    rho_box = box[2]
    g0, g1, *_ = model.branch_quantities(scenario, midpoint[0], midpoint[1], rho_box)
    eps = (g0.abs_upper() ** 2 + g1.abs_upper() ** 2).sqrt()
    distance = eps / model.lam
    assert midpoint[0] - distance > box[0].lower()
    assert midpoint[0] + distance < box[0].upper()
    assert midpoint[1] - distance > box[1].lower()
    assert midpoint[1] + distance < box[1].upper()

    contracted = [newton[i, 0].intersection(box[i]) for i in range(3)]
    report_rho_lo = arb(item["report_rho_lo"])
    report_rho_hi = arb(item["report_rho_hi"])
    assert report_rho_lo.upper() < contracted[2].lower()
    assert contracted[2].upper() < report_rho_hi.lower()

    _, _, gap_prime = model.path_derivative(
        scenario, contracted[0], contracted[1], contracted[2]
    )
    report_derivative_lo = arb(item["report_derivative_lo"])
    report_derivative_hi = arb(item["report_derivative_hi"])
    assert report_derivative_lo.upper() < gap_prime.lower()
    assert gap_prime.upper() < report_derivative_hi.lower()
    if int(item["derivative_sign"]) > 0:
        assert gap_prime.lower() > 0
    else:
        assert gap_prime.upper() < 0

    return {
        "name": item["name"],
        "scenario": scenario + 1,
        "rho_lower": str(contracted[2].lower()),
        "rho_upper": str(contracted[2].upper()),
        "beta_lower": str(contracted[0].lower()),
        "beta_upper": str(contracted[0].upper()),
        "b_lower": str(contracted[1].lower()),
        "b_upper": str(contracted[1].upper()),
        "report_rho_lo": item["report_rho_lo"],
        "report_rho_hi": item["report_rho_hi"],
        "report_derivative_lo": item["report_derivative_lo"],
        "report_derivative_hi": item["report_derivative_hi"],
        "gap_derivative_lower": str(gap_prime.lower()),
        "gap_derivative_upper": str(gap_prime.upper()),
    }


def verify_path_cell(model: Model, item: dict) -> dict:
    scenario = int(item["scenario"]) - 1
    rho_lo_decimal = Decimal(item["rho_lo"])
    rho_hi_decimal = Decimal(item["rho_hi"])
    assert 0 <= rho_lo_decimal < rho_hi_decimal
    assert Decimal(item["beta_radius"]) > 0
    assert Decimal(item["b_radius"]) > 0
    rho = interval(item["rho_lo"], item["rho_hi"])
    rho_mid = (arb(item["rho_lo"]) + arb(item["rho_hi"])) / 2
    half_width = (arb(item["rho_hi"]) - arb(item["rho_lo"])) / 2
    beta_mid = arb(item["beta_mid"])
    b_mid = arb(item["b_mid"])
    beta_box = ball(item["beta_mid"], item["beta_radius"])
    b_box = ball(item["b_mid"], item["b_radius"])
    assert beta_box.upper() < 0

    # The exact midpoint solution is enclosed by strong convexity.
    g0, g1, gap_mid, *_ = model.branch_quantities(scenario, beta_mid, b_mid, rho_mid)
    eps = (g0.abs_upper() ** 2 + g1.abs_upper() ** 2).sqrt()
    midpoint_distance = eps / model.lam

    beta_prime, b_prime, gap_prime = model.path_derivative(scenario, beta_box, b_box, rho)
    assert midpoint_distance + beta_prime.abs_upper() * half_width < arb(item["beta_radius"])
    assert midpoint_distance + b_prime.abs_upper() * half_width < arb(item["b_radius"])

    purpose = item["purpose"]
    if purpose == "derivative":
        wanted = int(item["sign"])
        if wanted > 0:
            assert gap_prime.lower() > 0
        else:
            assert gap_prime.upper() < 0
    elif purpose == "positive":
        # First certify H_s at the exact midpoint minimizer, then integrate
        # the validated derivative enclosure across the cell.
        k = model.gap_lipschitz(rho_mid)
        corrected_mid_lower = gap_mid.lower() - (k * midpoint_distance).upper()
        assert corrected_mid_lower - gap_prime.abs_upper() * half_width > 0
    elif purpose == "negative":
        # Directly enclose the branch gap on the whole validated path box.
        _, _, gap_box, *_ = model.branch_quantities(scenario, beta_box, b_box, rho)
        assert gap_box.upper() < 0
    else:
        raise AssertionError(f"unknown cell purpose {purpose}")

    return {
        "scenario": scenario + 1,
        "rho_lo": item["rho_lo"],
        "rho_hi": item["rho_hi"],
        "purpose": purpose,
        "gap_prime_lower": str(gap_prime.lower()),
        "gap_prime_upper": str(gap_prime.upper()),
    }


def verify_intercept_tail(model: Model, raw: dict) -> dict:
    assert Decimal(raw["b_radius"]) > 0
    assert Decimal(raw["prethreshold_lower"]) >= 0
    assert Decimal(raw["report_b_lo"]) < Decimal(raw["report_b_hi"])
    assert (
        Decimal(raw["report_rho_infty_lo"])
        < Decimal(raw["report_rho_infty_hi"])
    )
    b_box = ball(raw["b_mid"], raw["b_radius"])
    b_mid = b_box.mid()
    scenario = 1  # zero-based scenario 2

    def intercept_gradient(b: arb) -> arb:
        total = model.lam * b
        for i in range(3):
            a = sigmoid(-model.y[i] * b)
            total -= model.pi[scenario][i] * a * model.y[i]
        return total

    def intercept_hessian(b: arb) -> arb:
        total = model.lam
        for i in range(3):
            a = sigmoid(-model.y[i] * b)
            total += model.pi[scenario][i] * a * (1 - a)
        return total

    newton = b_mid - intercept_gradient(b_mid) / intercept_hessian(b_box)
    assert b_box.contains_interior(newton)
    b_exact = newton.intersection(b_box)

    c = arb(0)
    h = arb(0)
    gap = arb(0)
    for i in range(3):
        q = -model.y[i] * b_exact
        a = sigmoid(q)
        c -= model.pi[scenario][i] * a * model.y[i] * model.x[i]
        h += model.pi[scenario][i] * a
        gap += model.delta[i] * softplus(q)
    assert c.lower() > 0 and h.lower() > 0 and gap.upper() < 0
    threshold = c / h
    assert arb(raw["report_b_lo"]).upper() < b_exact.lower()
    assert b_exact.upper() < arb(raw["report_b_hi"]).lower()
    assert (
        arb(raw["report_rho_infty_lo"]).upper()
        < threshold.lower()
    )
    assert (
        threshold.upper()
        < arb(raw["report_rho_infty_hi"]).lower()
    )
    assert gap.upper() < arb(raw["report_tail_gap_upper"]).lower()

    # Just below the threshold, (0,b_infty) is not yet the exact minimizer,
    # but its smallest-norm subgradient is explicit.  Strong convexity and
    # the global gap Lipschitz bound transfer the strict scenario-2 gap to
    # the exact single-scenario minimizer on the whole remaining interval.
    rho_low = arb(raw["prethreshold_lower"])
    slope_residual = c - rho_low * h
    assert slope_residual.lower() > 0
    distance = slope_residual / model.lam
    correction = model.gap_lipschitz(threshold.upper()) * distance
    assert gap.upper() + correction.upper() < 0
    assert (
        gap.upper() + correction.upper()
        < arb(raw["report_prethreshold_gap_upper"]).lower()
    )

    return {
        "b_lower": str(b_exact.lower()),
        "b_upper": str(b_exact.upper()),
        "rho_infty_lower": str(threshold.lower()),
        "rho_infty_upper": str(threshold.upper()),
        "tail_gap_lower": str(gap.lower()),
        "tail_gap_upper": str(gap.upper()),
        "prethreshold_lower": raw["prethreshold_lower"],
        "prethreshold_gap_upper": str(gap.upper() + correction.upper()),
        "report_b_lo": raw["report_b_lo"],
        "report_b_hi": raw["report_b_hi"],
        "report_rho_infty_lo": raw["report_rho_infty_lo"],
        "report_rho_infty_hi": raw["report_rho_infty_hi"],
        "report_tail_gap_upper": raw["report_tail_gap_upper"],
        "report_prethreshold_gap_upper": raw[
            "report_prethreshold_gap_upper"
        ],
    }


def main() -> None:
    if not __debug__:
        raise RuntimeError(
            "certificate verification requires assertions; do not use "
            "python -O or PYTHONOPTIMIZE"
        )
    OUTPUT.write_text(
        json.dumps(
            {
                "status": "INCOMPLETE",
                "message": "verification did not finish successfully",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    input_bytes = INPUT.read_bytes()
    certificate_sha256 = hashlib.sha256(input_bytes).hexdigest()
    assert certificate_sha256 == EXPECTED_CERTIFICATE_SHA256
    raw = json.loads(input_bytes)
    assert raw["format"] == "n3-path-certificate-v1"
    assert int(raw["precision_bits"]) == 256
    assert raw["model"] == EXPECTED_MODEL
    assert len(raw["path_cells"]) == 23973
    ctx.prec = int(raw["precision_bits"])
    model = Model(raw["model"])

    roots = [verify_root_box(model, item) for item in raw["root_boxes"]]
    cells = [verify_path_cell(model, item) for item in raw["path_cells"]]
    tail = verify_intercept_tail(model, raw["intercept_root"])

    # The finitely many slabs must form genuine covers, not merely a list of
    # independently valid local statements.  Cells were generated in path
    # order separately for scenarios 1 and 2.
    coverage_blocks: dict[int, list[dict]] = {}
    for scenario, expected_end in ((1, "2.85"), (2, "3.18970")):
        group = [item for item in raw["path_cells"] if int(item["scenario"]) == scenario]
        assert Decimal(group[0]["rho_lo"]) == 0
        for left, right in zip(group, group[1:]):
            # These are decimal encodings of the same rational endpoint.  A
            # string comparison is intentional: constructing two Arb balls
            # from a non-dyadic decimal gives overlapping balls, for which
            # Arb's equality operator need not return true.
            assert left["rho_hi"] == right["rho_lo"]
        assert Decimal(group[-1]["rho_hi"]) >= Decimal(expected_end)

        blocks = []
        for item in group:
            key = (
                item["purpose"],
                int(item["sign"]) if item["purpose"] == "derivative" else None,
            )
            if not blocks or blocks[-1]["key"] != key:
                blocks.append(
                    {
                        "key": key,
                        "rho_lo": item["rho_lo"],
                        "rho_hi": item["rho_hi"],
                        "cells": 1,
                    }
                )
            else:
                assert blocks[-1]["rho_hi"] == item["rho_lo"]
                blocks[-1]["rho_hi"] = item["rho_hi"]
                blocks[-1]["cells"] += 1
        coverage_blocks[scenario] = blocks

    # The logical structure of the cover is part of the certificate.  These
    # assertions prevent a collection of individually true but insufficient
    # cells from being followed by an unsupported hard-coded conclusion.
    assert [b["key"] for b in coverage_blocks[1]] == [
        ("derivative", 1),
        ("positive", None),
        ("derivative", -1),
    ]
    assert [b["key"] for b in coverage_blocks[2]] == [
        ("derivative", 1),
        ("positive", None),
        ("derivative", -1),
        ("negative", None),
    ]

    expected_roots = [
        ("rho_1", 2, 1),
        ("rho_2", 1, 1),
        ("rho_3", 1, -1),
        ("rho_4", 2, -1),
    ]
    assert [
        (item["name"], int(item["scenario"]), int(item["derivative_sign"]))
        for item in raw["root_boxes"]
    ] == expected_roots

    def report_lo(index: int) -> Decimal:
        return Decimal(raw["root_boxes"][index]["report_rho_lo"])

    def report_hi(index: int) -> Decimal:
        return Decimal(raw["root_boxes"][index]["report_rho_hi"])

    b1 = coverage_blocks[1]
    b2 = coverage_blocks[2]
    # rho_2 lies in the increasing Delta_1 block and rho_3 in the decreasing
    # block; the latter block also reaches beyond rho_4.
    assert report_hi(1) < Decimal(b1[0]["rho_hi"])
    assert Decimal(b1[2]["rho_lo"]) < report_lo(2)
    assert report_hi(2) < Decimal(b1[2]["rho_hi"])
    assert report_hi(3) < Decimal(b1[2]["rho_hi"])
    # rho_1 lies in the increasing Delta_2 block and rho_4 in the decreasing
    # block.  That decreasing block contains rho_2 and rho_3 as well.
    assert report_hi(0) < Decimal(b2[0]["rho_hi"])
    assert Decimal(b2[2]["rho_lo"]) < report_lo(3)
    assert report_hi(3) < Decimal(b2[2]["rho_hi"])
    assert report_hi(1) < Decimal(b2[0]["rho_hi"])
    assert Decimal(b2[2]["rho_lo"]) < report_lo(2)
    assert report_hi(2) < report_lo(3)
    # The direct-negative cells overlap the residual-transfer tail.
    assert Decimal(b2[3]["rho_hi"]) >= Decimal(
        raw["intercept_root"]["prethreshold_lower"]
    )

    # Ordering is rigorous because the contracted root intervals are disjoint.
    for left, right in zip(roots, roots[1:]):
        assert arb(left["rho_upper"]) < arb(right["rho_lower"])

    summary = {
        "status": "PASS",
        "certificate_sha256": certificate_sha256,
        "arithmetic": "python-flint Arb outward-rounded balls",
        "precision_bits": ctx.prec,
        "root_enclosures": roots,
        "validated_path_cells": len(cells),
        "coverage_blocks": {
            str(scenario): [
                {
                    "purpose": key[0],
                    "sign": key[1],
                    "rho_lo": block["rho_lo"],
                    "rho_hi": block["rho_hi"],
                    "cells": block["cells"],
                }
                for block in blocks
                for key in [block["key"]]
            ]
            for scenario, blocks in coverage_blocks.items()
        },
        "tail": tail,
        "conclusion": {
            "active_set_path": [
                "{2} on [0,rho_1)",
                "{1,2} on [rho_1,rho_2]",
                "{1} on (rho_2,rho_3)",
                "{1,2} on [rho_3,rho_4]",
                "{2} on (rho_4,infinity)",
            ],
            "A1": "[rho_1,rho_4]",
            "A2": "[0,rho_2] union [rho_3,infinity)",
        },
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
