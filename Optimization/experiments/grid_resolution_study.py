#!/usr/bin/env python3
"""Quantify when uniform active-set sampling can miss narrow tie intervals.

The inputs are the 90-digit Decimal outputs produced by
``check_linear_family_decimal.py``.  This program performs no optimization:
it reads the recorded event radii, identifies the open intervals whose
recorded probe state is ``12``, and carries out the grid calculations with
Python ``Decimal`` arithmetic.

The resulting numbers are numerical evidence attached to finite family
members.  They are not interval certificates and do not replace the analytic
all-m lower-bound proof in the article.
"""

from __future__ import annotations

import argparse
import json
import re
from decimal import Decimal, getcontext
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_INPUTS = [
    HERE / "explicit_family_m4_L16_decimal90.json",
    HERE / "explicit_family_m5_L16_decimal90.json",
    HERE / "explicit_family_m6_L16_decimal90.json",
]
DEFAULT_OUTPUT = HERE / "grid_resolution_study.json"


def decimal_text(value: Decimal, digits: int = 40) -> str:
    """Return a stable significant-digit representation for JSON output."""

    return format(value, f".{digits}g")


def first_open_grid_node(left: Decimal, spacing: Decimal) -> tuple[int, Decimal]:
    """Return the first uniform-grid node strictly to the right of ``left``."""

    index = int(left / spacing) + 1
    return index, Decimal(index) * spacing


def study_instance(
    source: Path,
    domain_right: Decimal,
    reference_subintervals: int,
) -> dict:
    raw = json.loads(source.read_text(encoding="utf-8"))
    if raw.get("status") != "90_DIGIT_NUMERICAL_EVIDENCE_NOT_INTERVAL_CERTIFICATION":
        raise ValueError(f"unexpected status in {source}")

    match = re.search(r"_m(\d+)_L(\d+)_decimal", source.name)
    if match is None:
        raise ValueError(f"cannot infer m and L from {source.name}")
    m = int(match.group(1))
    scale_parameter = int(match.group(2))

    events = [Decimal(item["r"]) for item in raw["events"]]
    probes = raw["interval_probes"]
    if len(probes) != len(events) + 1:
        raise ValueError(f"{source.name}: expected one probe per open path interval")
    if any(right <= left for left, right in zip(events, events[1:])):
        raise ValueError(f"{source.name}: event radii are not strictly increasing")
    if not events or events[-1] >= domain_right:
        raise ValueError(f"{source.name}: events must lie in (0,{domain_right})")

    tie_intervals = []
    tie_widths: list[Decimal] = []
    for j, (left, right) in enumerate(zip(events, events[1:])):
        if probes[j + 1]["state"] != "12":
            continue
        width = right - left
        tie_widths.append(width)
        tie_intervals.append(
            {
                "left_event_index": j + 1,
                "right_event_index": j + 2,
                "left_r": decimal_text(left, 75),
                "right_r": decimal_text(right, 75),
                "width": decimal_text(width, 45),
            }
        )
    if not tie_intervals:
        raise ValueError(f"{source.name}: no recorded tie interval")

    adjacent_separations = [right - left for left, right in zip(events, events[1:])]
    minimum_separation = min(adjacent_separations)
    minimum_tie_width = min(tie_widths)
    reported_minimum = Decimal(raw["minimum_event_separation_r"])
    report_tolerance = Decimal("5e-40") * max(Decimal(1), abs(minimum_separation))
    if abs(reported_minimum - minimum_separation) > report_tolerance:
        raise ValueError(f"{source.name}: recorded and recomputed minimum separations disagree")

    reference_spacing = domain_right / Decimal(reference_subintervals)
    hit_count = 0
    for interval in tie_intervals:
        left = Decimal(interval["left_r"])
        right = Decimal(interval["right_r"])
        index, node = first_open_grid_node(left, reference_spacing)
        hit = index <= reference_subintervals and node < right
        interval["reference_grid_first_node_index"] = index
        interval["reference_grid_first_node_r"] = decimal_text(node, 45)
        interval["reference_grid_hits_open_interval"] = hit
        hit_count += int(hit)

    # If h is strictly smaller than every open tie-interval width, each such
    # interval contains a grid node, independently of its alignment with the
    # grid.  The strict inequality explains the floor-plus-one formula.
    sufficient_subintervals = int(domain_right / minimum_tie_width) + 1
    sufficient_spacing = domain_right / Decimal(sufficient_subintervals)
    if not sufficient_spacing < minimum_tie_width:
        raise AssertionError("strict grid-spacing guarantee was not achieved")

    return {
        "source": source.name,
        "m": m,
        "n": 2 * m,
        "L": scale_parameter,
        "event_count": len(events),
        "tie_interval_count": len(tie_intervals),
        "interval_state_sequence": [item["state"] for item in probes],
        "event_radii_scaled": [decimal_text(value, 75) for value in events],
        "tie_intervals": tie_intervals,
        "minimum_adjacent_event_separation": decimal_text(minimum_separation, 45),
        "minimum_tie_width": decimal_text(minimum_tie_width, 45),
        "reported_minimum_event_separation": raw["minimum_event_separation_r"],
        "reference_grid": {
            "subintervals": reference_subintervals,
            "points": reference_subintervals + 1,
            "spacing": decimal_text(reference_spacing, 45),
            "spacing_over_minimum_tie_width": decimal_text(
                reference_spacing / minimum_tie_width, 20
            ),
            "tie_intervals_hit": hit_count,
            "tie_intervals_missed": len(tie_intervals) - hit_count,
            "alignment_independent_hit_guarantee": reference_spacing < minimum_tie_width,
        },
        "strict_uniform_spacing_guarantee": {
            "condition": "grid spacing h < minimum tie-interval width",
            "sufficient_subintervals": sufficient_subintervals,
            "sufficient_grid_points": sufficient_subintervals + 1,
            "resulting_spacing": decimal_text(sufficient_spacing, 45),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--domain-right", default="0.5")
    parser.add_argument("--reference-subintervals", type=int, default=200_000)
    args = parser.parse_args()

    getcontext().prec = 90
    if args.reference_subintervals <= 0:
        raise ValueError("--reference-subintervals must be positive")
    domain_right = Decimal(args.domain_right)
    if domain_right <= 0:
        raise ValueError("--domain-right must be positive")

    instances = [
        study_instance(path.resolve(), domain_right, args.reference_subintervals)
        for path in args.inputs
    ]
    payload = {
        "status": "DECIMAL_POSTPROCESSING_OF_NUMERICAL_EVIDENCE_NOT_CERTIFICATION",
        "claim_boundary": (
            "Event radii come from 90-digit Decimal Newton solves. The grid calculations "
            "are exact Decimal postprocessing of their recorded digits; neither operation "
            "is outward-rounded interval certification."
        ),
        "scaled_radius_domain": {"left": "0", "right": decimal_text(domain_right)},
        "reference_uniform_grid": {
            "subintervals": args.reference_subintervals,
            "points": args.reference_subintervals + 1,
        },
        "instances": instances,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    for item in instances:
        grid = item["reference_grid"]
        sufficient = item["strict_uniform_spacing_guarantee"]
        print(
            f"m={item['m']} n={item['n']} K={item['event_count']} "
            f"ties={item['tie_interval_count']} min_width={item['minimum_tie_width']} "
            f"reference_hits={grid['tie_intervals_hit']}/{item['tie_interval_count']} "
            f"guarantee_points={sufficient['sufficient_grid_points']}"
        )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
