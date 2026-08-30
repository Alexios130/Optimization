#!/usr/bin/env python3
"""Render the certified n=3 active-set path and its numerical branch gaps.

The exact rational model and the event enclosures are read from the Arb
certificate files.  Smooth curves between the event markers are numerical
samples produced by the existing two-scenario solver in ``search_switches``.
Consequently, the figure distinguishes certified structure from numerical
rendering rather than presenting sampled curves as a proof.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from search_switches import Instance, TwoScenarioPath


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_CERTIFICATE = ROOT / "certificate" / "n3_certificate.json"
DEFAULT_SUMMARY = ROOT / "certificate" / "n3_certificate_summary.json"
DEFAULT_DATA = HERE / "n3_path_figure_data.json"
DEFAULT_PDF = HERE / "generated" / "n3_certified_path.pdf"


def rational(pair: list[int]) -> Fraction:
    return Fraction(int(pair[0]), int(pair[1]))


def load_exact_instance(model: dict) -> Instance:
    return Instance(
        x=[float(rational(value)) for value in model["x"]],
        y=[int(value) for value in model["y"]],
        pi=[
            [float(rational(value)) for value in scenario]
            for scenario in model["pi"]
        ],
        lam=float(rational(model["lambda"])),
    )


def certificate_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def event_records(summary: dict) -> list[dict]:
    records = []
    for item in summary["root_enclosures"]:
        lo = Fraction(item["report_rho_lo"])
        hi = Fraction(item["report_rho_hi"])
        records.append(
            {
                "name": item["name"],
                "scenario": int(item["scenario"]),
                "rho_lo": item["report_rho_lo"],
                "rho_hi": item["report_rho_hi"],
                "rho_mid": float((lo + hi) / 2),
                "derivative_lo": item["report_derivative_lo"],
                "derivative_hi": item["report_derivative_hi"],
            }
        )
    if any(b["rho_mid"] <= a["rho_mid"] for a, b in zip(records, records[1:])):
        raise ValueError("certificate event enclosures are not strictly ordered")
    return records


def certified_states(summary: dict, event_count: int) -> tuple[list[str], list[str]]:
    """Read, validate, and format the itinerary recorded by the Arb verifier."""
    entries = summary.get("conclusion", {}).get("active_set_path")
    if not isinstance(entries, list) or len(entries) != event_count + 1:
        raise ValueError("the passing summary has an inconsistent active-set itinerary")

    labels: list[str] = []
    solver_codes: list[str] = []
    for entry in entries:
        match = re.fullmatch(r"\{([12](?:,[12])?)\} on .+", str(entry))
        if match is None:
            raise ValueError(f"cannot parse active-set entry {entry!r}")
        members = match.group(1)
        labels.append(rf"$\{{{members}\}}$")
        solver_codes.append(members.replace(",", ""))
    return labels, solver_codes


def numerical_samples(path: TwoScenarioPath, rho_max: float, count: int) -> tuple[np.ndarray, np.ndarray]:
    rho = np.linspace(0.0, rho_max, count)
    gaps = np.empty((count, 2), dtype=float)
    for j, value in enumerate(rho):
        gaps[j], _ = path.gaps(float(value))
    if not np.all(np.isfinite(gaps)):
        raise FloatingPointError("non-finite branch gap in figure data")
    return rho, gaps


def write_data(
    output: Path,
    model: dict,
    certificate_hash: str,
    summary: dict,
    events: list[dict],
    rho: np.ndarray,
    gaps: np.ndarray,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": "NUMERICAL_CURVES_WITH_CERTIFIED_EVENT_MARKERS",
        "claim_boundary": (
            "The event intervals and active-set path come from the passing Arb summary. "
            "The plotted branch-gap curves are double-precision samples from "
            "search_switches.TwoScenarioPath and are included only for visualization."
        ),
        "certificate_sha256": certificate_hash,
        "certificate_status": summary["status"],
        "certificate_arithmetic": summary["arithmetic"],
        "certificate_precision_bits": summary["precision_bits"],
        "exact_model": model,
        "branch_gap_definition": "Delta_s(rho)=L_1(theta_s(rho);rho)-L_2(theta_s(rho);rho)",
        "events": events,
        "certified_active_set_path": summary["conclusion"]["active_set_path"],
        "sample_count": len(rho),
        "rho": [format(value, ".17g") for value in rho],
        "delta_1": [format(value, ".17g") for value in gaps[:, 0]],
        "delta_2": [format(value, ".17g") for value in gaps[:, 1]],
    }
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def render_figure(
    pdf_path: Path,
    rho: np.ndarray,
    gaps: np.ndarray,
    events: list[dict],
    state_labels: list[str],
    rho_max: float,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(7.05, 4.15), layout="constrained")
    grid = fig.add_gridspec(2, 1, height_ratios=(4.1, 1.0), hspace=0.04)
    ax = fig.add_subplot(grid[0])
    ribbon = fig.add_subplot(grid[1], sharex=ax)

    blue = "#0072B2"
    vermillion = "#D55E00"
    event_color = "#4D4D4D"
    ax.plot(rho, 1_000.0 * gaps[:, 0], color=blue, lw=1.8, label=r"$\Delta_1(\rho)$")
    ax.plot(rho, 1_000.0 * gaps[:, 1], color=vermillion, lw=1.8, label=r"$\Delta_2(\rho)$")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.grid(axis="y", color="#D9D9D9", lw=0.5)
    ax.set_ylabel(r"branch gap ($10^{-3}$)")
    ax.legend(loc="upper left", frameon=False, ncol=2, handlelength=2.4)
    ax.tick_params(axis="x", labelbottom=False)
    ax.text(
        0.988,
        0.96,
        "(a)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontweight="bold",
    )

    event_radii = [item["rho_mid"] for item in events]
    for j, value in enumerate(event_radii, start=1):
        ax.axvline(value, color=event_color, lw=0.8, ls=(0, (3, 2)), zorder=0)
        ribbon.axvline(value, color=event_color, lw=0.8, ls=(0, (3, 2)), zorder=3)
        ax.annotate(
            rf"$\rho_{j}$",
            xy=(value, 1.0),
            xycoords=("data", "axes fraction"),
            xytext=(0, -3 if j % 2 else -15),
            textcoords="offset points",
            ha="center",
            va="top",
            color=event_color,
        )

    edges = [0.0, *event_radii, rho_max]
    colors = ["#B9DDF2", "#E2E2E2", "#F4C6AE", "#E2E2E2", "#B9DDF2"]
    if len(state_labels) != len(colors):
        raise ValueError("the figure layout expects the certified five-cell itinerary")
    for j, (left, right, state, color) in enumerate(
        zip(edges, edges[1:], state_labels, colors)
    ):
        hatch = "////" if j in (1, 3) else None
        ribbon.add_patch(
            Rectangle(
                (left, 0.0),
                right - left,
                1.0,
                facecolor=color,
                edgecolor="#777777",
                linewidth=0.55,
                hatch=hatch,
            )
        )
        ribbon.text(
            0.5 * (left + right),
            0.5,
            state,
            ha="center",
            va="center",
            fontsize=6.2 if j == 3 else 8.5,
        )
    ribbon.set_ylim(0.0, 1.0)
    ribbon.set_xlim(0.0, rho_max)
    ribbon.set_yticks([])
    ribbon.set_ylabel(r"$I(\rho)$", rotation=0, labelpad=19)
    ribbon.set_xlabel(r"uncertainty radius $\rho$")
    ribbon.text(0.012, 0.88, "(b)", transform=ribbon.transAxes, va="top", fontweight="bold")
    for spine in ("left", "right", "top"):
        ribbon.spines[spine].set_visible(False)

    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fixed_time = datetime(2026, 8, 30, tzinfo=timezone.utc)
    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        metadata={
            "Title": "Certified n=3 active-set path",
            "Author": "Alexis Seferlis",
            "Creator": "experiments/make_n3_path_figure.py",
            "CreationDate": fixed_time,
            "ModDate": fixed_time,
        },
    )
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--samples", type=int, default=801)
    args = parser.parse_args()
    if args.samples < 101:
        raise ValueError("--samples must be at least 101")

    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    digest = certificate_sha256(args.certificate)
    if summary.get("status") != "PASS":
        raise ValueError("the n=3 certificate summary does not pass")
    if digest != summary.get("certificate_sha256"):
        raise ValueError("certificate SHA-256 does not match the passing summary")

    instance = load_exact_instance(certificate["model"])
    path = TwoScenarioPath(instance)
    events = event_records(summary)
    state_labels, certified_solver_codes = certified_states(summary, len(events))
    rho_max = float(Fraction(summary["tail"]["report_rho_infty_hi"]))
    rho, gaps = numerical_samples(path, rho_max, args.samples)

    # This numerical assertion checks that the rendering agrees with the
    # certified state labels away from the tiny event enclosures.
    edges = [0.0, *[item["rho_mid"] for item in events], rho_max]
    observed_states = []
    for left, right in zip(edges, edges[1:]):
        values, _ = path.gaps(0.5 * (left + right))
        observed_states.append(path.state(values, tol=1e-9))
    if observed_states != certified_solver_codes:
        raise ValueError(
            "numerical rendering disagrees with the certified itinerary: "
            f"observed {observed_states}, certified {certified_solver_codes}"
        )

    write_data(args.data, certificate["model"], digest, summary, events, rho, gaps)
    render_figure(args.pdf, rho, gaps, events, state_labels, rho_max)
    print(f"certificate_sha256={digest}")
    print(f"event_midpoints={[format(item['rho_mid'], '.16g') for item in events]}")
    print(f"state_sequence={observed_states}")
    print(f"wrote {args.data}")
    print(f"wrote {args.pdf}")


if __name__ == "__main__":
    main()
