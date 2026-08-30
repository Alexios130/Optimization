#!/usr/bin/env python3
"""Numerical search for many-switch p=1, S=2 robust-logistic paths.

This is a discovery program, not a proof.  For two scenarios it uses the
exact branch-gap classification from the thesis: if theta_s(rho) minimizes

    J_s = L_s + lambda/2 * (beta**2 + b**2),

then Delta_s = L_1(theta_s)-L_2(theta_s), Delta_1 <= Delta_2, and the
minimax active set is {1}, {2}, or {1,2} according to the signs of those two
gaps.  We solve the smooth single-scenario branches by damped Newton and use
Brent refinement for proposed gap roots.

All results printed by this program are floating-point candidates only.
They must be certified independently before entering a theorem.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar, root


def softplus(q: np.ndarray) -> np.ndarray:
    return np.logaddexp(0.0, q)


def sigmoid(q: np.ndarray) -> np.ndarray:
    out = np.empty_like(q, dtype=float)
    pos = q >= 0.0
    out[pos] = 1.0 / (1.0 + np.exp(-q[pos]))
    ez = np.exp(q[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


@dataclass
class Instance:
    x: list[float]
    y: list[int]
    pi: list[list[float]]
    lam: float

    def arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (
            np.asarray(self.x, dtype=float),
            np.asarray(self.y, dtype=float),
            np.asarray(self.pi, dtype=float),
        )


@dataclass
class PathSummary:
    K: int
    roots: list[float]
    root_branches: list[int]
    interval_states: list[str]
    rho_max: float
    min_event_separation: float
    min_probe_margin: float


class TwoScenarioPath:
    def __init__(self, instance: Instance):
        self.instance = instance
        self.x, self.y, self.pi = instance.arrays()
        self.lam = float(instance.lam)
        if self.pi.shape != (2, len(self.x)):
            raise ValueError("pi must have shape (2,n)")
        if not np.allclose(self.pi.sum(axis=1), 1.0, atol=1e-10):
            raise ValueError("scenario rows must sum to one")
        self.b0 = np.array([self._intercept_root(s) for s in range(2)])
        self.sign = np.empty(2, dtype=int)
        self.threshold = np.empty(2, dtype=float)
        for s in range(2):
            a = sigmoid(-self.y * self.b0[s])
            kappa = float(np.dot(self.pi[s], a))
            g = float(np.dot(self.pi[s], a * (-self.y * self.x)))
            self.sign[s] = -1 if g > 0.0 else 1
            self.threshold[s] = abs(g) / kappa if kappa > 0.0 else 0.0

    def _intercept_root(self, scenario: int) -> float:
        p = self.pi[scenario]

        def fun(b: float) -> float:
            return self.lam * b + float(np.dot(p, sigmoid(-self.y * b) * (-self.y)))

        # lambda*b dominates outside this deterministic bracket.
        B = 2.0 / self.lam + 2.0
        return float(brentq(fun, -B, B, xtol=2e-14, rtol=2e-14))

    def _objective_half(self, z: np.ndarray, scenario: int, rho: float, h: int) -> float:
        u, b = z
        q = -self.y * (h * u * self.x + b) + rho * u
        return float(np.dot(self.pi[scenario], softplus(q)) + 0.5 * self.lam * (u * u + b * b))

    def _grad_hess_half(
        self, z: np.ndarray, scenario: int, rho: float, h: int
    ) -> tuple[np.ndarray, np.ndarray]:
        u, b = z
        c = rho - self.y * h * self.x
        e = -self.y
        q = c * u + e * b
        a = sigmoid(q)
        w = self.pi[scenario] * a * (1.0 - a)
        grad = np.array(
            [
                self.lam * u + np.dot(self.pi[scenario] * a, c),
                self.lam * b + np.dot(self.pi[scenario] * a, e),
            ]
        )
        H = np.array(
            [
                [self.lam + np.dot(w, c * c), np.dot(w, c * e)],
                [np.dot(w, c * e), self.lam + np.dot(w, e * e)],
            ]
        )
        return grad, H

    def solve_branch(
        self, scenario: int, rho: float, guess: np.ndarray | None = None
    ) -> np.ndarray:
        """Return (beta,b) for the unique minimizer of J_s at rho."""
        if rho >= self.threshold[scenario] - 2e-13:
            return np.array([0.0, self.b0[scenario]])
        h = int(self.sign[scenario])
        if guess is None or h * guess[0] <= 0.0:
            frac = max(1e-3, 1.0 - rho / max(self.threshold[scenario], 1e-15))
            z0 = np.array([0.2 * frac, self.b0[scenario]])
        else:
            z0 = np.array([h * guess[0], guess[1]])

        ans = root(
            lambda z: self._grad_hess_half(z, scenario, rho, h)[0],
            z0,
            jac=lambda z: self._grad_hess_half(z, scenario, rho, h)[1],
            method="hybr",
            options={"xtol": 1e-11},
        )
        if ans.success and ans.x[0] > -1e-9:
            u = max(0.0, float(ans.x[0]))
            theta = np.array([h * u, float(ans.x[1])])
            grad, _ = self._grad_hess_half(np.array([u, theta[1]]), scenario, rho, h)
            if np.linalg.norm(grad, ord=np.inf) < 5e-8:
                return theta

        # Slower safety fallback.  Strong convexity makes the bounded-half
        # minimization unambiguous.
        opt = minimize(
            lambda z: self._objective_half(z, scenario, rho, h),
            x0=np.array([max(0.0, z0[0]), z0[1]]),
            jac=lambda z: self._grad_hess_half(z, scenario, rho, h)[0],
            bounds=[(0.0, None), (None, None)],
            method="L-BFGS-B",
            options={"ftol": 1e-14, "gtol": 1e-10, "maxiter": 300},
        )
        if not opt.success:
            raise RuntimeError(f"branch solve failed: s={scenario+1}, rho={rho}: {opt.message}")
        return np.array([h * opt.x[0], opt.x[1]])

    def losses(self, theta: np.ndarray, rho: float) -> np.ndarray:
        beta, b = theta
        q = -self.y * (beta * self.x + b) + rho * abs(beta)
        r = softplus(q)
        return self.pi @ r

    def gaps(self, rho: float, guesses: list[np.ndarray] | None = None) -> tuple[np.ndarray, list[np.ndarray]]:
        thetas: list[np.ndarray] = []
        vals = np.empty(2)
        for s in range(2):
            guess = None if guesses is None else guesses[s]
            theta = self.solve_branch(s, rho, guess)
            thetas.append(theta)
            loss = self.losses(theta, rho)
            vals[s] = loss[0] - loss[1]
        return vals, thetas

    def gap_derivative(self, scenario: int, rho: float, theta: np.ndarray | None = None) -> float:
        """Implicit derivative of Delta_s on a nonzero smooth branch."""
        if theta is None:
            theta = self.solve_branch(scenario, rho)
        beta, b = theta
        if abs(beta) < 1e-10:
            raise ValueError("gap derivative requested on the zero-slope tail")
        h = 1 if beta > 0 else -1
        u = abs(beta)
        c = rho - self.y * h * self.x
        e = -self.y
        q = c * u + e * b
        a = sigmoid(q)
        ap = a * (1.0 - a)
        _, H = self._grad_hess_half(np.array([u, b]), scenario, rho, h)
        grad_rho = np.array(
            [
                np.dot(self.pi[scenario], a + ap * u * c),
                np.dot(self.pi[scenario], ap * u * e),
            ]
        )
        zprime = -np.linalg.solve(H, grad_rho)
        dpi = self.pi[0] - self.pi[1]
        gap_z = np.array([np.dot(dpi * a, c), np.dot(dpi * a, e)])
        gap_rho = float(np.dot(dpi * a, np.full_like(a, u)))
        return float(gap_rho + np.dot(gap_z, zprime))

    @staticmethod
    def state(gaps: Iterable[float], tol: float = 2e-8) -> str:
        d1, d2 = gaps
        if d1 > tol:
            return "1"
        if d2 < -tol:
            return "2"
        return "12"

    def scan(self, grid_size: int = 301, refine: bool = True) -> PathSummary:
        rho_max = float(max(self.threshold))
        if rho_max <= 1e-12:
            gaps, _ = self.gaps(0.0)
            return PathSummary(0, [], [], [self.state(gaps)], rho_max, math.inf, float(np.min(np.abs(gaps))))

        # A mixed linear/Chebyshev grid resolves both early and late events.
        t = np.linspace(0.0, 1.0, grid_size)
        cheb = 0.5 * (1.0 - np.cos(np.pi * t))
        rho_grid = np.unique(np.r_[rho_max * t, rho_max * cheb, self.threshold])
        rho_grid.sort()
        gap_grid = np.empty((len(rho_grid), 2))
        guesses: list[np.ndarray] | None = None
        for k, rho in enumerate(rho_grid):
            gap_grid[k], guesses = self.gaps(float(rho), guesses)

        brackets: list[tuple[float, float, int]] = []
        for s in range(2):
            for k in range(len(rho_grid) - 1):
                a, b = gap_grid[k, s], gap_grid[k + 1, s]
                if a == 0.0:
                    continue
                if a * b < 0.0:
                    brackets.append((float(rho_grid[k]), float(rho_grid[k + 1]), s))

        roots: list[tuple[float, int]] = []
        for lo, hi, s in brackets:
            if refine:
                def f(r: float) -> float:
                    return float(self.gaps(r)[0][s])

                try:
                    rr = float(brentq(f, lo, hi, xtol=2e-12, rtol=2e-12, maxiter=100))
                except (ValueError, RuntimeError):
                    rr = 0.5 * (lo + hi)
            else:
                rr = 0.5 * (lo + hi)
            if not roots or abs(rr - roots[-1][0]) > 5e-8 * max(1.0, rho_max) or s != roots[-1][1]:
                roots.append((rr, s))
        roots.sort()

        # Collapse nearly coincident roots before classifying open intervals.
        unique_boundaries: list[float] = []
        for rr, _ in roots:
            if not unique_boundaries or rr - unique_boundaries[-1] > 5e-8 * max(1.0, rho_max):
                unique_boundaries.append(rr)

        edges = [0.0] + unique_boundaries + [rho_max]
        states: list[str] = []
        margins: list[float] = []
        for k in range(len(edges) - 1):
            lo, hi = edges[k], edges[k + 1]
            if hi - lo <= 1e-12 * max(1.0, rho_max):
                continue
            mid = 0.5 * (lo + hi)
            gg, _ = self.gaps(mid)
            states.append(self.state(gg, tol=1e-9))
            margins.append(float(max(abs(gg[0]), abs(gg[1]))))
        if not states:
            gg, _ = self.gaps(0.0)
            states = [self.state(gg)]
            margins = [float(max(abs(gg[0]), abs(gg[1])))]

        # Keep only boundaries across which the active set really changes.
        effective_roots: list[float] = []
        effective_branches: list[int] = []
        effective_states = [states[0]]
        for j, boundary in enumerate(unique_boundaries):
            if j + 1 >= len(states):
                break
            if states[j + 1] != effective_states[-1]:
                effective_roots.append(boundary)
                # record the branch whose refined root is nearest
                branch = min(roots, key=lambda z: abs(z[0] - boundary))[1]
                effective_branches.append(branch + 1)
                effective_states.append(states[j + 1])

        sep = min(np.diff(effective_roots)) if len(effective_roots) >= 2 else math.inf
        return PathSummary(
            K=len(effective_roots),
            roots=effective_roots,
            root_branches=effective_branches,
            interval_states=effective_states,
            rho_max=rho_max,
            min_event_separation=float(sep),
            min_probe_margin=float(min(margins)),
        )


def random_instance(rng: np.random.Generator, n: int, alpha: float = 0.35) -> Instance:
    # Random scales are deliberately heterogeneous: separated sigmoid
    # transition scales are the most plausible source of repeated changes.
    signs = rng.choice([-1.0, 1.0], size=n)
    magnitudes = np.exp(rng.uniform(math.log(0.35), math.log(12.0), size=n))
    x = signs * magnitudes
    y = rng.choice([-1, 1], size=n)
    if np.all(y == y[0]):
        y[rng.integers(n)] *= -1
    pi = rng.dirichlet(np.full(n, alpha), size=2)
    lam = float(np.exp(rng.uniform(math.log(0.03), math.log(2.0))))
    return Instance(x=x.tolist(), y=y.astype(int).tolist(), pi=pi.tolist(), lam=lam)


def canonical_n3() -> Instance:
    return Instance(
        x=[3.5, 3.0, 4.0],
        y=[-1, 1, -1],
        pi=[[12 / 19, 1 / 19, 6 / 19], [1 / 8, 1 / 16, 13 / 16]],
        lam=0.5,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--grid", type=int, default=181)
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument("--alpha", type=float, default=0.35)
    parser.add_argument("--keep", type=int, default=12)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-n3", action="store_true")
    args = parser.parse_args()

    if args.check_n3:
        path = TwoScenarioPath(canonical_n3())
        summary = path.scan(grid_size=max(args.grid, 301))
        print(json.dumps({"instance": asdict(canonical_n3()), "summary": asdict(summary)}, indent=2))
        return

    rng = np.random.default_rng(args.seed)
    kept: list[dict] = []
    failures = 0
    best_k = -1
    for trial in range(args.trials):
        inst = random_instance(rng, args.n, args.alpha)
        try:
            summary = TwoScenarioPath(inst).scan(grid_size=args.grid, refine=True)
        except (RuntimeError, ValueError, FloatingPointError):
            failures += 1
            continue
        record = {"trial": trial, "instance": asdict(inst), "summary": asdict(summary)}
        kept.append(record)
        kept.sort(
            key=lambda r: (
                r["summary"]["K"],
                r["summary"]["min_event_separation"],
                r["summary"]["min_probe_margin"],
            ),
            reverse=True,
        )
        kept = kept[: args.keep]
        if summary.K > best_k:
            best_k = summary.K
            print(
                f"trial={trial} new_best_K={best_k} states={summary.interval_states} "
                f"roots={[round(r, 8) for r in summary.roots]}",
                flush=True,
            )

    payload = {
        "status": "NUMERICAL_CANDIDATES_NOT_PROOFS",
        "seed": args.seed,
        "n": args.n,
        "trials": args.trials,
        "grid": args.grid,
        "alpha": args.alpha,
        "failures": failures,
        "top": kept,
    }
    text = json.dumps(payload, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print(f"wrote {args.output}")
    else:
        print(text)


if __name__ == "__main__":
    main()
