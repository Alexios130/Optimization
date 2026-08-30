#!/usr/bin/env python3
"""Independent exact-rational interval checker for the n=3 path.

Every interval endpoint is a fractions.Fraction.  The only transcendental
operation needed by the KKT equations is exp; ``exp_bounds`` encloses it by
an exact Taylor sum plus an exact geometric tail bound.  The code therefore
does not trust floating-point arithmetic when it returns a certificate.

SciPy is used only to propose affine predictors.  A predictor is accepted
solely when the rational interval checks succeed.

The checker validates branch-gap signs away from four rational event
brackets, strict derivative signs on containing event windows, and the
zero-slope residual bridge to the analytic stabilization tail.  This checker
is independent of the primary Arb implementation and serves as a rigorous
cross-check; the Arb verifier remains the compact primary certificate.
"""

from __future__ import annotations

from fractions import Fraction as F

import numpy as np
from scipy.optimize import root


# Fixed decimal lattice used only to control fraction growth.  All rounding
# is outward, hence it cannot invalidate an enclosure.
Q = 10**18
EXP_TERMS = 32


def flq(x: F) -> F:
    return F((x.numerator * Q) // x.denominator, Q)


def ceq(x: F) -> F:
    return F(-((-x.numerator * Q) // x.denominator), Q)


def as_f(x) -> F:
    if isinstance(x, F):
        return x
    return F(str(x))


class I:
    """Closed interval on the exact fixed-point lattice (1/Q) Z."""

    __slots__ = ("_lo", "_hi")

    def __init__(self, lo, hi=None):
        lo = as_f(lo)
        hi = lo if hi is None else as_f(hi)
        if lo > hi:
            raise ValueError((lo, hi))
        self._lo = (lo.numerator * Q) // lo.denominator
        self._hi = -((-hi.numerator * Q) // hi.denominator)

    @classmethod
    def raw(cls, lo: int, hi: int):
        if lo > hi:
            raise ValueError((lo, hi))
        out = object.__new__(cls)
        out._lo = lo
        out._hi = hi
        return out

    @property
    def lo(self):
        return F(self._lo, Q)

    @property
    def hi(self):
        return F(self._hi, Q)

    def __add__(self, other):
        other = ii(other)
        return I.raw(self._lo + other._lo, self._hi + other._hi)

    __radd__ = __add__

    def __neg__(self):
        return I.raw(-self._hi, -self._lo)

    def __sub__(self, other):
        return self + (-ii(other))

    def __rsub__(self, other):
        return ii(other) - self

    def __mul__(self, other):
        other = ii(other)
        vals = (
            self._lo * other._lo,
            self._lo * other._hi,
            self._hi * other._lo,
            self._hi * other._hi,
        )
        return I.raw(min(vals) // Q, -((-max(vals)) // Q))

    __rmul__ = __mul__

    def reciprocal(self):
        if self._lo <= 0 <= self._hi:
            raise ZeroDivisionError(self)
        scale = Q * Q
        if self._lo > 0:
            return I.raw(scale // self._hi, -((-scale) // self._lo))
        # Both endpoints are negative.  1/x is decreasing.
        return I.raw(scale // self._hi, -((-scale) // self._lo))

    def __truediv__(self, other):
        return self * ii(other).reciprocal()

    def __rtruediv__(self, other):
        return ii(other) / self

    def square(self):
        if self._lo >= 0:
            return I.raw((self._lo * self._lo) // Q, -((-(self._hi * self._hi)) // Q))
        if self._hi <= 0:
            return I.raw((self._hi * self._hi) // Q, -((-(self._lo * self._lo)) // Q))
        high = max(self._lo * self._lo, self._hi * self._hi)
        return I.raw(0, -((-high) // Q))

    def pow(self, n: int):
        if n < 0:
            return self.pow(-n).reciprocal()
        out = I(1)
        base = self
        while n:
            if n & 1:
                out = out * base
            n >>= 1
            if n:
                base = base * base
        return out

    def abs_upper(self) -> F:
        return F(max(abs(self._lo), abs(self._hi)), Q)

    def __repr__(self):
        return f"[{float(self.lo):.6g},{float(self.hi):.6g}]"


def fraction_text(value: F) -> str:
    """Return an exact, unambiguous representation of a Fraction."""
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def ii(x) -> I:
    return x if isinstance(x, I) else I(x)


def exp_positive_bounds(x: F) -> I:
    """Enclose exp(x), x >= 0, using positive Taylor terms."""
    if x < 0:
        raise ValueError(x)
    term = I(1)
    total = I(1)
    for k in range(1, EXP_TERMS + 1):
        term = term * I(x) / k
        total = total + term
    first_omitted = term * I(x) / (EXP_TERMS + 1)
    ratio = x / F(EXP_TERMS + 2)
    if ratio >= 1:
        raise ValueError("increase EXP_TERMS")
    tail_hi = first_omitted.hi / (1 - ratio)
    return I(total.lo, ceq(total.hi + tail_hi))


def exp_point_bounds(x: F) -> I:
    if x >= 0:
        return exp_positive_bounds(x)
    return exp_positive_bounds(-x).reciprocal()


def exp_i(x: I) -> I:
    low = exp_point_bounds(x.lo)
    high = exp_point_bounds(x.hi)
    return I(low.lo, high.hi)


def logistic_i(q: I) -> I:
    e = exp_i(q)
    return e / (I(1) + e)


x = (F(7, 2), F(3), F(4))
y = (F(-1), F(1), F(-1))
p = (
    (F(12, 19), F(1, 19), F(6, 19)),
    (F(1, 8), F(1, 16), F(13, 16)),
)
d = tuple(a - b for a, b in zip(p[0], p[1]))
lam = F(1, 2)


def qzat(beta: I, b: I, rho: I):
    # This checker is for the smooth beta<0 branches, so |beta|=-beta.
    q = tuple(-yi * (xi * beta + b) - rho * beta for xi, yi in zip(x, y))
    z = tuple(I(-yi * xi) - rho for xi, yi in zip(x, y))
    t = tuple(F(-yi) for yi in y)
    a = tuple(logistic_i(qi) for qi in q)
    h = tuple(ai * (I(1) - ai) for ai in a)
    return q, z, t, a, h


def gradient_box(beta: I, b: I, rho: I, scenario: int):
    _, z, t, a, _ = qzat(beta, b, rho)
    e1 = I(lam) * beta
    e2 = I(lam) * b
    for wi, ai, zi, ti in zip(p[scenario], a, z, t):
        e1 = e1 + wi * ai * zi
        e2 = e2 + wi * ai * ti
    return e1, e2


def inflate_from_residual(beta: I, b: I, rho: I, scenario: int):
    e1, e2 = gradient_box(beta, b, rho, scenario)
    # ||residual||_2 <= ||residual||_1; strong convexity then gives
    # distance <= (|e1|+|e2|)/lambda.
    delta = (e1.abs_upper() + e2.abs_upper()) / lam
    pad = I(-delta, delta)
    return beta + pad, b + pad, delta, (e1, e2)


def gap_sign(beta: I, b: I, rho: I) -> int:
    """Sign of G without logarithms.

    304 G = log((1+e^q1)^154 / ((1+e^q2)^3(1+e^q3)^151)).
    """
    q, _, _, _, _ = qzat(beta, b, rho)
    c = tuple(I(1) + exp_i(qi) for qi in q)
    numerator = c[0].pow(154)
    denominator = c[1].pow(3) * c[2].pow(151)
    if numerator.lo > denominator.hi:
        return 1
    if numerator.hi < denominator.lo:
        return -1
    return 0


def derivative_gap_box(beta: I, b: I, rho: I, scenario: int) -> I:
    q, z, t, a, h = qzat(beta, b, rho)
    w = p[scenario]
    h11 = I(lam)
    h12 = I(0)
    h22 = I(lam)
    r1 = I(0)
    r2 = I(0)
    for wi, ai, hi, zi, ti in zip(w, a, h, z, t):
        h11 = h11 + wi * hi * zi.square()
        h12 = h12 + wi * hi * zi * ti
        h22 = h22 + wi * hi
        # partial_rho grad, holding (beta,b) fixed
        r1 = r1 + wi * (hi * (-beta) * zi - ai)
        r2 = r2 + wi * hi * ti * (-beta)
    det = h11 * h22 - h12.square()
    if det.lo <= 0:
        raise ArithmeticError(f"nonpositive determinant enclosure {det}")
    beta_p = -(h22 * r1 - h12 * r2) / det
    b_p = -((-h12) * r1 + h11 * r2) / det
    gb = I(0)
    gi = I(0)
    gr = I(0)
    for di, ai, zi, ti in zip(d, a, z, t):
        gb = gb + di * ai * zi
        gi = gi + di * ai * ti
        gr = gr + di * ai * (-beta)
    return gb * beta_p + gi * b_p + gr


# Floating-point proposer.  None of these values are trusted by the checker.
xf = np.array([3.5, 3.0, 4.0])
yf = np.array([-1.0, 1.0, -1.0])
pf = np.array([[12 / 19, 1 / 19, 6 / 19], [1 / 8, 1 / 16, 13 / 16]])
df = pf[0] - pf[1]


def propose(rho: float, scenario: int):
    def equations(theta):
        beta, b = theta
        q = -yf * (xf * beta + b) - rho * beta
        a = 1 / (1 + np.exp(-q))
        z = -yf * xf - rho
        return 0.5 * theta + pf[scenario] @ (a[:, None] * np.c_[z, -yf])

    sol = root(equations, np.array([-0.35, -0.35]), tol=1e-12)
    if np.linalg.norm(equations(sol.x)) > 1e-9 or sol.x[0] >= 0:
        raise ArithmeticError((rho, scenario, sol.x, equations(sol.x)))
    beta, b = sol.x
    q = -yf * (xf * beta + b) - rho * beta
    a = 1 / (1 + np.exp(-q))
    h = a * (1 - a)
    z = -yf * xf - rho
    w = pf[scenario]
    hes = np.array(
        [
            [0.5 + np.sum(w * h * z * z), np.sum(w * h * z * (-yf))],
            [np.sum(w * h * z * (-yf)), 0.5 + np.sum(w * h)],
        ]
    )
    er = np.array(
        [np.sum(w * (h * (-beta) * z - a)), np.sum(w * h * (-yf) * (-beta))]
    )
    slope = np.linalg.solve(hes, -er)
    return sol.x, slope


def affine_predictor(lo: F, hi: F, scenario: int):
    mid = (lo + hi) / 2
    theta, slope = propose(float(mid), scenario)
    r = I(lo, hi)
    dm = r - I(mid)
    beta = I(F(format(theta[0], ".17g"))) + I(F(format(slope[0], ".17g"))) * dm
    b = I(F(format(theta[1], ".17g"))) + I(F(format(slope[1], ".17g"))) * dm
    return r, beta, b


def certify_sign_cell(lo: F, hi: F, scenario: int, wanted: int):
    rho, bp, ip = affine_predictor(lo, hi, scenario)
    beta, b, delta, residual = inflate_from_residual(bp, ip, rho, scenario)
    sign = gap_sign(beta, b, rho) if beta.hi < 0 else 0
    ok = sign == wanted
    return ok, {
        "rho": (lo, hi),
        "beta": beta,
        "b": b,
        "delta": delta,
        "residual": residual,
        "sign": sign,
    }


def cover(lo: F, hi: F, scenario: int, wanted: int, depth=0):
    ok, record = certify_sign_cell(lo, hi, scenario, wanted)
    if ok:
        return [record]
    if depth >= 25:
        raise ArithmeticError(("coverage failed", float(lo), float(hi), scenario, record))
    mid = (lo + hi) / 2
    return cover(lo, mid, scenario, wanted, depth + 1) + cover(
        mid, hi, scenario, wanted, depth + 1
    )


def point_gap_sign(rho: F, scenario: int):
    r, bp, ip = affine_predictor(rho, rho, scenario)
    beta, b, delta, residual = inflate_from_residual(bp, ip, r, scenario)
    return gap_sign(beta, b, r), (beta, b, delta, residual)


def event_derivative(lo: F, hi: F, scenario: int):
    rho, bp, ip = affine_predictor(lo, hi, scenario)
    beta, b, delta, residual = inflate_from_residual(bp, ip, rho, scenario)
    if beta.hi >= 0:
        raise ArithmeticError("event box reaches beta=0")
    return derivative_gap_box(beta, b, rho, scenario), (beta, b, delta, residual)


def certify_zero_slope_tail(lo: F, hi: F, bbar: F):
    """Certify Delta_2<0 using a beta=0 subgradient residual.

    The selected slope subgradient is the left endpoint c-rho*H of the
    exact interval [c-rho*H,c+rho*H].  This remains a valid selection on
    both sides of the stabilization threshold.
    """
    rho = I(lo, hi)
    b = I(bbar)
    a = tuple(logistic_i(I(-yi) * b) for yi in y)
    c = I(0)
    hsum = I(0)
    db = I(lam) * b
    for wi, ai, yi, xi in zip(p[1], a, y, x):
        c = c - wi * ai * yi * xi
        hsum = hsum + wi * ai
        db = db - wi * ai * yi
    d_beta = c - rho * hsum
    delta = (d_beta.abs_upper() + db.abs_upper()) / lam

    # At beta=0 the exact identity Delta=(3/304)b avoids logarithms.
    gap_at_candidate = I(F(3, 304)) * b
    # K <= sum |d_i| (|x_i|+rho+1), since sqrt(u^2+1)<=u+1.
    k_upper = sum(abs(di) * (abs(xi) + rho.hi + 1) for di, xi in zip(d, x))
    transferred_upper = gap_at_candidate.hi + k_upper * delta
    if transferred_upper >= 0:
        raise ArithmeticError(
            ("zero-slope tail failed", rho, d_beta, db, delta, transferred_upper)
        )
    return {
        "rho": (lo, hi),
        "bbar": bbar,
        "slope_residual": d_beta,
        "intercept_residual": db,
        "delta": delta,
        "gap_upper": transferred_upper,
    }


def certify_stabilization_threshold():
    """Bracket b_infty and the induced exact rho_infty."""
    blo = F("-0.5888140329295523")
    bhi = F("-0.5888140329295522")

    def intercept_equation(b0):
        b = I(b0)
        return I(lam) * b + logistic_i(b) - I(F(1, 16))

    flo = intercept_equation(blo)
    fhi = intercept_equation(bhi)
    if flo.hi >= 0 or fhi.lo <= 0:
        raise ArithmeticError(("b-infinity bracket failed", flo, fhi))
    # f'(b)=1/2+sigma(b)(1-sigma(b))>0 proves uniqueness.
    b = I(blo, bhi)
    aa = logistic_i(b)
    rho = (I(62) * aa - I(3)) / (I(14) * aa + I(1))
    if bhi >= 0 or rho.lo <= 0:
        raise ArithmeticError(("wrong stabilization signs", b, rho))
    return (blo, bhi), rho, (flo, fhi)


def main():
    print(
        "display_note: decimal interval and width values below are "
        "human-readable summaries only; every pass/fail decision uses the "
        "exact Fraction endpoints stored internally"
    )
    # scenario uses zero-based indexing; wanted is the sign of L1-L2.
    # (narrow root lower, narrow root upper, monotonic-window lower,
    #  monotonic-window upper, scenario, derivative sign)
    r1 = (
        F("0.524565"), F("0.524566"), F("0.5245"), F("0.5246"), 1, +1
    )
    r2 = (
        F("1.036949"), F("1.036950"), F("1.0369"), F("1.0370"), 0, +1
    )
    r3 = (
        F("2.696120"), F("2.696121"), F("2.6960"), F("2.6962"), 0, -1
    )
    r4 = (
        F("2.832429"), F("2.832430"), F("2.8323"), F("2.8325"), 1, -1
    )
    events = (r1, r2, r3, r4)
    for k, (lo, hi, wlo, whi, scenario, derivative_sign) in enumerate(events, 1):
        slo, _ = point_gap_sign(lo, scenario)
        shi, _ = point_gap_sign(hi, scenario)
        deriv, _ = event_derivative(wlo, whi, scenario)
        derivative_ok = deriv.lo > 0 if derivative_sign > 0 else deriv.hi < 0
        if slo * shi != -1 or not derivative_ok:
            raise ArithmeticError(("event failed", k, slo, shi, deriv))
        print(
            f"event {k}: exact rho bracket "
            f"[{fraction_text(lo)},{fraction_text(hi)}], "
            f"endpoint signs {slo:+d}/{shi:+d}; G' decimal interval "
            f"summary {deriv} on exact window "
            f"[{fraction_text(wlo)},{fraction_text(whi)}]"
        )

    domains = (
        # scenario 2 gap g2
        (F(0), r1[2], 1, -1, "g2 before r1"),
        (r1[3], r4[2], 1, +1, "g2 between r1 and r4"),
        (r4[3], F("3.1890"), 1, -1, "g2 after r4"),
        # scenario 1 gap g1 (only needed through r4)
        (F(0), r2[2], 0, -1, "g1 before r2"),
        (r2[3], r3[2], 0, +1, "g1 between r2 and r3"),
        (r3[3], r4[3], 0, -1, "g1 after r3 through r4"),
    )
    total = 0
    for lo, hi, scenario, wanted, label in domains:
        cells = cover(lo, hi, scenario, wanted)
        total += len(cells)
        max_width = max(float(c["rho"][1] - c["rho"][0]) for c in cells)
        min_width = min(float(c["rho"][1] - c["rho"][0]) for c in cells)
        print(
            f"{label}: {len(cells)} cells, decimal width summaries "
            f"{min_width:.3g}..{max_width:.3g}"
        )
    b_box, rho_inf, f_box = certify_stabilization_threshold()
    tail = certify_zero_slope_tail(
        F("3.1890"), rho_inf.hi, F("-0.58881403292955225")
    )
    print(
        "decimal interval summaries: b_infinity in "
        f"[{float(b_box[0]):.16f},{float(b_box[1]):.16f}], "
        f"rho_infinity in [{float(rho_inf.lo):.16f},{float(rho_inf.hi):.16f}]"
    )
    print(
        "zero-slope residual tail certified on [3.1890,rho_infinity]: "
        f"transferred-gap upper decimal summary "
        f"{float(tail['gap_upper']):.6g}"
    )
    print("status: PASS")
    print("arithmetic: exact Fraction intervals with proved exponential remainders")
    print(f"lattice_denominator: {Q}")
    print(f"exp_terms: {EXP_TERMS}")
    print(f"validated_rational_sign_cells: {total}")


if __name__ == "__main__":
    main()
