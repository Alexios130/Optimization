#!/usr/bin/env python3
"""High-precision exploratory calculations for the rational n=3 instance.

This is not, by itself, an interval certificate.  It computes accurate
centres for the interval-Newton boxes that a proof-producing checker must
validate.
"""

from decimal import Decimal as D, getcontext


getcontext().prec = 90

ZERO = D(0)
ONE = D(1)
HALF = D(1) / D(2)

x = (D(7) / D(2), D(3), D(4))
y = (D(-1), D(1), D(-1))
p1 = (D(12) / D(19), D(1) / D(19), D(6) / D(19))
p2 = (D(1) / D(8), D(1) / D(16), D(13) / D(16))
d = tuple(a - b for a, b in zip(p1, p2))
lam = HALF


def exp(z):
    return z.exp()


def logistic(z):
    if z >= ZERO:
        e = exp(-z)
        return ONE / (ONE + e)
    e = exp(z)
    return e / (ONE + e)


def softplus(z):
    if z >= ZERO:
        return z + (ONE + exp(-z)).ln()
    return (ONE + exp(z)).ln()


def qza(beta, b, rho):
    q = tuple(-yi * (xi * beta + b) - rho * beta for xi, yi in zip(x, y))
    z = tuple(-yi * xi - rho for xi, yi in zip(x, y))
    t = tuple(-yi for yi in y)
    a = tuple(logistic(qi) for qi in q)
    h = tuple(ai * (ONE - ai) for ai in a)
    return q, z, t, a, h


def gap(beta, b, rho):
    q, _, _, _, _ = qza(beta, b, rho)
    return sum(di * softplus(qi) for di, qi in zip(d, q))


def event_system(v, alpha):
    beta, b, rho = v
    q, z, t, a, h = qza(beta, b, rho)
    w = tuple(p2i + alpha * di for p2i, di in zip(p2, d))
    e1 = lam * beta + sum(wi * ai * zi for wi, ai, zi in zip(w, a, z))
    e2 = lam * b + sum(wi * ai * ti for wi, ai, ti in zip(w, a, t))
    e3 = sum(di * softplus(qi) for di, qi in zip(d, q))
    qr = tuple(-beta for _ in q)
    j11 = lam + sum(wi * hi * zi * zi for wi, hi, zi in zip(w, h, z))
    j12 = sum(wi * hi * zi * ti for wi, hi, zi, ti in zip(w, h, z, t))
    j13 = sum(
        wi * (hi * qri * zi - ai)
        for wi, hi, qri, zi, ai in zip(w, h, qr, z, a)
    )
    j21 = j12
    j22 = lam + sum(wi * hi for wi, hi in zip(w, h))
    j23 = sum(wi * hi * ti * qri for wi, hi, ti, qri in zip(w, h, t, qr))
    j31 = sum(di * ai * zi for di, ai, zi in zip(d, a, z))
    j32 = sum(di * ai * ti for di, ai, ti in zip(d, a, t))
    j33 = sum(di * ai * qri for di, ai, qri in zip(d, a, qr))
    return (e1, e2, e3), (
        (j11, j12, j13),
        (j21, j22, j23),
        (j31, j32, j33),
    )


def tie_system(v, rho):
    beta, b, alpha = v
    q, z, t, a, h = qza(beta, b, rho)
    w = tuple(p2i + alpha * di for p2i, di in zip(p2, d))
    e1 = lam * beta + sum(wi * ai * zi for wi, ai, zi in zip(w, a, z))
    e2 = lam * b + sum(wi * ai * ti for wi, ai, ti in zip(w, a, t))
    e3 = sum(di * softplus(qi) for di, qi in zip(d, q))
    j11 = lam + sum(wi * hi * zi * zi for wi, hi, zi in zip(w, h, z))
    j12 = sum(wi * hi * zi * ti for wi, hi, zi, ti in zip(w, h, z, t))
    j21 = j12
    j22 = lam + sum(wi * hi for wi, hi in zip(w, h))
    # The alpha column is the gradient of L1-L2.
    j13 = sum(di * ai * zi for di, ai, zi in zip(d, a, z))
    j23 = sum(di * ai * ti for di, ai, ti in zip(d, a, t))
    return (e1, e2, e3), (
        (j11, j12, j13),
        (j21, j22, j23),
        (j13, j23, ZERO),
    )


def solve_linear(a, b):
    a = [list(row) for row in a]
    b = list(b)
    n = len(b)
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(a[i][k]))
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            b[k], b[pivot] = b[pivot], b[k]
        piv = a[k][k]
        for i in range(k + 1, n):
            f = a[i][k] / piv
            if f == ZERO:
                continue
            for j in range(k, n):
                a[i][j] -= f * a[k][j]
            b[i] -= f * b[k]
    out = [ZERO] * n
    for i in range(n - 1, -1, -1):
        out[i] = (b[i] - sum(a[i][j] * out[j] for j in range(i + 1, n))) / a[i][i]
    return out


def newton(fun, initial, max_iter=30):
    v = tuple(D(str(z)) for z in initial)
    for _ in range(max_iter):
        f, j = fun(v)
        delta = solve_linear(j, tuple(-u for u in f))
        v = tuple(vi + di for vi, di in zip(v, delta))
        if max(abs(di) for di in delta) < D("1e-80"):
            break
    f, _ = fun(v)
    return v, f


def branch_derivative(beta, b, rho, alpha):
    """Derivative along a unique-scenario branch, plus dG/drho."""
    q, z, t, a, h = qza(beta, b, rho)
    w = tuple(p2i + alpha * di for p2i, di in zip(p2, d))
    qr = tuple(-beta for _ in q)
    h11 = lam + sum(wi * hi * zi * zi for wi, hi, zi in zip(w, h, z))
    h12 = sum(wi * hi * zi * ti for wi, hi, zi, ti in zip(w, h, z, t))
    h22 = lam + sum(wi * hi for wi, hi in zip(w, h))
    er1 = sum(
        wi * (hi * qri * zi - ai)
        for wi, hi, qri, zi, ai in zip(w, h, qr, z, a)
    )
    er2 = sum(wi * hi * ti * qri for wi, hi, ti, qri in zip(w, h, t, qr))
    beta_r, b_r = solve_linear(((h11, h12), (h12, h22)), (-er1, -er2))
    gb = sum(di * ai * zi for di, ai, zi in zip(d, a, z))
    gi = sum(di * ai * ti for di, ai, ti in zip(d, a, t))
    gr = sum(di * ai * qri for di, ai, qri in zip(d, a, qr))
    dg = gb * beta_r + gi * b_r + gr
    return beta_r, b_r, dg


def tie_derivative(beta, b, alpha, rho):
    """Derivative (beta', b', alpha') along the equality/KKT branch."""
    _, jac = tie_system((beta, b, alpha), rho)
    q, z, t, a, h = qza(beta, b, rho)
    w = tuple(p2i + alpha * di for p2i, di in zip(p2, d))
    qr = tuple(-beta for _ in q)
    fr1 = sum(
        wi * (hi * qri * zi - ai)
        for wi, hi, qri, zi, ai in zip(w, h, qr, z, a)
    )
    fr2 = sum(wi * hi * ti * qri for wi, hi, ti, qri in zip(w, h, t, qr))
    fr3 = sum(di * ai * qri for di, ai, qri in zip(d, a, qr))
    return tuple(solve_linear(jac, (-fr1, -fr2, -fr3)))


def rho_infty():
    def f(v):
        (b,) = v
        a = logistic(b)
        # Scenario 2: negative-label mass 15/16, positive-label mass 1/16.
        val = lam * b + a - D(1) / D(16)
        der = lam + a * (ONE - a)
        return (val,), ((der,),)

    (b,), residual = newton(f, ("-0.5888",))
    ai = tuple(logistic(-yi * b) for yi in y)
    c = -sum(wi * aa * yi * xi for wi, aa, yi, xi in zip(p2, ai, y, x))
    h = sum(wi * aa for wi, aa in zip(p2, ai))
    return b, c / h, residual[0]


def main():
    starts = (
        (ZERO, ("-.4947", "-.1498", ".5246")),
        (ONE, ("-.4931", "-.2412", "1.0369")),
        (ONE, ("-.1380", "-.5515", "2.6961")),
        (ZERO, ("-.1260", "-.5354", "2.8324")),
    )
    events = []
    for alpha, start in starts:
        v, f = newton(lambda z, a=alpha: event_system(z, a), start)
        deriv = branch_derivative(*v, alpha)
        events.append((alpha, v, f, deriv))
    for k, (alpha, v, f, deriv) in enumerate(events, 1):
        print(f"event {k}; alpha={alpha}")
        print("  beta =", v[0])
        print("  b    =", v[1])
        print("  rho  =", v[2])
        print("  residual =", tuple(str(u) for u in f))
        print("  beta', b', G' =", tuple(str(u) for u in deriv))
        td = tie_derivative(v[0], v[1], alpha, v[2])
        print("  tie beta', b', alpha' =", tuple(str(u) for u in td))
    b, rho, residual = rho_infty()
    print("rho infinity")
    print("  b    =", b)
    print("  rho  =", rho)
    print("  residual =", residual)

    # Useful interior tie points, for validated-continuation seed boxes.
    tie_starts = (
        ("0.75", ("-.495", "-.19", ".4")),
        ("0.90", ("-.495", "-.215", ".7")),
        ("2.74", ("-.145", "-.55", ".7")),
        ("2.78", ("-.135", "-.54", ".4")),
    )
    for rr, start in tie_starts:
        rho = D(rr)
        v, f = newton(lambda z, r=rho: tie_system(z, r), start)
        print(f"tie rho={rr}: beta={v[0]} b={v[1]} alpha={v[2]}")
        print("  residual =", tuple(str(u) for u in f))


if __name__ == "__main__":
    main()
