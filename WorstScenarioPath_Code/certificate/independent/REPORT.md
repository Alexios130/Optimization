# Rigorous `n=3` path certificate: findings and proof architecture

## Scope

This report concerns the exact rational instance stated in the separately
submitted manuscript:

\[
n=3,\quad p=1,\quad S=2,\quad \lambda=\frac12,
\]

\[
x=\left(\frac72,3,4\right),\qquad y=(-1,1,-1),
\]

\[
\pi_1=\left(\frac{12}{19},\frac1{19},\frac6{19}\right),
\qquad
\pi_2=\left(\frac18,\frac1{16},\frac{13}{16}\right).
\]

A key feature of the exact path is that it has **two nontrivial tie
intervals**, rather than merely two isolated tie radii.

## 1. Exact reduction to two scalar branch gaps

For each scenario define

\[
J_s(\theta;\rho)
=L_s(\theta;\rho)+\frac{\lambda}{2}\|\theta\|_2^2,
\qquad
\theta_s(\rho)=\arg\min_\theta J_s(\theta;\rho),
\]

and define the two branch diagnostics

\[
\Delta_s(\rho)
=L_1(\theta_s(\rho);\rho)-L_2(\theta_s(\rho);\rho).
\]

Strong convexity gives

\[
\begin{aligned}
\Delta_2-\Delta_1
&=[J_1(\theta_2)-J_1(\theta_1)]
 +[J_2(\theta_1)-J_2(\theta_2)]\\
&\geq \lambda\|\theta_2-\theta_1\|_2^2\geq0.
\end{aligned}
\]

Consequently, for every fixed radius,

\[
I(\rho)=
\begin{cases}
\{1\},&\Delta_1(\rho)>0,\\
\{2\},&\Delta_2(\rho)<0,\\
\{1,2\},&\Delta_1(\rho)\leq0\leq\Delta_2(\rho).
\end{cases}
\]

For example, if \(\Delta_1>0\), then

\[
F(\theta_1)=J_1(\theta_1)
\leq J_1(\theta)\leq F(\theta)
\]

for every \(\theta\).  Hence \(\theta_1\) is the unique minimax optimizer
and scenario 1 is strictly active.  The scenario-2 case is symmetric.  If
\(\Delta_1\leq0\leq\Delta_2\), neither scenario can be uniquely active;
because there are exactly two scenarios, they tie.

This lemma is the cleanest way to certify the whole path.  It removes the
need to assume or numerically trace a multiplier path on the tie intervals.

## 2. Smooth branch equations

The optimizer cannot have \(\beta>0\).  A short analytic proof is available:

1. Any KKT mixture of the two scenarios assigns at most \(1/16\) weight to
   the sole positive-label observation.  If \(\beta>0\) and \(b\geq0\), the
   intercept stationarity residual is strictly positive, so necessarily
   \(b<0\).
2. For \(0\leq\rho\leq3\), subtract \((3-\rho)\) times the intercept KKT
   equation from the slope KKT equation.  The positive-label contribution
   cancels, while

   \[
   \lambda[\beta-(3-\rho)b]
   +\sum_{y_i=-1}w_i\sigma(q_i)(x_i-3+2\rho)>0,
   \]

   contradicting stationarity.
3. If \(\rho\geq3\), every slope coefficient
   \(-y_ix_i+\rho\) is nonnegative and the negative-label contributions are
   strictly positive, again contradicting slope stationarity.

Thus before the zero-slope threshold one may use \(|\beta|=-\beta\).  Put

\[
q_i=-y_i(x_i\beta+b)-\rho\beta,\qquad
z_i=-y_ix_i-\rho,\qquad t_i=-y_i,
\]

\[
a_i=\sigma(q_i),\qquad h_i=a_i(1-a_i).
\]

The scenario-\(s\) branch is determined by

\[
E_{\beta,s}=\frac12\beta+\sum_i\pi_{si}a_i z_i=0,
\qquad
E_{b,s}=\frac12b+\sum_i\pi_{si}a_i t_i=0.
\]

Its Hessian is

\[
H_s=
\begin{pmatrix}
\frac12+\sum_i\pi_{si}h_i z_i^2&
\sum_i\pi_{si}h_i z_it_i\\
\sum_i\pi_{si}h_i z_it_i&
\frac12+\sum_i\pi_{si}h_i
\end{pmatrix}
\succeq\frac12I.
\]

Hence the branch is differentiable while \(\beta_s<0\), and

\[
\theta_s'=-H_s^{-1}\partial_\rho E_s.
\]

These formulas are evaluated outward in the event derivative checks.

The exact scenario difference is

\[
\pi_1-\pi_2=\frac1{304}(154,-3,-151),
\]

so

\[
304\Delta
=154\log(1+e^{q_1})-3\log(1+e^{q_2})
-151\log(1+e^{q_3}).
\]

The verifier does not need an interval logarithm.  It uses the equivalent
comparison

\[
(1+e^{q_1})^{154}
\mathrel{\gtrless}
(1+e^{q_2})^3(1+e^{q_3})^{151}.
\]

## 3. Four rigorously isolated events

The exact-rational checker proves that the four and only four relevant
branch zeros lie in the following rational intervals:

| event | branch equation | rigorous radius enclosure | sign of derivative |
|---|---|---:|---:|
| \(\rho_1\) | \(\Delta_2=0\) | \((0.524565,0.524566)\) | positive |
| \(\rho_2\) | \(\Delta_1=0\) | \((1.036949,1.036950)\) | positive |
| \(\rho_3\) | \(\Delta_1=0\) | \((2.696120,2.696121)\) | negative |
| \(\rho_4\) | \(\Delta_2=0\) | \((2.832429,2.832430)\) | negative |

On wider rational windows containing these narrow enclosures, the checker
obtained

\[
\begin{array}{c|c|c}
\text{event}&\text{monotonicity window}&\text{outward interval for }\Delta_s'\\ \hline
1&[0.5245,0.5246]&[0.00281407,0.00305716]\\
2&[1.0369,1.0370]&[0.00383586,0.00433250]\\
3&[2.6960,2.6962]&[-0.0168025,-0.0102149]\\
4&[2.8323,2.8325]&[-0.0174878,-0.0107319].
\end{array}
\]

The gap signs at the two endpoints of every narrow radius enclosure are
opposite.  The strict derivative enclosure on the containing window proves
existence and uniqueness of the zero and its crossing orientation.

For reference only, the high-precision centres are

\[
\begin{aligned}
\rho_1&\approx0.5245653271531518056719427285,\\
\rho_2&\approx1.036949195872801538478249203,\\
\rho_3&\approx2.696120891018807272210231760,\\
\rho_4&\approx2.832429831560884088881420461.
\end{aligned}
\]

The proof uses the rational enclosures, not these decimal centres.

## 4. Exhaustive no-other-zero certificate

For a rational radius cell \(R\), the checker constructs a rational affine
predictor \(\bar\theta_s(\rho)\).  It outward-evaluates the branch KKT
residual over the whole cell.  If

\[
\|E_s(\bar\theta_s(\rho),\rho)\|_1\leq\varepsilon_R,
\]

then strong convexity gives, uniformly on that cell,

\[
\|\theta_s(\rho)-\bar\theta_s(\rho)\|_2
\leq\frac{\varepsilon_R}{\lambda}.
\]

The checker inflates both coordinates by this amount and evaluates the
exponential comparison above on the resulting rational box.  A cell is
accepted only if zero is excluded.

The completed adaptive run accepted the following exhaustive covers.  The
infinite constant tail is established analytically in Section 5 below; it is
not counted among these finite branch cells.

| certified interval | assertion | accepted rational cells |
|---|---:|---:|
| \([0,0.5245]\) | \(\Delta_2<0\) | 3,828 |
| \([0.5246,2.8323]\) | \(\Delta_2>0\) | 16,987 |
| \([2.8325,3.1890]\) | \(\Delta_2<0\) | 10,507 |
| \([0,1.0369]\) | \(\Delta_1<0\) | 4,995 |
| \([1.0370,2.6960]\) | \(\Delta_1>0\) | 17,343 |
| \([2.6962,2.8325]\) | \(\Delta_1<0\) | 8,841 |

Total: **62,501 accepted exact-rational sign cells**.  Together with the
four strict-derivative event windows, this rules out every additional branch
zero in the part of the path where a switch could affect the active set.

The floating-point solver is used only to propose predictors.  Every
accepted statement is subsequently checked with exact rational interval
arithmetic.  A poor floating-point proposal can make the program fail or
subdivide further, but it cannot make an invalid cell pass.

## 5. Exact zero-slope tail

At \(\beta=0\) and \(b<0\),

\[
L_1(0,b;\rho)-L_2(0,b;\rho)=\frac{3b}{304}<0,
\]

so scenario 2 is strictly active.  Its intercept-only optimum is the unique
root of

\[
\frac12 b_\infty+\sigma(b_\infty)-\frac1{16}=0.
\]

Set \(a_\infty=\sigma(b_\infty)\).  The slope subdifferential of \(J_2\) at
\((0,b_\infty)\) is

\[
[c_\infty-\rho H_\infty,c_\infty+\rho H_\infty],
\]

where

\[
c_\infty=\frac{62a_\infty-3}{16},
\qquad
H_\infty=\frac{14a_\infty+1}{16}.
\]

Therefore the exact stabilization threshold is

\[
\rho_\infty
=\frac{c_\infty}{H_\infty}
=\frac{62a_\infty-3}{14a_\infty+1}.
\]

The rational checker proves

\[
b_\infty\in
[-0.5888140329295523,-0.5888140329295522]
\]

and

\[
\rho_\infty\in
[3.18979449782973638,3.18979449782973749].
\]

For every \(\rho\geq\rho_\infty\),

\[
\theta^\ast(\rho)=(0,b_\infty),\qquad I(\rho)=\{2\}.
\]

To bridge the smooth branch certificate to the exact threshold, the checker
uses the rational candidate \((0,-0.58881403292955225)\) and the selected
subgradient \(c-\rho H\).  On

\[
[3.1890,\rho_\infty]
\]

the residual-to-distance transfer gives the strict upper bound

\[
\Delta_2(\rho)<-0.0010258.
\]

Thus there is no unverified interval adjacent to the kink.

## 6. Complete active-set theorem

Let \(\rho_1<\rho_2<\rho_3<\rho_4\) denote the unique roots isolated above.
Then the certified path is

\[
I(\rho)=
\begin{cases}
\{2\},&0\leq\rho<\rho_1,\\
\{1,2\},&\rho_1\leq\rho\leq\rho_2,\\
\{1\},&\rho_2<\rho<\rho_3,\\
\{1,2\},&\rho_3\leq\rho\leq\rho_4,\\
\{2\},&\rho>\rho_4.
\end{cases}
\]

In particular, the full tie set and activity regions are

\[
T=[\rho_1,\rho_2]\cup[\rho_3,\rho_4],
\]

\[
A_1=[\rho_1,\rho_4],
\qquad
A_2=[0,\rho_2]\cup[\rho_3,\infty).
\]

Scenario 2 therefore has exactly two connected activity components.

## 7. Minimality corollary

The completed theorem for \(n\leq2\) says that every activity region is
connected in the declared class.  The certified rational \(n=3\) instance
has

\[
A_2=[0,\rho_2]\cup[\rho_3,\infty),
\qquad \rho_2<\rho_3,
\]

which is disconnected.  Hence

\[
\boxed{n_{\min}=3}
\]

within the precise class

\[
p=1,\quad S=2,\quad
\text{logistic loss, normalized nonnegative scenario rows,}
\]

\[
\text{observation-wise uncertainty, and full ridge regularization.}
\]

This must not be described as minimality over arbitrary feature dimension or
arbitrary number of scenarios.

## 8. Code and verification status

The independent implementation consists of:

- `compute_decimal.py`: 90-digit exploratory centres and derivatives.  This
  file is diagnostic and is **not** a proof.
- `rational_interval_check.py`: exact fixed-point/rational outward checker.

The publication package also contains a separate frozen machine-readable
certificate and a compact primary verifier in
`certificate/verify_n3_certificate.py`.  That verifier uses 256-bit
outward-rounded Arb arithmetic through `python-flint`, has no NumPy or SciPy
dependency, and validates 23,973 path cells.  The fixed-lattice program in
this directory was written independently and validates 62,501 cells using
exact rational endpoints and proved exponential remainders.  The two
implementations therefore provide distinct rigorous checks of the same four
events, sign chart, bridge, and stabilization tail.  Their recorded outputs
are under `verification/`, while the frozen certificate's SHA-256 digest is
listed in `REVIEWER_GUIDE.md`.
