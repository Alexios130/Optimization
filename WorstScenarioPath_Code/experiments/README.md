# Numerical experiments

These programs and JSON outputs document the computational research process.
They are deliberately separated from `certificate/` because dense grids,
Brent refinement, Decimal Newton iterations, small residuals, and nonzero
numerical derivatives do not exclude every hidden event.

The finite `L=16` family files illustrate the analytic lower-bound mechanism;
they do not prove the all-`m` theorem. Exact commands are listed in
`REPRODUCIBILITY_MANIFEST.md`.

`grid_resolution_study.py` postprocesses the recorded 90-digit event radii to
measure the open tie intervals and to test direct active-set sampling on a
uniform grid.  In particular, its statement that a grid misses a tie interval
means that no grid node lies strictly inside that interval; it does not mean
that the branch-gap continuation programs fail to locate the two boundary
events.  The generated evidence is `grid_resolution_study.json`.

`make_n3_path_figure.py` reads the exact rational model and certified event
enclosures from `certificate/`, samples the two branch gaps with the existing
two-scenario solver, and produces `n3_path_figure_data.json` together with the
vector figure `article/figures/n3_certified_path.pdf`.  The event markers and
state sequence in that figure are certified; the smooth plotted curves are a
double-precision visualization.
