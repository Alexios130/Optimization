# Reproducibility code for worst-scenario switch paths

This repository contains only the computational reproducibility material for
Alexis Seferlis's article *Switch Complexity of Worst-Scenario Paths in Robust
Logistic Regression*.

Repository: <https://github.com/Alexios130/Optimization>

The manuscript source, manuscript PDF, bibliography, and rendered figure files
are intentionally **not** stored here. They belong to the journal submission,
not to the code repository.

## Included material

- `certificate/`: the frozen exact-rational `n=3` certificate, the primary
  outward-rounded Arb verifier, an independent exact-rational checker, and the
  proposal/diagnostic programs explicitly discussed in the article;
- `experiments/`: the five programs and machine-readable JSON outputs used for
  the finite-family, grid-resolution, and path-plot calculations;
- `verification/`: recorded outputs from the two rigorous verification routes;
- `REPRODUCIBILITY.md`: the exact program-to-claim map and evidentiary status;
- `reproduce.sh`: one entry point for the clean reruns; and
- `SHA256SUMS`: byte-level fingerprints for every archived file other than the
  checksum file itself.

There are exactly nine Python programs. No exploratory `n=4` candidate,
manuscript copy, LaTeX build file, rendered figure, internal audit, or virtual
environment is included.

## Quick verification

From the repository root on a system with Bash, Python 3, `venv`, and
`sha256sum`, run:

```bash
sha256sum -c SHA256SUMS
./reproduce.sh --quick
```

The quick mode verifies the theorem-level Arb certificate and regenerates the
grid-resolution JSON. For every archived computation, including the slow
independent rational verifier and the numerical finite-family reruns, use:

```bash
./reproduce.sh --full
```

The recorded independent-verifier runtime was 611.93 seconds on one CPU core;
runtime is hardware-dependent.

The path-plot program is retained because the article's figure caption names
both the program and its JSON output. A reviewer may regenerate the rendered
PDF locally with `./reproduce.sh --figure`; that generated PDF is intentionally
ignored and is not part of the repository.

