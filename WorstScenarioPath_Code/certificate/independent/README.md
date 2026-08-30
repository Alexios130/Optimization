# Independent verification programs

`rational_interval_check.py` is a separately written rigorous cross-check. It
uses exact rational interval endpoints, lattice denominator `Q=10^18`, and 32
Taylor terms with proved geometric remainder bounds for the exponential. Its
successful final summary includes `status: PASS`, `lattice_denominator:
1000000000000000000`, `exp_terms: 32`, and
`validated_rational_sign_cells: 62501`. In the latest tested environment it
took 611.93 seconds as one CPU-bound process and printed the summary only after
completing the exhaustive cover; runtime is hardware-dependent.

The program prints event brackets exactly as fractions. Any decimal interval
or width printed for readability is labelled as a display summary. The proof
logic uses only the exact rational endpoints stored internally.

`compute_decimal.py` uses 90-digit decimal arithmetic to obtain exploratory
centres and derivatives. It is useful diagnostically, but it is not
outward-rounded and is not used as a proof.
