# Recorded verification

`arb_verification_output.txt` records the audited primary Arb run, and
`fixed_lattice_output.txt` records the independently implemented rational
interval run. Re-run both programs rather than trusting these text files; the
executable verifiers and frozen certificate are the authoritative
reproducibility objects.

The recorded independent output labels decimal intervals and widths as
display summaries. Its event brackets are exact fractions, and the checker
uses exact rational endpoints for every decision.
