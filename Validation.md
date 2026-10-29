# Validation

## Populations

### MP pool

A broad stress-test population using CIFs obtained from the Materials Project

### MAGNDATA

A practical population reproducing magnetic structures recorded in MAGNDATA

## Mathematical checks

[`verify_local_output.py`](verification/verify_local_output.py) · [`verify_web_text.py`](verification/verify_web_text.py)

Each APOSTRUCT and ISODISTORT output is checked for:

- Exact rational linear independence of the displayed mode basis
- Invariance of the displayed mode fields under the reported subgroup operations
- APOSTRUCT outputs only: each mode lies in the Bloch components of the star of its labelled k, the symmetry of the primary distortion is exactly the reported subgroup (proved, or indeterminate when the displacements do not resolve the order parameter), and modes of different irreps at one k star have orthogonal little-group orbits (no mixing of irreps). These take the parent operations from the atoms, without tables

## Output comparison

[`verify_comparison.py`](verification/verify_comparison.py)

APOSTRUCT and ISODISTORT outputs are compared.

- **Identical**: Display settings, child structures, mode types and counts, atom correspondence, normalization factors and vectors match, including conventions such as ordering.
- **Equivalent**: Equivalence of the child structures and mode spaces is established through setting or mode-basis transformations. Names (space groups, irreps, k points, Wyckoff letters, OPD labels) must match. Includes Identical results.
- **Unproven**: Equivalence could not be established.

APOSTRUCT aims at Equivalent, not Identical. Names follow ISODISTORT, but choices inside an equivalence class (which domain, the basis and origin within the standard setting, signs and bases of modes) follow APOSTRUCT's own short rules and may differ. Where the ISODISTORT output fails the mathematical checks, APOSTRUCT does not follow it; `web_math` in the records marks these branches.

Execution failures are recorded separately from comparison verdicts.

## Records

Calculation, check and comparison records for each branch (an input and OPD pair) are in [`branches/`](branches/).

- [`inputs.jsonl.gz`](branches/inputs.jsonl.gz): Inputs
- [`results.jsonl.gz`](branches/results.jsonl.gz): Verdicts and mathematical check statuses for each branch
- [`manifest.json`](branches/manifest.json): Branch and input counts, file hashes and verifier versions
