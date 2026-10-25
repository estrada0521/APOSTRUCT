# Validation Status

Validation uses comparisons with ISODISTORT output and independent mathematical checks.

## Mathematical checks

ISODISTORT outputs are used for comparison only after passing these independent checks:

- Exact rational linear independence of the complete displayed mode basis
- Invariance of the displayed mode fields under the reported subgroup operations

## Verdicts

- **Identical**: The OPD display setting, child structure, mode types and counts, atom correspondence, normalization factors, and vectors match ISODISTORT output exactly, including BYU-specific conventions such as ordering.
- **Equivalent**: The child structure and the spanned mode space are proven equivalent to ISODISTORT output through explicit crystallographic and mode-space equivalences (mode order, signs, basis within the same mode space, equivalent settings). Every Identical result is also Equivalent.
- **Unproven**: The comparator could not prove equivalence. This includes physically different outputs and cases where the comparator cannot establish equivalence. It does not mean that the APOSTRUCT output is wrong.

## Validation Populations

### MP pool

A population based on structures from the Materials Project, broadly covering fixed and parametric K, multi-K, and general OPDs.

#### Selection of the structure pool

The CIF pool was created from the Materials Project.
Structures in the same space group with the same combination of occupied Wyckoff site multiplicities and letters were treated as the same topology, and only one entry was acquired per `(space group, occupied Wyckoff pattern)`.
Structures differing only in elemental substitution or in the free coordinates of the same Wyckoff pattern were removed as duplicates at this stage.

### MAGNDATA

This population is based on experimentally realized magnetic structures recorded in MAGNDATA. Parent CIFs are obtained from the Materials Project, followed by COD, requiring the same parent space group and a unique match through chemical formulas, ICSD IDs, citation information or related metadata.

Inputs are generated from the primary irreps and child magnetic space group BNS number in the mcif. Each primary irrep is mapped to a unique Source label and its k point selected. Magnetic modes are enabled for an entire element if any of its sites carries a magnetic moment. If a primary irrep has an ordinary component, displacements are enabled for all elements and strain is disabled.

Records with missing required information, ambiguous Source labels, or inputs and OPDs that cannot be submitted consistently to ISODISTORT are excluded.

## Reproducibility

[`branches/`](branches/) contains each branch's canonical input, selected OPD, comparison result, and mathematical check status.
[`verification/`](verification/) contains the verification tools needed to reproduce the results.

See [`inputs.jsonl.gz`](branches/inputs.jsonl.gz) for inputs, [`results.jsonl.gz`](branches/results.jsonl.gz) for verdicts, and [`manifest.json`](branches/manifest.json) for counts and verifier versions.
