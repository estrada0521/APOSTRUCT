# Validation Status

Validation uses comparisons with ISODISTORT output and independent mathematical checks.
The figures in this document are a snapshot **as of October 23, 2026**.

## Definition of the Valid Population

Only Web outputs that pass both independent checks are included in the valid population:

- Exact rational linear independence of the complete displayed mode basis
- Invariance of the displayed mode fields under the reported subgroup operations

Currently, **215 outputs** are excluded. Most involve parametric K with a parameter greater than 1/2.

## Verdicts

Results are divided into three categories.

- **Identical**: The OPD display setting, child structure, mode types and counts, atom correspondence, normalization factors, and vectors match the Web output exactly, including BYU-specific conventions such as ordering.
- **Equivalent**: The child structure and the spanned mode space are proven equivalent to the Web output through explicit crystallographic and mode-space equivalences (mode order, signs, basis within the same mode space, equivalent settings). Every Identical result is also Equivalent.
- **Unproven**: The comparator could not prove equivalence. This includes physically different outputs and cases where the comparator cannot establish equivalence. It does not mean that the APOSTRUCT output is wrong.

## Validation Populations

### 1. MP pool: broad stress population

A population based on structures from the Materials Project, broadly covering fixed and parametric K, multi-K, and general OPDs.

#### Selection of the structure pool

The CIF pool was created from the Materials Project.
Structures in the same space group with the same combination of occupied Wyckoff site multiplicities and letters were treated as the same topology, and only one entry was acquired per `(space group, occupied Wyckoff pattern)`.
Structures differing only in elemental substitution or in the free coordinates of the same Wyckoff pattern were removed as duplicates at this stage.
The resulting pool contains **19,086 CIFs**.
Of these, **5,483 CIFs** have been used in the populations below so far.

#### Results

| Item | Count |
|---|---:|
| All | 9,634 |
| Identical | 8,844 |
| Equivalent | 9,630 |
| Unproven | 4 |

- **Equivalent: 9,630 / 9,634 = 99.96%**
- **Identical: 8,844 / 9,634 = 91.80%**

#### By K-signature

`F` denotes a fixed K and `P` denotes a parametric K.

| K | Validation population | Equivalent | Identical |
|---|---:|---:|---:|
| F | 6,233 | 6,233 (100.00%) | 5,913 (94.87%) |
| FF | 846 | 846 (100.00%) | 742 (87.71%) |
| FFF | 189 | 189 (100.00%) | 167 (88.36%) |
| FFFF | 202 | 202 (100.00%) | 186 (92.08%) |
| P | 1,287 | 1,287 (100.00%) | 1,169 (90.83%) |
| PP | 20 | 18 (90.00%) | 7 (35.00%) |
| FP | 126 | 126 (100.00%) | 108 (85.71%) |
| FFP | 453 | 453 (100.00%) | 364 (80.35%) |
| FFFP | 246 | 246 (100.00%) | 171 (69.51%) |
| FPP | 31 | 29 (93.55%) | 17 (54.84%) |
| FFPP | 1 | 1 (100.00%) | 0 (0.00%) |

### MAGNDATA: practical population

A population based on experimentally realized magnetic structures recorded in MAGNDATA.
Its purpose is to assess reliability on practical inputs.

#### Selection of structures and inputs

Rather than submitting MAGNDATA mcif files directly to Local, ordinary CIFs of the same parent structures were obtained from the Materials Project, followed by the Crystallography Open Database (COD).

Each parent CIF had to match the MAGNDATA parent space group. Only uniquely identified matches were accepted:

- Exactly one candidate matches the MAGNDATA/mcif chemical formula, ignoring element order: 495 from MP
- A unique match to the chemical formula on the MAGNDATA top page: 153 from MP
- A unique match through a recorded ICSD ID: 344 from MP
- A unique match using the chemical formula, DOI, citation, and element set together with the parent space group: 59 from COD

This identified **1,051 parent CIFs**.

Input generation uses only information read directly from the mcif:

- Map each primary IR to a unique IR label in Source and place its K in a slot
- Use the child magnetic space group's BNS number as recorded in the mcif
- Enable magnetic modes for an entire element if any of its sites carries a magnetic moment
- If the primary IR has an ordinary component, enable displacive modes for all elements and disable strain

After excluding missing primary IR or BNS information, non-unique Source labels, and invalid mcif records that cannot be submitted to the Web as the same input and OPD, the final population contains **951 entries**.
It consists of 862 one-K, 88 two-K, and 1 three-K entries.

| Item | Count |
|---|---:|
| All | 951 |
| Identical | 914 |
| Equivalent | 951 |
| Unproven | 0 |

- **Equivalent: 951 / 951 = 100.00%**
- **Identical: 914 / 951 = 96.11%**

#### By K-signature

| K | Population | Equivalent | Identical |
|---|---:|---:|---:|
| F | 821 | 821 (100.00%) | 788 (95.98%) |
| FF | 84 | 84 (100.00%) | 81 (96.43%) |
| P | 41 | 41 (100.00%) | 40 (97.56%) |
| PP | 3 | 3 (100.00%) | 3 (100.00%) |
| FP | 1 | 1 (100.00%) | 1 (100.00%) |
| FPP | 1 | 1 (100.00%) | 1 (100.00%) |

## Reproducibility

`branches/` contains each branch's canonical input, selected OPD, comparison result, and mathematical check status.
`verification/` contains the verification tools needed to reproduce the results above.
