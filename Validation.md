# Validation Status

Validation is carried out by comparison with the output of the public Web version of ISODISTORT, together with mathematical checks that are independent of the Web implementation.
The figures in this document are a snapshot **as of October 15, 2026**.
For each Web branch, they refer to the most recent Validation output produced by a local runtime that has been landed on main.

## Verdicts

The unit of validation is a pair consisting of one input structure and one order-parameter direction (OPD).

- **Strict pass**: The display setting of the OPD, the child structure, and the mode types, mode counts, atom correspondence, normalization factors, and vectors agree exactly with the Web output, including their order.
- **Physical pass**: Even when there are differences in the displayed mode order, in signs, in the choice of basis within the same mode space, or in equivalent crystallographic settings, the child structure and the physical mode space that is spanned are equivalent, and only purely conventional differences, such as the choice of basis, remain. Every Strict pass is also included in Physical pass.
- **Physical fail**: The OPD does not exist, the child structure is physically different, required modes are missing, or the mode vectors do not span the same physical space — that is, the result is physically different from the reference implementation. Potential candidate Physical passes, in which neither the reference implementation nor the local implementation violates mathematics but a non-trivial physical equivalence could not be decided, are also treated conservatively as fail.

## Definition of the Valid Population

Two necessary-condition checks are applied to each output.

- **Mode basis:** exact rational linear independence of the complete displayed mode basis
- **Group invariance:** invariance of the displayed mode fields under the reported subgroup operations, in the settings maintained by the verifier

The Validation population consists of the branches for which both Web and Local have been computed and for which the Web output was not refuted by either check.

## Reproducibility

The public checkout includes a `Branches/` bundle containing each canonical input, selected OPD, and recorded comparison and mathematics statuses.
The standalone `Verification/` package contains the exact parser, mathematical checks, Strict/Physical comparator, and result projection used by the private monorepo Validation operation.
A third party can select a branch, run APOSTRUCT, submit the same input to the public ISODISTORT Web application, and apply those checks to the two outputs without depending on the private sweep/store schema.

## Two Validation Populations

### MP pool: broad stress population

A population based on structures derived from the Materials Project, broadly including fixed K, parametric K, multi-K, and general OPDs.

#### Selection of the structure pool

The CIF pool was created from the Materials Project.
Structures in the same space group whose combinations of occupied Wyckoff site multiplicities and letters coincide are regarded as the same topology, and only one entry was acquired per `(space group, occupied Wyckoff pattern)`.
Structures that differ only by elemental substitution, and structures that differ only in the free coordinates of the same Wyckoff pattern, were removed as duplicates at this stage.
The pool obtained in this way contains **19,086 CIFs**.
Of these, **2,557 CIFs** have actually been used in the Web populations described below at this point.

#### Sampling of inputs

The current population is not a single simple random sampling, but the union of several reproducible campaigns.
Each campaign records its seed, filters, target count, and the ordered input IDs that were generated.

1. CIFs are filtered by space group, crystal system, centering, number of sites, number of atoms, and so on.
2. Depending on the campaign, a bucket of space group, crystal system, or centering is selected uniformly first, and then a CIF within that bucket is selected uniformly. For comparison, campaigns that select from all candidates without using buckets are also included.
3. Strain, and the enablement or disablement of displacive and magnetic modes for each element, are selected. When random selection is specified, each element is selected independently, but strata fixed to ordinary-only, magnetic-only, mixed, and so on are also included.
4. The number of K slots is selected in the range from 1 to 4, and the number and arrangement of fixed and parametric slots are selected. Multi-parametric configurations, different positions of the parametric slots, multi-dimensional parametric K, and multi-K, which were thin in the initial population, were explicitly supplemented in later campaigns.
5. For each slot, only the K/IR combinations that are valid in Source for that space group, for the Wyckoff rows of the selected elements, and for the mode kind are enumerated, and the selection is made from those candidates.
6. The values of parametric K are enumerated as irreducible rational numbers within the Source domain with a denominator of at most 6, and are selected from the set obtained after removing duplicates of the canonical values.
7. After the OPD list has been obtained from the Web, the 1-based ordinals explicitly specified by the campaign are collected. When several OPDs are collected from a single input, each is counted as a separate branch.

The pass rates below are therefore **the achievement rates of a stratified stress test**, centred on random sampling but with the difficult configurations that were lacking deliberately made heavier.

#### Results

| Item | Count |
|---|---:|
| Valid Validation population | 3,040 |
| Strict pass | 2,599 |
| Physical pass | 3,035 |
| Physical fail | 5 |

- **Physical pass: 3,035 / 3,040 = 99.84%**
- **Strict pass: 2,599 / 3,040 = 85.49%**

#### By K-signature

`F` denotes a fixed K and `P` denotes a parametric K. Every column uses the Validation population defined above.

| K | Validation population | Physical pass | Strict pass |
|---|---:|---:|---:|
| F | 788 | 788 (100.00%) | 734 (93.15%) |
| FF | 297 | 297 (100.00%) | 252 (84.85%) |
| FFF | 189 | 189 (100.00%) | 167 (88.36%) |
| FFFF | 202 | 202 (100.00%) | 186 (92.08%) |
| P | 687 | 686 (99.85%) | 599 (87.19%) |
| PP | 20 | 18 (90.00%) | 7 (35.00%) |
| FP | 126 | 126 (100.00%) | 106 (84.13%) |
| FFP | 453 | 453 (100.00%) | 363 (80.13%) |
| FFFP | 246 | 246 (100.00%) | 168 (68.29%) |
| FPP | 31 | 29 (93.55%) | 17 (54.84%) |
| FFPP | 1 | 1 (100.00%) | 0 (0.00%) |

### MAGNDATA: practical population

A population based on the experimentally realized magnetic structures recorded in MAGNDATA.
Its purpose is to examine the reliability against practical inputs.

#### Selection of structures and inputs

Rather than entering the MAGNDATA mcif directly into Local, an ordinary CIF of the same parent structure was obtained from the Materials Project, and then from the Crystallography Open Database (COD).

For the correspondence of the parent CIF, the parent SG of MAGNDATA was required, and only cases in which uniqueness such as the following could be confirmed were adopted.
No arbitrary single entry was chosen from among ambiguous candidates.

- Exactly one candidate matches the chemical formula of the MAGNDATA/mcif, ignoring the order of elements: 495 from MP
- Uniquely matches the chemical formula on the MAGNDATA top page: 153 from MP
- Uniquely corresponds by way of a recorded ICSD ID: 344 from MP
- Uniquely corresponds by combining the chemical formula, DOI, citation, and element set with the parent SG: 59 from COD

With the above, **1,051** parent CIFs could be matched.

For input generation, only information that can be read directly from the mcif, rather than inferred, is used.

- The primary IR is mapped to a unique IR label in Source, and the K belonging to that IR is placed in a slot
- For the child magnetic space group, the BNS number recorded in the mcif is used as it is
- If even one site of a given element carries a magnetic moment, that element as a whole is turned magnetic-on
- If the primary IR has an ordinary component, all elements are turned displacive-on, and strain is turned off

Excluding a missing primary IR or BNS, a non-unique Source label, and invalid mcif records that cannot be submitted to the Web as the same input and OPD, the final Validation population is **927** entries.
It consists of 458 zero-K, 28 nonzero Type I/III, 412 nonzero Type IV, 19 two-K, and 10 three-or-more-K entries.

| Item | Count |
|---|---:|
| Population / inputs / CIFs | 927 / 927 / 927 |
| Validated and judged | 927 |
| Strict pass | 892 |
| Physical pass | 927 |
| Physical fail | 0 |

- **Physical pass: 927 / 927 = 100.00%**
- **Strict pass: 892 / 927 = 96.22%**

#### By K-signature

| K | Population | Physical pass | Strict pass |
|---|---:|---:|---:|
| F | 798 | 798 (100.00%) | 767 (96.12%) |
| FF | 84 | 84 (100.00%) | 81 (96.43%) |
| P | 40 | 40 (100.00%) | 39 (97.50%) |
| PP | 3 | 3 (100.00%) | 3 (100.00%) |
| FP | 1 | 1 (100.00%) | 1 (100.00%) |
| FPP | 1 | 1 (100.00%) | 1 (100.00%) |

## Limits of Validation

Among the 4,141 successful Web branches, 4,139 have a recorded Local comparison (the other 2 are Local operational failures).
Of these, **171** are refuted by the current mathematical checks and **1** is indeterminate, so all 172 are excluded from the valid population.
Agreement with an excluded Web result is not counted as correctness evidence.

## Computational Performance

The performance population consists of the branches within the valid Validation population for which `wall_s` is recorded for both Web and Local.
The Web time is the sum of the successful OPD-list request and the selected Complete Mode Details request.
The Local time is the wall time of one call of the local executable `apo` per branch (including its start-up, measured during a sweep with six workers).
Because the execution environments differ, the figures below are reference values and not a direct comparison of speed.

| Performance population | Local median | Web median | Local p95 | Web p95 |
|---:|---:|---:|---:|---:|
| 2,081 / 3,967 | 0.03 s | 8.27 s | 1.11 s | 30.89 s |

```text
          0.01        0.1         1           10          100 s
Local         [=│=====]------------|
        p5=0.02  Q1=0.02  median=0.03  Q3=0.11  p95=1.11
Web                                         [│]-----|
        p5=7.31  Q1=7.48  median=8.27  Q3=10.44  p95=30.89
```

### Performance by K-signature

| K | Performance population | Local median | Web median | Local p95 | Web p95 |
|---|---:|---:|---:|---:|---:|
| F | 951 | 0.02 s | 7.48 s | 0.04 s | 7.94 s |
| FF | 164 | 0.02 s | 8.80 s | 0.09 s | 10.68 s |
| FFF | 53 | 0.04 s | 9.24 s | 0.31 s | 37.38 s |
| FFFF | 47 | 0.05 s | 9.30 s | 0.30 s | 20.05 s |
| P | 235 | 0.06 s | 8.48 s | 1.39 s | 32.51 s |
| PP | 23 | 0.25 s | 13.47 s | 1.70 s | 37.34 s |
| FP | 127 | 0.11 s | 10.58 s | 1.82 s | 41.10 s |
| FFP | 240 | 0.17 s | 11.54 s | 1.97 s | 47.92 s |
| FFFP | 208 | 0.24 s | 14.03 s | 2.42 s | 50.26 s |
| FPP | 32 | 0.19 s | 13.41 s | 1.94 s | 53.72 s |
| FFPP | 1 | 1.67 s | 66.64 s | 1.67 s | 66.64 s |

Even for the stress cases that include parametric K and multi-K, the Local p95 stays within 2.5 s.
