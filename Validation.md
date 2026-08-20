# Validation Status

Validation is carried out by comparison with the output of the public Web version of ISODISTORT, together with mathematical checks that are independent of the Web implementation.
The figures in this document are a snapshot **as of August 12, 2026**.
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

The public checkout includes a `Branches/` bundle containing each canonical
input, selected OPD, and recorded comparison and mathematics statuses. It does
not include saved Web or Local output payloads or a Web automation wrapper.
The standalone `Verification/` package contains the exact parser, mathematical
checks, Strict/Physical comparator, and result projection used by the private
monorepo Validation operation. A third party can select a branch, run APOSTRUCT,
submit the same input to the public ISODISTORT Web application, and apply those
checks to the two outputs without depending on the private sweep/store schema.

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
| Strict pass | 1,897 |
| Physical pass | 3,034 |
| Physical fail | 6 |

- **Physical pass: 3,034 / 3,040 = 99.80%**
- **Strict pass: 1,897 / 3,040 = 62.40%**

#### By K-signature

`F` denotes a fixed K and `P` denotes a parametric K. Every column uses the Validation population defined above.

| K | Validation population | Physical pass | Strict pass |
|---|---:|---:|---:|
| F | 788 | 788 (100.00%) | 678 (86.04%) |
| FF | 297 | 297 (100.00%) | 217 (73.06%) |
| FFF | 189 | 189 (100.00%) | 121 (64.02%) |
| FFFF | 202 | 202 (100.00%) | 138 (68.32%) |
| P | 687 | 686 (99.85%) | 439 (63.90%) |
| PP | 20 | 18 (90.00%) | 1 (5.00%) |
| FP | 126 | 126 (100.00%) | 62 (49.21%) |
| FFP | 453 | 452 (99.78%) | 154 (34.00%) |
| FFFP | 246 | 246 (100.00%) | 81 (32.93%) |
| FPP | 31 | 29 (93.55%) | 6 (19.35%) |
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
| Strict pass | 882 |
| Physical pass | 927 |
| Physical fail | 0 |

- **Physical pass: 927 / 927 = 100.00%**
- **Strict pass: 882 / 927 = 95.15%**

#### By K-signature

| K | Population | Physical pass | Strict pass |
|---|---:|---:|---:|
| F | 798 | 798 (100.00%) | 760 (95.24%) |
| FF | 84 | 84 (100.00%) | 81 (96.43%) |
| P | 40 | 40 (100.00%) | 36 (90.00%) |
| PP | 3 | 3 (100.00%) | 3 (100.00%) |
| FP | 1 | 1 (100.00%) | 1 (100.00%) |
| FPP | 1 | 1 (100.00%) | 1 (100.00%) |

## Limits of Validation

**65** cases have been confirmed so far that fall outside the valid population because of a mathematical violation on the Web side.
Among the corresponding Local outputs, 64 satisfied the two current checks and 1 was refuted.
Since this does not constitute a proof of Local, these cases are treated as subjects requiring verification by the stronger mathematical checks to come.

## Computational Performance

The performance population consists of the branches within the valid Validation population for which `wall_s` is recorded for both Web and Local.
The Web time measures the interval from the CIF upload to the conversion of the Complete Mode Details into text.
The Local time measures the local mode-details computation itself.
Because the execution environments differ, the figures below are reference values and not a direct comparison of speed.

| Performance population | Local median | Web median | Local p95 | Web p95 |
|---:|---:|---:|---:|---:|
| 2,075 / 3,967 | 0.69 s | 8.25 s | 42.04 s | 32.29 s |

```text
          0.1                 1                  10                100 s
Local        |----[=======│=============]-------------------|
        p5=0.14  Q1=0.28  median=0.69  Q3=3.78  p95=42.04
Web                                          |│=]---------|
        p5=7.31  Q1=7.48  median=8.25  Q3=10.52  p95=32.29
```

### Performance by K-signature

| K | Performance population | Local median | Web median | Local p95 | Web p95 |
|---|---:|---:|---:|---:|---:|
| F | 951 | 0.29 s | 7.48 s | 1.00 s | 7.94 s |
| FF | 164 | 0.45 s | 8.81 s | 2.51 s | 10.68 s |
| FFF | 53 | 1.29 s | 9.24 s | 15.33 s | 37.38 s |
| FFFF | 47 | 1.76 s | 9.30 s | 11.24 s | 20.05 s |
| P | 234 | 2.07 s | 8.48 s | 51.45 s | 30.86 s |
| PP | 21 | 12.17 s | 13.47 s | 88.07 s | 37.34 s |
| FP | 127 | 3.77 s | 10.58 s | 77.02 s | 41.10 s |
| FFP | 239 | 5.84 s | 11.84 s | 61.47 s | 53.99 s |
| FFFP | 208 | 10.09 s | 14.24 s | 84.17 s | 62.02 s |
| FPP | 30 | 10.55 s | 13.83 s | 78.43 s | 53.72 s |
| FFPP | 1 | 86.27 s | 66.64 s | 86.27 s | 66.64 s |

A clear long tail remains on the Local side for the stress cases that include parametric K and multi-K.
