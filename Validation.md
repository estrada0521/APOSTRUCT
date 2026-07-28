# Validation Status of Distortropy

Distortropy is validated through both comparison with the public Web version of ISODISTORT and mathematical checks independent of the Web implementation.
The figures in this document are a snapshot as of **July 27, 2026**.

## Verdicts

The unit of validation is a pair consisting of one input structure and one order-parameter direction (OPD).

- **Strict pass**: The OPD display setting, child structure, mode types, mode counts, atom correspondence, normalization factors, and vectors agree exactly with the Web output within the accepted numerical precision.
- **Physical pass**: Even if the displayed mode order, signs, basis choice within the same mode space, or equivalent crystallographic setting differ, the child structure and the physical mode space spanned by the modes are equivalent, and the differences are purely conventional.
  Every Strict pass is also included in Physical pass.
- **Physical fail**: The OPD does not exist, the child structure is physically different, required modes are missing, or the mode vectors do not span the same physical space and are physically different from the reference implementation (as discussed below, this definition has a minor qualification because some outputs from the reference implementation do not satisfy mathematical requirements).

## Two Validation Populations

### MP pool: Broad Stress Population

This population is based on structures from Materials Project and broadly includes fixed K points, parametric K points, multi-K inputs, and general OPDs.

#### Selection of the structure pool

The CIF pool was created from Materials Project.
Structures in the same space group with the same combination of occupied Wyckoff-site multiplicities and letters were treated as the same topology, and only one structure was acquired for each `(space group, occupied Wyckoff pattern)`.
Structures differing only by element substitution, or only by free coordinates within the same Wyckoff pattern, were removed as duplicates at this stage.
The resulting pool contains **19,086 CIFs**.
As of this snapshot, 2,665 of those CIFs are used in the Web population below.

#### Input sampling

The current population is not a single simple random sample, but the union of multiple reproducible campaigns.
Each campaign records its seed, filters, target count, and generated ordered input IDs.

1. CIFs are filtered by properties such as space group, crystal system, centering, site count, and atom count.
2. Depending on the campaign, a space-group, crystal-system, or centering bucket is first selected uniformly, followed by uniform selection of a CIF within that bucket.
   Some comparison campaigns instead select from all candidates without using buckets.
3. Strain and the enabled or disabled displacive and magnetic modes for each element are selected.
   Under random selection, each element is selected independently, while the population also includes strata fixed to ordinary-only, magnetic-only, mixed, and similar conditions.
4. Between one and four K slots are selected, together with the number and placement of fixed and parametric K points.
   Multi-parametric inputs, alternate parametric-slot positions, multidimensional parametric K points, and multi-K inputs that were sparse in the initial population were explicitly supplemented in later campaigns.
5. For each slot, only K/IR combinations that are valid in Source for the space group, the Wyckoff rows occupied by the selected elements, and the mode kind are enumerated, and a candidate is selected from that set.
6. Parametric K values are selected by enumerating reduced rational numbers within the Source domain with denominators no greater than 6 and deduplicating their canonical values.
7. After the Web OPD list is obtained, the 1-based ordinal explicitly specified by the campaign is collected.
   If multiple OPDs are collected for one input, each is counted as a separate branch.

The pass rates below therefore represent the achievement rate of a **stratified stress test** centered on random sampling while intentionally increasing the weight of difficult, underrepresented shapes.

| Item | Count |
|---|---:|
| Web population | 3,206 branches |
| Inputs / CIFs | 3,094 / 2,665 |
| Validated | 3,195 |
| Judged | 3,148 |
| Strict pass | 1,142 |
| Physical-only pass | 1,807 |
| Physical fail | 199 |
| Unjudged | 47 |
| Not yet validated | 11 |

- **Physical pass: 2,949 / 3,148 = 93.68%**
- **Strict pass: 1,142 / 3,148 = 36.28%**
- If Unjudged cases are conservatively counted as failures, the Physical pass rate is **2,949 / 3,195 = 92.30%**.

#### By K signature

`F` denotes a fixed K point and `P` denotes a parametric K point.

| K composition (order-independent) | Population | Judged | Physical pass | Strict pass | Unjudged |
|---|---:|---:|---:|---:|---:|
| F | 788 | 788 | 782 (99.24%) | 512 (64.97%) | 0 |
| FF | 297 | 294 | 293 (99.66%) | 186 (63.27%) | 3 |
| FFF | 189 | 185 | 185 (100.00%) | 66 (35.68%) | 4 |
| FFFF | 202 | 200 | 200 (100.00%) | 70 (35.00%) | 2 |
| P | 761 | 748 | 685 (91.58%) | 212 (28.34%) | 13 |
| PP | 24 | 22 | 16 (72.73%) | 2 (9.09%) | 2 |
| FP | 141 | 140 | 123 (87.86%) | 32 (22.86%) | 1 |
| FFP | 490 | 477 | 423 (88.68%) | 50 (10.48%) | 13 |
| FFFP | 273 | 258 | 221 (85.66%) | 11 (4.26%) | 15 |
| FPP | 40 | 36 | 21 (58.33%) | 1 (2.78%) | 4 |

### MAGNDATA: Practical Population

This population is based on experimentally realized magnetic structures recorded in MAGNDATA.
Its purpose is to evaluate reliability on practical inputs.

#### Selection of structures and inputs

Rather than passing the MAGNDATA mcif directly to Local, a normal CIF for the same parent structure was obtained first from Materials Project and then from the Crystallography Open Database (COD).

The MAGNDATA parent space group was required when matching parent CIFs, and a record was adopted only when uniqueness could be established under one of the following conditions.
No arbitrary candidate was selected from an ambiguous set.

- Exactly one MP candidate matched the MAGNDATA/mcif chemical formula, ignoring element order: 495 cases
- A unique MP candidate matched the chemical formula on the MAGNDATA top page: 153 cases
- A unique MP candidate was identified by the recorded ICSD ID: 344 cases
- A unique COD candidate was identified by combining the chemical formula, DOI, citation, or element set with the parent space group: 59 cases

This procedure matched parent CIFs for **1,051 cases**.

Input generation uses only information that can be read directly from the mcif rather than inferred.

- The primary IR is mapped to a unique Source IR label, and the K point belonging to that IR is entered into the slot
- The BNS number recorded in the mcif is used directly as the child magnetic space group
- If at least one site of an element has a magnetic moment, magnetic mode is enabled for that entire element
- If the primary IR contains an ordinary component, displacive mode is enabled for all elements, and strain is disabled

After excluding cases with missing primary IRs, non-unique Source labels, missing BNS numbers, or inconsistencies with the Web implementation, the final Validation population contains **927 cases**.
It consists of 458 zero-K cases, 28 nonzero Type I/III cases, 412 nonzero Type IV cases, 19 two-K cases, and 10 cases with three or more K points.

| Item | Count |
|---|---:|
| Population / inputs / CIFs | 927 / 927 / 927 |
| Validated and Judged | 927 |
| Strict pass | 783 |
| Physical-only pass | 144 |
| Physical fail | 0 |
| Unjudged | 0 |

- **Physical pass: 927 / 927 = 100.00%**
- **Strict pass: 783 / 927 = 84.47%**

#### By K signature

As in the MP pool, entries are aggregated by the number of F and P slots without distinguishing slot order.

| K composition (order-independent) | Population | Physical pass | Strict pass |
|---|---:|---:|---:|
| F | 798 | 798 (100.00%) | 679 (85.09%) |
| FF | 84 | 84 (100.00%) | 69 (82.14%) |
| P | 40 | 40 (100.00%) | 30 (75.00%) |
| PP | 3 | 3 (100.00%) | 3 (100.00%) |
| FP | 1 | 1 (100.00%) | 1 (100.00%) |
| FPP | 1 | 1 (100.00%) | 1 (100.00%) |

The result is clearly better than for the stress population because, even for the same K signature, the stress population intentionally contains many difficult and complex cases.

## Web Comparison and Mathematical Checks

Agreement with the Web implementation is an important compatibility metric, but the Web implementation is not treated as the sole source of truth.
Distortropy checks both Web and Local outputs against necessary conditions, including the exact rational rank of the displayed mode basis, without assuming that either output is correct.

External blind A/B audits confirmed that both Web and Local outputs include cases that do not satisfy necessary mathematical conditions on the displayed subgroup or mode basis.
The remaining long tail will therefore be assessed along two separate axes.

1. **Web compatibility**: Whether the output agrees with the Web implementation at the Strict or Physical level.
2. **Independent mathematics**: Whether it satisfies group-action invariance, basis independence, and related conditions.

## Reproducibility and Notes

The completed authorities for this snapshot are as follows.

- MP pool: runtime `6b9aa754cf70`, comparator `ae1ad1486dbb`
- MAGNDATA: runtime `000ce2ae1f98`, comparator `validation.v7+e276b8fbd237`

Computation time limits, Web acquisition failures, Local timeouts, and malformed inputs are recorded as operational failures separately from scientific Physical/Strict verdicts.
Particularly for heavy multi-K and parametric cases, there is a long tail in computation time as well as in output correctness.
