# Standalone verification

This directory is the single-branch verification boundary shared by the private monorepo Validation operation and the public distribution.
It does not fetch ISODISTORT Web output, run validation sweeps, or read private content-addressed stores.

The public `branches/` bundle contains the canonical input and OPD for each published result.
To inspect a branch, use its CIF and canonical input to run APOSTRUCT locally, submit the same request to the public ISODISTORT Web application yourself, and save the Web `Complete Mode Details` text.

Run the entry points from the root of this checkout on those explicit files:

```bash
python3 -m verification.verify_local_output local.json
python3 -m verification.verify_web_text complete-mode-details.txt
python3 -m verification.verify_comparison \
  local.json complete-mode-details.txt \
  --input case.json --opd P1
```

`verify_local_output.py` and `verify_web_text.py` apply the two independent necessary-condition checks under `verification/mathematics/`:

- exact rational linear independence of the complete displayed mode basis;
- invariance of displayed mode fields under the reported subgroup operations, allowing only the uncertainty implied by printed rounding cells.

`verify_comparison.py` is the same comparator called by the private Validation harness.
An identical verdict requires the displayed setting, structure, labels, normalizations, ordering, and vectors to agree.
An equivalent verdict admits only the explicit crystallographic and mode-space equivalences implemented under `verification/comparison/`.
An unproven verdict means those checks did not prove equivalence; it is deliberately fail-closed and does not say which side is wrong.
Records written before 2026-10-23 name these strict_pass, physical_pass and physical_fail; `verification.result.canonical_verdict` reads both.

The verifier needs these Python packages:

```bash
python3 -m pip install numpy sympy scipy spglib gemmi
```

To reproduce the published population counts without rerunning a branch:

```bash
python3 -m verification.reproduce
```
