# Standalone verification

This directory is the single-branch verification boundary shared by the
private monorepo Validation operation and the public distribution. It does
not fetch ISODISTORT Web output, run validation sweeps, or read private
content-addressed stores.

The public `Branches/` bundle contains the canonical input and OPD for each
published result. To inspect a branch, use its CIF and canonical input to run
APOSTRUCT locally, submit the same request to the public ISODISTORT Web
application yourself, and save the Web `Complete Mode Details` text.

The standalone entry points then operate on those explicit files:

```bash
apostruct-verify-local local.json
apostruct-verify-web complete-mode-details.txt
apostruct-compare \
  local.json complete-mode-details.txt \
  --input case.json --opd P1
```

`verify_local_output.py` and `verify_web_text.py` apply the two independent
necessary-condition checks under `Verification/mathematics/`:

- exact rational linear independence of the complete displayed mode basis;
- invariance of displayed mode fields under the reported subgroup operations,
  allowing only the uncertainty implied by printed rounding cells.

`verify_comparison.py` is the same comparator called by the private Validation
harness. A Strict pass requires the displayed setting, structure, labels,
normalizations, ordering, and vectors to agree. A Physical pass admits only
the explicit crystallographic and mode-space equivalences implemented under
`Verification/comparison/`. A Physical fail means those checks did not prove
equivalence; it is deliberately fail-closed.

The verifier depends on the public APOSTRUCT and Source modules included in
the same checkout. Install the optional comparison dependency with:

```bash
python3 -m pip install -e '.[verification]'
```

To reproduce the published population counts without rerunning a branch:

```bash
apostruct-validation-summary
```

The public bundle contains no saved Web or Local output payloads and no Web
automation wrapper.
