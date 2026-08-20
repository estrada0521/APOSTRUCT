# Published branches

`Branches/` is a generated data bundle, not a copy of the private Validation
operation. It contains no saved Web output, Local output, sweep definition,
worker plan, dashboard state, or Web automation.

## Files

- `manifest.json`: bundle schema, counts, comparator versions, and SHA-256
  digests of the two compressed JSONL files.
- `inputs.jsonl.gz`: one canonical input payload per content-addressed
  `input_id`.
- `results.jsonl.gz`: one landed result per `(input_id, opd)` branch.

An input row has this shape:

```json
{"input_id": "...", "input": {"cif": "...", "k": [], "sites": {}, "strain": false}}
```

A result row records the input and CIF identity, OPD, pool, K signature,
content IDs of the observed Web and Local outputs, Strict/Physical verdict,
comparator version and tolerance, Local runtime and landed commit, and the
available mathematics status and verifier version for each output. The output
IDs are provenance; their payloads are deliberately not distributed.

The generator accepts only a completed Validation set whose exact Local output
matches a dirty=false compute from a runtime observed on main. Candidate and
experimental sets are excluded. Each branch retains the comparator and
mathematics-verifier versions that actually produced its recorded result; a
refactor of the shipped tools does not relabel historical observations or
force a corpus-wide rejudgement.

To inspect a branch, extract its canonical input, use the referenced CIF under
`Assets/`, run APOSTRUCT locally, and submit the same request and OPD to the
public ISODISTORT Web application yourself. See `Verification/README.md` for
the single-branch comparison and mathematics commands.
