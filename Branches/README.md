# Published branches

`Branches/` is a generated data bundle, not a copy of the private Validation operation.
It contains no saved Web output, Local output, sweep definition, worker plan, dashboard state, or Web automation.

## Files

- `manifest.json`: bundle schema, counts, comparator versions, and SHA-256 digests of the two compressed JSONL files.
- `inputs.jsonl.gz`: one canonical input payload per content-addressed `input_id`.
- `results.jsonl.gz`: one landed result per `(input_id, opd)` branch.

An input row has this shape:

```json
{"input_id": "...", "input": {"cif": "...", "k": [], "sites": {}, "strain": false}}
```

A result row records the input and CIF identity, OPD, pool, K signature, content IDs of the observed Web and Local outputs, Strict/Physical verdict, comparator version and tolerance, Local runtime and landed commit, and the available mathematics status and verifier version for each output.
The output IDs are provenance; their payloads are deliberately not distributed.
MAGNDATA rows also record the public `magndata_id` and the database from which the ordinary parent CIF was obtained.
For Crystallography Open Database parents, `parent_cif_source_id` is the public COD identifier; Materials Project parent IDs are omitted because the available mapping contains internal synthetic identifiers rather than public Materials Project IDs.

The generator starts from the canonical successful Web observations and accepts only a completed Validation set whose exact Web and Local output IDs match a dirty=false compute from a runtime observed on main.
Candidate and experimental sets are excluded.
For each `(input_id, opd)` it selects the newest such landed runtime; append order resolves repeated judgements within that runtime.
Each branch retains the comparator and mathematics-verifier versions that actually produced its recorded result.
A refactor of the shipped tools does not relabel historical observations or force a corpus-wide rejudgement.

To inspect a branch, extract its canonical input, use the referenced CIF under `Assets/`, run APOSTRUCT locally, and submit the same request and OPD to the public ISODISTORT Web application yourself.
See `Verification/README.md` for the single-branch comparison and mathematics commands.
