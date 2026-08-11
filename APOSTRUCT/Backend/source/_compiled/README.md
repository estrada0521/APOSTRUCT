# Compiled Source Tables

These files are deterministic, read-only projections of the distributed
`Source/data_*.txt` tables. The text files remain the authority. The runtime
verifies their hashes from `manifest.json` before loading requested sections.

In a source checkout, regenerate the directory from the repository root with:

```bash
python3 -m Tests.APOSTRUCT.source_storage.build --replace
```

Do not edit the `.npy` files or manifest directly.
