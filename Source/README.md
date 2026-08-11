# Source

Bundled upstream ISOTROPY-suite material.

- `data_*.txt` — the authoritative runtime **source of truth (SoT)**.  The local
  reconstruction may read these at runtime; nothing else counts as SoT.
- `iso`, `findsym` — the reference **binaries**.  They are oracle *producers*
  (they generate comparison oracle by reading `data_*`), **not** SoT and **not**
  a local-runtime input.
- `Decompiled/`, decompilation / call-graph material — reverse-engineering
  reference only.

Web ISODISTORT output is never SoT and never a local-runtime input (see the
root `README.md` "Web 非注入" policy).
