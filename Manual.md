# Manual

## Build

```sh
sh build.sh
CC=x86_64-w64-mingw32-gcc sh build.sh
```

Executable: `bin/apo` (`bin/apo.exe` with MinGW-w64). Data: adjacent `bin/source.blob`.

## Commands

| Invocation | Output |
|---|---|
| `apo info CIF` | Parent cell and selectable sites |
| `apo kpoints CIF` | k-point labels, vectors, stars, parameter names |
| `apo irreps CIF SELECTION` | Available irreps |
| `apo opds CIF SELECTION` | OPD labels, subgroups, bases, origins, domains |
| `apo modes CIF SELECTION --opd LABEL` | Child structure and mode definitions |
| `apo run --case FILE --upto STAGE` | Calculation through the selected stage |
| `apo serve [--host HOST] [--port PORT] [--open-browser]` | GUI server |
| `apo show --case FILE [--server URL] [--amplitude ID=VALUE ...]` | GUI launch revision |

`SELECTION` consists of the site, strain, k-point, and irrep arguments below.
`STAGE`: `kpoints`, `irreps`, `opds`, `modes`; default: `modes`.

## Selection arguments

| Argument | Value |
|---|---|
| `--displacive SITE [SITE ...]` | Displacement sites |
| `--magnetic SITE [SITE ...]` | Magnetic sites |
| `--occupancy SITE [SITE ...]` | Occupancy sites |
| `--strain` | Homogeneous strain |
| `--k LABEL [PARAMS ...]` | k-point label and positional or named parameters |
| `--irrep LABEL` | Irrep label |
| `--opd LABEL` | OPD label |
| `--sg N` | Space group 1–230 in its default setting; replaces CIF for `kpoints`, `irreps`, `opds` |

| Selection | Definition |
|---|---|
| `SITE` | Site label or element type from `info`; a type selects all matching sites |
| k label / parameter names | Returned by `kpoints` |
| k parameters | Rational or decimal values; `1/3` or `b=1/3` |
| k slots | 1–4 ordered `--k` / `--irrep` pairs |
| Irrep label | Returned by `irreps` |
| OPD label | Returned by `opds`, including domain suffixes |
| Strain-only `irreps` | Default k: `GM` |
| `--sg` site | General-position site `A`; default mode: displacive |

```sh
apo info structure.cif
apo kpoints structure.cif
apo irreps structure.cif --k R --displacive O
apo opds structure.cif --k R --irrep R4- --displacive O
apo modes structure.cif --k R --irrep R4- --displacive O --opd P1
apo irreps structure.cif --k DT b=1/3 --displacive O
apo modes structure.cif --k R --irrep R4- --k M --irrep M3+ --displacive O --opd 'P1(1)P1(2)'
apo opds --sg 62 --k GM --irrep mGM1+ --magnetic A
```

## Case JSON

```json
{
  "structure": "structure.cif",
  "sites": {"O": {"displacive": true, "magnetic": false, "occupancy": false}},
  "strain": false,
  "k": [{"label": "R", "params": {}, "ir": "R4-"}],
  "opd": "P1"
}
```

| Field | Type / meaning |
|---|---|
| `structure` | CIF path; relative to the case file's directory |
| `cif` | Alternative to `structure`: 12-character ID in `assets/cif/` or `assets/magndata/cif/` |
| `sites` | Element → Boolean `displacive`, `magnetic`, optional `occupancy` |
| `strain` | Boolean; default `false` |
| `k` | Up to four ordered objects: `label`, optional `params`, `ir` |
| `opd` | OPD label |

Exactly one of `structure` and `cif`. Accepted fields: those listed above.
`ir` is required from `opds`; `opd` is required for `modes`.
Case filename extension: `.json`. `--case -`: stdin.
Case input and command-line selection arguments are mutually exclusive.

```sh
apo run --case case.json --upto opds -o opds.json
apo run --case case.json --upto modes -o modes.json
apo modes --case case.json --format text
```

## Output arguments

| Argument | Effect |
|---|---|
| `-o PATH`, `--output PATH` | Output file; default: stdout |
| `--indent N` | JSON indentation |
| `--full-state` | Complete calculation state |
| `--format json`, `--format text` | Mode output format; default: JSON |
| `--certify` | Exact certification; counters on stderr |

`--format` applies to `modes` and `run --upto modes`.
Errors: stderr, exit status 2.
Mode definition fields: `definition_id`, `role`, `slot`, `invariants`, `setting`.
`definition_id` is scoped to its result. `setting` contains `normfactor` and `vectors`.

## GUI

| Parameter | Default / meaning |
|---|---|
| `serve --host` | `127.0.0.1` |
| `serve --port` | `8300` |
| `serve --open-browser` | Open the input page |
| `show --server` | `http://127.0.0.1:8300` |
| `show --amplitude ID=VALUE` | Initial mode slider coefficient; repeatable |

| Case contents | Opened stage |
|---|---|
| k points | Irrep selection |
| k points + irreps | OPD selection |
| k points + irreps + OPD | Mode output |

Delivery: FIFO through the open input page. Pending capacity: 128 requests or 32 MiB.
Each mode result adds and selects an output tab. Tab state persists within the input-page session.

```sh
apo serve --open-browser
apo show --case case.json
apo show --case case.json --amplitude displacive-1=0.3 --amplitude displacive-2=-0.2
```
