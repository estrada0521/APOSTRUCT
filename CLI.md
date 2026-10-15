# APOSTRUCT CLI Guide

[日本語](CLI_ja.md) | [README](README.md)

The command-line interface exposes the same calculation services as the local graphical interface.
It is designed as a stateless pipeline: inspect one JSON result, choose identifiers from it, and pass those identifiers to the next command.

```text
parent -> k points -> irreps -> OPDs -> modes
```

Run `apo --version` to report the installed version and, in a Git checkout, its commit.
Run `apo COMMAND --help` for the complete option list.
This guide describes the scientific and machine-facing contracts shared by the commands.

## Commands

| Command | Result |
|---|---|
| `info` | normalized parent cell and selectable crystallographic sites |
| `kpoints` | k-point labels, expressions, stars, and parameter names |
| `irreps` | active ordinary or magnetic irreps for selected mode kinds |
| `opds` | OPD labels, subgroups, bases, origins, domains, and ferroic metadata |
| `modes` | atomic/strain mode definitions for one selected OPD |
| `run` | execute a saved JSON case to a requested stage |
| `serve` | start the local graphical interface |
| `show` | load a saved JSON case into a running graphical interface |

## Parent Input

Pass a CIF path as the positional argument:

```bash
apo info structure.cif
apo kpoints structure.cif
```

The parent setting, sites, occupancies, and coordinates come from that file.
Ordinary decimal, estimated-standard-deviation, and rational values such as `1/2` are accepted for numeric coordinates.
If a site cannot be assigned to a Wyckoff position, `info.sites[].wyckoff_mapping_error` retains the reason.
Selecting that site for a displacive or magnetic calculation is an error that names its label, coordinates, and reason.

## Site and Mode Selection

Read the available selectors from `info.sites`.
A unique `label` selects one crystallographic site.
A `type` selects every site carrying that type:

```bash
apo info structure.cif | jq '.sites[] | {label, type, wyckoff}'
apo irreps structure.cif --k GM --displacive O --magnetic Fe
```

`--strain` adds homogeneous strain.
Strain-only `irreps` infers `GM`; other atomic mode selections require explicit sites and k points.

## K Points and Parameters

`kpoints` is the authority for selectable labels and parameter names:

```bash
apo kpoints structure.cif | \
  jq '.kpoints[] | {label, kvector, parameter_names, miller_love_kvector}'
```

CLI input follows `kvector` and `parameter_names`, matching the graphical interface.
Values can be positional in that order or named:

```bash
apo irreps structure.cif --k DT 1/3 --displacive O
apo irreps structure.cif --k DT b=1/3 --displacive O
```

`miller_love_kvector` and resolved `miller_love_parameters` show the corresponding Miller-Love form.
CLI input always follows `kvector` and `parameter_names`.
Parametric values that collapse onto a special k point are rejected; select the corresponding fixed-point label instead.

One to four k/irrep factors can be coupled.
Repeated `--k` and `--irrep` options are paired by order:

```bash
apo opds structure.cif \
  --k R --irrep R4- \
  --k M --irrep M3+ \
  --displacive O
```

## OPDs and Domains

`opds` returns the exact labels accepted downstream.
Do not reconstruct a label from the displayed direction:

```bash
apo opds structure.cif \
  --k R --irrep R4- --displacive O \
  -o opds.json

jq '.opds[] | {label, opd, subgroup, basis, origin, index, cell_index}' \
  opds.json
```

For coupled selections, labels include one domain number per factor, for example `P1(1)P1(2)`.
Shell-quote labels containing parentheses:

```bash
apo modes structure.cif \
  --k R --irrep R4- \
  --k M --irrep M3+ \
  --displacive O --opd 'P1(1)P1(2)'
```

`index` is the group index and `cell_index` is the supercell factor.
A classified empty `ferroic_properties` list is distinguished from an unclassified result by `ferroic_classified`.

## Mode Details

`modes` emits compact JSON by default:

```bash
apo modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  -o modes.json
```

Atomic definitions carry payload-local `definition_id` values and structured mode identity: kind, k point, k vector, irrep, gid, direction, site, Wyckoff position, and site irrep.
`role` distinguishes primary from induced secondary definitions.

Use human-readable mode tables when needed:

```bash
apo modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  --format text
```

## Saved JSON Cases

A case preserves the ordered scientific selection.
It contains exactly one parent source: `structure` or a repository-local `cif` content ID.
The `cif` form resolves a 12-character content ID to `assets/cif/` or `assets/magndata/cif/` in this checkout (the CIFs of the published branches); otherwise use `structure`.

```json
{
  "structure": "structure.cif",
  "sites": {
    "O": {"displacive": true, "magnetic": false}
  },
  "strain": false,
  "k": [
    {"label": "R", "ir": "R4-"}
  ],
  "opd": "P1"
}
```

Run it to one stage:

```bash
apo run --case case.json --upto kpoints
apo run --case case.json --upto irreps
apo run --case case.json --upto opds
apo run --case case.json --upto modes
```

For saved or automated work, prefer `run` so that one process advances the pipeline without reloading the bundled tables for separate stage commands.
Use the individual commands when discovering or inspecting the identifiers to put in a case.

Use `--case -` for standard input.
Relative structure paths resolve from the case file's directory.
A k item uses `label`, optional exact `params`, and `ir` from the OPD stage onward.
JSON is the only saved-case format.

## JSON and Output Control

Compact results identify their shape with an unversioned schema name.
Consumers should branch on this field rather than infer the stage from optional keys:

```text
APOSTRUCT.cli.info
APOSTRUCT.cli.kpoints
APOSTRUCT.cli.irreps
APOSTRUCT.cli.opds
APOSTRUCT.cli.modes
```

Common controls:

```bash
--indent N       JSON indentation
-o, --output     write to a path instead of stdout
--full-state     complete calculation state, where supported
```

Errors are written to stderr and exit with status 2. Unknown k, irrep, and OPD selections fail rather than falling back to a different result.

## Local Interface

```bash
apo serve --host 127.0.0.1 --port 8300 --open-browser
```

The GUI provides the CIF workflow over the same calculations.

Send a case to an interface that is already open:

```bash
apo show --case case.json --server http://127.0.0.1:8300
```

The deepest selection in the case determines the displayed stage.
A case with K but no irrep opens the irrep step; selected irreps without `opd` open the OPD step; `opd` opens mode details and the structure viewer.
Initial mode-slider amplitudes use the payload-local definition order:

```bash
apo show --case case.json \
  --amplitude magnetic-1=1 --amplitude magnetic-2=19/20
```

Relative structure paths are resolved by the sending command, which transfers the CIF contents to the server.
The case and calculated state are held only in the running server's memory.
