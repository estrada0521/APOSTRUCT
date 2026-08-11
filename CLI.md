# APOSTRUCT CLI Guide

[日本語](CLI_ja.md) | [README](README.md)

The command-line interface exposes the same calculation services as the local
graphical interface. It is designed as a stateless pipeline: inspect one JSON
result, choose identifiers from it, and pass those identifiers to the next
command.

```text
parent -> k points -> irreps -> OPDs -> modes
                                  `----> invariants
```

Run `apo --version` to report the installed version and, in a Git
checkout, its commit. Run `apo COMMAND --help` for the complete option list. This guide
describes the scientific and machine-facing contracts shared by the commands.

## Commands

| Command | Result |
|---|---|
| `settings` | Source International settings for one space group |
| `directions` | ordinary or time-odd OP directions compatible with an exact G:H embedding |
| `info` | normalized parent cell and selectable crystallographic sites |
| `kpoints` | k-point labels, expressions, stars, and parameter names |
| `irreps` | active ordinary or magnetic irreps for selected mode kinds |
| `opds` | OPD labels, subgroups, bases, origins, domains, and ferroic metadata |
| `modes` | atomic/strain mode definitions for one selected OPD |
| `invariants` | invariant basis on one selected OPD domain |
| `combine-modes` | per-atom vectors from weighted compact mode definitions |
| `run` | execute a saved JSON case to a requested stage |
| `serve` | start the local graphical interface |
| `show` | load a saved JSON case into a running graphical interface |

## Subgroup-Compatible Directions

`directions` answers the structure-independent inverse question: which
ordinary or time-odd order-parameter directions are compatible with one exact
parent-to-subgroup embedding?

List the accepted International setting IDs first when either group is not in
its default setting:

```bash
apo settings --sg 14
```

```bash
apo directions \
  --parent-sg 14 --parent-setting 64 \
  --subgroup-sg 2 --subgroup-setting 2 \
  --basis="0,0,3;0,1,0;-1,0,0"
```

Omit the setting options to use the defaults:

```bash
apo directions \
  --parent-sg 221 --subgroup-sg 140 \
  --basis="-1,1,0;-1,-1,0;0,0,2" \
  --origin=0,0,0
```

For a magnetic child, give its BNS number instead of an ordinary subgroup
number. The parent is the paramagnetic gray extension of `--parent-sg`:

```bash
apo directions \
  --parent-sg 62 --subgroup-msg 62.448 \
  --basis="1,0,0;0,1,0;0,0,1"
```

The basis vectors and origin use the selected International settings of the
parent and subgroup. COPL uses the same IDs for the subset in its setting menu;
APOSTRUCT also exposes non-legacy axis relabelings stored by Source. Space-group
numbers alone are insufficient because the same subgroup type can have
inequivalent orientations, cells, origins, and domains inside a parent. No CIF,
Wyckoff site, or distortion-mode selection is used.

Each result gives the k point, irrep, OPD label, structured direction matrix,
stabilizer subgroup, group index, and cell index. Direction-matrix rows follow
the full irrep coordinates and columns follow `parameters`. Parametric rows
retain both public and Miller-Love parameter values. `opd` is null when no OPD
label is available for the returned direction. `primary` means that the
direction alone has exactly the requested subgroup embedding; more than one
returned direction may meet that condition. Its non-null `domain` can be reused
directly by `modes --from-directions`. `secondary` means its stabilizer is a
supergroup but the direction is allowed by the requested subgroup.

## Parent Inputs

### CIF

Pass a CIF path as the positional argument:

```bash
apo info structure.cif
apo kpoints structure.cif
```

The parent setting, sites, occupancies, and coordinates come from that file.
Ordinary decimal, estimated-standard-deviation, and rational values such as
`1/2` are accepted for numeric coordinates. If a site cannot be assigned to a
Wyckoff position, `info.sites[].wyckoff_mapping_error` retains the reason.
Selecting that site for a displacive or magnetic calculation is an error that
names its label, coordinates, and reason.

### Space group and Wyckoff sites

Use `--sg` when a concrete structure is unnecessary:

```bash
apo kpoints --sg 221
apo info --sg 221 --wyckoff 1a 1b 3c
```

The symbolic route uses the default International Tables setting. Wyckoff
letters and multiplicity-letter forms are both accepted. Input free coordinates
may remain absent through `invariants`:

```bash
apo opds \
  --sg 137 --wyckoff 2a 4d \
  --k M --irrep M1 --displacive a d
```

`info` marks such sites as unrealized and may display generic placeholder
coordinates for identification. Those values are not a physical structure and
are not used as geometry by the k-point, irrep, OPD, or invariant calculation.
`modes` realizes atomic geometry and therefore requires every free coordinate:

```bash
apo modes \
  --sg 137 --wyckoff 2a 4d:z=1/7 \
  --k M --irrep M1 --displacive a d --opd P1
```

The supplied values define a geometry; they are not generic placeholders.

## Site and Mode Selection

Read the available selectors from `info.sites`. A unique `label` selects one
crystallographic site. A `type` selects every site carrying that type:

```bash
apo info structure.cif | jq '.sites[] | {label, type, wyckoff}'
apo irreps structure.cif --k GM --displacive O --magnetic Fe
```

For repeated symbolic Wyckoff orbits, labels such as `i1` and `i2` select them
independently, while type `i` selects both. `--strain` adds homogeneous strain.
Site-free, strain-only `irreps --sg SG --strain` infers `GM`. Site-free
strain-only `modes` also requires no Wyckoff sites; other atomic mode selections
require explicit sites and k points.

## K Points and Parameters

`kpoints` is the authority for selectable labels and parameter names:

```bash
apo kpoints --sg 225 | \
  jq '.kpoints[] | {label, kvector, parameter_names, miller_love_kvector}'
```

CLI input follows `kvector` and `parameter_names`, matching the graphical
interface. Values can be positional in that order or named:

```bash
apo irreps structure.cif --k DT 1/3 --displacive O
apo irreps structure.cif --k DT b=1/3 --displacive O
```

`miller_love_kvector` and resolved `miller_love_parameters` show the
corresponding Miller-Love form. CLI input always follows `kvector` and
`parameter_names`. Parametric values that collapse onto a special k point are
rejected; select the corresponding fixed-point label instead.

One to four k/irrep factors can be coupled. Repeated `--k` and `--irrep` options
are paired by order:

```bash
apo opds structure.cif \
  --k R --irrep R4- \
  --k M --irrep M3+ \
  --displacive O
```

## OPDs and Domains

`opds` returns the exact labels accepted downstream. Do not reconstruct a
label from the displayed direction:

```bash
apo opds structure.cif \
  --k R --irrep R4- --displacive O \
  -o opds.json

jq '.opds[] | {label, opd, subgroup, basis, origin, index, cell_index}' \
  opds.json
```

For coupled selections, labels include one domain number per factor, for
example `P1(1)P1(2)`. Shell-quote labels containing parentheses:

```bash
apo modes structure.cif \
  --k R --irrep R4- \
  --k M --irrep M3+ \
  --displacive O --opd 'P1(1)P1(2)'
```

`index` is the group index and `cell_index` is the supercell factor. A
classified empty `ferroic_properties` list is distinguished from an
unclassified result by `ferroic_classified`.

## Mode Details

`modes` emits compact JSON by default:

```bash
apo modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  -o modes.json
```

When a reported subgroup embedding selects a symmetry-equivalent domain that
differs from the direct OPD result, reuse the exact `directions` result instead
of reducing it back to the OPD label:

```bash
apo directions \
  --parent-sg 225 --subgroup-sg 87 \
  --basis="1/2,-1/2,0;1/2,1/2,0;0,0,1" \
  -o directions.json

apo modes structure.cif \
  --from-directions directions.json --direction-row 3 \
  --displacive Br -o modes.json
```

`direction-row` is the 1-based `directions[].row` value and must select a
primary row. The calculation preserves the saved basis, origin, setting,
direction subspace, and domain. Direct `--k`, `--irrep`, and `--opd` arguments
are then unnecessary and unavailable.

Atomic definitions carry payload-local `definition_id` values and structured
mode identity: kind, k point, k vector, irrep, gid, direction, site, Wyckoff
position, and site irrep. `role` distinguishes primary from induced secondary
definitions. An unambiguous primary also carries its invariant factor slot and
global parameter names.

Use human-readable mode tables when needed:

```bash
apo modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  --format text
```

### Combining definitions

Inspect the definition IDs in one compact modes payload:

```bash
jq '.mode_details.displacive_definitions[],
    .mode_details.magnetic_definitions[] |
    {definition_id, normfactor, role, mode}' modes.json
```

`--weight` multiplies the published rows directly:

```bash
apo combine-modes modes.json \
  --weight magnetic-1=1 --weight magnetic-2=3/4
```

`--amplitude` first applies each definition's published mode normfactor:

```bash
apo combine-modes modes.json \
  --amplitude magnetic-1=1 --amplitude magnetic-2=1
```

These amplitudes follow each definition's mode normalization and do not
establish a common physical unit across unrelated definitions. Returned vector
components use the child crystallographic `dxyz` basis. The magnetic net is
the sum over the returned conventional child-cell atoms.

## Invariants

Direct `invariants` calculations use the selected primary factors and exact
OPD domains:

```bash
apo invariants structure.cif \
  --k GM --irrep GM4- \
  --k GM --irrep GM3+ \
  --displacive O --strain \
  --opd 'P1(1)P1(3)' \
  --minimum-degree 1 --maximum-degree 4
```

Degrees 1 through 12 are supported. Every requested degree is emitted; a
degree with no invariants has `count: 0`, `invariants: []`, and
`polynomials: []`. Display strings are retained for reading, while structured
polynomials are intended for machine use:

```json
{
  "terms": [
    {"coefficient": "1", "exponents": [3, 0]},
    {"coefficient": "-3", "exponents": [1, 2]}
  ]
}
```

Exponent positions follow the top-level `variables` array. Coefficients are
exact SymPy expression strings, not floating-point approximations.

### Including induced secondary factors

`modes` records ordered primary and secondary `invariant_factors`, including
their exact domains and global parameter offsets. Reuse that payload without
rerunning the parent, OPD, or mode calculation:

```bash
jq '.invariant_factors[] |
    {slot, role, label, opd, domain, parameter_offset}' modes.json

apo invariants \
  --from-modes modes.json \
  --secondary 4 \
  --minimum-degree 2 --maximum-degree 4
```

`--secondary` is repeatable and takes payload-local factor slots. The primary
factors are always included. There is no manual secondary-domain override; the
domain carried by the modes result remains authoritative.

## Saved JSON Cases

A case preserves the ordered scientific selection. It contains exactly one
parent source: `structure`, a repository-local `cif` content ID, or `sg` with
optional `wyckoff`. The `cif` form resolves an existing `Assets/cif` content ID
in a development checkout; standalone users normally use `structure`.

```json
{
  "sg": 221,
  "wyckoff": ["1a", "1b", "3c"],
  "sites": {
    "c": {"displacive": true, "magnetic": false}
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
apo run --case case.json --upto invariants \
  --minimum-degree 2 --maximum-degree 6
```

For saved or automated work, prefer `run` so that one process advances the
pipeline without reloading the bundled tables for separate stage commands.
Use the individual commands when discovering or inspecting the identifiers to
put in a case.

Use `--case -` for standard input. Relative structure paths resolve from the
case file's directory. A k item uses `label`, optional exact `params`, and `ir`
from the OPD stage onward. JSON is the only saved-case format.

## JSON and Output Control

Compact results identify their shape with an unversioned schema name. Consumers
should branch on this field rather than infer the stage from optional keys:

```text
APOSTRUCT.cli.settings
APOSTRUCT.cli.directions
APOSTRUCT.cli.info
APOSTRUCT.cli.kpoints
APOSTRUCT.cli.irreps
APOSTRUCT.cli.opds
APOSTRUCT.cli.modes
APOSTRUCT.cli.invariants
APOSTRUCT.cli.mode_combination
```

Common controls:

```bash
--indent N       JSON indentation
-o, --output     write to a path instead of stdout
--full-state     complete calculation state, where supported
```

`invariants` and `combine-modes` are compact-only. `run --upto invariants`
therefore rejects `--full-state`. Errors are written to stderr and exit with
status 2. Unknown k, irrep, and OPD selections fail rather than falling back to
a different result.

## Local Interface

```bash
apo serve --host 127.0.0.1 --port 8300 --open-browser
```

The GUI provides CIF and symbolic-parent workflows over the same calculations.
Debug display changes the amount of information shown, not the selected
calculation.

Send a case to an interface that is already open:

```bash
apo show --case case.json --server http://127.0.0.1:8300
```

The deepest selection in the case determines the displayed stage. A case with
K but no irrep opens the irrep step; selected irreps without `opd` open the OPD
step; `opd` opens mode details and the structure viewer. Initial mode-slider
amplitudes use the payload-local definition order:

```bash
apo show --case case.json \
  --amplitude magnetic-1=1 --amplitude magnetic-2=19/20
```

Relative structure paths are resolved by the sending command, which transfers
the CIF contents to the server. The case and calculated state are held only in
the running server's memory.
