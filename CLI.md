# CLI

[日本語](CLI_ja.md) · [README](README.md)

Each command is stateless. Pass a label returned in one stage's JSON to the next stage.

```text
info -> kpoints -> irreps -> opds -> modes
```

## Commands

| Command | Result |
| --- | --- |
| `info` | Normalized parent cell and selectable sites |
| `kpoints` | k-point labels, kvectors, stars and parameter names |
| `irreps` | Irreps active for the selected sites and strain |
| `opds` | OPD list, with the subgroup, basis, origin and domains of each OPD |
| `modes` | Mode details for one OPD |
| `run` | Runs a JSON case up to a given stage |
| `serve` | Starts the browser interface |
| `show` | Loads a case into a running interface |

## Parent

Pass the CIF path as the positional argument. Setting, sites, occupancies and coordinates follow the CIF.
Coordinates may be decimals, values with uncertainties (`0.1234(5)`) or rationals such as `1/2`.
A site whose Wyckoff position cannot be determined carries the reason in `info.sites[].wyckoff_mapping_error`; selecting it stops the calculation.

Without a structure, `kpoints`, `irreps` and `opds` take `--sg N` (1 to 230) instead of the CIF: the parent is space group N in its default setting, and k points, irreps and OPDs depend only on the group.
The parent then has one site `A` at a general position, selected as displacive unless `--displacive A` or `--magnetic A` is given; its mode totals in the output refer to that site.

```sh
apo opds --sg 221 --k R --irrep R5-
apo opds --sg 62 --k GM --irrep mGM1+ --magnetic A
```

## Sites

`--displacive`, `--magnetic` and `--occupancy` take a `label` or `type` from `info.sites`. A `label` selects one site, a `type` all sites of that type.
`--strain` adds homogeneous strain. For strain-only `irreps`, `--k` may be omitted and is taken as `GM`.

```sh
apo irreps structure.cif --k GM --displacive O --magnetic Fe
```

## k points

Labels and parameter names follow `label` and `parameter_names` in `kpoints`. Give parameters in order or by name.

```sh
apo irreps structure.cif --k DT 1/3 --displacive O
apo irreps structure.cif --k DT b=1/3 --displacive O
```

`miller_love_kvector` is for display in the output and is not accepted as input.
Parameter values at which a parametric k coincides with a special point are rejected; use the fixed-point label instead.

`--k` and `--irrep` can be given one to four times and pair in the order written.

```sh
apo opds structure.cif --k R --irrep R4- --k M --irrep M3+ --displacive O
```

## OPDs

Pass a `label` returned by `opds` to `--opd` as is. Do not rebuild it from the direction.
With several k points the label carries a domain number per slot (`P1(1)P1(2)`); quote it in the shell.

```sh
apo modes structure.cif --k R --irrep R4- --k M --irrep M3+ --displacive O --opd 'P1(1)P1(2)'
```

Each OPD row separates what does not depend on choices from what does. `invariants` holds the facts of the isotropy subgroup: `parameter_count` (free parameters, the dimension of its fixed space), `index` (the group index, also the number of domains), `cell_index` (the supercell factor), `k_active_count` (the number of active star arms) and `ferroic_properties`. `setting` holds one presentation of it: the direction `opd` in the BYU irrep basis, the subgroup `basis` and `origin`, and the active arms `k_active`. These change with the choice of basis, domain or arms; `subgroup` and the label prefix (P, C, S, nD) do not.

## Modes

Each mode has a `definition_id` (an ID within that payload), a `role` (primary or induced secondary) and the `slot` of its k point. `invariants` identifies it apart from any choice: kind, k point, irrep, site, Wyckoff position and site irrep. `setting` holds one presentation: the label, the direction in the BYU irrep basis, the active arm `kvector`, the parameter names, `normfactor` and the `vectors`: a displacement or moment per atom, or for an occupancy mode one occupancy value per atom, whose signs and, in a space with several modes, whose basis follow the choices of basis and domain. `subgroup_details` separates the `subgroup` from the `setting` (basis, origin, parent setting) of the child cell.
`--format text` prints the mode table (`modes` and `run --upto modes` only).

## Case

A selection collected in JSON.

```json
{
  "structure": "structure.cif",
  "sites": {"O": {"displacive": true, "magnetic": false}},
  "strain": false,
  "k": [{"label": "R", "ir": "R4-"}],
  "opd": "P1"
}
```

- Fields are limited to `structure`, `cif`, `sites`, `strain`, `k` and `opd`. Exactly one of `structure` and `cif`.
- Each site may also set `"occupancy": true`.
- A relative `structure` path is resolved from the case file's directory.
- `cif` is a 12-character content ID looked up in `assets/cif/` or `assets/magndata/cif/` of this checkout.
- Each `k` item has `label`, `params` (an object of names and values) and `ir`. `ir` is required from the OPD stage on.
- The file needs the `.json` extension. `--case -` reads from standard input.

```sh
apo run --case case.json --upto opds
apo modes --case case.json
```

`run` computes all stages in one process; use it for automated work. `--case` cannot be combined with the individual selection options.

## Output

Output is JSON, and `schema` (`APOSTRUCT.cli.info` and so on) names the stage.

| Option | Effect |
| --- | --- |
| `-o PATH` | Writes to a file instead of standard output |
| `--indent N` | JSON indentation |
| `--full-state` | Prints the complete calculation state |

Errors go to standard error with exit status 2. An unknown k, irrep or OPD is never replaced by another one.

## GUI

```sh
apo serve --host 127.0.0.1 --port 8300 --open-browser
apo show --case case.json --server http://127.0.0.1:8300
```

`show` opens the deepest stage in the case: irrep selection for k only, OPD selection up to irreps, and mode details with the structure view when `opd` is present.
`--amplitude DEFINITION_ID=VALUE` sets an initial mode slider value (`--amplitude magnetic-1=1`).
`show` reads the CIF itself and sends its text to the server. Each mode-details result adds an output tab without replacing existing tabs or changing the tab being viewed. Each tab retains its mode amplitudes, camera, display options and scroll positions until closed. Results live in the current GUI session, without disk persistence; the input window can restore a closed/reloaded output window while it remains open.

Consecutive `show` requests are delivered in order. Keep the input page open so it can receive them. If its pending queue is full (128 requests or 32 MiB), a new request reports an error instead of silently replacing an earlier result; wait for the input page to receive pending requests before retrying.
The Examples band above step 1 opens a published branch (a random MP example, or a MAGNDATA example by category) from `branches/examples.json.gz`, shows its input ID and the computation time, and copies the current selection as a case for `--case`.
