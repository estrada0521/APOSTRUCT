# APOSTRUCT

[日本語](README_ja.md) | [CLI guide](CLI.md)

**APOSTRUCT** is a local Python package for symmetry analysis of crystalline structures.
It computes k vectors, irreducible representations, order-parameter directions, isotropy subgroups, invariant bases, symmetry-adapted modes, and Landau invariants from the bundled BYU ISOTROPY tables.

<p align="center">
  <img src="media/APOSTRUCT.png" alt="APOSTRUCT graphical interface" width="100%">
</p>

On the name — apo- is the Greek prefix for "away from," used here as in apomorphy: a state derived away from its ancestor.
Where the iso- of ISODISTORT names the subgroup that leaves the order parameter invariant — the symmetry that remains — apo- names the departure from the parent.
struct is what departs: the structure itself, displacive, magnetic, and strain alike.

## Interfaces And Scope

- Parent input: a CIF, or a space group with ordered Wyckoff sites.
- Distortion input: displacive sites, magnetic sites, homogeneous strain, and
  one to four ordered k/irrep factors.
- Forward pipeline: parent information, k points, active irreps, OPDs,
  subgroups, invariant bases, and atomic or strain mode definitions.
- Embedding query: ordinary or time-odd OP directions compatible with a given
  parent/subgroup basis and origin.
- Output: compact JSON, optional full pipeline state, text mode tables, and
  reusable saved JSON cases.
- Interfaces: the `apo` command and a local browser interface.

The CLI and graphical interface are entry points to the same backend services.
They differ in workflow and presentation, not in the scientific calculation
path.

## Install

Clone the repository and install the package:

```bash
python -m pip install .
```

For development, use `python -m pip install -e .` instead. Both forms install
the `apo` command. The required group-theory tables are included in the
package.

## Quick Start

For a saved selection, run the pipeline to the required stage in one process.
For example, for a SrTiO3 parent written with Sr at 1a, Ti at 1b, and O at
3c, `case.json` may contain:

```json
{
  "structure": "structure.cif",
  "sites": {"O": {"displacive": true, "magnetic": false}},
  "strain": false,
  "k": [{"label": "R", "ir": "R5-"}],
  "opd": "P1"
}
```

```bash
apo run --case case.json --upto opds
apo run --case case.json --upto modes
apo run --case case.json --upto invariants \
  --minimum-degree 2 --maximum-degree 6
```

`run` avoids starting a new process and loading the bundled tables separately
for every pipeline stage. Use the individual commands when inspecting the
available identifiers or choosing a new case:

```bash
# Inspect the parent and its selectable crystallographic sites.
apo info structure.cif

# Discover reciprocal-space and representation choices.
apo kpoints structure.cif
apo irreps structure.cif --k R --displacive O

# Enumerate OPDs and inspect their exact labels.
apo opds structure.cif --k R --irrep R5- --displacive O
```

When no concrete structure is needed, use the space-group/Wyckoff route:

```bash
apo opds \
  --sg 221 --wyckoff 1a 1b 3c \
  --k R --irrep R5- --displacive c
```

This origin choice labels the antiphase octahedral rotation `R5-`. Translating
the parent origin by `(1/2,1/2,1/2)` places Ti at 1a and O at 3d, and labels the
same physical rotation by the commonly cited Howard-Stokes symbol `R4+`.

Free Wyckoff coordinates may remain unspecified through OPD and invariant
calculations. They are required only when atomic mode geometry must be
realized. Site-free strain modes do not require Wyckoff sites.

The [CLI guide](CLI.md) defines the complete command workflow, input
conventions, coupled selections, parametric k points, saved cases, exact
embedding reuse, secondary invariant factors, and machine-readable output.

## Graphical Interface

```bash
apo serve --host 127.0.0.1 --port 8300 --open-browser
```

Load an existing case into that running interface with
`apo show --case case.json`.

The interface supports the ordinary CIF workflow and a direct
space-group/Wyckoff workflow for symbolic calculations.

## Output

Pipeline commands emit compact JSON to standard output by default. `modes`
also supports `--format text`; `--full-state` exposes the complete calculation
state where available. Use `-o PATH` to write a file. No result file
is created unless one is requested.

## Validation

See [Validation methodology and results](Validation.md).

## Notice

The original ISODISTORT software and the group-theory tables are the work of
Harold T. Stokes, Dorian M. Hatch, and Branton J. Campbell at Brigham Young
University. Please acknowledge their work when using APOSTRUCT in research.
See [NOTICE](NOTICE) for details.
