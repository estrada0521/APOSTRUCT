# APOSTRUCT

[日本語](README_ja.md) | [CLI guide](CLI.md)

**APOSTRUCT** is a local Python package for symmetry analysis of crystalline structures.
It computes k vectors, irreducible representations, order-parameter directions, isotropy subgroups, and symmetry-adapted modes from the bundled BYU ISOTROPY tables.

<p align="center">
  <img src="media/APOSTRUCT.png" alt="APOSTRUCT graphical interface" width="100%">
</p>

On the name — apo- is the Greek prefix for "away from," used here as in apomorphy: a state derived away from its ancestor.
Where the iso- of ISODISTORT names the subgroup that leaves the order parameter invariant — the symmetry that remains — apo- names the departure from the parent. struct is what departs: the structure itself, displacive, magnetic, and strain alike.

## Interfaces And Scope

- Parent input: a CIF.
- Distortion input: displacive sites, magnetic sites, homogeneous strain, and one to four ordered k/irrep factors.
- Pipeline: parent information, k points, active irreps, OPDs, subgroups, and atomic or strain mode definitions.
- Output: compact JSON, optional full pipeline state, text mode tables, and reusable saved JSON cases.
- Interfaces: the `apo` command and a local browser interface.

The CLI and graphical interface are entry points to the same backend services.
They differ in workflow and presentation, not in the scientific calculation path.

## Install

Clone the repository and build the executable with a C compiler (macOS or Linux; only libm and pthreads are needed):

```bash
sh APOSTRUCT/Backend/native/build.sh
ln -s "$PWD/APOSTRUCT/Backend/native/apo" ~/.local/bin/apo
```

The build writes `APOSTRUCT/Backend/native/apo` and the group-theory table blob next to it.
Link the executable into any directory on `PATH`; it finds its tables and the browser interface through the link.
Running `apo` needs no Python. `python -m pip install ".[verification]"` installs only the verification tools.

## Quick Start

For a saved selection, run the pipeline to the required stage in one process.
For example, for a SrTiO3 parent written with Sr at 1a, Ti at 1b, and O at 3c, `case.json` may contain:

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
```

`run` avoids starting a new process and loading the bundled tables separately for every pipeline stage.
Use the individual commands when inspecting the available identifiers or choosing a new case:

```bash
# Inspect the parent and its selectable crystallographic sites.
apo info structure.cif

# Discover reciprocal-space and representation choices.
apo kpoints structure.cif
apo irreps structure.cif --k R --displacive O

# Enumerate OPDs and inspect their exact labels.
apo opds structure.cif --k R --irrep R5- --displacive O
```

The [CLI guide](CLI.md) defines the complete command workflow, input conventions, coupled selections, parametric k points, saved cases, and machine-readable output.

## Graphical Interface

```bash
apo serve --host 127.0.0.1 --port 8300 --open-browser
```

Load an existing case into that running interface with `apo show --case case.json`.

The interface supports the CIF workflow.

## Output

Pipeline commands emit compact JSON to standard output by default.
`modes` also supports `--format text`; `--full-state` exposes the complete calculation state where available.
Use `-o PATH` to write a file.
No result file is created unless one is requested.

## Validation

See [Validation methodology and results](Validation.md).
The source checkout also contains the standalone, single-branch procedures in [`Verification/`](Verification/README.md), separate from the `apo` command.

## Notice

The original ISODISTORT software and the group-theory tables are the work of Harold T. Stokes, Dorian M. Hatch, and Branton J. Campbell at Brigham Young University.
Please acknowledge their work when using APOSTRUCT in research.
See [NOTICE](NOTICE) for details.
