# Distortropy

[日本語](README_ja.md) | [CLI guide](CLI.md)

**Distortropy** is a local Python package for symmetry analysis of crystalline
distortions. It computes k vectors, irreducible representations,
order-parameter directions, isotropy subgroups, invariant bases, and
symmetry-adapted modes from the bundled BYU ISOTROPY tables. Normal calculations
do not call the BYU web service or an external executable.

<p align="center">
  <img src="media/distortropy.png" alt="Distortropy graphical interface" width="70%">
</p>

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
- Interfaces: the `distortropy` command and a local browser interface.

The CLI and graphical interface are entry points to the same backend services.
They differ in workflow and presentation, not in the scientific calculation
path.

## Install

Clone the repository and install the package:

```bash
python -m pip install .
```

For development, use `python -m pip install -e .` instead. Both forms install
the `distortropy` command. The required group-theory tables are included in the
package.

## Quick Start

Use identifiers returned by one stage as input to the next:

```bash
# Inspect the parent and its selectable crystallographic sites.
distortropy info structure.cif

# Discover reciprocal-space and representation choices.
distortropy kpoints structure.cif
distortropy irreps structure.cif --k R --displacive O

# Enumerate OPDs, then compute one exact returned label.
distortropy opds structure.cif --k R --irrep R4- --displacive O
distortropy modes structure.cif \
  --k R --irrep R4- --displacive O --opd P1

# Compute the selected-domain Landau invariant basis.
distortropy invariants structure.cif \
  --k R --irrep R4- --displacive O --opd P1 \
  --minimum-degree 2 --maximum-degree 6
```

When no concrete structure is needed, use the space-group/Wyckoff route:

```bash
distortropy opds \
  --sg 221 --wyckoff 1a 1b 3c \
  --k R --irrep R4- --displacive c
```

Free Wyckoff coordinates may remain unspecified through OPD and invariant
calculations. They are required only when atomic mode geometry must be
realized.

The [CLI guide](CLI.md) defines the complete command workflow, input
conventions, coupled selections, parametric k points, saved cases, exact
embedding reuse, secondary invariant factors, and machine-readable output.

## Graphical Interface

```bash
distortropy serve --host 127.0.0.1 --port 8300 --open-browser
```

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
University. Please acknowledge their work when using Distortropy in research.
See [NOTICE](NOTICE) for details.
