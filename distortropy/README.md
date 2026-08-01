# Distortropy Runtime

This package implements Distortropy's local crystal-symmetry pipeline in Python.
It accepts either a CIF or a symbolic space-group/Wyckoff parent and computes
the parent symmetry, reciprocal-space choices, isotropy subgroups, child
structures, and symmetry-adapted mode definitions:

```text
CIF input or symbolic SG/Wyckoff input
-> parent structure and setting
-> distortion types (strain and crystallographic-site masks)
-> k-point choices
-> irrep choices
-> OPD/subgroup choices
-> mode-detail object
```

## Local UI

The frontend is exposed as a first-class localhost tool:

```text
/tools/distortropy/
```

The UI accepts a CIF or ordered Wyckoff sites for a selected space group and
presents the calculation as a sequence of choices:

1. parent structure and input sites;
2. k-vector choices;
3. irrep choices for the selected k point;
4. OPD/subgroup choices for the selected irrep;
5. mode details after an OPD is explicitly selected.

These are UI states rather than backend layers. Local execution is
self-contained. Backend ownership follows the scientific data flow:

- `Backend/source/`: typed access to distributed Source tables, with no
  dependency on downstream domains;
- `Backend/parent/`: CIF and symbolic parent inputs plus displayed-setting state;
- `Backend/reciprocal/`: k-vector and irrep selection;
- `Backend/isotropy/`: OPD, subgroup, domain, and Source isotropy operations;
- `Backend/modes/`: subduction, projection, structure, and presentation;
- `Backend/pipeline.py`: composition only; it must not reimplement domain
  calculations.
- `Frontend/`: localhost UI assets served by the application shell.

Dependencies flow downward in that order. Shared Source decoding must not
import parent, reciprocal, isotropy, or mode behavior; domain-specific catalog
logic stays with the domain that consumes it.

All domains share the decoder in `Backend/source/iso_data.py`. A domain may
wrap those tables with typed records, but it must not carry a second Source
file parser or an independently evolving table cache.

Source tables used by the runtime belong under the shared root `Source/`. Runtime code
imports through `distortropy/Backend/`.

Runtime boundary:

- Distortropy does not invoke external crystallographic executables.
- Scientific calculations live inside the relevant backend domain rather than
  in the application shell or frontend.
- Generated validation data and diagnostic artifacts are not runtime
  dependencies.

## Command Line

The command-line interface emits compact JSON for structure, k-point, irrep,
OPD, and mode results. The complete pipeline state remains available explicitly:

```bash
distortropy kpoints structure.cif
distortropy irreps structure.cif --k DT b=1/3 --magnetic Fe
distortropy opds structure.cif --k R --irrep R4- --displacive O
distortropy modes structure.cif --k R --irrep R4- --displacive O --opd P1
distortropy modes structure.cif --k R --irrep R4- --displacive O --opd P1 --full-state
distortropy invariants structure.cif --k R --irrep R4- --displacive O --opd P2
distortropy opds --sg 221 --wyckoff 1a 1b 3c --k R --irrep R4- --displacive c
distortropy modes --sg 2 --wyckoff i:x=1/7,y=2/11,z=3/13 --k GM --irrep GM1+ --displacive i --opd P1
```

The symbolic route does not invent free Wyckoff coordinates or cell metrics
for K, irrep, OPD, or invariant calculations. Free site parameters may remain
absent through those stages. Mode details require them because atom positions
must then be realized.

Site-mode selectors accept either the unique `info.sites[].label` or its
Wyckoff/element `type`; a type selects every matching crystallographic site.

Case input is JSON. It may name a CIF path with `structure`, a repository asset
with its 12-character `cif` content ID, or `sg` with ordered `wyckoff` items.
Each Wyckoff item may be a label or an object containing `wyckoff` and optional
`parameters`. A case can be read from a file or standard input, and `run` stops
after one requested pipeline stage without invoking preceding commands as
subprocesses:

```bash
distortropy run --case case.json --upto opds
cat case.json | distortropy run --case - --upto modes
```

Mode results are JSON unless `--format text` is requested. The command layer
normalizes inputs and projects results; scientific calculations remain in the
backend domains listed above.

`invariants` uses the selected primary OPD domains and the same internal
invariant service as the graphical interface. Secondary irreps, full-IR
invariants, and gradient invariants are not implicit in this command.

## Current Scope

The active pipeline supports:

```text
one CIF or symbolic parent
-> strain on/off and displacive/magnetic crystallographic-site masks
-> ordered single or multiple K/IR selections
-> OPD table
-> selected mode details
```

Strain, magnetic, and coupled multi-IR paths share the same domain services.
