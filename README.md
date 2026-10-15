<div align="center">

# APOSTRUCT

Local implementation of ISODISTORT

[CLI Guide](CLI.md) · [Validation](Validation.md) · [日本語](README_ja.md)

</div>

<p align="center">
  <img src="media/APOSTRUCT.png" width="100%" alt="APOSTRUCT graphical interface">
</p>

---

## Overview

From a parent CIF, a choice of displacive, magnetic and strain degrees of freedom, and k points with irreps, it computes the OPD list and, for a selected OPD, the subgroup, child structure and mode details.

## How to use

Requires a C compiler (macOS, Linux). The build writes `bin/apo` and `bin/source.blob`, which holds the tables.

```sh
sh build.sh
ln -s "$PWD/bin/apo" ~/.local/bin/apo
```

CLI example (SrTiO3, Sr(1a), Ti(1b), O(3c)). See the [CLI Guide](CLI.md) for details.

```sh
apo info structure.cif
apo kpoints structure.cif
apo irreps structure.cif --k R --displacive O
apo opds structure.cif --k R --irrep R5- --displacive O
apo modes structure.cif --k R --irrep R5- --displacive O --opd P1
```

The GUI opens at `http://127.0.0.1:8300`.

```sh
apo serve --open-browser
```

## Layout

| Path | Content |
| --- | --- |
| `core/` | Computation (C) |
| `cli/` | Command line |
| `gui/` | `apo serve` and the browser interface |
| `source/` | BYU ISOTROPY tables and the blob generator |
| `verification/` | Verifier (Python): comparison with the ISODISTORT Web application and mathematical checks independent of the Web. Results in [Validation.md](Validation.md) |
| `branches/` `assets/` | Judgement records and their CIFs |

## Notice

The original ISODISTORT software and group-theory tables are the work of Harold T. Stokes, Dorian M. Hatch and Branton J. Campbell (Brigham Young University). Please acknowledge them when using this in research. See [NOTICE](NOTICE).
