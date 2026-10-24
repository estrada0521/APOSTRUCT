<div align="center">

# APOSTRUCT

An engine for crystal symmetry analysis.
Runs ISODISTORT-equivalent calculations locally.

[Manual](Manual.md) · [Validation](Validation.md) · [日本語](README.jp.md)

</div>

<p align="center">
  <img src="https://media.githubusercontent.com/media/estrada0521/APOSTRUCT/b76ddd9/media/APOSTRUCT.png" width="100%" alt="APOSTRUCT graphical interface">
</p>

---

## Overview

From a parent CIF, a choice of displacive, magnetic and strain degrees of freedom, and k points with irreps, it computes the OPD list and, for a selected OPD, the subgroup, child structure and mode details. Labels and the choices that mathematics leaves open (irrep bases, settings, representative directions, ...) follow the conventions of the standard crystallographic references and of the ISOTROPY Software Suite, as its data records them.

## How to use

Requires a C compiler (macOS, Linux). The build writes `bin/apo` and `bin/source.blob`, which holds the tables. For Windows, build with MinGW-w64 (`CC=x86_64-w64-mingw32-gcc sh build.sh`, from MSYS2 or as a cross build) and keep `bin/apo.exe` next to `bin/source.blob`.

```sh
sh build.sh
ln -s "$PWD/bin/apo" ~/.local/bin/apo
```

CLI example (SrTiO3, Sr(1a), Ti(1b), O(3c)). See the [Manual](Manual.md) for details.

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
| `core/` | Computation engine (C) |
| `cli/` | Command line |
| `gui/` | `apo serve` and the browser interface |
| `source/` | Basic static data and tables of BYU conventions |
| `verification/` | Verifier (Python): comparison with the ISODISTORT Web application and mathematical checks independent of the Web. Results in [Validation.md](Validation.md) |
| `branches/` `assets/` | Judgement records and their CIFs |

## Notice

APOSTRUCT output follows the conventions of Harold T. Stokes, Dorian M. Hatch and Branton J. Campbell (Brigham Young University, hereafter BYU). Please acknowledge their work when using APOSTRUCT in research. See [NOTICE](NOTICE) for details.
