<div align="center">

# APOSTRUCT

An engine for crystal symmetry analysis.
Runs ISODISTORT-equivalent calculations locally.

[Manual](Manual.md) · [Validation](Validation.md) · [NOTICE](NOTICE)

</div>

<p align="center">
  <img src="https://media.githubusercontent.com/media/estrada0521/APOSTRUCT/b76ddd9/media/APOSTRUCT.png" width="100%" alt="APOSTRUCT graphical interface">
</p>

---

## Overview

From a parent CIF, computes irreducible representations, order parameter directions (OPDs), subgroups, child structures and modes for displacements, magnetism, occupancy and strain.

## Layout

| Path | Content |
| --- | --- |
| `core/` | Computation engine (C) |
| `cli/` | Command line |
| `gui/` | `apo serve` and the browser interface |
| `source/` | Basic static data and convention tables |
| `verification/` | Verifier (Python): comparison with ISODISTORT output and independent mathematical checks. Results in [Validation.md](Validation.md) |
| `branches/` · `assets/` | Judgement records and their CIFs |

## Setup

- macOS and Linux: build with a C compiler.

  ```sh
  sh build.sh
  ln -s "$PWD/bin/apo" ~/.local/bin/apo
  ```

- Windows: use MinGW-w64 in MSYS2 or cross-compile.

  ```sh
  CC=x86_64-w64-mingw32-gcc sh build.sh
  ```

The build generates `bin/apo` (`bin/apo.exe` on Windows) and `bin/source.blob`, which contains the tables. Keep the executable and `source.blob` in the same directory.

## Example

SrTiO3 (Sr[1a], Ti[1b], O[3c]):

```sh
apo info structure.cif  # Inspect the parent structure and sites
apo kpoints structure.cif  # List k points
apo irreps structure.cif --k R --displacive O  # List irreps for oxygen displacements at R
apo opds structure.cif --k R --irrep R5- --displacive O  # List OPDs and subgroups for R5-
apo modes structure.cif --k R --irrep R5- --displacive O --opd P1  # Compute the child structure and displacement modes for P1
```

See the [Manual](Manual.md) for details.
