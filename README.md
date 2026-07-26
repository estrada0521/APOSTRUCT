# Distortropy

[日本語](README_ja.md)

**Distortropy** is a pure-Python, offline implementation of **ISODISTORT**.

<p align="center">
  <img src="media/distortropy-1.png" alt="Distortropy structure selection" width="49%">
  <img src="media/distortropy-2.png" alt="Distortropy mode viewer" width="49%">
</p>

## Setup

Clone the repository, then install its dependencies and the `isodistort`
command in editable mode:

```bash
python -m pip install -e .
```

Place the following nine group-theory table files in `Source/` at the
repository root. They are distributed with the BYU ISOTROPY Software Suite
<https://iso.byu.edu/>.

```text
const.dat  data_isotropy.txt  data_irreps.txt  data_images.txt  data_little.txt
data_magnetic.txt  data_space.txt  data_ssgmag.txt  data_wyckoff.txt
```

## How to Use

A structure is provided either as a CIF file, or — for quick queries without a
file — directly by space group and Wyckoff positions via `--sg` / `--wyckoff`,
accepted anywhere a CIF is.

```bash
# 0. Show the parsed structure: space group, sites, Wyckoff positions
isodistort info structure.cif
isodistort info --sg 205 --wyckoff a c        # on the fly, no CIF

# 1. Enumerate k points for the parent space group
isodistort kpoints structure.cif
isodistort kpoints --sg 205

# 2. Enumerate irreducible representations at selected k points
isodistort irreps structure.cif \
  --k L \
  --k GP 1/3 1/4 2/5 \
  --k B 1/3 2/5 \
  --k W 1/4
# without a CIF: isodistort irreps --sg 62 --wyckoff a c --k L --k GP 1/3 1/4 2/5 ...

# 3. Enumerate order-parameter directions (OPDs) for selected irreps
isodistort opds structure.cif \
  --k L --irrep L1 \
  --k GP 1/3 1/4 2/5 --irrep GP1GQ1 \
  --k B 1/3 2/5 --irrep mB2BA2 \
  --k W 1/4 --irrep W1WA1 \
  --displacive Sn Fe --magnetic O Fe --strain

# 4. Compute complete mode details for a selection
isodistort modes structure.cif \
  --k L --irrep L1 \
  --k GP 1/3 1/4 2/5 --irrep GP1GQ1 \
  --k B 1/3 2/5 --irrep mB2BA2 \
  --k W 1/4 --irrep W1WA1 \
  --displacive Sn Fe --magnetic O Fe --strain \
  --opd 'P1(1)P3(1)C2(1)P2(1)'
```

### JSON case input

```json
{ "structure": "structure.cif",
  "k": [
    { "label": "L", "ir": "L1" },
    { "label": "GP", "params": { "a": "1/3", "b": "1/4", "g": "2/5" },
      "ir": "GP1GQ1" },
    { "label": "B", "params": { "a": "1/3", "g": "2/5" }, "ir": "mB2BA2" },
    { "label": "W", "params": { "g": "1/4" }, "ir": "W1WA1" }
  ],
  "sites": {
    "Sn": { "displacive": true,  "magnetic": false },
    "O":  { "displacive": false, "magnetic": true },
    "Fe": { "displacive": true,  "magnetic": true },
    "Ba": { "displacive": false, "magnetic": false }
  },
  "strain": true,
  "opd": "P1(1)P3(1)C2(1)P2(1)"
}
```

```bash
isodistort modes --case case.json
```

### `.in` case input

```text
CIF structure.cif

K L
IR L1
K GP 1/3 1/4 2/5
IR GP1GQ1
K B 1/3 2/5
IR mB2BA2
K W 1/4
IR W1WA1

STRAIN
DISPLACIVE Sn Fe
MAGNETIC O Fe
OPD P1(1)P3(1)C2(1)P2(1)
```

```bash
isodistort modes --case case.in
```

## Output

`kpoints`, `irreps`, and `opds` write JSON to standard output. For `modes`, a
JSON case produces JSON, while a `.in` case or direct command-line selection
produces complete mode-details text. Use `modes --format json|text` to override
the format and `-o PATH` to save the result to a specific file. No output file
is created unless `-o` is provided.

## Server

Start the standalone local server to use the graphical interface:

```bash
isodistort serve --host 127.0.0.1 --port 8300
```

## Validation

Validation against the reference implementation is ongoing. For practical
inputs, physically non-equivalent output is uncommon, but exact equivalence is
not guaranteed. Differences are more likely for heavy, nontrivial cases such
as the multi-k example above.

## Notice

The original ISODISTORT software and the group-theory tables are the work of
Harold T. Stokes, Dorian M. Hatch, and Branton J. Campbell at Brigham Young
University. Please acknowledge their work when using Distortropy in research.
See [NOTICE](NOTICE) for details.
