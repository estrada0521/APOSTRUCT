# Source

The information apo cannot compute: the conventions its output is written in. They come from the data of the ISOTROPY
Software Suite (see [NOTICE](../NOTICE)), which gathers the standard crystallographic references (International Tables,
CDML, Kovalev, Bradley–Cracknell, OG and BNS) with the suite's own choices (irrep bases, settings, star arms, the image
catalog, the presentation of each isotropy subgroup), reorganized here as one table per kind of entity. What apo
computes from it (group elements and products, isotropy subgroups, irrep dimensions and types, the star of k, labels,
the pointers between tables, ...) is left out. `build.sh` compiles the tables into `bin/source.blob`
(`blob_build.c`), and `layout.txt` says how each table column becomes a blob section.

A blob section is named `directory.table.column` (with `_start` and `_count` for the pointers the build renders, and
`by_COLUMN_start` / `by_COLUMN_count` for the records grouped by a key); `layout.txt` also records each section's name in
the BYU tables (`byu=`), which the verifier uses. A 4x4 setting matrix is stored once; the build adds its inverse.

Each table file starts with a tab-separated header; every later line is one entity, with tab-separated fields. A
number field is a space-separated list, a text field is the text.

| Table | Rows | Holds |
| --- | --- | --- |
| `group/point_ops.txt` | 72 | point operations: matrix |
| `group/point_groups.txt` | 32 | point groups: symbol, generators |
| `group/lattices.txt` | 14 | Bravais lattices: symbol, the conventional-to-primitive matrices, number of k families common to the lattice |
| `group/space_groups.txt` | 230 | space groups: lattice, point group, BC symbol, generators, the BC-to-ITA setting matrix |
| `group/magnetic_point_groups.txt` | 122 | magnetic point groups: class code |
| `group/magnetic_groups.txt` | 1651 | magnetic space groups: space group, OG and BNS symbols, generators (translation and point operation, `'` for time reversal), and for type IV the reference space group with its basis and origin |
| `settings/inter_settings.txt` | 777 | International Tables settings: space group, whether it is the default, symbol, origin and cell choices, matrix |
| `settings/wyckoff_exceptions.txt` | 107 | Wyckoff formulas in other settings that differ from the ones apo derives from the default setting |
| `settings/wyckoff.txt` | 1731 | Wyckoff positions: space group, letter, site point group, formula |
| `settings/site_point_groups.txt` | 32 | site point groups: irreps and vector representations |
| `settings/site_settings.txt` | 152 | orientations of the site point groups: operations, vector basis |
| `k/slots.txt` | 212 | k families per lattice: name, base point and parameter directions (the number of parameters is the number of nonzero directions) |
| `k/space_group_slots.txt` | 4141 | k families per space group: Kovalev number, conventional-to-ML matrix row, and BYU's choices among equivalent star arms (arms stored with the opposite sign, a different arm order) |
| `k/conv2ml.txt` | 42 | the conventional-to-ML parameter matrices |
| `irreps/little.txt` | 10294 | little irreps: space group, k, partner, Kovalev label, real form, matrices at the seed positions (the irrep basis) |
| `irreps/old.txt` | 4777 | irreps at k points without parameters: space group, Mulliken label, image |
| `isotropy/ordinary.txt`, `isotropy/magnetic.txt` | 15239, 16721 | the presentation of each isotropy subgroup: subgroup basis and origin, representative direction (the OPD labels are those of image directions) |
| `images/images.txt`, `images/image_directions.txt` | 141, 689 | the image catalog: images with their directions, and each image with time reversal added |

`constants.txt` holds the constants that matrix and direction entries refer to: entry `k` means constant `k`, `-k` its
negative, and `0` no constant. Each has a decimal value and, when it lies in Q(√2, √3), an exact form; the ones without
a form are the coefficients of representative directions.

`bin/source_dump bin/source.blob FILE` writes what apo reads from the Source, independent of how the tables encode it:
every section, including those apo builds at startup, with constants written as their values, and the little irrep
matrices at every position.
