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
the BYU tables (`byu=`), which `byu_export` and the verifier use.

Each table file starts with a tab-separated header; every later line is one entity, with tab-separated fields. A
number field is a space-separated list, a text field is the text.

| Table | Rows | Holds |
| --- | --- | --- |
| `group/point_ops.txt` | 72 | point operations: matrix, class code |
| `group/point_groups.txt` | 32 | point groups: symbol, order, class code, generators |
| `group/lattices.txt` | 14 | Bravais lattices: symbol, the conventional-to-primitive matrices, special k points |
| `group/space_groups.txt` | 230 | space groups: lattice, point group, BC symbol, default setting, generators, setting matrices |
| `group/magnetic_ops.txt` | 144 | magnetic point operations: ordinary part, time reversal, class code |
| `group/magnetic_point_groups.txt` | 122 | magnetic point groups: class code |
| `group/magnetic_groups.txt` | 1651 | magnetic space groups: space group, OG and BNS symbols, generators, type-4 supergroup row |
| `group/type4_supergroups.txt` | 520 | type-4 reference supergroups: group, basis, origin |
| `settings/inter_settings.txt` | 800 | International Tables settings: number, symbols, origin and cell choices, matrices, Wyckoff formulas |
| `settings/wyckoff.txt` | 1731 | Wyckoff positions: space group, letter, site point group, formula |
| `settings/site_point_groups.txt` | 32 | site point groups: irreps and vector representations |
| `settings/site_settings.txt` | 152 | orientations of the site point groups: operations, vector basis |
| `k/slots.txt` | 378 | k families per lattice: name, number of parameters, vector |
| `k/space_group_slots.txt` | 6210 | k families per space group: Kovalev number, conventional-to-ML matrix row |
| `k/star_flips.txt`, `k/star_orders.txt`, `k/conv2ml.txt` | | the choices among equivalent star arms |
| `irreps/little.txt` | 10300 | little irreps: space group, k, partner, Kovalev label, real form, matrices at the seed positions (the irrep basis) |
| `irreps/old.txt` | 4777 | irreps at k points without parameters: space group, Mulliken label, image |
| `isotropy/ordinary.txt`, `isotropy/magnetic.txt` | 15239, 16721 | the presentation of each isotropy subgroup: subgroup basis and origin, representative direction, and the magnetic OPD label (the ordinary one is its image direction's) |
| `images/*.txt` | | the image catalog: images, their subgroups and directions |

`const.dat` names the constants that matrix and direction entries refer to by number.

`bin/byu_export bin/source.blob DIR [--full-matrices]` writes the tables back in the BYU format (`data_*.txt`): every
section apo reads, including those it builds at startup, and with `--full-matrices` the little irrep matrices at every
position. They equal the corresponding BYU tables.
