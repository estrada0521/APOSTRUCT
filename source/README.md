# Source

The information apo cannot compute: group data, names, the BYU irrep basis and BYU's presentation choices. It is the
data of the BYU ISOTROPY Software Suite (see [NOTICE](../NOTICE)), reorganized as one table per kind of entity; what apo
computes from it (group elements and products, isotropy subgroups, irrep dimensions and types, the star of k, labels,
the pointers between tables, ...) is left out. `build.sh` compiles the tables into `bin/source.blob`
(`blob_build.c`), and `layout.txt` says how each table column becomes a blob section.

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
| `k/star_flips.txt`, `k/star_orders.txt`, `k/conv2ml.txt` | | BYU's choices among equivalent star arms |
| `irreps/little.txt` | 10300 | little irreps: space group, k, partner, Kovalev label, real form, matrices at the seed positions (the BYU basis) |
| `irreps/old.txt` | 4777 | irreps at k points without parameters: space group, Mulliken label, image |
| `isotropy/ordinary.txt`, `isotropy/magnetic.txt` | 15239, 16721 | BYU's presentation of each isotropy subgroup: OPD label, subgroup basis and origin, representative direction |
| `images/*.txt` | | the image catalog: images, their subgroups and directions |

`const.dat` names the constants that matrix and direction entries refer to by number.
