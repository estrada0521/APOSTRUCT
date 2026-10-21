# Source

The information apo cannot compute, in two parts:

- `crystallography/`: the mathematics and the labels of the crystallographic literature, which do not depend on BYU's
  choices: point operations and groups, lattices, space-group operations, magnetic groups with their OG and BNS symbols,
  the International Tables settings and Wyckoff positions, the k families (CDML names, Kovalev numbers) and the Kovalev
  and Mulliken irrep labels.
- `byu/`: the conventions of the BYU ISOTROPY Software Suite (see [NOTICE](../NOTICE)), where one of several equivalent
  answers is chosen: the basis of each little irrep and its real form, the complex-conjugate partner numbering, the
  default setting of each space group, the orientation bases of the site point groups, the choices among equivalent star
  arms and ML matrices, the image catalog, and the presentation of each isotropy subgroup (OPD label, basis and origin,
  representative direction).

What apo computes from these (group elements and products, isotropy subgroups, irrep dimensions and types, the star of k,
labels, the pointers between tables, ...) is left out. Nothing here is a cache: derivations run at startup. A cache of
something derivable but too slow to compute on every run would be written at build time by a generator in this
directory.

One mixture remains: the subgroup identification (`core/group/id_subgroup.c`) tries the generators of
`group/space_groups.txt` and `group/magnetic_groups.txt` in their listed order and keeps the first basis that fits, so
the choice and order of the generators also decide which of the equivalent settings of a child structure is shown (with
the space-group generators reversed, 15 of 4,141 outputs change). Beyond the group they generate, they carry a BYU
choice.

`build.sh` compiles the tables into `bin/source.blob` (`blob_build.c`), and `layout.txt` says how each table column
becomes a blob section. A table found under both parts holds different columns of the same records, with the same key in
its first column. A blob section is named `directory.table.column` (with `_start` and `_count` for the pointers the build
renders, and `by_COLUMN_start` / `by_COLUMN_count` for the records grouped by a key); `layout.txt` also records each
section's name in the BYU tables (`byu=`), which `byu_export` and the verifier use.

Each table file starts with a tab-separated header; every later line is one entity, with tab-separated fields. A
number field is a space-separated list, a text field is the text.

| Table | Rows | `crystallography/` | `byu/` |
| --- | --- | --- | --- |
| `group/point_ops.txt` | 72 | matrix, class code | |
| `group/point_groups.txt` | 32 | symbol, order, class code, generators | |
| `group/lattices.txt` | 14 | symbol, conventional-to-primitive matrices, special k points | |
| `group/space_groups.txt` | 230 | lattice, point group, BC symbol, generators, setting matrices | default setting |
| `group/magnetic_ops.txt` | 144 | ordinary part, time reversal, class code | |
| `group/magnetic_point_groups.txt` | 122 | class code | |
| `group/magnetic_groups.txt` | 1651 | space group, OG and BNS symbols, generators, type-4 reference row | |
| `group/type4_supergroups.txt` | 520 | OG reference supergroup, basis, origin | |
| `settings/inter_settings.txt` | 800 | number, symbols, origin and cell choices, matrices, Wyckoff formulas | |
| `settings/wyckoff.txt` | 1731 | space group, letter, site point group, formula | |
| `settings/site_point_groups.txt` | 32 | irreps, vector representations | |
| `settings/site_settings.txt` | 152 | point group, operations of each orientation | vector basis |
| `k/slots.txt` | 378 | k families per lattice: name, number of parameters, vector | |
| `k/space_group_slots.txt` | 6210 | Kovalev number | conventional-to-ML matrix row |
| `k/star_flips.txt`, `k/star_orders.txt`, `k/conv2ml.txt` | | | choices among equivalent star arms, ML matrices |
| `irreps/little.txt` | 10300 | space group, k, Kovalev label | partner, real form, matrices at the seed positions |
| `irreps/old.txt` | 4777 | space group, Mulliken label | image |
| `isotropy/ordinary.txt`, `isotropy/magnetic.txt` | 15239, 16721 | | OPD label, subgroup basis and origin, representative direction |
| `images/*.txt` | | | the image catalog: images, their subgroups and directions |

`const.dat` names the constants that matrix and direction entries refer to by number.

`bin/byu_export bin/source.blob DIR [--full-matrices]` writes the tables back in the BYU format (`data_*.txt`): every
section apo reads, including those it builds at startup, and with `--full-matrices` the little irrep matrices at every
position. They equal the corresponding BYU tables.
