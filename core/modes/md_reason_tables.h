/* Python exception types and texts of the native error codes (generated once from native/__init__.py). */

typedef struct {
    int status;
    const char *type, *text;
} md_reason_entry;

static const md_reason_entry child_layout_errors[] = {
    {1, "MemoryError", "child layout ran out of memory"},
    {2, "ArithmeticError", "child layout value out of range"},
    {3, "ValueError", "Formula15 child sites are inconsistent"},
    {4, "ValueError", "Formula15 occurrence ordinals are not contiguous"},
    {5, "ValueError", "Formula15 does not cover the occurrence product"},
    {6, "ValueError", "Formula15 flat occurrence order is inconsistent"},
    {7, "ValueError", "Formula15 parent coset identity is inconsistent"},
    {8, "ValueError", "Formula15 child-site membership differs from its multiplicity"},
    {9, "ValueError", "Formula15 produced duplicate child atom identity"},
    {10, "ValueError", "Source atom table does not bijectively cover canonical child atoms"},
    {11, "ValueError", "Source atom table contains a nonfinite position"},
    {12, "ValueError", "invalid atom operation record"},
    {13, "ValueError", "child layout capacity exceeded"},
    {14, "ValueError", "Source atom positions and operation records differ"},
    {15, "ValueError", "magnetic presentation atoms share one position"},
    {16, "ValueError", "magnetic child site does not partition canonical atoms"},
    {17, "ValueError", "magnetic child-site representative is not one canonical atom"},
    {0, 0, 0},
};

static const md_reason_entry wyckoff_split_errors[] = {
    {1, "ArithmeticError", "Wyckoff split value out of range"},
    {2, "KeyError", "child operation image not found"},
    {3, "ValueError", "child Wyckoff orbit collision"},
    {4, "KeyError", "no child Wyckoff candidate"},
    {5, "ValueError", "could not solve child formula"},
    {6, "MemoryError", "Wyckoff split capacity exceeded"},
    {7, "TypeError", "child Wyckoff orbit has no site point group"},
    {0, 0, 0},
};

static const md_reason_entry child_site_errors[] = {
    {1, "ValueError", "Source Wyckoff geometry is inconsistent"},
    {2, "ArithmeticError", "child site arithmetic out of range"},
    {3, "ValueError", "parent Wyckoff row id does not identify one Source row"},
    {4, "ValueError", "Source Wyckoff split produced no Formula15 rows"},
    {5, "ValueError", "Formula15 did not partition the presentation atoms"},
    {6, "MemoryError", "child site expansion ran out of memory"},
    {7, "ValueError", "child site output capacity exceeded"},
    {8, "ValueError", "zero denominator"},
    {0, 0, 0},
};

static const md_reason_entry child_setting_errors[] = {
    {1, "ValueError", "child setting Source tables are inconsistent"},
    {2, "ArithmeticError", "child setting value out of range"},
    {3, "ValueError", "child presentation transform is singular"},
    {4, "ValueError", "child presentation transform does not define a conventional lattice"},
    {5, "ValueError", "child space-group action contains duplicate affine operations"},
    {6, "ValueError", "child setting transform collapses affine operations"},
    {7, "ValueError", "displayed child action matches no Source inter setting"},
    {8, "ValueError", "Source Wyckoff parameters are linearly dependent"},
    {9, "ValueError", "Formula15 child site contains duplicate points"},
    {10, "ValueError", "Source child-setting Wyckoff identity is not unique"},
    {11, "ValueError", "equivalent Source child settings disagree on one Wyckoff identity"},
    {12, "ValueError", "child setting capacity exceeded"},
    {0, 0, 0},
};

static const md_reason_entry site_atom_errors[] = {
    {1, "ValueError", "site atom table Source tables are inconsistent"},
    {2, "ArithmeticError", "site atom table value out of range"},
    {3, "ValueError", "singular basis matrix"},
    {4, "ValueError", "parent representative has no Source orbit transport"},
    {5, "ValueError", "Source site cosets overlap"},
    {6, "ValueError", "Source atom operations are not unique in the supercell"},
    {7, "ValueError", "Source representative order is not an atom permutation"},
    {8, "KeyError", "rebased atom operation has no Source atom"},
    {0, 0, 0},
};

static const md_reason_entry row_assignment_errors[] = {
    {1, "ValueError", "invalid atom operation record"},
    {2, "ValueError", "child centering ordinal out of range"},
    {3, "ValueError", "presentation basis is singular"},
    {4, "ValueError", "canonical child atoms share one affine occurrence identity"},
    {5, "ValueError", "mode row did not identify one canonical child atom"},
    {6, "MemoryError", "row assignment ran out of memory"},
    {7, "ArithmeticError", "row assignment value out of range"},
    {0, 0, 0},
};

static const md_reason_entry magnetic_errors[] = {
    {1, "ValueError", "magnetic Source tables are inconsistent"},
    {2, "IndexError", "magnetic value out of range"},
    {3, "MemoryError", "magnetic kernel ran out of memory"},
    {4, "KeyError", "magnetic point operation not found"},
    {5, "ValueError", "dependent magnetic Wyckoff parameter slots"},
    {6, "KeyError", "no magnetic Wyckoff formula"},
    {7, "ValueError", "canonical atoms share one child-cinter point"},
    {8, "ValueError", "ordinary child sites do not partition magnetic atoms"},
    {9, "KeyError", "magnetic Source operation leaves the canonical atom population"},
    {10, "ValueError", "magnetic operation is not an atom permutation"},
    {11, "ValueError", "magnetic atom component splits an ordinary child site"},
    {12, "ValueError", "canonical magnetic orbit does not equal atom membership"},
    {13, "ValueError", "magnetic representative does not identify one canonical atom"},
    {14, "ValueError", "singular magnetic basis"},
    {15, "ValueError", "magnetic kernel capacity exceeded"},
    {0, 0, 0},
};

static const md_reason_entry orderparam_errors[] = {
    {101, "KeyError", "operation record not found in generated space-group records"},
    {102, "ValueError", "orderparam_to_subgroup_ has no Source answer"},
    {103, "ArithmeticError", "orderparam_to_subgroup_ value out of range"},
    {104, "ZeroDivisionError", "operation record has a zero denominator"},
    {105, "ValueError", "Source little irrep tables are inconsistent"},
    {106, "ValueError", "the little irrep has no native irrep tables"},
    {0, 0, 0},
};

static const md_reason_entry recordwise_errors[] = {
    {1, "MemoryError", "recordwise projection ran out of memory"},
    {2, "ValueError", "bridge Source family has no active basis values"},
    {3, "ValueError", "recordwise projection Source tables are missing"},
    {4, "IndexError", "tuple index out of range"},
    {5, "ValueError", "recordwise projection capacity exceeded"},
    {6, "ValueError", "recordwise projection has no irrep matrices"},
    {7, "RuntimeError", "native bridge source modes failed"},
    {0, 0, 0},
};

static const md_reason_entry projection_entry_errors[] = {
    {1, "ValueError", "projection Source tables are inconsistent"},
    {2, "KeyError", "projection Source entry not found"},
    {3, "ValueError", "projection Source entry is not integral"},
    {4, "ValueError", "projection entry capacity exceeded"},
    {0, 0, 0},
};

static const md_reason_entry mode_block_errors[] = {
    {1, "MemoryError", "mode block projection ran out of memory"},
    {2, "ValueError", "recordwise projection produced no Source modes"},
    {3, "TypeError", "int() argument must be a string, a bytes-like object or a real number, not 'NoneType'"},
    {0, 0, 0},
};

static const md_reason_entry kernel_errors[] = {
    {1, "MemoryError", "native isotropy traversal ran out of memory"},
    {2, "ValueError", "orderparam_to_eqs_: OPD row has no pivot"},
    {3, "ValueError", "orderparam_to_eqs_: reduced row has no pivot"},
    {4, "ValueError", "orderparam_to_eqs_: duplicate pivot"},
    {5, "ValueError", "orderparam_to_eqs_: invalid column state"},
    {6, "ValueError", "orderparam_to_eqs_: inconsistent equation count"},
    {7, "ValueError", "eqs_to_orderparam_: equation row has no pivot"},
    {8, "ValueError", "eqs_to_orderparam_: zero pivot"},
    {9, "ValueError", "eqs_to_orderparam_: inconsistent free dimension"},
    {10, "ValueError", "fixed-space intersection changed free dimension"},
    {11, "KeyError", "kernel fraction coset not found"},
    {0, 0, 0},
};

static const md_reason_entry table_errors[] = {
    {1, "ArithmeticError", "traversal table value out of range"},
    {2, "KeyError", "kernel fraction coset not found"},
    {3, "KeyError", "operation record not represented by ispace/kernel pair"},
    {4, "KeyError", "point-operation product not found"},
    {5, "ValueError", "kernel fractions are not distinct modulo the kernel lattice"},
    {0, 0, 0},
};

static const md_reason_entry id_subgroup_errors[] = {
    {-7, "ZeroDivisionError", "id_subgroup_ operation with zero denominator"},
    {-6, "ValueError", "id_subgroup_ input exceeds capacity"},
    {-5, "KeyError", "id_subgroup_ operation has no Source record"},
    {-4, "ValueError", "id_subgroup_ Source tables are inconsistent"},
    {-3, "ValueError", "id_subgroup_ matrix did not reduce to integer"},
    {-2, "ArithmeticError", "id_subgroup_ value out of range"},
    {-1, "ValueError", "rowop matrix dump is too narrow"},
    {0, 0, 0},
};

static const md_reason_entry bridge_weight_errors[] = {
    {1, "KeyError", "operation record not found in generated space-group records"},
    {2, "ValueError", "phase count does not divide the irrep dimension"},
    {3, "ValueError", "non-real scalar phase"},
    {4, "ValueError", "odd real block cannot hold complex phase"},
    {5, "ArithmeticError", "bridge weight phase out of range"},
    {0, 0, 0},
};

static const md_reason_entry subgroup_ops_errors[] = {
    {1, "KeyError", "operation record not represented in the parent space group"},
    {2, "ValueError", "singular subgroup basis"},
    {3, "RuntimeError", "generate_domains2 point-coset scan did not close"},
    {4, "ArithmeticError", "subgroup operation value out of range"},
    {5, "ValueError", "subgroup operation capacity exceeded"},
    {6, "ZeroDivisionError", "operation record has a zero denominator"},
    {7, "ValueError", "subgroup operation Source tables are inconsistent"},
    {0, 0, 0},
};

static const char *const block_plan_errors[] = {
    0,
    "mode topology requires Source rows and centerings",
    "mode topology does not cover every canonical child atom",
    "mode presentation basis is singular",
    "mode rows do not cover the canonical presentation layout",
    "mode rows do not bijectively cover canonical child atoms",
};
