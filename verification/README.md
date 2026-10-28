# verification

- `web/` reads the parent, subgroup and undistorted-structure sections of the Web "Complete Mode Details" text.
- `compare/` compares an apo output with that text: identical, equivalent or unproven (`compare_local_to_web`). The group-theoretic frame (OPD, subgroup, lattices, atoms) is Python; the mode definitions, from both texts as read in place, are compared in C (`definitions.py`, `native/definitions.c`).
- `mathematics/` checks an output on its own: mode basis, group invariance, fixed-space dimension (`fast.py` runs them in C, `native/checks.c`; the Python modules are the reference).
- `symmetry/` holds the group, setting and Source-table machinery both use; `apo_kernels.py` is the few apo C kernels it calls.
- `json.h` is the JSON cursor both C parts read apo outputs with.
- `verify_comparison.py`, `verify_local_output.py`, `verify_web_text.py` are the commands; `version.py` names the comparator by its sources.
