/* Little-irrep sections built from other tables: little_irr_full_matrices_pointer, where each stored matrix starts in
 * the sparse triples (each matrix lists its entries by column, then row, so a new one starts where that order drops);
 * little_irr_full_dim, the number of entries of the stored identity matrix (0 for an irrep without matrices);
 * little_irr_old, the old (data_irreps) id, which the irreps at k points without parameters take in gid order. */

#include <stdlib.h>
#include <string.h>

#include "tables.h"

int section_find_stored(const char *name, section_view *view);

static struct {
    int ready, failed;
    int64_t *dims, *old, *pointer;
    blob_section entries[3];
} D;

static int build(void)
{
    section_view sg, kslot, k_dim, lattice, irr, position, triples;
    if (section_find_stored("little.little_irr_space_group", &sg) || section_find_stored("little.little_irr_k", &kslot) ||
        section_find_stored("little.little_k_dim", &k_dim) || section_find_stored("space.ispace_lattice", &lattice) ||
        section_find_stored("little.little_irr_full_matrices_irr_pointer", &irr) ||
        section_find_stored("little.little_irr_full_matrices_position", &position) ||
        section_find_stored("little.little_irr_full_matrices", &triples))
        return 1;
    int64_t count = section_count(&sg), irrs = section_count(&irr), positions = section_count(&position), olds = 0;
    int64_t entries = section_count(&triples) / 3, matrices = 0;
    D.dims = calloc((size_t)count + 1, sizeof(int64_t));
    D.old = calloc((size_t)count + 1, sizeof(int64_t));
    D.pointer = calloc((size_t)positions + 2, sizeof(int64_t));
    if (!D.dims || !D.old || !D.pointer) return 1;
    for (int64_t t = 0; t < entries; t++) {
        int64_t row = section_int(&triples, 3 * t), col = section_int(&triples, 3 * t + 1);
        int64_t prow = t ? section_int(&triples, 3 * t - 3) : 0, pcol = t ? section_int(&triples, 3 * t - 2) : 0;
        if (!t || col < pcol || (col == pcol && row <= prow)) {
            if (matrices == positions) return 1;
            D.pointer[matrices++] = t + 1;
        }
    }
    if (matrices != positions) return 1;
    D.pointer[matrices] = entries + 1;
    for (int64_t g = 0; g < count; g++) {
        int64_t s = section_int(&sg, g), k = section_int(&kslot, g);
        if (s < 1 || s > section_count(&lattice) || k < 1 || k > 27) continue;
        if (section_int(&k_dim, (section_int(&lattice, s - 1) - 1) * 27 + k - 1) == 0) D.old[g] = ++olds;
    }
    for (int64_t g = 0; g + 1 < irrs && g < count; g++) {
        int64_t first = section_int(&irr, g), end = section_int(&irr, g + 1);
        if (end <= first || first < 1 || first > positions || section_int(&position, first - 1) != 1) continue;
        D.dims[g] = D.pointer[first] - D.pointer[first - 1];
    }
    const char *names[3] = {"little.little_irr_full_dim", "little.little_irr_old", "little.little_irr_full_matrices_pointer"};
    int64_t counts[3] = {count, count, positions + 1};
    for (int i = 0; i < 3; i++) {
        strncpy(D.entries[i].name, names[i], BLOB_NAME_SIZE - 1);
        D.entries[i].kind = BLOB_I64;
        D.entries[i].count = (uint64_t)counts[i];
        D.entries[i].bytes = (uint64_t)counts[i] * sizeof(int64_t);
    }
    return 0;
}

/* section_find for the sections above; *handled 0 for any other name */
int little_derived_section(const char *name, section_view *view, int *handled)
{
    int which = -1;
    for (int i = 0; i < 3 && which < 0; i++)
        if (!strncmp(name, i == 0 ? "little.little_irr_full_dim" : i == 1 ? "little.little_irr_old" : "little.little_irr_full_matrices_pointer",
                     BLOB_NAME_SIZE))
            which = i;
    *handled = which >= 0;
    if (!*handled) return 1;
    if (!D.ready && !D.failed) {
        D.failed = build();
        D.ready = !D.failed;
    }
    if (!D.ready) return 1;
    view->entry = &D.entries[which];
    view->data = (const unsigned char *)(which == 0 ? D.dims : which == 1 ? D.old : D.pointer);
    return 0;
}

/* the dimension of little irrep gid, 0 when unknown */
int little_full_dim(int gid)
{
    section_view view;
    if (section_find("little.little_irr_full_dim", &view) || gid < 1 || gid > section_count(&view)) return 0;
    return (int)section_int(&view, gid - 1);
}

/* the old (data_irreps) id of little irrep gid, 0 when it has none */
int little_old_id(int gid)
{
    section_view view;
    if (section_find("little.little_irr_old", &view) || gid < 1 || gid > section_count(&view)) return 0;
    return (int)section_int(&view, gid - 1);
}
