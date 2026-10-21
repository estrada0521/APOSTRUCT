/* Little-irrep sections built from other tables: little_irr_full_dim, the number of entries of the stored identity
 * matrix (0 for an irrep without matrices); little_irr_old, the old (data_irreps) id, which the irreps at k points
 * without parameters take in gid order. */

#include <stdlib.h>
#include <string.h>

#include "tables.h"

int section_find_stored(const char *name, section_view *view);

static struct {
    int ready, failed;
    int64_t *dims, *old;
    blob_section entries[2];
} D;

static int build(void)
{
    section_view sg, kslot, k_dim, lattice, irr, position, pointer;
    if (section_find_stored("little.little_irr_space_group", &sg) || section_find_stored("little.little_irr_k", &kslot) ||
        section_find_stored("little.little_k_dim", &k_dim) || section_find_stored("space.ispace_lattice", &lattice) ||
        section_find_stored("little.little_irr_full_matrices_irr_pointer", &irr) ||
        section_find_stored("little.little_irr_full_matrices_position", &position) ||
        section_find_stored("little.little_irr_full_matrices_pointer", &pointer))
        return 1;
    int64_t count = section_count(&sg), irrs = section_count(&irr), positions = section_count(&position), olds = 0;
    D.dims = calloc((size_t)count + 1, sizeof(int64_t));
    D.old = calloc((size_t)count + 1, sizeof(int64_t));
    if (!D.dims || !D.old) return 1;
    for (int64_t g = 0; g < count; g++) {
        int64_t s = section_int(&sg, g), k = section_int(&kslot, g);
        if (s < 1 || s > section_count(&lattice) || k < 1 || k > 27) continue;
        if (section_int(&k_dim, (section_int(&lattice, s - 1) - 1) * 27 + k - 1) == 0) D.old[g] = ++olds;
    }
    for (int64_t g = 0; g + 1 < irrs && g < count; g++) {
        int64_t first = section_int(&irr, g), end = section_int(&irr, g + 1);
        if (end <= first || first < 1 || first > positions || section_int(&position, first - 1) != 1) continue;
        D.dims[g] = section_int(&pointer, first) - section_int(&pointer, first - 1);
    }
    const char *names[2] = {"little.little_irr_full_dim", "little.little_irr_old"};
    for (int i = 0; i < 2; i++) {
        strncpy(D.entries[i].name, names[i], BLOB_NAME_SIZE - 1);
        D.entries[i].kind = BLOB_I64;
        D.entries[i].count = (uint64_t)count;
        D.entries[i].bytes = (uint64_t)count * sizeof(int64_t);
    }
    return 0;
}

/* section_find for little_irr_full_dim and little_irr_old; *handled 0 for any other name */
int little_derived_section(const char *name, section_view *view, int *handled)
{
    int which = !strncmp(name, "little.little_irr_full_dim", BLOB_NAME_SIZE) ? 0 : !strncmp(name, "little.little_irr_old", BLOB_NAME_SIZE) ? 1 : -1;
    *handled = which >= 0;
    if (!*handled) return 1;
    if (!D.ready && !D.failed) {
        D.failed = build();
        D.ready = !D.failed;
    }
    if (!D.ready) return 1;
    view->entry = &D.entries[which];
    view->data = (const unsigned char *)(which ? D.old : D.dims);
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
