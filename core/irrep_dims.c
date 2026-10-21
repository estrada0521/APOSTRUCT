/* little_irr_full_dim built from the stored little-irrep matrices: an irrep's dimension is the number of entries of its
 * identity matrix (the first stored position), and 0 for an irrep without matrices. */

#include <stdlib.h>
#include <string.h>

#include "tables.h"

int section_find_stored(const char *name, section_view *view);

static struct {
    int ready, failed;
    int64_t *dims;
    blob_section entry;
} D;

static int build(void)
{
    section_view sg, irr, position, pointer;
    if (section_find_stored("little.little_irr_space_group", &sg) || section_find_stored("little.little_irr_full_matrices_irr_pointer", &irr) ||
        section_find_stored("little.little_irr_full_matrices_position", &position) ||
        section_find_stored("little.little_irr_full_matrices_pointer", &pointer))
        return 1;
    int64_t count = section_count(&sg), irrs = section_count(&irr), positions = section_count(&position);
    D.dims = calloc((size_t)count + 1, sizeof(int64_t));
    if (!D.dims) return 1;
    for (int64_t g = 0; g + 1 < irrs && g < count; g++) {
        int64_t first = section_int(&irr, g), end = section_int(&irr, g + 1);
        if (end <= first || first < 1 || first > positions || section_int(&position, first - 1) != 1) continue;
        D.dims[g] = section_int(&pointer, first) - section_int(&pointer, first - 1);
    }
    strncpy(D.entry.name, "little.little_irr_full_dim", BLOB_NAME_SIZE - 1);
    D.entry.kind = BLOB_I64;
    D.entry.count = (uint64_t)count;
    D.entry.bytes = (uint64_t)count * sizeof(int64_t);
    return 0;
}

/* section_find for little_irr_full_dim; *handled 0 for any other name */
int irrep_dims_section(const char *name, section_view *view, int *handled)
{
    *handled = !strncmp(name, "little.little_irr_full_dim", BLOB_NAME_SIZE);
    if (!*handled) return 1;
    if (!D.ready && !D.failed) {
        D.failed = build();
        D.ready = !D.failed;
    }
    if (!D.ready) return 1;
    view->entry = &D.entry;
    view->data = (const unsigned char *)D.dims;
    return 0;
}

/* the dimension of little irrep gid, 0 when unknown */
int little_full_dim(int gid)
{
    section_view view;
    if (section_find("little.little_irr_full_dim", &view) || gid < 1 || gid > section_count(&view)) return 0;
    return (int)section_int(&view, gid - 1);
}
