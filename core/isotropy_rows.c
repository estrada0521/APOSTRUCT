/* Sizes of an isotropy row read off its order parameter: the irrep's dimension, and the number of direction rows
 * (free parameters) its codes span. Codes 0 only pad the end of the order-parameter array. */

#include "tables.h"

static struct {
    section_view irrep_pointer, op_pointer, op, little_old, full_dim;
    int ready;
    int64_t last_code;
} R;

static int load(void)
{
    if (R.ready) return 0;
    if (section_find("isotropy.isotropy_irrep_pointer", &R.irrep_pointer) ||
        section_find("isotropy.isotropy_orderparam_pointer", &R.op_pointer) || section_find("isotropy.isotropy_orderparam", &R.op) ||
        section_find("little.little_irr_old", &R.little_old) || section_find("little.little_irr_full_dim", &R.full_dim))
        return 1;
    R.last_code = section_count(&R.op);
    while (R.last_code > 0 && section_int(&R.op, R.last_code - 1) == 0) R.last_code--;
    R.ready = 1;
    return 0;
}

int64_t isotropy_row_count(void) { return load() ? 0 : section_count(&R.op_pointer); }

/* the dimension of the irrep of an isotropy row (1-based); 0 when unknown */
int isotropy_row_dim(int64_t row_id)
{
    if (load() || row_id < 1 || row_id > section_count(&R.op_pointer)) return 0;
    int64_t lo = 1, hi = section_count(&R.irrep_pointer) - 1, old = 0;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (section_int(&R.irrep_pointer, mid - 1) <= row_id) old = mid, lo = mid + 1;
        else hi = mid - 1;
    }
    for (int64_t g = 0; old && g < section_count(&R.little_old); g++)
        if (section_int(&R.little_old, g) == old) return (int)section_int(&R.full_dim, g);
    return 0;
}

/* the number of direction rows of an isotropy row */
int isotropy_row_free(int64_t row_id)
{
    int dim = isotropy_row_dim(row_id);
    if (!dim) return 0;
    int64_t start = section_int(&R.op_pointer, row_id - 1) - 1;
    int64_t end = row_id < section_count(&R.op_pointer) ? section_int(&R.op_pointer, row_id) - 1 : R.last_code;
    return (int)((end - start) / dim);
}
