/* Sizes of an isotropy row read off its order parameter: the irrep's dimension, and the number of direction rows
 * (free parameters) its codes span. Codes 0 only pad the end of the order-parameter array. The ordinary and the
 * magnetic isotropy tables share the layout (rows grouped by the ordinary old irrep). */

#include "tables.h"

typedef struct {
    section_view irrep_pointer, op_pointer, op, little_old, full_dim;
    int ready;
    int64_t last_code;
} family;

static family ordinary, magnetic;

static int load(family *f, const char *irrep_pointer, const char *op_pointer, const char *op)
{
    if (f->ready) return 0;
    if (section_find(irrep_pointer, &f->irrep_pointer) || section_find(op_pointer, &f->op_pointer) || section_find(op, &f->op) ||
        section_find("irreps.little.old", &f->little_old) || section_find("irreps.little.dim", &f->full_dim))
        return 1;
    f->last_code = section_count(&f->op);
    while (f->last_code > 0 && section_int(&f->op, f->last_code - 1) == 0) f->last_code--;
    f->ready = 1;
    return 0;
}

static family *of(int is_magnetic)
{
    if (is_magnetic)
        return load(&magnetic, "isotropy.magnetic.by_old_start", "isotropy.magnetic.direction_codes_start", "isotropy.magnetic.direction_codes")
                   ? 0
                   : &magnetic;
    return load(&ordinary, "isotropy.ordinary.by_old_start", "isotropy.ordinary.direction_codes_start", "isotropy.ordinary.direction_codes")
               ? 0
               : &ordinary;
}

static int row_dim(family *f, int64_t row_id)
{
    if (!f || row_id < 1 || row_id > section_count(&f->op_pointer)) return 0;
    int64_t lo = 1, hi = section_count(&f->irrep_pointer) - 1, old = 0;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (section_int(&f->irrep_pointer, mid - 1) <= row_id) old = mid, lo = mid + 1;
        else hi = mid - 1;
    }
    for (int64_t g = 0; old && g < section_count(&f->little_old); g++)
        if (section_int(&f->little_old, g) == old) return (int)section_int(&f->full_dim, g);
    return 0;
}

static int row_free(family *f, int64_t row_id)
{
    int dim = row_dim(f, row_id);
    if (!dim) return 0;
    int64_t start = section_int(&f->op_pointer, row_id - 1) - 1;
    int64_t end = row_id < section_count(&f->op_pointer) ? section_int(&f->op_pointer, row_id) - 1 : f->last_code;
    return (int)((end - start) / dim);
}

int64_t isotropy_row_count(void)
{
    family *f = of(0);
    return f ? section_count(&f->op_pointer) : 0;
}

/* the dimension of the irrep of an isotropy row (1-based); 0 when unknown */
int isotropy_row_dim(int64_t row_id) { return row_dim(of(0), row_id); }
/* the number of direction rows of an isotropy row */
int isotropy_row_free(int64_t row_id) { return row_free(of(0), row_id); }
int mag_isotropy_row_dim(int64_t row_id) { return row_dim(of(1), row_id); }
int mag_isotropy_row_free(int64_t row_id) { return row_free(of(1), row_id); }
