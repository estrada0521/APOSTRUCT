/* The group of an isotropy row, computed from its irrep and order-parameter direction: the parent operations that fix
 * the direction, identified as a subgroup (a space group for the ordinary table, a magnetic group for the magnetic one). */

#include <stdlib.h>

#include "geometry.h"
#include "magnetic.h"
#include "irreps.h"
#include "tables.h"

int orderparam_to_subgroup_rows(int gid, int param_count, const int64_t *kparam, const double *values, int64_t value_count, int row_count,
                                int candidate_count, const int64_t *candidates, const int *candidate_matrix_ops, const int *candidate_signs,
                                int max_iter, int64_t *basis_out, int64_t *kept_out, int *kept_count);
int id_subgroup_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t *out);
int id_subgroup_magnetic_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t capacity,
                                  int64_t *out);
int isotropy_row_dim(int64_t row_id);
int isotropy_row_free(int64_t row_id);
int mag_isotropy_row_dim(int64_t row_id);
int mag_isotropy_row_free(int64_t row_id);

typedef struct {
    section_view irrep_pointer, op_pointer, op;
    int ready, magnetic, *subgroup;
    int64_t rows, olds;
} family;

static family ordinary, magnetic = {.magnetic = 1};
static section_view values, little_old;
static int *gid_of_old;

static int load(family *f)
{
    if (f->ready) return f->ready < 0;
    f->ready = -1;
    const char *names[2][3] = {{"isotropy.ordinary.by_old_start", "isotropy.ordinary.direction_codes_start", "isotropy.ordinary.direction_codes"},
                               {"isotropy.magnetic.by_old_start", "isotropy.magnetic.direction_codes_start", "isotropy.magnetic.direction_codes"}};
    if (geo_load() || (f->magnetic && mag_tables()) || section_find(names[f->magnetic][0], &f->irrep_pointer) ||
        section_find(names[f->magnetic][1], &f->op_pointer) || section_find(names[f->magnetic][2], &f->op) ||
        section_find("const.values", &values) ||
        section_find("irreps.little.old", &little_old))
        return 1;
    f->rows = section_count(&f->op_pointer);
    f->olds = section_count(&f->irrep_pointer);
    f->subgroup = calloc((size_t)f->rows + 1, sizeof(int));
    if (!f->subgroup) return 1;
    if (!gid_of_old) {
        int64_t olds = old_irrep_count();
        if (!(gid_of_old = calloc((size_t)olds + 1, sizeof(int)))) return 1;
        for (int64_t g = section_count(&little_old) - 1; g >= 0; g--) {
            int64_t old = section_int(&little_old, g);
            if (old > 0 && old <= olds) gid_of_old[old] = (int)g + 1;
        }
    }
    f->ready = 1;
    return 0;
}

static int64_t old_of_row(const family *f, int64_t row_id)
{
    int64_t lo = 1, hi = f->olds - 1, old = 0;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (section_int(&f->irrep_pointer, mid - 1) <= row_id) old = mid, lo = mid + 1;
        else hi = mid - 1;
    }
    return old;
}

static int compute(const family *f, int64_t row_id)
{
    int64_t old = old_of_row(f, row_id);
    int gid = old > 0 ? gid_of_old[old] : 0;
    int dim = f->magnetic ? mag_isotropy_row_dim(row_id) : isotropy_row_dim(row_id);
    int nf = f->magnetic ? mag_isotropy_row_free(row_id) : isotropy_row_free(row_id);
    if (!gid || dim <= 0 || dim > 48 || nf <= 0 || nf > 48) return -1;
    int parent = old_irrep_space_group((int)old);
    static double v[48 * 48];
    static int64_t candidates[192][5], kept[1536][5], out[8192];
    static int ops[192], signs[192], nc, cached = -1;
    int nk = 0;
    if (cached != parent * 2 + f->magnetic) {
        cached = -1;
        nc = f->magnetic ? magnetic_records_iso(mag_orderparam_group(parent), 192, candidates[0]) : geo_sg_records(parent, candidates);
        if (nc <= 0) return -1;
        for (int j = 0; j < nc; j++) {
            if (f->magnetic) ops[j] = mag_point_action(1, (int)candidates[j][4], &signs[j]);
            else ops[j] = (int)candidates[j][4], signs[j] = 1;
        }
        cached = parent * 2 + f->magnetic;
    }
    const double *constants = (const double *)values.data;
    int64_t p = section_int(&f->op_pointer, row_id - 1) - 1;
    for (int i = 0; i < 48 * 48; i++) v[i] = 0;
    for (int r = 0; r < nf; r++)
        for (int c = 0; c < dim; c++) v[r * 48 + c] = constants[section_int(&f->op, p + r * dim + c)];
    int64_t none = 0, basis[9] = {0};
    if (orderparam_to_subgroup_rows(gid, 0, &none, v, (int64_t)nf * 48, nf, nc, candidates[0], ops, signs, 256, basis, kept[0], &nk))
        return -1;
    if (f->magnetic) return id_subgroup_magnetic_identify(parent, basis, nk, kept[0], 0, 8192, out) > 0 ? (int)out[0] : -1;
    return id_subgroup_identify(parent, basis, nk, kept[0], 0, out) ? (int)out[0] : -1;
}

static int row_subgroup(family *f, int64_t row_id)
{
    if (load(f) || row_id < 1 || row_id > f->rows) return 0;
    if (!f->subgroup[row_id]) f->subgroup[row_id] = compute(f, row_id);
    return f->subgroup[row_id] > 0 ? f->subgroup[row_id] : 0;
}

/* the space group number of ordinary isotropy row row_id (1-based); 0 when it cannot be computed */
int isotropy_row_subgroup(int64_t row_id) { return row_subgroup(&ordinary, row_id); }
/* the magnetic group of magnetic isotropy row row_id */
int mag_isotropy_row_subgroup(int64_t row_id) { return row_subgroup(&magnetic, row_id); }
