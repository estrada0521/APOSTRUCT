/* The OPD-table presentation of a subgroup basis (FINDSYM triclinic_lattparam_ and monoclinic_lattparam_,
 * the fixed-coupled trigonal basal turn), admitted for a magnetic subgroup
 * only while the selected BNS operation set stays fixed.  Floats keep the
 * operation order of the former Python kernels (Neumaier sums, radians/degrees, libm). */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"
#include "pynum.h"

enum {
    OB_RANGE = 1,
    OB_ZERO = 2,
    OB_TABLES = 3,
    OB_DECLINED = 4,
};

#define OB_ATOL 1e-5
#define OB_NCMP 1e-6
#define OB_SNAP 1.7e-5

static double clamp_unit(double x) { return (x < 1.0 ? x : 1.0) > -1.0 ? (x < 1.0 ? x : 1.0) : -1.0; }
static int ncmp(double x) { return fabs(x) > OB_NCMP; }
static int64_t nint(double x) { return x >= 0.0 ? (int64_t)floor(x + 0.5) : (int64_t)ceil(x - 0.5); }

static void get_gmat(const double *cell, double g[3][3])
{
    double a = cell[0], b = cell[1], c = cell[2];
    g[0][0] = a * a;
    g[0][1] = a * b * geo_cos_deg(cell[5]);
    g[0][2] = a * c * geo_cos_deg(cell[4]);
    g[1][0] = a * b * geo_cos_deg(cell[5]);
    g[1][1] = b * b;
    g[1][2] = b * c * geo_cos_deg(cell[3]);
    g[2][0] = a * c * geo_cos_deg(cell[4]);
    g[2][1] = b * c * geo_cos_deg(cell[3]);
    g[2][2] = c * c;
}

static void congruence(double t[3][3], double g[3][3], double out[3][3])
{
    double r[3][3], terms[9];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            int n = 0;
            for (int p = 0; p < 3; p++)
                for (int q = 0; q < 3; q++) terms[n++] = t[row][p] * g[p][q] * t[col][q];
            r[row][col] = py_sum(terms, 9);
        }
    memcpy(out, r, sizeof r);
}

static void float_matmul(double l[3][3], double r[3][3], double out[3][3])
{
    double o[3][3], terms[3];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            for (int k = 0; k < 3; k++) terms[k] = l[row][k] * r[k][col];
            o[row][col] = py_sum(terms, 3);
        }
    memcpy(out, o, sizeof o);
}

static void gram_lattparam(double g[3][3], double *p)
{
    double a = sqrt(g[0][0]), b = sqrt(g[1][1]), c = sqrt(g[2][2]);
    p[0] = a;
    p[1] = b;
    p[2] = c;
    p[3] = geo_acos_upstream(clamp_unit(g[1][2] / b / c));
    p[4] = geo_acos_upstream(clamp_unit(g[0][2] / a / c));
    p[5] = geo_acos_upstream(clamp_unit(g[0][1] / a / b));
}

static void cell_cartesian(const double *cell, double v[3][3])
{
    double a = cell[0], b = cell[1], c = cell[2];
    double ca = geo_cos_deg(cell[3]), cb = geo_cos_deg(cell[4]), cg = geo_cos_deg(cell[5]), sg = geo_sin_deg(cell[5]);
    double cx = c * cb, cy = c * (ca - cb * cg) / sg;
    double rest = c * c - cx * cx - cy * cy;
    double cz = sqrt(0.0 < rest ? rest : 0.0);
    double rows[3][3] = {{a, 0.0, 0.0}, {b * cg, b * sg, 0.0}, {cx, cy, cz}};
    memcpy(v, rows, sizeof rows);
}

/* ---- the magnetic BNS operation set ---- */

typedef struct {
    rat rotation[9], translation[3];
    int antiunitary;
} mag_op;

static int op_equal(const mag_op *a, const mag_op *b)
{
    if (a->antiunitary != b->antiunitary) return 0;
    for (int i = 0; i < 9; i++)
        if (!rat_eq(a->rotation[i], b->rotation[i])) return 0;
    for (int i = 0; i < 3; i++)
        if (!rat_eq(a->translation[i], b->translation[i])) return 0;
    return 1;
}

static int add_unique(mag_op *set, int *count, const mag_op *op)
{
    for (int i = 0; i < *count; i++)
        if (op_equal(set + i, op)) return 0;
    set[(*count)++] = *op;
    return 1;
}

static int magnetic_set(int group, mag_op *set, int *count)
{
    int64_t records[192][5];
    int n = magnetic_records(group, 2, 192, &records[0][0]);
    if (n < 0) return OB_TABLES;
    *count = 0;
    for (int r = 0; r < n; r++) {
        mag_op op;
        int sign, point_op = mag_point_action(1, (int)records[r][4], &sign);
        int64_t m[9];
        geo_point_op_matrix(point_op, m);
        for (int i = 0; i < 9; i++) op.rotation[i] = rat_int(m[i]);
        for (int axis = 0; axis < 3; axis++) op.translation[axis] = rat_mod1(rat_make(records[r][axis], records[r][3]));
        op.antiunitary = sign < 0;
        add_unique(set, count, &op);
    }
    return mag_take_error() || geo_take_error() ? OB_TABLES : 0;
}

/* _magnetic_basis_candidate_preserves_group */
static int preserves_group(int group, rat source[3][3], rat candidate[3][3], int *preserved)
{
    static mag_op set[192], moved[192];
    int count, moved_count = 0;
    if (magnetic_set(group, set, &count)) return OB_TABLES;
    rat source_inverse[3][3], transform[3][3], inverse[3][3];
    if (rat_inverse3((const rat (*)[3])source, source_inverse)) return OB_ZERO;
    for (int row = 0; row < 3; row++) rat_row_mul(candidate[row], (const rat (*)[3])source_inverse, transform[row]);
    if (rat_inverse3((const rat (*)[3])transform, inverse)) return OB_ZERO;
    for (int i = 0; i < count; i++) {
        rat r[3][3], left[3][3], out[3][3];
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) r[row][col] = set[i].rotation[3 * row + col];
        for (int row = 0; row < 3; row++) rat_row_mul(transform[row], (const rat (*)[3])r, left[row]);
        for (int row = 0; row < 3; row++) rat_row_mul(left[row], (const rat (*)[3])inverse, out[row]);
        mag_op op;
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) op.rotation[3 * row + col] = out[row][col];
        rat t[3];
        rat_row_mul(set[i].translation, (const rat (*)[3])inverse, t);
        for (int axis = 0; axis < 3; axis++) op.translation[axis] = rat_mod1(t[axis]);
        op.antiunitary = set[i].antiunitary;
        add_unique(moved, &moved_count, &op);
    }
    int same = moved_count == count;
    for (int i = 0; i < moved_count && same; i++) {
        int found = 0;
        for (int j = 0; j < count && !found; j++) found = op_equal(moved + i, set + j);
        same = found;
    }
    *preserved = same;
    return rat_overflow ? OB_RANGE : 0;
}

static void apply_int(const int64_t t[3][3], rat rows[3][3], rat out[3][3])
{
    rat o[3][3];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            rat total = rat_int(0);
            for (int k = 0; k < 3; k++) total = rat_add(total, rat_mul(rat_int(t[row][k]), rows[k][col]));
            o[row][col] = total;
        }
    memcpy(out, o, sizeof o);
}

typedef struct {
    int group;
    rat (*rows)[3];
} filter_ctx;

/* transform_filter: whether the cumulative integer transform keeps the group. */
static int admits(const filter_ctx *f, const int64_t t[3][3], int *admitted)
{
    if (!f->group) {
        *admitted = 1;
        return 0;
    }
    rat candidate[3][3];
    apply_int(t, f->rows, candidate);
    return preserves_group(f->group, f->rows, candidate, admitted);
}

/* ---- triclinic_lattparam_ ---- */

static int lattparam_sort(double g[3][3], double perm[3][3])
{
    int used[3] = {0, 0, 0}, picked = 0, changed = 0;
    double maxdiag = g[0][0];
    if (g[1][1] > maxdiag) maxdiag = g[1][1];
    if (g[2][2] > maxdiag) maxdiag = g[2][2];
    memset(perm, 0, sizeof(double) * 9);
    for (int position = 0; position < 3; position++) {
        double best = maxdiag + 1.0;
        for (int index = 0; index < 3; index++) {
            if (used[index]) continue;
            double value = g[index][index];
            if (value < best && ncmp(value - best)) {
                picked = index;
                best = value;
            }
        }
        used[picked] = 1;
        perm[position][picked] = 1.0;
        if (picked != position) changed = 1;
    }
    double det = perm[0][0] * (perm[1][1] * perm[2][2] - perm[1][2] * perm[2][1]) -
                 perm[0][1] * (perm[1][0] * perm[2][2] - perm[1][2] * perm[2][0]) +
                 perm[0][2] * (perm[1][0] * perm[2][1] - perm[1][1] * perm[2][0]);
    if (det < 0)
        for (int col = 0; col < 3; col++) perm[2][col] = -perm[2][col];
    return changed;
}

static int lattparam_sign(double g[3][3], int zero_default, double flips[3][3])
{
    static const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    for (int i = 0; i < 3; i++) {
        int row = pairs[i][0], col = pairs[i][1];
        if (fabs(g[row][col]) < OB_SNAP) g[row][col] = g[col][row] = 0.0;
    }
    int signs[3];
    for (int i = 0; i < 3; i++) {
        double value = g[pairs[i][0]][pairs[i][1]];
        signs[i] = !ncmp(value) ? zero_default : value > 0 ? 1 : -1;
    }
    int zeta = signs[0], eta = signs[1], xi = signs[2], total = zeta + eta + xi;
    if (abs(total) == 3) return 0;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) flips[row][col] = row == col ? 1.0 : 0.0;
    if (total == xi) flips[0][0] = -1.0;
    if (total == eta) flips[1][1] = -1.0;
    if (total == zeta) flips[2][2] = -1.0;
    return 1;
}

typedef struct {
    double total[3][3];
    const filter_ctx *filter;
} reducer;

static int reduce_apply(reducer *r, double step[3][3], double g[3][3], int *accepted)
{
    double candidate[3][3];
    float_matmul(step, r->total, candidate);
    *accepted = 1;
    if (r->filter->group) {
        int64_t exact[3][3];
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) exact[row][col] = nint(candidate[row][col]);
        int status = admits(r->filter, (const int64_t (*)[3])exact, accepted);
        if (status) return status;
        if (!*accepted) return 0;
    }
    memcpy(r->total, candidate, sizeof candidate);
    congruence(step, g, g);
    return 0;
}

static int sort_and_sign(reducer *r, double g[3][3], int zero_default)
{
    double step[3][3];
    int accepted, status;
    if (lattparam_sort(g, step) && (status = reduce_apply(r, step, g, &accepted))) return status;
    if (lattparam_sign(g, zero_default, step) && (status = reduce_apply(r, step, g, &accepted))) return status;
    return 0;
}

static int triclinic_reduce(double g[3][3], const filter_ctx *filter, double total[3][3])
{
    reducer r = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, filter};
    double previous[6], params[6];
    int status;
    gram_lattparam(g, previous);
    for (int outer = 0; outer < 10; outer++) {
        int order[3] = {0, 1, 2}, failures = 0;
        for (;;) {
            int pivot = order[0], second = order[1], third = order[2];
            if (g[pivot][pivot] == 0.0) return OB_ZERO;
            int64_t shear_second = nint(-(g[second][pivot] / g[pivot][pivot]));
            int64_t shear_third = nint(-(g[third][pivot] / g[pivot][pivot]));
            if (shear_second == 0 && shear_third == 0) {
                failures++;
            } else {
                double step[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, candidate[3][3];
                step[second][pivot] = (double)shear_second;
                step[third][pivot] = (double)shear_third;
                congruence(step, g, candidate);
                int changes = 0;
                for (int d = 0; d < 3; d++) changes |= ncmp(candidate[d][d] - g[d][d]);
                if (changes) {
                    int accepted;
                    if ((status = reduce_apply(&r, step, g, &accepted))) return status;
                    if (!accepted) failures++;
                } else {
                    failures++;
                }
            }
            if (failures > 2) break;
            for (int i = 0; i < 3; i++) order[i] = (order[i] + 1) % 3;
        }
        if ((status = sort_and_sign(&r, g, -1))) return status;
        double product = g[0][1] * g[0][2] * g[1][2];
        double slack = g[0][0] + g[1][1] - fabs(2.0 * g[0][1]) - fabs(2.0 * g[0][2]) - fabs(2.0 * g[1][2]);
        if (product < 0.0 && slack < 0.0 && ncmp(product) && ncmp(slack)) {
            double step[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            step[2][0] = g[0][2] > 0.0 ? -1.0 : 1.0;
            step[2][1] = g[1][2] > 0.0 ? -1.0 : 1.0;
            int accepted;
            if ((status = reduce_apply(&r, step, g, &accepted))) return status;
            if ((status = sort_and_sign(&r, g, 1))) return status;
        }
        gram_lattparam(g, params);
        int fixed = 1;
        for (int i = 0; i < 6; i++) fixed &= !ncmp(params[i] - previous[i]);
        if (fixed) break;
        memcpy(previous, params, sizeof params);
    }
    memcpy(total, r.total, sizeof r.total);
    return 0;
}

static int triclinic_candidate(int subgroup, rat rows[3][3], const double *cell, int group, rat out[3][3])
{
    memcpy(out, rows, sizeof(rat) * 9);
    if (subgroup > 2 || !cell) return 0;
    double parent[3][3], basis[3][3], child[3][3], g[3][3], terms[3], lattice[6];
    cell_cartesian(cell, parent);
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) basis[row][col] = rat_to_double(rows[row][col]);
    float_matmul(basis, parent, child);
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            for (int axis = 0; axis < 3; axis++) terms[axis] = child[row][axis] * child[col][axis];
            g[row][col] = py_sum(terms, 3);
        }
    gram_lattparam(g, lattice);
    for (int i = 3; i < 6; i++)
        if (fabs(lattice[i] - 90.0) < 0.001) lattice[i] = 90.0;
    get_gmat(lattice, g);
    filter_ctx filter = {group, rows};
    double total[3][3], reduced[3][3], start[3][3], final_params[6];
    int status = triclinic_reduce(g, &filter, total);
    if (status) return status;
    get_gmat(lattice, start);
    congruence(total, start, reduced);
    gram_lattparam(reduced, final_params);
    int unchanged = 1;
    for (int i = 0; i < 6; i++) unchanged &= !ncmp(final_params[i] - lattice[i]);
    if (unchanged) return 0;
    int64_t step[3][3];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) step[row][col] = nint(total[row][col]);
    apply_int((const int64_t (*)[3])step, rows, out);
    return 0;
}

/* ---- monoclinic_lattparam_ ---- */

/* vector_change_setting_record: the reduced (x, y, z, den) of a record
 * moved between the cinter, cml and pml settings (codes 0, 1, 2). */
static int to_cinter(int sg, int setting, rat m[3][3])
{
    rat origin[3];
    if (!setting) {
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) m[row][col] = rat_int(row == col);
        return 0;
    }
    return geo_setting_to_cinter(sg, setting == 2, 0, m, origin) ? OB_TABLES : 0;
}

static int change_setting(int sg, int from, int to, const int64_t *record, int64_t *out)
{
    rat a[3][3], b[3][3], inverse[3][3], m[3][3];
    if (to_cinter(sg, from, a) || to_cinter(sg, to, b)) return OB_TABLES;
    if (rat_inverse3((const rat (*)[3])b, inverse)) return OB_ZERO;
    for (int row = 0; row < 3; row++) rat_row_mul(a[row], (const rat (*)[3])inverse, m[row]);
    if (!record[3]) return OB_ZERO;
    rat v[3] = {rat_make(record[0], record[3]), rat_make(record[1], record[3]), rat_make(record[2], record[3])}, w[3];
    rat_row_mul(v, (const rat (*)[3])m, w);
    i128 den = 1;
    for (int axis = 0; axis < 3; axis++) den = rat_mul_i(den / rat_gcd(den, w[axis].d), w[axis].d);
    i128 n[4], g = 0;
    for (int axis = 0; axis < 3; axis++) n[axis] = rat_mul_i(w[axis].n, den / w[axis].d);
    n[3] = den;
    for (int i = 0; i < 4; i++)
        if (n[i]) g = g ? rat_gcd(g, n[i]) : (n[i] < 0 ? -n[i] : n[i]);
    if (!g) g = 1;
    for (int i = 0; i < 4; i++) {
        i128 value = n[i] / g;
        if (value > INT64_MAX || value < INT64_MIN) return OB_RANGE;
        out[i] = (int64_t)value;
    }
    return 0;
}

static void vmlt(const int64_t *m, const int64_t *v, int64_t *out)
{
    int64_t x = v[0], y = v[1], z = v[2], den = v[3];
    out[0] = x * m[0] + y * m[3] + z * m[6];
    out[1] = x * m[1] + y * m[4] + z * m[7];
    out[2] = x * m[2] + y * m[5] + z * m[8];
    out[3] = den;
}

static int64_t gcd64(int64_t a, int64_t b) { return (int64_t)gcd_u64((uint64_t)(a < 0 ? -a : a), (uint64_t)(b < 0 ? -b : b)); }

static int preserves_translations(int sg, const int64_t t[3][3], int *ok)
{
    *ok = 0;
    int64_t flat[9], adj[9];
    for (int i = 0; i < 9; i++) flat[i] = t[i / 3][i % 3];
    adj[0] = flat[4] * flat[8] - flat[5] * flat[7];
    adj[1] = flat[2] * flat[7] - flat[1] * flat[8];
    adj[2] = flat[1] * flat[5] - flat[2] * flat[4];
    adj[3] = flat[5] * flat[6] - flat[3] * flat[8];
    adj[4] = flat[0] * flat[8] - flat[2] * flat[6];
    adj[5] = flat[2] * flat[3] - flat[0] * flat[5];
    adj[6] = flat[3] * flat[7] - flat[4] * flat[6];
    adj[7] = flat[1] * flat[6] - flat[0] * flat[7];
    adj[8] = flat[0] * flat[4] - flat[1] * flat[3];
    int64_t det = flat[0] * adj[0] + flat[1] * adj[3] + flat[2] * adj[6];
    if (det != 1 && det != -1) return OB_DECLINED;
    int64_t inverse[9], pml[9];
    for (int i = 0; i < 9; i++) inverse[i] = adj[i] / det;
    for (int unit = 0; unit < 3; unit++) {
        int64_t vector[4] = {unit == 0, unit == 1, unit == 2, 1}, a[4], b[4], c[4];
        int status = change_setting(sg, 2, 0, vector, a);
        if (status) return status;
        vmlt(inverse, a, b);
        if ((status = change_setting(sg, 0, 2, b, c))) return status;
        if (c[3] != 1) return 0;
        for (int k = 0; k < 3; k++) pml[3 * unit + k] = c[k];
    }
    int64_t records[GEO_MAX_OPS][5];
    int count = geo_sg_records(sg, records);
    if (geo_take_error()) return OB_TABLES;
    for (int r = 1; r < count; r++) {
        int64_t moved[4];
        vmlt(pml, records[r], moved);
        int64_t den = records[r][3] / gcd64(records[r][3], moved[3]) * moved[3];
        for (int axis = 0; axis < 3; axis++) {
            int64_t delta = records[r][axis] * (den / records[r][3]) - moved[axis] * (den / moved[3]);
            if (delta % den) return 0;
        }
    }
    *ok = 1;
    return 0;
}

static void matmlt(const int64_t l[3][3], const int64_t r[3][3], int64_t out[3][3])
{
    const int64_t *lf = &l[0][0], *rf = &r[0][0];
    int64_t o[9];
    for (int outer = 0; outer < 3; outer++)
        for (int group = 0; group < 3; group++) {
            int64_t total = 0;
            for (int k = 0; k < 3; k++) total += lf[outer + 3 * k] * rf[3 * group + k];
            o[outer + 3 * group] = total;
        }
    memcpy(out, o, sizeof o);
}

static void transform_gram(const int64_t t[3][3], double g[3][3], double out[3][3])
{
    double r[3][3], terms[9];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            int n = 0;
            for (int p = 0; p < 3; p++)
                for (int q = 0; q < 3; q++) terms[n++] = (double)t[row][p] * g[p][q] * (double)t[col][q];
            r[row][col] = py_sum(terms, 9);
        }
    memcpy(out, r, sizeof r);
}

typedef struct {
    int sg;
    int64_t transform[3][3];
    const filter_ctx *filter;
} mono;

static int mono_accept(mono *m, const int64_t candidate[3][3], int *accepted)
{
    int status = preserves_translations(m->sg, candidate, accepted);
    if (status || !*accepted) return status;
    int64_t combined[3][3];
    matmlt(m->transform, candidate, combined);
    return admits(m->filter, (const int64_t (*)[3])combined, accepted);
}

static void mono_accumulate(mono *m, const int64_t candidate[3][3]) { matmlt(m->transform, candidate, m->transform); }

static int monoclinic_transform(int sg, const double *lattice, const filter_ctx *filter, int64_t out[3][3])
{
    double g[3][3];
    get_gmat(lattice, g);
    mono m = {sg, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, filter};
    int64_t axis_vector[4] = {0, 0, 1, 1}, v[4];
    int status = change_setting(sg, 1, 0, axis_vector, v);
    if (status) return status;
    int unique, axis_a, axis_c;
    if (v[0]) unique = 0, axis_a = 1, axis_c = 2;
    else if (v[1]) unique = 1, axis_a = 0, axis_c = 2;
    else if (v[2]) unique = 2, axis_a = 0, axis_c = 1;
    else return OB_DECLINED;
    int pair[2] = {axis_a, axis_c};
    for (int pass = 0; pass < 2; pass++) {
        int other = pair[0], target = pair[1];
        int64_t step = g[target][other] > 0 ? -1 : 1, offset = 0;
        for (;;) {
            offset += step;
            if (offset > 1000000 || offset < -1000000) return OB_DECLINED;
            int64_t candidate[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            candidate[target][other] = offset;
            int accepted;
            if ((status = mono_accept(&m, (const int64_t (*)[3])candidate, &accepted))) return status;
            if (!accepted) continue;
            double next[3][3];
            transform_gram((const int64_t (*)[3])candidate, g, next);
            if (next[target][target] < g[target][target] - OB_ATOL) {
                memcpy(g, next, sizeof next);
                mono_accumulate(&m, (const int64_t (*)[3])candidate);
                offset = 0;
                continue;
            }
            break;
        }
        int swap = pair[0];
        pair[0] = pair[1];
        pair[1] = swap;
    }
    int low = axis_a < axis_c ? axis_a : axis_c, high = axis_a < axis_c ? axis_c : axis_a, accepted;
    if (g[high][high] < g[low][low] - OB_ATOL) {
        int64_t candidate[3][3] = {{0}};
        candidate[unique][unique] = 1;
        candidate[low][high] = -1;
        candidate[high][low] = 1;
        if ((status = mono_accept(&m, (const int64_t (*)[3])candidate, &accepted))) return status;
        if (accepted) {
            double next[3][3];
            transform_gram((const int64_t (*)[3])candidate, g, next);
            memcpy(g, next, sizeof next);
            mono_accumulate(&m, (const int64_t (*)[3])candidate);
        }
    }
    if (g[high][low] > OB_ATOL) {
        int64_t candidate[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        candidate[unique][unique] = -1;
        candidate[high][high] = -1;
        if ((status = mono_accept(&m, (const int64_t (*)[3])candidate, &accepted))) return status;
        if (accepted) mono_accumulate(&m, (const int64_t (*)[3])candidate);
    }
    memcpy(out, m.transform, sizeof m.transform);
    return 0;
}

static double norm3(const double *v)
{
    double terms[3] = {v[0] * v[0], v[1] * v[1], v[2] * v[2]};
    return sqrt(py_sum(terms, 3));
}

static double angle3(const double *a, const double *b)
{
    double denom = norm3(a) * norm3(b);
    if (denom <= 0) return 90.0;
    double terms[3] = {a[0] * b[0], a[1] * b[1], a[2] * b[2]};
    return geo_acos_upstream(clamp_unit(py_sum(terms, 3) / denom));
}

static int monoclinic_candidate(int subgroup, rat rows[3][3], const double *cell, int group, rat out[3][3])
{
    memcpy(out, rows, sizeof(rat) * 9);
    if (subgroup < 3 || subgroup > 15 || !cell) return 0;
    double parent[3][3], vectors[3][3], terms[3], lattice[6];
    cell_cartesian(cell, parent);
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            for (int idx = 0; idx < 3; idx++) terms[idx] = rat_to_double(rows[row][idx]) * parent[idx][col];
            vectors[row][col] = py_sum(terms, 3);
        }
    lattice[0] = norm3(vectors[0]);
    lattice[1] = norm3(vectors[1]);
    lattice[2] = norm3(vectors[2]);
    lattice[3] = angle3(vectors[1], vectors[2]);
    lattice[4] = angle3(vectors[0], vectors[2]);
    lattice[5] = angle3(vectors[0], vectors[1]);
    filter_ctx filter = {group, rows};
    int64_t transform[3][3];
    int prior = rat_overflow;
    rat_overflow = 0;
    int failed = monoclinic_transform(subgroup, lattice, &filter, transform) || rat_overflow;
    rat_overflow = prior;
    if (!failed) apply_int((const int64_t (*)[3])transform, rows, out);
    return 0;
}

/* basis: three rows of rationals (int128 numerator, denominator words);
 * cell: the parent lattice parameters or 0; magnetic_group 0 for an
 * ordinary subgroup.  out: the presented rows in the same words. */
int opd_basis_presentation(int subgroup, const int64_t *basis, const double *cell, int parametric, int coupled,
                           int magnetic_group, int64_t *out)
{
    rat source[3][3], rows[3][3];
    if (geo_load()) return OB_TABLES;
    rat_overflow = 0;
    for (int row = 0; row < 3; row++) read_rat3(basis + 12 * row, source[row]);
    if (rat_overflow) return OB_RANGE;
    if (magnetic_group && mag_tables()) return OB_TABLES;
    int status = triclinic_candidate(subgroup, source, cell, magnetic_group, rows);
    if (status) return status;
    if (coupled && !parametric && subgroup == 160) {
        static const int64_t turn[3][3] = {{0, 1, 0}, {-1, -1, 0}, {0, 0, 1}};
        apply_int(turn, rows, rows);
    }
    rat presented[3][3];
    if ((status = monoclinic_candidate(subgroup, rows, cell, magnetic_group, presented))) return status;
    if (magnetic_group) {
        int preserved;
        if ((status = preserves_group(magnetic_group, source, presented, &preserved))) return status;
        if (!preserved) memcpy(presented, source, sizeof source);
    }
    for (int row = 0; row < 3; row++) write_rat3(out + 12 * row, presented[row]);
    return rat_overflow ? OB_RANGE : 0;
}
