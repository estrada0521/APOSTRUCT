/* Source magnetic space groups (data_magnetic): the generated operation
 * records in their coordinate settings, magnetic Wyckoff rows and branches
 * (get_wyckoff_all_magnetic_), point identification (id_wyckoff2_magnetic_),
 * centered orbits, and the magnetic regrouping of a child atom population. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"
#include "tables.h"

enum {
    MG_OK = 0,
    MG_TABLES = 1,
    MG_RANGE = 2,
    MG_MEMORY = 3,
    MG_POINT_OP = 4,
    MG_DEPENDENT = 5,
    MG_NO_FORMULA = 6,
    MG_SHARED_POINT = 7,
    MG_PARTITION = 8,
    MG_LEAVES = 9,
    MG_PERMUTATION = 10,
    MG_SPLIT_SITE = 11,
    MG_ORBIT = 12,
    MG_REPRESENTATIVE = 13,
    MG_SINGULAR = 14,
    MG_CAPACITY = 15,
};

enum {
    SET_CLOSURE = 0,
    SET_BINARY = 1,
    SET_CINTER = 2,
};

#define MG_MAX_OPS 192
#define MG_MAX_ROWS 32

static struct {
    section_view m2n, r, nlabel, type, sg, t4_pointer, t4_sg, t4_basis, t4_origin, gen_count, gen, op_code,
        point_group, point_group_code;
    section_view psettings, point_op, lattice, label_bc, sg_pointer;
    int ready;
    short products[144 * 144];
} M;

static int mag_bad;

static void fail(int status)
{
    if (!mag_bad) mag_bad = status;
}

static int load(void)
{
    if (M.ready) return 0;
    if (geo_load()) return MG_TABLES;
    struct {
        const char *name;
        section_view *view;
        int text;
    } items[] = {
        {"magnetic.mag_point_op_mag2nonmag", &M.m2n, 0},
        {"magnetic.mag_point_op_r", &M.r, 0},
        {"magnetic.mag_nlabel", &M.nlabel, 1},
        {"magnetic.mag_type", &M.type, 0},
        {"magnetic.mag_space_group", &M.sg, 0},
        {"magnetic.mag_type4_pointer", &M.t4_pointer, 0},
        {"magnetic.mag_type4_supergroup", &M.t4_sg, 0},
        {"magnetic.mag_type4_supergroup_basis", &M.t4_basis, 0},
        {"magnetic.mag_type4_supergroup_origin", &M.t4_origin, 0},
        {"magnetic.mag_gen_count", &M.gen_count, 0},
        {"magnetic.mag_gen", &M.gen, 0},
        {"magnetic.mag_point_op_code", &M.op_code, 0},
        {"magnetic.mag_point_group", &M.point_group, 0},
        {"magnetic.mag_point_group_code", &M.point_group_code, 0},
        {"space.ipoint_op_psettings", &M.psettings, 0},
        {"space.ipoint_op", &M.point_op, 0},
        {"space.ispace_lattice", &M.lattice, 0},
        {"space.space_label_bc", &M.label_bc, 1},
        {"magnetic.mag_space_group_pointer", &M.sg_pointer, 0},
    };
    for (size_t i = 0; i < sizeof items / sizeof items[0]; i++) {
        if (section_find(items[i].name, items[i].view)) return MG_TABLES;
        if (!items[i].text && !section_is_int(items[i].view)) return MG_TABLES;
    }
    int64_t count = section_count(&M.m2n);
    if (count != 144 || section_count(&M.r) != 144) return MG_TABLES;
    short of[73][2] = {{0}};
    for (int i = 0; i < 144; i++) {
        int64_t ordinary = section_int(&M.m2n, i);
        if (ordinary < 1 || ordinary > 72) return MG_TABLES;
        of[ordinary][section_int(&M.r, i) != 0] = (short)(i + 1);
    }
    /* products by the ordinary parts and the time reversals */
    for (int left = 0; left < 144; left++)
        for (int right = 0; right < 144; right++) {
            int ordinary = (int)geo_point_product((int)section_int(&M.m2n, right), (int)section_int(&M.m2n, left));
            int reversed = (section_int(&M.r, left) != 0) != (section_int(&M.r, right) != 0);
            M.products[right * 144 + left] = ordinary ? of[ordinary][reversed] : 0;
        }
    M.ready = 1;
    return 0;
}

static int64_t at(const section_view *view, int64_t index)
{
    if (index < 0 || index >= section_count(view)) {
        fail(MG_TABLES);
        return 0;
    }
    return section_int(view, index);
}

static void values_of(const int64_t *record, rat *out)
{
    for (int axis = 0; axis < 3; axis++) out[axis] = record[3] ? rat_make(record[axis], record[3]) : rat_int(0);
}

static void fold3(rat *v)
{
    for (int axis = 0; axis < 3; axis++) v[axis] = rat_mod1(v[axis]);
}

static int eq3(const rat *a, const rat *b) { return rat_eq(a[0], b[0]) && rat_eq(a[1], b[1]) && rat_eq(a[2], b[2]); }
static int zero3(const rat *a) { return rat_zero(a[0]) && rat_zero(a[1]) && rat_zero(a[2]); }

/* fraction_record3 */
static void record_of(const rat *v, int64_t *out)
{
    i128 den = 1;
    for (int axis = 0; axis < 3; axis++) den = rat_mul_i(den / rat_gcd(den, v[axis].d), v[axis].d);
    i128 n[3], g = den;
    for (int axis = 0; axis < 3; axis++) {
        n[axis] = rat_mul_i(v[axis].n, den / v[axis].d);
        g = rat_gcd(g, n[axis]);
    }
    if (g > 1) {
        for (int axis = 0; axis < 3; axis++) n[axis] /= g;
        den /= g;
    }
    for (int axis = 0; axis < 3; axis++) {
        if (n[axis] > INT64_MAX || n[axis] < INT64_MIN) rat_overflow = 1;
        out[axis] = (int64_t)n[axis];
    }
    if (den > INT64_MAX) rat_overflow = 1;
    out[3] = (int64_t)den;
}

static void mat_mul(const rat (*a)[3], const rat (*b)[3], rat (*out)[3])
{
    rat t[3][3];
    for (int row = 0; row < 3; row++) rat_row_mul(a[row], b, t[row]);
    memcpy(out, t, sizeof t);
}

typedef struct {
    int group, type, sg, reference_sg;
    int64_t basis[9], origin[4];
    rat to_cinter[3][3], to_pml[3][3], origin_cinter[3];
    int centering_count;
    rat centering[4][3];
} group_setting;

static int64_t centering_count_of(int sg)
{
    const char *text;
    int64_t length;
    if (section_string(&M.label_bc, sg - 1, &text, &length) || !length) return 1;
    switch (text[0]) {
    case 'A': case 'B': case 'C': case 'I': return 2;
    case 'F': return 4;
    case 'R': return 3;
    default: return 1;
    }
}

static int load_setting(int group, group_setting *s)
{
    memset(s, 0, sizeof *s);
    if (group < 1 || group > section_count(&M.type)) return MG_RANGE;
    s->group = group;
    s->type = (int)at(&M.type, group - 1);
    s->sg = (int)at(&M.sg, group - 1);
    s->reference_sg = s->sg;
    for (int i = 0; i < 9; i++) s->basis[i] = i % 4 == 0;
    s->origin[3] = 1;
    if (s->type == 4) {
        int64_t pointer = at(&M.t4_pointer, group - 1);
        if (pointer <= 0) return MG_TABLES;
        s->reference_sg = (int)at(&M.t4_sg, pointer - 1);
        for (int i = 0; i < 9; i++) s->basis[i] = at(&M.t4_basis, (pointer - 1) * 9 + i);
        for (int i = 0; i < 4; i++) s->origin[i] = at(&M.t4_origin, (pointer - 1) * 4 + i);
    }
    if (geo_setting_to_cinter(s->sg, 1, 0, s->to_cinter, s->origin_cinter)) return MG_TABLES;
    if (rat_inverse3((const rat (*)[3])s->to_cinter, s->to_pml)) return MG_SINGULAR;
    int64_t count = centering_count_of(s->sg);
    for (int x = 0; x < count && s->centering_count < count; x++)
        for (int y = 0; y < count && s->centering_count < count; y++)
            for (int z = 0; z < count && s->centering_count < count; z++) {
                rat v[3] = {rat_int(x), rat_int(y), rat_int(z)}, image[3];
                rat_row_mul(v, (const rat (*)[3])s->to_cinter, image);
                fold3(image);
                int seen = 0;
                for (int i = 0; i < s->centering_count && !seen; i++) seen = eq3(s->centering[i], image);
                if (!seen) memcpy(s->centering[s->centering_count++], image, sizeof image);
            }
    if (s->centering_count != count) return MG_TABLES;
    return mag_bad;
}

static int pml_rotation(int sg, int op, rat (*out)[3])
{
    int64_t lattice = at(&M.lattice, sg - 1);
    int64_t code = at(&M.psettings, (lattice - 1) * 72 + op - 1);
    if (!code) return 1;
    for (int i = 0; i < 9; i++) {
        out[i / 3][i % 3] = rat_int(code % 3 - 1);
        code /= 3;
    }
    return 0;
}

static int cinter_point_operation(const rat (*matrix)[3])
{
    int64_t target[9];
    for (int i = 0; i < 9; i++) {
        rat value = matrix[i / 3][i % 3];
        if (value.d != 1) return 0;
        target[i] = (int64_t)value.n;
    }
    for (int op = 1; op <= 72; op++) {
        int same = 1;
        for (int i = 0; i < 9 && same; i++) same = at(&M.point_op, (op - 1) * 9 + i) == target[i];
        if (same) return op;
    }
    return 0;
}

static int pml_point_operation(int sg, const rat (*matrix)[3])
{
    for (int op = 1; op <= 72; op++) {
        rat candidate[3][3];
        if (pml_rotation(sg, op, candidate)) continue;
        int same = 1;
        for (int i = 0; i < 9 && same; i++) same = rat_eq(candidate[i / 3][i % 3], matrix[i / 3][i % 3]);
        if (same) return op;
    }
    return 0;
}

static int magnetic_point_operation(int ordinary, int antiunitary)
{
    int64_t count = section_count(&M.m2n);
    for (int64_t i = 0; i < count; i++)
        if (at(&M.m2n, i) == ordinary && (at(&M.r, i) != 0) == (antiunitary != 0)) return (int)i + 1;
    return 0;
}

static int ordinary_of(int magnetic_op) { return (int)at(&M.m2n, magnetic_op - 1); }

int mag_tables(void)
{
    mag_bad = 0;
    return load();
}

int mag_take_error(void)
{
    int status = mag_bad;
    mag_bad = 0;
    return status;
}

int64_t mag_group_total(void) { return section_count(&M.sg); }
int mag_ordinary_group(int group) { return (int)at(&M.sg, group - 1); }
int64_t mag_group_code(int group) { return at(&M.point_group_code, at(&M.point_group, group - 1) - 1); }
int64_t mag_op_code(int op) { return at(&M.op_code, op - 1); }
int mag_nonmagnetic_op(int op) { return ordinary_of(op); }
int mag_time_reversed(int op) { return at(&M.r, op - 1) != 0; }

int mag_point_action(int mode, int op, int *sign)
{
    *sign = mode == 1 && mag_time_reversed(op) ? -1 : 1;
    return mode ? mag_nonmagnetic_op(op) : op;
}
int mag_point_product(int left, int right)
{
    if (load() || left < 1 || left > 144 || right < 1 || right > 144) {
        fail(MG_TABLES);
        return 0;
    }
    return M.products[(right - 1) * 144 + left - 1];
}
int64_t mag_point_product_count(void) { return 144 * 144; }
int mag_orderparam_group(int sg) { return (int)at(&M.sg_pointer, sg - 1) + 1; }

int mag_generators(int group, int64_t (*out)[5])
{
    int count = (int)at(&M.gen_count, group - 1);
    if (count < 0 || count > 4) {
        fail(MG_TABLES);
        return 0;
    }
    for (int i = 0; i < count; i++)
        for (int k = 0; k < 5; k++) out[i][k] = at(&M.gen, (int64_t)(group - 1) * 20 + i * 5 + k);
    return count;
}

static void cinter_rotation_of(int ordinary, rat (*out)[3])
{
    for (int i = 0; i < 9; i++) out[i / 3][i % 3] = rat_int(at(&M.point_op, (ordinary - 1) * 9 + i));
}

static void pml_record_to_cinter(const group_setting *s, const int64_t *record, int64_t *out)
{
    rat rotation[3][3], cinter_rotation[3][3];
    if (pml_rotation(s->sg, ordinary_of((int)record[4]), rotation)) {
        fail(MG_POINT_OP);
        return;
    }
    mat_mul((const rat (*)[3])s->to_pml, (const rat (*)[3])rotation, cinter_rotation);
    mat_mul((const rat (*)[3])cinter_rotation, (const rat (*)[3])s->to_cinter, cinter_rotation);
    rat translation[3], transformed[3], rotated_origin[3];
    values_of(record, translation);
    rat_row_mul(translation, (const rat (*)[3])s->to_cinter, transformed);
    rat_row_mul(s->origin_cinter, (const rat (*)[3])cinter_rotation, rotated_origin);
    for (int axis = 0; axis < 3; axis++)
        transformed[axis] = rat_add(rat_sub(transformed[axis], rotated_origin[axis]), s->origin_cinter[axis]);
    int ordinary = cinter_point_operation((const rat (*)[3])cinter_rotation);
    int magnetic = ordinary ? magnetic_point_operation(ordinary, (int)at(&M.r, record[4] - 1)) : 0;
    if (!magnetic) {
        fail(MG_POINT_OP);
        return;
    }
    fold3(transformed);
    record_of(transformed, out);
    out[4] = magnetic;
}

static void cinter_record_to_pml(const group_setting *s, const int64_t *record, int64_t *out)
{
    rat cinter_rotation[3][3], pml[3][3];
    cinter_rotation_of(ordinary_of((int)record[4]), cinter_rotation);
    mat_mul((const rat (*)[3])s->to_cinter, (const rat (*)[3])cinter_rotation, pml);
    mat_mul((const rat (*)[3])pml, (const rat (*)[3])s->to_pml, pml);
    rat translation[3], rotated_origin[3], shifted[3], pml_translation[3];
    values_of(record, translation);
    rat_row_mul(s->origin_cinter, (const rat (*)[3])cinter_rotation, rotated_origin);
    for (int axis = 0; axis < 3; axis++)
        shifted[axis] = rat_sub(rat_add(translation[axis], rotated_origin[axis]), s->origin_cinter[axis]);
    rat_row_mul(shifted, (const rat (*)[3])s->to_pml, pml_translation);
    int ordinary = pml_point_operation(s->sg, (const rat (*)[3])pml);
    int magnetic = ordinary ? magnetic_point_operation(ordinary, (int)at(&M.r, record[4] - 1)) : 0;
    if (!magnetic) {
        fail(MG_POINT_OP);
        return;
    }
    record_of(pml_translation, out);
    out[4] = magnetic;
}

static void compose(const group_setting *s, const int64_t *left, const int64_t *right, int64_t *out)
{
    rat r[3], rotated[3], l[3], t[3];
    values_of(right, r);
    geo_vrot(s->sg, ordinary_of((int)left[4]), r, rotated);
    values_of(left, l);
    for (int axis = 0; axis < 3; axis++) t[axis] = rat_add(l[axis], rotated[axis]);
    fold3(t);
    record_of(t, out);
    out[4] = mag_point_product((int)left[4], (int)right[4]);
}

typedef struct {
    rat f[4][3];
} formula;

typedef struct {
    int ordinal;
    formula bns, cinter;
} wyckoff_row;

typedef struct {
    int row;
    formula f;
    int fixed;
    rat folded_base[3];
    int constraint_count;
    int64_t vector[3][3], base[3][3], base_scale[3];
    rat reduced[3][3], transform[3][3];
    int pivots[3], pivot_count;
} candidate;

typedef struct {
    int ready;
    group_setting setting;
    int count;
    int64_t closure[MG_MAX_OPS][5], binary[MG_MAX_OPS][5], cinter[MG_MAX_OPS][5];
    int rows_ready, row_count;
    wyckoff_row rows[MG_MAX_ROWS];
    int candidate_count;
    candidate *candidates;
} group_cache;

static group_cache **groups;

static int by_point_op(const void *a, const void *b)
{
    const int64_t *x = a, *y = b;
    return (x[4] > y[4]) - (x[4] < y[4]);
}

static group_cache *group_for(int group)
{
    if (load()) {
        fail(MG_TABLES);
        return 0;
    }
    int64_t total = section_count(&M.type);
    if (group < 1 || group > total) {
        fail(MG_RANGE);
        return 0;
    }
    if (!groups) {
        groups = calloc((size_t)total + 1, sizeof *groups);
        if (!groups) {
            fail(MG_MEMORY);
            return 0;
        }
    }
    if (groups[group]) return groups[group];
    group_cache *g = calloc(1, sizeof *g);
    if (!g) {
        fail(MG_MEMORY);
        return 0;
    }
    int status = load_setting(group, &g->setting);
    if (status) {
        free(g);
        fail(status);
        return 0;
    }
    int count = (int)at(&M.gen_count, group - 1);
    for (int i = 0; i < count; i++)
        for (int k = 0; k < 5; k++) g->closure[i][k] = at(&M.gen, (int64_t)(group - 1) * 20 + i * 5 + k);
    int seen[145] = {0};
    for (int i = 0; i < count; i++)
        if (g->closure[i][4] >= 0 && g->closure[i][4] <= 144) seen[g->closure[i][4]] = 1;
    for (int changed = 1; changed && !mag_bad;) {
        changed = 0;
        int outer = count;
        for (int li = 0; li < outer && !mag_bad; li++) {
            int inner = count;
            for (int ri = 0; ri < inner && !mag_bad; ri++) {
                int64_t product[5];
                compose(&g->setting, g->closure[li], g->closure[ri], product);
                if (product[4] < 1 || product[4] > 144) {
                    fail(MG_POINT_OP);
                    break;
                }
                if (seen[product[4]]) continue;
                if (count >= MG_MAX_OPS) {
                    fail(MG_CAPACITY);
                    break;
                }
                seen[product[4]] = 1;
                memcpy(g->closure[count++], product, sizeof product);
                changed = 1;
            }
        }
    }
    qsort(g->closure, (size_t)count, sizeof g->closure[0], by_point_op);
    g->count = count;
    for (int i = 0; i < count && !mag_bad; i++) {
        pml_record_to_cinter(&g->setting, g->closure[i], g->cinter[i]);
        cinter_record_to_pml(&g->setting, g->cinter[i], g->binary[i]);
    }
    if (mag_bad || rat_overflow) {
        if (!mag_bad) mag_bad = MG_RANGE;
        free(g);
        return 0;
    }
    groups[group] = g;
    return g;
}

/* nice_wyckoff_: the canonical spatial formula of a PML formula. */
static int nice_formula(const group_setting *s, const formula *in, formula *out)
{
    rat params[3][3];
    int active[3], active_count = 0;
    for (int i = 0; i < 3; i++) {
        rat_row_mul(in->f[i + 1], (const rat (*)[3])s->to_cinter, params[i]);
        if (!zero3(params[i])) active[active_count++] = i;
    }
    rat m[3][3];
    memcpy(m, params, sizeof m);
    int pivots[3], rank = 0;
    for (int col = 0; col < 3 && rank < 3; col++) {
        int pivot = -1;
        for (int row = rank; row < 3 && pivot < 0; row++)
            if (!rat_zero(m[row][col])) pivot = row;
        if (pivot < 0) continue;
        for (int j = 0; j < 3; j++) {
            rat t = m[rank][j];
            m[rank][j] = m[pivot][j];
            m[pivot][j] = t;
        }
        rat scale = m[rank][col];
        for (int j = 0; j < 3; j++) m[rank][j] = rat_div(m[rank][j], scale);
        for (int row = 0; row < 3; row++) {
            if (row == rank || rat_zero(m[row][col])) continue;
            rat factor = m[row][col];
            for (int j = 0; j < 3; j++) m[row][j] = rat_sub(m[row][j], rat_mul(factor, m[rank][j]));
        }
        pivots[rank++] = col;
    }
    if (rank != active_count) return MG_DEPENDENT;
    rat base[3];
    rat_row_mul(in->f[0], (const rat (*)[3])s->to_cinter, base);
    for (int axis = 0; axis < 3; axis++) base[axis] = rat_add(base[axis], s->origin_cinter[axis]);
    fold3(base);
    for (int r = 0; r < rank; r++) {
        rat coefficient = base[pivots[r]];
        for (int axis = 0; axis < 3; axis++) base[axis] = rat_sub(base[axis], rat_mul(coefficient, m[r][axis]));
    }
    for (int axis = 0; axis < 3; axis++) base[axis] = rat_sub(base[axis], s->origin_cinter[axis]);
    rat_row_mul(base, (const rat (*)[3])s->to_pml, out->f[0]);
    for (int i = 1; i < 4; i++)
        for (int axis = 0; axis < 3; axis++) out->f[i][axis] = rat_int(0);
    for (int r = 0; r < rank; r++) rat_row_mul(m[r], (const rat (*)[3])s->to_pml, out->f[active[r] + 1]);
    return 0;
}

static void formula_to_cinter(const group_setting *s, const formula *in, formula *out)
{
    formula t;
    rat_row_mul(in->f[0], (const rat (*)[3])s->to_cinter, t.f[0]);
    for (int axis = 0; axis < 3; axis++) t.f[0][axis] = rat_add(t.f[0][axis], s->origin_cinter[axis]);
    for (int i = 1; i < 4; i++) rat_row_mul(in->f[i], (const rat (*)[3])s->to_cinter, t.f[i]);
    *out = t;
}

static void reduce_bns_base(const group_setting *s, const rat *base, rat *out)
{
    rat t[3];
    rat_row_mul(base, (const rat (*)[3])s->to_cinter, t);
    for (int axis = 0; axis < 3; axis++) t[axis] = rat_add(t[axis], s->origin_cinter[axis]);
    fold3(t);
    for (int axis = 0; axis < 3; axis++) t[axis] = rat_sub(t[axis], s->origin_cinter[axis]);
    rat_row_mul(t, (const rat (*)[3])s->to_pml, out);
}

/* _formula_solver: Gauss-Jordan on the axis x parameter matrix. */
static void formula_solver(const formula *f, candidate *c)
{
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            c->reduced[row][col] = f->f[col + 1][row];
            c->transform[row][col] = rat_int(row == col);
        }
    c->pivot_count = 0;
    int pivot_row = 0;
    for (int col = 0; col < 3 && pivot_row < 3; col++) {
        int source = -1;
        for (int row = pivot_row; row < 3 && source < 0; row++)
            if (!rat_zero(c->reduced[row][col])) source = row;
        if (source < 0) continue;
        if (source != pivot_row)
            for (int j = 0; j < 3; j++) {
                rat t = c->reduced[pivot_row][j];
                c->reduced[pivot_row][j] = c->reduced[source][j];
                c->reduced[source][j] = t;
                t = c->transform[pivot_row][j];
                c->transform[pivot_row][j] = c->transform[source][j];
                c->transform[source][j] = t;
            }
        rat scale = c->reduced[pivot_row][col];
        for (int j = 0; j < 3; j++) {
            c->reduced[pivot_row][j] = rat_div(c->reduced[pivot_row][j], scale);
            c->transform[pivot_row][j] = rat_div(c->transform[pivot_row][j], scale);
        }
        for (int row = 0; row < 3; row++) {
            if (row == pivot_row) continue;
            rat factor = c->reduced[row][col];
            if (rat_zero(factor)) continue;
            for (int j = 0; j < 3; j++) {
                c->reduced[row][j] = rat_sub(c->reduced[row][j], rat_mul(factor, c->reduced[pivot_row][j]));
                c->transform[row][j] = rat_sub(c->transform[row][j], rat_mul(factor, c->transform[pivot_row][j]));
            }
        }
        c->pivots[c->pivot_count++] = col;
        pivot_row++;
    }
}

static i128 lcm_den(const rat *v, int n)
{
    i128 den = 1;
    for (int i = 0; i < n; i++) den = rat_mul_i(den / rat_gcd(den, v[i].d), v[i].d);
    return den;
}

static int load_rows(group_cache *g)
{
    if (g->rows_ready) return 0;
    const group_setting *s = &g->setting;
    int start, count;
    geo_wyckoff_row_range(s->reference_sg, &start, &count);
    if (count > MG_MAX_ROWS) return MG_CAPACITY;
    rat inverse[3][3], origin[3];
    if (s->type == 4) {
        rat basis[3][3];
        for (int i = 0; i < 9; i++) basis[i / 3][i % 3] = rat_int(s->basis[i]);
        if (rat_inverse3((const rat (*)[3])basis, inverse)) return MG_SINGULAR;
        values_of(s->origin, origin);
    }
    for (int r = 0; r < count; r++) {
        formula raw, nice;
        geo_wyckoff_vectors(start + r, raw.f);
        if (s->type == 4) {
            rat shifted[3];
            for (int axis = 0; axis < 3; axis++) shifted[axis] = rat_sub(raw.f[0][axis], origin[axis]);
            rat_row_mul(shifted, (const rat (*)[3])inverse, raw.f[0]);
            for (int i = 1; i < 4; i++) rat_row_mul(raw.f[i], (const rat (*)[3])inverse, raw.f[i]);
        }
        int status = nice_formula(s, &raw, &nice);
        if (status) return status;
        g->rows[r].ordinal = r + 1;
        g->rows[r].bns = nice;
        formula_to_cinter(s, &nice, &g->rows[r].cinter);
    }
    g->row_count = count;

    int capacity = 0, used = 0;
    for (int r = 0; r < count && !mag_bad; r++) {
        formula branches[MG_MAX_OPS];
        int branch_count = 1;
        branches[0] = g->rows[r].bns;
        for (int op = 1; op < g->count; op++) {
            const int64_t *record = g->binary[op];
            int ordinary = ordinary_of((int)record[4]);
            rat translation[3];
            values_of(record, translation);
            formula f;
            geo_vrot(s->sg, ordinary, g->rows[r].bns.f[0], f.f[0]);
            for (int axis = 0; axis < 3; axis++) f.f[0][axis] = rat_add(f.f[0][axis], translation[axis]);
            for (int i = 1; i < 4; i++) geo_vrot(s->sg, ordinary, g->rows[r].bns.f[i], f.f[i]);
            int same = 0;
            for (int b = 0; b < branch_count && !same; b++) {
                int integral = 1;
                for (int axis = 0; axis < 3 && integral; axis++)
                    integral = rat_sub(f.f[0][axis], branches[b].f[0][axis]).d == 1;
                same = integral;
                for (int i = 1; i < 4 && same; i++) same = eq3(f.f[i], branches[b].f[i]);
            }
            if (same) continue;
            reduce_bns_base(s, f.f[0], f.f[0]);
            branches[branch_count++] = f;
        }
        for (int b = 0; b < branch_count; b++) {
            formula cinter;
            formula_to_cinter(s, &branches[b], &cinter);
            for (int c = 0; c < s->centering_count; c++) {
                if (used >= capacity) {
                    capacity = capacity ? 2 * capacity : 64;
                    candidate *grown = realloc(g->candidates, sizeof *grown * (size_t)capacity);
                    if (!grown) return MG_MEMORY;
                    g->candidates = grown;
                }
                candidate *cand = &g->candidates[used++];
                memset(cand, 0, sizeof *cand);
                cand->row = r;
                cand->f = cinter;
                for (int axis = 0; axis < 3; axis++) cand->f.f[0][axis] = rat_add(cinter.f[0][axis], s->centering[c][axis]);
                if (zero3(cand->f.f[1]) && zero3(cand->f.f[2]) && zero3(cand->f.f[3])) {
                    cand->fixed = 1;
                    memcpy(cand->folded_base, cand->f.f[0], sizeof cand->folded_base);
                    fold3(cand->folded_base);
                    continue;
                }
                formula_solver(&cand->f, cand);
                i128 base_scale = lcm_den(cand->f.f[0], 3);
                if (base_scale > INT64_MAX) rat_overflow = 1;
                for (int row = 0; row < 3; row++) {
                    if (!rat_zero(cand->reduced[row][0]) || !rat_zero(cand->reduced[row][1]) ||
                        !rat_zero(cand->reduced[row][2]))
                        continue;
                    int k = cand->constraint_count++;
                    i128 scale = lcm_den(cand->transform[row], 3);
                    for (int axis = 0; axis < 3; axis++) {
                        cand->vector[k][axis] = (int64_t)(cand->transform[row][axis].n * (scale / cand->transform[row][axis].d));
                        cand->base[k][axis] = (int64_t)(cand->f.f[0][axis].n * (base_scale / cand->f.f[0][axis].d));
                    }
                    cand->base_scale[k] = (int64_t)base_scale;
                }
            }
        }
    }
    g->candidate_count = used;
    g->rows_ready = 1;
    return mag_bad;
}

static i128 floor_rat(rat a)
{
    i128 q = a.n / a.d;
    if ((a.n % a.d) && a.n < 0) q--;
    return q;
}

static i128 ceil_rat(rat a) { return -floor_rat(rat_neg(a)); }

static int formula_solution(const candidate *c, const rat *point, rat *out)
{
    i128 start[3], stop[3];
    for (int axis = 0; axis < 3; axis++) {
        rat d = rat_sub(c->f.f[0][axis], point[axis]);
        start[axis] = floor_rat(d) - 1;
        stop[axis] = ceil_rat(d) + 2;
    }
    i128 shift[3] = {start[0], start[1], start[2]};
    if (c->constraint_count) {
        rat delta[3];
        for (int axis = 0; axis < 3; axis++) delta[axis] = rat_sub(point[axis], c->f.f[0][axis]);
        int constraint_rows[3], n = 0;
        for (int row = 0; row < 3; row++)
            if (rat_zero(c->reduced[row][0]) && rat_zero(c->reduced[row][1]) && rat_zero(c->reduced[row][2]))
                constraint_rows[n++] = row;
        rat required[3];
        for (int k = 0; k < n; k++) {
            rat acc = rat_int(0);
            for (int axis = 0; axis < 3; axis++) acc = rat_add(acc, rat_mul(c->transform[constraint_rows[k]][axis], delta[axis]));
            required[k] = rat_neg(acc);
        }
        int found = 0;
        for (i128 x = start[0]; x < stop[0] && !found; x++)
            for (i128 y = start[1]; y < stop[1] && !found; y++)
                for (i128 z = start[2]; z < stop[2] && !found; z++) {
                    i128 s[3] = {x, y, z};
                    int ok = 1;
                    for (int k = 0; k < n && ok; k++) {
                        rat acc = rat_int(0);
                        for (int axis = 0; axis < 3; axis++)
                            acc = rat_add(acc, rat_mul(c->transform[constraint_rows[k]][axis], rat_int(s[axis])));
                        ok = rat_eq(acc, required[k]);
                    }
                    if (ok) {
                        memcpy(shift, s, sizeof s);
                        found = 1;
                    }
                }
        if (!found) return 0;
    }
    rat target[3];
    for (int axis = 0; axis < 3; axis++) target[axis] = rat_sub(rat_add(point[axis], rat_int(shift[axis])), c->f.f[0][axis]);
    for (int i = 0; i < 3; i++) out[i] = rat_int(0);
    for (int row = 0; row < c->pivot_count; row++) {
        rat acc = rat_int(0);
        for (int axis = 0; axis < 3; axis++) acc = rat_add(acc, rat_mul(c->transform[row][axis], target[axis]));
        out[c->pivots[row]] = acc;
    }
    return 1;
}

/* id_wyckoff2_magnetic_ in cinter: the row ordinal and folded
 * representative of a folded cinter point, or 0. */
static int identify(group_cache *g, const rat *point, rat *representative)
{
    if (load_rows(g)) return 0;
    rat folded[3];
    memcpy(folded, point, sizeof folded);
    fold3(folded);
    i128 den = lcm_den(folded, 3), numerators[3];
    for (int axis = 0; axis < 3; axis++) numerators[axis] = rat_mul_i(folded[axis].n, den / folded[axis].d);
    for (int i = 0; i < g->candidate_count; i++) {
        const candidate *c = &g->candidates[i];
        rat solution[3] = {rat_int(0), rat_int(0), rat_int(0)};
        if (c->fixed) {
            if (!eq3(c->folded_base, folded)) continue;
        } else {
            int admits = 1;
            for (int k = 0; k < c->constraint_count && admits; k++) {
                i128 common = rat_mul_i(den / rat_gcd(den, c->base_scale[k]), c->base_scale[k]);
                i128 total = 0;
                for (int axis = 0; axis < 3; axis++)
                    total = rat_add_i(total, rat_mul_i(c->vector[k][axis],
                                                       rat_add_i(rat_mul_i(numerators[axis], common / den),
                                                                 -rat_mul_i(c->base[k][axis], common / c->base_scale[k]))));
                admits = total % common == 0;
            }
            if (!admits || !formula_solution(c, folded, solution)) continue;
        }
        for (int axis = 0; axis < 3; axis++) {
            rat value = c->f.f[0][axis];
            for (int p = 0; p < 3; p++) value = rat_add(value, rat_mul(solution[p], c->f.f[p + 1][axis]));
            representative[axis] = value;
        }
        fold3(representative);
        return g->rows[c->row].ordinal;
    }
    return 0;
}

/* The full centered orbit of a cinter point in Source operation order. */
static int orbit_points(group_cache *g, const rat *point, rat (*out)[3], int capacity)
{
    int count = 0;
    for (int op = 0; op < g->count; op++) {
        rat rotation[3][3], translation[3], base[3];
        cinter_rotation_of(ordinary_of((int)g->cinter[op][4]), rotation);
        values_of(g->cinter[op], translation);
        rat_row_mul(point, (const rat (*)[3])rotation, base);
        for (int axis = 0; axis < 3; axis++) base[axis] = rat_add(base[axis], translation[axis]);
        for (int c = 0; c < g->setting.centering_count; c++) {
            rat image[3];
            for (int axis = 0; axis < 3; axis++) image[axis] = rat_mod1(rat_add(base[axis], g->setting.centering[c][axis]));
            int seen = 0;
            for (int i = 0; i < count && !seen; i++) seen = eq3(out[i], image);
            if (seen) continue;
            if (count >= capacity) {
                fail(MG_CAPACITY);
                return count;
            }
            memcpy(out[count++], image, sizeof image);
        }
    }
    return count;
}

/* records_out: the group's records in setting (closure, binary or cinter);
 * returns the count or a negative error. */
int magnetic_records(int group, int setting, int capacity, int64_t *records_out)
{
    mag_bad = 0;
    rat_overflow = 0;
    group_cache *g = group_for(group);
    if (!g) return -(mag_bad ? mag_bad : MG_TABLES);
    if (g->count > capacity) return -MG_CAPACITY;
    const int64_t (*source)[5] = setting == SET_CLOSURE ? g->closure : setting == SET_BINARY ? g->binary : g->cinter;
    memcpy(records_out, source, sizeof(int64_t) * 5 * (size_t)g->count);
    return g->count;
}

static int find(int *parent, int i)
{
    while (parent[i] != i) i = parent[i] = parent[parent[i]];
    return i;
}

static void unite(int *parent, int a, int b)
{
    a = find(parent, a);
    b = find(parent, b);
    if (a < b)
        parent[b] = a;
    else if (b < a)
        parent[a] = b;
}

/* The magnetic regrouping of one child atom population.  Points are int128
 * (num, den) words per coordinate: parent_points (selected-basis parent
 * PML) and displayed (presentation cinter).  atom_order: (source raw index,
 * centering ordinal); sites: CSR atom indices of the ordinary child sites.
 * components_out per component: representative atom, multiplicity, Source
 * row ordinal, presentation row ordinal, first site, site count; sites_out
 * the component sites.  labels_ok: 1 when the presentation rows form a
 * bijective relabeling. */
int magnetic_site_groups(int group, int child_sg, int atom_count, const int64_t *parent_points,
                         const int64_t *atom_order, int site_count, const int64_t *site_offsets,
                         const int64_t *site_atoms, const int64_t *displayed, const int64_t *basis,
                         const int64_t *origin, int capacity, int64_t *components_out, int64_t *sites_out,
                         int *component_count, int *labels_ok)
{
    mag_bad = 0;
    rat_overflow = 0;
    group_cache *g = group_for(group);
    if (!g) return mag_bad ? mag_bad : MG_TABLES;
    int rows_status = load_rows(g);
    if (rows_status) return rows_status;
    if (atom_count <= 0) return MG_PARTITION;
    if (!origin[3]) return MG_RANGE;
    rat result_basis[3][3], inverse[3][3], result_origin[3];
    for (int i = 0; i < 9; i++) result_basis[i / 3][i % 3] = rat_int(basis[i]);
    if (rat_inverse3((const rat (*)[3])result_basis, inverse)) return MG_SINGULAR;
    for (int axis = 0; axis < 3; axis++) result_origin[axis] = rat_make(origin[axis], origin[3]);
    rat child_to_cinter[3][3], child_origin[3];
    if (geo_setting_to_cinter(child_sg, 1, 0, child_to_cinter, child_origin)) return MG_TABLES;

    int status = MG_OK;
    rat (*pml)[3] = malloc(sizeof(rat[3]) * (size_t)atom_count);
    rat (*cinter)[3] = malloc(sizeof(rat[3]) * (size_t)atom_count);
    rat (*shown)[3] = malloc(sizeof(rat[3]) * (size_t)atom_count);
    int *parent = malloc(sizeof(int) * (size_t)atom_count);
    int *site_of = malloc(sizeof(int) * (size_t)atom_count);
    rat (*orbit)[3] = malloc(sizeof(rat[3]) * (size_t)(atom_count + 1));
    if (!pml || !cinter || !shown || !parent || !site_of || !orbit) {
        status = MG_MEMORY;
        goto done;
    }
    for (int a = 0; a < atom_count; a++) {
        rat p[3];
        read_rat3(parent_points + 12 * a, p);
        read_rat3(displayed + 12 * a, shown[a]);
        for (int axis = 0; axis < 3; axis++) p[axis] = rat_sub(p[axis], result_origin[axis]);
        rat_row_mul(p, (const rat (*)[3])inverse, pml[a]);
        rat_row_mul(pml[a], (const rat (*)[3])child_to_cinter, cinter[a]);
        for (int axis = 0; axis < 3; axis++) cinter[a][axis] = rat_add(cinter[a][axis], child_origin[axis]);
        fold3(cinter[a]);
        for (int b = 0; b < a && !status; b++)
            if (eq3(cinter[a], cinter[b])) status = MG_SHARED_POINT;
        parent[a] = a;
        site_of[a] = -1;
    }
    if (status) goto done;
    for (int s = 0; s < site_count && !status; s++) {
        if (site_offsets[s + 1] <= site_offsets[s]) status = MG_PARTITION;
        for (int64_t i = site_offsets[s]; i < site_offsets[s + 1] && !status; i++) {
            int64_t a = site_atoms[i];
            if (a < 0 || a >= atom_count || site_of[a] >= 0) {
                status = MG_PARTITION;
                break;
            }
            site_of[a] = s;
            unite(parent, (int)site_atoms[site_offsets[s]], (int)a);
        }
    }
    for (int a = 0; a < atom_count && !status; a++)
        if (site_of[a] < 0) status = MG_PARTITION;
    if (status) goto done;

    for (int op = 0; op < g->count && !status; op++) {
        rat rotation[3][3], linear[3][3], translation[3], shift[3];
        if (pml_rotation(child_sg, ordinary_of((int)g->binary[op][4]), rotation)) {
            status = MG_POINT_OP;
            break;
        }
        mat_mul((const rat (*)[3])rotation, (const rat (*)[3])child_to_cinter, linear);
        rat t[3];
        for (int axis = 0; axis < 3; axis++) t[axis] = rat_make(g->binary[op][axis], g->binary[op][3]);
        rat_row_mul(t, (const rat (*)[3])child_to_cinter, translation);
        for (int axis = 0; axis < 3; axis++) shift[axis] = rat_add(translation[axis], child_origin[axis]);
        unsigned char *hit = calloc((size_t)atom_count, 1);
        if (!hit) {
            status = MG_MEMORY;
            break;
        }
        for (int a = 0; a < atom_count && !status; a++) {
            rat image[3];
            rat_row_mul(pml[a], (const rat (*)[3])linear, image);
            for (int axis = 0; axis < 3; axis++) image[axis] = rat_add(image[axis], shift[axis]);
            fold3(image);
            int target = -1;
            for (int b = 0; b < atom_count && target < 0; b++)
                if (eq3(cinter[b], image)) target = b;
            if (target < 0)
                status = MG_LEAVES;
            else if (hit[target])
                status = MG_PERMUTATION;
            else {
                hit[target] = 1;
                unite(parent, a, target);
            }
        }
        free(hit);
    }
    if (status) goto done;

    int components = 0, written_sites = 0, presentation_ok = 1;
    int64_t *canonical_label = malloc(sizeof(int64_t) * (size_t)(capacity + 1));
    int64_t *presentation_label = malloc(sizeof(int64_t) * (size_t)(capacity + 1));
    int label_pairs = 0;
    if (!canonical_label || !presentation_label) status = MG_MEMORY;
    for (int root = 0; root < atom_count && !status; root++) {
        if (find(parent, root) != root) continue;
        if (components >= capacity) {
            status = MG_CAPACITY;
            break;
        }
        int64_t *component = components_out + 6 * components++;
        int size = 0, seed = -1;
        component[4] = written_sites;
        for (int a = 0; a < atom_count; a++) {
            if (find(parent, a) != root) continue;
            size++;
            if (seed < 0 || atom_order[2 * a] < atom_order[2 * seed] ||
                (atom_order[2 * a] == atom_order[2 * seed] && atom_order[2 * a + 1] < atom_order[2 * seed + 1]))
                seed = a;
            int s = site_of[a], seen = 0;
            for (int i = (int)component[4]; i < written_sites && !seen; i++) seen = sites_out[i] == s;
            if (!seen) {
                if (written_sites >= capacity) {
                    status = MG_CAPACITY;
                    break;
                }
                sites_out[written_sites++] = s;
            }
        }
        if (status) break;
        component[5] = written_sites - component[4];
        rat standard[3];
        int ordinal = identify(g, cinter[seed], standard);
        if (!ordinal) {
            status = mag_bad ? mag_bad : MG_NO_FORMULA;
            break;
        }
        int orbit_count = orbit_points(g, cinter[seed], orbit, atom_count + 1);
        if (orbit_count != size) {
            status = MG_ORBIT;
            break;
        }
        for (int i = 0; i < orbit_count && !status; i++) {
            int member = 0;
            for (int a = 0; a < atom_count && !member; a++) member = find(parent, a) == root && eq3(cinter[a], orbit[i]);
            if (!member) status = MG_ORBIT;
        }
        if (status) break;
        int match = -1, matches = 0;
        for (int a = 0; a < atom_count; a++)
            if (find(parent, a) == root && eq3(cinter[a], standard)) {
                if (match < 0) match = a;
                matches++;
            }
        if (matches != 1) {
            status = MG_REPRESENTATIVE;
            break;
        }
        component[0] = match;
        component[1] = size;
        component[2] = ordinal;
        component[3] = 0;
        if (presentation_ok) {
            rat shown_standard[3];
            int shown_ordinal = identify(g, shown[match], shown_standard);
            mag_bad = 0;
            int shown_count = shown_ordinal ? orbit_points(g, shown[match], orbit, atom_count + 1) : 0;
            mag_bad = 0;
            int consistent = shown_ordinal && shown_count == size;
            for (int i = 0; i < label_pairs && consistent; i++) {
                if (canonical_label[i] == ordinal && presentation_label[i] != shown_ordinal) consistent = 0;
                if (presentation_label[i] == shown_ordinal && canonical_label[i] != ordinal) consistent = 0;
            }
            if (!consistent) {
                presentation_ok = 0;
            } else {
                canonical_label[label_pairs] = ordinal;
                presentation_label[label_pairs++] = shown_ordinal;
                component[3] = shown_ordinal;
            }
        }
    }
    if (!status) {
        for (int i = 0; i < components; i++) {
            int64_t *component = components_out + 6 * i;
            int64_t size = 0;
            for (int64_t k = component[4]; k < component[4] + component[5]; k++)
                size += site_offsets[sites_out[k] + 1] - site_offsets[sites_out[k]];
            if (size != component[1]) status = MG_SPLIT_SITE;
        }
    }
    free(canonical_label);
    free(presentation_label);
    *component_count = components;
    *labels_ok = presentation_ok;
    if (!status && (mag_bad || rat_overflow)) status = mag_bad ? mag_bad : MG_RANGE;
done:
    free(pml);
    free(cinter);
    free(shown);
    free(parent);
    free(site_of);
    free(orbit);
    return status;
}

/* _reduce_fraction_vector: positive denominator, common factor of the
 * nonzero entries divided out. */
static int reduce_vector(int64_t *v)
{
    if (!v[3]) return 1;
    if (v[3] < 0)
        for (int i = 0; i < 4; i++) v[i] = -v[i];
    int64_t g = 0;
    for (int i = 0; i < 4; i++)
        if (v[i]) g = (int64_t)gcd_i128(g, v[i]);
    if (g > 1)
        for (int i = 0; i < 4; i++) v[i] /= g;
    return 0;
}

static int64_t floor_mod(int64_t a, int64_t m)
{
    int64_t r = a % m;
    if (r && ((r < 0) != (m < 0))) r += m;
    return r;
}

/* SourceData._reduce_operation_record */
static int reduce_operation(int64_t *record)
{
    if (!record[3]) return 0;
    int64_t v[4] = {floor_mod(record[0], record[3]), floor_mod(record[1], record[3]), floor_mod(record[2], record[3]),
                    record[3]};
    if (reduce_vector(v)) return 1;
    for (int i = 0; i < 3; i++) record[i] = floor_mod(v[i], v[3]);
    record[3] = v[3];
    return 0;
}

/* SourceData._compose_magnetic_operation_records */
static int compose_iso(int sg, const int64_t *left, const int64_t *right, int64_t *out)
{
    int64_t lattice = at(&M.lattice, sg - 1);
    int64_t backward[9], forward[9], point[9];
    section_view lattice_ml, lattice_den;
    if (section_find("space.lattice_ml", &lattice_ml) || section_find("space.lattice_ml_denom", &lattice_den)) return 1;
    for (int i = 0; i < 9; i++) {
        backward[i] = at(&lattice_ml, (lattice - 1) * 36 + i);
        forward[i] = at(&lattice_ml, (lattice - 1) * 36 + 9 + i);
        point[i] = at(&M.point_op, (int64_t)(ordinary_of((int)left[4]) - 1) * 9 + i);
    }
    int64_t setting_den = at(&lattice_den, (lattice - 1) * 2);
    if (!right[3] || !left[3]) return 1;
    int64_t first[3], second[3], rotated[4];
    for (int c = 0; c < 3; c++) {
        first[c] = 0;
        for (int r = 0; r < 3; r++) first[c] += right[r] * forward[c + 3 * r];
    }
    for (int c = 0; c < 3; c++) {
        second[c] = 0;
        for (int r = 0; r < 3; r++) second[c] += first[r] * point[c + 3 * r];
    }
    for (int c = 0; c < 3; c++) {
        rotated[c] = 0;
        for (int r = 0; r < 3; r++) rotated[c] += second[r] * backward[c + 3 * r];
    }
    rotated[3] = right[3] * setting_den;
    if (reduce_vector(rotated)) return 1;
    int64_t den = (int64_t)(left[3] / gcd_i128(left[3], rotated[3]) * rotated[3]);
    if (den < 0) den = -den;
    int64_t v[4];
    for (int a = 0; a < 3; a++) v[a] = left[a] * (den / left[3]) + rotated[a] * (den / rotated[3]);
    v[3] = den;
    if (reduce_vector(v)) return 1;
    for (int i = 0; i < 4; i++) out[i] = v[i];
    out[4] = mag_point_product((int)left[4], (int)right[4]);
    return reduce_operation(out);
}

static int64_t (*iso_closures)[MG_MAX_OPS][5];
static int *iso_counts;

/* SourceData.generate_magnetic_space_group_records: the closure with each
 * pass scanning the records present when it began. */
int magnetic_records_iso(int group, int capacity, int64_t *records_out)
{
    mag_bad = 0;
    if (load()) return -MG_TABLES;
    int64_t total = section_count(&M.type);
    if (group < 1 || group > total) return -MG_RANGE;
    if (!iso_closures) {
        iso_closures = calloc((size_t)total + 1, sizeof *iso_closures);
        iso_counts = calloc((size_t)total + 1, sizeof *iso_counts);
        if (!iso_closures || !iso_counts) return -MG_MEMORY;
    }
    if (!iso_counts[group]) {
        int64_t (*records)[5] = iso_closures[group];
        int sg = (int)at(&M.sg, group - 1);
        int count = (int)at(&M.gen_count, group - 1);
        int seen[145] = {0};
        for (int i = 0; i < count; i++) {
            for (int k = 0; k < 5; k++) records[i][k] = at(&M.gen, (int64_t)(group - 1) * 20 + i * 5 + k);
            if (reduce_operation(records[i])) return -MG_RANGE;
            if (records[i][4] >= 1 && records[i][4] <= 144) seen[records[i][4]] = 1;
        }
        for (int changed = 1; changed;) {
            changed = 0;
            int scan = count;
            for (int li = 0; li < scan; li++)
                for (int ri = 0; ri < scan; ri++) {
                    int64_t product[5];
                    if (compose_iso(sg, records[li], records[ri], product)) return -MG_RANGE;
                    if (product[4] < 1 || product[4] > 144) return -MG_POINT_OP;
                    if (seen[product[4]]) continue;
                    if (count >= MG_MAX_OPS) return -MG_CAPACITY;
                    seen[product[4]] = 1;
                    memcpy(records[count++], product, sizeof product);
                    changed = 1;
                }
        }
        qsort(records, (size_t)count, sizeof records[0], by_point_op);
        if (mag_bad) return -mag_bad;
        iso_counts[group] = count;
    }
    if (iso_counts[group] > capacity) return -MG_CAPACITY;
    memcpy(records_out, iso_closures[group], sizeof(int64_t) * 5 * (size_t)iso_counts[group]);
    return iso_counts[group];
}
