/* get_new_wyckoff_ Formula15 provenance for child_sites.c: transform the
 * parent Wyckoff branches into the child cell, group them into child orbits,
 * match each orbit to its first Source child row and solve Formula15.  All
 * arithmetic is on exact i128 rationals; overflow is a range error. */

#include <stdint.h>

#include "intmath.h"
#include <stdlib.h>
#include <string.h>


enum {
    WS_OK = 0,
    WS_RANGE = 1,
    WS_IMAGE = 2,
    WS_COLLISION = 3,
    WS_NO_ROW = 4,
    WS_FORMULA = 5,
    WS_CAPACITY = 6,
    WS_SITE_PG = 7,
};

typedef struct {
    i128 n, d;
} rat;

static int ws_overflow;

static i128 ws_gcd(i128 a, i128 b) { return gcd_i128(a, b); }

static i128 ws_mul_i(i128 a, i128 b)
{
    i128 out;
    if (__builtin_mul_overflow(a, b, &out)) {
        ws_overflow = 1;
        return 0;
    }
    return out;
}

static i128 ws_add_i(i128 a, i128 b)
{
    i128 out;
    if (__builtin_add_overflow(a, b, &out)) {
        ws_overflow = 1;
        return 0;
    }
    return out;
}

static rat ws_rat(i128 n, i128 d)
{
    if (!d) {
        ws_overflow = 1;
        return (rat){0, 1};
    }
    if (d < 0) {
        n = -n;
        d = -d;
    }
    i128 g = ws_gcd(n, d);
    if (g > 1) {
        n /= g;
        d /= g;
    }
    return (rat){n, d};
}

static rat ws_add(rat a, rat b)
{
    i128 g = ws_gcd(a.d, b.d);
    return ws_rat(ws_add_i(ws_mul_i(a.n, b.d / g), ws_mul_i(b.n, a.d / g)), ws_mul_i(a.d / g, b.d));
}

static rat ws_neg(rat a) { return (rat){-a.n, a.d}; }
static rat ws_sub(rat a, rat b) { return ws_add(a, ws_neg(b)); }
static rat ws_mul(rat a, rat b) { return ws_rat(ws_mul_i(a.n, b.n), ws_mul_i(a.d, b.d)); }
static rat ws_div(rat a, rat b) { return ws_rat(ws_mul_i(a.n, b.d), ws_mul_i(a.d, b.n)); }
static rat ws_int(i128 v) { return (rat){v, 1}; }
static int ws_zero(rat a) { return a.n == 0; }
static int ws_eq(rat a, rat b) { return a.n == b.n && a.d == b.d; }

static rat ws_mod1(rat a)
{
    i128 r = a.n % a.d;
    if (r < 0) r += a.d;
    return (rat){r, a.d};
}

/* A 4-int record (x, y, z, den) in _fraction_record normal form. */
typedef struct {
    int64_t v[4];
} record;

static record ws_record(const rat *values)
{
    i128 den = 1;
    for (int i = 0; i < 3; i++) den = ws_mul_i(den / ws_gcd(den, values[i].d), values[i].d);
    i128 n[3], g = den;
    for (int i = 0; i < 3; i++) {
        n[i] = ws_mul_i(values[i].n, den / values[i].d);
        g = ws_gcd(g, n[i]);
    }
    if (!g) g = 1;
    record out;
    for (int i = 0; i < 3; i++) {
        i128 value = n[i] / g;
        if (value > INT32_MAX || value < INT32_MIN) ws_overflow = 1;
        out.v[i] = (int64_t)value;
    }
    i128 d = den / g;
    if (d > INT32_MAX) ws_overflow = 1;
    out.v[3] = (int64_t)d;
    return out;
}

static void ws_values(const int64_t *rec, rat *out)
{
    for (int i = 0; i < 3; i++) out[i] = ws_rat(rec[i], rec[3]);
}

static int ws_record_zero(const int64_t *rec) { return !rec[0] && !rec[1] && !rec[2]; }

typedef struct {
    record r[4];
} branch;

static int ws_branch_eq(const branch *a, const branch *b) { return !memcmp(a, b, sizeof *a); }

/* row-vector product v . m for a 3x3 rational matrix. */
static void ws_row_mul(const rat *v, const rat *m, rat *out)
{
    for (int col = 0; col < 3; col++) {
        rat total = ws_int(0);
        for (int row = 0; row < 3; row++) total = ws_add(total, ws_mul(v[row], m[row * 3 + col]));
        out[col] = total;
    }
}

static int ws_rank(const rat *columns, int count)
{
    rat m[3][3];
    for (int row = 0; row < 3; row++)
        for (int c = 0; c < count; c++) m[row][c] = columns[c * 3 + row];
    int rank = 0;
    for (int c = 0; c < count && rank < 3; c++) {
        int pivot = -1;
        for (int row = rank; row < 3; row++)
            if (!ws_zero(m[row][c])) {
                pivot = row;
                break;
            }
        if (pivot < 0) continue;
        for (int k = 0; k < count; k++) {
            rat t = m[rank][k];
            m[rank][k] = m[pivot][k];
            m[pivot][k] = t;
        }
        rat scale = m[rank][c];
        for (int k = 0; k < count; k++) m[rank][k] = ws_div(m[rank][k], scale);
        for (int row = 0; row < 3; row++) {
            if (row == rank || ws_zero(m[row][c])) continue;
            rat factor = m[row][c];
            for (int k = 0; k < count; k++) m[row][k] = ws_sub(m[row][k], ws_mul(factor, m[rank][k]));
        }
        rank++;
    }
    return rank;
}

static void ws_cross(const rat *a, const rat *b, rat *out)
{
    out[0] = ws_sub(ws_mul(a[1], b[2]), ws_mul(a[2], b[1]));
    out[1] = ws_sub(ws_mul(a[2], b[0]), ws_mul(a[0], b[2]));
    out[2] = ws_sub(ws_mul(a[0], b[1]), ws_mul(a[1], b[0]));
}

typedef struct {
    int rank;
    rat basis[3], normal[3];
} span;

static void ws_span(const rat *params, int count, span *s)
{
    s->rank = count ? ws_rank(params, count) : 0;
    for (int i = 0; i < 3; i++) s->basis[i] = s->normal[i] = ws_int(0);
    for (int c = 0; c < count; c++)
        if (!ws_zero(params[3 * c]) || !ws_zero(params[3 * c + 1]) || !ws_zero(params[3 * c + 2])) {
            memcpy(s->basis, params + 3 * c, sizeof s->basis);
            break;
        }
    for (int c = 0; c < count; c++) {
        rat cand[3];
        ws_cross(s->basis, params + 3 * c, cand);
        if (!ws_zero(cand[0]) || !ws_zero(cand[1]) || !ws_zero(cand[2])) {
            memcpy(s->normal, cand, sizeof s->normal);
            break;
        }
    }
}

static int ws_in_span(const span *s, const rat *v)
{
    if (s->rank == 0) return ws_zero(v[0]) && ws_zero(v[1]) && ws_zero(v[2]);
    if (s->rank == 1) {
        rat c[3];
        ws_cross(s->basis, v, c);
        return ws_zero(c[0]) && ws_zero(c[1]) && ws_zero(c[2]);
    }
    if (s->rank == 2) {
        rat total = ws_int(0);
        for (int i = 0; i < 3; i++) total = ws_add(total, ws_mul(s->normal[i], v[i]));
        return ws_zero(total);
    }
    return 1;
}

static int ws_affine_subset(const branch *group, const branch *child)
{
    rat child_params[9], base[3], child_base[3];
    int count = 0;
    for (int i = 1; i < 4; i++)
        if (!ws_record_zero(child->r[i].v)) ws_values(child->r[i].v, child_params + 3 * count++);
    span s;
    ws_span(child_params, count, &s);
    for (int i = 1; i < 4; i++) {
        if (ws_record_zero(group->r[i].v)) continue;
        rat p[3];
        ws_values(group->r[i].v, p);
        if (!ws_in_span(&s, p)) return 0;
    }
    ws_values(group->r[0].v, base);
    ws_values(child->r[0].v, child_base);
    rat delta[3];
    for (int i = 0; i < 3; i++) delta[i] = ws_sub(base[i], child_base[i]);
    for (int a = -2; a <= 2; a++)
        for (int b = -2; b <= 2; b++)
            for (int c = -2; c <= 2; c++) {
                rat shifted[3] = {ws_sub(delta[0], ws_int(a)), ws_sub(delta[1], ws_int(b)), ws_sub(delta[2], ws_int(c))};
                if (ws_in_span(&s, shifted)) return 1;
            }
    return 0;
}

static void ws_factor_column(i128 rows[][3], int nrows, int col)
{
    i128 g = 0;
    for (int r = 0; r < nrows; r++) g = ws_gcd(g, rows[r][col]);
    if (g <= 1) return;
    for (int r = 0; r < nrows; r++) rows[r][col] /= g;
}

static void ws_rowop2(i128 rows[][3], int nrows)
{
    if (!nrows) return;
    int skipped = 0;
    for (int col = 0; col < 3;) {
        int pivot_row = skipped + col;
        if (pivot_row >= nrows) break;
        while (rows[pivot_row][col] == 0) {
            int swap = -1;
            for (int c = col + 1; c < 3; c++)
                if (rows[pivot_row][c] != 0) {
                    swap = c;
                    break;
                }
            if (swap >= 0) {
                for (int r = 0; r < nrows; r++) {
                    i128 t = rows[r][col];
                    rows[r][col] = rows[r][swap];
                    rows[r][swap] = t;
                }
                break;
            }
            skipped++;
            pivot_row = skipped + col;
            if (pivot_row >= nrows) break;
        }
        if (pivot_row >= nrows) break;
        i128 pivot = rows[pivot_row][col];
        for (int target = 0; target < 3; target++) {
            if (target == col || rows[pivot_row][target] == 0) continue;
            i128 value = rows[pivot_row][target];
            for (int r = 0; r < nrows; r++)
                if (r != pivot_row)
                    rows[r][target] = ws_add_i(ws_mul_i(pivot, rows[r][target]), -ws_mul_i(rows[r][col], value));
            rows[pivot_row][target] = 0;
            ws_factor_column(rows, nrows, target);
        }
        col++;
    }
    for (int c = 0; c < 3; c++) ws_factor_column(rows, nrows, c);
}

static int ws_rowop_accepts(const branch *child, const branch *group, const int *parent_active,
                            int parent_active_count, const int itry[3])
{
    int child_active[3], child_count = 0;
    for (int i = 1; i < 4; i++)
        if (!ws_record_zero(child->r[i].v)) child_active[child_count++] = i;
    i128 product = 1;
    for (int i = 0; i < 4; i++) product = ws_mul_i(product, child->r[i].v[3]);
    for (int i = 0; i < 4; i++) product = ws_mul_i(product, group->r[i].v[3]);
    i128 child_base[4];
    for (int i = 0; i < 4; i++) child_base[i] = child->r[0].v[i];
    for (int axis = 0; axis < 3; axis++) child_base[axis] -= child_base[3] * itry[axis];
    i128 rows[8][3];
    int nrows = 0;
    for (int k = 0; k < child_count; k++) {
        const int64_t *v = child->r[child_active[k]].v;
        for (int axis = 0; axis < 3; axis++) rows[nrows][axis] = ws_mul_i(product, v[axis]) / v[3];
        nrows++;
    }
    for (int axis = 0; axis < 3; axis++)
        rows[nrows][axis] = ws_mul_i(product, group->r[0].v[axis]) / group->r[0].v[3]
                          - ws_mul_i(product, child_base[axis]) / child_base[3];
    nrows++;
    for (int k = 0; k < parent_active_count; k++) {
        const int64_t *v = group->r[parent_active[k]].v;
        for (int axis = 0; axis < 3; axis++) rows[nrows][axis] = ws_mul_i(product, v[axis]) / v[3];
        nrows++;
    }
    ws_rowop2(rows, nrows);
    for (int r = 0; r < nrows; r++)
        for (int c = child_count; c < 3; c++)
            if (rows[r][c] != 0) return 0;
    return 1;
}

/* _solve_child_parameter_column */
static int ws_solve_column(const rat *child, const rat *target, int allow_shift, rat *out)
{
    static const int primary[3] = {0, 1, -1};
    int offsets[125][3], count = 0;
    if (allow_shift) {
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                for (int c = 0; c < 3; c++) {
                    offsets[count][0] = primary[a];
                    offsets[count][1] = primary[b];
                    offsets[count][2] = primary[c];
                    count++;
                }
        for (int a = -2; a <= 2; a++)
            for (int b = -2; b <= 2; b++)
                for (int c = -2; c <= 2; c++) {
                    if (a >= -1 && a <= 1 && b >= -1 && b <= 1 && c >= -1 && c <= 1) continue;
                    offsets[count][0] = a;
                    offsets[count][1] = b;
                    offsets[count][2] = c;
                    count++;
                }
    } else {
        offsets[0][0] = offsets[0][1] = offsets[0][2] = 0;
        count = 1;
    }
    for (int o = 0; o < count; o++) {
        rat m[3][4];
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) m[row][col] = child[(col + 1) * 3 + row];
            m[row][3] = ws_sub(target[row], ws_int(offsets[o][row]));
        }
        rat shifted[3] = {m[0][3], m[1][3], m[2][3]};
        int pivots[3], pivot_count = 0, pivot_row = 0;
        for (int col = 0; col < 3; col++) {
            int pivot = -1;
            for (int row = pivot_row; row < 3; row++)
                if (!ws_zero(m[row][col])) {
                    pivot = row;
                    break;
                }
            if (pivot < 0) continue;
            for (int k = 0; k < 4; k++) {
                rat t = m[pivot_row][k];
                m[pivot_row][k] = m[pivot][k];
                m[pivot][k] = t;
            }
            rat scale = m[pivot_row][col];
            for (int k = 0; k < 4; k++) m[pivot_row][k] = ws_div(m[pivot_row][k], scale);
            for (int row = 0; row < 3; row++) {
                if (row == pivot_row || ws_zero(m[row][col])) continue;
                rat factor = m[row][col];
                for (int k = 0; k < 4; k++) m[row][k] = ws_sub(m[row][k], ws_mul(factor, m[pivot_row][k]));
            }
            pivots[pivot_count++] = col;
            pivot_row++;
            if (pivot_row == 3) break;
        }
        int inconsistent = 0;
        for (int row = 0; row < 3; row++)
            if (ws_zero(m[row][0]) && ws_zero(m[row][1]) && ws_zero(m[row][2]) && !ws_zero(m[row][3])) inconsistent = 1;
        if (inconsistent) continue;
        rat v[3] = {ws_int(0), ws_int(0), ws_int(0)};
        for (int k = 0; k < pivot_count; k++) v[pivots[k]] = m[k][3];
        int ok = 1;
        for (int row = 0; row < 3 && ok; row++) {
            rat total = ws_int(0);
            for (int col = 0; col < 3; col++) total = ws_add(total, ws_mul(child[(col + 1) * 3 + row], v[col]));
            ok = ws_eq(total, shifted[row]);
        }
        if (ok) {
            memcpy(out, v, sizeof v);
            return 1;
        }
    }
    return 0;
}

static void ws_formula_block(rat constant, const rat *coefficients, int force, int64_t *out)
{
    rat values[4] = {ws_mod1(constant), coefficients[0], coefficients[1], coefficients[2]};
    if (!force && ws_zero(values[0]) && ws_zero(values[1]) && ws_zero(values[2]) && ws_zero(values[3])) {
        for (int i = 0; i < 5; i++) out[i] = 0;
        return;
    }
    i128 den = 1;
    for (int i = 0; i < 4; i++) den = ws_mul_i(den / ws_gcd(den, values[i].d), values[i].d);
    for (int i = 0; i < 4; i++) {
        i128 value = ws_mul_i(values[i].n, den / values[i].d);
        if (value > INT32_MAX || value < INT32_MIN) ws_overflow = 1;
        out[i] = (int64_t)value;
    }
    if (den > INT32_MAX) ws_overflow = 1;
    out[4] = (int64_t)den;
}

#define WS_MAX_FLAT 65536

/* parent_branches: branch_count x 16 (four 4-int records); basis 9; origin 4;
 * kernel: kernel_count x 4; child ops: op_count x 5 with rotation op_count x 9
 * (row-major integer matrix over rotation_dens[i]); child rows: row_count x
 * (site_pg, offset0) and row_count x 16 normalized records; point_op_codes,
 * point_group_codes: Source tables.  Output per group g (in Source order):
 * out_rows[g] = (child row offset + 1, representative flat index, orbit
 * start in out_orbits, orbit size, matched flat index), out_formula[g] 15 ints. */
int wyckoff_split_provenance(int branch_count, const int64_t *parent_branches, const int64_t *basis,
                             const int64_t *origin, int kernel_count, const int64_t *kernel,
                             int op_count, const int64_t *ops, const int64_t *rotations,
                             const int64_t *rotation_dens, int row_count, const int64_t *row_info,
                             const int64_t *row_records, int point_op_code_count, const int64_t *point_op_codes,
                             int point_group_code_count, const int64_t *point_group_codes, int capacity,
                             int64_t *out_rows, int64_t *out_formula, int64_t *out_orbits, int *group_count)
{
    ws_overflow = 0;
    int flat_count = kernel_count * branch_count;
    if (flat_count <= 0 || flat_count > WS_MAX_FLAT) return WS_CAPACITY;
    branch *flat = malloc(sizeof *flat * (size_t)flat_count);
    int *visited = calloc((size_t)flat_count, sizeof *visited);
    int *orbit_mark = calloc((size_t)flat_count, sizeof *orbit_mark);
    int *affine = malloc(sizeof *affine * (size_t)flat_count);
    int status = WS_OK;
    if (!flat || !visited || !orbit_mark || !affine) {
        status = WS_CAPACITY;
        goto done;
    }
    /* transform_wyckoff_branches_to_child */
    rat b[9], inv[9];
    for (int i = 0; i < 9; i++) b[i] = ws_int(basis[i]);
    {
        rat det = ws_sub(ws_add(ws_mul(b[0], ws_sub(ws_mul(b[4], b[8]), ws_mul(b[5], b[7]))),
                                ws_mul(b[2], ws_sub(ws_mul(b[3], b[7]), ws_mul(b[4], b[6])))),
                         ws_mul(b[1], ws_sub(ws_mul(b[3], b[8]), ws_mul(b[5], b[6]))));
        if (ws_zero(det)) {
            status = WS_RANGE;
            goto done;
        }
        rat adj[9] = {
            ws_sub(ws_mul(b[4], b[8]), ws_mul(b[5], b[7])), ws_sub(ws_mul(b[2], b[7]), ws_mul(b[1], b[8])),
            ws_sub(ws_mul(b[1], b[5]), ws_mul(b[2], b[4])), ws_sub(ws_mul(b[5], b[6]), ws_mul(b[3], b[8])),
            ws_sub(ws_mul(b[0], b[8]), ws_mul(b[2], b[6])), ws_sub(ws_mul(b[2], b[3]), ws_mul(b[0], b[5])),
            ws_sub(ws_mul(b[3], b[7]), ws_mul(b[4], b[6])), ws_sub(ws_mul(b[1], b[6]), ws_mul(b[0], b[7])),
            ws_sub(ws_mul(b[0], b[4]), ws_mul(b[1], b[3])),
        };
        for (int i = 0; i < 9; i++) inv[i] = ws_div(adj[i], det);
    }
    rat origin_values[3], origin_child[3];
    ws_values(origin, origin_values);
    ws_row_mul(origin_values, inv, origin_child);
    for (int k = 0; k < kernel_count; k++) {
        rat kv[3];
        ws_values(kernel + 4 * k, kv);
        for (int br = 0; br < branch_count; br++) {
            const int64_t *src = parent_branches + 16 * br;
            branch *out = flat + k * branch_count + br;
            rat base[3], child_base[3];
            ws_values(src, base);
            for (int i = 0; i < 3; i++) base[i] = ws_add(base[i], kv[i]);
            ws_row_mul(base, inv, child_base);
            for (int i = 0; i < 3; i++) child_base[i] = ws_mod1(ws_sub(child_base[i], origin_child[i]));
            out->r[0] = ws_record(child_base);
            {
                rat again[3];
                ws_values(out->r[0].v, again);
                for (int i = 0; i < 3; i++) again[i] = ws_mod1(again[i]);
                out->r[0] = ws_record(again);
            }
            for (int p = 1; p < 4; p++) {
                rat param[3], child_param[3];
                ws_values(src + 4 * p, param);
                ws_row_mul(param, inv, child_param);
                out->r[p] = ws_record(child_param);
            }
        }
    }
    if (ws_overflow) {
        status = WS_RANGE;
        goto done;
    }
    /* group_child_orbits; index_by_branch keeps the last flat index of equal branches */
    int groups = 0, orbit_fill = 0;
    for (int index = 0; index < flat_count; index++) {
        if (visited[index]) continue;
        if (groups >= capacity) {
            status = WS_CAPACITY;
            goto done;
        }
        int group_id = groups + 1;
        visited[index] = group_id;
        const branch *br = flat + index;
        int active[3], active_count = 0;
        for (int p = 1; p < 4; p++)
            if (!ws_record_zero(br->r[p].v)) active[active_count++] = p;
        memset(orbit_mark, 0, sizeof *orbit_mark * (size_t)flat_count);
        orbit_mark[index] = 1;
        int64_t code = 0;
        for (int o = 0; o < op_count; o++) {
            const int64_t *op = ops + 5 * o, *rot = rotations + 9 * o;
            branch image;
            for (int p = 0; p < 4; p++) {
                const int64_t *v = br->r[p].v;
                rat rotated[3];
                for (int a = 0; a < 3; a++) {
                    i128 total = 0;
                    for (int c = 0; c < 3; c++) total = ws_add_i(total, ws_mul_i(rot[a * 3 + c], v[c]));
                    rotated[a] = ws_rat(total, ws_mul_i(v[3], rotation_dens[o]));
                }
                if (p == 0) {
                    record r = ws_record(rotated);
                    rat values[3], translation[3];
                    ws_values(r.v, values);
                    ws_values(op, translation);
                    for (int a = 0; a < 3; a++) values[a] = ws_add(translation[a], values[a]);
                    record added = ws_record(values);
                    ws_values(added.v, values);
                    for (int a = 0; a < 3; a++) values[a] = ws_mod1(values[a]);
                    image.r[0] = ws_record(values);
                } else {
                    image.r[p] = ws_record(rotated);
                }
            }
            if (ws_overflow) {
                status = WS_RANGE;
                goto done;
            }
            int target = -1;
            for (int t = flat_count - 1; t >= 0; t--)
                if (ws_branch_eq(flat + t, &image)) {
                    target = t;
                    break;
                }
            if (target < 0) {
                status = WS_IMAGE;
                goto done;
            }
            orbit_mark[target] = 1;
            if (visited[target] != 0 && visited[target] != group_id) {
                status = WS_COLLISION;
                goto done;
            }
            visited[target] = group_id;
            if (target == index) {
                int op_index = (int)op[4];
                if (op_index < 1 || op_index > point_op_code_count) {
                    status = WS_RANGE;
                    goto done;
                }
                code += point_op_codes[op_index - 1];
            }
        }
        int site_pg = 0;
        for (int i = 0; i < point_group_code_count; i++)
            if (point_group_codes[i] == code) {
                site_pg = i + 1;
                break;
            }
        int orbit_start = orbit_fill, orbit_size = 0;
        for (int t = 0; t < flat_count; t++)
            if (orbit_mark[t]) {
                if (orbit_fill >= flat_count) {
                    status = WS_CAPACITY;
                    goto done;
                }
                out_orbits[orbit_fill++] = t;
                orbit_size++;
            }
        /* child_row_matches_for_group: first match only */
        if (!site_pg && row_count > 0) {
            status = WS_SITE_PG;
            goto done;
        }
        int matched_row = -1, matched_flat = -1;
        int first_affine_row = -1, first_affine_flat = -1;
        static const int tries[3] = {0, 1, -1};
        for (int r = 0; r < row_count && matched_row < 0; r++) {
            if (!site_pg || row_info[2 * r] != site_pg) continue;
            branch child;
            for (int p = 0; p < 4; p++)
                for (int i = 0; i < 4; i++) child.r[p].v[i] = row_records[16 * r + 4 * p + i];
            int affine_count = 0;
            for (int k = 0; k < orbit_size; k++) {
                int t = (int)out_orbits[orbit_start + k];
                if (ws_affine_subset(flat + t, &child)) affine[affine_count++] = t;
            }
            if (affine_count && first_affine_row < 0) {
                first_affine_row = r;
                first_affine_flat = affine[0];
            }
            for (int a = 0; a < 3 && matched_row < 0; a++)
                for (int bb = 0; bb < 3 && matched_row < 0; bb++)
                    for (int c = 0; c < 3 && matched_row < 0; c++) {
                        int itry[3] = {tries[a], tries[bb], tries[c]};
                        for (int k = 0; k < affine_count; k++)
                            if (ws_rowop_accepts(&child, flat + affine[k], active, active_count, itry)) {
                                matched_row = r;
                                matched_flat = affine[k];
                                break;
                            }
                    }
            if (ws_overflow) {
                status = WS_RANGE;
                goto done;
            }
        }
        if (matched_row < 0) {
            matched_row = first_affine_row;
            matched_flat = first_affine_flat;
        }
        if (matched_row < 0) {
            status = WS_NO_ROW;
            goto done;
        }
        /* formula15_for_child_row */
        rat child_vectors[12], branch_vectors[12];
        for (int p = 0; p < 4; p++) {
            ws_values(row_records + 16 * matched_row + 4 * p, child_vectors + 3 * p);
            ws_values(flat[matched_flat].r[p].v, branch_vectors + 3 * p);
        }
        rat solved[4][3];
        for (int column = 0; column < 4; column++) {
            rat target[3];
            for (int a = 0; a < 3; a++)
                target[a] = column == 0 ? ws_sub(branch_vectors[a], child_vectors[a]) : branch_vectors[3 * column + a];
            if (!ws_solve_column(child_vectors, target, column == 0, solved[column])) {
                status = WS_FORMULA;
                goto done;
            }
        }
        for (int axis = 0; axis < 3; axis++) {
            int active_axis = !ws_zero(child_vectors[3 * (axis + 1)]) || !ws_zero(child_vectors[3 * (axis + 1) + 1])
                           || !ws_zero(child_vectors[3 * (axis + 1) + 2]);
            rat coefficients[3] = {solved[1][axis], solved[2][axis], solved[3][axis]};
            ws_formula_block(solved[0][axis], coefficients, active_axis, out_formula + 15 * groups + 5 * axis);
        }
        if (ws_overflow) {
            status = WS_RANGE;
            goto done;
        }
        out_rows[5 * groups] = row_info[2 * matched_row + 1] + 1;
        out_rows[5 * groups + 1] = index;
        out_rows[5 * groups + 2] = orbit_start;
        out_rows[5 * groups + 3] = orbit_size;
        out_rows[5 * groups + 4] = matched_flat;
        groups++;
    }
    *group_count = groups;
done:
    free(flat);
    free(visited);
    free(orbit_mark);
    free(affine);
    return status;
}
