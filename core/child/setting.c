/* The displayed child setting: the exact coordinate map from the Source
 * default child setting to the presentation, the Source inter settings whose
 * space-group action it reproduces, and the Wyckoff identity of a child
 * orbit in those settings. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"

enum {
    CH_TABLES = 1,
    CH_RANGE = 2,
    CH_SINGULAR = 3,
    CH_LATTICE = 4,
    CH_DUPLICATE = 5,
    CH_COLLAPSE = 6,
    CH_NO_SETTING = 7,
    CH_DEPENDENT = 8,
    CH_DUPLICATE_POINTS = 9,
    CH_NOT_UNIQUE = 10,
    CH_DISAGREE = 11,
    CH_CAPACITY = 12,
};

#define CH_MAX_SETTINGS 32
#define CH_MAX_ACTION (GEO_MAX_OPS * 8)

typedef struct {
    rat r[3][3], t[3];
} affine;

typedef struct {
    int count;
    affine ops[CH_MAX_ACTION];
} action;

static int ch_bad;

static void mat_mul(const rat (*a)[3], const rat (*b)[3], rat (*out)[3])
{
    rat t[3][3];
    for (int row = 0; row < 3; row++) rat_row_mul(a[row], b, t[row]);
    memcpy(out, t, sizeof t);
}

static void mat4_mul(const i128 *left, const i128 *right, i128 *out)
{
    i128 t[16];
    for (int row = 0; row < 4; row++)
        for (int col = 0; col < 4; col++) {
            i128 total = 0;
            for (int k = 0; k < 4; k++) total = rat_add_i(total, rat_mul_i(left[4 * row + k], right[4 * k + col]));
            t[4 * row + col] = total;
        }
    memcpy(out, t, sizeof t);
}

/* trans4_(first, middle, last) = last . middle . first */
static void trans4(const int64_t *first, i128 *middle, const int64_t *last)
{
    i128 f[16], l[16];
    for (int i = 0; i < 16; i++) {
        f[i] = first[i];
        l[i] = last[i];
    }
    mat4_mul(middle, f, middle);
    mat4_mul(l, middle, middle);
}

/* subgroup_change_setting_ to cinter as the exact affine map and its inverse. */
static int change_setting(int parent_sg, int subgroup, const int64_t *basis, const int64_t *origin, int parent_setting,
                          int subgroup_setting, rat (*matrix)[3], rat *shift, rat (*inverse)[3], rat *inverse_shift)
{
    i128 den = origin[3] ? origin[3] : 1, m[16];
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) m[4 * row + col] = rat_mul_i(basis[3 * row + col], den);
        m[4 * row + 3] = 0;
    }
    for (int col = 0; col < 3; col++) m[12 + col] = origin[col];
    m[15] = den;
    int64_t first[16], last[16];
    geo_cinter_base(parent_sg, 1, first);
    geo_cinter_base(subgroup, 0, last);
    trans4(first, m, last);
    geo_inter_block(parent_setting ? parent_setting : geo_default_setting(parent_sg), 0, first);
    geo_inter_block(subgroup_setting ? subgroup_setting : geo_default_setting(subgroup), 1, last);
    trans4(first, m, last);
    if (!m[15]) return CH_TABLES;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) matrix[row][col] = rat_make(m[4 * row + col], m[15]);
    for (int col = 0; col < 3; col++) shift[col] = rat_make(m[12 + col], m[15]);
    if (rat_inverse3((const rat (*)[3])matrix, inverse)) return CH_SINGULAR;
    rat_row_mul(shift, (const rat (*)[3])inverse, inverse_shift);
    for (int col = 0; col < 3; col++) inverse_shift[col] = rat_neg(inverse_shift[col]);
    return 0;
}

static rat rat_det3(const rat (*m)[3])
{
    rat det = rat_int(0);
    for (int col = 0; col < 3; col++)
        det = rat_add(det, rat_mul(m[0][col], rat_sub(rat_mul(m[1][(col + 1) % 3], m[2][(col + 2) % 3]),
                                                       rat_mul(m[1][(col + 2) % 3], m[2][(col + 1) % 3]))));
    return det;
}

static int same3(const rat *a, const rat *b) { return rat_eq(a[0], b[0]) && rat_eq(a[1], b[1]) && rat_eq(a[2], b[2]); }

static int same_affine(const affine *a, const affine *b)
{
    return same3(a->r[0], b->r[0]) && same3(a->r[1], b->r[1]) && same3(a->r[2], b->r[2]) && same3(a->t, b->t);
}

static int action_has(const action *a, const affine *op)
{
    for (int i = 0; i < a->count; i++)
        if (same_affine(&a->ops[i], op)) return 1;
    return 0;
}

static int action_add(action *a, const affine *op)
{
    if (action_has(a, op)) return 0;
    if (a->count == CH_MAX_ACTION) {
        if (!ch_bad) ch_bad = CH_CAPACITY;
        return 0;
    }
    a->ops[a->count++] = *op;
    return 1;
}

static int same_action(const action *a, const action *b)
{
    if (a->count != b->count) return 0;
    for (int i = 0; i < a->count; i++)
        if (!action_has(b, &a->ops[i])) return 0;
    return 1;
}

static int lattice_offsets(const rat (*matrix)[3], rat (*out)[3], int capacity)
{
    i128 den = 1;
    for (int i = 0; i < 9; i++) den = rat_mul_i(den / rat_gcd(den, matrix[i / 3][i % 3].d), matrix[i / 3][i % 3].d);
    int count = 0;
    for (i128 x = 0; x < den; x++)
        for (i128 y = 0; y < den; y++)
            for (i128 z = 0; z < den; z++) {
                rat v[3] = {rat_int(x), rat_int(y), rat_int(z)}, image[3];
                rat_row_mul(v, matrix, image);
                for (int axis = 0; axis < 3; axis++) image[axis] = rat_mod1(image[axis]);
                int seen = 0;
                for (int i = 0; i < count && !seen; i++) seen = same3(out[i], image);
                if (seen) continue;
                if (count == capacity) return -CH_LATTICE;
                memcpy(out[count++], image, sizeof image);
            }
    rat det = rat_det3(matrix);
    if (rat_zero(det)) return -CH_SINGULAR;
    rat expected = rat_div(rat_int(1), det);
    if (expected.n < 0) expected.n = -expected.n;
    if (expected.d != 1 || expected.n != count) return -CH_LATTICE;
    return count;
}

/* The space-group action of sg in the affine setting (matrix, origin). */
static int setting_action(int sg, const rat (*matrix)[3], const rat *origin, action *out)
{
    rat inverse[3][3], offsets[8][3];
    if (rat_inverse3(matrix, inverse)) return CH_SINGULAR;
    int offset_count = lattice_offsets(matrix, offsets, 8);
    if (offset_count < 0) return -offset_count;
    int64_t records[GEO_MAX_OPS][5];
    int record_count = geo_sg_records(sg, records);
    out->count = 0;
    for (int g = 0; g < record_count; g++) {
        rat rotation[3][3], t[3], translated[3], rotated_origin[3];
        for (int col = 0; col < 3; col++) {
            rat unit[3] = {rat_int(col == 0), rat_int(col == 1), rat_int(col == 2)};
            geo_vrot(sg, (int)records[g][4], unit, rotation[col]);
        }
        affine op;
        mat_mul((const rat (*)[3])rotation, matrix, op.r);
        mat_mul((const rat (*)[3])inverse, (const rat (*)[3])op.r, op.r);
        for (int axis = 0; axis < 3; axis++) t[axis] = rat_make(records[g][axis], records[g][3]);
        rat_row_mul(t, matrix, translated);
        rat_row_mul(origin, (const rat (*)[3])op.r, rotated_origin);
        rat base[3];
        for (int axis = 0; axis < 3; axis++)
            base[axis] = rat_mod1(rat_sub(rat_add(translated[axis], origin[axis]), rotated_origin[axis]));
        for (int o = 0; o < offset_count; o++) {
            for (int axis = 0; axis < 3; axis++) op.t[axis] = rat_mod1(rat_add(base[axis], offsets[o][axis]));
            action_add(out, &op);
        }
    }
    if (out->count != record_count * offset_count) return CH_DUPLICATE;
    return ch_bad;
}

static int transform_action(const action *in, const rat (*matrix)[3], const rat *origin, action *out)
{
    rat inverse[3][3];
    if (rat_inverse3(matrix, inverse)) return CH_SINGULAR;
    out->count = 0;
    for (int i = 0; i < in->count; i++) {
        affine op;
        rat translated[3], rotated_origin[3];
        mat_mul((const rat (*)[3])in->ops[i].r, matrix, op.r);
        mat_mul((const rat (*)[3])inverse, (const rat (*)[3])op.r, op.r);
        rat_row_mul(in->ops[i].t, matrix, translated);
        rat_row_mul(origin, (const rat (*)[3])op.r, rotated_origin);
        for (int axis = 0; axis < 3; axis++)
            op.t[axis] = rat_mod1(rat_sub(rat_add(translated[axis], origin[axis]), rotated_origin[axis]));
        action_add(out, &op);
    }
    if (out->count != in->count) return CH_COLLAPSE;
    return ch_bad;
}

static action *cached_action(int sg, int setting_id)
{
    static action *cache[1024];
    static int keys[1024];
    if (setting_id < 1 || setting_id >= 1024) {
        if (!ch_bad) ch_bad = CH_TABLES;
        return 0;
    }
    if (keys[setting_id] == sg) return cache[setting_id];
    action *a = cache[setting_id] ? cache[setting_id] : malloc(sizeof *a);
    if (!a) {
        if (!ch_bad) ch_bad = CH_RANGE;
        return 0;
    }
    cache[setting_id] = a;
    keys[setting_id] = 0;
    rat matrix[3][3], origin[3], cml[3][3], cml_origin[3];
    geo_setting_to_cinter(sg, 1, setting_id, matrix, origin);
    geo_setting_to_cinter(sg, 0, setting_id, cml, cml_origin);
    int status = setting_action(sg, (const rat (*)[3])matrix, cml_origin, a);
    if (status) {
        if (!ch_bad) ch_bad = status;
        return 0;
    }
    keys[setting_id] = sg;
    return a;
}

static int finish(int status)
{
    if (!status) status = ch_bad;
    if (!status && geo_take_error()) status = CH_TABLES;
    if (!status && rat_overflow) status = CH_RANGE;
    return status;
}

/* origin_words: the presentation origin (3 rationals as int128 words).
 * Writes the child-to-display matrix (9) and origin (3) as rational words
 * and the matching inter setting ids; returns their count or a negative
 * status. */
int child_presentation_setting(int parent_sg, int child_sg, int parent_setting, const int64_t *basis,
                               const int64_t *origin, const int64_t *presentation_pml, const int64_t *origin_words,
                               int capacity, int *setting_ids, int64_t *matrix_words, int64_t *origin_out_words)
{
    ch_bad = 0;
    rat_overflow = 0;
    if (geo_load()) return -CH_TABLES;
    geo_take_error();
    int default_id = geo_default_setting(child_sg);
    rat c2p[3][3], c2p_origin[3], unused[3][3], unused_origin[3];
    int status = change_setting(parent_sg, child_sg, basis, origin, parent_setting, default_id, c2p, c2p_origin,
                                unused, unused_origin);
    if (status) return -finish(status);
    rat presentation[3][3], to_cinter[3][3], to_origin[3], display[3][3], display_inverse[3][3];
    for (int i = 0; i < 9; i++) presentation[i / 3][i % 3] = rat_int(presentation_pml[i]);
    geo_setting_to_cinter(parent_sg, 1, parent_setting, to_cinter, to_origin);
    mat_mul((const rat (*)[3])presentation, (const rat (*)[3])to_cinter, display);
    if (rat_inverse3((const rat (*)[3])display, display_inverse)) return -finish(CH_SINGULAR);
    rat c2d[3][3], c2d_origin[3], display_origin[3], shifted[3];
    mat_mul((const rat (*)[3])c2p, (const rat (*)[3])display_inverse, c2d);
    read_rat3(origin_words, display_origin);
    for (int axis = 0; axis < 3; axis++) shifted[axis] = rat_sub(c2p_origin[axis], display_origin[axis]);
    rat_row_mul(shifted, (const rat (*)[3])display_inverse, c2d_origin);
    int ids[CH_MAX_SETTINGS], candidates[CH_MAX_SETTINGS];
    int id_count = geo_inter_setting_ids(child_sg, CH_MAX_SETTINGS, ids), candidate_count = 0;
    const int64_t identity[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}, zero[4] = {0, 0, 0, 1};
    for (int i = 0; i < id_count; i++) {
        rat m[3][3], o[3], inverse[3][3], inverse_origin[3];
        status = change_setting(child_sg, child_sg, identity, zero, default_id, ids[i], m, o, inverse, inverse_origin);
        if (status) return -finish(status);
        int match = same3(inverse[0], c2d[0]) && same3(inverse[1], c2d[1]) && same3(inverse[2], c2d[2]);
        for (int axis = 0; axis < 3 && match; axis++)
            match = rat_zero(rat_mod1(rat_sub(inverse_origin[axis], c2d_origin[axis])));
        if (match) candidates[candidate_count++] = ids[i];
    }
    if (!candidate_count) {
        memcpy(candidates, ids, sizeof(int) * (size_t)id_count);
        candidate_count = id_count;
    }
    action *base = cached_action(child_sg, default_id);
    static action target;
    if (!base) return -finish(0);
    status = transform_action(base, (const rat (*)[3])c2d, c2d_origin, &target);
    if (status) return -finish(status);
    int selected = 0;
    if (same_action(base, &target))
        setting_ids[selected++] = default_id;
    else
        for (int i = 0; i < candidate_count; i++) {
            action *a = cached_action(child_sg, candidates[i]);
            if (!a) return -finish(0);
            if (!same_action(a, &target)) continue;
            if (selected == capacity) return -finish(CH_CAPACITY);
            setting_ids[selected++] = candidates[i];
        }
    if (!selected) return -finish(CH_NO_SETTING);
    for (int row = 0; row < 3; row++) write_rat3(matrix_words + 12 * row, c2d[row]);
    write_rat3(origin_out_words, c2d_origin);
    status = finish(0);
    return status ? -status : selected;
}

static i128 floor_div(rat a)
{
    i128 q = a.n / a.d;
    return (a.n % a.d && a.n < 0) ? q - 1 : q;
}

/* One Source Wyckoff representative solved exactly modulo the lattice:
 * 1 with solution, 0 when none, negative on dependent parameters. */
static int solve_parameters(const rat (*vectors)[3], const rat *target, rat *solution, int *count_out)
{
    rat parameters[3][3];
    int count = 0;
    for (int v = 1; v < 4; v++)
        if (!rat_zero(vectors[v][0]) || !rat_zero(vectors[v][1]) || !rat_zero(vectors[v][2]))
            memcpy(parameters[count++], vectors[v], sizeof parameters[0]);
    *count_out = count;
    if (!count) {
        for (int axis = 0; axis < 3; axis++)
            if (!rat_zero(rat_mod1(rat_sub(target[axis], vectors[0][axis])))) return 0;
        return 1;
    }
    i128 low[3], high[3];
    for (int axis = 0; axis < 3; axis++) {
        rat lower = rat_int(0), upper = rat_int(0);
        for (int c = 0; c < count; c++) {
            if (parameters[c][axis].n < 0) lower = rat_add(lower, parameters[c][axis]);
            if (parameters[c][axis].n > 0) upper = rat_add(upper, parameters[c][axis]);
        }
        rat delta = rat_sub(target[axis], vectors[0][axis]);
        low[axis] = -floor_div(rat_neg(rat_sub(lower, delta)));
        high[axis] = floor_div(rat_sub(upper, delta));
    }
    i128 shift[3];
    for (shift[0] = low[0]; shift[0] <= high[0]; shift[0]++)
        for (shift[1] = low[1]; shift[1] <= high[1]; shift[1]++)
            for (shift[2] = low[2]; shift[2] <= high[2]; shift[2]++) {
                rat aug[3][4];
                for (int row = 0; row < 3; row++) {
                    for (int c = 0; c < count; c++) aug[row][c] = parameters[c][row];
                    aug[row][count] = rat_sub(rat_add(target[row], rat_int(shift[row])), vectors[0][row]);
                }
                int pivot_row = 0, pivots[3], pivot_count = 0;
                for (int c = 0; c < count; c++) {
                    int pivot = -1;
                    for (int row = pivot_row; row < 3 && pivot < 0; row++)
                        if (!rat_zero(aug[row][c])) pivot = row;
                    if (pivot < 0) continue;
                    rat swap[4];
                    memcpy(swap, aug[pivot_row], sizeof swap);
                    memcpy(aug[pivot_row], aug[pivot], sizeof swap);
                    memcpy(aug[pivot], swap, sizeof swap);
                    rat divisor = aug[pivot_row][c];
                    for (int k = 0; k <= count; k++) aug[pivot_row][k] = rat_div(aug[pivot_row][k], divisor);
                    for (int row = 0; row < 3; row++) {
                        if (row == pivot_row || rat_zero(aug[row][c])) continue;
                        rat factor = aug[row][c];
                        for (int k = 0; k <= count; k++)
                            aug[row][k] = rat_sub(aug[row][k], rat_mul(factor, aug[pivot_row][k]));
                    }
                    pivots[pivot_count++] = c;
                    pivot_row++;
                }
                if (pivot_count != count) return -CH_DEPENDENT;
                int inconsistent = 0;
                for (int row = 0; row < 3 && !inconsistent; row++) {
                    int empty = 1;
                    for (int c = 0; c < count; c++) empty &= rat_zero(aug[row][c]);
                    inconsistent = empty && !rat_zero(aug[row][count]);
                }
                if (inconsistent) continue;
                for (int c = 0; c < count; c++) solution[c] = rat_int(0);
                for (int row = 0; row < pivot_count; row++) solution[pivots[row]] = aug[row][count];
                int in_cell = 1;
                for (int c = 0; c < count && in_cell; c++)
                    in_cell = solution[c].n >= 0 && rat_cmp(solution[c], rat_int(1)) < 0;
                if (!in_cell) continue;
                int exact = 1;
                for (int axis = 0; axis < 3 && exact; axis++) {
                    rat total = rat_int(0);
                    for (int c = 0; c < count; c++) total = rat_add(total, rat_mul(parameters[c][axis], solution[c]));
                    exact = rat_eq(total, rat_sub(rat_add(target[axis], rat_int(shift[axis])), vectors[0][axis]));
                }
                if (exact) return 1;
            }
    return 0;
}

static int orbit_matches(const rat (*vectors)[3], const action *a, const rat (*points)[3], int point_count)
{
    rat active[3][3];
    int active_count = 0;
    for (int v = 1; v < 4; v++)
        if (!rat_zero(vectors[v][0]) || !rat_zero(vectors[v][1]) || !rat_zero(vectors[v][2]))
            memcpy(active[active_count++], vectors[v], sizeof active[0]);
    for (int p = 0; p < point_count; p++) {
        rat solution[3];
        int count;
        int solved = solve_parameters(vectors, points[p], solution, &count);
        if (solved < 0) return solved;
        if (!solved) continue;
        rat representative[3];
        for (int axis = 0; axis < 3; axis++) {
            rat value = vectors[0][axis];
            for (int c = 0; c < count; c++) value = rat_add(value, rat_mul(active[c][axis], solution[c]));
            representative[axis] = rat_mod1(value);
        }
        int orbit_count = 0, equal = 1;
        static rat orbit[CH_MAX_ACTION][3];
        for (int i = 0; i < a->count && equal; i++) {
            rat image[3];
            rat_row_mul(representative, (const rat (*)[3])a->ops[i].r, image);
            for (int axis = 0; axis < 3; axis++) image[axis] = rat_mod1(rat_add(image[axis], a->ops[i].t[axis]));
            int seen = 0;
            for (int j = 0; j < orbit_count && !seen; j++) seen = same3(orbit[j], image);
            if (seen) continue;
            int member = 0;
            for (int j = 0; j < point_count && !member; j++) member = same3(points[j], image);
            if (!member) equal = 0;
            memcpy(orbit[orbit_count++], image, sizeof image);
        }
        if (equal && orbit_count == point_count) return 1;
    }
    return 0;
}

/* The child Wyckoff row (0-based offset) of one Formula15 orbit in every
 * displayed setting: xyz_words holds the atoms' structure coordinates. */
int child_public_wyckoff(int child_sg, int setting_count, const int *setting_ids, const int64_t *matrix_words,
                         const int64_t *origin_words, int point_count, const int64_t *xyz_words, int site_pg)
{
    ch_bad = 0;
    rat_overflow = 0;
    if (geo_load()) return -CH_TABLES;
    geo_take_error();
    rat matrix[3][3], origin[3];
    for (int row = 0; row < 3; row++) read_rat3(matrix_words + 12 * row, matrix[row]);
    read_rat3(origin_words, origin);
    rat (*points)[3] = malloc(sizeof *points * (size_t)(point_count + 1));
    if (!points) return -CH_RANGE;
    for (int p = 0; p < point_count; p++) {
        rat xyz[3];
        read_rat3(xyz_words + 12 * p, xyz);
        rat_row_mul(xyz, (const rat (*)[3])matrix, points[p]);
        for (int axis = 0; axis < 3; axis++) points[p][axis] = rat_mod1(rat_add(points[p][axis], origin[axis]));
        for (int q = 0; q < p; q++)
            if (same3(points[q], points[p])) {
                free(points);
                return -finish(CH_DUPLICATE_POINTS);
            }
    }
    int start, row_count, agreed = -1, status = 0;
    geo_wyckoff_row_range(child_sg, &start, &row_count);
    for (int s = 0; s < setting_count && !status; s++) {
        action *a = cached_action(child_sg, setting_ids[s]);
        if (!a) {
            status = CH_TABLES;
            break;
        }
        rat to_cinter[3][3], unused_origin[3], cml[3][3], cml_origin[3];
        geo_setting_to_cinter(child_sg, 1, setting_ids[s], to_cinter, unused_origin);
        geo_setting_to_cinter(child_sg, 0, setting_ids[s], cml, cml_origin);
        int found = -1, matches = 0;
        for (int offset = 0; offset < row_count; offset++) {
            int row_id = start + offset, row_pg = geo_wyckoff_site_pg(row_id);
            if (row_pg != site_pg || geo_multiplicity(child_sg, row_pg) != point_count) continue;
            rat pml[4][3], vectors[4][3];
            geo_inter_wyckoff_vectors(setting_ids[s], offset, pml);
            for (int v = 0; v < 4; v++) rat_row_mul(pml[v], (const rat (*)[3])to_cinter, vectors[v]);
            for (int axis = 0; axis < 3; axis++) vectors[0][axis] = rat_add(vectors[0][axis], cml_origin[axis]);
            int hit = orbit_matches((const rat (*)[3])vectors, a, (const rat (*)[3])points, point_count);
            if (hit < 0) {
                status = -hit;
                break;
            }
            if (hit) {
                matches++;
                found = offset;
            }
        }
        if (status) break;
        if (matches != 1) status = CH_NOT_UNIQUE;
        else if (agreed < 0) agreed = found;
    }
    free(points);
    status = finish(status);
    return status ? -status : agreed;
}
