/* reciprocal.catalog.k_coordinate_map_for_slot (the affine map behind the
 * displayed K-family expression, display conventions included) and
 * request_context._pml_vector_to_case_k_params. */

#include <stdint.h>
#include <string.h>

#include "geometry.h"
#include "tables.h"

enum {
    KC_TABLES = 1,
    KC_MAP = 2,
    KC_RANGE = 3,
};

static const char NAMES[3] = {'a', 'b', 'g'};

typedef struct {
    int dim, count;
    char names[3];
    rat origin[3], columns[3][3];
} kmap;

#include "k_display.h"

static struct {
    section_view k, k_dim, k_label, label_bc, irr_sg, irr_k;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (geo_load() || section_find("little.little_k", &T.k) || section_find("little.little_k_dim", &T.k_dim) ||
        section_find("little.little_k_label", &T.k_label) || section_find("space.space_label_bc", &T.label_bc) ||
        section_find("little.little_irr_space_group", &T.irr_sg) || section_find("little.little_irr_k", &T.irr_k))
        return KC_TABLES;
    T.ready = 1;
    return 0;
}

static void trimmed(const section_view *view, int64_t index, const char **text, int64_t *length)
{
    if (section_string(view, index, text, length)) {
        *text = "";
        *length = 0;
        return;
    }
    while (*length && (**text == ' ' || **text == '\t')) {
        (*text)++;
        (*length)--;
    }
}

static int label_is(const char *text, int64_t length, const char *want)
{
    return (int64_t)strlen(want) == length && !memcmp(text, want, (size_t)length);
}

static int reciprocal_to_cinter(int sg, rat (*out)[3])
{
    rat matrix[3][3], origin[3], inverse[3][3];
    if (geo_setting_to_cinter(sg, 1, 0, matrix, origin) || rat_inverse3((const rat (*)[3])matrix, inverse)) return KC_TABLES;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) out[row][col] = inverse[col][row];
    return 0;
}

static void apply_k_matrix(rat (*records)[3], int count, const int (*m)[3])
{
    for (int index = 0; index < count; index++) {
        rat t[3];
        for (int row = 0; row < 3; row++) {
            t[row] = rat_int(0);
            for (int col = 0; col < 3; col++) t[row] = rat_add(t[row], rat_mul(rat_int(m[row][col]), records[index][col]));
        }
        if (index > 0) {
            i128 scale = 0;
            for (int axis = 0; axis < 3; axis++)
                if (t[axis].n && t[axis].d == 1) scale = scale ? rat_gcd(scale, t[axis].n < 0 ? -t[axis].n : t[axis].n) : (t[axis].n < 0 ? -t[axis].n : t[axis].n);
            if (scale > 1)
                for (int axis = 0; axis < 3; axis++) t[axis] = rat_div(t[axis], rat_int(scale));
        }
        memcpy(records[index], t, sizeof t);
    }
}

/* fraction_rref over the nonzero rows; returns the pivot count. */
static int rref_rows(rat (*rows)[3], int count, int *pivots)
{
    int pivot_row = 0, found = 0;
    for (int col = 0; col < 3; col++) {
        int pivot = -1;
        for (int row = pivot_row; row < count && pivot < 0; row++)
            if (!rat_zero(rows[row][col])) pivot = row;
        if (pivot < 0) continue;
        rat swap[3];
        memcpy(swap, rows[pivot_row], sizeof swap);
        memcpy(rows[pivot_row], rows[pivot], sizeof swap);
        memcpy(rows[pivot], swap, sizeof swap);
        rat scale = rows[pivot_row][col];
        for (int axis = 0; axis < 3; axis++) rows[pivot_row][axis] = rat_div(rows[pivot_row][axis], scale);
        for (int row = 0; row < count; row++) {
            if (row == pivot_row || rat_zero(rows[row][col])) continue;
            rat factor = rows[row][col];
            for (int axis = 0; axis < 3; axis++) rows[row][axis] = rat_sub(rows[row][axis], rat_mul(factor, rows[pivot_row][axis]));
        }
        pivots[found++] = col;
        pivot_row++;
    }
    return found;
}

static void fixed_map(kmap *map, const char *expression)
{
    map->count = 0;
    for (int axis = 0; axis < 3; axis++) map->origin[axis] = rat_int(0);
    for (int axis = 0; axis < 3; axis++) {
        char c = expression[axis];
        if (c == '0') continue;
        if (c == 'h') {
            map->origin[axis] = rat_make(1, 2);
            continue;
        }
        map->names[map->count] = c;
        for (int k = 0; k < 3; k++) map->columns[map->count][k] = rat_int(k == axis);
        map->count++;
    }
}

static int name_rank(char c) { return c == 'a' ? 0 : c == 'b' ? 1 : 2; }

/* _fmt_display_k_vector before string formatting: the display records,
 * their parameter names (0 when unnamed), and a fixed expression ("abg",
 * "a0g", ... with 'h' for 1/2) when a whole-vector override applies. */
/* _display_records_for_star_arm: the lattice's display transform of Source
 * records (for lattices 9, 10, 13, 14 the pml reciprocal-to-cinter matrix). */
int kc_display_records(int sg, int lattice, const char *symbol, int64_t symbol_length, const char *label, int64_t label_length,
                       const rat (*records)[3], int count, rat (*d)[3])
{
    if (load()) return KC_TABLES;
    rat_overflow = 0;
    int n = count;
    for (int i = 0; i < n; i++) memcpy(d[i], records[i], sizeof d[i]);
    if (lattice == 9 || lattice == 10 || lattice == 13 || lattice == 14) {
        rat recip[3][3], t[4][3];
        if (reciprocal_to_cinter(sg, recip)) return KC_TABLES;
        for (int index = 0; index < n; index++) rat_row_mul(d[index], (const rat (*)[3])recip, t[index]);
        memcpy(d, t, sizeof(rat) * 3 * (size_t)n);
    } else if (lattice == 2) {
        for (int index = 0; index < n; index++) {
            rat t[3] = {d[index][1], d[index][2], d[index][0]};
            memcpy(d[index], t, sizeof t);
        }
    } else if (lattice == 3) {
        static const int m[3][3] = {{0, 1, -1}, {0, 1, 0}, {1, 0, 0}};
        apply_k_matrix(d, n, m);
        if (label_is(label, label_length, "V") && rat_eq(d[0][0], rat_make(-1, 2)) && rat_zero(d[0][1]) && rat_zero(d[0][2])) {
            d[0][0] = rat_make(-1, 2);
            d[0][1] = rat_make(1, 2);
            d[0][2] = rat_int(0);
        }
        if ((label_is(label, label_length, "Y") || label_is(label, label_length, "M")) && rat_zero(d[0][0]) &&
            rat_eq(d[0][1], rat_make(1, 2)))
            d[0][1] = rat_int(1);
    } else if (lattice == 5 && symbol_length && symbol[0] == 'C') {
        static const int m[3][3] = {{1, 1, 0}, {-1, 1, 0}, {0, 0, 1}};
        apply_k_matrix(d, n, m);
    } else if (lattice == 5 && symbol_length && symbol[0] == 'A') {
        static const int m[3][3] = {{0, 0, 1}, {-1, 1, 0}, {1, 1, 0}};
        apply_k_matrix(d, n, m);
        static const char *const flipped[] = {"Y", "R", "S", "T", "D", "H"};
        int flip = 0;
        for (int i = 0; i < 6; i++) flip |= label_is(label, label_length, flipped[i]);
        if (flip && rat_cmp(d[0][2], rat_int(0)) > 0) d[0][2] = rat_neg(d[0][2]);
    } else if (lattice == 6) {
        static const int m[3][3] = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
        apply_k_matrix(d, n, m);
    } else if (lattice == 7) {
        static const int m[3][3] = {{-1, 1, 1}, {1, -1, 1}, {1, 1, -1}};
        apply_k_matrix(d, n, m);
    }
    return rat_overflow ? KC_RANGE : 0;
}

int kc_display_form(int sg, int lattice, const char *symbol, int64_t symbol_length, const char *label, int64_t label_length,
                    const rat (*records)[3], int count, kc_form *form)
{
    if (load()) return KC_TABLES;
    int n = count;
    rat (*d)[3] = form->d;
    char *param_names = form->names;
    memset(param_names, 0, 4);
    form->fixed = 0;
    int status = kc_display_records(sg, lattice, symbol, symbol_length, label, label_length, records, count, d);
    if (status) return status;
    rat_overflow = 0;
    if (lattice == 9 || lattice == 10 || lattice == 13 || lattice == 14) {
        rat t[4][3];
        memcpy(t, d, sizeof(rat) * 3 * (size_t)n);
        rat rows[3][3];
        int nonzero = 0;
        for (int index = 1; index < n; index++)
            if (!rat_zero(t[index][0]) || !rat_zero(t[index][1]) || !rat_zero(t[index][2])) memcpy(rows[nonzero++], t[index], sizeof rows[0]);
        int pivots[3];
        int found = rref_rows(rows, nonzero, pivots);
        rat constant[3];
        memcpy(constant, t[0], sizeof constant);
        for (int i = 0; i < found; i++) {
            rat shift = constant[pivots[i]];
            if (rat_zero(shift)) continue;
            for (int axis = 0; axis < 3; axis++) constant[axis] = rat_sub(constant[axis], rat_mul(shift, rows[i][axis]));
        }
        memcpy(d[0], constant, sizeof constant);
        for (int i = 0; i < found; i++) {
            memcpy(d[i + 1], rows[i], sizeof rows[i]);
            param_names[i + 1] = NAMES[pivots[i]];
        }
        n = 1 + found;
    }
    form->count = n;
    if (rat_overflow) return KC_RANGE;
    if (n > 1 && rat_zero(d[0][0]) && rat_zero(d[0][1]) && rat_zero(d[0][2])) {
        if (label_is(label, label_length, "GP")) return form->fixed = "abg", 0;
        if (lattice == 3 && label_is(label, label_length, "B")) return form->fixed = "a0g", 0;
        if (n == 3) {
            int z2 = 1, z0 = 1, z1 = 1;
            for (int i = 1; i < 3; i++) {
                z2 &= rat_zero(d[i][2]);
                z0 &= rat_zero(d[i][0]);
                z1 &= rat_zero(d[i][1]);
            }
            if (z2) return form->fixed = "ab0", 0;
            if (z0) return form->fixed = "0bg", 0;
            if (z1) return form->fixed = "a0g", 0;
        }
    }
    if (n == 3) {
        if (rat_zero(d[0][0]) && rat_zero(d[0][1]) && rat_eq(d[0][2], rat_make(1, 2)) && rat_zero(d[1][2]) && rat_zero(d[2][2]))
            return form->fixed = "abh", 0;
        if (rat_eq(d[0][0], rat_make(1, 2)) && rat_zero(d[0][1]) && rat_zero(d[0][2]) && rat_zero(d[1][0]) && rat_zero(d[2][0]))
            return form->fixed = "hbg", 0;
    }
    int named = 0;
    for (int i = 1; i < 4; i++) named |= param_names[i] != 0;
    if (!named)
        for (int index = 1; index < n; index++)
            for (int axis = 0; axis < 3 && !param_names[index]; axis++)
                if (!rat_zero(d[index][axis])) param_names[index] = NAMES[axis];
    return 0;
}

static int kmap_build(int sg, int kslot, kmap *map)
{
    if (load()) return KC_TABLES;
    rat_overflow = 0;
    int lattice = geo_sg_lattice(sg);
    int64_t lslot = (int64_t)(lattice - 1) * 27 + kslot - 1;
    int dim = (int)section_int(&T.k_dim, lslot);
    map->dim = dim;
    map->count = 0;
    if (dim < 0 || dim > 3) return KC_TABLES;
    rat records[4][3];
    for (int index = 0; index < dim + 1; index++) {
        int64_t den = section_int(&T.k, lslot * 16 + 4 * index + 3);
        for (int axis = 0; axis < 3; axis++)
            records[index][axis] = den ? rat_make(section_int(&T.k, lslot * 16 + 4 * index + axis), den) : rat_int(0);
    }
    const char *symbol, *label;
    int64_t symbol_length, label_length;
    trimmed(&T.label_bc, sg - 1, &symbol, &symbol_length);
    trimmed(&T.k_label, lslot, &label, &label_length);
    kc_form form;
    int status = kc_display_form(sg, lattice, symbol, symbol_length, label, label_length, (const rat (*)[3])records, dim + 1, &form);
    if (status) return status;
    if (form.fixed) {
        fixed_map(map, form.fixed);
        return 0;
    }
    int n = form.count;
    int present[3] = {0, 0, 0};
    rat columns[3][3];
    for (int i = 0; i < 3; i++)
        for (int axis = 0; axis < 3; axis++) columns[i][axis] = rat_int(0);
    for (int axis = 0; axis < 3; axis++) {
        rat constant = form.d[0][axis];
        int parts = 0;
        for (int index = 1; index < n; index++) {
            rat value = form.d[index][axis];
            if (rat_zero(value)) continue;
            char name = form.names[index] ? form.names[index] : NAMES[index - 1 < 2 ? index - 1 : 2];
            int rank = name_rank(name);
            present[rank] = 1;
            columns[rank][axis] = rat_add(columns[rank][axis], value);
            parts++;
        }
        if (parts && constant.d == 1 && rat_cmp(constant, rat_int(1)) > 0) constant = rat_int(0);
        map->origin[axis] = constant;
    }
    for (int rank = 0; rank < 3; rank++) {
        if (!present[rank]) continue;
        map->names[map->count] = NAMES[rank];
        memcpy(map->columns[map->count], columns[rank], sizeof columns[rank]);
        map->count++;
    }
    return rat_overflow ? KC_RANGE : 0;
}

static int solve(const rat (*eq)[3], int width, const rat *rhs, rat *out)
{
    rat m[3][4];
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < width; col++) m[row][col] = eq[row][col];
        m[row][width] = rhs[row];
    }
    int pivot_row = 0, pivots[3], found = 0;
    for (int col = 0; col < width; col++) {
        int pivot = -1;
        for (int row = pivot_row; row < 3 && pivot < 0; row++)
            if (!rat_zero(m[row][col])) pivot = row;
        if (pivot < 0) continue;
        rat swap[4];
        memcpy(swap, m[pivot_row], sizeof swap);
        memcpy(m[pivot_row], m[pivot], sizeof swap);
        memcpy(m[pivot], swap, sizeof swap);
        rat factor = m[pivot_row][col];
        for (int i = 0; i <= width; i++) m[pivot_row][i] = rat_div(m[pivot_row][i], factor);
        for (int row = 0; row < 3; row++) {
            if (row == pivot_row || rat_zero(m[row][col])) continue;
            rat f = m[row][col];
            for (int i = 0; i <= width; i++) m[row][i] = rat_sub(m[row][i], rat_mul(f, m[pivot_row][i]));
        }
        pivots[found++] = col;
        pivot_row++;
    }
    for (int row = 0; row < 3; row++) {
        int zero = 1;
        for (int col = 0; col < width; col++) zero &= rat_zero(m[row][col]);
        if (zero && !rat_zero(m[row][width])) return 0;
    }
    if (found < width) return 0;
    for (int col = 0; col < width; col++) out[col] = rat_int(0);
    for (int i = 0; i < found; i++) out[pivots[i]] = m[i][width];
    return 1;
}

/* The visible K-family parameters of a PML reciprocal point: count >= 0,
 * -1 for none (Python None), < -1 for a failure. */
int case_kparams_of_pml(int gid, const rat *pml, rat *out)
{
    if (load()) return -KC_TABLES - 1;
    int sg = (int)section_int(&T.irr_sg, gid - 1), kslot = (int)section_int(&T.irr_k, gid - 1);
    int lattice = geo_sg_lattice(sg);
    int dim = (int)section_int(&T.k_dim, (int64_t)(lattice - 1) * 27 + kslot - 1);
    if (dim <= 0) return 0;
    kmap map;
    int status = kmap_build(sg, kslot, &map);
    if (status == KC_TABLES) return -KC_TABLES - 1;
    if (status || map.count != dim) return -1;
    rat recip[3][3], target[3];
    if (reciprocal_to_cinter(sg, recip)) return -KC_TABLES - 1;
    rat_overflow = 0;
    rat_row_mul(pml, (const rat (*)[3])recip, target);
    rat eq[3][3];
    for (int axis = 0; axis < 3; axis++)
        for (int p = 0; p < dim; p++) eq[axis][p] = map.columns[p][axis];
    int have = 0;
    rat best[3];
    for (int sx = -1; sx <= 1; sx++)
        for (int sy = -1; sy <= 1; sy++)
            for (int sz = -1; sz <= 1; sz++) {
                int shift[3] = {sx, sy, sz};
                rat rhs[3], solution[3];
                for (int axis = 0; axis < 3; axis++)
                    rhs[axis] = rat_sub(rat_add(target[axis], rat_int(shift[axis])), map.origin[axis]);
                if (!solve((const rat (*)[3])eq, dim, rhs, solution)) continue;
                int inside = 1;
                for (int p = 0; p < dim; p++) inside &= rat_cmp(solution[p], rat_int(0)) >= 0 && rat_cmp(solution[p], rat_int(1)) < 0;
                if (!inside) continue;
                int better = !have;
                for (int p = 0; p < dim && have && !better; p++) {
                    int c = rat_cmp(solution[p], best[p]);
                    if (c) {
                        better = c < 0;
                        break;
                    }
                }
                if (better) {
                    memcpy(best, solution, sizeof best);
                    have = 1;
                }
            }
    if (rat_overflow) return -KC_RANGE - 1;
    if (!have) return -1;
    memcpy(out, best, sizeof(rat) * (size_t)dim);
    return dim;
}

/* out: dim, count, names (char codes) x3, origin (n, d) x3, columns (n, d) x9. */
int k_coordinate_map_words(int sg, int kslot, int64_t *out)
{
    kmap map;
    int status = kmap_build(sg, kslot, &map);
    if (status) return status;
    out[0] = map.dim;
    out[1] = map.count;
    for (int i = 0; i < 3; i++) out[2 + i] = i < map.count ? map.names[i] : 0;
    for (int axis = 0; axis < 3; axis++) {
        out[5 + 2 * axis] = (int64_t)map.origin[axis].n;
        out[6 + 2 * axis] = (int64_t)map.origin[axis].d;
    }
    for (int i = 0; i < 3; i++)
        for (int axis = 0; axis < 3; axis++) {
            rat value = i < map.count ? map.columns[i][axis] : rat_int(0);
            out[11 + 6 * i + 2 * axis] = (int64_t)value.n;
            out[12 + 6 * i + 2 * axis] = (int64_t)value.d;
        }
    return 0;
}
