/* Mixed-radix basis-transform candidates of id_subgroup_: non-pivot columns take
 * trial values from a fixed list, pivot columns are back-substituted, and only
 * determinant +1 integer matrices qualify. */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "intmath.h"
#include "magnetic.h"
#include "rational.h"


static const int SEARCH_VALUES[11] = {0, 1, -1, 2, -2, 3, -3, 4, -4, 6, -6};

enum {
    ID_NARROW = -1,
    ID_OVERFLOW = -2,
};

/* rowop: column-stride-50 matrix of the given length; cols is fixed at 9.
 * Scans candidate indices start..limit and returns the first that yields a
 * determinant +1 matrix (written row-major into matrix), 0 when none, or a
 * negative error. */
static int64_t id_subgroup_next_candidate(const int64_t *rowop, int64_t length, int rows, int base,
                                          int64_t start, int64_t limit, int64_t *matrix)
{
    const int cols = 9;
    int pivot_by_row[64];
    int pivot_flags[10] = {0};
    if (rows > 64) return ID_NARROW;
    for (int row = 1; row <= rows; row++) {
        pivot_by_row[row - 1] = 0;
        for (int col = row; col <= cols; col++) {
            int64_t index = (row - 1) + (int64_t)(col - 1) * 50;
            if (index >= length) return ID_NARROW;
            if (rowop[index] != 0) {
                pivot_by_row[row - 1] = col;
                pivot_flags[col] = 1;
                break;
            }
        }
    }
    for (int64_t candidate = start; candidate <= limit; candidate++) {
        i128 values[10] = {0};
        int64_t residue = candidate - 1;
        for (int col = 1; col <= cols; col++)
            if (!pivot_flags[col]) {
                values[col] = SEARCH_VALUES[residue % base];
                residue /= base;
            }
        int ok = 1;
        for (int row = 1; row <= rows && ok; row++) {
            int pivot_col = pivot_by_row[row - 1];
            if (!pivot_col) continue;
            i128 pivot = rowop[(row - 1) + (int64_t)(pivot_col - 1) * 50];
            i128 total = 0;
            for (int col = pivot_col + 1; col <= cols; col++) {
                int64_t index = (row - 1) + (int64_t)(col - 1) * 50;
                if (index >= length) return ID_NARROW;
                total += (i128)rowop[index] * values[col];
            }
            if (total > ((i128)1 << 100) || total < -((i128)1 << 100)) return ID_OVERFLOW;
            if (pivot == 0 || total % pivot != 0) {
                ok = 0;
                break;
            }
            values[pivot_col] = -(total / pivot);
        }
        if (!ok) continue;
        i128 m[3][3];
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) m[row][col] = values[1 + row + col * 3];
        i128 det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                 - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                 + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        if (det != 1) continue;
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) {
                if (m[row][col] > INT64_MAX || m[row][col] < INT64_MIN) return ID_OVERFLOW;
                matrix[row * 3 + col] = (int64_t)m[row][col];
            }
        return candidate;
    }
    return 0;
}

/* Per-candidate origin solve of id_subgroup_ (generator translation
 * differences, the origin congruences and solve_eqs_mod_int_ with ifirst=1):
 * every operation is exact; i128 overflow or a nonconverging Smith form is a
 * range error. */

enum {
    OS_NONE = 0,
    OS_FOUND = 1,
    OS_FALLBACK = 2,
};

#define OS_MAX_GEN 16
#define OS_MAX_ROWS (3 * OS_MAX_GEN)

static int os_mul(i128 a, i128 b, i128 *out) { return __builtin_mul_overflow(a, b, out); }
static int os_add(i128 a, i128 b, i128 *out) { return __builtin_add_overflow(a, b, out); }
static int os_sub(i128 a, i128 b, i128 *out) { return __builtin_sub_overflow(a, b, out); }

static i128 os_abs(i128 a) { return a < 0 ? -a : a; }

static i128 os_gcd(i128 a, i128 b) { return gcd_i128(a, b); }

static i128 os_mod(i128 a, i128 m)
{
    i128 r = a % m;
    return r < 0 ? r + m : r;
}

/* Python a % m != 0 for m > 0 is the same as C a % m != 0. */

static int os_row_echelon(i128 *aug, int rows, int cols, int width)
{
    int total = cols + width;
    for (int pivot_index = 0; pivot_index < rows; pivot_index++) {
        for (;;) {
            if (pivot_index >= cols) break;
            int nonzero = 0, best = -1;
            i128 best_abs = 0;
            for (int row = pivot_index; row < rows; row++) {
                i128 value = aug[row * total + pivot_index];
                if (!value) continue;
                nonzero++;
                if (best < 0 || os_abs(value) < best_abs) {
                    best = row;
                    best_abs = os_abs(value);
                }
            }
            if (!nonzero) break;
            if (best != pivot_index)
                for (int col = 0; col < total; col++) {
                    i128 t = aug[best * total + col];
                    aug[best * total + col] = aug[pivot_index * total + col];
                    aug[pivot_index * total + col] = t;
                }
            if (aug[pivot_index * total + pivot_index] < 0)
                for (int col = 0; col < total; col++)
                    aug[pivot_index * total + col] = -aug[pivot_index * total + col];
            if (nonzero == 1) break;
            i128 pivot = aug[pivot_index * total + pivot_index];
            for (int row = pivot_index + 1; row < rows; row++) {
                i128 value = aug[row * total + pivot_index];
                if (!value) continue;
                i128 quotient = value / pivot;
                for (int col = pivot_index; col < total; col++) {
                    i128 product, next;
                    if (os_mul(quotient, aug[pivot_index * total + col], &product)) return 1;
                    if (os_sub(aug[row * total + col], product, &next)) return 1;
                    aug[row * total + col] = next;
                }
            }
        }
    }
    return 0;
}

static int os_has_offdiagonal(const i128 *m, int rows, int cols)
{
    for (int col = 0; col < cols; col++)
        for (int row = 0; row < rows; row++)
            if (row != col && m[row * cols + col]) return 1;
    return 0;
}

/* Smith form of rows x 3 matrix a; writes diagonal into a, left (rows x rows)
 * and right (3 x 3, already transposed as upstream returns it). */
static int os_smith(i128 *a, int rows, i128 *left, i128 *right)
{
    enum { COLS = 3 };
    i128 aug[OS_MAX_ROWS * (OS_MAX_ROWS + COLS)];
    i128 taug[COLS * (OS_MAX_ROWS + COLS)];
    i128 rt[COLS * COLS];
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < rows; j++) left[i * rows + j] = i == j;
    for (int i = 0; i < COLS; i++)
        for (int j = 0; j < COLS; j++) rt[i * COLS + j] = i == j;
    for (int iteration = 0; iteration < 100; iteration++) {
        int total = COLS + rows;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < COLS; j++) aug[i * total + j] = a[i * COLS + j];
            for (int j = 0; j < rows; j++) aug[i * total + COLS + j] = left[i * rows + j];
        }
        if (os_row_echelon(aug, rows, COLS, rows)) return 1;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < COLS; j++) a[i * COLS + j] = aug[i * total + j];
            for (int j = 0; j < rows; j++) left[i * rows + j] = aug[i * total + COLS + j];
        }
        if (!os_has_offdiagonal(a, rows, COLS)) {
            for (int i = 0; i < COLS; i++)
                for (int j = 0; j < COLS; j++) right[i * COLS + j] = rt[j * COLS + i];
            return 0;
        }
        int ttotal = rows + COLS;
        for (int i = 0; i < COLS; i++) {
            for (int j = 0; j < rows; j++) taug[i * ttotal + j] = a[j * COLS + i];
            for (int j = 0; j < COLS; j++) taug[i * ttotal + rows + j] = rt[i * COLS + j];
        }
        if (os_row_echelon(taug, COLS, rows, COLS)) return 1;
        for (int i = 0; i < COLS; i++) {
            for (int j = 0; j < rows; j++) a[j * COLS + i] = taug[i * ttotal + j];
            for (int j = 0; j < COLS; j++) rt[i * COLS + j] = taug[i * ttotal + rows + j];
        }
        if (!os_has_offdiagonal(a, rows, COLS)) {
            for (int i = 0; i < COLS; i++)
                for (int j = 0; j < COLS; j++) right[i * COLS + j] = rt[j * COLS + i];
            return 0;
        }
    }
    return 1;
}

/* Returns OS_FOUND (solution written), OS_NONE (no solution) or OS_FALLBACK. */
static int os_solve_candidate(const int64_t *basis, int k, const int64_t *inputs,
                              const int64_t *cands, const int64_t *mats, int64_t *solution)
{
    i128 m[9];
    for (int i = 0; i < 9; i++) m[i] = basis[i];
    i128 inv[9] = {
        m[4] * m[8] - m[5] * m[7], m[2] * m[7] - m[1] * m[8], m[1] * m[5] - m[2] * m[4],
        m[5] * m[6] - m[3] * m[8], m[0] * m[8] - m[2] * m[6], m[2] * m[3] - m[0] * m[5],
        m[3] * m[7] - m[4] * m[6], m[1] * m[6] - m[0] * m[7], m[0] * m[4] - m[1] * m[3],
    };
    i128 diff[OS_MAX_GEN][4];
    i128 modulus = 1;
    for (int g = 0; g < k; g++) {
        const int64_t *in = inputs + 4 * g, *cd = cands + 4 * g;
        i128 din = in[3], dc = cd[3], den, n[3];
        if (!din || !dc) return OS_FALLBACK;
        if (os_mul(din, dc, &den)) return OS_FALLBACK;
        for (int c = 0; c < 3; c++) {
            i128 t = 0, p, q;
            for (int r = 0; r < 3; r++) {
                if (os_mul((i128)in[r], inv[r * 3 + c], &p) || os_add(t, p, &t)) return OS_FALLBACK;
            }
            if (os_mul(t, dc, &p) || os_mul((i128)cd[c], din, &q) || os_sub(p, q, &n[c])) return OS_FALLBACK;
        }
        i128 g0 = os_gcd(os_gcd(os_gcd(n[0], n[1]), n[2]), den);
        if (!g0) return OS_FALLBACK;
        if (den < 0) g0 = -g0;
        den /= g0;
        for (int c = 0; c < 3; c++) diff[g][c] = os_mod(n[c] / g0, den);
        diff[g][3] = den;
        if (os_mul(modulus, den, &modulus)) return OS_FALLBACK;
    }
    int rows = 3 * k;
    i128 a[OS_MAX_ROWS * 3], rhs[OS_MAX_ROWS];
    for (int g = 0; g < k; g++)
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) a[(3 * g + r) * 3 + c] = (r == c) - (i128)mats[9 * g + c * 3 + r];
            i128 p;
            if (os_mul(modulus, diff[g][r], &p)) return OS_FALLBACK;
            rhs[3 * g + r] = p / diff[g][3];
        }
    i128 left[OS_MAX_ROWS * OS_MAX_ROWS], right[9];
    if (os_smith(a, rows, left, right)) return OS_FALLBACK;
    i128 trhs[OS_MAX_ROWS];
    for (int i = 0; i < rows; i++) {
        i128 t = 0, p;
        for (int j = 0; j < rows; j++)
            if (os_mul(left[i * rows + j], rhs[j], &p) || os_add(t, p, &t)) return OS_FALLBACK;
        trhs[i] = (t != 0 && t % modulus == 0) ? 0 : t;
    }
    int diag_count = rows < 3 ? rows : 3;
    i128 lcm = 1;
    int any = 0;
    for (int i = 0; i < diag_count; i++) {
        i128 d = os_abs(a[i * 3 + i]);
        if (!d) continue;
        any = 1;
        i128 g = os_gcd(lcm, d);
        if (os_mul(lcm / g, d, &lcm)) return OS_FALLBACK;
    }
    if (!any) {
        for (int i = 0; i < rows; i++)
            if (trhs[i] % modulus) return OS_NONE;
        solution[0] = solution[1] = solution[2] = 0;
        solution[3] = 1;
        return OS_FOUND;
    }
    i128 product;
    if (os_mul(lcm, modulus, &product)) return OS_FALLBACK;
    i128 y[3] = {0, 0, 0};
    for (int i = 0; i < diag_count; i++) {
        i128 d = a[i * 3 + i], p;
        if (!d) {
            if (trhs[i] % modulus) return OS_NONE;
            continue;
        }
        if (os_mul(product, trhs[i], &p)) return OS_FALLBACK;
        y[i] = (p / modulus) / d;
    }
    for (int i = diag_count; i < rows; i++)
        if (trhs[i] % modulus) return OS_NONE;
    i128 sol[3], common = product;
    for (int i = 0; i < 3; i++) {
        i128 t = 0, p;
        for (int j = 0; j < 3; j++)
            if (os_mul(right[i * 3 + j], y[j], &p) || os_add(t, p, &t)) return OS_FALLBACK;
        sol[i] = os_mod(t, product);
        common = os_gcd(common, sol[i]);
    }
    if (!common) common = 1;
    for (int i = 0; i < 3; i++) {
        i128 v = sol[i] / common;
        if (v > INT64_MAX || v < INT64_MIN) return OS_FALLBACK;
        solution[i] = (int64_t)v;
    }
    i128 den = product / common;
    if (den > INT64_MAX) return OS_FALLBACK;
    solution[3] = (int64_t)den;
    return OS_FOUND;
}

/* Scans candidates start..limit like id_subgroup_next_candidate and solves
 * the origin system for each.  Returns the candidate index and writes
 * status (OS_FOUND or OS_FALLBACK), matrix and solution; returns 0 when the
 * range holds no solvable candidate.  *evaluated counts determinant +1
 * candidates reached (the caller's lattice constraint depends on it). */
static int64_t id_subgroup_origin_scan(const int64_t *rowop, int64_t length, int rows, int base,
                                int64_t start, int64_t limit,
                                int k, const int64_t *inputs, const int64_t *cands,
                                const int64_t *mats, int64_t *status, int64_t *matrix,
                                int64_t *solution, int64_t *evaluated)
{
    *evaluated = 0;
    *status = OS_NONE;
    if (k <= 0 || k > OS_MAX_GEN) return ID_NARROW;
    while (start <= limit) {
        int64_t index = id_subgroup_next_candidate(rowop, length, rows, base, start, limit, matrix);
        if (index <= 0) return index;
        (*evaluated)++;
        int result = os_solve_candidate(matrix, k, inputs, cands, mats, solution);
        if (result != OS_NONE) {
            *status = result;
            return index;
        }
        start = index + 1;
    }
    return 0;
}

/* id_subgroup_'s rowop_ block for one generator mapping: the equations
 * C T - T R = 0 over the nine column-major unknowns of T (column stride 50,
 * as upstream), reduced with rowop_ rules (row factors divided out).
 * input, candidate: count x 9 row-major integer matrices.  out holds
 * length = max(rows, 1) + 8 * 50 values; returns the non-pivot column
 * count, or -1 when a value leaves int64. */
static int id_subgroup_rowop_block(int count, const int64_t *input, const int64_t *candidate, int64_t length, int64_t *out)
{
    const int leading = 50, cols = 9;
    int rows = count * 9;
    for (int64_t i = 0; i < length; i++) out[i] = 0;
#define AT(row, col) out[((row) - 1) + ((int64_t)(col) - 1) * leading]
    int equation = 1;
    for (int g = 0; g < count; g++) {
        const int64_t *r = input + 9 * g, *c = candidate + 9 * g;
        for (int col = 0; col < 3; col++)
            for (int row = 0; row < 3; row++) {
                int64_t coefficients[9] = {0};
                for (int inner = 0; inner < 3; inner++) {
                    coefficients[inner + col * 3] += c[row * 3 + inner];
                    coefficients[row + inner * 3] -= r[inner * 3 + col];
                }
                for (int variable = 1; variable <= 9; variable++) AT(equation, variable) = coefficients[variable - 1];
                equation++;
            }
    }
    int skipped = 0, pivot_row = 1;
    while (pivot_row <= rows) {
        int pivot_col = skipped + pivot_row;
        if (cols < pivot_col) break;
        int done = 0;
        while (AT(pivot_row, pivot_col) == 0) {
            int swap_row = 0;
            for (int row = pivot_row + 1; row <= rows && !swap_row; row++)
                if (AT(row, pivot_col) != 0) swap_row = row;
            if (swap_row) {
                for (int col = 1; col <= cols; col++) {
                    int64_t t = AT(pivot_row, col);
                    AT(pivot_row, col) = AT(swap_row, col);
                    AT(swap_row, col) = t;
                }
                break;
            }
            skipped++;
            pivot_col++;
            if (cols < pivot_col) {
                done = 1;
                break;
            }
        }
        if (done) break;
        int64_t pivot = AT(pivot_row, pivot_col);
        for (int row = 1; row <= rows; row++) {
            int64_t elim = AT(row, pivot_col);
            if (row == pivot_row || elim == 0) continue;
            for (int col = 1; col <= cols; col++) {
                if (col == pivot_col) continue;
                i128 value = (i128)pivot * AT(row, col) - (i128)AT(pivot_row, col) * elim;
                if (value > INT64_MAX || value < INT64_MIN) return -1;
                AT(row, col) = (int64_t)value;
            }
            AT(row, pivot_col) = 0;
            int64_t g = 0;
            for (int col = 1; col <= cols; col++) g = (int64_t)gcd_i128(g, AT(row, col));
            if (g > 1)
                for (int col = 1; col <= cols; col++) AT(row, col) /= g;
        }
        pivot_row++;
    }
    int flags[10] = {0}, free_count = 0;
    for (int row = 1; row <= rows; row++)
        for (int col = row; col <= cols; col++) {
            int64_t index = (row - 1) + (int64_t)(col - 1) * leading;
            if (index < length && out[index] != 0) {
                flags[col] = 1;
                break;
            }
        }
    for (int col = 1; col <= cols; col++) free_count += !flags[col];
#undef AT
    return free_count;
}

/* id_subgroup_ and id_subgroup_magnetic_: normalize the basis (nice_lattice2_),
 * bring the operations into it, and scan Source candidate groups for the
 * first generator mapping whose basis transform and origin solve. */

enum {
    ID_VALUE = -3,
    ID_TABLES = -4,
    ID_KEY = -5,
    ID_CAPACITY = -6,
    ID_ZERO = -7,
};

#define ID_MAX_INPUT 192
#define ID_MAX_GEN 5
#define ID_MAX_OP 144

static int id_bad;

static void id_fail(int status)
{
    if (!id_bad) id_bad = status;
}

/* matmlt_ on flat 3x3 buffers. */
static void matmlt(const i128 *left, const i128 *right, i128 *out)
{
    i128 t[9];
    for (int outer = 0; outer < 3; outer++)
        for (int group = 0; group < 3; group++) {
            i128 total = 0;
            for (int k = 0; k < 3; k++) total = rat_add_i(total, rat_mul_i(left[outer + 3 * k], right[3 * group + k]));
            t[outer + 3 * group] = total;
        }
    memcpy(out, t, sizeof t);
}

static void widen(const int64_t *in, i128 *out, int count)
{
    for (int i = 0; i < count; i++) out[i] = in[i];
}

/* reduc3_: divide out the common factor of the matrix and its denominator. */
static i128 reduce_matrix(i128 *m, i128 den)
{
    i128 g = den < 0 ? -den : den;
    for (int i = 0; i < 9; i++) g = gcd_i128(g, m[i]);
    if (g <= 0) g = 1;
    den /= g;
    for (int i = 0; i < 9; i++) m[i] /= g;
    if (den < 0) {
        den = -den;
        for (int i = 0; i < 9; i++) m[i] = -m[i];
    }
    return den;
}

static i128 floor_mod_i(i128 a, i128 m)
{
    i128 r = a % m;
    return r && (r < 0) != (m < 0) ? r + m : r;
}

/* A rational vector over its least common denominator, wrapped into [0, 1) when asked. */
static void record_of(const rat *v, int wrap, int64_t *out)
{
    i128 den = 1, n[3];
    for (int a = 0; a < 3; a++) den = rat_mul_i(den / gcd_i128(den, v[a].d), v[a].d);
    for (int a = 0; a < 3; a++) {
        n[a] = rat_mul_i(v[a].n, den / v[a].d);
        if (wrap) n[a] = floor_mod_i(n[a], den);
    }
    for (int a = 0; a < 3; a++) {
        if (n[a] > INT64_MAX || n[a] < INT64_MIN) rat_overflow = 1;
        out[a] = (int64_t)n[a];
    }
    if (den > INT64_MAX) rat_overflow = 1;
    out[3] = (int64_t)den;
}

static int rat_inverse_of(const i128 *flat, rat (*inv)[3])
{
    rat rows[3][3];
    for (int i = 0; i < 9; i++) rows[i / 3][i % 3] = rat_int(flat[i]);
    if (rat_inverse3((const rat (*)[3])rows, inv)) {
        id_fail(ID_VALUE);
        return 1;
    }
    return 0;
}

/* matinv_ as an integer numerator matrix over one reduced denominator. */
static i128 inverse_numerators(const i128 *flat, i128 *nums)
{
    rat inv[3][3];
    if (rat_inverse_of(flat, inv)) return 1;
    i128 den = 1;
    for (int i = 0; i < 9; i++) den = rat_mul_i(den / gcd_i128(den, inv[i / 3][i % 3].d), inv[i / 3][i % 3].d);
    i128 g = den;
    for (int i = 0; i < 9; i++) {
        nums[i] = rat_mul_i(inv[i / 3][i % 3].n, den / inv[i / 3][i % 3].d);
        g = gcd_i128(g, nums[i]);
    }
    if (g > 1) {
        for (int i = 0; i < 9; i++) nums[i] /= g;
        den /= g;
    }
    return den;
}

static int code_of(int magnetic, int op) { return (int)(magnetic ? mag_op_code(op) : geo_point_op_code(op)); }

static int product_of(int magnetic, int left, int right)
{
    return magnetic ? mag_point_product(left, right) : (int)geo_point_product(right, left);
}

static int spatial_of(int magnetic, int op) { return magnetic ? mag_nonmagnetic_op(op) : op; }

static void lattice_setting(int lattice, i128 *left, i128 *right, i128 *den)
{
    int64_t record[36];
    *den = geo_ml_setting(lattice, record);
    widen(record, left, 9);
    widen(record + 9, right, 9);
}

/* The candidate generator matrix S0 * ipoint_op * S1 of a space group. */
static int candidate_matrix(int sg, int op, i128 *out)
{
    i128 left[9], right[9], den, point[9];
    int64_t raw[9];
    lattice_setting(geo_sg_lattice(sg), left, right, &den);
    geo_point_op_matrix(op, raw);
    widen(raw, point, 9);
    matmlt(left, point, out);
    matmlt(out, right, out);
    if (reduce_matrix(out, den) != 1) {
        id_fail(ID_VALUE);
        return 1;
    }
    return 0;
}

static int nice_lattice2(int lattice, const i128 *basis, i128 *out)
{
    i128 left[9], right[9], den, b[9];
    lattice_setting(lattice, left, right, &den);
    matmlt(right, basis, b);
    for (int changed = 1; changed;) {
        changed = 0;
        for (int target = 0; target < 3; target++) {
            i128 length = 0;
            for (int a = 0; a < 3; a++) length = rat_add_i(length, rat_mul_i(b[3 * target + a], b[3 * target + a]));
            int replaced = 0;
            for (int other = 0; other < 3 && !replaced; other++) {
                if (other == target) continue;
                for (int sign = -1; sign <= 1 && !replaced; sign += 2) {
                    i128 c[3], square = 0;
                    for (int a = 0; a < 3; a++) {
                        c[a] = sign * b[3 * other + a] + b[3 * target + a];
                        square = rat_add_i(square, rat_mul_i(c[a], c[a]));
                    }
                    if (square < length) {
                        memcpy(b + 3 * target, c, sizeof c);
                        replaced = changed = 1;
                    }
                }
            }
        }
    }
    matmlt(left, b, out);
    if (reduce_matrix(out, den) != 1) {
        id_fail(ID_VALUE);
        return 1;
    }
    return 0;
}

typedef struct {
    int magnetic, count;
    int64_t active[ID_MAX_INPUT][5], local[ID_MAX_INPUT][5];
    int64_t matrices[ID_MAX_INPUT][9];
    int slot_by_op[ID_MAX_OP + 1];
} inputs_t;

static int prepare_inputs(int magnetic, int parent, const i128 *basis, int op_count, const int64_t *ops, inputs_t *in)
{
    in->magnetic = magnetic;
    in->count = 0;
    for (int i = 0; i < op_count; i++) {
        if (!ops[5 * i + 4]) continue;
        if (in->count == ID_MAX_INPUT) return ID_CAPACITY;
        memcpy(in->active[in->count++], ops + 5 * i, sizeof in->active[0]);
    }
    i128 nums[9], left[9], right[9], den0;
    i128 den = inverse_numerators(basis, nums);
    if (id_bad) return id_bad;
    lattice_setting(geo_sg_lattice(parent), left, right, &den0);
    memset(in->slot_by_op, 0, sizeof in->slot_by_op);
    for (int i = 0; i < in->count; i++) {
        const int64_t *record = in->active[i];
        if (!record[3]) return ID_ZERO;
        rat v[3];
        for (int c = 0; c < 3; c++) {
            i128 t = 0;
            for (int r = 0; r < 3; r++) t = rat_add_i(t, rat_mul_i(record[r], nums[r * 3 + c]));
            v[c] = rat_make(t, rat_mul_i(den, record[3]));
        }
        record_of(v, 1, in->local[i]);
        in->local[i][4] = record[4];
        int op = (int)record[4];
        if (op < 1 || op > ID_MAX_OP) return ID_TABLES;
        in->slot_by_op[op] = i + 1;
        int64_t raw[9];
        i128 point[9], m[9];
        geo_point_op_matrix(spatial_of(magnetic, op), raw);
        widen(raw, point, 9);
        matmlt(nums, left, m);
        matmlt(m, point, m);
        matmlt(m, right, m);
        matmlt(m, basis, m);
        if (reduce_matrix(m, rat_mul_i(den0, den)) != 1) return ID_VALUE;
        for (int k = 0; k < 9; k++) in->matrices[i][k] = (int64_t)m[k];
    }
    return 0;
}

/* Close a generator correspondence under point-operation products; on
 * success canonical (when given) maps each input slot to its candidate op. */
static int closure_accepts(const inputs_t *in, int gen_count, const int *selected, const int *candidate_ops,
                           int *canonical)
{
    int index_by_slot[ID_MAX_INPUT + 1] = {0}, index_by_op[ID_MAX_OP + 1] = {0};
    int slots[ID_MAX_OP + ID_MAX_INPUT], cops[ID_MAX_OP + ID_MAX_INPUT], mapped = 0;
    for (int g = 0; g < gen_count; g++) {
        if (index_by_slot[selected[g]]) return 0;
        mapped++;
        index_by_slot[selected[g]] = mapped;
        index_by_op[candidate_ops[g]] = mapped;
        slots[mapped - 1] = selected[g];
        cops[mapped - 1] = candidate_ops[g];
    }
    for (int changed = 1; changed;) {
        changed = 0;
        int fixed = mapped;
        for (int li = 0; li < fixed; li++)
            for (int ri = 0; ri < fixed; ri++) {
                int cp = product_of(in->magnetic, cops[li], cops[ri]);
                int ip = product_of(in->magnetic, (int)in->local[slots[li] - 1][4], (int)in->local[slots[ri] - 1][4]);
                if (cp < 1 || cp > ID_MAX_OP || ip < 1 || ip > ID_MAX_OP) {
                    id_fail(ID_TABLES);
                    return 0;
                }
                int islot = in->slot_by_op[ip];
                if (!islot) return 0;
                int existing_c = index_by_op[cp], existing_i = index_by_slot[islot];
                if (!existing_c) {
                    if (existing_i) return 0;
                    if (mapped == ID_MAX_OP + ID_MAX_INPUT) {
                        id_fail(ID_CAPACITY);
                        return 0;
                    }
                    mapped++;
                    index_by_op[cp] = mapped;
                    index_by_slot[islot] = mapped;
                    cops[mapped - 1] = cp;
                    slots[mapped - 1] = islot;
                    changed = 1;
                    continue;
                }
                if (slots[existing_c - 1] != islot || existing_i != existing_c) return 0;
            }
    }
    if (mapped < in->count) return 0;
    if (canonical) {
        for (int s = 1; s <= in->count; s++) canonical[s] = 0;
        for (int i = 0; i < mapped; i++) canonical[slots[i]] = cops[i];
    }
    return 1;
}

typedef struct {
    int group, sg, base, scan_ordinal, mapping_ordinal, gen_count;
    int64_t candidate_index, candidate_basis[9], solution[4];
    int64_t gens[ID_MAX_GEN][5];
    int match_count[ID_MAX_GEN], match[ID_MAX_GEN][ID_MAX_INPUT], selected[ID_MAX_GEN];
    i128 basis[9];
} found_t;

/* One candidate group: enumerate accepted mappings in mixed-radix order and
 * solve each.  Returns 1 when found, 0 otherwise; *reached reports whether a
 * determinant +1 basis was evaluated for some mapping. */
static int try_candidate(const inputs_t *in, found_t *f, int *reached)
{
    *reached = 0;
    int ops[ID_MAX_GEN];
    int64_t total = 1;
    for (int g = 0; g < f->gen_count; g++) {
        ops[g] = (int)f->gens[g][4];
        if (ops[g] < 1 || ops[g] > ID_MAX_OP) {
            id_fail(ID_TABLES);
            return 0;
        }
        int code = code_of(in->magnetic, ops[g]);
        f->match_count[g] = 0;
        for (int s = 0; s < in->count; s++)
            if (code_of(in->magnetic, (int)in->local[s][4]) == code) f->match[g][f->match_count[g]++] = s + 1;
        if (!f->match_count[g]) return 0;
        total *= f->match_count[g];
    }
    int have_matrices = 0;
    int64_t cmats[ID_MAX_GEN][9], cands[ID_MAX_GEN][4], inputs[ID_MAX_GEN][4], imats[ID_MAX_GEN][9];
    for (int g = 0; g < f->gen_count; g++) memcpy(cands[g], f->gens[g], sizeof cands[g]);
    f->mapping_ordinal = 0;
    for (int64_t ordinal = 0; ordinal < total; ordinal++) {
        int64_t value = ordinal;
        int selected[ID_MAX_GEN], duplicate = 0;
        for (int g = 0; g < f->gen_count; g++) {
            selected[g] = f->match[g][value % f->match_count[g]];
            value /= f->match_count[g];
            for (int h = 0; h < g; h++) duplicate |= selected[h] == selected[g];
        }
        if (duplicate || !closure_accepts(in, f->gen_count, selected, ops, 0)) {
            if (id_bad) return 0;
            continue;
        }
        f->mapping_ordinal++;
        if (!have_matrices) {
            for (int g = 0; g < f->gen_count; g++) {
                i128 m[9];
                if (candidate_matrix(f->sg, spatial_of(in->magnetic, ops[g]), m)) return 0;
                for (int k = 0; k < 9; k++) cmats[g][k] = (int64_t)m[k];
            }
            have_matrices = 1;
        }
        for (int g = 0; g < f->gen_count; g++) {
            memcpy(imats[g], in->matrices[selected[g] - 1], sizeof imats[g]);
            memcpy(inputs[g], in->local[selected[g] - 1], sizeof inputs[g]);
        }
        int64_t length = (f->gen_count * 9 > 1 ? f->gen_count * 9 : 1) + 8 * 50, rowop[45 + 400];
        int free_count = id_subgroup_rowop_block(f->gen_count, imats[0], cmats[0], length, rowop);
        if (free_count < 0) {
            id_fail(ID_OVERFLOW);
            return 0;
        }
        int64_t limit = 1, status, evaluated;
        for (int i = 0; i < free_count; i++) limit *= f->base;
        int64_t index = id_subgroup_origin_scan(rowop, length, f->gen_count * 9, f->base, 1, limit, f->gen_count,
                                                inputs[0], cands[0], cmats[0], &status, f->candidate_basis,
                                                f->solution, &evaluated);
        if (index < 0) {
            id_fail(index == ID_NARROW ? ID_NARROW : ID_OVERFLOW);
            return 0;
        }
        if (evaluated) *reached = 1;
        if (index > 0) {
            if (status == OS_FALLBACK) {
                id_fail(ID_OVERFLOW);
                return 0;
            }
            f->candidate_index = index;
            memcpy(f->selected, selected, sizeof selected);
            return 1;
        }
    }
    return 0;
}

static int write_result(const i128 *basis_flat, found_t *f, int64_t *out)
{
    i128 cb[9];
    widen(f->candidate_basis, cb, 9);
    matmlt(basis_flat, cb, f->basis);
    rat origin[3];
    for (int c = 0; c < 3; c++) {
        i128 t = 0;
        for (int r = 0; r < 3; r++) t = rat_add_i(t, rat_mul_i(f->solution[r], f->basis[3 * r + c]));
        origin[c] = rat_make(t, f->solution[3]);
    }
    for (int i = 0; i < 9; i++) {
        if (f->basis[i] > INT64_MAX || f->basis[i] < INT64_MIN) rat_overflow = 1;
        out[i] = (int64_t)f->basis[i];
    }
    record_of(origin, 0, out + 9);
    return 13;
}

static int finish(int status)
{
    if (!status && id_bad) status = id_bad;
    int table = geo_take_error();
    if (!status && table) status = table == GEO_RANGE ? ID_OVERFLOW : ID_TABLES;
    if (!status && rat_overflow) status = ID_OVERFLOW;
    return status;
}

static int start_search(int magnetic, int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag,
                        i128 *basis_flat, inputs_t *in)
{
    id_bad = 0;
    rat_overflow = 0;
    if (geo_load() || (magnetic && mag_tables())) return ID_TABLES;
    geo_take_error();
    i128 b[9];
    widen(basis, b, 9);
    if (flag >= 0) {
        if (nice_lattice2(geo_sg_lattice(parent), b, basis_flat)) return id_bad;
    } else
        memcpy(basis_flat, b, sizeof b);
    return prepare_inputs(magnetic, parent, basis_flat, op_count, ops, in);
}

static const int BASE_VALUE_COUNTS[4] = {3, 5, 7, 11};

/* out: subgroup, basis (9), origin (4).  Returns 1 when identified, 0 when
 * not, or a negative error. */
int id_subgroup_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t *out)
{
    static inputs_t in;
    static found_t f;
    i128 basis_flat[9];
    int status = start_search(0, parent, basis, op_count, ops, flag, basis_flat, &in);
    if (status) return finish(status);
    int64_t code = 0, group_count;
    for (int i = 0; i < in.count; i++) code += geo_point_op_code((int)in.local[i][4]);
    const int64_t *codes = geo_point_group_codes(&group_count);
    int point_group = 0;
    for (int64_t i = 0; i < group_count && !point_group; i++)
        if (codes[i] == code) point_group = (int)i + 1;
    if (!point_group) return finish(0);
    int first = flag ? abs(flag) : 1, last = flag ? abs(flag) : 230, constraint = 0;
    for (int b = 0; b < 4; b++)
        for (int sg = first; sg <= last; sg++) {
            if (geo_sg_point_group(sg) != point_group) continue;
            int lattice = geo_sg_lattice(sg);
            if (constraint && lattice != constraint) continue;
            f.group = 0;
            f.sg = sg;
            f.base = BASE_VALUE_COUNTS[b];
            f.gen_count = geo_pml_generators(sg, f.gens);
            if (f.gen_count < 0) return finish(ID_KEY);
            int reached;
            if (try_candidate(&in, &f, &reached)) {
                out[0] = sg;
                write_result(basis_flat, &f, out + 1);
                status = finish(0);
                return status ? status : 1;
            }
            if (id_bad) return finish(0);
            if (reached) constraint = lattice;
        }
    return finish(0);
}

/* The lifted full-Seitz correspondence of the accepted magnetic embedding;
 * returns the record count written (0 when the lift is incomplete). */
static int lift_correspondence(const inputs_t *in, const found_t *f, const int *canonical, int64_t *out)
{
    static int64_t records[ID_MAX_INPUT][5];
    int record_count = magnetic_records_iso(f->group, ID_MAX_INPUT, records[0]);
    if (record_count < 0) {
        id_fail(ID_TABLES);
        return 0;
    }
    i128 cb[9];
    widen(f->candidate_basis, cb, 9);
    rat inv[3][3], origin[3];
    if (rat_inverse_of(cb, inv)) return 0;
    for (int a = 0; a < 3; a++) origin[a] = rat_make(f->solution[a], f->solution[3]);
    for (int s = 1; s <= in->count; s++) {
        int op = canonical[s], found = -1;
        for (int r = 0; r < record_count; r++)
            if (records[r][4] == op) found = r;
        if (found < 0) return 0;
        i128 m[9];
        if (candidate_matrix(f->sg, mag_nonmagnetic_op(op), m)) return 0;
        rat matrix[3][3], local[3], rotated[3], transformed[3], translation[3];
        for (int i = 0; i < 9; i++) matrix[i / 3][i % 3] = rat_int(m[i]);
        rat_row_mul(origin, (const rat (*)[3])matrix, rotated);
        for (int a = 0; a < 3; a++) local[a] = rat_make(in->local[s - 1][a], in->local[s - 1][3]);
        rat_row_mul(local, (const rat (*)[3])inv, transformed);
        for (int a = 0; a < 3; a++) {
            translation[a] = rat_add(rat_sub(transformed[a], origin[a]), rotated[a]);
            if (!rat_eq(rat_mod1(translation[a]), rat_mod1(rat_make(records[found][a], records[found][3])))) return 0;
        }
        int64_t *row = out + 21 * (s - 1);
        row[0] = s;
        memcpy(row + 1, in->active[s - 1], sizeof(int64_t) * 5);
        memcpy(row + 6, in->local[s - 1], sizeof(int64_t) * 5);
        record_of(translation, 0, row + 11);
        row[15] = op;
        memcpy(row + 16, records[found], sizeof(int64_t) * 5);
    }
    return in->count;
}

/* Returns the length written to out, 0 when none is identified, or a
 * negative error.  out, in order: magnetic group, space group, basis (9), origin (4), base
 * value count, candidate scan ordinal, mapping ordinal, basis candidate
 * ordinal, candidate basis (9), candidate origin (4), generator count g,
 * match counts (g), match slots (their sum), selected slots (g), candidate
 * ops (g), input count n, n x (input op, canonical op, selected), lifted
 * count m, m x (slot, input, local, transformed, canonical records). */
int id_subgroup_magnetic_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag,
                                  int64_t capacity, int64_t *out)
{
    static inputs_t in;
    static found_t f;
    i128 basis_flat[9];
    int status = start_search(1, parent, basis, op_count, ops, flag, basis_flat, &in);
    if (status) return finish(status);
    int64_t code = 0, groups = mag_group_total();
    for (int i = 0; i < in.count; i++) code += mag_op_code((int)in.local[i][4]);
    int constraint = 0;
    for (int b = 0; b < 4; b++) {
        int ordinal = 0;
        for (int group = 1; group <= groups; group++) {
            if (mag_group_code(group) != code) continue;
            ordinal++;
            int sg = mag_ordinary_group(group), lattice = geo_sg_lattice(sg);
            if (constraint && lattice != constraint) continue;
            f.group = group;
            f.sg = sg;
            f.base = BASE_VALUE_COUNTS[b];
            f.scan_ordinal = ordinal;
            f.gen_count = mag_generators(group, f.gens);
            int reached;
            if (try_candidate(&in, &f, &reached)) {
                if (capacity < 64 + 6 * ID_MAX_INPUT + 21 * ID_MAX_INPUT) return finish(ID_CAPACITY);
                int canonical[ID_MAX_INPUT + 1];
                int cops[ID_MAX_GEN];
                for (int g = 0; g < f.gen_count; g++) cops[g] = (int)f.gens[g][4];
                closure_accepts(&in, f.gen_count, f.selected, cops, canonical);
                int64_t *p = out;
                *p++ = group;
                *p++ = sg;
                p += write_result(basis_flat, &f, p);
                *p++ = f.base;
                *p++ = f.scan_ordinal;
                *p++ = f.mapping_ordinal;
                *p++ = f.candidate_index;
                memcpy(p, f.candidate_basis, sizeof f.candidate_basis);
                p += 9;
                memcpy(p, f.solution, sizeof f.solution);
                p += 4;
                *p++ = f.gen_count;
                for (int g = 0; g < f.gen_count; g++) *p++ = f.match_count[g];
                for (int g = 0; g < f.gen_count; g++)
                    for (int i = 0; i < f.match_count[g]; i++) *p++ = f.match[g][i];
                for (int g = 0; g < f.gen_count; g++) *p++ = f.selected[g];
                for (int g = 0; g < f.gen_count; g++) *p++ = cops[g];
                *p++ = in.count;
                for (int s = 1; s <= in.count; s++) {
                    if (!canonical[s]) return finish(ID_KEY);
                    int selected = 0;
                    for (int g = 0; g < f.gen_count; g++) selected |= f.selected[g] == s;
                    *p++ = in.local[s - 1][4];
                    *p++ = canonical[s];
                    *p++ = selected;
                }
                int64_t *count_at = p++;
                *count_at = lift_correspondence(&in, &f, canonical, p);
                p += 21 * *count_at;
                status = finish(0);
                return status ? status : (int)(p - out);
            }
            if (id_bad) return finish(0);
            if (reached) constraint = lattice;
        }
    }
    return finish(0);
}
