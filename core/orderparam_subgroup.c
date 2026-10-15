/* orderparam_to_subgroup_ pieces: the newlat_order_ lattice search and the
 * per-row invariance filter of operation records, and the subduction
 * frequency of a little irrep over a subgroup.  Matrices come from
 * bw_record_matrix. */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "intmath.h"
#include "irreps.h"
#include "magnetic.h"


int bw_record_matrix(int dim, int arms, const int64_t *k, int subtract_generated,
                     int element_count, const int *elements, int position,
                     const double *sparse, const int *record, double *w);
int bw_find_position(int element_count, const int *elements, const int *record);
int bw_record_phase_coefficients(int dim, int arms, const int64_t *k, int subtract_generated, int element_count,
                                 const int *elements, const int *record, double *coefficients);
void bw_record_row(int dim, int arms, const double *coefficients, const double *sparse, int i, double *out);

enum {
    OP_OK = 0,
    OP_KEY = 101,   /* an operation record is not a generated operation */
    OP_VALUE = 102, /* the Source recurrence has no answer (phase, lattice, quotient) */
    OP_RANGE = 103, /* a value or buffer exceeds the fixed-width limits */
    OP_ZERO = 104,  /* a zero record denominator */
    OP_TABLES = 105, /* the Source little tables are inconsistent */
    OP_NO_TABLES = 106, /* the little irrep has no native matrix tables */
};

enum { BW_NO_POSITION = 1, BW_OVERFLOW = 5 };

typedef struct {
    int dim, arms;
    const int64_t *k;
    int subtract;
    int element_count;
    const int *elements;
    const double *sparse;
    const int *available;
} op_context;

/* While a caller holds one irrep context fixed (op_translation_cache_begin
 * .. _end), the phase coefficients of the pure translations newlat_search
 * tries are kept and reused across rows; they depend on the context and the
 * translation only. */
#define OP_CACHE_SLOTS 65536

static struct {
    int active;
    int keys[OP_CACHE_SLOTS][5];
    unsigned char used[OP_CACHE_SLOTS];
    int status[OP_CACHE_SLOTS];
    double *coefficients[OP_CACHE_SLOTS];
} tcache;

void op_translation_cache_begin(void)
{
    tcache.active = 1;
}

void op_translation_cache_end(void)
{
    for (int i = 0; i < OP_CACHE_SLOTS; i++) {
        if (tcache.used[i]) free(tcache.coefficients[i]);
        tcache.used[i] = 0;
        tcache.coefficients[i] = 0;
    }
    tcache.active = 0;
}

/* op_matrix's position and phase steps for a record: the position and the
 * phase coefficients, or the op_matrix status. */
static int record_phases(const op_context *c, const int *record, int *position, const double **coefficients,
                         double *scratch)
{
    *position = bw_find_position(c->element_count, c->elements, record);
    if (!*position) return OP_KEY;
    if (!c->available[*position - 1]) return OP_VALUE;
    size_t slot = 0;
    int cache = tcache.active;
    if (cache) {
        uint64_t h = (uint64_t)(uint32_t)record[0] * 0x9E3779B97F4A7C15ull ^
                     (uint64_t)(uint32_t)record[1] * 0xC2B2AE3D27D4EB4Full ^
                     (uint64_t)(uint32_t)record[2] * 0x165667B19E3779F9ull ^ (uint64_t)(uint32_t)record[4] ^
                     (uint64_t)(uint32_t)record[3] * 0x27D4EB2F165667C5ull;
        slot = (size_t)(h ^ (h >> 29)) & (OP_CACHE_SLOTS - 1);
        for (int probe = 0; probe < 16; probe++, slot = (slot + 1) & (OP_CACHE_SLOTS - 1)) {
            if (!tcache.used[slot]) break;
            const int *key = tcache.keys[slot];
            if (key[0] == record[0] && key[1] == record[1] && key[2] == record[2] && key[3] == record[3] &&
                key[4] == record[4]) {
                *coefficients = tcache.coefficients[slot] ? tcache.coefficients[slot] : scratch;
                return tcache.status[slot];
            }
        }
        cache = !tcache.used[slot];
    }
    int status = bw_record_phase_coefficients(c->dim, c->arms, c->k, c->subtract, c->element_count, c->elements,
                                              record, scratch);
    if (status) status = status == BW_NO_POSITION ? OP_KEY : status == BW_OVERFLOW ? OP_RANGE : OP_VALUE;
    *coefficients = scratch;
    if (!cache) return status;
    double *copy = 0;
    if (!status) {
        copy = malloc(sizeof(double) * (2 * (size_t)c->arms + 1));
        if (!copy) return status;
        memcpy(copy, scratch, sizeof(double) * 2 * (size_t)c->arms);
    }
    tcache.used[slot] = 1;
    tcache.keys[slot][0] = record[0];
    tcache.keys[slot][1] = record[1];
    tcache.keys[slot][2] = record[2];
    tcache.keys[slot][3] = record[3];
    tcache.keys[slot][4] = record[4];
    tcache.status[slot] = status;
    tcache.coefficients[slot] = copy;
    return status;
}

/* While set, the invariance test of a record is decided by this hook (an
 * exact one for rows known over a cyclotomic field) after the position and
 * phase checks. */
int (*op_invariant_hook)(void *ctx, const int *record, int position, double sign, const double *row, int *result);
void *op_invariant_ctx;

/* op_matrix then op_invariant, building W one row at a time and stopping at
 * the first row that fails (the rows it builds are op_matrix's bit for bit). */
static int op_record_invariant(const op_context *c, const int *record, double sign, const double *row, int *result)
{
    int position;
    double scratch[96], w[48];
    const double *coefficients;
    int status = record_phases(c, record, &position, &coefficients, scratch);
    if (status) return status;
    if (op_invariant_hook) return op_invariant_hook(op_invariant_ctx, record, position, sign, row, result);
    int dim = c->dim;
    const double *sparse = c->sparse + (size_t)(position - 1) * dim * dim;
    for (int col = 0; col < dim; col++) {
        bw_record_row(dim, c->arms, coefficients, sparse, col, w);
        double transformed = 0.0;
        for (int source = 0; source < dim; source++) {
            double value = w[source];
            if (sign < 0) value = -value;
            transformed += value * row[source];
        }
        if (1e-6 < fabs(transformed - row[col])) {
            *result = 0;
            return OP_OK;
        }
    }
    *result = 1;
    return OP_OK;
}

static int op_allowed(const op_context *c, int identity_op, const double *row,
                      int has_contains, const int64_t *contains, int64_t contains_den,
                      const int v[3], int *result)
{
    if (has_contains)
        for (int col = 0; col < 3; col++) {
            i128 t = (i128)v[0] * contains[col] + (i128)v[1] * contains[3 + col]
                   + (i128)v[2] * contains[6 + col];
            if (t % contains_den) {
                *result = 0;
                return OP_OK;
            }
        }
    int record[5] = {v[0], v[1], v[2], 1, identity_op};
    return op_record_invariant(c, record, 1.0, row, result);
}

static int newlat_search(int dim, int arms, const int64_t *k, int subtract,
                             int element_count, const int *elements, const double *sparse,
                             const int *available, int identity_op, const double *row,
                             int has_contains, const int64_t *contains, int64_t contains_den,
                             int max_iter, int64_t *rows_out)
{
    op_context c = {dim, arms, k, subtract, element_count, elements, sparse, available};
    if (dim <= 0 || dim > 48 || max_iter <= 0 || max_iter > 1024) return OP_RANGE;
    if (has_contains && contains_den <= 0) return OP_VALUE;
    int periods[3] = {0, 0, 0};
    for (int axis = 0; axis < 3; axis++) {
        int trial[3] = {0, 0, 0};
        for (int value = 1; value <= max_iter; value++) {
            trial[axis] = value;
            int ok, status = op_allowed(&c, identity_op, row, has_contains, contains, contains_den, trial, &ok);
            if (status) return status;
            if (ok) {
                periods[axis] = value;
                break;
            }
        }
        if (!periods[axis]) return OP_VALUE;
    }
    int64_t guard = (int64_t)max_iter * max_iter * max_iter;
    for (int axis = 0; axis < 3; axis++) {
        int trial[3] = {0, 0, 0};
        trial[axis] = 1;
        int found = 0;
        for (int64_t step = 0; step < guard; step++) {
            int ok, status = op_allowed(&c, identity_op, row, has_contains, contains, contains_den, trial, &ok);
            if (status) return status;
            if (ok) {
                for (int i = 0; i < 3; i++) rows_out[axis * 3 + i] = trial[i];
                found = 1;
                break;
            }
            trial[0] += 1;
            if (trial[0] != periods[0] || axis == 0) continue;
            trial[0] = 0;
            trial[1] += 1;
            if (trial[1] != periods[1] || axis == 1) continue;
            trial[1] = 0;
            trial[2] += 1;
            if (trial[2] != periods[2] || axis == 2) continue;
            return OP_VALUE;
        }
        if (!found) return OP_VALUE;
    }
    return OP_OK;
}

static i128 op_gcd(i128 a, i128 b) { return gcd_i128(a, b); }

/* _vadd_fraction_operation */
static int64_t gcd_s64(int64_t a, int64_t b) { return (int64_t)gcd_u64((uint64_t)(a < 0 ? -a : a), (uint64_t)(b < 0 ? -b : b)); }

int op_vadd(const int64_t *fraction, const int64_t *operation, int64_t *out)
{
    int small = 1;
    for (int i = 0; i < 4 && small; i++)
        small = fraction[i] > -(1 << 30) && fraction[i] < (1 << 30) && operation[i] > -(1 << 30) && operation[i] < (1 << 30);
    if (small) {
        /* the i128 recurrence below, exactly, while every product fits 64 bits */
        int64_t fd = fraction[3], od = operation[3];
        if (!fd || !od) return OP_ZERO;
        int64_t common = fd * od, n[3], den = 1;
        int64_t abs_common = common < 0 ? -common : common;
        for (int axis = 0; axis < 3; axis++) {
            n[axis] = fraction[axis] * od + operation[axis] * fd;
            int64_t g = gcd_s64(n[axis], common);
            int64_t part = abs_common / g;
            den = den / gcd_s64(den, part) * part;
        }
        for (int axis = 0; axis < 3; axis++) {
            __int128 value = (__int128)n[axis] * den;
            if (value % common) return OP_VALUE;
            value /= common;
            if (value > INT32_MAX || value < INT32_MIN) return OP_RANGE;
            out[axis] = (int64_t)value;
        }
        if (den > INT32_MAX) return OP_RANGE;
        out[3] = den;
        out[4] = operation[4];
        return OP_OK;
    }
    i128 fd = fraction[3], od = operation[3];
    if (!fd || !od) return OP_ZERO;
    i128 common = fd * od, n[3], den = 1;
    i128 abs_common = common < 0 ? -common : common;
    for (int axis = 0; axis < 3; axis++) {
        n[axis] = (i128)fraction[axis] * od + (i128)operation[axis] * fd;
        i128 part = abs_common / op_gcd(n[axis], common);
        den = den / op_gcd(den, part) * part;
    }
    for (int axis = 0; axis < 3; axis++) {
        i128 value = n[axis] * den;
        if (value % common) return OP_VALUE;
        value /= common;
        if (value > INT32_MAX || value < INT32_MIN) return OP_RANGE;
        out[axis] = (int64_t)value;
    }
    if (den > INT32_MAX) return OP_RANGE;
    out[3] = (int64_t)den;
    out[4] = operation[4];
    return OP_OK;
}

static int op_filter(const op_context *c, const double *row, int candidate_count,
                     const int64_t *candidates, const int *matrix_ops, const int *signs,
                     int fraction_count, const int64_t *fractions, int64_t *kept, int *kept_from,
                     int *kept_count)
{
    int count = 0;
    for (int i = 0; i < candidate_count; i++) {
        for (int f = 0; f < fraction_count; f++) {
            int64_t shifted[5];
            int status = op_vadd(fractions + 4 * f, candidates + 5 * i, shifted);
            if (status) return status;
            int record[5] = {(int)shifted[0], (int)shifted[1], (int)shifted[2], (int)shifted[3],
                             matrix_ops[i]};
            int invariant;
            status = op_record_invariant(c, record, signs[i] < 0 ? -1.0 : 1.0, row, &invariant);
            if (status) return status;
            if (invariant) {
                for (int j = 0; j < 5; j++) kept[5 * count + j] = shifted[j];
                kept_from[count] = i;
                count++;
                break;
            }
        }
    }
    *kept_count = count;
    return OP_OK;
}

static i128 op_abs(i128 a) { return a < 0 ? -a : a; }

static i128 op_det(const i128 *m)
{
    return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
         + m[2] * (m[3] * m[7] - m[4] * m[6]);
}

static void op_adjugate(const i128 *m, i128 *a)
{
    a[0] = m[4] * m[8] - m[5] * m[7];
    a[1] = m[2] * m[7] - m[1] * m[8];
    a[2] = m[1] * m[5] - m[2] * m[4];
    a[3] = m[5] * m[6] - m[3] * m[8];
    a[4] = m[0] * m[8] - m[2] * m[6];
    a[5] = m[2] * m[3] - m[0] * m[5];
    a[6] = m[3] * m[7] - m[4] * m[6];
    a[7] = m[1] * m[6] - m[0] * m[7];
    a[8] = m[0] * m[4] - m[1] * m[3];
}

/* SourceData._reduce_basis_rows */
static void op_reduce_rows(int64_t *rows)
{
    for (int pass = 0; pass < 3; pass++)
        for (int r = 0; r < 3; r++) {
            int64_t size = 0;
            for (int a = 0; a < 3; a++) size = llabs(rows[r * 3 + a]) > size ? llabs(rows[r * 3 + a]) : size;
            for (int o = 0; o < 3; o++) {
                if (o == r) continue;
                for (int sign = -1; sign <= 1; sign += 2) {
                    int64_t cand[3], csize = 0;
                    for (int a = 0; a < 3; a++) {
                        cand[a] = rows[r * 3 + a] + sign * rows[o * 3 + a];
                        if (llabs(cand[a]) > csize) csize = llabs(cand[a]);
                    }
                    if (csize < size) {
                        for (int a = 0; a < 3; a++) rows[r * 3 + a] = cand[a];
                        size = csize;
                    }
                }
            }
        }
}

static int op_cmp3(const void *left, const void *right)
{
    const int64_t *a = left, *b = right;
    for (int i = 0; i < 3; i++)
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

/* get_new_fractionals_: the x in [0, m)^3, sorted, with x M = 0 (mod m), as
 * (x M / m, 1).  The solutions are the row lattice of m M^-1 plus m Z^3; its
 * Hermite form enumerates them directly. */
int op_new_fractionals(const int64_t *matrix, int64_t *out, int capacity, int *count)
{
    i128 m[9], adj[9];
    for (int i = 0; i < 9; i++) m[i] = matrix[i];
    i128 det = op_det(m);
    if (!det) return OP_VALUE;
    i128 abs_det = det < 0 ? -det : det;
    op_adjugate(m, adj);
    i128 factor = abs_det;
    for (int i = 0; i < 9; i++) factor = op_gcd(factor, adj[i]);
    i128 modulus = abs_det / factor;
    if (abs_det == 1) {
        if (capacity < 1) return OP_RANGE;
        out[0] = out[1] = out[2] = 0;
        out[3] = 1;
        *count = 1;
        return OP_OK;
    }
    if (abs_det > capacity || modulus > ((i128)1 << 40)) return OP_RANGE;
    i128 g[6][3];
    for (int r = 0; r < 3; r++)
        for (int col = 0; col < 3; col++) {
            i128 value = modulus * adj[r * 3 + col];
            if (value % det) return OP_VALUE;
            g[r][col] = value / det;
            g[3 + r][col] = r == col ? modulus : 0;
        }
    int rows = 6;
    i128 h[3][3] = {{0}};
    for (int col = 0, top = 0; col < 3; col++, top++) {
        for (;;) {
            int pivot = -1;
            for (int r = top; r < rows; r++)
                if (g[r][col] && (pivot < 0 || op_abs(g[r][col]) < op_abs(g[pivot][col]))) pivot = r;
            if (pivot < 0) return OP_VALUE;
            if (pivot != top)
                for (int k = 0; k < 3; k++) {
                    i128 t = g[pivot][k];
                    g[pivot][k] = g[top][k];
                    g[top][k] = t;
                }
            int done = 1;
            for (int r = top + 1; r < rows; r++) {
                if (!g[r][col]) continue;
                i128 q = g[r][col] / g[top][col];
                for (int k = 0; k < 3; k++) g[r][k] -= q * g[top][k];
                if (g[r][col]) done = 0;
            }
            if (done) break;
        }
        if (g[top][col] < 0)
            for (int k = 0; k < 3; k++) g[top][k] = -g[top][k];
        for (int k = 0; k < 3; k++) h[col][k] = g[top][k];
    }
    i128 n0 = modulus / h[0][0], n1 = modulus / h[1][1], n2 = modulus / h[2][2];
    if (modulus % h[0][0] || modulus % h[1][1] || modulus % h[2][2] || n0 * n1 * n2 != abs_det)
        return OP_VALUE;
    int total = 0;
    for (i128 i = 0; i < n0; i++)
        for (i128 j = 0; j < n1; j++)
            for (i128 k = 0; k < n2; k++) {
                for (int col = 0; col < 3; col++) {
                    i128 value = i * h[0][col] + j * h[1][col] + k * h[2][col];
                    value %= modulus;
                    if (value < 0) value += modulus;
                    out[4 * total + col] = (int64_t)value;
                }
                total++;
            }
    int64_t (*items)[4] = (int64_t (*)[4])out;
    qsort(items, (size_t)total, sizeof items[0], op_cmp3);
    for (int t = 1; t < total; t++)
        if (!op_cmp3(items[t - 1], items[t])) return OP_VALUE;
    for (int t = 0; t < total; t++) {
        i128 x = items[t][0], y = items[t][1], z = items[t][2];
        for (int col = 0; col < 3; col++) {
            i128 value = x * m[col] + y * m[3 + col] + z * m[6 + col];
            if (value % modulus) return OP_VALUE;
            i128 q = value / modulus;
            if (q > INT32_MAX || q < INT32_MIN) return OP_VALUE;
            items[t][col] = (int64_t)q;
        }
        items[t][3] = 1;
    }
    *count = total;
    return OP_OK;
}

#define OP_MAX_FRACTIONS 4096
#define OP_MAX_CANDIDATES 4096

typedef struct {
    int64_t fractions[OP_MAX_FRACTIONS * 4], raw[OP_MAX_FRACTIONS * 4];
    int64_t current[OP_MAX_CANDIDATES * 5], kept[OP_MAX_CANDIDATES * 5];
    int ops[OP_MAX_CANDIDATES], signs[OP_MAX_CANDIDATES], from[OP_MAX_CANDIDATES];
} op_work;

/* The orderparam_to_subgroup_ recurrence (SourceData._orderparam_to_subgroup_rows):
 * values holds row_count rows at stride 48.  Writes the basis and the kept
 * candidate records. */
int subgroup_rows_in_context(int dim, int arms, const int64_t *k, int subtract,
                                int element_count, const int *elements, const double *sparse,
                                const int *available, int identity_op,
                                const double *values, int64_t value_count, int row_count,
                                int candidate_count, const int64_t *candidates,
                                const int *candidate_matrix_ops, const int *candidate_signs,
                                int max_iter, int64_t *basis_out, int64_t *kept_out, int *kept_count)
{
    op_context c = {dim, arms, k, subtract, element_count, elements, sparse, available};
    if (dim <= 0 || dim > 48 || row_count < 0 || candidate_count > OP_MAX_CANDIDATES) return OP_RANGE;
    if (row_count > 0 && (int64_t)(row_count - 1) * 48 + dim > value_count) return OP_RANGE;
    static op_work *work;
    if (!work) work = malloc(sizeof *work);
    if (!work) return OP_RANGE;
    int status = OP_OK;
    int remaining = candidate_count;
    for (int i = 0; i < remaining; i++) {
        for (int j = 0; j < 5; j++) work->current[5 * i + j] = candidates[5 * i + j];
        work->ops[i] = candidate_matrix_ops[i];
        work->signs[i] = candidate_signs[i];
    }
    int64_t basis[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (int r = 0; r < row_count; r++) {
        const double *row = values + (int64_t)r * 48;
        int fraction_count = 0;
        int64_t next[9];
        if (r == 0) {
            status = newlat_search(dim, arms, k, subtract, element_count, elements, sparse,
                                              available, identity_op, row, 0, 0, 1, max_iter, next);
            if (status) goto done;
            op_reduce_rows(next);
            status = op_new_fractionals(next, work->fractions, OP_MAX_FRACTIONS, &fraction_count);
            if (status) goto done;
        } else {
            i128 b[9], adj[9];
            for (int i = 0; i < 9; i++) b[i] = basis[i];
            i128 det = op_det(b);
            if (!det) {
                status = OP_VALUE;
                goto done;
            }
            op_adjugate(b, adj);
            int64_t contains[9];
            for (int i = 0; i < 9; i++) contains[i] = (int64_t)adj[i];
            int64_t abs_det = (int64_t)(det < 0 ? -det : det);
            status = newlat_search(dim, arms, k, subtract, element_count, elements, sparse,
                                              available, identity_op, row, 1, contains, abs_det, max_iter, next);
            if (status) goto done;
            op_reduce_rows(next);
            int64_t transform[9];
            for (int row_index = 0; row_index < 3; row_index++)
                for (int col = 0; col < 3; col++) {
                    i128 value = 0;
                    for (int t = 0; t < 3; t++) value += (i128)next[row_index * 3 + t] * adj[t * 3 + col];
                    if (value % det) {
                        status = OP_VALUE;
                        goto done;
                    }
                    value /= det;
                    if (value > INT32_MAX || value < INT32_MIN) {
                        status = OP_RANGE;
                        goto done;
                    }
                    transform[row_index * 3 + col] = (int64_t)value;
                }
            status = op_new_fractionals(transform, work->raw, OP_MAX_FRACTIONS, &fraction_count);
            if (status) goto done;
            for (int f = 0; f < fraction_count; f++) {
                const int64_t *item = work->raw + 4 * f;
                i128 v[3], g = 0;
                for (int a = 0; a < 3; a++) {
                    v[a] = (i128)item[0] * basis[a] + (i128)item[1] * basis[3 + a] + (i128)item[2] * basis[6 + a];
                    g = op_gcd(g, v[a]);
                }
                g = op_gcd(g, item[3]);
                if (!g) g = 1;
                for (int a = 0; a < 3; a++) work->fractions[4 * f + a] = (int64_t)(v[a] / g);
                work->fractions[4 * f + 3] = (int64_t)(item[3] / g);
            }
        }
        for (int i = 0; i < 9; i++) basis[i] = next[i];
        int count = 0;
        status = op_filter(&c, row, remaining, work->current, work->ops, work->signs, fraction_count,
                           work->fractions, work->kept, work->from, &count);
        if (status) goto done;
        for (int i = 0; i < count; i++) {
            for (int j = 0; j < 5; j++) work->current[5 * i + j] = work->kept[5 * i + j];
            work->ops[i] = work->ops[work->from[i]];
            work->signs[i] = work->signs[work->from[i]];
        }
        remaining = count;
    }
    for (int i = 0; i < 9; i++) basis_out[i] = basis[i];
    for (int i = 0; i < remaining; i++)
        for (int j = 0; j < 5; j++) kept_out[5 * i + j] = work->current[5 * i + j];
    *kept_count = remaining;
done:
    return status;
}

static int context_status(int status)
{
    return status == IR_RANGE ? OP_RANGE : status == IR_NO_TABLES ? OP_NO_TABLES : OP_TABLES;
}

/* newlat_order_ for one row of a little irrep at a k parameter. */
int orderparam_newlat_search(int gid, int param_count, const int64_t *kparam, const double *row, int has_contains,
                             const int64_t *contains, int64_t contains_den, int max_iter, int64_t *rows_out)
{
    irrep_context c;
    int status = irrep_context_for(gid, param_count, kparam, &c);
    if (status) return context_status(status);
    return newlat_search(c.dim, c.arms, c.k, c.subtract, c.element_count, c.elements, c.sparse, c.available,
                         c.identity_op, row, has_contains, contains, contains_den, max_iter, rows_out);
}

/* orderparam_to_subgroup_ of a little irrep at a k parameter. */
int orderparam_to_subgroup_rows(int gid, int param_count, const int64_t *kparam, const double *values,
                                int64_t value_count, int row_count, int candidate_count, const int64_t *candidates,
                                const int *candidate_matrix_ops, const int *candidate_signs, int max_iter,
                                int64_t *basis_out, int64_t *kept_out, int *kept_count)
{
    irrep_context c;
    int status = irrep_context_for(gid, param_count, kparam, &c);
    if (status) return context_status(status);
    return subgroup_rows_in_context(c.dim, c.arms, c.k, c.subtract, c.element_count, c.elements, c.sparse,
                                    c.available, c.identity_op, values, value_count, row_count, candidate_count,
                                    candidates, candidate_matrix_ops, candidate_signs, max_iter, basis_out, kept_out,
                                    kept_count);
}

int bridge_weight_matrices(int dim, int arms, const int64_t *k, int subtract_generated, int element_count,
                           const int *elements, int op_limit, const int *canonical_position, const double *sparse,
                           const int *available, int record_count, const int *records, double *out);

/* get_irreps_ matrices of a little irrep at a k parameter for many records
 * (record_count x dim x dim); returns the bridge-weight status, or
 * 100 + an orderparam status for the context. */
int irrep_matrix_stack(int gid, int param_count, const int64_t *kparam, int record_count, const int *records,
                       int capacity, double *out)
{
    irrep_context c;
    int status = irrep_context_for(gid, param_count, kparam, &c);
    if (status) return 100 + context_status(status);
    if ((int64_t)record_count * c.dim * c.dim > capacity) return 100 + OP_RANGE;
    const int none = 0;
    return bridge_weight_matrices(c.dim, c.arms, c.k, c.subtract, c.element_count, c.elements, 0, &none, c.sparse,
                                  c.available, record_count, records, out);
}

/* numpy's pairwise summation of n strided doubles (ndarray.trace). */
static double numpy_sum(const double *a, int n, int stride)
{
    if (n < 8) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += a[(size_t)i * stride];
        return sum;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; j++) r[j] = a[(size_t)j * stride];
        int i = 8;
        for (; i < n - n % 8; i += 8)
            for (int j = 0; j < 8; j++) r[j] += a[(size_t)(i + j) * stride];
        double sum = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; i++) sum += a[(size_t)i * stride];
        return sum;
    }
    int half = n / 2;
    half -= half % 8;
    return numpy_sum(a, half, stride) + numpy_sum(a + (size_t)half * stride, n - half, stride);
}

/* coupled _subduction_frequency: the average character of the little irrep
 * over the operations times the kernel cosets inside the subgroup lattice;
 * mode: the point action of the records (mag_point_action). */
int subduction_frequency(int gid, int param_count, const int64_t *kparam, int mode, const int64_t *kernel_basis,
                         const int64_t *subgroup_basis, int op_count, const int64_t *records, int *frequency)
{
    *frequency = 0;
    if (!op_count) return OP_OK;
    if (mode && mag_tables()) return 100 + OP_TABLES;
    irrep_context c;
    int status = irrep_context_for(gid, param_count, kparam, &c);
    if (status) return 100 + context_status(status);
    int dim = c.dim;
    i128 b[9], adj[9];
    for (int i = 0; i < 9; i++) b[i] = kernel_basis[i];
    i128 kernel_det = op_abs(op_det(b));
    if (kernel_det > INT32_MAX) return 100 + OP_RANGE;
    int64_t *fractions = malloc(sizeof(int64_t) * 4 * (size_t)(kernel_det + 1));
    int *shifted = malloc(sizeof(int) * 5 * (size_t)(kernel_det + 1));
    if (!fractions || !shifted) {
        free(fractions);
        free(shifted);
        return 100 + OP_RANGE;
    }
    int fraction_count = 0;
    status = op_new_fractionals(kernel_basis, fractions, (int)kernel_det + 1, &fraction_count);
    if (status) status += 100;
    for (int i = 0; i < 9; i++) b[i] = subgroup_basis[i];
    i128 det = op_det(b);
    op_adjugate(b, adj);
    int kept = 0;
    for (int f = 0; f < fraction_count && det && !status; f++) {
        const int64_t *x = fractions + 4 * f;
        i128 scale = (i128)(x[3] ? x[3] : 1) * det;
        int inside = 1;
        for (int col = 0; col < 3 && inside; col++)
            inside = ((i128)x[0] * adj[col] + (i128)x[1] * adj[3 + col] + (i128)x[2] * adj[6 + col]) % scale == 0;
        if (inside) memmove(fractions + 4 * kept++, x, sizeof(int64_t) * 4);
    }
    double sum = 0.0, *stack = malloc(sizeof(double) * ((size_t)kept * dim * dim + 1));
    if (!stack && !status) status = 100 + OP_RANGE;
    for (int i = 0; i < op_count && kept && !status; i++) {
        int sign;
        int matrix_op = mag_point_action(mode, (int)records[5 * i + 4], &sign);
        if (mode && mag_take_error()) status = 100 + OP_TABLES;
        for (int f = 0; f < kept && !status; f++) {
            int64_t record[5];
            status = op_vadd(fractions + 4 * f, records + 5 * i, record);
            if (status) status += 100;
            for (int j = 0; j < 4; j++) shifted[5 * f + j] = (int)record[j];
            shifted[5 * f + 4] = matrix_op;
        }
        if (!status)
            status = irrep_matrix_stack(gid, param_count, kparam, kept, shifted, kept * dim * dim, stack);
        for (int f = 0; f < kept && !status; f++) {
            double trace = numpy_sum(stack + (size_t)f * dim * dim, dim, dim + 1);
            sum += sign < 0 ? -trace : trace;
        }
    }
    free(fractions);
    free(shifted);
    free(stack);
    if (!status && det && kept) *frequency = (int)rint(sum / ((double)op_count * kept));
    return status;
}
