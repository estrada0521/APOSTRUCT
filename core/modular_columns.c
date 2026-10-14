/* The column structure of a sparse integer matrix over GF(p): the pivot
 * columns of its reduced row echelon form (the first maximal independent
 * columns in order) and, for every other column, the support of its unique
 * expression through the pivot columns (the nullspace basis of the rref). */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum {
    MC_MEMORY = 1,
    MC_RANGE = 2,
};

static uint64_t mc_pow(uint64_t base, uint64_t exp, uint64_t p)
{
    uint64_t result = 1;
    base %= p;
    while (exp) {
        if (exp & 1) result = result * base % p;
        base = base * base % p;
        exp >>= 1;
    }
    return result;
}

typedef struct {
    int pivot_row;
    int count;
    int *rows;
    uint32_t *values;
} basis_vector;

/* columns: col_count sparse columns (offsets, rows, values).  pivots_out:
 * the pivot columns ascending, *pivot_count of them.  support_offsets /
 * supports: per non-pivot column (ascending), the ascending column indices
 * of its relation (itself included).  Returns 0 or a status. */
int modular_columns(int row_count, int col_count, const int64_t *offsets, const int *rows, const int64_t *values,
                    int64_t prime, int *pivots_out, int *pivot_count, int64_t support_capacity,
                    int64_t *support_offsets, int *supports)
{
    uint64_t p = (uint64_t)prime;
    if (p < 3 || p > (1u << 31)) return MC_RANGE;
    uint32_t *dense = calloc((size_t)row_count + 1, sizeof(uint32_t));
    unsigned char *touched = calloc((size_t)row_count + 1, 1);
    int *touched_rows = malloc(sizeof(int) * ((size_t)row_count + 1));
    basis_vector *basis = calloc((size_t)col_count + 1, sizeof(basis_vector));
    uint32_t *transform = calloc((size_t)col_count * (size_t)col_count + 1, sizeof(uint32_t));
    uint32_t *coefficients = calloc((size_t)col_count + 1, sizeof(uint32_t));
    int *row_pivot = malloc(sizeof(int) * ((size_t)row_count + 1));
    int status = !dense || !touched || !touched_rows || !basis || !transform || !coefficients || !row_pivot ? MC_MEMORY : 0;
    int basis_count = 0, free_count = 0;
    int64_t support_used = 0;
    for (int r = 0; r < row_count && !status; r++) row_pivot[r] = -1;
    if (!status) support_offsets[0] = 0;
    for (int col = 0; col < col_count && !status; col++) {
        int touched_count = 0;
        for (int64_t k = offsets[col]; k < offsets[col + 1]; k++) {
            int r = rows[k];
            if (r < 0 || r >= row_count) {
                status = MC_RANGE;
                break;
            }
            int64_t v = values[k] % (int64_t)p;
            if (v < 0) v += (int64_t)p;
            dense[r] = (uint32_t)(((uint64_t)dense[r] + (uint64_t)v) % p);
            if (!touched[r]) {
                touched[r] = 1;
                touched_rows[touched_count++] = r;
            }
        }
        if (status) break;
        memset(coefficients, 0, sizeof(uint32_t) * (size_t)(basis_count + 1));
        /* Reduce against the basis in insertion order; each basis vector is
         * zero on the pivot rows of the vectors inserted before it. */
        for (int b = 0; b < basis_count; b++) {
            uint32_t a = dense[basis[b].pivot_row];
            if (!a) continue;
            for (int k = 0; k < basis[b].count; k++) {
                int r = basis[b].rows[k];
                uint64_t product = (uint64_t)a * basis[b].values[k] % p;
                dense[r] = (uint32_t)(((uint64_t)dense[r] + p - product) % p);
                if (!touched[r]) {
                    touched[r] = 1;
                    touched_rows[touched_count++] = r;
                }
            }
            const uint32_t *t = transform + (size_t)b * col_count;
            for (int m = 0; m <= b; m++)
                if (t[m]) coefficients[m] = (uint32_t)(((uint64_t)coefficients[m] + (uint64_t)a * t[m]) % p);
        }
        int lead = -1, nonzero = 0;
        for (int i = 0; i < touched_count; i++) {
            int r = touched_rows[i];
            if (dense[r]) {
                nonzero++;
                if (lead < 0 || r < lead) lead = r;
            }
        }
        if (lead >= 0) {
            uint64_t inverse = mc_pow(dense[lead], p - 2, p);
            basis_vector *v = basis + basis_count;
            v->pivot_row = lead;
            v->count = 0;
            v->rows = malloc(sizeof(int) * (size_t)nonzero);
            v->values = malloc(sizeof(uint32_t) * (size_t)nonzero);
            if (!v->rows || !v->values) {
                status = MC_MEMORY;
                break;
            }
            for (int i = 0; i < touched_count; i++) {
                int r = touched_rows[i];
                if (!dense[r]) continue;
                v->rows[v->count] = r;
                v->values[v->count++] = (uint32_t)((uint64_t)dense[r] * inverse % p);
            }
            uint32_t *t = transform + (size_t)basis_count * col_count;
            for (int m = 0; m < basis_count; m++)
                if (coefficients[m]) t[m] = (uint32_t)((p - coefficients[m]) * inverse % p);
            t[basis_count] = (uint32_t)inverse;
            row_pivot[lead] = basis_count;
            pivots_out[basis_count++] = col;
        } else {
            int count = 1;
            for (int m = 0; m < basis_count; m++) count += coefficients[m] != 0;
            if (support_used + count > support_capacity) {
                status = MC_RANGE;
                break;
            }
            int *out = supports + support_used;
            int n = 0;
            for (int m = 0; m < basis_count; m++)
                if (coefficients[m]) out[n++] = pivots_out[m];
            out[n++] = col;
            support_used += n;
            support_offsets[++free_count] = support_used;
        }
        for (int i = 0; i < touched_count; i++) {
            dense[touched_rows[i]] = 0;
            touched[touched_rows[i]] = 0;
        }
    }
    *pivot_count = basis_count;
    for (int b = 0; b < basis_count; b++) {
        free(basis[b].rows);
        free(basis[b].values);
    }
    free(dense);
    free(touched);
    free(touched_rows);
    free(basis);
    free(transform);
    free(coefficients);
    free(row_pivot);
    return status;
}

#include <math.h>

static int is_prime32(uint64_t n)
{
    if (n < 2) return 0;
    static const uint64_t small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (size_t i = 0; i < sizeof small / sizeof small[0]; i++)
        if (n % small[i] == 0) return n == small[i];
    uint64_t d = n - 1;
    int r = 0;
    while (!(d & 1)) {
        d >>= 1;
        r++;
    }
    static const uint64_t bases[] = {2, 7, 61};
    for (size_t i = 0; i < 3; i++) {
        uint64_t x = mc_pow(bases[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (int j = 1; j < r && composite; j++) {
            x = x * x % n;
            composite = x != n - 1;
        }
        if (composite) return 0;
    }
    return 1;
}

static uint64_t prime_below(uint64_t n)
{
    do n--;
    while (!is_prime32(n));
    return n;
}

/* Rank over GF(p) of the sparse integer columns cols[0..k); -1 for memory. */
static int rank_mod(const int64_t *offsets, const int *rows, const int64_t *values, const int *cols, int k, uint64_t p,
                    int *row_of)
{
    int m = 0;
    for (int j = 0; j < k; j++)
        for (int64_t t = offsets[cols[j]]; t < offsets[cols[j] + 1]; t++)
            if (row_of[rows[t]] < 0) row_of[rows[t]] = m++;
    uint64_t *a = calloc((size_t)m * (size_t)k + 1, sizeof(uint64_t));
    if (!a) return -1;
    for (int j = 0; j < k; j++)
        for (int64_t t = offsets[cols[j]]; t < offsets[cols[j] + 1]; t++) {
            int64_t v = values[t] % (int64_t)p;
            if (v < 0) v += (int64_t)p;
            uint64_t *cell = a + (size_t)row_of[rows[t]] * k + j;
            *cell = (*cell + (uint64_t)v) % p;
        }
    for (int j = 0; j < k; j++)
        for (int64_t t = offsets[cols[j]]; t < offsets[cols[j] + 1]; t++) row_of[rows[t]] = -1;
    int rank = 0;
    for (int j = 0; j < k && rank < m; j++) {
        int pivot = -1;
        for (int r = rank; r < m && pivot < 0; r++)
            if (a[(size_t)r * k + j]) pivot = r;
        if (pivot < 0) continue;
        for (int c = 0; c < k; c++) {
            uint64_t swap = a[(size_t)pivot * k + c];
            a[(size_t)pivot * k + c] = a[(size_t)rank * k + c];
            a[(size_t)rank * k + c] = swap;
        }
        uint64_t inverse = mc_pow(a[(size_t)rank * k + j], p - 2, p);
        for (int r = rank + 1; r < m; r++) {
            uint64_t f = a[(size_t)r * k + j] * inverse % p;
            if (!f) continue;
            for (int c = j; c < k; c++) a[(size_t)r * k + c] = (a[(size_t)r * k + c] + p - f * a[(size_t)rank * k + c] % p) % p;
        }
        rank++;
    }
    free(a);
    return rank;
}

/* Whether the integer columns cols[0..k) are linearly dependent over Q:
 * every k x k minor is bounded by the product of the column norms
 * (Hadamard), so ranks below k modulo primes whose product exceeds twice
 * that bound prove each minor zero, while one prime of rank k proves
 * independence.  1 dependent, 0 independent, -1 memory. */
static int dependent_over_q(const int64_t *offsets, const int *rows, const int64_t *values, const int *cols, int k,
                            int *row_of)
{
    double bound_bits = 2.0;
    for (int j = 0; j < k; j++) {
        double norm = 0.0;
        for (int64_t t = offsets[cols[j]]; t < offsets[cols[j] + 1]; t++) norm += (double)values[t] * (double)values[t];
        bound_bits += 0.5 * log2(norm > 1.0 ? norm : 1.0);
    }
    double covered = 0.0;
    uint64_t p = (uint64_t)1 << 31;
    while (covered <= bound_bits) {
        p = prime_below(p);
        int rank = rank_mod(offsets, rows, values, cols, k, p, row_of);
        if (rank < 0) return -1;
        if (rank == k) return 0;
        covered += log2((double)p);
    }
    return 1;
}

/* The exact (over Q) rref pivot columns of a sparse integer matrix.  The
 * GF(p) pivots are exact once every GF(p) relation is proved dependent over
 * Q; a false relation retries with the next prime.  Returns 0, MC_MEMORY,
 * or MC_RANGE after 64 primes. */
int independent_columns(int row_count, int col_count, const int64_t *offsets, const int *rows, const int64_t *values,
                        int *pivots_out, int *pivot_count)
{
    int64_t capacity = ((int64_t)col_count + 1) * ((int64_t)col_count + 1);
    int64_t *support_offsets = malloc(sizeof(int64_t) * ((size_t)col_count + 1));
    int *supports = malloc(sizeof(int) * (size_t)capacity);
    int *row_of = malloc(sizeof(int) * ((size_t)row_count + 1));
    if (!support_offsets || !supports || !row_of) {
        free(support_offsets);
        free(supports);
        free(row_of);
        return MC_MEMORY;
    }
    for (int r = 0; r < row_count; r++) row_of[r] = -1;
    int status = MC_RANGE;
    uint64_t p = 1u << 31;
    for (int attempt = 0; attempt < 64; attempt++) {
        p = prime_below(p);
        int count = 0;
        int s = modular_columns(row_count, col_count, offsets, rows, values, (int64_t)p, pivots_out, &count, capacity, support_offsets,
                                supports);
        if (s) {
            status = s;
            break;
        }
        int proved = 1;
        for (int f = 0; f < col_count - count && proved > 0; f++)
            proved = dependent_over_q(offsets, rows, values, supports + support_offsets[f],
                                      (int)(support_offsets[f + 1] - support_offsets[f]), row_of);
        if (proved < 0) {
            status = MC_MEMORY;
            break;
        }
        if (proved) {
            *pivot_count = count;
            status = 0;
            break;
        }
    }
    free(support_offsets);
    free(supports);
    free(row_of);
    return status;
}
