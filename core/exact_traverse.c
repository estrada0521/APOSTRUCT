/* The isotropy-subspace traversal of isotropy_traverse.c over Q(zeta_m).
 *
 * The same breadth-first walk, equation keys, orbit representatives,
 * projective classes and child flags, with every test exact: a value is
 * nonzero when it is not 0 in the field, signs come from a certified
 * evaluation, and keys are the reduced equations themselves. Elements are
 * stored flat as phi coefficients and a denominator (xs = phi + 1 words). */

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "cyclo.h"

enum {
    XT_OK = 0,
    XT_NO_MEMORY = 1,
    XT_EQS_NO_PIVOT = 2,
    XT_EQS_REDUCED_NO_PIVOT = 3,
    XT_EQS_DUPLICATE_PIVOT = 4,
    XT_EQS_INVALID_COLUMN_STATE = 5,
    XT_EQS_INCONSISTENT_COUNT = 6,
    XT_E2O_NO_PIVOT = 7,
    XT_E2O_ZERO_PIVOT = 8,
    XT_E2O_INCONSISTENT = 9,
    XT_FIXED_FREE_CHANGED = 10,
    XT_FRACTION_NOT_FOUND = 11,
    XT_ARITHMETIC = 12,
};

#define TRY(expr)                        \
    do {                                 \
        int status_ = (expr);            \
        if (status_) return status_;     \
    } while (0)

static const cy_field *F;
static int XS; /* words per element */

typedef i128 xe; /* an element is XS consecutive words */

static xe *at(xe *base, size_t index) { return base + index * (size_t)XS; }
static const xe *cat(const xe *base, size_t index) { return base + index * (size_t)XS; }

static void ld(const xe *e, cy *t)
{
    memcpy(t->c, e, sizeof(i128) * (size_t)F->phi);
    t->den = e[F->phi];
}

static void st(const cy *t, xe *e)
{
    memcpy(e, t->c, sizeof(i128) * (size_t)F->phi);
    e[F->phi] = t->den;
}

static int ez(const xe *e)
{
    for (int i = 0; i < F->phi; i++)
        if (e[i]) return 0;
    return 1;
}

static void eint(xe *e, int64_t n)
{
    memset(e, 0, sizeof(i128) * (size_t)XS);
    e[0] = n;
    e[F->phi] = 1;
}

static void ecopy(xe *dst, const xe *src)
{
    if (dst != src) memcpy(dst, src, sizeof(i128) * (size_t)XS);
}

static int eeq(const xe *a, const xe *b) { return !memcmp(a, b, sizeof(i128) * (size_t)XS); }

static int erational(const xe *a)
{
    for (int i = 1; i < F->phi; i++)
        if (a[i]) return 0;
    return 1;
}

/* b scaled by the rational n / d */
static int escale(const xe *b, i128 n, i128 d, xe *out)
{
    cy x, q, r;
    ld(b, &x);
    cy_rat(F, n, d, &q);
    cy_mul(F, &x, &q, &r);
    st(&r, out);
    return cy_overflow ? XT_ARITHMETIC : 0;
}

static int eadd(const xe *a, const xe *b, xe *out)
{
    if (ez(b)) {
        ecopy(out, a);
        return 0;
    }
    if (ez(a)) {
        ecopy(out, b);
        return 0;
    }
    cy x, y, r;
    ld(a, &x);
    ld(b, &y);
    cy_add(F, &x, &y, &r);
    st(&r, out);
    return cy_overflow ? XT_ARITHMETIC : 0;
}

static int eneg(const xe *a, xe *out);

static int esub(const xe *a, const xe *b, xe *out)
{
    if (ez(b)) {
        ecopy(out, a);
        return 0;
    }
    cy x, y, r;
    ld(a, &x);
    ld(b, &y);
    cy_sub(F, &x, &y, &r);
    st(&r, out);
    return cy_overflow ? XT_ARITHMETIC : 0;
}

static int emul(const xe *a, const xe *b, xe *out)
{
    if (ez(a) || ez(b)) {
        eint(out, 0);
        return 0;
    }
    if (erational(a)) return escale(b, a[0], a[F->phi], out);
    if (erational(b)) return escale(a, b[0], b[F->phi], out);
    cy x, y, r;
    ld(a, &x);
    ld(b, &y);
    cy_mul(F, &x, &y, &r);
    st(&r, out);
    return cy_overflow ? XT_ARITHMETIC : 0;
}

static uint64_t mix64(uint64_t hash);
static uint64_t words_hash(uint64_t hash, const xe *values, size_t words);

/* inverses seen in this traversal, keyed by the element's words */
static struct {
    xe *keys, *values;
    int *slots;
    int count, capacity, slot_count;
} INV;

static void inverse_cache_reset(void)
{
    free(INV.keys);
    free(INV.values);
    free(INV.slots);
    memset(&INV, 0, sizeof INV);
}

static int inverse_of(const xe *b, cy *inv)
{
    uint64_t h = mix64(words_hash(1469598103934665603ULL, b, (size_t)XS));
    if (INV.slot_count) {
        for (int slot = (int)(h & (uint64_t)(INV.slot_count - 1)); INV.slots[slot] >= 0; slot = (slot + 1) & (INV.slot_count - 1)) {
            int id = INV.slots[slot];
            if (eeq(INV.keys + (size_t)id * XS, b)) {
                ld(INV.values + (size_t)id * XS, inv);
                return 0;
            }
        }
    }
    cy y;
    ld(b, &y);
    if (exact_inverse(F, &y, inv)) return XT_ARITHMETIC;
    if ((INV.count + 1) * 2 > INV.slot_count) {
        int slot_count = INV.slot_count ? INV.slot_count * 2 : 256;
        int *slots = malloc(sizeof(int) * (size_t)slot_count);
        if (!slots) return 0;
        for (int i = 0; i < slot_count; i++) slots[i] = -1;
        for (int id = 0; id < INV.count; id++) {
            uint64_t k = mix64(words_hash(1469598103934665603ULL, INV.keys + (size_t)id * XS, (size_t)XS));
            int slot = (int)(k & (uint64_t)(slot_count - 1));
            while (slots[slot] >= 0) slot = (slot + 1) & (slot_count - 1);
            slots[slot] = id;
        }
        free(INV.slots);
        INV.slots = slots;
        INV.slot_count = slot_count;
    }
    if (INV.count == INV.capacity) {
        int capacity = INV.capacity ? INV.capacity * 2 : 128;
        xe *keys = realloc(INV.keys, sizeof(xe) * (size_t)capacity * XS), *values = keys ? realloc(INV.values, sizeof(xe) * (size_t)capacity * XS) : 0;
        if (keys) INV.keys = keys;
        if (values) INV.values = values;
        if (!keys || !values) return 0;
        INV.capacity = capacity;
    }
    ecopy(INV.keys + (size_t)INV.count * XS, b);
    st(inv, INV.values + (size_t)INV.count * XS);
    int slot = (int)(h & (uint64_t)(INV.slot_count - 1));
    while (INV.slots[slot] >= 0) slot = (slot + 1) & (INV.slot_count - 1);
    INV.slots[slot] = INV.count++;
    return 0;
}

static int ediv(const xe *a, const xe *b, xe *out)
{
    if (ez(a)) {
        eint(out, 0);
        return 0;
    }
    if (erational(b)) {
        if (!b[0]) return XT_ARITHMETIC;
        return escale(a, b[F->phi], b[0], out);
    }
    cy x, inv, r;
    ld(a, &x);
    if (inverse_of(b, &inv)) return XT_ARITHMETIC;
    cy_mul(F, &x, &inv, &r);
    st(&r, out);
    return cy_overflow ? XT_ARITHMETIC : 0;
}

static int eneg(const xe *a, xe *out)
{
    ecopy(out, a);
    for (int i = 0; i < F->phi; i++) out[i] = -out[i];
    return 0;
}

static int esign(const xe *a, int *sign)
{
    cy x;
    ld(a, &x);
    return cy_sign(F, &x, sign) ? XT_ARITHMETIC : 0;
}

/* |a| == 1 */
static int eunit(const xe *a)
{
    if (a[F->phi] != 1 || (a[0] != 1 && a[0] != -1)) return 0;
    for (int i = 1; i < F->phi; i++)
        if (a[i]) return 0;
    return 1;
}

typedef struct {
    int *data;
    size_t used, capacity;
} IntList;

static int list_push(IntList *list, int value)
{
    if (list->used == list->capacity) {
        size_t capacity = list->capacity ? list->capacity * 2 : 64;
        int *data = realloc(list->data, capacity * sizeof(int));
        if (!data) return XT_NO_MEMORY;
        list->data = data;
        list->capacity = capacity;
    }
    list->data[list->used++] = value;
    return XT_OK;
}

/* Equation keys: the reduced equations, column by column. */
typedef struct {
    xe *values;
    size_t used_words, capacity_words;
    size_t *offset;
    int *eq_count, *length;
    uint64_t *hash;
    int count, capacity;
    int *slots;
    size_t slot_count;
} KeyTable;

static uint64_t mix64(uint64_t hash)
{
    hash ^= hash >> 30;
    hash *= 0xbf58476d1ce4e5b9ULL;
    hash ^= hash >> 27;
    hash *= 0x94d049bb133111ebULL;
    return hash ^ (hash >> 31);
}

static uint64_t words_hash(uint64_t hash, const xe *values, size_t words)
{
    const uint64_t *w = (const uint64_t *)values;
    for (size_t i = 0; i < 2 * words; i++) hash = (hash ^ w[i]) * 1099511628211ULL;
    return hash;
}

static int key_grow_slots(KeyTable *table)
{
    size_t slot_count = table->slot_count ? table->slot_count * 2 : 1024;
    int *slots = malloc(slot_count * sizeof(int));
    if (!slots) return XT_NO_MEMORY;
    for (size_t index = 0; index < slot_count; index++) slots[index] = -1;
    for (int id = 0; id < table->count; id++) {
        size_t slot = table->hash[id] & (slot_count - 1);
        while (slots[slot] >= 0) slot = (slot + 1) & (slot_count - 1);
        slots[slot] = id;
    }
    free(table->slots);
    table->slots = slots;
    table->slot_count = slot_count;
    return XT_OK;
}

static int key_intern(KeyTable *table, int eq_count, const xe *values, int length, int *id_out)
{
    if ((size_t)(table->count + 1) * 2 > table->slot_count) TRY(key_grow_slots(table));
    size_t words = (size_t)length * XS;
    uint64_t hash = mix64(words_hash(1469598103934665603ULL ^ (uint64_t)eq_count, values, words));
    size_t slot = hash & (table->slot_count - 1);
    while (table->slots[slot] >= 0) {
        int id = table->slots[slot];
        if (table->hash[id] == hash && table->eq_count[id] == eq_count && table->length[id] == length &&
            !memcmp(table->values + table->offset[id], values, words * sizeof(xe))) {
            *id_out = id;
            return XT_OK;
        }
        slot = (slot + 1) & (table->slot_count - 1);
    }
    if (table->count == table->capacity) {
        int capacity = table->capacity ? table->capacity * 2 : 1024;
        size_t *offset = realloc(table->offset, (size_t)capacity * sizeof(size_t));
        int *counts = realloc(table->eq_count, (size_t)capacity * sizeof(int));
        int *lengths = realloc(table->length, (size_t)capacity * sizeof(int));
        uint64_t *hashes = realloc(table->hash, (size_t)capacity * sizeof(uint64_t));
        if (offset) table->offset = offset;
        if (counts) table->eq_count = counts;
        if (lengths) table->length = lengths;
        if (hashes) table->hash = hashes;
        if (!offset || !counts || !lengths || !hashes) return XT_NO_MEMORY;
        table->capacity = capacity;
    }
    if (table->used_words + words > table->capacity_words) {
        size_t capacity = table->capacity_words ? table->capacity_words : 4096;
        while (capacity < table->used_words + words) capacity *= 2;
        xe *grown = realloc(table->values, capacity * sizeof(xe));
        if (!grown) return XT_NO_MEMORY;
        table->values = grown;
        table->capacity_words = capacity;
    }
    int id = table->count++;
    table->offset[id] = table->used_words;
    memcpy(table->values + table->used_words, values, words * sizeof(xe));
    table->used_words += words;
    table->eq_count[id] = eq_count;
    table->length[id] = length;
    table->hash[id] = hash;
    table->slots[slot] = id;
    *id_out = id;
    return XT_OK;
}

/* xrowop2 on a rows x cols matrix, in place. */
static int xrowop2(xe *matrix, int rows, int cols, xe *swap, xe *tmp)
{
    int pivot_col = 0;
    for (int pivot_row = 0; pivot_row < rows; pivot_row++) {
        int reduced = 0;
        while (pivot_col < cols) {
            int swap_row = -1;
            for (int row = pivot_row; row < rows; row++)
                if (!ez(at(matrix, (size_t)row * cols + pivot_col))) {
                    swap_row = row;
                    break;
                }
            if (swap_row < 0) {
                pivot_col += 1;
                continue;
            }
            if (swap_row != pivot_row) {
                size_t words = (size_t)cols * XS;
                memcpy(swap, at(matrix, (size_t)pivot_row * cols), words * sizeof(xe));
                memcpy(at(matrix, (size_t)pivot_row * cols), at(matrix, (size_t)swap_row * cols), words * sizeof(xe));
                memcpy(at(matrix, (size_t)swap_row * cols), swap, words * sizeof(xe));
            }
            xe *pivot_values = at(matrix, (size_t)pivot_row * cols);
            ecopy(tmp, at(pivot_values, (size_t)pivot_col));
            for (int col = pivot_col; col < cols; col++) TRY(ediv(at(pivot_values, (size_t)col), tmp, at(pivot_values, (size_t)col)));
            for (int row = 0; row < rows; row++) {
                if (row == pivot_row) continue;
                xe *row_values = at(matrix, (size_t)row * cols);
                xe *factor = at(row_values, (size_t)pivot_col);
                if (ez(factor)) continue;
                xe *product = tmp + XS;
                for (int col = 0; col < cols; col++) {
                    if (col == pivot_col || ez(at(pivot_values, (size_t)col))) continue;
                    TRY(emul(factor, at(pivot_values, (size_t)col), product));
                    TRY(esub(at(row_values, (size_t)col), product, at(row_values, (size_t)col)));
                }
                eint(factor, 0);
            }
            pivot_col += 1;
            reduced = 1;
            break;
        }
        if (!reduced) break;
    }
    return XT_OK;
}

typedef struct {
    int dim, space_count, fraction_count, pair_count, matrix_count;
    const xe *matrices; /* matrix_count x dim x dim, [source_row][col] */
    const int *pair_matrix;
    const int *space_product;
    const int *space_fraction;
    const int *fraction_residue;
    int residue_modulus;
    int *residue_slots;
    size_t residue_slot_count;
    int *projective;
    int *orbit_pairs;
    int orbit_pair_count;
    KeyTable keys;
    unsigned char *seen;
    size_t seen_capacity;
    xe *work, *equations, *equations2, *orderparam, *orderparam2, *transformed, *key_values, *swap, *terms, *tmp;
    int *column_state, *assign_eq, *assign_col, *pivots;
    xe *assign_value;
    unsigned char *is_pivot;
    int *matrix_key, *orbit_matrix_key, *check_cache, *transformed_key_cache;
    signed char *flags;
    int *minus_list;
    int *nz_start, *nz_row; /* per matrix, dim + 1 column starts into nz_row (source rows) */
} Context;

static const xe *matrix_of_pair(const Context *ctx, int pair)
{
    return cat(ctx->matrices, (size_t)ctx->pair_matrix[pair] * ctx->dim * ctx->dim);
}

static int matrix_index(const Context *ctx, const xe *matrix)
{
    return (int)((matrix - ctx->matrices) / ((ptrdiff_t)ctx->dim * ctx->dim * XS));
}

static int orderparam_to_eqs(Context *ctx, int free_rows, const xe *orderparam, xe *equations, int *eq_count_out)
{
    int cols = ctx->dim;
    xe *reduced = ctx->work, *t = ctx->tmp;
    memcpy(reduced, orderparam, (size_t)free_rows * cols * XS * sizeof(xe));
    for (int row = 0; row < free_rows; row++) {
        int pivot = -1;
        for (int col = 0; col < cols; col++)
            if (!ez(at(reduced, (size_t)row * cols + col))) {
                pivot = col;
                break;
            }
        if (pivot < 0) return XT_EQS_NO_PIVOT;
        ecopy(t, at(reduced, (size_t)row * cols + pivot));
        for (int col = pivot; col < cols; col++) TRY(ediv(at(reduced, (size_t)row * cols + col), t, at(reduced, (size_t)row * cols + col)));
        for (int other = 0; other < free_rows; other++) {
            if (other == row) continue;
            ecopy(t, at(reduced, (size_t)other * cols + pivot));
            if (ez(t)) continue;
            for (int col = 0; col < cols; col++) {
                if (ez(at(reduced, (size_t)row * cols + col))) continue;
                TRY(emul(t, at(reduced, (size_t)row * cols + col), t + XS));
                TRY(esub(at(reduced, (size_t)other * cols + col), t + XS, at(reduced, (size_t)other * cols + col)));
            }
        }
    }
    int *column_state = ctx->column_state;
    memset(column_state, 0, (size_t)cols * sizeof(int));
    int eq_count = 0, assigned = 0;
    for (int row = 0; row < free_rows; row++) {
        int pivot = -1;
        for (int col = 0; col < cols; col++)
            if (!ez(at(reduced, (size_t)row * cols + col))) {
                pivot = col;
                break;
            }
        if (pivot < 0) return XT_EQS_REDUCED_NO_PIVOT;
        if (column_state[pivot] != 0) return XT_EQS_DUPLICATE_PIVOT;
        column_state[pivot] = -1;
        for (int col = pivot + 1; col < cols; col++) {
            const xe *value = at(reduced, (size_t)row * cols + col);
            if (ez(value)) continue;
            int existing = column_state[col];
            if (existing < 1) {
                if (existing != 0) return XT_EQS_INVALID_COLUMN_STATE;
                eq_count += 1;
                column_state[col] = eq_count;
                ctx->assign_eq[assigned] = eq_count - 1;
                ctx->assign_col[assigned] = col;
                eint(at(ctx->assign_value, (size_t)assigned++), 1);
            }
            ctx->assign_eq[assigned] = column_state[col] - 1;
            ctx->assign_col[assigned] = pivot;
            eneg(value, at(ctx->assign_value, (size_t)assigned++));
        }
    }
    for (int col = 0; col < cols; col++)
        if (column_state[col] == 0) {
            eq_count += 1;
            ctx->assign_eq[assigned] = eq_count - 1;
            ctx->assign_col[assigned] = col;
            eint(at(ctx->assign_value, (size_t)assigned++), 1);
        }
    if (cols - free_rows != eq_count) return XT_EQS_INCONSISTENT_COUNT;
    for (size_t i = 0; i < (size_t)eq_count * cols; i++) eint(at(equations, i), 0);
    for (int index = 0; index < assigned; index++)
        ecopy(at(equations, (size_t)ctx->assign_eq[index] * cols + (cols - ctx->assign_col[index] - 1)), at(ctx->assign_value, (size_t)index));
    *eq_count_out = eq_count;
    return XT_OK;
}

static int eqs_to_orderparam(Context *ctx, int cols, int eq_count, const xe *equations, xe *orderparam, int *free_out)
{
    memset(ctx->is_pivot, 0, (size_t)cols);
    for (int row = 0; row < eq_count; row++) {
        int pivot = -1;
        for (int col = 0; col < cols; col++)
            if (!ez(cat(equations, (size_t)row * cols + col))) {
                pivot = col;
                break;
            }
        if (pivot < 0) return XT_E2O_NO_PIVOT;
        ctx->pivots[row] = pivot;
        ctx->is_pivot[pivot] = 1;
    }
    int free_count = 0;
    for (int source_col = cols - 1; source_col >= 0; source_col--) {
        if (ctx->is_pivot[source_col]) continue;
        xe *out = at(orderparam, (size_t)free_count * cols);
        free_count += 1;
        for (int col = 0; col < cols; col++) eint(at(out, (size_t)col), 0);
        eint(at(out, (size_t)(cols - source_col - 1)), 1);
        for (int eq_row = 0; eq_row < eq_count; eq_row++) {
            int pivot_col = ctx->pivots[eq_row];
            const xe *pivot_value = cat(equations, (size_t)eq_row * cols + pivot_col);
            if (ez(pivot_value)) return XT_E2O_ZERO_PIVOT;
            xe *target = at(out, (size_t)(cols - pivot_col - 1));
            TRY(ediv(cat(equations, (size_t)eq_row * cols + source_col), pivot_value, target));
            eneg(target, target);
        }
    }
    if (cols - eq_count != free_count) return XT_E2O_INCONSISTENT;
    *free_out = free_count;
    return XT_OK;
}

static int equations_are_reduced(int cols, int eq_count, const xe *equations)
{
    int previous = -1;
    for (int row = 0; row < eq_count; row++) {
        int pivot = -1;
        for (int col = 0; col < cols; col++) {
            const xe *value = cat(equations, (size_t)row * cols + col);
            if (!ez(value)) {
                pivot = col;
                if (!(value[0] == 1 && value[F->phi] == 1 && eunit(value))) return 0;
                break;
            }
        }
        if (pivot < 0) return 0;
        if (previous >= 0 && pivot <= previous) return 0;
        previous = pivot;
        for (int other = 0; other < eq_count; other++)
            if (other != row && !ez(cat(equations, (size_t)other * cols + pivot))) return 0;
    }
    return 1;
}

static uint64_t residue_hash(const int *residue)
{
    uint64_t hash = 1469598103934665603ULL;
    for (int axis = 0; axis < 3; axis++) hash = (hash ^ (uint64_t)(uint32_t)residue[axis]) * 1099511628211ULL;
    return hash;
}

static int fraction_of_residue(const Context *ctx, const int *residue, int *fraction_out)
{
    size_t slot = residue_hash(residue) & (ctx->residue_slot_count - 1);
    while (ctx->residue_slots[slot] >= 0) {
        int fraction = ctx->residue_slots[slot];
        const int *stored = ctx->fraction_residue + fraction * 3;
        if (stored[0] == residue[0] && stored[1] == residue[1] && stored[2] == residue[2]) {
            *fraction_out = fraction;
            return XT_OK;
        }
        slot = (slot + 1) & (ctx->residue_slot_count - 1);
    }
    return XT_FRACTION_NOT_FOUND;
}

static int pair_product(const Context *ctx, int left, int right, int *product)
{
    int fraction_count = ctx->fraction_count;
    int left_space = left / fraction_count, left_fraction = left % fraction_count;
    int right_space = right / fraction_count, right_fraction = right % fraction_count;
    const int *space = ctx->space_product + (left_space * ctx->space_count + right_space) * 2;
    int correction = space[1] - 1;
    int rotated = ctx->space_fraction[left_space * fraction_count + right_fraction] - 1;
    int residue[3];
    for (int axis = 0; axis < 3; axis++) {
        long long sum = (long long)ctx->fraction_residue[left_fraction * 3 + axis] + ctx->fraction_residue[correction * 3 + axis] +
                        ctx->fraction_residue[rotated * 3 + axis];
        residue[axis] = (int)(sum % ctx->residue_modulus);
    }
    int fraction;
    TRY(fraction_of_residue(ctx, residue, &fraction));
    *product = (space[0] - 1) * fraction_count + fraction;
    return XT_OK;
}

/* out[row][col] = sum over source rows k of matrix[k][col] orderparam[row][k] */
static int transform(Context *ctx, int free_rows, const xe *orderparam, const xe *matrix, xe *out)
{
    int dim = ctx->dim, mi = matrix_index(ctx, matrix);
    const int *start = ctx->nz_start + (size_t)mi * (dim + 1);
    for (int row = 0; row < free_rows; row++)
        for (int col = 0; col < dim; col++) {
            xe *total = at(out, (size_t)row * dim + col);
            eint(total, 0);
            for (int k = start[col]; k < start[col + 1]; k++) {
                int source = ctx->nz_row[k];
                const xe *x = cat(orderparam, (size_t)row * dim + source);
                if (ez(x)) continue;
                TRY(emul(cat(matrix, (size_t)source * dim + col), x, ctx->tmp));
                TRY(eadd(total, ctx->tmp, total));
            }
        }
    return XT_OK;
}

static int orderparam_check(Context *ctx, int free_rows, const xe *orderparam, const xe *matrix, int *preserves)
{
    TRY(transform(ctx, free_rows, orderparam, matrix, ctx->transformed));
    *preserves = !memcmp(ctx->transformed, orderparam, (size_t)free_rows * ctx->dim * XS * sizeof(xe));
    return XT_OK;
}

static int equation_data(Context *ctx, int free_rows, const xe *orderparam, xe *equations, int *eq_count_out, int *key_out)
{
    int dim = ctx->dim, eq_count;
    TRY(orderparam_to_eqs(ctx, free_rows, orderparam, equations, &eq_count));
    xe *reduced = ctx->equations2;
    memcpy(reduced, equations, (size_t)eq_count * dim * XS * sizeof(xe));
    if (!equations_are_reduced(dim, eq_count, reduced)) TRY(xrowop2(reduced, eq_count, dim, ctx->swap, ctx->tmp));
    int length = 0;
    for (int col = 0; col < dim; col++)
        for (int row = 0; row < eq_count; row++) ecopy(at(ctx->key_values, (size_t)length++), at(reduced, (size_t)row * dim + col));
    TRY(key_intern(&ctx->keys, eq_count, ctx->key_values, length, key_out));
    if ((size_t)ctx->keys.count > ctx->seen_capacity) {
        size_t capacity = ctx->seen_capacity ? ctx->seen_capacity * 2 : 4096;
        while (capacity < (size_t)ctx->keys.count) capacity *= 2;
        unsigned char *seen = realloc(ctx->seen, capacity);
        if (!seen) return XT_NO_MEMORY;
        memset(seen + ctx->seen_capacity, 0, capacity - ctx->seen_capacity);
        ctx->seen = seen;
        ctx->seen_capacity = capacity;
    }
    *eq_count_out = eq_count;
    return XT_OK;
}

static int nice_score(int dim, int free_rows, const xe *orderparam, int score[5])
{
    memset(score, 0, 5 * sizeof(int));
    for (int component = 0; component < dim; component++) {
        int seen_nonzero = 0;
        for (int row = 0; row < free_rows; row++) {
            const xe *value = cat(orderparam, (size_t)row * dim + component);
            if (ez(value)) continue;
            if (!eunit(value)) score[1] += 1;
            if (seen_nonzero)
                score[0] += 1;
            else
                score[2] += 1;
            int sign;
            TRY(esign(value, &sign));
            if (sign < 0) score[3] += 1;
            score[4] += component + 1;
            seen_nonzero = 1;
        }
    }
    return XT_OK;
}

static int nice_prefers_second(int dim, int free_rows, const xe *first, const xe *second, int *prefers)
{
    int first_nonzero = !ez(first), second_nonzero = !ez(second);
    *prefers = 0;
    if (first_nonzero && !second_nonzero) return XT_OK;
    if (!first_nonzero && second_nonzero) {
        *prefers = 1;
        return XT_OK;
    }
    int a[5], b[5];
    TRY(nice_score(dim, free_rows, first, a));
    TRY(nice_score(dim, free_rows, second, b));
    for (int index = 0; index < 5; index++) {
        if (a[index] < b[index]) return XT_OK;
        if (b[index] < a[index]) {
            *prefers = 1;
            return XT_OK;
        }
    }
    return XT_OK;
}

static int is_full(int dim, int free_rows, const xe *orderparam)
{
    if (free_rows != dim) return 0;
    for (int row = 0; row < dim; row++)
        for (int col = 0; col < dim; col++) {
            const xe *v = cat(orderparam, (size_t)row * dim + col);
            if (row == col ? !(eunit(v) && v[0] == 1) : !ez(v)) return 0;
        }
    return 1;
}

static int orbit(Context *ctx, int free_rows, const xe *orderparam, int use_projective, IntList *keys, xe *representative)
{
    int dim = ctx->dim, eq_count, key;
    size_t words = (size_t)free_rows * dim * XS;
    TRY(equation_data(ctx, free_rows, orderparam, ctx->equations, &eq_count, &key));
    keys->used = 0;
    TRY(list_push(keys, key));
    memcpy(representative, orderparam, words * sizeof(xe));
    for (int matrix = 0; matrix < ctx->matrix_count; matrix++) ctx->matrix_key[matrix] = ctx->orbit_matrix_key[matrix] = -1;
    if (is_full(dim, free_rows, orderparam)) return XT_OK;
    int pair_count = use_projective ? ctx->orbit_pair_count : ctx->pair_count;
    for (int index = 0; index < pair_count; index++) {
        int pair = use_projective ? ctx->orbit_pairs[index] : index;
        int matrix = ctx->pair_matrix[pair];
        int matrix_key = ctx->matrix_key[matrix];
        if (matrix_key < 0) {
            TRY(transform(ctx, free_rows, orderparam, matrix_of_pair(ctx, pair), ctx->orderparam2));
            TRY(equation_data(ctx, free_rows, ctx->orderparam2, ctx->equations, &eq_count, &matrix_key));
            int known = 0;
            for (size_t item = 0; item < keys->used; item++)
                if (keys->data[item] == matrix_key) {
                    known = 1;
                    break;
                }
            if (!known) {
                TRY(list_push(keys, matrix_key));
                memcpy(ctx->equations2, ctx->equations, (size_t)eq_count * dim * XS * sizeof(xe));
                TRY(xrowop2(ctx->equations2, eq_count, dim, ctx->swap, ctx->tmp));
                int next_free, prefers;
                TRY(eqs_to_orderparam(ctx, dim, eq_count, ctx->equations2, ctx->orderparam, &next_free));
                if (next_free == free_rows) {
                    TRY(nice_prefers_second(dim, free_rows, representative, ctx->orderparam, &prefers));
                    if (prefers) memcpy(representative, ctx->orderparam, words * sizeof(xe));
                }
            }
            ctx->matrix_key[matrix] = matrix_key;
        }
        ctx->orbit_matrix_key[matrix] = matrix_key;
    }
    if (use_projective)
        for (int matrix = 0; matrix < ctx->matrix_count; matrix++) {
            int representative_key = ctx->orbit_matrix_key[ctx->projective[matrix]];
            if (representative_key >= 0) ctx->orbit_matrix_key[matrix] = representative_key;
        }
    return XT_OK;
}

static int fixed_space_intersection(Context *ctx, int free_rows, const xe *orderparam, const xe *matrix, int *found, int *next_free_out,
                                    xe *next_orderparam, int *key_out)
{
    int dim = ctx->dim;
    *found = 0;
    TRY(transform(ctx, free_rows, orderparam, matrix, ctx->transformed));
    xe *coefficients = ctx->terms; /* dim x free_rows */
    int any = 0;
    for (int component = 0; component < dim; component++)
        for (int source_row = 0; source_row < free_rows; source_row++) {
            xe *target = at(coefficients, (size_t)component * free_rows + (free_rows - source_row - 1));
            TRY(esub(at(ctx->transformed, (size_t)source_row * dim + component), cat(orderparam, (size_t)source_row * dim + component), target));
            any = any || !ez(target);
        }
    if (!any) return XT_OK;
    TRY(xrowop2(coefficients, dim, free_rows, ctx->swap, ctx->tmp));
    int rank = dim;
    while (rank > 0) {
        int row_nonzero = 0;
        for (int col = 0; col < free_rows; col++)
            if (!ez(at(coefficients, (size_t)(rank - 1) * free_rows + col))) {
                row_nonzero = 1;
                break;
            }
        if (row_nonzero) break;
        rank -= 1;
    }
    if (rank <= 0 || rank >= free_rows) return XT_OK;
    int next_free;
    xe *basis = ctx->orderparam2; /* next_free x free_rows */
    TRY(eqs_to_orderparam(ctx, free_rows, rank, coefficients, basis, &next_free));
    xe *raw = ctx->orderparam;
    for (int row = 0; row < next_free; row++)
        for (int component = 0; component < dim; component++) {
            xe *total = at(raw, (size_t)row * dim + component);
            eint(total, 0);
            for (int source_row = 0; source_row < free_rows; source_row++) {
                const xe *b = at(basis, (size_t)row * free_rows + source_row), *x = cat(orderparam, (size_t)source_row * dim + component);
                if (ez(b) || ez(x)) continue;
                TRY(emul(b, x, ctx->tmp));
                TRY(eadd(total, ctx->tmp, total));
            }
        }
    int eq_count, key, canonical_free;
    TRY(equation_data(ctx, next_free, raw, ctx->equations, &eq_count, &key));
    TRY(eqs_to_orderparam(ctx, dim, eq_count, ctx->equations, next_orderparam, &canonical_free));
    if (canonical_free != next_free) return XT_FIXED_FREE_CHANGED;
    *found = 1;
    *next_free_out = next_free;
    *key_out = key;
    return XT_OK;
}

static void mark_minus(Context *ctx, int pair, int *minus_count)
{
    if (ctx->flags[pair] == -1) return;
    ctx->flags[pair] = -1;
    ctx->minus_list[(*minus_count)++] = pair;
}

static int close_child_flags(Context *ctx, int free_rows, const xe *orderparam, const signed char *parent_flags, int selected_pair)
{
    int pair_count = ctx->pair_count, minus_count = 0, product;
    memset(ctx->flags, 0, (size_t)pair_count);
    for (int matrix = 0; matrix < ctx->matrix_count; matrix++) ctx->check_cache[matrix] = ctx->transformed_key_cache[matrix] = -1;
    for (int marked = 0; marked < pair_count; marked++) {
        if (parent_flags[marked] != -1) continue;
        mark_minus(ctx, marked, &minus_count);
        TRY(pair_product(ctx, selected_pair, marked, &product));
        mark_minus(ctx, product, &minus_count);
        TRY(pair_product(ctx, marked, selected_pair, &product));
        mark_minus(ctx, product, &minus_count);
    }
    for (int pair = 0; pair < pair_count; pair++) {
        if (ctx->flags[pair] != 0) continue;
        int matrix = ctx->pair_matrix[pair];
        int preserves = ctx->check_cache[matrix];
        if (preserves < 0) {
            TRY(orderparam_check(ctx, free_rows, orderparam, matrix_of_pair(ctx, pair), &preserves));
            ctx->check_cache[matrix] = preserves;
        }
        if (preserves) {
            int snapshot = minus_count;
            for (int index = 0; index < snapshot; index++) {
                int marked = ctx->minus_list[index];
                TRY(pair_product(ctx, pair, marked, &product));
                mark_minus(ctx, product, &minus_count);
                TRY(pair_product(ctx, marked, pair, &product));
                mark_minus(ctx, product, &minus_count);
            }
            continue;
        }
        int key = ctx->orbit_matrix_key[matrix];
        if (key < 0) key = ctx->transformed_key_cache[matrix];
        if (key < 0) {
            int eq_count;
            TRY(transform(ctx, free_rows, orderparam, matrix_of_pair(ctx, pair), ctx->orderparam2));
            TRY(equation_data(ctx, free_rows, ctx->orderparam2, ctx->equations, &eq_count, &key));
            ctx->transformed_key_cache[matrix] = key;
        }
        if (!ctx->seen[key]) continue;
        for (int index = 0; index < minus_count; index++) {
            TRY(pair_product(ctx, pair, ctx->minus_list[index], &product));
            ctx->flags[product] = 1;
        }
    }
    return XT_OK;
}

/* The sign of the first nonzero entry of a matrix (1 when all are zero). */
static int first_sign(const xe *values, int cells, int *sign)
{
    *sign = 1;
    for (int index = 0; index < cells; index++)
        if (!ez(cat(values, (size_t)index))) return esign(cat(values, (size_t)index), sign);
    return XT_OK;
}

static int projective_pairs(Context *ctx)
{
    int dim = ctx->dim, cells = dim * dim;
    int *sign_of = malloc(sizeof(int) * ((size_t)ctx->matrix_count + 1));
    if (!sign_of) return XT_NO_MEMORY;
    for (int matrix = 0; matrix < ctx->matrix_count; matrix++) {
        ctx->projective[matrix] = -1;
        int status = first_sign(cat(ctx->matrices, (size_t)matrix * cells), cells, &sign_of[matrix]);
        if (status) {
            free(sign_of);
            return status;
        }
    }
    ctx->orbit_pair_count = 0;
    for (int pair = 0; pair < ctx->pair_count; pair++) {
        int matrix = ctx->pair_matrix[pair];
        if (ctx->projective[matrix] >= 0) continue;
        const xe *values = cat(ctx->matrices, (size_t)matrix * cells);
        int representative = -1;
        for (int earlier = 0; earlier < ctx->orbit_pair_count && representative < 0; earlier++) {
            int other_matrix = ctx->pair_matrix[ctx->orbit_pairs[earlier]];
            const xe *other = cat(ctx->matrices, (size_t)other_matrix * cells);
            int same_sign = sign_of[matrix] == sign_of[other_matrix], equal = 1;
            for (int cell = 0; cell < cells && equal; cell++) {
                const xe *x = cat(values, (size_t)cell), *y = cat(other, (size_t)cell);
                if (same_sign) {
                    equal = eeq(x, y);
                } else {
                    eneg(y, ctx->tmp);
                    equal = eeq(x, ctx->tmp);
                }
            }
            if (equal) representative = other_matrix;
        }
        if (representative < 0) {
            ctx->projective[matrix] = matrix;
            ctx->orbit_pairs[ctx->orbit_pair_count++] = pair;
        } else {
            ctx->projective[matrix] = representative;
        }
    }
    free(sign_of);
    return XT_OK;
}

typedef struct {
    int free_rows;
    xe *node_orderparam, *row_orderparam;
    signed char *flags;
} Node;

static void context_free(Context *ctx)
{
    free(ctx->residue_slots);
    free(ctx->projective);
    free(ctx->orbit_pairs);
    free(ctx->keys.values);
    free(ctx->keys.offset);
    free(ctx->keys.eq_count);
    free(ctx->keys.length);
    free(ctx->keys.hash);
    free(ctx->keys.slots);
    free(ctx->seen);
    free(ctx->work);
    free(ctx->equations);
    free(ctx->equations2);
    free(ctx->orderparam);
    free(ctx->orderparam2);
    free(ctx->transformed);
    free(ctx->key_values);
    free(ctx->swap);
    free(ctx->terms);
    free(ctx->tmp);
    free(ctx->column_state);
    free(ctx->assign_eq);
    free(ctx->assign_col);
    free(ctx->assign_value);
    free(ctx->pivots);
    free(ctx->is_pivot);
    free(ctx->matrix_key);
    free(ctx->orbit_matrix_key);
    free(ctx->check_cache);
    free(ctx->transformed_key_cache);
    free(ctx->flags);
    free(ctx->minus_list);
    free(ctx->nz_start);
    free(ctx->nz_row);
}

static xe *elements(size_t count) { return malloc((count + 2) * (size_t)XS * sizeof(xe)); }

static int context_init(Context *ctx)
{
    size_t dim = (size_t)ctx->dim, square = dim * dim + 1, assign = dim * (2 * dim + 1) + dim + 1;
    size_t matrices = (size_t)ctx->matrix_count + 1, pairs = (size_t)ctx->pair_count + 1;
    ctx->work = elements(square);
    ctx->equations = elements(square);
    ctx->equations2 = elements(square);
    ctx->orderparam = elements(square);
    ctx->orderparam2 = elements(square);
    ctx->transformed = elements(square);
    ctx->key_values = elements(square);
    ctx->swap = elements(square);
    ctx->terms = elements(square + dim);
    ctx->tmp = elements(4);
    ctx->column_state = malloc((dim + 1) * sizeof(int));
    ctx->assign_eq = malloc(assign * sizeof(int));
    ctx->assign_col = malloc(assign * sizeof(int));
    ctx->assign_value = elements(assign);
    ctx->pivots = malloc((dim + 1) * sizeof(int));
    ctx->is_pivot = malloc(dim + 1);
    ctx->matrix_key = malloc(matrices * sizeof(int));
    ctx->orbit_matrix_key = malloc(matrices * sizeof(int));
    ctx->check_cache = malloc(matrices * sizeof(int));
    ctx->transformed_key_cache = malloc(matrices * sizeof(int));
    ctx->projective = malloc(matrices * sizeof(int));
    ctx->flags = malloc(pairs);
    ctx->minus_list = malloc(pairs * sizeof(int));
    ctx->orbit_pairs = malloc(pairs * sizeof(int));
    ctx->nz_start = malloc(matrices * (dim + 1) * sizeof(int));
    ctx->nz_row = malloc(matrices * square * sizeof(int));
    size_t slot_count = 16;
    while (slot_count < (size_t)ctx->fraction_count * 2) slot_count *= 2;
    ctx->residue_slots = malloc(slot_count * sizeof(int));
    ctx->residue_slot_count = slot_count;
    if (!ctx->work || !ctx->equations || !ctx->equations2 || !ctx->orderparam || !ctx->orderparam2 || !ctx->transformed ||
        !ctx->key_values || !ctx->swap || !ctx->terms || !ctx->tmp || !ctx->column_state || !ctx->assign_eq || !ctx->assign_col ||
        !ctx->assign_value || !ctx->pivots || !ctx->is_pivot || !ctx->matrix_key || !ctx->orbit_matrix_key || !ctx->check_cache ||
        !ctx->transformed_key_cache || !ctx->projective || !ctx->flags || !ctx->minus_list || !ctx->orbit_pairs || !ctx->nz_start ||
        !ctx->nz_row || !ctx->residue_slots)
        return XT_NO_MEMORY;
    for (size_t index = 0; index < slot_count; index++) ctx->residue_slots[index] = -1;
    for (int fraction = 0; fraction < ctx->fraction_count; fraction++) {
        size_t slot = residue_hash(ctx->fraction_residue + fraction * 3) & (slot_count - 1);
        while (ctx->residue_slots[slot] >= 0) slot = (slot + 1) & (slot_count - 1);
        ctx->residue_slots[slot] = fraction;
    }
    for (int matrix = 0; matrix < ctx->matrix_count; matrix++) {
        int *start = ctx->nz_start + (size_t)matrix * (dim + 1), *rows = ctx->nz_row + (size_t)matrix * square, n = 0;
        const xe *m = cat(ctx->matrices, (size_t)matrix * dim * dim);
        for (size_t col = 0; col < dim; col++) {
            start[col] = (int)((size_t)matrix * square) + n;
            for (size_t source = 0; source < dim; source++)
                if (!ez(cat(m, source * dim + col))) rows[n++] = (int)source;
        }
        start[dim] = (int)((size_t)matrix * square) + n;
    }
    return XT_OK;
}

static int push_node(Node **queue, size_t *capacity, size_t *tail, int free_rows, int dim, int pair_count, const xe *node_orderparam,
                     const xe *row_orderparam, const signed char *flags)
{
    if (*tail == *capacity) {
        size_t grown_capacity = *capacity ? *capacity * 2 : 64;
        Node *grown = realloc(*queue, grown_capacity * sizeof(Node));
        if (!grown) return XT_NO_MEMORY;
        *queue = grown;
        *capacity = grown_capacity;
    }
    size_t size = (size_t)free_rows * dim;
    Node *node = &(*queue)[*tail];
    node->free_rows = free_rows;
    node->node_orderparam = elements(size);
    node->row_orderparam = elements(size);
    node->flags = malloc((size_t)pair_count + 1);
    *tail += 1;
    if (!node->node_orderparam || !node->row_orderparam || !node->flags) return XT_NO_MEMORY;
    memcpy(node->node_orderparam, node_orderparam, size * XS * sizeof(xe));
    memcpy(node->row_orderparam, row_orderparam, size * XS * sizeof(xe));
    memcpy(node->flags, flags, (size_t)pair_count);
    return XT_OK;
}

typedef int (*xt_node_fn)(void *user, int free_rows, const i128 *row_orderparam);

/* matrices: matrix_count x dim x dim elements of f ([source_row][col]), the
 * rest as in isotropy_traverse.  on_node gets each node's representative. */
int exact_traverse(const cy_field *f, int dim, int space_count, int fraction_count, const i128 *matrices, int matrix_count,
                   const int *pair_matrix, const int *space_product, const int *space_fraction, const int *fraction_residue,
                   int residue_modulus, xt_node_fn on_node, void *user)
{
    F = f;
    XS = f->phi + 1;
    inverse_cache_reset();
    Context ctx;
    memset(&ctx, 0, sizeof ctx);
    ctx.dim = dim;
    ctx.space_count = space_count;
    ctx.fraction_count = fraction_count;
    ctx.pair_count = space_count * fraction_count;
    ctx.matrix_count = matrix_count;
    ctx.matrices = matrices;
    ctx.pair_matrix = pair_matrix;
    ctx.space_product = space_product;
    ctx.space_fraction = space_fraction;
    ctx.fraction_residue = fraction_residue;
    ctx.residue_modulus = residue_modulus;
    Node *queue = NULL;
    size_t queue_capacity = 0, head = 0, tail = 0;
    IntList keys = {0}, candidates = {0};
    signed char *candidate_flags = malloc((size_t)ctx.pair_count + 1);
    xe *representative = elements((size_t)dim * dim), *next_orderparam = elements((size_t)dim * dim), *first = elements((size_t)dim * dim);
    unsigned char *matrix_done = malloc((size_t)matrix_count + 1);
    int status = context_init(&ctx);
    if (!status && (!candidate_flags || !representative || !next_orderparam || !first || !matrix_done)) status = XT_NO_MEMORY;
    if (status) goto done;
    for (int row = 0; row < dim; row++)
        for (int col = 0; col < dim; col++) eint(at(first, (size_t)row * dim + col), row == col);
    if ((status = orbit(&ctx, dim, first, 0, &keys, representative))) goto done;
    for (size_t index = 0; index < keys.used; index++) ctx.seen[keys.data[index]] = 1;
    for (int pair = 0; pair < ctx.pair_count; pair++) {
        int preserves;
        if ((status = orderparam_check(&ctx, dim, representative, matrix_of_pair(&ctx, pair), &preserves))) goto done;
        ctx.flags[pair] = preserves ? -1 : 0;
    }
    if ((status = push_node(&queue, &queue_capacity, &tail, dim, dim, ctx.pair_count, first, representative, ctx.flags))) goto done;
    int projective_ready = 0;
    while (head < tail) {
        Node node = queue[head++];
        int decision = on_node(user, node.free_rows, node.row_orderparam);
        if (decision < 0) {
            status = -decision;
            goto done;
        }
        if (decision > 0) goto done;
        memcpy(candidate_flags, node.flags, (size_t)ctx.pair_count);
        candidates.used = 0;
        int marked_count = 0;
        for (int pair = 0; pair < ctx.pair_count; pair++)
            if (candidate_flags[pair] == -1) ctx.minus_list[marked_count++] = pair;
        for (int pair = 0; pair < ctx.pair_count; pair++) {
            if (candidate_flags[pair] != 0) continue;
            for (int index = 0; index < marked_count; index++) {
                int product;
                if ((status = pair_product(&ctx, pair, ctx.minus_list[index], &product))) goto done;
                candidate_flags[product] = 1;
            }
            if ((status = list_push(&candidates, pair))) goto done;
        }
        memset(matrix_done, 0, (size_t)ctx.matrix_count);
        for (size_t index = 0; index < candidates.used; index++) {
            int pair = candidates.data[index];
            int found, next_free, key;
            if (matrix_done[ctx.pair_matrix[pair]]) continue;
            matrix_done[ctx.pair_matrix[pair]] = 1;
            if ((status = fixed_space_intersection(&ctx, node.free_rows, node.node_orderparam, matrix_of_pair(&ctx, pair), &found, &next_free,
                                                   next_orderparam, &key)))
                goto done;
            if (!found || ctx.seen[key]) continue;
            if (!projective_ready) {
                if ((status = projective_pairs(&ctx))) goto done;
                projective_ready = 1;
            }
            if ((status = orbit(&ctx, next_free, next_orderparam, 1, &keys, representative))) goto done;
            int overlaps = 0;
            for (size_t item = 0; item < keys.used; item++)
                if (ctx.seen[keys.data[item]]) {
                    overlaps = 1;
                    break;
                }
            if (overlaps) continue;
            if ((status = close_child_flags(&ctx, next_free, next_orderparam, node.flags, pair))) goto done;
            for (size_t item = 0; item < keys.used; item++) ctx.seen[keys.data[item]] = 1;
            if ((status = push_node(&queue, &queue_capacity, &tail, next_free, dim, ctx.pair_count, next_orderparam, representative,
                                    ctx.flags)))
                goto done;
        }
    }
done:
    for (size_t index = 0; index < tail; index++) {
        free(queue[index].node_orderparam);
        free(queue[index].row_orderparam);
        free(queue[index].flags);
    }
    free(queue);
    free(keys.data);
    free(candidates.data);
    free(candidate_flags);
    free(representative);
    free(next_orderparam);
    free(first);
    free(matrix_done);
    context_free(&ctx);
    return status;
}

/* A node's representative as doubles, and its orthogonalized form
 * (orthogonal_orderparam_: kept when already orthogonal, else one reverse
 * Gram-Schmidt pass), computed exactly and correctly rounded. */
int exact_rows_to_doubles(const cy_field *f, int rows, int cols, const i128 *in, double *plain, double *orthogonal)
{
    F = f;
    XS = f->phi + 1;
    size_t cells = (size_t)rows * cols;
    xe *v = elements(cells), *t = elements(4);
    int status = v && t ? XT_OK : XT_NO_MEMORY;
    if (!status) memcpy(v, in, cells * XS * sizeof(xe));
    for (size_t i = 0; i < cells && !status; i++) {
        cy x;
        ld(cat(in, i), &x);
        if (cy_to_double(F, &x, &plain[i])) status = XT_ARITHMETIC;
    }
    /* dot(a, b) into t[0] */
#define DOT(a, b)                                                                     \
    do {                                                                              \
        eint(t, 0);                                                                   \
        for (int col = 0; col < cols && !status; col++) {                             \
            const xe *x_ = cat(v, (size_t)(a) * cols + col), *y_ = cat(v, (size_t)(b) * cols + col); \
            if (ez(x_) || ez(y_)) continue;                                           \
            status = emul(x_, y_, at(t, 1));                                          \
            if (!status) status = eadd(t, at(t, 1), t);                               \
        }                                                                             \
    } while (0)
    int orthogonal_already = 1;
    for (int row = 1; row < rows && orthogonal_already && !status; row++)
        for (int prior = 0; prior < row && orthogonal_already && !status; prior++) {
            DOT(row, prior);
            if (!ez(t)) orthogonal_already = 0;
        }
    if (!orthogonal_already)
        for (int row = rows - 2; row >= 0 && !status; row--)
            for (int other = row + 1; other < rows && !status; other++) {
                DOT(other, row);
                ecopy(at(t, 2), t);
                DOT(other, other);
                if (status || ez(t)) continue;
                status = ediv(at(t, 2), t, at(t, 3));
                for (int col = 0; col < cols && !status; col++) {
                    const xe *y = cat(v, (size_t)other * cols + col);
                    if (ez(y)) continue;
                    status = emul(at(t, 3), y, at(t, 1));
                    if (!status) status = esub(at(v, (size_t)row * cols + col), at(t, 1), at(v, (size_t)row * cols + col));
                }
            }
#undef DOT
    for (size_t i = 0; i < cells && !status; i++) {
        cy x;
        ld(cat(v, i), &x);
        if (cy_to_double(F, &x, &orthogonal[i])) status = XT_ARITHMETIC;
    }
    free(v);
    free(t);
    return status;
}
