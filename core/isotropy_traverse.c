/* Breadth-first traversal of the isotropy-subspace lattice for dynamic rows.
 *
 * Every floating-point operation follows the Python kernels it replaced in
 * the same order, so the emitted subspaces are bit-identical. Build with
 * -ffp-contract=off.
 * Subspaces are dense free x dim row-major matrices; equation systems are
 * dense eq_count x dim row-major matrices. */

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "irreps.h"
#include "pynum.h"
#include "kernels.h"
#include "magnetic.h"

enum {
    ISO_OK = 0,
    ISO_NO_MEMORY = 1,
    ISO_EQS_NO_PIVOT = 2,
    ISO_EQS_REDUCED_NO_PIVOT = 3,
    ISO_EQS_DUPLICATE_PIVOT = 4,
    ISO_EQS_INVALID_COLUMN_STATE = 5,
    ISO_EQS_INCONSISTENT_COUNT = 6,
    ISO_E2O_NO_PIVOT = 7,
    ISO_E2O_ZERO_PIVOT = 8,
    ISO_E2O_INCONSISTENT = 9,
    ISO_FIXED_FREE_CHANGED = 10,
    ISO_FRACTION_NOT_FOUND = 11,
};

typedef int (*NodeCallback)(void *user, int free_rows, const double *row_orderparam);

#define TRY(expr) do { int status_ = (expr); if (status_) return status_; } while (0)

static int nonzero(double value) { return 1e-6 < fabs(value); }

typedef struct { double *data; size_t used, capacity; } DoubleArena;
typedef struct { int *data; size_t used, capacity; } IntList;

static int arena_reserve(DoubleArena *arena, size_t extra)
{
    if (arena->used + extra <= arena->capacity) return ISO_OK;
    size_t capacity = arena->capacity ? arena->capacity : 4096;
    while (capacity < arena->used + extra) capacity *= 2;
    double *data = realloc(arena->data, capacity * sizeof(double));
    if (!data) return ISO_NO_MEMORY;
    arena->data = data;
    arena->capacity = capacity;
    return ISO_OK;
}

static int list_push(IntList *list, int value)
{
    if (list->used == list->capacity) {
        size_t capacity = list->capacity ? list->capacity * 2 : 64;
        int *data = realloc(list->data, capacity * sizeof(int));
        if (!data) return ISO_NO_MEMORY;
        list->data = data;
        list->capacity = capacity;
    }
    list->data[list->used++] = value;
    return ISO_OK;
}

/* Equation keys, interned so sets become flags over key ids. Python compares
 * key tuples with float ==, so -0.0 and 0.0 hash and compare equal. */
typedef struct {
    DoubleArena values;
    size_t *offset;
    int *eq_count, *length;
    uint64_t *hash;
    int count, capacity;
    int *slots;
    size_t slot_count;
} KeyTable;

/* Doubles of small rationals have zero low mantissa bits; spread the high
 * bits into the low ones the slot masks keep. */
static uint64_t mix64(uint64_t hash)
{
    hash ^= hash >> 30;
    hash *= 0xbf58476d1ce4e5b9ULL;
    hash ^= hash >> 27;
    hash *= 0x94d049bb133111ebULL;
    return hash ^ (hash >> 31);
}

static uint64_t key_hash(int eq_count, const double *values, int length)
{
    uint64_t hash = 1469598103934665603ULL ^ (uint64_t)eq_count;
    for (int index = 0; index < length; index++) {
        double value = values[index] == 0.0 ? 0.0 : values[index];
        uint64_t bits;
        memcpy(&bits, &value, sizeof bits);
        hash = (hash ^ bits) * 1099511628211ULL;
    }
    return mix64(hash);
}

static int key_equal(const KeyTable *table, int id, int eq_count, const double *values, int length)
{
    if (table->eq_count[id] != eq_count || table->length[id] != length) return 0;
    const double *stored = table->values.data + table->offset[id];
    for (int index = 0; index < length; index++) if (!(stored[index] == values[index])) return 0;
    return 1;
}

static int key_grow_slots(KeyTable *table)
{
    size_t slot_count = table->slot_count ? table->slot_count * 2 : 1024;
    int *slots = malloc(slot_count * sizeof(int));
    if (!slots) return ISO_NO_MEMORY;
    for (size_t index = 0; index < slot_count; index++) slots[index] = -1;
    for (int id = 0; id < table->count; id++) {
        size_t slot = table->hash[id] & (slot_count - 1);
        while (slots[slot] >= 0) slot = (slot + 1) & (slot_count - 1);
        slots[slot] = id;
    }
    free(table->slots);
    table->slots = slots;
    table->slot_count = slot_count;
    return ISO_OK;
}

static int key_intern(KeyTable *table, int eq_count, const double *values, int length, int *id_out)
{
    if ((size_t)(table->count + 1) * 2 > table->slot_count) TRY(key_grow_slots(table));
    uint64_t hash = key_hash(eq_count, values, length);
    size_t slot = hash & (table->slot_count - 1);
    while (table->slots[slot] >= 0) {
        int id = table->slots[slot];
        if (table->hash[id] == hash && key_equal(table, id, eq_count, values, length)) { *id_out = id; return ISO_OK; }
        slot = (slot + 1) & (table->slot_count - 1);
    }
    if (table->count == table->capacity) {
        int capacity = table->capacity ? table->capacity * 2 : 1024;
        size_t *offset = realloc(table->offset, (size_t)capacity * sizeof(size_t));
        if (!offset) return ISO_NO_MEMORY;
        table->offset = offset;
        int *counts = realloc(table->eq_count, (size_t)capacity * sizeof(int));
        if (!counts) return ISO_NO_MEMORY;
        table->eq_count = counts;
        int *lengths = realloc(table->length, (size_t)capacity * sizeof(int));
        if (!lengths) return ISO_NO_MEMORY;
        table->length = lengths;
        uint64_t *hashes = realloc(table->hash, (size_t)capacity * sizeof(uint64_t));
        if (!hashes) return ISO_NO_MEMORY;
        table->hash = hashes;
        table->capacity = capacity;
    }
    TRY(arena_reserve(&table->values, (size_t)length));
    int id = table->count++;
    table->offset[id] = table->values.used;
    if (length) memcpy(table->values.data + table->values.used, values, (size_t)length * sizeof(double));
    table->values.used += (size_t)length;
    table->eq_count[id] = eq_count;
    table->length[id] = length;
    table->hash[id] = hash;
    table->slots[slot] = id;
    *id_out = id;
    return ISO_OK;
}

/* xrowop2 on a dense rows x cols row-major matrix, in place. */
static void xrowop2(double *matrix, int rows, int cols, double *swap)
{
    int pivot_col = 0;
    for (int pivot_row = 0; pivot_row < rows; pivot_row++) {
        int reduced = 0;
        while (pivot_col < cols) {
            int swap_row = -1;
            for (int row = pivot_row; row < rows; row++) {
                if (fabs(matrix[row * cols + pivot_col]) > 1e-6) { swap_row = row; break; }
            }
            if (swap_row < 0) { pivot_col += 1; continue; }
            if (swap_row != pivot_row) {
                memcpy(swap, matrix + pivot_row * cols, (size_t)cols * sizeof(double));
                memcpy(matrix + pivot_row * cols, matrix + swap_row * cols, (size_t)cols * sizeof(double));
                memcpy(matrix + swap_row * cols, swap, (size_t)cols * sizeof(double));
            }
            double *pivot_values = matrix + pivot_row * cols;
            double pivot = pivot_values[pivot_col];
            for (int col = pivot_col; col < cols; col++) pivot_values[col] /= pivot;
            for (int row = 0; row < rows; row++) {
                if (row == pivot_row) continue;
                double *row_values = matrix + row * cols;
                double factor = row_values[pivot_col];
                if (fabs(factor) <= 1e-6) continue;
                for (int col = 0; col < cols; col++) {
                    if (col == pivot_col) continue;
                    row_values[col] -= factor * pivot_values[col];
                }
                row_values[pivot_col] = 0.0;
            }
            pivot_col += 1;
            reduced = 1;
            break;
        }
        if (!reduced) break;
    }
}

typedef struct {
    int dim, space_count, fraction_count, pair_count, matrix_count;
    const double *matrices;          /* matrix_count x dim x dim, [source_row][col] */
    const int *pair_matrix;          /* pair_count; matrices numbered by first appearance */
    const int *space_product;        /* space_count^2 x 2, 1-based (space, fraction) */
    const int *space_fraction;       /* space_count x fraction_count, 1-based */
    const int *fraction_residue;     /* fraction_count x 3, modulo residue_modulus */
    int residue_modulus;
    int *residue_slots;
    size_t residue_slot_count;
    int *projective;                 /* matrix -> first sign-equivalent matrix */
    int *orbit_pairs;
    int orbit_pair_count;
    KeyTable keys;
    unsigned char *seen;
    size_t seen_capacity;
    double *work, *equations, *equations2, *orderparam, *orderparam2, *transformed, *key_values, *swap, *terms;
    int *column_state, *assign_eq, *assign_col, *pivots;
    double *assign_value;
    unsigned char *is_pivot;
    int *matrix_key, *orbit_matrix_key, *check_cache, *transformed_key_cache;
    signed char *flags;
    int *minus_list;
    int *nz_offset;                  /* matrix -> first column start in nz_starts, -1 until built */
    IntList nz_starts, nz_rows;
    DoubleArena nz_values;
} Context;

static const double *matrix_of_pair(const Context *ctx, int pair)
{
    return ctx->matrices + (size_t)ctx->pair_matrix[pair] * ctx->dim * ctx->dim;
}

/* orderparam_to_active_eqs: equations[eq][display_col], eq_count x cols. */
static int orderparam_to_eqs(Context *ctx, int free_rows, const double *orderparam, double *equations, int *eq_count_out)
{
    int cols = ctx->dim;
    double *reduced = ctx->work;
    memcpy(reduced, orderparam, (size_t)(free_rows * cols) * sizeof(double));
    for (int row = 0; row < free_rows; row++) {
        int pivot = -1;
        for (int col = 0; col < cols; col++) if (nonzero(reduced[row * cols + col])) { pivot = col; break; }
        if (pivot < 0) return ISO_EQS_NO_PIVOT;
        double pivot_value = reduced[row * cols + pivot];
        for (int col = pivot; col < cols; col++) reduced[row * cols + col] /= pivot_value;
        for (int other = 0; other < free_rows; other++) {
            if (other == row) continue;
            double factor = reduced[other * cols + pivot];
            if (nonzero(factor)) {
                for (int col = 0; col < cols; col++) reduced[other * cols + col] -= factor * reduced[row * cols + col];
            }
        }
    }
    int *column_state = ctx->column_state;
    memset(column_state, 0, (size_t)cols * sizeof(int));
    int eq_count = 0, assigned = 0;
    for (int row = 0; row < free_rows; row++) {
        int pivot = -1;
        for (int col = 0; col < cols; col++) if (nonzero(reduced[row * cols + col])) { pivot = col; break; }
        if (pivot < 0) return ISO_EQS_REDUCED_NO_PIVOT;
        if (column_state[pivot] != 0) return ISO_EQS_DUPLICATE_PIVOT;
        column_state[pivot] = -1;
        for (int col = pivot + 1; col < cols; col++) {
            double value = reduced[row * cols + col];
            if (!nonzero(value)) continue;
            int existing = column_state[col];
            if (existing < 1) {
                if (existing != 0) return ISO_EQS_INVALID_COLUMN_STATE;
                eq_count += 1;
                column_state[col] = eq_count;
                ctx->assign_eq[assigned] = eq_count - 1; ctx->assign_col[assigned] = col; ctx->assign_value[assigned++] = 1.0;
            }
            ctx->assign_eq[assigned] = column_state[col] - 1; ctx->assign_col[assigned] = pivot; ctx->assign_value[assigned++] = -value;
        }
    }
    for (int col = 0; col < cols; col++) {
        if (column_state[col] == 0) {
            eq_count += 1;
            ctx->assign_eq[assigned] = eq_count - 1; ctx->assign_col[assigned] = col; ctx->assign_value[assigned++] = 1.0;
        }
    }
    if (cols - free_rows != eq_count) return ISO_EQS_INCONSISTENT_COUNT;
    memset(equations, 0, (size_t)(eq_count * cols) * sizeof(double));
    for (int index = 0; index < assigned; index++) {
        equations[ctx->assign_eq[index] * cols + (cols - ctx->assign_col[index] - 1)] = ctx->assign_value[index];
    }
    *eq_count_out = eq_count;
    return ISO_OK;
}

/* eqs_to_orderparam: orderparam[free][cols] from eq_count x cols equations. */
static int eqs_to_orderparam(Context *ctx, int cols, int eq_count, const double *equations, double *orderparam, int *free_out)
{
    memset(ctx->is_pivot, 0, (size_t)cols);
    for (int row = 0; row < eq_count; row++) {
        int pivot = -1;
        for (int col = 0; col < cols; col++) if (nonzero(equations[row * cols + col])) { pivot = col; break; }
        if (pivot < 0) return ISO_E2O_NO_PIVOT;
        ctx->pivots[row] = pivot;
        ctx->is_pivot[pivot] = 1;
    }
    int free_count = 0;
    for (int source_col = cols - 1; source_col >= 0; source_col--) {
        if (ctx->is_pivot[source_col]) continue;
        double *out = orderparam + free_count * cols;
        free_count += 1;
        memset(out, 0, (size_t)cols * sizeof(double));
        out[cols - source_col - 1] = 1.0;
        for (int eq_row = 0; eq_row < eq_count; eq_row++) {
            int pivot_col = ctx->pivots[eq_row];
            double pivot_value = equations[eq_row * cols + pivot_col];
            if (!nonzero(pivot_value)) return ISO_E2O_ZERO_PIVOT;
            out[cols - pivot_col - 1] = -(equations[eq_row * cols + source_col] / pivot_value);
        }
    }
    if (cols - eq_count != free_count) return ISO_E2O_INCONSISTENT;
    *free_out = free_count;
    return ISO_OK;
}

static int equations_are_reduced(int cols, int eq_count, const double *equations)
{
    int previous = -1;
    for (int row = 0; row < eq_count; row++) {
        int pivot = -1;
        for (int col = 0; col < cols; col++) {
            double value = equations[row * cols + col];
            if (fabs(value) > 1e-6) {
                pivot = col;
                if (fabs(value - 1.0) > 1e-6) return 0;
                break;
            }
        }
        if (pivot < 0) return 0;
        if (previous >= 0 && pivot <= previous) return 0;
        previous = pivot;
        for (int other = 0; other < eq_count; other++) {
            if (other != row && fabs(equations[other * cols + pivot]) > 1e-6) return 0;
        }
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
            return ISO_OK;
        }
        slot = (slot + 1) & (ctx->residue_slot_count - 1);
    }
    return ISO_FRACTION_NOT_FOUND;
}

/* Pair p = space * fraction_count + fraction, both 0-based. */
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
        long long sum = (long long)ctx->fraction_residue[left_fraction * 3 + axis]
            + ctx->fraction_residue[correction * 3 + axis]
            + ctx->fraction_residue[rotated * 3 + axis];
        residue[axis] = (int)(sum % ctx->residue_modulus);
    }
    int fraction;
    TRY(fraction_of_residue(ctx, residue, &fraction));
    *product = (space[0] - 1) * fraction_count + fraction;
    return ISO_OK;
}

/* Per matrix column, its nonzero (source_row, value) entries in ascending
 * source_row order as absolute ranges of nz_rows/nz_values; cached for the
 * traversal's matrices, rebuilt (and dropped by the caller) otherwise. */
static int column_nonzeros(Context *ctx, const double *matrix, const int **starts)
{
    int dim = ctx->dim;
    ptrdiff_t index = -1;
    if (ctx->nz_offset && matrix >= ctx->matrices && matrix < ctx->matrices + (size_t)ctx->matrix_count * dim * dim)
        index = (matrix - ctx->matrices) / (dim * dim);
    if (index >= 0 && ctx->nz_offset[index] >= 0) {
        *starts = ctx->nz_starts.data + ctx->nz_offset[index];
        return ISO_OK;
    }
    int offset = (int)ctx->nz_starts.used;
    for (int col = 0; col < dim; col++) {
        TRY(list_push(&ctx->nz_starts, (int)ctx->nz_rows.used));
        for (int source_row = 0; source_row < dim; source_row++) {
            double value = matrix[source_row * dim + col];
            if (value == 0.0) continue;
            TRY(list_push(&ctx->nz_rows, source_row));
            TRY(arena_reserve(&ctx->nz_values, 1));
            ctx->nz_values.data[ctx->nz_values.used++] = value;
        }
    }
    TRY(list_push(&ctx->nz_starts, (int)ctx->nz_rows.used));
    if (index >= 0) ctx->nz_offset[index] = offset;
    *starts = ctx->nz_starts.data + offset;
    return ISO_OK;
}

static void drop_uncached(Context *ctx, const double *matrix, const int *starts)
{
    int dim = ctx->dim;
    if (ctx->nz_offset && matrix >= ctx->matrices && matrix < ctx->matrices + (size_t)ctx->matrix_count * dim * dim) return;
    ctx->nz_rows.used = (size_t)starts[0];
    ctx->nz_values.used = (size_t)starts[0];
    ctx->nz_starts.used = (size_t)(starts - ctx->nz_starts.data);
}

/* _transform_orderparam_sparse over the nonzero matrix entries. */
static int transform(Context *ctx, int free_rows, const double *orderparam, const double *matrix, double *out)
{
    int dim = ctx->dim;
    const int *starts;
    TRY(column_nonzeros(ctx, matrix, &starts));
    const int *rows = ctx->nz_rows.data;
    const double *values = ctx->nz_values.data;
    for (int row = 0; row < free_rows; row++) {
        for (int col = 0; col < dim; col++) {
            double total = 0.0;
            for (int k = starts[col]; k < starts[col + 1]; k++) total += values[k] * orderparam[row * dim + rows[k]];
            out[row * dim + col] = total;
        }
    }
    drop_uncached(ctx, matrix, starts);
    return ISO_OK;
}

static int orderparam_check(Context *ctx, int free_rows, const double *orderparam, const double *matrix, int *preserves)
{
    int dim = ctx->dim;
    const int *starts;
    TRY(column_nonzeros(ctx, matrix, &starts));
    const int *rows = ctx->nz_rows.data;
    const double *values = ctx->nz_values.data;
    *preserves = 1;
    for (int row = 0; row < free_rows && *preserves; row++) {
        for (int col = 0; col < dim; col++) {
            double transformed = 0.0;
            for (int k = starts[col]; k < starts[col + 1]; k++) transformed += values[k] * orderparam[row * dim + rows[k]];
            if (nonzero(transformed - orderparam[row * dim + col])) {
                *preserves = 0;
                break;
            }
        }
    }
    drop_uncached(ctx, matrix, starts);
    return ISO_OK;
}

/* _orderparam_equation_data: unrounded equations plus the interned key. */
static int equation_data(Context *ctx, int free_rows, const double *orderparam, double *equations, int *eq_count_out, int *key_out)
{
    int dim = ctx->dim, eq_count;
    TRY(orderparam_to_eqs(ctx, free_rows, orderparam, equations, &eq_count));
    double *rounded = ctx->equations2;
    for (int index = 0; index < eq_count * dim; index++) rounded[index] = python_round(equations[index], 12);
    if (!equations_are_reduced(dim, eq_count, rounded)) xrowop2(rounded, eq_count, dim, ctx->swap);
    int length = 0;
    for (int col = 0; col < dim; col++) {
        for (int row = 0; row < eq_count; row++) ctx->key_values[length++] = python_round(rounded[row * dim + col], 10);
    }
    TRY(key_intern(&ctx->keys, eq_count, ctx->key_values, length, key_out));
    if ((size_t)ctx->keys.count > ctx->seen_capacity) {
        size_t capacity = ctx->seen_capacity ? ctx->seen_capacity * 2 : 4096;
        while (capacity < (size_t)ctx->keys.count) capacity *= 2;
        unsigned char *seen = realloc(ctx->seen, capacity);
        if (!seen) return ISO_NO_MEMORY;
        memset(seen + ctx->seen_capacity, 0, capacity - ctx->seen_capacity);
        ctx->seen = seen;
        ctx->seen_capacity = capacity;
    }
    *eq_count_out = eq_count;
    return ISO_OK;
}

/* nice_orderparam(dim, free, a, b) == 2. */
static void nice_score(int dim, int free_rows, const double *orderparam, int score[5])
{
    memset(score, 0, 5 * sizeof(int));
    for (int component = 0; component < dim; component++) {
        int seen_nonzero = 0;
        for (int row = 0; row < free_rows; row++) {
            double value = orderparam[row * dim + component];
            if (!nonzero(value)) continue;
            if (nonzero(fabs(value) - 1.0)) score[1] += 1;
            if (seen_nonzero) score[0] += 1; else score[2] += 1;
            if (value < 0.0) score[3] += 1;
            score[4] += component + 1;
            seen_nonzero = 1;
        }
    }
}

static int nice_prefers_second(int dim, int free_rows, const double *first, const double *second)
{
    int first_nonzero = nonzero(first[0]), second_nonzero = nonzero(second[0]);
    if (first_nonzero && !second_nonzero) return 0;
    if (!first_nonzero && second_nonzero) return 1;
    int a[5], b[5];
    nice_score(dim, free_rows, first, a);
    nice_score(dim, free_rows, second, b);
    for (int index = 0; index < 5; index++) {
        if (a[index] < b[index]) return 0;
        if (b[index] < a[index]) return 1;
    }
    return 0;
}

static int is_full(int dim, int free_rows, const double *orderparam)
{
    if (free_rows != dim) return 0;
    for (int row = 0; row < dim; row++) {
        for (int col = 0; col < dim; col++) {
            if (fabs(orderparam[row * dim + col] - (row == col ? 1.0 : 0.0)) > 1e-10) return 0;
        }
    }
    return 1;
}

/* _orbit_keys_and_representative; orbit_matrix_key feeds child flags. */
static int orbit(Context *ctx, int free_rows, const double *orderparam, int use_projective, IntList *keys, double *representative)
{
    int dim = ctx->dim, eq_count, key;
    size_t size = (size_t)free_rows * dim;
    TRY(equation_data(ctx, free_rows, orderparam, ctx->equations, &eq_count, &key));
    keys->used = 0;
    TRY(list_push(keys, key));
    memcpy(representative, orderparam, size * sizeof(double));
    for (int matrix = 0; matrix < ctx->matrix_count; matrix++) ctx->matrix_key[matrix] = ctx->orbit_matrix_key[matrix] = -1;
    if (is_full(dim, free_rows, orderparam)) return ISO_OK;
    int pair_count = use_projective ? ctx->orbit_pair_count : ctx->pair_count;
    for (int index = 0; index < pair_count; index++) {
        int pair = use_projective ? ctx->orbit_pairs[index] : index;
        int matrix = ctx->pair_matrix[pair];
        int matrix_key = ctx->matrix_key[matrix];
        if (matrix_key < 0) {
            TRY(transform(ctx, free_rows, orderparam, matrix_of_pair(ctx, pair), ctx->orderparam2));
            TRY(equation_data(ctx, free_rows, ctx->orderparam2, ctx->equations, &eq_count, &matrix_key));
            int known = 0;
            for (size_t item = 0; item < keys->used; item++) if (keys->data[item] == matrix_key) { known = 1; break; }
            if (!known) {
                TRY(list_push(keys, matrix_key));
                memcpy(ctx->equations2, ctx->equations, (size_t)(eq_count * dim) * sizeof(double));
                xrowop2(ctx->equations2, eq_count, dim, ctx->swap);
                int next_free;
                TRY(eqs_to_orderparam(ctx, dim, eq_count, ctx->equations2, ctx->orderparam, &next_free));
                if (next_free == free_rows && nice_prefers_second(dim, free_rows, representative, ctx->orderparam)) {
                    memcpy(representative, ctx->orderparam, size * sizeof(double));
                }
            }
            ctx->matrix_key[matrix] = matrix_key;
        }
        ctx->orbit_matrix_key[matrix] = matrix_key;
    }
    if (use_projective) {
        for (int matrix = 0; matrix < ctx->matrix_count; matrix++) {
            int representative_key = ctx->orbit_matrix_key[ctx->projective[matrix]];
            if (representative_key >= 0) ctx->orbit_matrix_key[matrix] = representative_key;
        }
    }
    return ISO_OK;
}

/* _fixed_space_intersection; *found is 0 when no child results. */
static int fixed_space_intersection(Context *ctx, int free_rows, const double *orderparam, const double *matrix,
                                    int *found, int *next_free_out, double *next_orderparam, int *key_out)
{
    int dim = ctx->dim;
    *found = 0;
    TRY(transform(ctx, free_rows, orderparam, matrix, ctx->transformed));
    double *coefficients = ctx->terms + dim;   /* dim x free_rows */
    int any = 0;
    for (int component = 0; component < dim; component++) {
        for (int source_row = 0; source_row < free_rows; source_row++) {
            double value = ctx->transformed[source_row * dim + component] - orderparam[source_row * dim + component];
            any = any || nonzero(value);
            coefficients[component * free_rows + (free_rows - source_row - 1)] = value;
        }
    }
    if (!any) return ISO_OK;
    xrowop2(coefficients, dim, free_rows, ctx->swap);
    int rank = dim;
    while (rank > 0) {
        int row_nonzero = 0;
        for (int col = 0; col < free_rows; col++) if (nonzero(coefficients[(rank - 1) * free_rows + col])) { row_nonzero = 1; break; }
        if (row_nonzero) break;
        rank -= 1;
    }
    if (rank <= 0 || rank >= free_rows) return ISO_OK;
    int next_free;
    double *basis = ctx->orderparam2;   /* next_free x free_rows */
    TRY(eqs_to_orderparam(ctx, free_rows, rank, coefficients, basis, &next_free));
    double *raw = ctx->orderparam;
    for (int row = 0; row < next_free; row++) {
        for (int component = 0; component < dim; component++) {
            for (int source_row = 0; source_row < free_rows; source_row++) {
                ctx->terms[source_row] = basis[row * free_rows + source_row] * orderparam[source_row * dim + component];
            }
            raw[row * dim + component] = py_sum(ctx->terms, free_rows);
        }
    }
    int eq_count, key, canonical_free;
    TRY(equation_data(ctx, next_free, raw, ctx->equations, &eq_count, &key));
    TRY(eqs_to_orderparam(ctx, dim, eq_count, ctx->equations, next_orderparam, &canonical_free));
    if (canonical_free != next_free) return ISO_FIXED_FREE_CHANGED;
    *found = 1;
    *next_free_out = next_free;
    *key_out = key;
    return ISO_OK;
}

static void mark_minus(Context *ctx, int pair, int *minus_count)
{
    if (ctx->flags[pair] == -1) return;
    ctx->flags[pair] = -1;
    ctx->minus_list[(*minus_count)++] = pair;
}

/* _close_child_pair_flags into ctx->flags. */
static int close_child_flags(Context *ctx, int free_rows, const double *orderparam, const signed char *parent_flags, int selected_pair)
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
    return ISO_OK;
}

/* _projective_orbit_pairs: the first pair of each matrix class up to sign. */
static double first_sign(const double *values, int cells)
{
    for (int index = 0; index < cells; index++) if (values[index] != 0.0) return values[index] < 0.0 ? -1.0 : 1.0;
    return 1.0;
}

static uint64_t projective_hash(const double *values, int cells, double sign)
{
    uint64_t hash = 1469598103934665603ull;
    for (int index = 0; index < cells; index++) {
        double value = values[index] * sign;
        if (value == 0.0) value = 0.0;
        uint64_t bits;
        memcpy(&bits, &value, sizeof bits);
        hash = (hash ^ bits) * 1099511628211ull;
    }
    return mix64(hash);
}

/* Representative of each matrix up to a global sign, first in pair order.
 * A hash of the sign-normalised entries (zeros unsigned) replaces the scan
 * over earlier representatives; equal hashes are confirmed entrywise. */
static int projective_pairs(Context *ctx)
{
    int dim = ctx->dim, cells = dim * dim;
    size_t size = 1;
    while (size < (size_t)ctx->matrix_count * 2 + 2) size <<= 1;
    int *table = malloc(size * sizeof(int));
    if (!table) return ISO_NO_MEMORY;
    for (size_t slot = 0; slot < size; slot++) table[slot] = -1;
    for (int matrix = 0; matrix < ctx->matrix_count; matrix++) ctx->projective[matrix] = -1;
    ctx->orbit_pair_count = 0;
    for (int pair = 0; pair < ctx->pair_count; pair++) {
        int matrix = ctx->pair_matrix[pair];
        if (ctx->projective[matrix] >= 0) continue;
        const double *values = ctx->matrices + (size_t)matrix * cells;
        double sign = first_sign(values, cells);
        size_t slot = projective_hash(values, cells, sign) & (size - 1);
        int representative = -1;
        for (; table[slot] >= 0; slot = (slot + 1) & (size - 1)) {
            const double *other = ctx->matrices + (size_t)table[slot] * cells;
            double other_sign = first_sign(other, cells);
            int equal = 1;
            for (int cell = 0; cell < cells; cell++) if (!(values[cell] * sign == other[cell] * other_sign)) { equal = 0; break; }
            if (equal) { representative = table[slot]; break; }
        }
        if (representative < 0) {
            table[slot] = matrix;
            ctx->projective[matrix] = matrix;
            ctx->orbit_pairs[ctx->orbit_pair_count++] = pair;
        } else {
            ctx->projective[matrix] = representative;
        }
    }
    free(table);
    return ISO_OK;
}

typedef struct { int free_rows; double *node_orderparam, *row_orderparam; signed char *flags; } Node;

static void context_free(Context *ctx)
{
    free(ctx->residue_slots); free(ctx->projective); free(ctx->orbit_pairs);
    free(ctx->keys.values.data); free(ctx->keys.offset); free(ctx->keys.eq_count); free(ctx->keys.length); free(ctx->keys.hash); free(ctx->keys.slots);
    free(ctx->seen);
    free(ctx->work); free(ctx->equations); free(ctx->equations2); free(ctx->orderparam); free(ctx->orderparam2);
    free(ctx->transformed); free(ctx->key_values); free(ctx->swap); free(ctx->terms);
    free(ctx->column_state); free(ctx->assign_eq); free(ctx->assign_col); free(ctx->assign_value);
    free(ctx->pivots); free(ctx->is_pivot);
    free(ctx->matrix_key); free(ctx->orbit_matrix_key); free(ctx->check_cache); free(ctx->transformed_key_cache);
    free(ctx->flags); free(ctx->minus_list);
    free(ctx->nz_offset); free(ctx->nz_starts.data); free(ctx->nz_rows.data); free(ctx->nz_values.data);
}

static int context_init(Context *ctx)
{
    size_t dim = (size_t)ctx->dim, square = dim * dim + 1, assign = dim * (2 * dim + 1) + dim + 1;
    size_t matrices = (size_t)ctx->matrix_count + 1, pairs = (size_t)ctx->pair_count + 1;
    ctx->work = malloc(square * sizeof(double));
    ctx->equations = malloc(square * sizeof(double));
    ctx->equations2 = malloc(square * sizeof(double));
    ctx->orderparam = malloc(square * sizeof(double));
    ctx->orderparam2 = malloc(square * sizeof(double));
    ctx->transformed = malloc(square * sizeof(double));
    ctx->key_values = malloc(square * sizeof(double));
    ctx->swap = malloc(square * sizeof(double));
    ctx->terms = malloc((square + dim) * sizeof(double));
    ctx->column_state = malloc((dim + 1) * sizeof(int));
    ctx->assign_eq = malloc(assign * sizeof(int));
    ctx->assign_col = malloc(assign * sizeof(int));
    ctx->assign_value = malloc(assign * sizeof(double));
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
    ctx->nz_offset = malloc(matrices * sizeof(int));
    if (ctx->nz_offset)
        for (size_t index = 0; index < matrices; index++) ctx->nz_offset[index] = -1;
    size_t slot_count = 16;
    while (slot_count < (size_t)ctx->fraction_count * 2) slot_count *= 2;
    ctx->residue_slots = malloc(slot_count * sizeof(int));
    ctx->residue_slot_count = slot_count;
    if (!ctx->work || !ctx->equations || !ctx->equations2 || !ctx->orderparam || !ctx->orderparam2 ||
        !ctx->transformed || !ctx->key_values || !ctx->swap || !ctx->terms || !ctx->column_state ||
        !ctx->assign_eq || !ctx->assign_col || !ctx->assign_value || !ctx->pivots || !ctx->is_pivot ||
        !ctx->matrix_key || !ctx->orbit_matrix_key || !ctx->check_cache || !ctx->transformed_key_cache ||
        !ctx->projective || !ctx->flags || !ctx->minus_list || !ctx->orbit_pairs || !ctx->residue_slots) {
        return ISO_NO_MEMORY;
    }
    for (size_t index = 0; index < slot_count; index++) ctx->residue_slots[index] = -1;
    for (int fraction = 0; fraction < ctx->fraction_count; fraction++) {
        size_t slot = residue_hash(ctx->fraction_residue + fraction * 3) & (slot_count - 1);
        while (ctx->residue_slots[slot] >= 0) slot = (slot + 1) & (slot_count - 1);
        ctx->residue_slots[slot] = fraction;
    }
    return ISO_OK;
}

static int push_node(Node **queue, size_t *capacity, size_t *tail, int free_rows, int dim, int pair_count,
                     const double *node_orderparam, const double *row_orderparam, const signed char *flags)
{
    if (*tail == *capacity) {
        size_t grown_capacity = *capacity ? *capacity * 2 : 64;
        Node *grown = realloc(*queue, grown_capacity * sizeof(Node));
        if (!grown) return ISO_NO_MEMORY;
        *queue = grown;
        *capacity = grown_capacity;
    }
    size_t size = (size_t)free_rows * dim;
    Node *node = &(*queue)[*tail];
    node->free_rows = free_rows;
    node->node_orderparam = malloc((size + 1) * sizeof(double));
    node->row_orderparam = malloc((size + 1) * sizeof(double));
    node->flags = malloc((size_t)pair_count + 1);
    *tail += 1;
    if (!node->node_orderparam || !node->row_orderparam || !node->flags) return ISO_NO_MEMORY;
    memcpy(node->node_orderparam, node_orderparam, size * sizeof(double));
    memcpy(node->row_orderparam, row_orderparam, size * sizeof(double));
    memcpy(node->flags, flags, (size_t)pair_count);
    return ISO_OK;
}

/* Walk the subspace lattice; on_node receives each node in queue order and
 * returns 1 to stop, 0 to continue, or a positive status ORed with
 * NODE_FAILED. */
#define NODE_FAILED (1 << 20)

static int isotropy_traverse(
    int dim, int space_count, int fraction_count,
    const double *matrices, int matrix_count, const int *pair_matrix,
    const int *space_product, const int *space_fraction,
    const int *fraction_residue, int residue_modulus,
    NodeCallback on_node, void *user)
{
    Context ctx;
    memset(&ctx, 0, sizeof ctx);
    ctx.dim = dim; ctx.space_count = space_count; ctx.fraction_count = fraction_count;
    ctx.pair_count = space_count * fraction_count; ctx.matrix_count = matrix_count;
    ctx.matrices = matrices; ctx.pair_matrix = pair_matrix;
    ctx.space_product = space_product; ctx.space_fraction = space_fraction;
    ctx.fraction_residue = fraction_residue; ctx.residue_modulus = residue_modulus;
    Node *queue = NULL;
    size_t queue_capacity = 0, head = 0, tail = 0;
    IntList keys = {0}, candidates = {0};
    signed char *candidate_flags = malloc((size_t)ctx.pair_count + 1);
    double *representative = malloc(((size_t)dim * dim + 1) * sizeof(double));
    double *next_orderparam = malloc(((size_t)dim * dim + 1) * sizeof(double));
    double *first = calloc((size_t)dim * dim + 1, sizeof(double));
    unsigned char *matrix_done = malloc((size_t)matrix_count + 1);
    int status = context_init(&ctx);
    if (!status && (!candidate_flags || !representative || !next_orderparam || !first || !matrix_done)) status = ISO_NO_MEMORY;
    if (status) goto done;

    for (int row = 0; row < dim; row++) first[row * dim + row] = 1.0;
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
        if (decision & NODE_FAILED) { status = decision & ~NODE_FAILED; goto done; }
        if (decision > 0) goto done;

        memcpy(candidate_flags, node.flags, (size_t)ctx.pair_count);
        candidates.used = 0;
        int marked_count = 0;
        for (int pair = 0; pair < ctx.pair_count; pair++) if (candidate_flags[pair] == -1) ctx.minus_list[marked_count++] = pair;
        for (int pair = 0; pair < ctx.pair_count; pair++) {
            if (candidate_flags[pair] != 0) continue;
            for (int index = 0; index < marked_count; index++) {
                int product;
                if ((status = pair_product(&ctx, pair, ctx.minus_list[index], &product))) goto done;
                candidate_flags[product] = 1;
            }
            if ((status = list_push(&candidates, pair))) goto done;
        }

        /* Only the first candidate of each matrix can add a child: later ones
         * reproduce its intersection, whose key is then already seen. */
        memset(matrix_done, 0, (size_t)ctx.matrix_count);
        for (size_t index = 0; index < candidates.used; index++) {
            int pair = candidates.data[index];
            int found, next_free, key;
            if (matrix_done[ctx.pair_matrix[pair]]) continue;
            matrix_done[ctx.pair_matrix[pair]] = 1;
            if ((status = fixed_space_intersection(&ctx, node.free_rows, node.node_orderparam, matrix_of_pair(&ctx, pair),
                                                   &found, &next_free, next_orderparam, &key))) goto done;
            if (!found || ctx.seen[key]) continue;
            if (!projective_ready) {
                if ((status = projective_pairs(&ctx))) goto done;
                projective_ready = 1;
            }
            if ((status = orbit(&ctx, next_free, next_orderparam, 1, &keys, representative))) goto done;
            int overlaps = 0;
            for (size_t item = 0; item < keys.used; item++) if (ctx.seen[keys.data[item]]) { overlaps = 1; break; }
            if (overlaps) continue;
            if ((status = close_child_flags(&ctx, next_free, next_orderparam, node.flags, pair))) goto done;
            for (size_t item = 0; item < keys.used; item++) ctx.seen[keys.data[item]] = 1;
            if ((status = push_node(&queue, &queue_capacity, &tail, next_free, dim, ctx.pair_count,
                                    next_orderparam, representative, ctx.flags))) goto done;
        }
    }
done:
    for (size_t index = 0; index < tail; index++) {
        free(queue[index].node_orderparam); free(queue[index].row_orderparam); free(queue[index].flags);
    }
    free(queue); free(keys.data); free(candidates.data); free(candidate_flags);
    free(representative); free(next_orderparam); free(first); free(matrix_done);
    context_free(&ctx);
    return status;
}

/* Rank of a rows x cols matrix: singular values (one-sided Jacobi) above tol. */
static int matrix_rank(const double *a, int rows, int cols, double tol, double *work)
{
    int n = cols, m = rows;
    double *u = work;
    memcpy(u, a, (size_t)m * n * sizeof(double));
    for (int sweep = 0; sweep < 60; sweep++) {
        double off = 0.0;
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++) {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (int i = 0; i < m; i++) {
                    double up = u[i * n + p], uq = u[i * n + q];
                    alpha += up * up;
                    beta += uq * uq;
                    gamma += up * uq;
                }
                if (gamma == 0.0) continue;
                double scale = sqrt(alpha * beta);
                if (scale > 0.0 && fabs(gamma) / scale > off) off = fabs(gamma) / scale;
                if (scale == 0.0 || fabs(gamma) <= 1e-15 * scale) continue;
                double zeta = (beta - alpha) / (2.0 * gamma);
                double t = (zeta >= 0.0 ? 1.0 : -1.0) / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
                double c = 1.0 / sqrt(1.0 + t * t), s = c * t;
                for (int i = 0; i < m; i++) {
                    double up = u[i * n + p], uq = u[i * n + q];
                    u[i * n + p] = c * up - s * uq;
                    u[i * n + q] = s * up + c * uq;
                }
            }
        if (off <= 1e-15) break;
    }
    int rank = 0;
    for (int j = 0; j < n; j++) {
        double norm = 0.0;
        for (int i = 0; i < m; i++) norm += u[i * n + j] * u[i * n + j];
        if (sqrt(norm) > tol) rank++;
    }
    return rank;
}

/* get_isotropy_ dynamic rows: each traversal node becomes a subgroup through
 * orderparam_to_subgroup_ and id_subgroup_ (or their magnetic variants);
 * nodes that identify none are skipped. */

typedef struct {
    int dim, arms, subtract, element_count, identity_op, magnetic, parent, occurrence;
    const int64_t *k;
    const int *elements, *available;
    const double *sparse;
    int candidate_count;
    const int64_t *candidates;
    const int *candidate_ops, *candidate_signs;
    int counts[64];
    int64_t *ints;
    size_t ints_used, ints_cap;
    double *doubles;
    size_t doubles_used, doubles_cap;
} row_builder;

static int64_t *row_ints(row_builder *b, size_t count)
{
    if (b->ints_used + count > b->ints_cap) {
        size_t cap = b->ints_cap ? b->ints_cap : 4096;
        while (b->ints_used + count > cap) cap *= 2;
        int64_t *grown = realloc(b->ints, cap * sizeof *grown);
        if (!grown) return 0;
        b->ints = grown;
        b->ints_cap = cap;
    }
    int64_t *at = b->ints + b->ints_used;
    b->ints_used += count;
    return at;
}

static double *row_doubles(row_builder *b, size_t count)
{
    if (b->doubles_used + count > b->doubles_cap) {
        size_t cap = b->doubles_cap ? b->doubles_cap : 4096;
        while (b->doubles_used + count > cap) cap *= 2;
        double *grown = realloc(b->doubles, cap * sizeof *grown);
        if (!grown) return 0;
        b->doubles = grown;
        b->doubles_cap = cap;
    }
    double *at = b->doubles + b->doubles_used;
    b->doubles_used += count;
    return at;
}

#define SELECTION_CAPACITY 8192

extern int (*op_invariant_hook)(void *ctx, const int *record, int position, double sign, const double *row, int *result);
extern void *op_invariant_ctx;
int xt_invariant(void *user, const int *record, int position, double sign, const double *row, int *result);

typedef struct {
    const double *base;
    const uint64_t *rows;
    int gid, dim, arms, subtract, element_count;
    const int64_t *k;
    const int *elements;
    int64_t *words;
} xt_invariance;

static const uint64_t *exact_rows; /* the current node's rows in every embedding */
static int exact_gid;

static int build_row(void *user, int free_rows, const double *row_orderparam, const double *orthogonal)
{
    row_builder *b = user;
    if (b->occurrence && free_rows != 1) return 0;
    if (free_rows >= 64) return NODE_FAILED | ISO_NO_MEMORY;
    int dim = b->dim;
    size_t stride_values = (size_t)free_rows * 48;
    double *values = calloc(stride_values + 1, sizeof(double));
    int64_t *kept = malloc(sizeof(int64_t) * 5 * (size_t)(b->candidate_count + 1));
    int64_t *selection = b->magnetic ? malloc(sizeof(int64_t) * SELECTION_CAPACITY) : 0;
    int64_t basis[9], ordinary[14];
    int kept_count = 0, status = 0, found;
    if (!values || !kept || (b->magnetic && !selection)) {
        status = NODE_FAILED | ISO_NO_MEMORY;
        goto done;
    }
    for (int row = 0; row < free_rows; row++)
        for (int col = 0; col < dim; col++) values[row * 48 + col] = row_orderparam[row * dim + col];
    int64_t *words = malloc(sizeof(int64_t) * 8 * (size_t)dim * dim + 8);
    xt_invariance invariance = {values, exact_rows, exact_gid, dim, b->arms, b->subtract, b->element_count, b->k, b->elements, words};
    if (!words) {
        status = NODE_FAILED | ISO_NO_MEMORY;
        goto done;
    }
    op_invariant_hook = xt_invariant;
    op_invariant_ctx = &invariance;
    status = subgroup_rows_in_context(dim, b->arms, b->k, b->subtract, b->element_count, b->elements, b->sparse,
                                         b->available, b->identity_op, values, (int64_t)stride_values, free_rows,
                                         b->candidate_count, b->candidates, b->candidate_ops, b->candidate_signs, 256,
                                         basis, kept, &kept_count);
    op_invariant_hook = 0;
    op_invariant_ctx = 0;
    free(words);
    if (status) {
        status |= NODE_FAILED;
        goto done;
    }
    found = b->magnetic
                ? id_subgroup_magnetic_identify(b->parent, basis, kept_count, kept, 0, SELECTION_CAPACITY, selection)
                : id_subgroup_identify(b->parent, basis, kept_count, kept, 0, ordinary);
    if (found < 0) {
        status = NODE_FAILED | (200 - found);
        goto done;
    }
    if (!found) goto done;
    b->counts[free_rows]++;
    size_t extra = b->magnetic ? 1 + 5 * (size_t)kept_count + (size_t)found : 0;
    int64_t *out = row_ints(b, 2 + (b->magnetic ? extra : 14));
    double *matrix = row_doubles(b, (size_t)free_rows * dim);
    if (!out || !matrix) {
        status = NODE_FAILED | ISO_NO_MEMORY;
        goto done;
    }
    *out++ = free_rows;
    if (b->magnetic) {
        *out++ = found;
        memcpy(out, selection, sizeof(int64_t) * (size_t)found);
        out += found;
        *out++ = kept_count;
        memcpy(out, kept, sizeof(int64_t) * 5 * (size_t)kept_count);
    } else {
        *out++ = 14;
        memcpy(out, ordinary, sizeof ordinary);
    }
    memcpy(matrix, orthogonal, sizeof(double) * (size_t)free_rows * dim);
    if (b->occurrence && b->counts[1] == b->occurrence) status = 1;
done:
    free(values);
    free(kept);
    free(selection);
    return status;
}
static row_builder rows;

#include "cyclo.h"

int exact_image_walk(int gid, int dim, int arms, const int64_t *k, int subtract, int element_count, const int *elements, int distinct_count,
                     const int *records, const int *signs, const int *operation_distinct, size_t pairs, int space_count,
                     int fraction_count, const int *space_product, const int *space_fraction, const int *fraction_residue,
                     int residue_modulus, int (*on_node)(void *, int, const uint64_t *), void *user);
int exact_rows_to_doubles(const cy_field *f, int rows, int cols, const uint64_t *in, double *plain, double *orthogonal);
static const cy_field *rows_field;

/* exact_traverse's node callback: the representative and its orthogonal
 * form as correctly rounded doubles, then build_row. */
static int exact_build_row(void *user, int free_rows, const uint64_t *row)
{
    int dim = ((row_builder *)user)->dim;
    double *plain = malloc(sizeof(double) * ((size_t)free_rows * dim + 1)), *orthogonal = malloc(sizeof(double) * ((size_t)free_rows * dim + 1));
    int status = plain && orthogonal ? exact_rows_to_doubles(rows_field, free_rows, dim, row, plain, orthogonal) : ISO_NO_MEMORY;
    if (status) status = 400 + status;
    else {
        exact_rows = row;
        status = build_row(user, free_rows, plain, orthogonal);
        exact_rows = 0;
        if (status & NODE_FAILED) status = status & ~NODE_FAILED;
        else if (status > 0) status = -1;
    }
    free(plain);
    free(orthogonal);
    /* exact_traverse: negative stops with that status, positive stops cleanly */
    return status == -1 ? 1 : -status;
}

/* The dynamic rows of one irrep: per row, the free count, the payload length
 * n and payload (subgroup, basis, origin; or the magnetic identification
 * record), and for magnetic rows the kept operation records; the
 * orthogonalized free x dim matrices follow in isotropy_rows_doubles.
 * occurrence > 0 stops at the occurrence-th one-parameter row. */
static int isotropy_rows(int dim, int space_count, int fraction_count, int distinct_count, const int *distinct_records,
                         const int *distinct_signs, const int *operation_distinct, size_t pairs, const int *space_product,
                         const int *space_fraction, const int *fraction_residue, int residue_modulus, int gid, int param_count,
                         const int64_t *kparam, int candidate_count, const int64_t *candidates, const int *candidate_ops,
                         const int *candidate_signs, int magnetic, int parent, int occurrence)
{
    irrep_context c;
    int context_status = irrep_context_for(gid, param_count, kparam, &c);
    if (context_status) return context_status == IR_RANGE ? 103 : context_status == IR_NO_TABLES ? 106 : 105;
    rows.dim = dim;
    rows.arms = c.arms;
    rows.k = c.k;
    rows.subtract = c.subtract;
    rows.element_count = c.element_count;
    rows.elements = c.elements;
    rows.sparse = c.sparse;
    rows.available = c.available;
    rows.identity_op = c.identity_op;
    rows.candidate_count = candidate_count;
    rows.candidates = candidates;
    rows.candidate_ops = candidate_ops;
    rows.candidate_signs = candidate_signs;
    rows.magnetic = magnetic;
    rows.parent = parent;
    rows.occurrence = occurrence;
    memset(rows.counts, 0, sizeof rows.counts);
    rows.ints_used = rows.doubles_used = 0;
    exact_gid = gid;
    return exact_image_walk(gid, dim, c.arms, c.k, c.subtract, c.element_count, c.elements, distinct_count, distinct_records, distinct_signs,
                            operation_distinct, pairs, space_count, fraction_count, space_product, space_fraction, fraction_residue,
                            residue_modulus, exact_build_row, &rows);
}

enum {
    IF_TABLES = 105,
    IF_RANGE = 103,
    IF_TRAVERSAL_TABLES = 1000,
    IF_STACK = 2000,
    IF_MAGNETIC = 3000,
};

static uint64_t words_hash(const void *data, size_t size)
{
    const unsigned char *bytes = data;
    uint64_t hash = 1469598103934665603ULL;
    for (size_t i = 0; i < size; i++) hash = (hash ^ bytes[i]) * 1099511628211ULL;
    return mix64(hash);
}

/* The slot of key in an open-addressed table of ids (-1 empty) whose keys
 * live at keys + id * size; *found is the id or -1. */
static size_t hashed_slot(const int *slots, size_t mask, const unsigned char *keys, size_t size, const void *key,
                          int *found)
{
    size_t slot = words_hash(key, size) & mask;
    while (slots[slot] >= 0 && memcmp(keys + (size_t)slots[slot] * size, key, size)) slot = (slot + 1) & mask;
    *found = slots[slot];
    return slot;
}

static int *magnetic_product_table(void)
{
    static int *table;
    if (table) return table;
    int64_t count = mag_point_product_count();
    int *grown = malloc(sizeof(int) * ((size_t)count + 1));
    if (!grown) return 0;
    for (int64_t k = 0; k < count; k++) grown[k] = mag_point_product((int)(k % 144) + 1, (int)(k / 144) + 1);
    if (mag_take_error()) {
        free(grown);
        return 0;
    }
    table = grown;
    return table;
}

/* get_isotropy_ dynamic rows of a little irrep at a k parameter from Source
 * alone: the generated (magnetic) space group, the kernel of the full order
 * parameter and its get_new_fractionals_ cosets, every coset pair's
 * get_irreps_ matrix (distinct matrices numbered by first appearance) and
 * the gmlt_ product tables, then isotropy_rows. */
int isotropy_rows_for(int gid, int param_count, const int64_t *kparam, int magnetic, int occurrence)
{
    irrep_context c;
    int status = irrep_context_for(gid, param_count, kparam, &c);
    if (status) return status == IR_RANGE ? IF_RANGE : status == IR_NO_TABLES ? 106 : IF_TABLES;
    int sg = irrep_space_group(gid), dim = c.dim, lattice = geo_sg_lattice(sg);
    if (!sg || (magnetic && mag_tables())) return IF_TABLES;
    int64_t records[192][5];
    int space_count = magnetic ? magnetic_records_iso(mag_orderparam_group(sg), 192, &records[0][0])
                               : geo_sg_records(sg, records);
    if (space_count <= 0 || geo_take_error() || (magnetic && mag_take_error())) return magnetic ? IF_MAGNETIC : IF_TABLES;
    int candidate_ops[192], candidate_signs[192];
    for (int i = 0; i < space_count; i++) candidate_ops[i] = mag_point_action(magnetic, (int)records[i][4], candidate_signs + i);
    if (magnetic && mag_take_error()) return IF_MAGNETIC;
    double *full = calloc((size_t)dim * 48 + 1, sizeof(double));
    int64_t *kept = malloc(sizeof(int64_t) * (5 * (size_t)space_count + 1));
    if (!full || !kept) {
        free(full);
        free(kept);
        return ISO_NO_MEMORY;
    }
    for (int row = 0; row < dim; row++) full[row * 48 + row] = 1.0;
    int64_t basis[9];
    int kept_count = 0;
    status = subgroup_rows_in_context(dim, c.arms, c.k, c.subtract, c.element_count, c.elements, c.sparse, c.available,
                                      c.identity_op, full, (int64_t)dim * 48, dim, space_count, &records[0][0],
                                      candidate_ops, candidate_signs, 256, basis, kept, &kept_count);
    free(full);
    free(kept);
    if (status) return status;
    i128 b[9];
    for (int i = 0; i < 9; i++) b[i] = basis[i];
    i128 det128 = b[0] * (b[4] * b[8] - b[5] * b[7]) - b[1] * (b[3] * b[8] - b[5] * b[6]) + b[2] * (b[3] * b[7] - b[4] * b[6]);
    i128 abs_det = det128 < 0 ? -det128 : det128;
    if (!abs_det || abs_det > (1 << 20)) return IF_RANGE;
    int64_t *fractions = malloc(sizeof(int64_t) * 4 * (size_t)(abs_det + 1));
    if (!fractions) return ISO_NO_MEMORY;
    int fraction_count = 0;
    status = op_new_fractionals(basis, fractions, (int)abs_det + 1, &fraction_count);
    if (status) {
        free(fractions);
        return status;
    }
    size_t pairs = (size_t)space_count * fraction_count, square = (size_t)dim * dim;
    int64_t *operations = malloc(sizeof(int64_t) * 5 * (pairs + 1));
    int *distinct = malloc(sizeof(int) * 5 * (pairs + 1)), *distinct_signs = malloc(sizeof(int) * (pairs + 1));
    int *operation_distinct = malloc(sizeof(int) * (pairs + 1));
    int *space_product = malloc(sizeof(int) * 2 * ((size_t)space_count * space_count + 1));
    int *space_fraction = malloc(sizeof(int) * (pairs + 1)), *residues32 = malloc(sizeof(int) * 3 * ((size_t)fraction_count + 1));
    int64_t *residues = malloc(sizeof(int64_t) * 3 * ((size_t)fraction_count + 1));
    status = !operations || !distinct || !distinct_signs || !operation_distinct ||
                     !space_product || !space_fraction || !residues32 || !residues
                 ? ISO_NO_MEMORY
                 : 0;
    int distinct_count = 0;
    size_t slot_count = 16;
    while (slot_count < 2 * pairs + 2) slot_count *= 2;
    int *slots = malloc(sizeof(int) * slot_count);
    int64_t *distinct_keys = malloc(sizeof(int64_t) * 5 * (pairs + 1));
    if (!status && (!slots || !distinct_keys)) status = ISO_NO_MEMORY;
    for (size_t i = 0; i < slot_count && slots; i++) slots[i] = -1;
    for (size_t p = 0; p < pairs && !status; p++) {
        int64_t *op = operations + 5 * p;
        status = op_vadd(fractions + 4 * (p % fraction_count), records[p / fraction_count], op);
        if (status) break;
        int found;
        size_t slot = hashed_slot(slots, slot_count - 1, (const unsigned char *)distinct_keys, sizeof(int64_t) * 5, op,
                                  &found);
        if (found < 0) {
            found = distinct_count++;
            distinct[5 * found] = (int)p;
            memcpy(distinct_keys + 5 * (size_t)found, op, sizeof(int64_t) * 5);
            slots[slot] = found;
        }
        operation_distinct[p] = found;
    }
    for (size_t i = 0; i < slot_count && slots; i++) slots[i] = -1;
    int *stack_records = malloc(sizeof(int) * 5 * ((size_t)distinct_count + 1));
    if (!status && !stack_records) status = ISO_NO_MEMORY;
    for (int d = 0; d < distinct_count && !status; d++) {
        const int64_t *op = operations + 5 * (size_t)distinct[5 * d];
        for (int k = 0; k < 4; k++) stack_records[5 * d + k] = (int)op[k];
        stack_records[5 * d + 4] = mag_point_action(magnetic, (int)op[4], distinct_signs + d);
    }
    if (!status && magnetic && mag_take_error()) status = IF_MAGNETIC;
    int64_t setting[36], translations[192][4], rotations[192][9], point_matrices[192][9];
    int ops[192];
    int64_t rotation_den = 0;
    if (!status) {
        rotation_den = geo_ml_setting(lattice, setting);
        for (int i = 0; i < space_count; i++) {
            int sign;
            ops[i] = (int)records[i][4];
            int rotation_op = magnetic ? mag_point_action(2, ops[i], &sign) : ops[i];
            int64_t m[9], t[9];
            geo_point_op_matrix(rotation_op, m);
            for (int k = 0; k < 4; k++) translations[i][k] = records[i][k];
            for (int r = 0; r < 3; r++)
                for (int col = 0; col < 3; col++) {
                    int64_t v = 0;
                    for (int k = 0; k < 3; k++) v += setting[9 + 3 * r + k] * m[3 * k + col];
                    t[3 * r + col] = v;
                }
            for (int r = 0; r < 3; r++)
                for (int col = 0; col < 3; col++) {
                    int64_t v = 0;
                    for (int k = 0; k < 3; k++) v += t[3 * r + k] * setting[3 * k + col];
                    rotations[i][3 * r + col] = v;
                }
            if (magnetic) memset(point_matrices[i], 0, sizeof point_matrices[i]);
            else geo_point_op_matrix(ops[i], point_matrices[i]);
        }
        if (geo_take_error() || (magnetic && mag_take_error())) status = IF_TABLES;
    }
    const int none = 0;
    const int *mlt = magnetic ? magnetic_product_table() : &none;
    if (!status && !mlt) status = IF_MAGNETIC;
    int64_t modulus = 0;
    if (!status) {
        int64_t adjugate[9] = {
            basis[4] * basis[8] - basis[5] * basis[7], basis[2] * basis[7] - basis[1] * basis[8],
            basis[1] * basis[5] - basis[2] * basis[4], basis[5] * basis[6] - basis[3] * basis[8],
            basis[0] * basis[8] - basis[2] * basis[6], basis[2] * basis[3] - basis[0] * basis[5],
            basis[3] * basis[7] - basis[4] * basis[6], basis[1] * basis[6] - basis[0] * basis[7],
            basis[0] * basis[4] - basis[1] * basis[3],
        };
        int tables = traversal_tables(space_count, &translations[0][0], ops, &rotations[0][0], rotation_den, magnetic,
                                      &point_matrices[0][0], mlt, fraction_count, fractions, adjugate,
                                      (int64_t)det128, residues, &modulus, space_product, space_fraction);
        if (tables) status = IF_TRAVERSAL_TABLES + tables;
        else if (modulus >= (int64_t)1 << 31) status = IF_TRAVERSAL_TABLES + 1;
    }
    for (int i = 0; i < 3 * fraction_count && !status; i++) residues32[i] = (int)residues[i];
    if (!status)
        status = isotropy_rows(dim, space_count, fraction_count, distinct_count, stack_records, distinct_signs, operation_distinct, pairs,
                               space_product,
                               space_fraction, residues32, (int)modulus, gid, param_count, kparam, space_count,
                               &records[0][0], candidate_ops, candidate_signs, magnetic, sg, occurrence);
    free(fractions);
    free(operations);
    free(distinct);
    free(distinct_signs);
    free(operation_distinct);
    free(space_product);
    free(space_fraction);
    free(residues32);
    free(residues);
    free(stack_records);
    free(slots);
    free(distinct_keys);
    return status;
}

int64_t isotropy_rows_ints(const int64_t **data)
{
    *data = rows.ints;
    return (int64_t)rows.ints_used;
}

int64_t isotropy_rows_doubles(const double **data)
{
    *data = rows.doubles;
    return (int64_t)rows.doubles_used;
}

enum {
    DA_UNMATCHED = 20,
    DO_MAGNETIC = 21,
    DO_STACK = 1000,
    DO_SAME = 2000,
};

typedef struct {
    Context ctx;
    irrep_context irrep;
    int dim, free_rows, candidate_count, numbered_count;
    size_t cells;
    const double *values;
    const int64_t *candidates;
    int *candidate_ops, *candidate_signs;
    const double *numbered;
    int map_capacity;
    int *domain_of, *first_number;
    int domain_count, domain_capacity;
    double *domains, *stride;
    int *numbers;
    size_t *subgroup_at;
    int64_t *subgroups, *kept;
    size_t subgroups_used, subgroups_capacity;
} orbit_state;

/* Records' get_irreps_ matrices under a point action, as [source_row][col]. */
static int action_matrices(int gid, int param_count, const int64_t *kparam, int mode, int dim, int count,
                           const int64_t *records, double **out)
{
    size_t square = (size_t)dim * dim;
    double *stack = malloc(sizeof(double) * (count * square + 1));
    int *matrix_records = malloc(sizeof(int) * (5 * (size_t)count + 1)), *signs = malloc(sizeof(int) * ((size_t)count + 1));
    *out = malloc(sizeof(double) * (count * square + 1));
    int status = !stack || !matrix_records || !signs || !*out ? ISO_NO_MEMORY : ISO_OK;
    for (int r = 0; r < count && !status; r++) {
        for (int j = 0; j < 4; j++) matrix_records[5 * r + j] = (int)records[5 * r + j];
        matrix_records[5 * r + 4] = mag_point_action(mode, (int)records[5 * r + 4], signs + r);
    }
    if (!status && mode && mag_take_error()) status = DO_MAGNETIC;
    if (!status && count) {
        status = irrep_matrix_stack(gid, param_count, kparam, count, matrix_records, (int)(count * square), stack);
        if (status) status += DO_STACK;
    }
    if (!status)
        for (int r = 0; r < count; r++)
            for (int source = 0; source < dim; source++)
                for (int col = 0; col < dim; col++) {
                    double value = stack[r * square + (size_t)col * dim + source];
                    (*out)[r * square + (size_t)source * dim + col] = signs[r] < 0 ? -value : value;
                }
    free(stack);
    free(matrix_records);
    free(signs);
    return status;
}

static int orbit_key(orbit_state *s, const double *matrix, double *transformed, int *key)
{
    int eq_count;
    TRY(transform(&s->ctx, s->free_rows, s->values, matrix, transformed));
    TRY(equation_data(&s->ctx, s->free_rows, transformed, s->ctx.equations, &eq_count, key));
    if (*key < s->map_capacity) return ISO_OK;
    int capacity = s->map_capacity ? s->map_capacity : 64;
    while (capacity <= *key) capacity *= 2;
    int *domain_of = realloc(s->domain_of, sizeof(int) * capacity);
    if (!domain_of) return ISO_NO_MEMORY;
    s->domain_of = domain_of;
    int *first_number = realloc(s->first_number, sizeof(int) * capacity);
    if (!first_number) return ISO_NO_MEMORY;
    s->first_number = first_number;
    for (int i = s->map_capacity; i < capacity; i++) s->domain_of[i] = s->first_number[i] = -1;
    s->map_capacity = capacity;
    return ISO_OK;
}

static int keep_domain(orbit_state *s, const double *transformed, int number, const int64_t *basis, int kept_count)
{
    if (s->domain_count == s->domain_capacity) {
        int capacity = s->domain_capacity ? s->domain_capacity * 2 : 16;
        double *domains = realloc(s->domains, sizeof(double) * (s->cells * capacity + 1));
        int *numbers = realloc(s->numbers, sizeof(int) * capacity);
        size_t *subgroup_at = realloc(s->subgroup_at, sizeof(size_t) * capacity);
        if (domains) s->domains = domains;
        if (numbers) s->numbers = numbers;
        if (subgroup_at) s->subgroup_at = subgroup_at;
        if (!domains || !numbers || !subgroup_at) return ISO_NO_MEMORY;
        s->domain_capacity = capacity;
    }
    size_t need = s->subgroups_used + 10 + 5 * (size_t)kept_count;
    if (need > s->subgroups_capacity) {
        size_t capacity = s->subgroups_capacity ? s->subgroups_capacity : 1024;
        while (capacity < need) capacity *= 2;
        int64_t *grown = realloc(s->subgroups, sizeof(int64_t) * capacity);
        if (!grown) return ISO_NO_MEMORY;
        s->subgroups = grown;
        s->subgroups_capacity = capacity;
    }
    int d = s->domain_count++;
    memcpy(s->domains + d * s->cells, transformed, sizeof(double) * s->cells);
    s->numbers[d] = number;
    s->subgroup_at[d] = s->subgroups_used;
    int64_t *at = s->subgroups + s->subgroups_used;
    memcpy(at, basis, sizeof(int64_t) * 9);
    at[9] = kept_count;
    memcpy(at + 10, s->kept, sizeof(int64_t) * 5 * (size_t)kept_count);
    s->subgroups_used = need;
    return ISO_OK;
}

/* _domain_orbit: the distinct subspaces the operations carry the source to,
 * one per subgroup embedding; prepass numbers them by the first numbered
 * operation reaching each subspace. */
static int orbit_pass(orbit_state *s, int op_count, const double *matrices, int prepass)
{
    size_t square = (size_t)s->dim * s->dim;
    double *transformed = s->ctx.transformed;
    s->domain_count = 0;
    s->subgroups_used = 0;
    for (int i = 0; i < s->map_capacity; i++) s->domain_of[i] = s->first_number[i] = -1;
    int key;
    if (prepass)
        for (int n = 0; n < s->numbered_count; n++) {
            TRY(orbit_key(s, s->numbered + n * square, transformed, &key));
            if (s->first_number[key] < 0) s->first_number[key] = n + 1;
        }
    for (int n = 0; n < op_count; n++) {
        TRY(orbit_key(s, matrices + n * square, transformed, &key));
        if (s->domain_of[key] >= 0) continue;
        memset(s->stride, 0, sizeof(double) * 48 * (size_t)s->free_rows);
        for (int row = 0; row < s->free_rows; row++)
            memcpy(s->stride + 48 * row, transformed + (size_t)row * s->dim, sizeof(double) * s->dim);
        int64_t basis[9];
        int kept_count = 0;
        const irrep_context *c = &s->irrep;
        int status = subgroup_rows_in_context(c->dim, c->arms, c->k, c->subtract, c->element_count, c->elements,
                                              c->sparse, c->available, c->identity_op, s->stride,
                                              48 * (int64_t)s->free_rows, s->free_rows, s->candidate_count,
                                              s->candidates, s->candidate_ops, s->candidate_signs, 256, basis,
                                              s->kept, &kept_count);
        if (status) return status;
        int equivalent = -1;
        for (int d = 0; d < s->domain_count && equivalent < 0; d++) {
            const int64_t *other = s->subgroups + s->subgroup_at[d];
            int same = same_subgroup(basis, kept_count, s->kept, other, (int)other[9], other + 10);
            if (same < 0) return DO_SAME - same;
            if (same) equivalent = d;
        }
        if (equivalent >= 0) {
            s->domain_of[key] = equivalent;
            continue;
        }
        s->domain_of[key] = s->domain_count;
        int number = s->first_number[key] > 0 ? s->first_number[key] : n + 1;
        TRY(keep_domain(s, transformed, number, basis, kept_count));
    }
    return ISO_OK;
}

/* Whether [left; other] (each rows x cols) provably has rank above `rank`
 * as matrix_rank counts it: `rank` rows of left (pivoted Gram-Schmidt) plus
 * one row of other whose Gram matrix exceeds 1e-8 I (Cholesky) give the
 * stacked matrix a (rank + 1)-th singular value of at least 1e-4, and any
 * orthogonal column rotation then keeps more than `rank` of its at most 48
 * columns above 1e-6.  work: certify_size(rows, cols) doubles. */
static size_t certify_size(int rows, int cols) { return (size_t)(3 * rows + 3) * cols + (size_t)(rows + 2) * (rows + 2); }

static int certified_rank_increase(const double *left, const double *other, int rows, int cols, int rank, double *work)
{
    if (rank < 0 || rank > rows || cols > 48 || rows > 48) return 0;
    double *q = work, *chosen = q + (size_t)rows * cols, *residual = chosen + (size_t)(rows + 1) * cols;
    double *gram = residual + (size_t)rows * cols;
    int used[48] = {0};
    memcpy(residual, left, sizeof(double) * (size_t)rows * cols);
    for (int k = 0; k < rank; k++) {
        int best = -1;
        double best_norm = 0.0;
        for (int i = 0; i < rows; i++) {
            if (used[i]) continue;
            double norm = 0.0;
            for (int c = 0; c < cols; c++) norm += residual[i * cols + c] * residual[i * cols + c];
            if (norm > best_norm) {
                best_norm = norm;
                best = i;
            }
        }
        if (best < 0) return 0;
        used[best] = 1;
        memcpy(chosen + (size_t)k * cols, left + (size_t)best * cols, sizeof(double) * cols);
        double scale = 1.0 / sqrt(best_norm);
        for (int c = 0; c < cols; c++) q[k * cols + c] = residual[best * cols + c] * scale;
        for (int i = 0; i < rows; i++) {
            if (used[i]) continue;
            double dot = 0.0;
            for (int c = 0; c < cols; c++) dot += residual[i * cols + c] * q[k * cols + c];
            for (int c = 0; c < cols; c++) residual[i * cols + c] -= dot * q[k * cols + c];
        }
    }
    int n = rank + 1;
    for (int i = 0; i < rows; i++) {
        memcpy(chosen + (size_t)rank * cols, other + (size_t)i * cols, sizeof(double) * cols);
        for (int a = 0; a < n; a++)
            for (int b = 0; b <= a; b++) {
                double dot = 0.0;
                for (int c = 0; c < cols; c++) dot += chosen[a * cols + c] * chosen[b * cols + c];
                gram[a * n + b] = dot - (a == b ? 1e-8 : 0.0);
            }
        int positive = 1;
        for (int a = 0; a < n && positive; a++) {
            double diagonal = gram[a * n + a];
            for (int t = 0; t < a; t++) diagonal -= gram[a * n + t] * gram[a * n + t];
            if (!(diagonal > 0.0)) {
                positive = 0;
                break;
            }
            double pivot = sqrt(diagonal);
            gram[a * n + a] = pivot;
            for (int b = a + 1; b < n; b++) {
                double value = gram[b * n + a];
                for (int t = 0; t < a; t++) value -= gram[b * n + t] * gram[a * n + t];
                gram[b * n + a] = value / pivot;
            }
        }
        if (positive) return 1;
    }
    return 0;
}

/* _domain_actions: per operation, the domain every domain is carried to;
 * subspaces outside the known keys match a domain of equal span. */
static int domain_actions(orbit_state *s, int op_count, const double *matrices, int *out)
{
    Context *ctx = &s->ctx;
    int dim = s->dim, free_rows = s->free_rows, domain_count = s->domain_count, status = ISO_OK;
    size_t cells = s->cells;
    int *ranks = 0;
    double *transformed = malloc((cells + 1) * sizeof(double));
    double *stacked = malloc((2 * cells + 1) * sizeof(double));
    double *work = malloc((2 * cells + 1) * sizeof(double));
    double *certify = malloc(sizeof(double) * certify_size(free_rows, dim));
    if (!transformed || !stacked || !work || !certify) status = ISO_NO_MEMORY;
    for (int op = 0; op < op_count && !status; op++) {
        const double *matrix = matrices + (size_t)op * dim * dim;
        int repeated = -1;
        for (int prev = 0; prev < op && repeated < 0; prev++) {
            const double *other = matrices + (size_t)prev * dim * dim;
            int same = 1;
            for (int i = 0; i < dim * dim && same; i++) same = other[i] == matrix[i];
            if (same) repeated = prev;
        }
        if (repeated >= 0) {
            memcpy(out + (size_t)op * domain_count, out + (size_t)repeated * domain_count, sizeof(int) * (size_t)domain_count);
            continue;
        }
        for (int d = 0; d < domain_count && !status; d++) {
            status = transform(ctx, free_rows, s->domains + (size_t)d * cells, matrix, transformed);
            if (status) break;
            int eq_count, key;
            status = equation_data(ctx, free_rows, transformed, ctx->equations, &eq_count, &key);
            if (status) break;
            int target = key < s->map_capacity ? s->domain_of[key] : -1;
            if (target < 0) {
                if (!ranks) {
                    ranks = malloc(sizeof(int) * (size_t)(domain_count + 1));
                    if (!ranks) {
                        status = ISO_NO_MEMORY;
                        break;
                    }
                    for (int e = 0; e < domain_count; e++)
                        ranks[e] = matrix_rank(s->domains + (size_t)e * cells, free_rows, dim, 1e-6, work);
                }
                int left_rank = matrix_rank(transformed, free_rows, dim, 1e-6, work);
                for (int e = 0; e < domain_count && target < 0; e++) {
                    if (ranks[e] != left_rank) continue;
                    if (certified_rank_increase(transformed, s->domains + (size_t)e * cells, free_rows, dim, left_rank,
                                                certify))
                        continue;
                    memcpy(stacked, transformed, cells * sizeof(double));
                    memcpy(stacked + cells, s->domains + (size_t)e * cells, cells * sizeof(double));
                    if (matrix_rank(stacked, 2 * free_rows, dim, 1e-6, work) == left_rank) target = e;
                }
                if (target < 0) {
                    status = DA_UNMATCHED;
                    break;
                }
            }
            out[(size_t)op * domain_count + d] = target;
        }
    }
    free(transformed);
    free(stacked);
    free(work);
    free(certify);
    free(ranks);
    return status;
}

static row_builder orbit_rows;

static int orbit_emit(orbit_state *s, int generator_count, const double *generators)
{
    int d_count = s->domain_count;
    int *actions = malloc(sizeof(int) * ((size_t)generator_count * d_count + 1));
    if (!actions) return ISO_NO_MEMORY;
    int status = domain_actions(s, generator_count, generators, actions);
    if (!status) {
        orbit_rows.ints_used = orbit_rows.doubles_used = 0;
        int64_t *ints = row_ints(&orbit_rows, 2 + (size_t)d_count * (1 + generator_count) + s->subgroups_used);
        double *doubles = row_doubles(&orbit_rows, s->cells * d_count + 1);
        if (!ints || !doubles) {
            status = ISO_NO_MEMORY;
        } else {
            *ints++ = d_count;
            *ints++ = generator_count;
            for (int d = 0; d < d_count; d++) *ints++ = s->numbers[d];
            for (size_t i = 0; i < (size_t)generator_count * d_count; i++) *ints++ = actions[i];
            memcpy(ints, s->subgroups, sizeof(int64_t) * s->subgroups_used);
            memcpy(doubles, s->domains, sizeof(double) * s->cells * d_count);
            orbit_rows.doubles_used = s->cells * d_count;
        }
    }
    free(actions);
    return status;
}

/* One slot of coupled_opd_rows: the domain orbit of a source subspace
 * (values: free_rows x dim) and the action of the generating operations on
 * it.  mode: the point action of the records (mag_point_action); candidates
 * are the orderparam_to_subgroup_ records.  numbered_count < 0: no numbered
 * operations.  When the numbered order is used and some action leaves its
 * domains, the orbit falls back to the operations' order.
 * Results: coupled_domain_orbit_ints (domain count, generator count, display
 * numbers, actions [generator][domain], per domain basis, kept count and
 * kept records) and coupled_domain_orbit_doubles (the domains). */
void op_translation_cache_begin(void);
void op_translation_cache_end(void);

int coupled_domain_orbit(int gid, int param_count, const int64_t *kparam, int mode, int free_rows,
                         const double *values, int candidate_count, const int64_t *candidates, int op_count,
                         const int64_t *op_records, int numbered_count, const int64_t *numbered_records,
                         int use_numbered, int generator_count, const int64_t *generator_records)
{
    orbit_state s;
    memset(&s, 0, sizeof s);
    if (mode && mag_tables()) return DO_MAGNETIC;
    int status = irrep_context_for(gid, param_count, kparam, &s.irrep);
    if (status) return status == IR_RANGE ? 103 : status == IR_NO_TABLES ? 106 : 105;
    s.dim = s.irrep.dim;
    s.free_rows = free_rows;
    s.cells = (size_t)free_rows * s.dim;
    s.values = values;
    s.candidate_count = candidate_count;
    s.candidates = candidates;
    s.candidate_ops = malloc(sizeof(int) * ((size_t)candidate_count + 1));
    s.candidate_signs = malloc(sizeof(int) * ((size_t)candidate_count + 1));
    s.numbered_count = numbered_count;
    s.ctx.dim = s.dim;
    double *ops = 0, *numbered = 0, *generators = 0;
    status = context_init(&s.ctx);
    if (!status && (!s.candidate_ops || !s.candidate_signs)) status = ISO_NO_MEMORY;
    for (int i = 0; i < candidate_count && !status; i++)
        s.candidate_ops[i] = mag_point_action(mode, (int)candidates[5 * i + 4], s.candidate_signs + i);
    if (!status && mode && mag_take_error()) status = DO_MAGNETIC;
    if (!status) status = action_matrices(gid, param_count, kparam, mode, s.dim, op_count, op_records, &ops);
    if (!status && numbered_count > 0)
        status = action_matrices(gid, param_count, kparam, mode, s.dim, numbered_count, numbered_records, &numbered);
    if (!status)
        status = action_matrices(gid, param_count, kparam, mode, s.dim, generator_count, generator_records, &generators);
    if (!status) {
        int context = irrep_context_for(gid, param_count, kparam, &s.irrep);
        if (context) status = context == IR_RANGE ? 103 : context == IR_NO_TABLES ? 106 : 105;
    }
    s.numbered = numbered;
    s.stride = malloc(sizeof(double) * (48 * (size_t)free_rows + 1));
    s.kept = malloc(sizeof(int64_t) * (5 * (size_t)candidate_count + 1));
    if (!status && (!s.stride || !s.kept)) status = ISO_NO_MEMORY;
    int numbered_order = use_numbered && numbered_count >= 0;
    op_translation_cache_begin();
    if (!status) {
        status = numbered_order ? orbit_pass(&s, numbered_count, numbered, 0)
                                : orbit_pass(&s, op_count, ops, numbered_count >= 0);
        if (!status) status = orbit_emit(&s, generator_count, generators);
        if (status == DA_UNMATCHED && numbered_order) {
            status = orbit_pass(&s, op_count, ops, 1);
            if (!status) status = orbit_emit(&s, generator_count, generators);
        }
    }
    op_translation_cache_end();
    free(ops);
    free(numbered);
    free(generators);
    free(s.candidate_ops);
    free(s.candidate_signs);
    free(s.stride);
    free(s.kept);
    free(s.domain_of);
    free(s.first_number);
    free(s.domains);
    free(s.numbers);
    free(s.subgroup_at);
    free(s.subgroups);
    context_free(&s.ctx);
    return status;
}

int64_t coupled_domain_orbit_ints(const int64_t **data)
{
    *data = orbit_rows.ints;
    return (int64_t)orbit_rows.ints_used;
}

int64_t coupled_domain_orbit_doubles(const double **data)
{
    *data = orbit_rows.doubles;
    return (int64_t)orbit_rows.doubles_used;
}
