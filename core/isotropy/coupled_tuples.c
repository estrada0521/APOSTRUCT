/* The per-domain-tuple part of coupled_opd_rows for one source combination:
 * for every relative domain tuple in order, the combined stabilizer (newlat
 * chain over the slots' domain rows, then the intersection of the slots'
 * domain subgroups in the retained stored order) and the per-slot subduction
 * frequencies over the slot kernels and the common kernel.  Tuples whose
 * common-kernel frequencies equal the expected ones are kept for the caller;
 * the first failure stops the scan and is reported after the kept tuples. */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "intmath.h"

int orderparam_newlat_search(int gid, int param_count, const int64_t *kparam, const double *row, int has_contains,
                             const int64_t *contains, int64_t contains_den, int max_iter, int64_t *rows_out);
int subgroup_intersection(int sg, int group, const int64_t *basis, int subgroup_count, const int64_t *bases,
                          const int *operation_counts, const int64_t *operations, int64_t capacity, int64_t *out);
int subduction_frequency(int gid, int param_count, const int64_t *kparam, int mode, const int64_t *kernel_basis,
                         const int64_t *subgroup_basis, int op_count, const int64_t *records, int *frequency);

enum {
    CT_OK = 0,
    CT_NEWLAT = 1,     /* status: the orderparam error code */
    CT_SINGULAR = 2,   /* the previous chain basis is singular */
    CT_SUBGROUP = 3,   /* status: the subgroup operation error code */
    CT_FREQUENCY = 4,  /* status: the subduction frequency return code */
    CT_MEMORY = 5,
    CT_RANGE = 6,
};

enum { CT_MAX_SLOTS = 16, CT_MAX_OPS = 192 };


static int64_t *out_ints;
static size_t out_used, out_cap;
static int64_t failure[3];

static int64_t *grow(size_t count)
{
    if (out_used + count > out_cap) {
        size_t next = out_cap ? out_cap : 4096;
        while (next < out_used + count) next *= 2;
        int64_t *grown = realloc(out_ints, sizeof(int64_t) * next);
        if (!grown) return 0;
        out_ints = grown;
        out_cap = next;
    }
    int64_t *at = out_ints + out_used;
    out_used += count;
    return at;
}

static void reduce_rows(int64_t *rows)
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

/* SourceData._search_newlat_basis's contains argument: the exact inverse of
 * basis scaled by the least common denominator of its reduced entries. */
static int inverse_contains(const int64_t *basis, int64_t *contains, int64_t *scale_out)
{
    i128 b[9];
    for (int i = 0; i < 9; i++) b[i] = basis[i];
    i128 adj[9] = {
        b[4] * b[8] - b[5] * b[7], b[2] * b[7] - b[1] * b[8], b[1] * b[5] - b[2] * b[4],
        b[5] * b[6] - b[3] * b[8], b[0] * b[8] - b[2] * b[6], b[2] * b[3] - b[0] * b[5],
        b[3] * b[7] - b[4] * b[6], b[1] * b[6] - b[0] * b[7], b[0] * b[4] - b[1] * b[3],
    };
    i128 det = b[0] * adj[0] + b[1] * adj[3] + b[2] * adj[6];
    if (!det) return CT_SINGULAR;
    i128 num[9], den[9], scale = 1;
    for (int i = 0; i < 9; i++) {
        i128 n = adj[i], d = det;
        if (d < 0) {
            n = -n;
            d = -d;
        }
        i128 g = gcd_i128(n, d);
        if (g < 0) return CT_RANGE;
        if (g > 1) {
            n /= g;
            d /= g;
        }
        if (!n) d = 1;
        num[i] = n;
        den[i] = d;
        i128 h = gcd_i128(scale, d);
        if (h <= 0) return CT_RANGE;
        scale = scale / h * d;
        if (scale > INT64_MAX) return CT_RANGE;
    }
    for (int i = 0; i < 9; i++) {
        i128 value = num[i] * (scale / den[i]);
        if (value > INT64_MAX || value < INT64_MIN) return CT_RANGE;
        contains[i] = (int64_t)value;
    }
    *scale_out = (int64_t)scale;
    return CT_OK;
}

/* Open-addressing maps keyed by int64 words: fixed width (prefix chains)
 * or variable length (frequency sets) with keys kept in an arena. */
typedef struct {
    size_t cap, used;
    uint64_t *hashes;
    int64_t *key_at, *key_len, *value_at;
    unsigned char *full;
    int64_t *arena;
    size_t arena_used, arena_cap;
} map64;

static int map_init(map64 *m)
{
    memset(m, 0, sizeof *m);
    m->cap = 256;
    m->hashes = malloc(sizeof(uint64_t) * m->cap);
    m->key_at = malloc(sizeof(int64_t) * m->cap);
    m->key_len = malloc(sizeof(int64_t) * m->cap);
    m->value_at = malloc(sizeof(int64_t) * m->cap);
    m->full = calloc(m->cap, 1);
    return m->hashes && m->key_at && m->key_len && m->value_at && m->full ? 0 : CT_MEMORY;
}

static void map_free(map64 *m)
{
    free(m->hashes);
    free(m->key_at);
    free(m->key_len);
    free(m->value_at);
    free(m->full);
    free(m->arena);
    memset(m, 0, sizeof *m);
}

static uint64_t hash_words(const int64_t *v, size_t n)
{
    uint64_t h = 1469598103934665603ull ^ n;
    for (size_t i = 0; i < n; i++) {
        h ^= (uint64_t)v[i];
        h *= 1099511628211ull;
        h ^= h >> 29;
    }
    return h;
}

static int64_t *map_find(map64 *m, const int64_t *key, size_t length)
{
    uint64_t h = hash_words(key, length);
    size_t mask = m->cap - 1, i = h & mask;
    while (m->full[i]) {
        if (m->hashes[i] == h && (size_t)m->key_len[i] == length &&
            !memcmp(m->arena + m->key_at[i], key, sizeof(int64_t) * length))
            return m->arena + m->value_at[i];
        i = (i + 1) & mask;
    }
    return 0;
}

static int arena_reserve(map64 *m, size_t count)
{
    if (m->arena_used + count <= m->arena_cap) return 0;
    size_t next = m->arena_cap ? m->arena_cap : 4096;
    while (next < m->arena_used + count) next *= 2;
    int64_t *grown = realloc(m->arena, sizeof(int64_t) * next);
    if (!grown) return CT_MEMORY;
    m->arena = grown;
    m->arena_cap = next;
    return 0;
}

static int map_grow(map64 *m)
{
    size_t cap = m->cap * 2, mask = cap - 1;
    uint64_t *hashes = malloc(sizeof(uint64_t) * cap);
    int64_t *key_at = malloc(sizeof(int64_t) * cap), *key_len = malloc(sizeof(int64_t) * cap);
    int64_t *value_at = malloc(sizeof(int64_t) * cap);
    unsigned char *full = calloc(cap, 1);
    if (!hashes || !key_at || !key_len || !value_at || !full) {
        free(hashes);
        free(key_at);
        free(key_len);
        free(value_at);
        free(full);
        return CT_MEMORY;
    }
    for (size_t i = 0; i < m->cap; i++) {
        if (!m->full[i]) continue;
        size_t j = m->hashes[i] & mask;
        while (full[j]) j = (j + 1) & mask;
        full[j] = 1;
        hashes[j] = m->hashes[i];
        key_at[j] = m->key_at[i];
        key_len[j] = m->key_len[i];
        value_at[j] = m->value_at[i];
    }
    free(m->hashes);
    free(m->key_at);
    free(m->key_len);
    free(m->value_at);
    free(m->full);
    m->hashes = hashes;
    m->key_at = key_at;
    m->key_len = key_len;
    m->value_at = value_at;
    m->full = full;
    m->cap = cap;
    return 0;
}

/* Inserts key (not present) with value_length words; returns the value. */
static int64_t *map_insert(map64 *m, const int64_t *key, size_t length, size_t value_length)
{
    if (2 * (m->used + 1) > m->cap && map_grow(m)) return 0;
    if (arena_reserve(m, length + value_length)) return 0;
    uint64_t h = hash_words(key, length);
    size_t mask = m->cap - 1, i = h & mask;
    while (m->full[i]) i = (i + 1) & mask;
    m->full[i] = 1;
    m->hashes[i] = h;
    m->key_at[i] = (int64_t)m->arena_used;
    m->key_len[i] = (int64_t)length;
    memcpy(m->arena + m->arena_used, key, sizeof(int64_t) * length);
    m->arena_used += length;
    m->value_at[i] = (int64_t)m->arena_used;
    m->arena_used += value_length;
    m->used++;
    return m->arena + m->value_at[i];
}

/* The subduction frequencies keyed by (basis, records) hold every slot's
 * values; the slots and kernels stay fixed between coupled_tuples_reset
 * calls (one coupled_opd_rows request), so the cache spans its scans. */
static map64 frequencies, steps;
static int frequencies_ready, steps_ready;

void coupled_tuples_reset(void)
{
    if (frequencies_ready) map_free(&frequencies);
    if (steps_ready) map_free(&steps);
    frequencies_ready = steps_ready = 0;
}

/* One newlat step (SourceData.newlat_order / newlat_order3 and the row
 * reduction), cached for the request by (gid, kparam, previous basis, row)
 * with -0.0 read as 0.0, as Python's float equality does. */
static int newlat_step(const int64_t *slot, const double *row, int has_previous, const int64_t *previous,
                       int64_t *basis_out, int *kind_out, int *code_out)
{
    int dim = (int)slot[8];
    int64_t key[2 + 6 + 1 + 9 + 48];
    size_t length = 0;
    key[length++] = slot[0];
    key[length++] = slot[1];
    for (int i = 0; i < 6; i++) key[length++] = slot[2 + i];
    key[length++] = has_previous;
    for (int i = 0; i < 9; i++) key[length++] = has_previous ? previous[i] : 0;
    for (int i = 0; i < dim && i < 48; i++) {
        double value = row[i] == 0.0 ? 0.0 : row[i];
        int64_t bits;
        memcpy(&bits, &value, sizeof bits);
        key[length++] = bits;
    }
    if (!steps_ready) {
        if (map_init(&steps)) return CT_MEMORY;
        steps_ready = 1;
    }
    int64_t *known = map_find(&steps, key, length);
    if (known) {
        *kind_out = (int)known[0];
        *code_out = (int)known[1];
        memcpy(basis_out, known + 2, sizeof(int64_t) * 9);
        return CT_OK;
    }
    int64_t contains[9] = {0}, scale = 1, found[9] = {0};
    int kind = CT_OK, code = 0;
    if (has_previous) kind = inverse_contains(previous, contains, &scale);
    if (!kind) {
        code = orderparam_newlat_search((int)slot[0], (int)slot[1], slot + 2, row, has_previous, contains, scale, 256,
                                        found);
        if (code) kind = CT_NEWLAT;
        else reduce_rows(found);
    }
    int64_t *stored = map_insert(&steps, key, length, 11);
    if (!stored) return CT_MEMORY;
    stored[0] = kind;
    stored[1] = code;
    memcpy(stored + 2, found, sizeof found);
    *kind_out = kind;
    *code_out = code;
    memcpy(basis_out, found, sizeof found);
    return CT_OK;
}

/* newlat_step for the selection stage: slot words gid, param count,
 * kparam[6], full dim; returns a CT_ kind with *code. */
int coupled_newlat_step(const int64_t *slot, const double *row, int has_previous, const int64_t *previous, int64_t *basis_out,
                        int *code_out)
{
    int kind = CT_OK;
    int status = newlat_step(slot, row, has_previous, previous, basis_out, &kind, code_out);
    return status ? status : kind;
}

/* slots (SLOT_WIDTH words each): gid, param count, kparam[6], full dim,
 * free rows (the expected common-kernel frequency), point mode, slot kernel
 * basis[9], domain count, stored order length (-1 none).  values: per slot
 * per domain free x 48 doubles; bases / op_counts / ops: per slot per domain
 * the subgroup (basis[9], records x 5); orders: the stored orders' point
 * operations in slot order; tuples: tuple_count x slot_count domains.
 * Results (coupled_tuples_ints): kept count, then per kept tuple its index,
 * op count, basis[9], records, binary and coupled frequencies per slot.
 * Failure (coupled_tuples_failure): tuple index (-1 none), kind, status. */
enum { SLOT_WIDTH = 22 };

int coupled_tuples(int sg, int group, int slot_count, const int64_t *slots, const int64_t *kernel_basis,
                   const double *values, const int64_t *bases, const int *op_counts, const int64_t *ops,
                   const int64_t *orders, int magnetic, int mag_count, const int *mag2nonmag, int tuple_count,
                   const int32_t *tuples)
{
    out_used = 0;
    failure[0] = -1;
    failure[1] = failure[2] = 0;
    if (slot_count <= 0 || slot_count > CT_MAX_SLOTS) return CT_RANGE;
    const double *slot_values[CT_MAX_SLOTS];
    const int64_t *slot_bases[CT_MAX_SLOTS], *slot_ops[CT_MAX_SLOTS], *slot_order[CT_MAX_SLOTS];
    const int *slot_op_counts[CT_MAX_SLOTS];
    int64_t *op_offsets[CT_MAX_SLOTS] = {0};
    const double *v = values;
    const int64_t *b = bases, *o = ops, *ord = orders;
    const int *c = op_counts;
    int status = CT_OK;
    for (int s = 0; s < slot_count && !status; s++) {
        const int64_t *slot = slots + SLOT_WIDTH * s;
        int domains = (int)slot[20], free_rows = (int)slot[9];
        slot_values[s] = v;
        v += (size_t)domains * free_rows * 48;
        slot_bases[s] = b;
        b += 9 * (size_t)domains;
        slot_op_counts[s] = c;
        slot_ops[s] = o;
        op_offsets[s] = malloc(sizeof(int64_t) * ((size_t)domains + 1));
        if (!op_offsets[s]) {
            status = CT_MEMORY;
            break;
        }
        int64_t offset = 0;
        for (int d = 0; d < domains; d++) {
            op_offsets[s][d] = offset;
            offset += c[d];
        }
        o += 5 * offset;
        c += domains;
        slot_order[s] = slot[21] >= 0 ? ord : 0;
        if (slot[21] > 0) ord += slot[21];
    }
    map64 prefixes;
    int prefix_ready = !status && !map_init(&prefixes);
    if (!frequencies_ready && !status) {
        if (map_init(&frequencies)) status = CT_MEMORY;
        else frequencies_ready = 1;
    }
    int frequency_ready = frequencies_ready;
    int64_t *subgroup_bases = malloc(sizeof(int64_t) * 9 * (size_t)slot_count);
    int *subgroup_counts = malloc(sizeof(int) * (size_t)slot_count);
    int64_t *subgroup_ops = 0, *head = grow(1);
    size_t subgroup_ops_cap = 0;
    int64_t candidates[CT_MAX_OPS * 5], key[12 + 5 * CT_MAX_OPS];
    if (!prefix_ready || !frequency_ready || !subgroup_bases || !subgroup_counts || !head) status = CT_MEMORY;
    if (!status) *head = 0;
    for (int t = 0; t < tuple_count && !status; t++) {
        const int32_t *domain = tuples + (size_t)t * slot_count;
        /* newlat chain over the leading domains */
        int64_t chain[9];
        int has_chain = 0;
        for (int s = 0; s < slot_count && !status; s++) {
            int64_t prefix[CT_MAX_SLOTS + 1];
            prefix[0] = s;
            for (int i = 0; i <= s; i++) prefix[i + 1] = domain[i];
            int64_t *known = map_find(&prefixes, prefix, (size_t)s + 2);
            if (known) {
                has_chain = (int)known[0];
                memcpy(chain, known + 1, sizeof chain);
                continue;
            }
            const int64_t *slot = slots + SLOT_WIDTH * s;
            int dim = (int)slot[8], free_rows = (int)slot[9];
            const double *rows = slot_values[s] + (size_t)domain[s] * free_rows * 48;
            for (int r = 0; r < free_rows && !status; r++) {
                double row[48] = {0};
                for (int i = 0; i < dim && i < 48; i++) row[i] = rows[r * 48 + i];
                int64_t found[9];
                int kind, code;
                if (newlat_step(slot, row, has_chain, chain, found, &kind, &code)) {
                    status = CT_MEMORY;
                    break;
                }
                if (kind) {
                    failure[0] = t;
                    failure[1] = kind;
                    failure[2] = code;
                    status = kind;
                    break;
                }
                memcpy(chain, found, sizeof chain);
                has_chain = 1;
            }
            if (status) break;
            int64_t *stored = map_insert(&prefixes, prefix, (size_t)s + 2, 10);
            if (!stored) {
                status = CT_MEMORY;
                break;
            }
            stored[0] = has_chain;
            memcpy(stored + 1, chain, sizeof chain);
        }
        if (status) break;
        int64_t final_basis[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        if (has_chain) memcpy(final_basis, chain, sizeof final_basis);
        /* intersection of the slots' domain subgroups */
        size_t total_ops = 0;
        for (int s = 0; s < slot_count; s++) total_ops += (size_t)slot_op_counts[s][domain[s]];
        if (total_ops + 1 > subgroup_ops_cap) {
            int64_t *grown = realloc(subgroup_ops, sizeof(int64_t) * 5 * (total_ops + 1));
            if (!grown) {
                status = CT_MEMORY;
                break;
            }
            subgroup_ops = grown;
            subgroup_ops_cap = total_ops + 1;
        }
        size_t at = 0;
        for (int s = 0; s < slot_count; s++) {
            int d = domain[s], count = slot_op_counts[s][d];
            memcpy(subgroup_bases + 9 * s, slot_bases[s] + 9 * (size_t)d, sizeof(int64_t) * 9);
            subgroup_counts[s] = count;
            memcpy(subgroup_ops + 5 * at, slot_ops[s] + 5 * (size_t)op_offsets[s][d], sizeof(int64_t) * 5 * (size_t)count);
            at += (size_t)count;
        }
        int candidate_count = subgroup_intersection(sg, group, final_basis, slot_count, subgroup_bases, subgroup_counts,
                                                    subgroup_ops, CT_MAX_OPS, candidates);
        if (candidate_count < 0) {
            failure[0] = t;
            failure[1] = CT_SUBGROUP;
            failure[2] = -candidate_count;
            status = CT_SUBGROUP;
            break;
        }
        /* the retained stored order: magnetic keeps the final slot's,
         * otherwise the last nonempty one; a slot's order applies only to
         * its domain 0 */
        const int64_t *retained = 0;
        int retained_count = 0;
        for (int s = slot_count - 1; s >= 0 && !retained; s--) {
            const int64_t *slot = slots + SLOT_WIDTH * s;
            if (domain[s] == 0 && slot[21] > 0) {
                retained = slot_order[s];
                retained_count = (int)slot[21];
            }
            if (magnetic) break;
        }
        if (retained) {
            int position[256], distinct = 0;
            for (int i = 0; i < 256; i++) position[i] = -1;
            for (int i = 0; i < retained_count && !status; i++) {
                int64_t op = retained[i];
                if (op < 0 || op > 255) status = CT_RANGE;
                else {
                    if (position[op] < 0) distinct++;
                    position[op] = i;
                }
            }
            if (status) break;
            int order_keys[CT_MAX_OPS];
            for (int i = 0; i < candidate_count; i++) {
                int64_t op = candidates[5 * i + 4];
                if (magnetic) op = op >= 1 && op <= mag_count ? mag2nonmag[op - 1] : -1;
                order_keys[i] = op >= 0 && op < 256 && position[op] >= 0 ? position[op] : distinct;
            }
            for (int i = 1; i < candidate_count; i++) {
                int k = order_keys[i];
                int64_t record[5];
                memcpy(record, candidates + 5 * i, sizeof record);
                int j = i - 1;
                while (j >= 0 && order_keys[j] > k) {
                    order_keys[j + 1] = order_keys[j];
                    memcpy(candidates + 5 * (j + 1), candidates + 5 * j, sizeof record);
                    j--;
                }
                order_keys[j + 1] = k;
                memcpy(candidates + 5 * (j + 1), record, sizeof record);
            }
        }
        /* subduction frequencies over the slot kernels, then the common
         * kernel; one cache entry per (basis, records) holds every slot's */
        size_t key_length = 10 + 5 * (size_t)candidate_count;
        key[0] = candidate_count;
        memcpy(key + 1, final_basis, sizeof final_basis);
        memcpy(key + 10, candidates, sizeof(int64_t) * 5 * (size_t)candidate_count);
        int64_t *entry = map_find(&frequencies, key, key_length);
        if (!entry) {
            entry = map_insert(&frequencies, key, key_length, 4 * (size_t)slot_count);
            if (!entry) {
                status = CT_MEMORY;
                break;
            }
            for (int i = 0; i < 4 * slot_count; i += 2) entry[i] = -1;
        }
        int binary[CT_MAX_SLOTS], coupled[CT_MAX_SLOTS], matches = 1;
        for (int pass = 0; pass < 2 && !status; pass++)
            for (int s = 0; s < slot_count && !status; s++) {
                const int64_t *slot = slots + SLOT_WIDTH * s;
                int64_t *cell = entry + 2 * (pass * slot_count + s);
                if (cell[0] < 0) {
                    int frequency = 0;
                    int code = subduction_frequency((int)slot[0], (int)slot[1], slot + 2, (int)slot[10],
                                                    pass ? kernel_basis : slot + 11, final_basis, candidate_count,
                                                    candidates, &frequency);
                    cell[0] = code;
                    cell[1] = frequency;
                }
                if (cell[0]) {
                    failure[0] = t;
                    failure[1] = CT_FREQUENCY;
                    failure[2] = cell[0];
                    status = CT_FREQUENCY;
                    break;
                }
                (pass ? coupled : binary)[s] = (int)cell[1];
                if (pass && cell[1] != slot[9]) matches = 0;
            }
        if (status || !matches) continue;
        int64_t *kept = grow(2 + 9 + 5 * (size_t)candidate_count + 2 * (size_t)slot_count);
        if (!kept) {
            status = CT_MEMORY;
            break;
        }
        out_ints[0]++;
        *kept++ = t;
        *kept++ = candidate_count;
        memcpy(kept, final_basis, sizeof final_basis);
        kept += 9;
        memcpy(kept, candidates, sizeof(int64_t) * 5 * (size_t)candidate_count);
        kept += 5 * candidate_count;
        for (int s = 0; s < slot_count; s++) *kept++ = binary[s];
        for (int s = 0; s < slot_count; s++) *kept++ = coupled[s];
    }
    if (prefix_ready) map_free(&prefixes);
    free(subgroup_bases);
    free(subgroup_counts);
    free(subgroup_ops);
    for (int s = 0; s < slot_count; s++) free(op_offsets[s]);
    if (failure[0] >= 0) return CT_OK;
    return status;
}

int64_t coupled_tuples_ints(const int64_t **data)
{
    *data = out_ints;
    return (int64_t)out_used;
}

void coupled_tuples_failure(int64_t *out)
{
    out[0] = failure[0];
    out[1] = failure[1];
    out[2] = failure[2];
}
