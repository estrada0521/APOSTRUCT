/* Parent atom rows presented in a child setting (presentation.py
 * present_mode_rows): transformed by the inverse basis, folded into the
 * cell, deduplicated up to the centering translations, then expanded by
 * them, with the same periodic position buckets and float semantics as the
 * original (Python % and sum). */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "pynum.h"

#include "pysum.h"

enum {
    PR_OK = 0,
    PR_SINGULAR = 1,
    PR_MEMORY = 2,
    PR_CAPACITY = 3,
};

static int periodic_close3(const double *left, const double *right, double tol)
{
    for (int i = 0; i < 3; i++)
        if (!(fabs(pymod1(left[i] - right[i] + 0.5) - 0.5) <= tol)) return 0;
    return 1;
}

static int same_vector(const double *left, const double *right, double tol)
{
    for (int i = 0; i < 3; i++)
        if (!(fabs(left[i] - right[i]) <= tol)) return 0;
    return 1;
}

typedef struct {
    int key[3];
    int head;
} bucket;

typedef struct {
    bucket *slots;
    size_t cap, used;
    int *next, *item;
    size_t entries, entry_cap;
} bucket_map;

static uint64_t key_hash(const int *key)
{
    uint64_t h = (uint64_t)(uint32_t)key[0] * 0x9e3779b97f4a7c15ull;
    h ^= (uint64_t)(uint32_t)key[1] + 0x7f4a7c159e3779b9ull + (h << 6) + (h >> 2);
    h ^= (uint64_t)(uint32_t)key[2] * 0xbf58476d1ce4e5b9ull;
    h ^= h >> 31;
    return h;
}

static int map_init(bucket_map *m)
{
    memset(m, 0, sizeof *m);
    m->cap = 256;
    m->slots = malloc(sizeof(bucket) * m->cap);
    m->entry_cap = 256;
    m->next = malloc(sizeof(int) * m->entry_cap);
    m->item = malloc(sizeof(int) * m->entry_cap);
    if (!m->slots || !m->next || !m->item) return PR_MEMORY;
    for (size_t i = 0; i < m->cap; i++) m->slots[i].head = -2;
    return PR_OK;
}

static void map_free(bucket_map *m)
{
    free(m->slots);
    free(m->next);
    free(m->item);
}

static bucket *map_find(const bucket_map *m, const int *key)
{
    size_t mask = m->cap - 1, i = key_hash(key) & mask;
    while (m->slots[i].head != -2) {
        if (!memcmp(m->slots[i].key, key, sizeof(int) * 3)) return &m->slots[i];
        i = (i + 1) & mask;
    }
    return 0;
}

static int map_add(bucket_map *m, const int *key, int value)
{
    if (2 * (m->used + 1) > m->cap) {
        bucket *old = m->slots;
        size_t old_cap = m->cap;
        m->cap *= 2;
        m->slots = malloc(sizeof(bucket) * m->cap);
        if (!m->slots) return PR_MEMORY;
        for (size_t i = 0; i < m->cap; i++) m->slots[i].head = -2;
        for (size_t j = 0; j < old_cap; j++) {
            if (old[j].head == -2) continue;
            size_t i = key_hash(old[j].key) & (m->cap - 1);
            while (m->slots[i].head != -2) i = (i + 1) & (m->cap - 1);
            m->slots[i] = old[j];
        }
        free(old);
    }
    bucket *b = map_find(m, key);
    if (!b) {
        size_t i = key_hash(key) & (m->cap - 1);
        while (m->slots[i].head != -2) i = (i + 1) & (m->cap - 1);
        b = &m->slots[i];
        memcpy(b->key, key, sizeof(int) * 3);
        b->head = -1;
        m->used++;
    }
    if (m->entries == m->entry_cap) {
        m->entry_cap *= 2;
        int *next = realloc(m->next, sizeof(int) * m->entry_cap);
        if (!next) return PR_MEMORY;
        m->next = next;
        int *item = realloc(m->item, sizeof(int) * m->entry_cap);
        if (!item) return PR_MEMORY;
        m->item = item;
    }
    m->item[m->entries] = value;
    m->next[m->entries] = b->head;
    b->head = (int)m->entries++;
    return PR_OK;
}

/* _periodic_bucket_key; 0 when a coordinate is not finite. */
static int bucket_key(const double *values, int count, int *key)
{
    for (int axis = 0; axis < 3; axis++) {
        double folded = pymod1(values[axis]);
        if (!isfinite(folded)) return 0;
        int k = (int)floor(folded * count);
        key[axis] = k < count - 1 ? k : count - 1;
    }
    return 1;
}

/* Visit every index in the buckets around center (with neighbours within
 * tol), each once; returns 1 as soon as visit returns 1. */
typedef int (*visitor)(void *context, int index);

static int visit_candidates(const bucket_map *m, const double *values, const int *center, int count, double tol,
                            unsigned char *mark, int *touched, visitor visit, void *context)
{
    int axis_keys[3][3], axis_count[3];
    for (int axis = 0; axis < 3; axis++) {
        double folded = pymod1(values[axis]);
        int key = center[axis];
        axis_keys[axis][0] = key;
        axis_count[axis] = 1;
        double lower = (double)key / count, upper = (double)(key + 1) / count;
        if (folded - lower <= tol) {
            int k = ((key - 1) % count + count) % count;
            int seen = 0;
            for (int i = 0; i < axis_count[axis]; i++) seen |= axis_keys[axis][i] == k;
            if (!seen) axis_keys[axis][axis_count[axis]++] = k;
        }
        if (upper - folded <= tol) {
            int k = (key + 1) % count;
            int seen = 0;
            for (int i = 0; i < axis_count[axis]; i++) seen |= axis_keys[axis][i] == k;
            if (!seen) axis_keys[axis][axis_count[axis]++] = k;
        }
    }
    int touched_count = 0, found = 0;
    for (int a = 0; a < axis_count[0] && !found; a++)
        for (int b = 0; b < axis_count[1] && !found; b++)
            for (int c = 0; c < axis_count[2] && !found; c++) {
                int key[3] = {axis_keys[0][a], axis_keys[1][b], axis_keys[2][c]};
                const bucket *bk = map_find(m, key);
                if (!bk) continue;
                for (int e = bk->head; e >= 0 && !found; e = m->next[e]) {
                    int index = m->item[e];
                    if (mark[index]) continue;
                    mark[index] = 1;
                    touched[touched_count++] = index;
                    if (visit(context, index)) found = 1;
                }
            }
    for (int i = 0; i < touched_count; i++) mark[touched[i]] = 0;
    return found;
}

typedef struct {
    const double *xyz, *dxyz;
    const double (*rows_xyz)[3];
    const double (*rows_dxyz)[3];
    const double (*translations)[3];
    int translation_count;
    double tol;
} quotient_probe;

static int quotient_matches(void *context, int index)
{
    const quotient_probe *p = context;
    if (!same_vector(p->dxyz, p->rows_dxyz[index], p->tol)) return 0;
    for (int t = 0; t < p->translation_count; t++) {
        double shifted[3];
        for (int i = 0; i < 3; i++) shifted[i] = p->rows_xyz[index][i] + p->translations[t][i];
        if (periodic_close3(p->xyz, shifted, p->tol)) return 1;
    }
    return 0;
}

static int conventional_matches(void *context, int index)
{
    const quotient_probe *p = context;
    return periodic_close3(p->xyz, p->rows_xyz[index], p->tol) && same_vector(p->dxyz, p->rows_dxyz[index], p->tol);
}

/* rows: count parent (xyz, dxyz); out: conventional rows as (source row,
 * centering ordinal) with their xyz and dxyz. */
int present_rows(int count, const double *parent_xyz, const double *parent_dxyz, const double *basis,
                 const double *origin, int centering_letter, double tol, int capacity, int *out_source,
                 int *out_centering, double *out_xyz, double *out_dxyz, int *out_count)
{
    const double *a = basis;
#define A(r, c) a[(r) * 3 + (c)]
    double det = A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) - A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
                 A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
    if (fabs(det) <= 1e-15) return PR_SINGULAR;
    double inv[3][3] = {
        {(A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) / det, (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2)) / det,
         (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) / det},
        {(A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)) / det, (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) / det,
         (A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2)) / det},
        {(A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) / det, (A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1)) / det,
         (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) / det},
    };
#undef A
    double translations[4][3] = {{0.0, 0.0, 0.0}};
    int translation_count = 1;
    switch (centering_letter) {
    case 'A': translations[1][0] = 0.0; translations[1][1] = 0.5; translations[1][2] = 0.5; translation_count = 2; break;
    case 'B': translations[1][0] = 0.5; translations[1][1] = 0.0; translations[1][2] = 0.5; translation_count = 2; break;
    case 'C': translations[1][0] = 0.5; translations[1][1] = 0.5; translations[1][2] = 0.0; translation_count = 2; break;
    case 'I': translations[1][0] = 0.5; translations[1][1] = 0.5; translations[1][2] = 0.5; translation_count = 2; break;
    case 'F':
        translations[1][0] = 0.0; translations[1][1] = 0.5; translations[1][2] = 0.5;
        translations[2][0] = 0.5; translations[2][1] = 0.0; translations[2][2] = 0.5;
        translations[3][0] = 0.5; translations[3][1] = 0.5; translations[3][2] = 0.0;
        translation_count = 4;
        break;
    case 'R':
        translations[1][0] = 2.0 / 3.0; translations[1][1] = 1.0 / 3.0; translations[1][2] = 1.0 / 3.0;
        translations[2][0] = 1.0 / 3.0; translations[2][1] = 2.0 / 3.0; translations[2][2] = 2.0 / 3.0;
        translation_count = 3;
        break;
    default: break;
    }
    int status = PR_OK;
    double (*xyz)[3] = malloc(sizeof(double[3]) * (size_t)(count + 1));
    double (*dxyz)[3] = malloc(sizeof(double[3]) * (size_t)(count + 1));
    int *quotient = malloc(sizeof(int) * (size_t)(count + 1));
    double (*q_xyz)[3] = malloc(sizeof(double[3]) * (size_t)(count + 1));
    double (*q_dxyz)[3] = malloc(sizeof(double[3]) * (size_t)(count + 1));
    unsigned char *mark = calloc((size_t)(count * 4 + 1), 1);
    int *touched = malloc(sizeof(int) * (size_t)(count * 4 + 1));
    bucket_map orbit_buckets, position_buckets;
    int maps = !map_init(&orbit_buckets) & !map_init(&position_buckets);
    if (!xyz || !dxyz || !quotient || !q_xyz || !q_dxyz || !mark || !touched || !maps) {
        status = PR_MEMORY;
        goto done;
    }
    double scale = 0.0;
    for (int r = 0; r < count; r++) {
        double shifted[3];
        for (int i = 0; i < 3; i++) shifted[i] = parent_xyz[3 * r + i] - origin[i];
        for (int axis = 0; axis < 3; axis++) {
            xyz[r][axis] = pymod1(apostruct_builtin_sum3(shifted[0] * inv[0][axis], shifted[1] * inv[1][axis],
                                                         shifted[2] * inv[2][axis]));
            dxyz[r][axis] = apostruct_builtin_sum3(parent_dxyz[3 * r] * inv[0][axis],
                                                   parent_dxyz[3 * r + 1] * inv[1][axis],
                                                   parent_dxyz[3 * r + 2] * inv[2][axis]);
        }
        for (int axis = 0; axis < 3; axis++)
            if (fabs(dxyz[r][axis]) > scale) scale = fabs(dxyz[r][axis]);
    }
    if (scale > 0.0)
        for (int r = 0; r < count; r++)
            for (int axis = 0; axis < 3; axis++) dxyz[r][axis] = dxyz[r][axis] / scale;

    double width = 0.0;
    int bucketed = isfinite(tol) && tol > 0.0;
    if (bucketed) {
        width = 2.0 * tol > sqrt(tol) ? 2.0 * tol : sqrt(tol);
        bucketed = width > 0.0 && isfinite(1.0 / width);
    }
    int periodic_count = bucketed ? (int)floor(1.0 / width) : 1;
    if (periodic_count < 1) periodic_count = 1;
    int quotient_count = 0, index_complete = 1;
    for (int r = 0; r < count && !status; r++) {
        quotient_probe probe = {xyz[r], dxyz[r], (const double (*)[3])q_xyz, (const double (*)[3])q_dxyz,
                                (const double (*)[3])translations, translation_count, tol};
        int key[3], found = 0;
        if (bucketed && index_complete && bucket_key(xyz[r], periodic_count, key)) {
            found = visit_candidates(&orbit_buckets, xyz[r], key, periodic_count, tol, mark, touched, quotient_matches,
                                     &probe);
        } else {
            for (int q = 0; q < quotient_count && !found; q++) found = quotient_matches(&probe, q);
        }
        if (found) continue;
        int q = quotient_count++;
        quotient[q] = r;
        memcpy(q_xyz[q], xyz[r], sizeof q_xyz[q]);
        memcpy(q_dxyz[q], dxyz[r], sizeof q_dxyz[q]);
        if (bucketed) {
            int keys[4][3], all = 1;
            for (int t = 0; t < translation_count; t++) {
                double shifted[3];
                for (int i = 0; i < 3; i++) shifted[i] = pymod1(xyz[r][i] + translations[t][i]);
                all &= bucket_key(shifted, periodic_count, keys[t]);
            }
            if (all) {
                for (int t = 0; t < translation_count && !status; t++) status = map_add(&orbit_buckets, keys[t], q);
            } else {
                index_complete = 0;
            }
        }
    }
    int conventional = 0;
    double (*c_xyz)[3] = malloc(sizeof(double[3]) * (size_t)(quotient_count * translation_count + 1));
    double (*c_dxyz)[3] = malloc(sizeof(double[3]) * (size_t)(quotient_count * translation_count + 1));
    if (!status && (!c_xyz || !c_dxyz)) status = PR_MEMORY;
    for (int q = 0; q < quotient_count && !status; q++)
        for (int t = 0; t < translation_count && !status; t++) {
            double candidate[3];
            for (int i = 0; i < 3; i++) candidate[i] = pymod1(q_xyz[q][i] + translations[t][i]);
            quotient_probe probe = {candidate, q_dxyz[q], (const double (*)[3])c_xyz, (const double (*)[3])c_dxyz,
                                    0, 0, tol};
            int key[3], found = 0, keyed = bucketed && bucket_key(candidate, periodic_count, key);
            if (keyed) {
                found = visit_candidates(&position_buckets, candidate, key, periodic_count, tol, mark, touched,
                                         conventional_matches, &probe);
            } else {
                for (int c = 0; c < conventional && !found; c++) found = conventional_matches(&probe, c);
            }
            if (found) continue;
            if (conventional >= capacity) {
                status = PR_CAPACITY;
                break;
            }
            memcpy(c_xyz[conventional], candidate, sizeof candidate);
            memcpy(c_dxyz[conventional], q_dxyz[q], sizeof c_dxyz[conventional]);
            out_source[conventional] = quotient[q];
            out_centering[conventional] = t;
            memcpy(out_xyz + 3 * conventional, candidate, sizeof candidate);
            memcpy(out_dxyz + 3 * conventional, q_dxyz[q], sizeof candidate);
            if (keyed) status = map_add(&position_buckets, key, conventional);
            conventional++;
        }
    free(c_xyz);
    free(c_dxyz);
    *out_count = conventional;
done:
    free(xyz);
    free(dxyz);
    free(quotient);
    free(q_xyz);
    free(q_dxyz);
    free(mark);
    free(touched);
    map_free(&orbit_buckets);
    map_free(&position_buckets);
    return status;
}
