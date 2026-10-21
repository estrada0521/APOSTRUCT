/* Assign presented atom rows to canonical child atoms by their affine
 * occurrence in the presentation cell: each operation record (with its
 * child centering carried into the parent PML frame) is keyed by its point
 * operation and its translation modulo the presentation lattice. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"

enum {
    RA_OK = 0,
    RA_RECORD = 1,
    RA_CENTERING = 2,
    RA_SINGULAR = 3,
    RA_SHARED = 4,
    RA_UNMATCHED = 5,
    RA_MEMORY = 6,
    RA_RANGE = 7,
};

typedef struct {
    int64_t key[5];
    int index;
} keyed;

static int64_t gcd64(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        int64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/* add_translation_to_operation_record */
static void add_translation(const int64_t *t, const int64_t *o, int64_t *out)
{
    int64_t den = t[3] * o[3];
    int64_t x = t[0] * o[3] + o[0] * t[3], y = t[1] * o[3] + o[1] * t[3], z = t[2] * o[3] + o[2] * t[3];
    int64_t common = gcd64(gcd64(x, y), gcd64(z, den));
    if (den < 0) common = -common;
    out[0] = x / common;
    out[1] = y / common;
    out[2] = z / common;
    out[3] = den / common;
    out[4] = o[4];
}

static int key_order(const void *a, const void *b)
{
    const keyed *x = a, *y = b;
    for (int i = 0; i < 5; i++)
        if (x->key[i] != y->key[i]) return x->key[i] < y->key[i] ? -1 : 1;
    return x->index - y->index;
}

typedef struct {
    int64_t det, adj[9];
} lattice;

static int occurrence_key(const lattice *l, const int64_t *record, int64_t *key)
{
    i128 den = (i128)record[3] * l->det;
    if (den > INT64_MAX) return RA_RANGE;
    int64_t c[3];
    for (int col = 0; col < 3; col++) {
        i128 value = 0;
        for (int row = 0; row < 3; row++) value += (i128)record[row] * l->adj[3 * row + col];
        value %= den;
        if (value < 0) value += den;
        c[col] = (int64_t)value;
    }
    int64_t common = gcd64(gcd64(c[0], c[1]), gcd64(c[2], (int64_t)den));
    key[0] = record[4];
    for (int col = 0; col < 3; col++) key[1 + col] = c[col] / common;
    key[4] = (int64_t)den / common;
    return 0;
}

/* atoms: atom_count x (kernel[4], coset[5], centering); rows: row_count x
 * (record[5], centering or -1).  out[row] = canonical atom index. */
int row_assignments(int child_sg, const int64_t *basis, int point_operation_count, int atom_count,
                    const int64_t *atoms, int row_count, const int64_t *rows, int *out)
{
    if (geo_load()) return RA_CENTERING;
    rat centering[4][3];
    int centering_count = geo_centering(child_sg, centering);
    int64_t centering_records[4][4];
    for (int c = 0; c < centering_count; c++) {
        rat pml[3];
        for (int col = 0; col < 3; col++) {
            rat acc = rat_int(0);
            for (int row = 0; row < 3; row++) acc = rat_add(acc, rat_mul(centering[c][row], rat_int(basis[3 * row + col])));
            pml[col] = acc;
        }
        i128 den = 1;
        for (int col = 0; col < 3; col++) den = den / rat_gcd(den, pml[col].d) * pml[col].d;
        for (int col = 0; col < 3; col++) centering_records[c][col] = (int64_t)(pml[col].n * (den / pml[col].d));
        centering_records[c][3] = (int64_t)den;
    }
    lattice l;
    const int64_t *b = basis;
    int64_t adj[9] = {
        b[4] * b[8] - b[5] * b[7], b[2] * b[7] - b[1] * b[8], b[1] * b[5] - b[2] * b[4],
        b[5] * b[6] - b[3] * b[8], b[0] * b[8] - b[2] * b[6], b[2] * b[3] - b[0] * b[5],
        b[3] * b[7] - b[4] * b[6], b[1] * b[6] - b[0] * b[7], b[0] * b[4] - b[1] * b[3],
    };
    l.det = b[0] * adj[0] + b[1] * adj[3] + b[2] * adj[6];
    if (!l.det) return RA_SINGULAR;
    for (int i = 0; i < 9; i++) l.adj[i] = l.det < 0 ? -adj[i] : adj[i];
    if (l.det < 0) l.det = -l.det;

    keyed *canonical = malloc(sizeof *canonical * (size_t)(atom_count + 1));
    if (!canonical) return RA_MEMORY;
    int status = RA_OK;
    for (int a = 0; a < atom_count && !status; a++) {
        const int64_t *atom = atoms + 10 * a;
        int64_t centering_ordinal = atom[9];
        if (centering_ordinal < 0 || centering_ordinal >= centering_count) {
            status = RA_CENTERING;
            break;
        }
        if (!atom[3] || !atom[7]) {
            status = RA_RECORD;
            break;
        }
        int64_t source[5], record[5];
        add_translation(atom, atom + 4, source);
        add_translation(centering_records[centering_ordinal], source, record);
        if (record[3] <= 0 || record[4] < 1 || record[4] > point_operation_count) {
            status = RA_RECORD;
            break;
        }
        status = occurrence_key(&l, record, canonical[a].key);
        canonical[a].index = a;
    }
    if (!status) qsort(canonical, (size_t)atom_count, sizeof *canonical, key_order);
    for (int a = 1; a < atom_count && !status; a++)
        if (!memcmp(canonical[a - 1].key, canonical[a].key, sizeof canonical[a].key)) status = RA_SHARED;
    for (int r = 0; r < row_count && !status; r++) {
        const int64_t *row = rows + 6 * r;
        int64_t record[5];
        memcpy(record, row, sizeof record);
        if (record[3] <= 0 || record[4] < 1 || record[4] > point_operation_count) {
            status = RA_RECORD;
            break;
        }
        if (row[5] >= 0) {
            if (row[5] >= centering_count) {
                status = RA_CENTERING;
                break;
            }
            int64_t shifted[5];
            add_translation(centering_records[row[5]], record, shifted);
            memcpy(record, shifted, sizeof shifted);
            if (record[3] <= 0) {
                status = RA_RECORD;
                break;
            }
        }
        keyed probe;
        status = occurrence_key(&l, record, probe.key);
        if (status) break;
        int lo = 0, hi = atom_count - 1, found = -1;
        while (lo <= hi && found < 0) {
            int mid = (lo + hi) / 2;
            int c = memcmp(canonical[mid].key, probe.key, sizeof probe.key) ? 0 : 1;
            if (c) {
                found = canonical[mid].index;
                break;
            }
            probe.index = -1;
            if (key_order(&canonical[mid], &probe) < 0)
                lo = mid + 1;
            else
                hi = mid - 1;
        }
        if (found < 0) {
            status = RA_UNMATCHED;
            break;
        }
        out[r] = found;
    }
    free(canonical);
    return status;
}
