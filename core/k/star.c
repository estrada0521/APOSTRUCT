/* The star of k of every (space group, k slot) with irreps, built from the k point's formula (little_k) and the parent's
 * point operations, in place of BYU's star tables: the arms are the images of the first arm in Source operation order,
 * equal arms (constant part mod 1, and up to sign when -k is in the star) kept once, constants reduced to [0,1).
 * BYU's own choices among equivalent arms are kept as conventions in little_k_star_flip (arms stored with the opposite
 * sign) and little_k_star_order (a different arm order); the conventional-to-ML parameter matrices are little_k_conv2ml.
 * section_find serves the star sections from here. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "rational.h"
#include "tables.h"

int section_find_stored(const char *name, section_view *view);

enum { SLOTS = 230 * 27, MAX_ARMS = 48 };

typedef rat arm_t[4][3];

static struct {
    int ready, failed;
    int64_t *star, *special, pointer[SLOTS], count[SLOTS], minusk[SLOTS], conv2ml[SLOTS];
    int64_t nstar, nspecial;
    blob_section entries[6];
} K;

/* an arm in integers over one denominator: rows constant and three parameter directions */
typedef int64_t iarm_t[4][3];

static void transform(const iarm_t v, const int64_t (*m)[3], iarm_t out)
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 3; c++) out[r][c] = v[r][0] * m[0][c] + v[r][1] * m[1][c] + v[r][2] * m[2][c];
}

/* equal arms over denominator den: constants equal mod 1, directions equal */
static int same_arm(const iarm_t a, const iarm_t b, int sign, int64_t den)
{
    for (int c = 0; c < 3; c++)
        if ((sign * a[0][c] - b[0][c]) % den) return 0;
    for (int r = 1; r < 4; r++)
        for (int c = 0; c < 3; c++)
            if (sign * a[r][c] != b[r][c]) return 0;
    return 1;
}

/* the reciprocal action of point operation op of sg on k row vectors: transpose(inverse(R)), R's rows the images of the
 * unit vectors; integral */
static int reciprocal(int sg, int op, int64_t (*out)[3])
{
    rat r[3][3], inverse[3][3];
    for (int i = 0; i < 3; i++) {
        rat e[3] = {rat_int(i == 0), rat_int(i == 1), rat_int(i == 2)};
        geo_vrot(sg, op, e, r[i]);
    }
    if (rat_inverse3((const rat (*)[3])r, inverse)) return 1;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            if (inverse[j][i].d != 1) return 1;
            out[i][j] = (int64_t)inverse[j][i].n;
        }
    return 0;
}

/* the generated arms over the denominator *den_out, before BYU's conventions; *minus whether -k is in the star */
static int generated(int sg, int64_t dim, const int64_t *first, iarm_t *kept, int *count, int *minus, int64_t *den_out)
{
    int64_t den = 1;
    for (int r = 0; r < 4; r++) {
        int64_t d = (dim > 0 || r == 0) ? first[4 * r + 3] : 0;
        if (d > 0) den = den / (int64_t)rat_gcd(den, d) * d;
    }
    iarm_t seed;
    for (int r = 0; r < 4; r++) {
        int64_t d = (dim > 0 || r == 0) ? first[4 * r + 3] : 0;
        for (int c = 0; c < 3; c++) seed[r][c] = d > 0 ? first[4 * r + c] * (den / d) : 0;
    }
    static int cached_sg, cached_n;
    static int64_t cached[GEO_MAX_OPS][3][3];
    if (sg != cached_sg) {
        int64_t records[GEO_MAX_OPS][5];
        cached_sg = 0;
        cached_n = geo_sg_records(sg, records);
        if (cached_n <= 0) return 1;
        for (int i = 0; i < cached_n; i++)
            if (reciprocal(sg, (int)records[i][4], cached[i])) return 1;
        cached_sg = sg;
    }
    int n = cached_n;
    iarm_t images[GEO_MAX_OPS];
    for (int i = 0; i < n; i++) transform((const int64_t (*)[3])seed, (const int64_t (*)[3])cached[i], images[i]);
    *minus = 0;
    if (!same_arm((const int64_t (*)[3])seed, (const int64_t (*)[3])seed, -1, den))
        for (int i = 0; i < n && !*minus; i++) *minus = same_arm((const int64_t (*)[3])images[i], (const int64_t (*)[3])seed, -1, den);
    *count = 0;
    for (int i = 0; i < n; i++) {
        int known = 0;
        for (int j = 0; j < *count && !known; j++)
            known = same_arm((const int64_t (*)[3])images[i], (const int64_t (*)[3])kept[j], 1, den) ||
                    (*minus && same_arm((const int64_t (*)[3])images[i], (const int64_t (*)[3])kept[j], -1, den));
        if (!known) {
            if (*count >= MAX_ARMS) return 1;
            memcpy(kept[(*count)++], images[i], sizeof(iarm_t));
        }
    }
    *den_out = den;
    return 0;
}

/* generated() as rationals, for the one-off convention generator */
int star_rule_generated(int sg, int kslot, int64_t dim, const int64_t *first, arm_t *out, int *count, int *minus)
{
    (void)kslot;
    iarm_t kept[MAX_ARMS];
    int64_t den;
    if (generated(sg, dim, first, kept, count, minus, &den)) return 1;
    for (int a = 0; a < *count; a++)
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 3; c++) out[a][r][c] = rat_make(kept[a][r][c], den);
    return 0;
}

/* a row over den as reduced numerators and their denominator */
static void encode_row(const int64_t *row, int64_t den, int64_t *w)
{
    int64_t g = den;
    for (int c = 0; c < 3; c++) g = (int64_t)rat_gcd(g, row[c] < 0 ? -row[c] : row[c]);
    for (int c = 0; c < 3; c++) w[c] = row[c] / g;
    w[3] = den / g;
}

static int build(void)
{
    section_view little_k, k_dim, irr_sg, irr_k, flips, orders, matrices, matrix_index;
    if (geo_load() || section_find_stored("k.slots.vector", &little_k) || section_find_stored("k.slots.dim", &k_dim) ||
        section_find_stored("irreps.little.sg", &irr_sg) || section_find_stored("irreps.little.k", &irr_k) ||
        section_find_stored("k.star_flips.flip", &flips) || section_find_stored("k.star_orders.order", &orders) ||
        section_find_stored("k.conv2ml.matrix", &matrices) || section_find_stored("k.space_group_slots.conv2ml", &matrix_index))
        return 1;
    static unsigned char present[SLOTS];
    for (int64_t g = 0; g < section_count(&irr_sg); g++) {
        int64_t sg = section_int(&irr_sg, g), kslot = section_int(&irr_k, g);
        if (sg >= 1 && sg <= 230 && kslot >= 1 && kslot <= 27) present[(sg - 1) * 27 + kslot - 1] = 1;
    }
    int64_t nmatrices = section_count(&matrices) / 16, star_cap = 1024, special_cap = 1024;
    K.star = malloc(sizeof(int64_t) * 16 * (size_t)star_cap);
    K.special = malloc(sizeof(int64_t) * 4 * (size_t)special_cap);
    if (!K.star || !K.special) return 1;
    rat_overflow = 0;
    int64_t next_flip = 0, next_order = 0;
    for (int slot = 0; slot < SLOTS; slot++) {
        if (!present[slot]) continue;
        int sg = slot / 27 + 1, kslot = slot % 27 + 1;
        int64_t lslot = (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1, dim = section_int(&k_dim, lslot), first[16];
        for (int i = 0; i < 16; i++) first[i] = section_int(&little_k, 16 * lslot + i);
        iarm_t arms[MAX_ARMS], ordered[MAX_ARMS];
        int n, minus;
        int64_t den;
        if (generated(sg, dim, first, arms, &n, &minus, &den)) return 1;
        memcpy(ordered, arms, sizeof(iarm_t) * (size_t)n);
        /* both convention lists are in (sg, k slot) order, as this loop */
        while (next_order + 2 < section_count(&orders) &&
               (section_int(&orders, next_order) - 1) * 27 + section_int(&orders, next_order + 1) - 1 < slot)
            next_order += 3 + section_int(&orders, next_order + 2);
        if (next_order + 2 < section_count(&orders) && section_int(&orders, next_order) == sg &&
            section_int(&orders, next_order + 1) == kslot) {
            if (section_int(&orders, next_order + 2) != n) return 1;
            for (int a = 0; a < n; a++) memcpy(ordered[a], arms[section_int(&orders, next_order + 3 + a) - 1], sizeof(iarm_t));
        }
        while (next_flip + 2 < section_count(&flips) &&
               (section_int(&flips, next_flip) - 1) * 27 + section_int(&flips, next_flip + 1) - 1 <= slot) {
            if (section_int(&flips, next_flip) == sg && section_int(&flips, next_flip + 1) == kslot) {
                int a = (int)section_int(&flips, next_flip + 2) - 1;
                if (a < 0 || a >= n) return 1;
                for (int r = 0; r < 4; r++)
                    for (int c = 0; c < 3; c++) ordered[a][r][c] = -ordered[a][r][c];
            }
            next_flip += 3;
        }
        int64_t base = dim > 0 ? K.nstar : K.nspecial;
        if (dim > 0 && base + n + nmatrices > star_cap) {
            while (base + n + nmatrices > star_cap) star_cap *= 2;
            int64_t *grown = realloc(K.star, sizeof(int64_t) * 16 * (size_t)star_cap);
            if (!grown) return 1;
            K.star = grown;
        }
        if (dim <= 0 && base + n > special_cap) {
            while (base + n > special_cap) special_cap *= 2;
            int64_t *grown = realloc(K.special, sizeof(int64_t) * 4 * (size_t)special_cap);
            if (!grown) return 1;
            K.special = grown;
        }
        K.pointer[slot] = base + 1;
        K.count[slot] = n;
        K.minusk[slot] = minus;
        for (int a = 0; a < n; a++) {
            for (int c = 0; c < 3; c++) ordered[a][0][c] = ((ordered[a][0][c] % den) + den) % den;
            if (dim > 0) {
                for (int r = 0; r < 4; r++) encode_row(ordered[a][r], den, K.star + 16 * (base + a) + 4 * r);
            } else {
                encode_row(ordered[a][0], den, K.special + 4 * (base + a));
            }
        }
        if (dim > 0) K.nstar += n;
        else K.nspecial += n;
    }
    if (rat_overflow) return 1;
    if (K.nstar + nmatrices > star_cap) {
        int64_t *grown = realloc(K.star, sizeof(int64_t) * 16 * (size_t)(K.nstar + nmatrices));
        if (!grown) return 1;
        K.star = grown;
    }
    /* the conventional-to-ML parameter matrices follow the arms in the star array, as in BYU's */
    for (int64_t j = 0; j < 16 * nmatrices; j++) K.star[16 * K.nstar + j] = section_int(&matrices, j);
    for (int slot = 0; slot < SLOTS; slot++) {
        int64_t index = slot < section_count(&matrix_index) ? section_int(&matrix_index, slot) : 0;
        K.conv2ml[slot] = index > 0 ? K.nstar + index : 0;
    }
    K.nstar += nmatrices;
    struct {
        const char *name;
        int64_t *data, count;
    } made[6] = {{"k.star.arms", K.star, K.nstar * 16},
                 {"k.star.special_arms", K.special, K.nspecial * 4},
                 {"k.star.arm_start", K.pointer, SLOTS},
                 {"k.star.arm_count", K.count, SLOTS},
                 {"k.star.minus_k", K.minusk, SLOTS},
                 {"k.star.conv2ml_start", K.conv2ml, SLOTS}};
    for (int i = 0; i < 6; i++) {
        blob_section *e = &K.entries[i];
        strncpy(e->name, made[i].name, BLOB_NAME_SIZE - 1);
        e->kind = BLOB_I64;
        e->count = (uint64_t)made[i].count;
        e->bytes = (uint64_t)made[i].count * sizeof(int64_t);
    }
    return 0;
}

/* section_find for the star sections; 1 when name is not one of them */
int star_rule_section(const char *name, section_view *view, int *handled)
{
    static const char *const names[6] = {"k.star.arms",        "k.star.special_arms",
                                         "k.star.arm_start", "k.star.arm_count",
                                         "k.star.minus_k", "k.star.conv2ml_start"};
    int64_t *data[6];
    *handled = 0;
    for (int i = 0; i < 6; i++)
        if (!strncmp(name, names[i], BLOB_NAME_SIZE)) {
            *handled = 1;
            if (!K.ready && !K.failed) {
                K.failed = build();
                K.ready = !K.failed;
            }
            if (!K.ready) return 1;
            data[0] = K.star, data[1] = K.special, data[2] = K.pointer, data[3] = K.count, data[4] = K.minusk, data[5] = K.conv2ml;
            view->entry = &K.entries[i];
            view->data = (const unsigned char *)data[i];
            return 0;
        }
    return 1;
}
