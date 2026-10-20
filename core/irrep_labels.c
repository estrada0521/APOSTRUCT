/* Irrep labels assembled from parts: the k name, the local number counted within the inversion parity, the parity sign
 * from the exact matrix of the inversion, and for a pair the partner under the k name or its conjugate alias. */

#include <stdlib.h>
#include <string.h>

#include "cyclo.h"
#include "geometry.h"
#include "irreps.h"
#include "tables.h"

int exact_record_matrix(int gid, int dim, int arms, const int64_t *k, int subtract_generated, int element_count,
                        const int *elements, int position, const int *record, int want_m, const cy_field **field_out, cy *w);

enum { LABEL_CAP = 24 };

static struct {
    section_view sg, kslot, number, type, partner, old, k_label;
    int ready, count;
    signed char *parity; /* 2 not computed */
    char (*full)[LABEL_CAP], (*old_label)[LABEL_CAP];
    int *gid_of_old, old_count;
} L;

static int load(void)
{
    if (L.ready) return 0;
    if (geo_load() || section_find("little.little_irr_space_group", &L.sg) || section_find("little.little_irr_k", &L.kslot) ||
        section_find("little.little_irr_irrnum", &L.number) || section_find("little.little_irr_type", &L.type) ||
        section_find("little.little_irr_cc", &L.partner) || section_find("little.little_irr_old", &L.old) ||
        section_find("little.little_k_label", &L.k_label))
        return 1;
    L.count = (int)section_count(&L.number);
    for (int g = 0; g < L.count; g++)
        if (section_int(&L.old, g) > L.old_count) L.old_count = (int)section_int(&L.old, g);
    L.parity = malloc((size_t)L.count + 1);
    L.full = calloc((size_t)L.count + 1, LABEL_CAP);
    L.old_label = calloc((size_t)L.old_count + 1, LABEL_CAP);
    L.gid_of_old = calloc((size_t)L.old_count + 1, sizeof(int));
    if (!L.parity || !L.full || !L.old_label || !L.gid_of_old) return 1;
    memset(L.parity, 2, (size_t)L.count + 1);
    for (int g = L.count - 1; g >= 0; g--) {
        int old = (int)section_int(&L.old, g);
        if (old > 0) L.gid_of_old[old] = g + 1;
    }
    L.ready = 1;
    return 0;
}

/* +1 or -1 when the inversion's first-arm block is the identity or its negative, 0 otherwise */
static int parity(int gid)
{
    if (L.parity[gid] != 2) return L.parity[gid];
    int result = 0;
    const int64_t none[1] = {0};
    irrep_context c;
    if (!irrep_context_for(gid, 0, none, &c)) {
        int position = 0;
        for (int q = 0; q < c.element_count && !position; q++) {
            int64_t m[9];
            geo_point_op_matrix(c.elements[5 * q + 4], m);
            int inversion = 1;
            for (int i = 0; i < 9; i++) inversion &= m[i] == (i % 4 == 0 ? -1 : 0);
            if (inversion) position = q + 1;
        }
        cy *w = position ? calloc((size_t)c.dim * c.dim, sizeof(cy)) : 0;
        const cy_field *f = 0;
        cy_overflow = 0;
        if (w && !exact_record_matrix(gid, c.dim, c.arms, c.k, c.subtract, c.element_count, c.elements, position,
                                      c.elements + 5 * (position - 1), 24, &f, w) && !cy_overflow) {
            int block = c.arms ? c.dim / c.arms : c.dim;
            cy trace, ref, t;
            cy_zero(f, &trace);
            for (int i = 0; i < block; i++) cy_add(f, &trace, w + i * c.dim + i, &trace);
            cy_rat(f, block, 1, &ref);
            cy_sub(f, &trace, &ref, &t);
            if (cy_is_zero(f, &t)) result = 1;
            cy_add(f, &trace, &ref, &t);
            if (cy_is_zero(f, &t)) result = -1;
        }
        free(w);
    }
    L.parity[gid] = (signed char)result;
    return result;
}

static void put(char *out, const char *s, int64_t n)
{
    size_t used = strlen(out);
    while (n > 0 && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    while (n > 0 && (*s == ' ' || *s == '\t')) s++, n--;
    if (used + (size_t)n < LABEL_CAP) memcpy(out + used, s, (size_t)n), out[used + n] = 0;
}

static void k_name(int sg, int kslot, char *out)
{
    const char *s;
    int64_t n;
    out[0] = 0;
    if (!section_string(&L.k_label, (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1, &s, &n)) put(out, s, n);
}

/* the local label of constituent j of (sg, kslot): its rank among the constituents of the same parity, then the sign */
static void local_label(int sg, int kslot, int j, char *out)
{
    signed char par[64];
    memset(par, 3, sizeof par);
    for (int g = 0; g < L.count; g++) {
        if (section_int(&L.sg, g) != sg || section_int(&L.kslot, g) != kslot) continue;
        int number = (int)section_int(&L.number, g), partner = (int)section_int(&L.partner, g);
        if (number < 1 || number >= 64) continue;
        int p = parity(g + 1);
        par[number] = (signed char)p;
        if (partner > 0 && partner < 64) par[partner] = (signed char)p;
    }
    int rank = 1;
    for (int i = 1; i < j && i < 64; i++) rank += par[i] == par[j];
    char text[8];
    int n = 0;
    if (rank >= 10) text[n++] = (char)('0' + rank / 10);
    text[n++] = (char)('0' + rank % 10);
    if (par[j] == 1) text[n++] = '+';
    if (par[j] == -1) text[n++] = '-';
    put(out, text, n);
}

static const char *conjugate_k(int sg, const char *k, char *out)
{
    static const char *const pairs[][2] = {{"GP", "GQ"}, {"LD", "LE"}, {"SM", "SN"}, {"DT", "DU"}};
    for (size_t i = 0; i < sizeof pairs / sizeof pairs[0]; i++)
        if (!strcmp(k, pairs[i][0])) return strcpy(out, pairs[i][1]);
    /* BYU names these conjugate points with C instead of A */
    static const struct {
        int sg;
        const char *k;
    } c_names[] = {{143, "PBCDE"}, {144, "PBCDE"}, {145, "PBCDE"}, {156, "D"}, {157, "PC"}, {158, "D"}, {159, "PC"}, {174, "BE"}};
    char suffix = 'A';
    for (size_t i = 0; i < sizeof c_names / sizeof c_names[0]; i++)
        if (c_names[i].sg == sg && strlen(k) == 1 && strchr(c_names[i].k, k[0])) suffix = 'C';
    strcpy(out, k);
    size_t n = strlen(out);
    out[n] = suffix, out[n + 1] = 0;
    return out;
}

/* the Miller-Love label of a little irrep, e.g. GM1+, W1W2, P1PC1 */
const char *irrep_full_label(int gid)
{
    if (load() || gid < 1 || gid > L.count) return "";
    char *out = L.full[gid];
    if (out[0]) return out;
    int g = gid - 1, sg = (int)section_int(&L.sg, g), kslot = (int)section_int(&L.kslot, g), type = (int)section_int(&L.type, g);
    int j = (int)section_int(&L.number, g), partner = (int)section_int(&L.partner, g);
    if (!partner) partner = j;
    char k[8], other[8], local[8];
    k_name(sg, kslot, k);
    local[0] = 0;
    local_label(sg, kslot, j, local);
    put(out, k, (int64_t)strlen(k));
    put(out, local, (int64_t)strlen(local));
    if (type == 2) {
        put(out, k, (int64_t)strlen(k));
        put(out, local, (int64_t)strlen(local));
    } else if (type == 3) {
        local[0] = 0;
        local_label(sg, kslot, partner, local);
        const char *name = partner == j ? conjugate_k(sg, k, other) : k;
        put(out, name, (int64_t)strlen(name));
        put(out, local, (int64_t)strlen(local));
    }
    return out;
}

/* the label of an old (data_irreps) irrep: its little irrep's parts under the original k name for both constituents */
const char *irrep_old_label(int old_id)
{
    if (load() || old_id < 1 || old_id > L.old_count || !L.gid_of_old[old_id]) return "";
    char *out = L.old_label[old_id];
    if (out[0]) return out;
    int g = L.gid_of_old[old_id] - 1, sg = (int)section_int(&L.sg, g), kslot = (int)section_int(&L.kslot, g);
    int type = (int)section_int(&L.type, g), j = (int)section_int(&L.number, g), partner = (int)section_int(&L.partner, g);
    if (!partner) partner = j;
    char k[8], local[8];
    k_name(sg, kslot, k);
    local[0] = 0;
    local_label(sg, kslot, j, local);
    put(out, k, (int64_t)strlen(k));
    put(out, local, (int64_t)strlen(local));
    if (type > 1) {
        local[0] = 0;
        local_label(sg, kslot, type == 2 ? j : partner, local);
        put(out, k, (int64_t)strlen(k));
        put(out, local, (int64_t)strlen(local));
    }
    return out;
}
