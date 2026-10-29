/* The Local theorems on one apo mode-details payload: mode basis, group invariance and fixed-space dimension at the
 * display precision of the Python checks in verification/mathematics; the k-star sector, the isotropy index and the
 * isotypic orthogonality of different irreps, which take the parent structure's operations from the child cell's atoms
 * and read no tables. Shares no code with core/. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../json.h"

#define PSCALE 100000

enum { SAT = 0, REF = 1, IND = 2, NA = 3 };

typedef struct {
    int n, cap;
    char **s;
    int *len;
} strtab;

typedef struct {
    int n;
    int *atom, *pos;
} block_t;

typedef struct {
    int kind, label_atom, nf_zero, block, nrows, has_k, primary, irrep;
    int *vec;
    double k[3];
} def_t;

typedef struct {
    int label, multiplicity, rep[3], nexp;
    int *exp;
} orbit_t;

typedef struct {
    char err[256];
    int natoms, nblocks, ndefs, nstrain, norbits, failed;
    int *atom_marker, *atom_pos, *atom_label;
    block_t *blocks;
    def_t *defs;
    int64_t (*strain)[6];
    int *strain_zero;
    orbit_t *orbits;
    strtab markers, label_atoms, irreps;
    char *subgroup;
    int has_basis, has_origin, has_cell, index;
    double basis[9], origin[3], cell[6];
} mc_t;

typedef struct {
    int r[9];
    int64_t tn[3], td[3];
    int tr;
} op_t;

static int fail(mc_t *m, const char *msg)
{
    if (!m->failed) snprintf(m->err, sizeof m->err, "%s", msg);
    m->failed = 1;
    return -1;
}

/* round-to-display as Python format(x, f".{digits}f"), read back as an integer */
static int quantize(double x, int digits, int64_t *out)
{
    char buf[400];
    snprintf(buf, sizeof buf, "%.*f", digits, x);
    int64_t v = 0;
    int neg = 0;
    for (const char *s = buf; *s; s++) {
        if (*s == '-') neg = 1;
        else if (*s >= '0' && *s <= '9') {
            if (v > INT64_MAX / 10 - 10) return -1;
            v = 10 * v + (*s - '0');
        }
    }
    *out = neg ? -v : v;
    return 0;
}

static int qpos(cur *c, int *out)
{
    double x;
    int64_t q;
    if (number(c, &x) || quantize(x, 5, &q)) return -1;
    q %= PSCALE;
    if (q < 0) q += PSCALE;
    *out = (int)q;
    return 0;
}

static int qvec(cur *c, int *out)
{
    double x;
    int64_t q;
    if (number(c, &x) || quantize(x, 4, &q) || q > INT32_MAX || q < -INT32_MAX) return -1;
    *out = (int)q;
    return 0;
}


/* a rational "p", "p/q" or decimal at s, advancing past it; -1 when there is none */
static int rational(const char **s, const char *end, double *out)
{
    const char *p = *s;
    while (p < end && *p == ' ') p++;
    char *stop;
    double a = strtod(p, &stop);
    if (stop == p) return -1;
    p = stop;
    if (p < end && *p == '/') {
        double b = strtod(p + 1, &stop);
        if (stop == p + 1 || b == 0) return -1;
        a /= b;
        p = stop;
    }
    *out = a;
    *s = p;
    return 0;
}

/* three rationals separated by commas, starting after `open` at s; the k of a label "[k1,k2,k3]" or an origin "(x,y,z)" */
static int triple(const char *s, const char *end, char open, double *out)
{
    while (s < end && *s != open) s++;
    if (s == end) return -1;
    s++;
    for (int a = 0; a < 3; a++) {
        if (rational(&s, end, out + a)) return -1;
        while (s < end && *s == ' ') s++;
        if (a < 2 && (s == end || *s++ != ',')) return -1;
    }
    return 0;
}

static int intern(strtab *t, const char *s, int len)
{
    for (int i = t->n - 1; i >= 0; i--)
        if (t->len[i] == len && !memcmp(t->s[i], s, len)) return i;
    if (t->n == t->cap) {
        t->cap = t->cap ? 2 * t->cap : 64;
        t->s = realloc(t->s, sizeof(char *) * t->cap);
        t->len = realloc(t->len, sizeof(int) * t->cap);
    }
    t->s[t->n] = malloc(len + 1);
    memcpy(t->s[t->n], s, len);
    t->s[t->n][len] = 0;
    t->len[t->n] = len;
    return t->n++;
}

/* the one [atom:site:kind] block of a mode label, as re.finditer(r"\[([^:\]]+):([^:\]]+):(dsp|mag|occ)]") */
static int label_atom(const char *s, int len, int kind, int *atom_start, int *atom_len)
{
    static const char *kinds[] = {"dsp", "mag", "occ"};
    int found = 0, i = 0;
    while (i < len) {
        if (s[i] != '[') {
            i++;
            continue;
        }
        int a = i + 1, j = a;
        while (j < len && s[j] != ':' && s[j] != ']') j++;
        if (j == a || j >= len || s[j] != ':') {
            i++;
            continue;
        }
        int b = j + 1, k = b;
        while (k < len && s[k] != ':' && s[k] != ']') k++;
        if (k == b || k >= len || s[k] != ':' || k + 4 >= len + 0 || s[k + 4] != ']') {
            i++;
            continue;
        }
        int which = -1;
        for (int w = 0; w < 3; w++)
            if (!memcmp(s + k + 1, kinds[w], 3)) which = w;
        if (which < 0) {
            i++;
            continue;
        }
        if (found++) return -1;
        if (which != kind) return -1;
        *atom_start = a;
        *atom_len = j - a;
        i = k + 5;
    }
    return found == 1 ? 0 : -1;
}

typedef struct {
    int64_t *keys;
    int *vals;
    int cap;
} hmap;

static uint64_t mix(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    return x ^ (x >> 33);
}

static int64_t pkey(int atom, const int *p) { return (((int64_t)atom * PSCALE + p[0]) * PSCALE + p[1]) * PSCALE + p[2]; }

static void hm_init(hmap *h, int n)
{
    h->cap = 16;
    while (h->cap < 2 * n + 2) h->cap *= 2;
    h->keys = malloc(sizeof(int64_t) * h->cap);
    h->vals = malloc(sizeof(int) * h->cap);
    for (int i = 0; i < h->cap; i++) h->vals[i] = -1;
}

static int hm_put(hmap *h, int64_t k, int v)
{
    uint64_t i = mix((uint64_t)k) & (h->cap - 1);
    while (h->vals[i] >= 0) {
        if (h->keys[i] == k) return h->vals[i];
        i = (i + 1) & (h->cap - 1);
    }
    h->keys[i] = k;
    h->vals[i] = v;
    return -1;
}

static int hm_get(const hmap *h, int64_t k)
{
    uint64_t i = mix((uint64_t)k) & (h->cap - 1);
    while (h->vals[i] >= 0) {
        if (h->keys[i] == k) return h->vals[i];
        i = (i + 1) & (h->cap - 1);
    }
    return -1;
}

static void hm_free(hmap *h)
{
    free(h->keys);
    free(h->vals);
}

static int same_block(const block_t *a, int n, const int *atom, const int *pos);

static int add_block(mc_t *m, int n, int *atom, int *pos)
{
    for (int b = 0; b < m->nblocks; b++)
        if (same_block(&m->blocks[b], n, atom, pos)) {
            free(atom);
            free(pos);
            return b;
        }
    hmap h;
    hm_init(&h, n);
    for (int i = 0; i < n; i++)
        if (hm_put(&h, pkey(atom[i], pos + 3 * i), i) >= 0) {
            hm_free(&h);
            free(atom);
            free(pos);
            return fail(m, "mode definition repeats an atom-position row");
        }
    hm_free(&h);
    m->blocks = realloc(m->blocks, sizeof(block_t) * (m->nblocks + 1));
    m->blocks[m->nblocks] = (block_t){n, atom, pos};
    return m->nblocks++;
}

/* carried atoms of a table block: a row without its own marker continues the previous row's atom */
static int table_block(mc_t *m, cur *c)
{
    int first = 1, n = 0, cap = 16, carried = -1;
    int *atom = malloc(sizeof(int) * cap), *pos = malloc(sizeof(int) * 3 * cap);
    if (arr_begin(c)) return fail(m, "atom block is not a list");
    for (int r; (r = arr_next(c, &first)) == 1;) {
        double x;
        if (number(c, &x) || x < 0 || x >= m->natoms || x != floor(x)) return fail(m, "atom block index");
        int i = (int)x;
        if (m->atom_marker[i] >= 0) carried = m->atom_marker[i];
        if (carried < 0) return fail(m, "mode definition starts with an unlabeled atom row");
        if (n == cap) {
            cap *= 2;
            atom = realloc(atom, sizeof(int) * cap);
            pos = realloc(pos, sizeof(int) * 3 * cap);
        }
        atom[n] = carried;
        memcpy(pos + 3 * n, m->atom_pos + 3 * i, sizeof(int) * 3);
        n++;
    }
    return add_block(m, n, atom, pos);
}

static int read_atoms(mc_t *m, cur *c)
{
    int first = 1, cap = 64;
    m->atom_marker = malloc(sizeof(int) * cap);
    m->atom_label = malloc(sizeof(int) * cap);
    m->atom_pos = malloc(sizeof(int) * 3 * cap);
    if (arr_begin(c)) return fail(m, "atoms is not a list");
    for (int r; (r = arr_next(c, &first)) == 1;) {
        if (m->natoms == cap) {
            cap *= 2;
            m->atom_marker = realloc(m->atom_marker, sizeof(int) * cap);
            m->atom_label = realloc(m->atom_label, sizeof(int) * cap);
            m->atom_pos = realloc(m->atom_pos, sizeof(int) * 3 * cap);
        }
        int i = m->natoms++, f = 1, seen = 0;
        m->atom_marker[i] = m->atom_label[i] = -1;
        const char *key;
        int klen;
        if (obj_begin(c)) return fail(m, "atom row");
        for (int k; (k = obj_next(c, &key, &klen, &f)) == 1;) {
            if (keyis(key, klen, "atom") && !is_null(c)) {
                const char *s;
                int len;
                if (string_value(c, &s, &len) || !len) return fail(m, "atom marker");
                m->atom_marker[i] = intern(&m->markers, s, len);
            } else if (keyis(key, klen, "label") && !is_null(c)) {
                const char *s;
                int len;
                if (string_value(c, &s, &len) || !len) return fail(m, "atom label");
                m->atom_label[i] = intern(&m->markers, s, len);
            } else if (keyis(key, klen, "xyz")) {
                int af = 1, a = 0;
                if (arr_begin(c)) return fail(m, "atom xyz");
                for (int q; (q = arr_next(c, &af)) == 1; a++)
                    if (a >= 3 || qpos(c, m->atom_pos + 3 * i + a)) return fail(m, "atom xyz");
                if (a != 3) return fail(m, "atom xyz");
                seen = 1;
            } else if (skip(c)) return fail(m, "atom row");
        }
        if (!seen) return fail(m, "atom without xyz");
    }
    return 0;
}

static int read_vectors(mc_t *m, cur *c, int kind, int **out, int *count)
{
    int first = 1, n = 0, cap = 16;
    int *v = malloc(sizeof(int) * 3 * cap);
    if (arr_begin(c)) return fail(m, "definition vectors");
    for (int r; (r = arr_next(c, &first)) == 1; n++) {
        if (n == cap) {
            cap *= 2;
            v = realloc(v, sizeof(int) * 3 * cap);
        }
        ws(c);
        if (*c->p == '[') {
            int af = 1, a = 0;
            arr_begin(c);
            for (int q; (q = arr_next(c, &af)) == 1; a++)
                if (a >= 3 || qvec(c, v + 3 * n + a)) return fail(m, "definition vector");
            if (a != 3) return fail(m, "definition vector");
        } else {
            if (kind != 2 || qvec(c, v + 3 * n)) return fail(m, "definition vector");
            v[3 * n + 1] = v[3 * n + 2] = 0;
        }
    }
    *out = v;
    *count = n;
    return 0;
}

static int same_block(const block_t *a, int n, const int *atom, const int *pos)
{
    return a->n == n && !memcmp(a->atom, atom, sizeof(int) * n) && !memcmp(a->pos, pos, sizeof(int) * 3 * n);
}

/* the older per-definition form: rows [{atom, xyz, dxyz}] */
static int row_block(mc_t *m, cur *c, int kind, def_t *d)
{
    int first = 1, n = 0, cap = 16, carried = -1;
    int *atom = malloc(sizeof(int) * cap), *pos = malloc(sizeof(int) * 3 * cap), *vec = malloc(sizeof(int) * 3 * cap);
    if (arr_begin(c)) return fail(m, "definition rows");
    for (int r; (r = arr_next(c, &first)) == 1; n++) {
        if (n == cap) {
            cap *= 2;
            atom = realloc(atom, sizeof(int) * cap);
            pos = realloc(pos, sizeof(int) * 3 * cap);
            vec = realloc(vec, sizeof(int) * 3 * cap);
        }
        int f = 1, seen = 0;
        const char *key;
        int klen;
        if (obj_begin(c)) return fail(m, "definition row");
        for (int k; (k = obj_next(c, &key, &klen, &f)) == 1;) {
            if (keyis(key, klen, "atom") && !is_null(c)) {
                const char *s;
                int len;
                if (string_value(c, &s, &len) || !len) return fail(m, "row atom");
                carried = intern(&m->markers, s, len);
            } else if (keyis(key, klen, "xyz") || keyis(key, klen, "dxyz")) {
                int is_pos = klen == 3, af = 1, a = 0;
                if (arr_begin(c)) return fail(m, "row vector");
                for (int q; (q = arr_next(c, &af)) == 1; a++)
                    if (a >= 3 || (is_pos ? qpos(c, pos + 3 * n + a) : qvec(c, vec + 3 * n + a))) return fail(m, "row vector");
                if (a != 3) return fail(m, "row vector");
                seen |= is_pos ? 1 : 2;
            } else if (skip(c)) return fail(m, "definition row");
        }
        if (seen != 3) return fail(m, "row needs xyz and dxyz");
        if (carried < 0) return fail(m, "mode definition starts with an unlabeled atom row");
        atom[n] = carried;
    }
    (void)kind;
    d->vec = vec;
    d->nrows = n;
    d->block = add_block(m, n, atom, pos);
    return d->block < 0 ? -1 : 0;
}

static int nf_is_zero(double nf)
{
    int64_t q;
    return !quantize(nf, 5, &q) && q == 0;
}

static int read_definitions(mc_t *m, cur *c, int kind, int *blocks_raw, int nblocks_raw)
{
    int first = 1;
    if (arr_begin(c)) return fail(m, "definitions is not a list");
    for (int r; (r = arr_next(c, &first)) == 1;) {
        m->defs = realloc(m->defs, sizeof(def_t) * (m->ndefs + 1));
        def_t *d = &m->defs[m->ndefs++];
        memset(d, 0, sizeof *d);
        d->kind = kind;
        d->block = -1;
        int f = 1, has_nf = 0, has_label = 0, raw_block = -1;
        const char *key;
        int klen;
        cur vectors = {0};
        if (obj_begin(c)) return fail(m, "definition");
        for (int k; (k = obj_next(c, &key, &klen, &f)) == 1;) {
            if (keyis(key, klen, "label")) {
                const char *s;
                int len, a, alen;
                if (string_value(c, &s, &len) || !len) return fail(m, "definition label");
                if (label_atom(s, len, kind, &a, &alen) || !alen) return fail(m, "mode label must contain one atom/site/kind block");
                d->label_atom = intern(&m->label_atoms, s + a, alen);
                d->has_k = !triple(s, s + a, '[', d->k);
                d->irrep = -1;
                for (int q = 0; q < a; q++)
                    if (s[q] == ']') {
                        int e = q + 1;
                        while (e < a && s[e] != '(') e++;
                        if (e > q + 1) d->irrep = intern(&m->irreps, s + q + 1, e - q - 1);
                        break;
                    }
                has_label = 1;
            } else if (keyis(key, klen, "role")) {
                const char *s;
                int len;
                if (string_value(c, &s, &len)) return fail(m, "definition role");
                d->primary = len == 7 && !memcmp(s, "primary", 7);
            } else if (keyis(key, klen, "normfactor")) {
                double nf;
                if (number(c, &nf)) return fail(m, "normfactor");
                d->nf_zero = nf_is_zero(nf);
                has_nf = 1;
            } else if (keyis(key, klen, "atom_block")) {
                double x;
                if (number(c, &x) || x < 0 || x >= nblocks_raw || x != floor(x)) return fail(m, "atom_block");
                raw_block = (int)x;
            } else if (keyis(key, klen, "vectors")) {
                vectors = *c;
                if (skip(c)) return fail(m, "vectors");
            } else if (keyis(key, klen, "rows")) {
                if (row_block(m, c, kind, d)) return -1;
            } else if (skip(c)) return fail(m, "definition");
        }
        if (!has_label || !has_nf) return fail(m, "definition needs label and normfactor");
        if (d->block < 0) {
            if (raw_block < 0 || !vectors.p) return fail(m, "definition needs rows or atom_block and vectors");
            d->block = blocks_raw[raw_block];
            if (read_vectors(m, &vectors, kind, &d->vec, &d->nrows)) return -1;
            if (d->nrows != m->blocks[d->block].n) return fail(m, "vectors and atom block differ in length");
        }
        if (!d->nrows) return fail(m, "mode definition has no vector rows");
    }
    return 0;
}

static int read_strain(mc_t *m, cur *c)
{
    int first = 1;
    if (arr_begin(c)) return fail(m, "strain definitions");
    for (int r; (r = arr_next(c, &first)) == 1;) {
        m->strain = realloc(m->strain, sizeof(int64_t[6]) * (m->nstrain + 1));
        m->strain_zero = realloc(m->strain_zero, sizeof(int) * (m->nstrain + 1));
        int i = m->nstrain++, f = 1, seen = 0;
        const char *key;
        int klen;
        if (obj_begin(c)) return fail(m, "strain definition");
        for (int k; (k = obj_next(c, &key, &klen, &f)) == 1;) {
            if (keyis(key, klen, "components")) {
                int af = 1, a = 0;
                if (arr_begin(c)) return fail(m, "strain components");
                for (int q; (q = arr_next(c, &af)) == 1; a++) {
                    double x;
                    if (a >= 6 || number(c, &x) || quantize(x, 4, &m->strain[i][a])) return fail(m, "strain components");
                }
                if (a != 6) return fail(m, "strain definition needs six components");
                seen |= 1;
            } else if (keyis(key, klen, "normfactor")) {
                double nf;
                if (number(c, &nf)) return fail(m, "strain normfactor");
                m->strain_zero[i] = nf == 0;
                seen |= 2;
            } else if (skip(c)) return fail(m, "strain definition");
        }
        if (seen != 3) return fail(m, "strain definition needs components and normfactor");
    }
    return 0;
}

static int read_orbits(mc_t *m, cur *c)
{
    int first = 1;
    if (arr_begin(c)) return fail(m, "undistorted atoms");
    for (int r; (r = arr_next(c, &first)) == 1;) {
        m->orbits = realloc(m->orbits, sizeof(orbit_t) * (m->norbits + 1));
        orbit_t *o = &m->orbits[m->norbits++];
        memset(o, 0, sizeof *o);
        int f = 1, seen = 0;
        const char *key;
        int klen;
        if (obj_begin(c)) return fail(m, "undistorted atom orbit");
        for (int k; (k = obj_next(c, &key, &klen, &f)) == 1;) {
            if (keyis(key, klen, "label")) {
                const char *s;
                int len;
                if (string_value(c, &s, &len) || !len) return fail(m, "orbit label");
                o->label = intern(&m->markers, s, len);
                seen |= 1;
            } else if (keyis(key, klen, "multiplicity")) {
                double x;
                if (number(c, &x) || x <= 0 || x != floor(x)) return fail(m, "orbit multiplicity");
                o->multiplicity = (int)x;
                seen |= 2;
            } else if (keyis(key, klen, "xyz")) {
                int af = 1, a = 0;
                if (arr_begin(c)) return fail(m, "orbit xyz");
                for (int q; (q = arr_next(c, &af)) == 1; a++)
                    if (a >= 3 || qpos(c, o->rep + a)) return fail(m, "orbit xyz");
                if (a != 3) return fail(m, "orbit xyz");
                seen |= 4;
            } else if (keyis(key, klen, "site")) {
                const char *s;
                int len;
                if (string_value(c, &s, &len) || !len) return fail(m, "orbit site");
                seen |= 8;
            } else if (keyis(key, klen, "_presentation_orbit_points") && !is_null(c)) {
                int pf = 1, cap = 16;
                o->exp = malloc(sizeof(int) * 3 * cap);
                if (arr_begin(c)) return fail(m, "orbit points");
                for (int q; (q = arr_next(c, &pf)) == 1; o->nexp++) {
                    if (o->nexp == cap) {
                        cap *= 2;
                        o->exp = realloc(o->exp, sizeof(int) * 3 * cap);
                    }
                    int af = 1, a = 0;
                    if (arr_begin(c)) return fail(m, "orbit point");
                    for (int s; (s = arr_next(c, &af)) == 1; a++)
                        if (a >= 3 || qpos(c, o->exp + 3 * o->nexp + a)) return fail(m, "orbit point");
                    if (a != 3) return fail(m, "orbit point");
                }
            } else if (skip(c)) return fail(m, "undistorted atom orbit");
        }
        if (seen != 15) return fail(m, "undistorted atom orbit needs label, site, multiplicity and xyz");
    }
    if (!m->norbits) return fail(m, "no undistorted atom orbits");
    return 0;
}

static int find_mode_details(cur *c, const char *const *path, int depth, cur *out)
{
    const char *key;
    int klen, first = 1;
    if (obj_begin(c)) return -1;
    for (int k; (k = obj_next(c, &key, &klen, &first)) == 1;) {
        if (keyis(key, klen, path[0])) {
            if (depth == 1) {
                *out = *c;
                return 0;
            }
            return find_mode_details(c, path + 1, depth - 1, out);
        }
        if (skip(c)) return -1;
    }
    return -1;
}

static int load(mc_t *m, const char *text, int64_t len)
{
    static const char *const preview[] = {"preview", "selected", "mode_details"}, *const state[] = {"selected", "mode_details"};
    static const char *const preview_i[] = {"preview", "selected", "orderparam", "isotropy", "i"},
                             *const state_i[] = {"selected", "orderparam", "isotropy", "i"};
    cur root = {text, text + len}, md, iv;
    double index = 0;
    if (!find_mode_details(&root, preview_i, 5, &iv) || (root = (cur){text, text + len}, !find_mode_details(&root, state_i, 4, &iv)))
        if (!number(&iv, &index) && index >= 1 && index == floor(index)) m->index = (int)index;
    root = (cur){text, text + len};
    if (find_mode_details(&root, preview, 3, &md)) {
        root = (cur){text, text + len};
        if (find_mode_details(&root, state, 2, &md)) return fail(m, "no mode_details");
    }
    cur atoms = {0}, blocks = {0}, defs[3] = {{0}}, strain = {0}, orbits = {0}, subgroup = {0};
    const char *key;
    int klen, first = 1;
    cur c = md;
    if (obj_begin(&c)) return fail(m, "mode_details is not an object");
    for (int k; (k = obj_next(&c, &key, &klen, &first)) == 1;) {
        cur *slot = keyis(key, klen, "atoms") ? &atoms
                    : keyis(key, klen, "atom_blocks") ? &blocks
                    : keyis(key, klen, "displacive_definitions") ? &defs[0]
                    : keyis(key, klen, "magnetic_definitions") ? &defs[1]
                    : keyis(key, klen, "occupancy_definitions") ? &defs[2]
                    : keyis(key, klen, "strain_definitions") ? &strain
                    : keyis(key, klen, "undistorted_atoms") ? &orbits
                    : 0;
        if (keyis(key, klen, "subgroup_details")) {
            const char *sk;
            int sklen, sf = 1;
            cur s = c;
            if (obj_begin(&s)) return fail(m, "subgroup_details");
            for (int q; (q = obj_next(&s, &sk, &sklen, &sf)) == 1;) {
                if (keyis(sk, sklen, "magnetic_subgroup")) subgroup = s;
                if (keyis(sk, sklen, "presentation_basis") && !is_null(&s)) {
                    cur b = s;
                    int rf = 1, row = 0;
                    if (arr_begin(&b)) return fail(m, "presentation_basis");
                    for (int r; (r = arr_next(&b, &rf)) == 1; row++) {
                        int af = 1, a = 0;
                        if (row >= 3 || arr_begin(&b)) return fail(m, "presentation_basis");
                        for (int q; (q = arr_next(&b, &af)) == 1; a++)
                            if (a >= 3 || number(&b, m->basis + 3 * row + a)) return fail(m, "presentation_basis");
                        if (a != 3) return fail(m, "presentation_basis");
                    }
                    m->has_basis = row == 3;
                }
                if (keyis(sk, sklen, "origin") && !is_null(&s)) {
                    cur o = s;
                    const char *os;
                    int olen;
                    if (!string_value(&o, &os, &olen)) m->has_origin = !triple(os, os + olen, '(', m->origin);
                }
                if (skip(&s)) return fail(m, "subgroup_details");
            }
        }
        if (keyis(key, klen, "lattice") && !is_null(&c)) {
            static const char *const names[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
            cur l = c;
            const char *lk;
            int lklen, lf = 1, seen = 0;
            if (obj_begin(&l)) return fail(m, "lattice");
            for (int q; (q = obj_next(&l, &lk, &lklen, &lf)) == 1;) {
                int which = -1;
                for (int i = 0; i < 6; i++)
                    if (keyis(lk, lklen, names[i])) which = i;
                if (which >= 0) {
                    if (number(&l, m->cell + which)) return fail(m, "lattice");
                    seen |= 1 << which;
                } else if (skip(&l)) return fail(m, "lattice");
            }
            m->has_cell = seen == 63;
        }
        if (slot && !is_null(&c)) *slot = c;
        if (skip(&c)) return fail(m, "mode_details");
    }
    if (!subgroup.p) return fail(m, "no magnetic subgroup");
    cur s = subgroup;
    skip(&s);
    m->subgroup = malloc(s.p - subgroup.p + 1);
    memcpy(m->subgroup, subgroup.p, s.p - subgroup.p);
    m->subgroup[s.p - subgroup.p] = 0;
    if (!defs[0].p || !defs[1].p || !strain.p || !orbits.p) return fail(m, "mode_details lacks a definition list");
    if (atoms.p && read_atoms(m, &atoms)) return -1;
    int *raw = 0, nraw = 0;
    if (blocks.p) {
        int bf = 1;
        if (arr_begin(&blocks)) return fail(m, "atom_blocks");
        for (int r; (r = arr_next(&blocks, &bf)) == 1;) {
            raw = realloc(raw, sizeof(int) * (nraw + 1));
            if ((raw[nraw++] = table_block(m, &blocks)) < 0) return -1;
        }
    }
    for (int k = 0; k < 3; k++)
        if (defs[k].p && read_definitions(m, &defs[k], k, raw, nraw)) return -1;
    free(raw);
    if (read_strain(m, &strain) || read_orbits(m, &orbits)) return -1;
    return 0;
}

void *mc_open(const char *text, int64_t len)
{
    mc_t *m = calloc(1, sizeof *m);
    load(m, text, len);
    return m;
}

const char *mc_error(void *h) { return ((mc_t *)h)->failed ? ((mc_t *)h)->err : 0; }

const char *mc_subgroup(void *h) { return ((mc_t *)h)->subgroup; }

int mc_definitions(void *h) { return ((mc_t *)h)->ndefs + ((mc_t *)h)->nstrain; }

void mc_close(void *h)
{
    mc_t *m = h;
    for (int i = 0; i < m->nblocks; i++) {
        free(m->blocks[i].atom);
        free(m->blocks[i].pos);
    }
    for (int i = 0; i < m->ndefs; i++) free(m->defs[i].vec);
    for (int i = 0; i < m->norbits; i++) free(m->orbits[i].exp);
    for (strtab *t = &m->markers; t <= &m->irreps; t++) {
        for (int i = 0; i < t->n; i++) free(t->s[i]);
        free(t->s);
        free(t->len);
    }
    free(m->blocks);
    free(m->defs);
    free(m->orbits);
    free(m->strain);
    free(m->strain_zero);
    free(m->atom_marker);
    free(m->atom_label);
    free(m->atom_pos);
    free(m->subgroup);
    free(m);
}

/* ---- basis: rank modulo the primes of mode_basis.py, per connected component of the columns ---- */

static const int64_t primes[] = {1000003, 1000033, 1000037, 1000039, 1000081, 1000099, 1000117, 1000121};

static int64_t inv_mod(int64_t a, int64_t p)
{
    int64_t r = 1, e = p - 2;
    a %= p;
    while (e) {
        if (e & 1) r = r * a % p;
        a = a * a % p;
        e >>= 1;
    }
    return r;
}

/* rank of the rows x cols matrix (column major, entries already reduced) modulo p; destroys a */
static int rank_mod(int64_t *a, int rows, int cols, int64_t p)
{
    int rank = 0;
    for (int j = 0; j < cols && rank < rows; j++) {
        int64_t *cj = a + (int64_t)j * rows;
        int piv = -1;
        for (int i = rank; i < rows; i++)
            if (cj[i]) {
                piv = i;
                break;
            }
        if (piv < 0) continue;
        if (piv != rank)
            for (int k = j; k < cols; k++) {
                int64_t *ck = a + (int64_t)k * rows, t = ck[piv];
                ck[piv] = ck[rank];
                ck[rank] = t;
            }
        int64_t inv = inv_mod(cj[rank], p);
        for (int k = j + 1; k < cols; k++) {
            int64_t *ck = a + (int64_t)k * rows;
            if (!ck[rank]) continue;
            int64_t f = ck[rank] * inv % p;
            for (int i = rank + 1; i < rows; i++)
                if (cj[i]) ck[i] = ((ck[i] - f * cj[i]) % p + p) % p;
            ck[rank] = 0;
        }
        rank++;
    }
    return rank;
}

typedef struct {
    int n;
    int *coord;
    int64_t *val;
} column_t;

static int uf_find(int *parent, int x)
{
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
}

/* full column rank modulo some prime of the list, as _group_certificate; 1 full, 0 deficient */
static int full_rank(column_t *cols, int ncols, int ncoords)
{
    if (!ncols) return 1;
    int *parent = malloc(sizeof(int) * ncols), *owner = malloc(sizeof(int) * ncoords);
    for (int j = 0; j < ncols; j++) parent[j] = j;
    for (int i = 0; i < ncoords; i++) owner[i] = -1;
    for (int j = 0; j < ncols; j++)
        for (int e = 0; e < cols[j].n; e++) {
            int i = cols[j].coord[e];
            if (owner[i] < 0) owner[i] = j;
            else {
                int a = uf_find(parent, owner[i]), b = uf_find(parent, j);
                if (a != b) parent[b] = a;
            }
        }
    int *start = calloc(ncols + 1, sizeof(int)), *order = malloc(sizeof(int) * ncols), *local = malloc(sizeof(int) * ncoords);
    for (int j = 0; j < ncols; j++) start[uf_find(parent, j) + 1]++;
    for (int j = 0; j < ncols; j++) start[j + 1] += start[j];
    int *fill = malloc(sizeof(int) * ncols);
    memcpy(fill, start, sizeof(int) * ncols);
    for (int j = 0; j < ncols; j++) order[fill[uf_find(parent, j)]++] = j;
    for (int i = 0; i < ncoords; i++) local[i] = -1;
    int result = 0;
    for (int q = 0; q < 8 && !result; q++) {
        int64_t p = primes[q];
        int total = 0;
        for (int r = 0; r < ncols; r++) {
            int a = start[r], b = start[r + 1], nc = b - a, nr = 0;
            if (!nc) continue;
            for (int t = a; t < b; t++)
                for (int e = 0; e < cols[order[t]].n; e++)
                    if (local[cols[order[t]].coord[e]] < 0) local[cols[order[t]].coord[e]] = nr++;
            int64_t *mat = calloc((size_t)nr * nc + 1, sizeof(int64_t));
            for (int t = a; t < b; t++)
                for (int e = 0; e < cols[order[t]].n; e++) {
                    int64_t v = cols[order[t]].val[e] % p;
                    mat[(int64_t)(t - a) * nr + local[cols[order[t]].coord[e]]] = v < 0 ? v + p : v;
                }
            total += rank_mod(mat, nr, nc, p);
            free(mat);
            for (int t = a; t < b; t++)
                for (int e = 0; e < cols[order[t]].n; e++) local[cols[order[t]].coord[e]] = -1;
        }
        result = total == ncols;
    }
    free(parent);
    free(owner);
    free(start);
    free(order);
    free(fill);
    free(local);
    return result;
}

int mc_basis(void *h)
{
    mc_t *m = h;
    if (m->failed) return -1;
    int ngroups = 0, *group_kind = 0, *group_atom = 0, *def_group = malloc(sizeof(int) * (m->ndefs + 1));
    for (int i = 0; i < m->ndefs; i++) {
        int g = 0;
        while (g < ngroups && !(group_kind[g] == m->defs[i].kind && group_atom[g] == m->defs[i].label_atom)) g++;
        if (g == ngroups) {
            group_kind = realloc(group_kind, sizeof(int) * (ngroups + 1));
            group_atom = realloc(group_atom, sizeof(int) * (ngroups + 1));
            group_kind[g] = m->defs[i].kind;
            group_atom[g] = m->defs[i].label_atom;
            ngroups++;
        }
        def_group[i] = g;
    }
    int status = SAT;
    for (int g = 0; g < ngroups && status == SAT; g++)
        for (int i = 0; i < m->ndefs; i++) {
            if (def_group[i] != g) continue;
            const def_t *d = &m->defs[i];
            int nonzero = 0;
            for (int r = 0; r < 3 * d->nrows && !nonzero; r++) nonzero = d->vec[r] != 0;
            if (d->nf_zero || !nonzero) {
                status = REF;
                break;
            }
        }
    for (int g = 0; g < ngroups && status == SAT; g++) {
        int ncols = 0, ncoords = 0, total = 0;
        for (int i = 0; i < m->ndefs; i++)
            if (def_group[i] == g) total += m->defs[i].nrows;
        column_t *cols = malloc(sizeof(column_t) * m->ndefs);
        hmap coords, seen;
        hm_init(&coords, 3 * total);
        for (int i = 0; i < m->ndefs && status == SAT; i++) {
            if (def_group[i] != g) continue;
            const def_t *d = &m->defs[i];
            const block_t *b = &m->blocks[d->block];
            column_t *col = &cols[ncols++];
            col->n = 0;
            col->coord = malloc(sizeof(int) * 3 * d->nrows);
            col->val = malloc(sizeof(int64_t) * 3 * d->nrows);
            hm_init(&seen, d->nrows);
            for (int r = 0; r < d->nrows; r++) {
                if (hm_put(&seen, pkey(0, b->pos + 3 * r), r) >= 0) {
                    status = fail(m, "mode definition repeats a periodic position");
                    break;
                }
                for (int a = 0; a < 3; a++) {
                    if (!d->vec[3 * r + a]) continue;
                    int64_t key = pkey(0, b->pos + 3 * r) * 3 + a;
                    int idx = hm_get(&coords, key);
                    if (idx < 0) hm_put(&coords, key, idx = ncoords++);
                    col->coord[col->n] = idx;
                    col->val[col->n++] = d->vec[3 * r + a];
                }
            }
            hm_free(&seen);
        }
        if (status == SAT && !full_rank(cols, ncols, ncoords)) status = REF;
        for (int j = 0; j < ncols; j++) {
            free(cols[j].coord);
            free(cols[j].val);
        }
        free(cols);
        hm_free(&coords);
    }
    if (status == SAT && m->nstrain) {
        column_t *cols = malloc(sizeof(column_t) * m->nstrain);
        for (int j = 0; j < m->nstrain; j++) {
            cols[j].n = 0;
            cols[j].coord = malloc(sizeof(int) * 6);
            cols[j].val = malloc(sizeof(int64_t) * 6);
            for (int a = 0; a < 6; a++)
                if (m->strain[j][a] && !m->strain_zero[j]) {
                    cols[j].coord[cols[j].n] = a;
                    cols[j].val[cols[j].n++] = m->strain[j][a];
                }
        }
        if (!full_rank(cols, m->nstrain, 6)) status = REF;
        for (int j = 0; j < m->nstrain; j++) {
            free(cols[j].coord);
            free(cols[j].val);
        }
        free(cols);
    }
    free(group_kind);
    free(group_atom);
    free(def_group);
    return m->failed ? -1 : status;
}

/* ---- invariance and fixed space ---- */

typedef struct {
    int n, bijective;
    int *start, *target, *match;
} graph_t;

static int64_t floordiv(int64_t a, int64_t b) { return a / b - ((a % b) && ((a < 0) != (b < 0))); }

static int cmp_int(const void *a, const void *b) { return (*(const int *)a > *(const int *)b) - (*(const int *)a < *(const int *)b); }

/* the sorted residues of the integers within the rounding cell of one transformed coordinate */
static int axis_candidates(const op_t *o, int axis, const int *p, int *out)
{
    int64_t td = o->td[axis], sum = 0, norm = 1;
    for (int s = 0; s < 3; s++) {
        sum += (int64_t)o->r[3 * axis + s] * p[s];
        norm += llabs(o->r[3 * axis + s]);
    }
    int64_t c2 = 2 * td * sum + 2 * (int64_t)PSCALE * o->tn[axis], b2 = td * norm;
    int64_t lo = -floordiv(-(c2 - b2), 2 * td), hi = floordiv(c2 + b2, 2 * td);
    int n = 0;
    for (int64_t v = lo; v <= hi && n < 16; v++) {
        int64_t r = v % PSCALE;
        out[n++] = (int)(r < 0 ? r + PSCALE : r);
    }
    qsort(out, n, sizeof(int), cmp_int);
    int u = 0;
    for (int i = 0; i < n; i++)
        if (!u || out[u - 1] != out[i]) out[u++] = out[i];
    return u;
}

static int cmp_i64(const void *a, const void *b) { return (*(const int64_t *)a > *(const int64_t *)b) - (*(const int64_t *)a < *(const int64_t *)b); }

/* _maximum_matching, in the same order */
static void max_matching(graph_t *g)
{
    int n = g->n;
    int *ts = malloc(sizeof(int) * n), *order = malloc(sizeof(int) * n), *queue = malloc(sizeof(int) * n);
    int *seen_s = calloc(n, sizeof(int)), *seen_t = calloc(n, sizeof(int)), *pred_s = malloc(sizeof(int) * n), *pred_t = malloc(sizeof(int) * n);
    int stamp = 0;
    for (int i = 0; i < n; i++) {
        g->match[i] = ts[i] = -1;
        order[i] = i;
    }
    int64_t *keys = malloc(sizeof(int64_t) * (n + 1));
    for (int i = 0; i < n; i++) keys[i] = (int64_t)(g->start[i + 1] - g->start[i]) << 32 | i;
    qsort(keys, n, sizeof(int64_t), cmp_i64);
    for (int i = 0; i < n; i++) order[i] = (int)(keys[i] & 0xffffffff);
    free(keys);
    for (int k = 0; k < n; k++) {
        int root = order[k];
        if (g->match[root] >= 0) continue;
        stamp++;
        int head = 0, tail = 0, end_s = -1, end_t = -1;
        queue[tail++] = root;
        seen_s[root] = stamp;
        while (head < tail && end_s < 0) {
            int s = queue[head++];
            for (int e = g->start[s]; e < g->start[s + 1]; e++) {
                int t = g->target[e];
                if (seen_t[t] == stamp) continue;
                seen_t[t] = stamp;
                int owner = ts[t];
                if (owner < 0) {
                    end_s = s;
                    end_t = t;
                    break;
                }
                if (seen_s[owner] != stamp) {
                    seen_s[owner] = stamp;
                    pred_s[owner] = s;
                    pred_t[owner] = t;
                    queue[tail++] = owner;
                }
            }
        }
        if (end_s < 0) continue;
        int s = end_s, t = end_t;
        for (;;) {
            g->match[s] = t;
            ts[t] = s;
            if (s == root) break;
            int ps = pred_s[s], pt = pred_t[s];
            s = ps;
            t = pt;
        }
    }
    free(ts);
    free(order);
    free(queue);
    free(seen_s);
    free(seen_t);
    free(pred_s);
    free(pred_t);
}

static void position_graph(const int *atom, const int *pos, int n, const op_t *o, graph_t *g)
{
    hmap h;
    hm_init(&h, n);
    for (int i = 0; i < n; i++) hm_put(&h, pkey(atom[i], pos + 3 * i), i);
    int cap = 2 * n + 4, count = 0;
    g->n = n;
    g->start = malloc(sizeof(int) * (n + 1));
    g->target = malloc(sizeof(int) * cap);
    g->match = malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) {
        int ax[3][16], na[3];
        for (int a = 0; a < 3; a++) na[a] = axis_candidates(o, a, pos + 3 * i, ax[a]);
        g->start[i] = count;
        for (int x = 0; x < na[0]; x++)
            for (int y = 0; y < na[1]; y++)
                for (int z = 0; z < na[2]; z++) {
                    int q[3] = {ax[0][x], ax[1][y], ax[2][z]}, idx = hm_get(&h, pkey(atom[i], q));
                    if (idx < 0) continue;
                    if (count == cap) g->target = realloc(g->target, sizeof(int) * (cap *= 2));
                    g->target[count++] = idx;
                }
    }
    g->start[n] = count;
    hm_free(&h);
    g->bijective = count == n;
    if (g->bijective) {
        char *hit = calloc(n, 1);
        for (int i = 0; i < n && g->bijective; i++) {
            if (g->start[i + 1] - g->start[i] != 1 || hit[g->target[g->start[i]]]) g->bijective = 0;
            else hit[g->target[g->start[i]]] = 1;
        }
        free(hit);
    }
    if (g->bijective)
        for (int i = 0; i < n; i++) g->match[i] = g->target[g->start[i]];
    else max_matching(g);
}

static void graph_free(graph_t *g)
{
    free(g->start);
    free(g->target);
    free(g->match);
}

static int perfect(const graph_t *g)
{
    for (int i = 0; i < g->n; i++)
        if (g->match[i] < 0) return 0;
    return 1;
}

static int det3(const int *r) { return r[0] * (r[4] * r[8] - r[5] * r[7]) - r[1] * (r[3] * r[8] - r[5] * r[6]) + r[2] * (r[3] * r[7] - r[4] * r[6]); }

static void vector_action(const op_t *o, int kind, int *t)
{
    int f = kind == 1 ? o->tr * det3(o->r) : 1;
    for (int i = 0; i < 9; i++) t[i] = kind == 2 ? (i % 4 == 0) : f * o->r[i];
}

static int64_t lcm64(int64_t a, int64_t b)
{
    int64_t x = a, y = b;
    while (y) {
        int64_t t = x % y;
        x = y;
        y = t;
    }
    return a / x * b;
}

static int representative_count(const orbit_t *orb, const op_t *ops, int nops)
{
    int n = nops + 1;
    int64_t d = 1, mod;
    for (int k = 0; k < nops; k++)
        for (int a = 0; a < 3; a++) d = lcm64(d, ops[k].td[a]);
    mod = (int64_t)PSCALE * d;
    int64_t (*center)[3] = malloc(sizeof(int64_t[3]) * n);
    const int ident[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    int *parent = malloc(sizeof(int) * n);
    for (int k = 0; k < n; k++) {
        const int *r = k ? ops[k - 1].r : ident;
        parent[k] = k;
        for (int a = 0; a < 3; a++) {
            int64_t s = 0;
            for (int b = 0; b < 3; b++) s += (int64_t)r[3 * a + b] * orb->rep[b];
            int64_t t = k ? (int64_t)PSCALE * ops[k - 1].tn[a] * (d / ops[k - 1].td[a]) : 0;
            int64_t v = (d * s + t) % mod;
            center[k][a] = v < 0 ? v + mod : v;
        }
    }
    for (int l = 0; l < n; l++)
        for (int r = 0; r < l; r++) {
            const int *rl = l ? ops[l - 1].r : ident, *rr = r ? ops[r - 1].r : ident;
            int close = 1;
            for (int a = 0; a < 3 && close; a++) {
                int64_t diff = ((center[l][a] - center[r][a]) % mod + mod) % mod, dist = diff < mod - diff ? diff : mod - diff, bound = 0;
                for (int b = 0; b < 3; b++) bound += llabs(rl[3 * a + b] - rr[3 * a + b]);
                close = 2 * dist <= d * bound;
            }
            if (close) {
                int a = uf_find(parent, l), b = uf_find(parent, r);
                if (a != b) parent[b] = a;
            }
        }
    int count = 0;
    for (int k = 0; k < n; k++) count += uf_find(parent, k) == k;
    free(center);
    free(parent);
    return count;
}

static int structure_ok(const mc_t *m, const op_t *ops, int nops)
{
    for (int i = 0; i < m->norbits; i++) {
        const orbit_t *o = &m->orbits[i];
        int observed;
        if (o->nexp) {
            int *atom = malloc(sizeof(int) * o->nexp);
            for (int j = 0; j < o->nexp; j++) atom[j] = o->label;
            int ok = 1;
            for (int k = 0; k < nops && ok; k++) {
                graph_t g;
                position_graph(atom, o->exp, o->nexp, &ops[k], &g);
                ok = perfect(&g);
                graph_free(&g);
            }
            free(atom);
            if (!ok) return 0;
            observed = o->nexp;
        } else observed = representative_count(o, ops, nops);
        if (observed != o->multiplicity) return 0;
    }
    return 1;
}

static int read_ops(const int64_t *w, op_t **out)
{
    int n = (int)w[0];
    op_t *ops = malloc(sizeof(op_t) * (n + 1));
    for (int k = 0; k < n; k++) {
        const int64_t *x = w + 1 + 16 * k;
        for (int i = 0; i < 9; i++) ops[k].r[i] = (int)x[i];
        for (int a = 0; a < 3; a++) {
            ops[k].tn[a] = x[9 + a];
            ops[k].td[a] = x[12 + a];
        }
        ops[k].tr = (int)x[15];
    }
    *out = ops;
    return n;
}

/* graphs[b * nops + k] for every block used by an atomic definition */
static graph_t *block_graphs(const mc_t *m, const op_t *ops, int nops)
{
    graph_t *gs = calloc((size_t)m->nblocks * (nops ? nops : 1), sizeof(graph_t));
    char *used = calloc(m->nblocks + 1, 1);
    for (int i = 0; i < m->ndefs; i++) used[m->defs[i].block] = 1;
    for (int b = 0; b < m->nblocks; b++)
        if (used[b])
            for (int k = 0; k < nops; k++) position_graph(m->blocks[b].atom, m->blocks[b].pos, m->blocks[b].n, &ops[k], &gs[(int64_t)b * nops + k]);
    free(used);
    return gs;
}

static void block_graphs_free(const mc_t *m, graph_t *gs, int nops)
{
    for (int64_t i = 0; i < (int64_t)m->nblocks * nops; i++)
        if (gs[i].start) graph_free(&gs[i]);
    free(gs);
}

static int vector_ok(const def_t *d, const block_t *b, const graph_t *g, const op_t *o)
{
    int t[9], bound2[3];
    vector_action(o, d->kind, t);
    for (int a = 0; a < 3; a++) bound2[a] = 1 + abs(t[3 * a]) + abs(t[3 * a + 1]) + abs(t[3 * a + 2]);
    graph_t v = {g->n, 0, malloc(sizeof(int) * (g->n + 1)), malloc(sizeof(int) * (g->start[g->n] + 1)), malloc(sizeof(int) * g->n)};
    int count = 0;
    for (int s = 0; s < g->n; s++) {
        int64_t e[3];
        for (int a = 0; a < 3; a++)
            e[a] = (int64_t)t[3 * a] * d->vec[3 * s] + (int64_t)t[3 * a + 1] * d->vec[3 * s + 1] + (int64_t)t[3 * a + 2] * d->vec[3 * s + 2];
        v.start[s] = count;
        for (int x = g->start[s]; x < g->start[s + 1]; x++) {
            int tgt = g->target[x], ok = 1;
            for (int a = 0; a < 3 && ok; a++) ok = 2 * llabs(d->vec[3 * tgt + a] - e[a]) <= bound2[a];
            if (ok) v.target[count++] = tgt;
        }
    }
    v.start[g->n] = count;
    int result;
    if (g->bijective) result = count == g->n;
    else {
        max_matching(&v);
        result = perfect(&v);
    }
    (void)b;
    free(v.start);
    free(v.target);
    free(v.match);
    return result;
}

/* settings: [nsettings, per setting nops then 16 words per operation]; returns the index of the first setting where
 * every atomic definition is invariant, -1 when none, -2 when there are no atomic definitions, -3 on error */
int mc_invariance(void *h, const int64_t *settings)
{
    mc_t *m = h;
    if (m->failed) return -3;
    if (!m->ndefs) return -2;
    int ns = (int)settings[0];
    const int64_t *w = settings + 1;
    for (int s = 0; s < ns; s++) {
        op_t *ops;
        int nops = read_ops(w, &ops);
        w += 1 + 16 * nops;
        int ok = structure_ok(m, ops, nops);
        graph_t *gs = ok ? block_graphs(m, ops, nops) : 0;
        for (int i = 0; ok && i < m->ndefs; i++)
            for (int k = 0; k < nops && ok; k++) ok = perfect(&gs[(int64_t)m->defs[i].block * nops + k]);
        for (int i = 0; ok && i < m->ndefs; i++)
            for (int k = 0; k < nops && ok; k++)
                ok = vector_ok(&m->defs[i], &m->blocks[m->defs[i].block], &gs[(int64_t)m->defs[i].block * nops + k], &ops[k]);
        if (gs) block_graphs_free(m, gs, nops);
        free(ops);
        if (ok) return s;
    }
    return -1;
}

/* the first setting whose operations permute every block's atoms: its index, and in *refuted whether some block's
 * definition count differs from its fixed-space dimension; -1 when none, -2 without atomic definitions */
int mc_fixed(void *h, const int64_t *settings, int *refuted)
{
    mc_t *m = h;
    if (m->failed) return -3;
    if (!m->ndefs) return -2;
    int nb = 0, *bkind = 0, *bblock = 0, *bcount = 0;
    for (int i = 0; i < m->ndefs; i++) {
        int b = 0;
        while (b < nb && !(bkind[b] == m->defs[i].kind && bblock[b] == m->defs[i].block)) b++;
        if (b == nb) {
            bkind = realloc(bkind, sizeof(int) * (nb + 1));
            bblock = realloc(bblock, sizeof(int) * (nb + 1));
            bcount = realloc(bcount, sizeof(int) * (nb + 1));
            bkind[b] = m->defs[i].kind;
            bblock[b] = m->defs[i].block;
            bcount[b] = 0;
            nb++;
        }
        bcount[b]++;
    }
    int ns = (int)settings[0], found = -1;
    const int64_t *w = settings + 1;
    for (int s = 0; s < ns && found < 0; s++) {
        op_t *ops;
        int nops = read_ops(w, &ops);
        w += 1 + 16 * nops;
        if (structure_ok(m, ops, nops)) {
            graph_t *gs = block_graphs(m, ops, nops);
            int ok = 1, bad = 0;
            for (int b = 0; b < nb && ok; b++) {
                const block_t *blk = &m->blocks[bblock[b]];
                int64_t total = bkind[b] == 2 ? blk->n : 3 * blk->n;
                for (int k = 0; k < nops && ok; k++) {
                    const graph_t *g = &gs[(int64_t)bblock[b] * nops + k];
                    if (!(ok = perfect(g))) break;
                    int t[9], trace, fixed = 0;
                    vector_action(&ops[k], bkind[b], t);
                    trace = bkind[b] == 2 ? 1 : t[0] + t[4] + t[8];
                    for (int i = 0; i < g->n; i++) fixed += g->match[i] == i;
                    total += (int64_t)trace * fixed;
                }
                if (ok && total != (int64_t)bcount[b] * (nops + 1)) bad = 1;
            }
            block_graphs_free(m, gs, nops);
            if (ok) {
                found = s;
                *refuted = bad;
            }
        }
        free(ops);
    }
    free(bkind);
    free(bblock);
    free(bcount);
    return found;
}

/* ---- the k-star sector and the isotropy index (no tables) ----
 * The parent structure's operations come from the child cell's atoms: every child-lattice translate is filled into one
 * parent cell, and an operation is a metric-preserving (x -> x Q + t, Q with entries -1..1) map of that atom set onto
 * itself, atoms matched by parent site and position. */

#define FTOL 2e-3

typedef struct {
    int natoms, nkinds, nparent, nrot, ncent, nreps, ops_n;
    int *kind;                 /* per child atom: its parent site */
    double *y;                 /* child coordinates */
    double *px;                /* parent-cell atoms */
    int *pkind;
    double B[9], Binv[9], M[9];
    double (*rot)[9];          /* point rotations, row form */
    double (*cent)[3];         /* pure translations of the parent structure */
    double (*reps)[3];         /* parent lattice translations modulo the child lattice, child coordinates */
    double (*ops)[12];         /* one operation (Q, t) per rotation */
    int *perm;                 /* nreps x natoms: atom i moved by rep r */
    hmap cell;                 /* (kind, position) -> child atom */
    strtab kinds;              /* parent site names */
} sg_t;

static double wrap1(double x)
{
    double w = x - floor(x);
    return w >= 1 - 1e-9 ? 0 : w;
}

static int near_mod1(const double *a, const double *b, double tol)
{
    for (int i = 0; i < 3; i++) {
        double d = fabs(wrap1(a[i] - b[i]));
        if (fmin(d, 1 - d) > tol) return 0;
    }
    return 1;
}

static void row_mul(const double *x, const double *A, double *out)
{
    for (int j = 0; j < 3; j++) out[j] = x[0] * A[j] + x[1] * A[3 + j] + x[2] * A[6 + j];
}

static void mat_mul(const double *A, const double *Bm, double *out)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) out[3 * i + j] = A[3 * i] * Bm[j] + A[3 * i + 1] * Bm[3 + j] + A[3 * i + 2] * Bm[6 + j];
}

static int inverse3(const double *A, double *out)
{
    double d = A[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (A[3] * A[8] - A[5] * A[6]) + A[2] * (A[3] * A[7] - A[4] * A[6]);
    if (fabs(d) < 1e-12) return -1;
    out[0] = (A[4] * A[8] - A[5] * A[7]) / d, out[1] = (A[2] * A[7] - A[1] * A[8]) / d, out[2] = (A[1] * A[5] - A[2] * A[4]) / d;
    out[3] = (A[5] * A[6] - A[3] * A[8]) / d, out[4] = (A[0] * A[8] - A[2] * A[6]) / d, out[5] = (A[2] * A[3] - A[0] * A[5]) / d;
    out[6] = (A[3] * A[7] - A[4] * A[6]) / d, out[7] = (A[1] * A[6] - A[0] * A[7]) / d, out[8] = (A[0] * A[4] - A[1] * A[3]) / d;
    return 0;
}

static int64_t cell_key(int kind, const double *y, const int *shift)
{
    int64_t k = kind;
    for (int a = 0; a < 3; a++) k = k * 1000 + (((int64_t)floor(wrap1(y[a]) * 1000) + shift[a] + 1000) % 1000);
    return k;
}

/* the child atom of this parent site at y (mod 1), or -1 */
static int atom_at(const sg_t *g, int kind, const double *y)
{
    static const int none[3] = {0, 0, 0};
    int i = hm_get(&g->cell, cell_key(kind, y, none));
    if (i >= 0 && near_mod1(g->y + 3 * i, y, FTOL)) return i;
    for (int s = 0; s < 27; s++) {
        int sh[3] = {s % 3 - 1, s / 3 % 3 - 1, s / 9 - 1};
        i = hm_get(&g->cell, cell_key(kind, y, sh));
        if (i >= 0 && near_mod1(g->y + 3 * i, y, FTOL)) return i;
    }
    return -1;
}

static int parent_atom_index(const sg_t *g, int kind, const double *x)
{
    for (int j = 0; j < g->nparent; j++)
        if (g->pkind[j] == kind && near_mod1(g->px + 3 * j, x, FTOL)) return j;
    return -1;
}

static int maps_parent(const sg_t *g, const double *Q, const double *t)
{
    for (int j = 0; j < g->nparent; j++) {
        double x[3];
        row_mul(g->px + 3 * j, Q, x);
        for (int a = 0; a < 3; a++) x[a] += t[a];
        if (parent_atom_index(g, g->pkind[j], x) < 0) return 0;
    }
    return 1;
}

static void sg_free(sg_t *g)
{
    free(g->kind), free(g->y), free(g->px), free(g->pkind), free(g->rot), free(g->cent), free(g->reps), free(g->ops), free(g->perm);
    if (g->cell.keys) hm_free(&g->cell);
    for (int i = 0; i < g->kinds.n; i++) free(g->kinds.s[i]);
    free(g->kinds.s), free(g->kinds.len);
}

/* the parent site of a child site marker ("Rb0_1" -> "Rb0"), interned in t */
static int site_kind(strtab *t, const strtab *markers, int marker)
{
    const char *s = markers->s[marker];
    int len = markers->len[marker], cut = len;
    while (cut > 0 && s[cut - 1] >= '0' && s[cut - 1] <= '9') cut--;
    if (cut > 0 && cut < len && s[cut - 1] == '_') len = cut - 1;
    return intern(t, s, len);
}

static void cell_metric(const double *c, double *G)
{
    double ca = cos(c[3] * M_PI / 180), cb = cos(c[4] * M_PI / 180), cg = cos(c[5] * M_PI / 180);
    G[0] = c[0] * c[0], G[4] = c[1] * c[1], G[8] = c[2] * c[2];
    G[1] = G[3] = c[0] * c[1] * cg, G[2] = G[6] = c[0] * c[2] * cb, G[5] = G[7] = c[1] * c[2] * ca;
}

/* the parent structure group of the payload's child cell; 0, or NA when the payload lacks the frame */
static int structure(const mc_t *m, sg_t *g)
{
    memset(g, 0, sizeof *g);
    if (!m->has_basis || !m->has_origin || !m->has_cell) return NA;
    memcpy(g->B, m->basis, sizeof g->B);
    if (inverse3(g->B, g->Binv)) return NA;
    double G[9], T1[9], BinvT[9];
    cell_metric(m->cell, G);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) BinvT[3 * i + j] = g->Binv[3 * j + i];
    mat_mul(g->Binv, G, T1);
    mat_mul(T1, BinvT, g->M);
    /* child atoms: the atoms list, or the definition rows of an older payload */
    int cap = 64;
    g->kind = malloc(sizeof(int) * cap);
    g->y = malloc(sizeof(double) * 3 * cap);
    hm_init(&g->cell, m->natoms + 1024);
    for (int b = -1, carried = -1; b < m->nblocks; b++) {
        int n = b < 0 ? m->natoms : m->blocks[b].n;
        if (b >= 0 && m->natoms) break;
        for (int r = 0; r < n; r++) {
            int marker = b < 0 ? (m->atom_label[r] >= 0 ? m->atom_label[r] : m->atom_marker[r]) : m->blocks[b].atom[r];
            const int *q = b < 0 ? m->atom_pos + 3 * r : m->blocks[b].pos + 3 * r;
            if (marker < 0) marker = carried;
            if (marker < 0) continue;
            carried = marker;
            int k = site_kind(&g->kinds, &m->markers, marker);
            double y[3] = {q[0] / (double)PSCALE, q[1] / (double)PSCALE, q[2] / (double)PSCALE};
            if (atom_at(g, k, y) >= 0) continue;
            if (g->natoms == cap) {
                cap *= 2;
                g->kind = realloc(g->kind, sizeof(int) * cap);
                g->y = realloc(g->y, sizeof(double) * 3 * cap);
            }
            if (g->natoms >= g->cell.cap / 2) {
                hmap bigger;
                hm_init(&bigger, 4 * g->natoms);
                static const int none[3] = {0, 0, 0};
                for (int i = 0; i < g->natoms; i++) hm_put(&bigger, cell_key(g->kind[i], g->y + 3 * i, none), i);
                hm_free(&g->cell);
                g->cell = bigger;
            }
            static const int none[3] = {0, 0, 0};
            g->kind[g->natoms] = k;
            memcpy(g->y + 3 * g->natoms, y, sizeof y);
            hm_put(&g->cell, cell_key(k, y, none), g->natoms);
            g->natoms++;
        }
    }
    g->nkinds = g->kinds.n;
    if (!g->natoms) return NA;
    /* one parent cell: child atoms in parent coordinates, closed under the child lattice vectors */
    int pcap = 64, *todo_k = malloc(sizeof(int) * pcap), ntodo = 0;
    double *todo = malloc(sizeof(double) * 3 * pcap);
    g->px = malloc(sizeof(double) * 3 * pcap);
    g->pkind = malloc(sizeof(int) * pcap);
    int pxcap = pcap;
    for (int i = 0; i < g->natoms; i++) {
        if (ntodo == pcap) pcap *= 2, todo = realloc(todo, sizeof(double) * 3 * pcap), todo_k = realloc(todo_k, sizeof(int) * pcap);
        row_mul(g->y + 3 * i, g->B, todo + 3 * ntodo);
        for (int a = 0; a < 3; a++) todo[3 * ntodo + a] = wrap1(todo[3 * ntodo + a] + m->origin[a]);
        todo_k[ntodo++] = g->kind[i];
    }
    while (ntodo) {
        ntodo--;
        double x[3];
        memcpy(x, todo + 3 * ntodo, sizeof x);
        int k = todo_k[ntodo];
        if (parent_atom_index(g, k, x) >= 0) continue;
        if (g->nparent == pxcap) pxcap *= 2, g->px = realloc(g->px, sizeof(double) * 3 * pxcap), g->pkind = realloc(g->pkind, sizeof(int) * pxcap);
        memcpy(g->px + 3 * g->nparent, x, sizeof x);
        g->pkind[g->nparent++] = k;
        if (g->nparent > 20000) break;
        for (int r = 0; r < 3; r++) {
            if (ntodo == pcap) pcap *= 2, todo = realloc(todo, sizeof(double) * 3 * pcap), todo_k = realloc(todo_k, sizeof(int) * pcap);
            for (int a = 0; a < 3; a++) todo[3 * ntodo + a] = wrap1(x[a] + g->B[3 * r + a]);
            todo_k[ntodo++] = k;
        }
    }
    free(todo), free(todo_k);
    /* rotations and pure translations */
    int *count = calloc(g->nkinds, sizeof(int)), rare = 0;
    for (int j = 0; j < g->nparent; j++) count[g->pkind[j]]++;
    for (int k = 0; k < g->nkinds; k++)
        if (count[k] && (!count[rare] || count[k] < count[rare])) rare = k;
    free(count);
    int p0 = -1;
    for (int j = 0; j < g->nparent && p0 < 0; j++)
        if (g->pkind[j] == rare) p0 = j;
    g->rot = malloc(sizeof(double[9]) * 48);
    g->ops = malloc(sizeof(double[12]) * 48);
    g->cent = malloc(sizeof(double[3]) * 64);
    for (int code = 0; code < 19683; code++) {
        double Q[9];
        for (int i = 0, c = code; i < 9; i++, c /= 3) Q[i] = c % 3 - 1;
        double QM[9], QMQ[9], QT[9];
        double det = Q[0] * (Q[4] * Q[8] - Q[5] * Q[7]) - Q[1] * (Q[3] * Q[8] - Q[5] * Q[6]) + Q[2] * (Q[3] * Q[7] - Q[4] * Q[6]);
        if (fabs(det) != 1) continue;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) QT[3 * i + j] = Q[3 * j + i];
        mat_mul(Q, g->M, QM);
        mat_mul(QM, QT, QMQ);
        int metric = 1;
        for (int i = 0; i < 9 && metric; i++) metric = fabs(QMQ[i] - g->M[i]) <= 1e-3 * sqrt(g->M[4 * (i / 3)] * g->M[4 * (i % 3)]);
        if (!metric) continue;
        int identity = Q[0] == 1 && Q[4] == 1 && Q[8] == 1 && !Q[1] && !Q[2] && !Q[3] && !Q[5] && !Q[6] && !Q[7], found = 0;
        for (int j = 0; j < g->nparent; j++) {
            if (g->pkind[j] != rare) continue;
            double moved[3], t[3];
            row_mul(g->px + 3 * p0, Q, moved);
            for (int a = 0; a < 3; a++) t[a] = wrap1(g->px[3 * j + a] - moved[a]);
            if (!maps_parent(g, Q, t)) continue;
            if (!found && g->nrot < 48) {
                memcpy(g->rot[g->nrot], Q, sizeof Q);
                memcpy(g->ops[g->nrot], Q, sizeof Q);
                memcpy(g->ops[g->nrot] + 9, t, sizeof t);
                g->nrot++;
            }
            found = 1;
            if (!identity) break;
            /* a pure translation of a parent lattice is a centring vector, a multiple of 1/6 */
            double snapped[3];
            for (int a = 0; a < 3; a++) snapped[a] = wrap1(round(t[a] * 6) / 6);
            if (g->ncent < 64 && near_mod1(snapped, t, FTOL)) memcpy(g->cent[g->ncent++], snapped, sizeof snapped);
        }
    }
    /* parent lattice translations (unit vectors and the pure translations) modulo the child lattice */
    int rcap = 64;
    g->reps = malloc(sizeof(double[3]) * rcap);
    g->reps[0][0] = g->reps[0][1] = g->reps[0][2] = 0;
    g->nreps = 1;
    for (int i = 0; i < g->nreps; i++)
        for (int gen = 0; gen < 3 + g->ncent; gen++) {
            double v[3] = {gen == 0, gen == 1, gen == 2}, s[3];
            if (gen >= 3) memcpy(v, g->cent[gen - 3], sizeof v);
            row_mul(v, g->Binv, s);
            for (int a = 0; a < 3; a++) s[a] = wrap1(s[a] + g->reps[i][a]);
            int known = 0;
            for (int j = 0; j < g->nreps && !known; j++) known = near_mod1(g->reps[j], s, 1e-6);
            if (known) continue;
            if (g->nreps == 4096) return IND;
            if (g->nreps == rcap) rcap *= 2, g->reps = realloc(g->reps, sizeof(double[3]) * rcap);
            memcpy(g->reps[g->nreps++], s, sizeof s);
        }
    g->perm = malloc(sizeof(int) * (size_t)g->nreps * g->natoms);
    for (int r = 0; r < g->nreps; r++)
        for (int i = 0; i < g->natoms; i++) {
            double y[3] = {g->y[3 * i] + g->reps[r][0], g->y[3 * i + 1] + g->reps[r][1], g->y[3 * i + 2] + g->reps[r][2]};
            if ((g->perm[(size_t)r * g->natoms + i] = atom_at(g, g->kind[i], y)) < 0) return IND;
        }
    return 0;
}

/* the field of one definition on the child atoms (vectors at display precision), zero elsewhere */
static int def_field(const mc_t *m, sg_t *g, const def_t *d, double *u)
{
    memset(u, 0, sizeof(double) * 3 * g->natoms);
    const block_t *b = &m->blocks[d->block];
    int kind = intern(&g->kinds, m->label_atoms.s[d->label_atom], m->label_atoms.len[d->label_atom]);
    for (int r = 0; r < d->nrows; r++) {
        double y[3] = {b->pos[3 * r] / (double)PSCALE, b->pos[3 * r + 1] / (double)PSCALE, b->pos[3 * r + 2] / (double)PSCALE};
        int at = atom_at(g, kind, y);
        if (at < 0) return -1;
        for (int a = 0; a < 3; a++) u[3 * at + a] += d->vec[3 * r + a] / 1e4;
    }
    return 0;
}

/* every mode lies in the Bloch sector of the arms of its label's k star that are periodic in the child cell (the
 * labelled arm itself need not be: a secondary irrep's representative) */
int mc_kstar(void *h)
{
    mc_t *m = h;
    sg_t g;
    int status = structure(m, &g);
    if (status) {
        sg_free(&g);
        return status;
    }
    int N = g.nreps, verdict = NA;
    double *u = malloc(sizeof(double) * 3 * g.natoms), *P = malloc(sizeof(double) * 3 * g.natoms);
    double *phase = malloc(sizeof(double) * 2 * N * 2 * (g.nrot + 1)), *c = malloc(sizeof(double) * N), *tp = malloc(sizeof(double) * 3 * N);
    for (int r = 0; r < N; r++) row_mul(g.reps[r], g.B, tp + 3 * r);
    for (int di = 0; di < m->ndefs && verdict != REF; di++) {
        const def_t *d = &m->defs[di];
        if (!d->has_k) continue;
        if (def_field(m, &g, d, u)) {
            verdict = IND;
            continue;
        }
        double norm = 0;
        for (int i = 0; i < 3 * g.natoms; i++) norm += u[i] * u[i];
        if (norm == 0) continue;
        int arms = 0;
        for (int q = 0; q < g.nrot; q++) {
            double Qinv[9];
            if (inverse3(g.rot[q], Qinv)) continue;
            for (int sign = 1; sign >= -1; sign -= 2) {
                double kk[3];
                for (int j = 0; j < 3; j++) kk[j] = sign * (d->k[0] * Qinv[3 * j] + d->k[1] * Qinv[3 * j + 1] + d->k[2] * Qinv[3 * j + 2]);
                int periodic = 1;
                for (int r = 0; r < 3 && periodic; r++) {
                    double v = g.B[3 * r] * kk[0] + g.B[3 * r + 1] * kk[1] + g.B[3 * r + 2] * kk[2];
                    periodic = fabs(v - round(v)) <= 1e-6;
                }
                if (!periodic) continue;
                double *ph = phase + 2 * N * arms;
                for (int r = 0; r < N; r++) {
                    double a = 2 * M_PI * (kk[0] * tp[3 * r] + kk[1] * tp[3 * r + 1] + kk[2] * tp[3 * r + 2]);
                    ph[2 * r] = cos(a), ph[2 * r + 1] = sin(a);
                }
                int known = 0;
                for (int e = 0; e < arms && !known; e++) {
                    known = 1;
                    for (int r = 0; r < 2 * N && known; r++) known = fabs(phase[2 * N * e + r] - ph[r]) < 1e-6;
                }
                if (!known) arms++;
            }
        }
        if (!arms) {
            verdict = REF;
            break;
        }
        for (int r = 0; r < N; r++) {
            c[r] = 0;
            for (int e = 0; e < arms; e++) c[r] += phase[2 * N * e + 2 * r];
        }
        memset(P, 0, sizeof(double) * 3 * g.natoms);
        for (int r = 0; r < N; r++) {
            if (fabs(c[r]) < 1e-12) continue;
            const int *perm = g.perm + (size_t)r * g.natoms;
            for (int i = 0; i < g.natoms; i++)
                for (int a = 0; a < 3; a++) P[3 * perm[i] + a] += c[r] * u[3 * i + a];
        }
        double res = 0;
        for (int i = 0; i < 3 * g.natoms; i++) res += (P[i] / N - u[i]) * (P[i] / N - u[i]);
        if (verdict == NA) verdict = SAT;
        if (sqrt(res / norm) > 1e-3) verdict = REF;
    }
    free(u), free(P), free(phase), free(c), free(tp);
    sg_free(&g);
    return verdict;
}

/* the stabilizer, in the parent structure group (with time reversal when a primary mode is magnetic), of a generic
 * combination of the primary modes has the declared index: satisfied; a larger stabilizer (the field does not resolve
 * the order parameter) is indeterminate; a smaller one contradicts the declared group, which leaves the field fixed */
int mc_isotropy(void *h)
{
    mc_t *m = h;
    if (m->index <= 0) return NA;
    sg_t g;
    int status = structure(m, &g);
    if (status) {
        sg_free(&g);
        return status;
    }
    double *disp = calloc(3 * g.natoms, sizeof(double)), *mag = calloc(3 * g.natoms, sizeof(double)), *u = malloc(sizeof(double) * 3 * g.natoms);
    uint64_t seed = 0x9e3779b97f4a7c15ULL;
    int primaries = 0, magnetic = 0, verdict = NA;
    for (int di = 0; di < m->ndefs; di++) {
        const def_t *d = &m->defs[di];
        if (!d->primary || d->kind > 1) continue;
        if (def_field(m, &g, d, u)) {
            verdict = IND;
            goto done;
        }
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        double coefficient = 0.5 + (double)(seed >> 11) / 9007199254740992.0;
        double *f = d->kind ? mag : disp;
        for (int i = 0; i < 3 * g.natoms; i++) f[i] += coefficient * u[i];
        primaries++;
        magnetic |= d->kind == 1;
    }
    if (!primaries) goto done;
    double scale = 0;
    for (int i = 0; i < 3 * g.natoms; i++) scale = fmax(scale, fmax(fabs(disp[i]), fabs(mag[i])));
    if (scale == 0) goto done;
    int64_t total = 0, stab = 0;
    for (int q = 0; q < g.nrot; q++) {
        const double *Q = g.ops[q], *t = g.ops[q] + 9;
        double Qc[9], T1[9];
        mat_mul(g.B, Q, T1);
        mat_mul(T1, g.Binv, Qc);
        double det = Q[0] * (Q[4] * Q[8] - Q[5] * Q[7]) - Q[1] * (Q[3] * Q[8] - Q[5] * Q[6]) + Q[2] * (Q[3] * Q[7] - Q[4] * Q[6]);
        /* an operation that does not keep the child lattice cannot fix a field of that period (the declared group keeps it) */
        int keeps = 1;
        for (int e = 0; e < 9 && keeps; e++) keeps = fabs(Qc[e] - round(Qc[e])) <= 1e-6;
        for (int r = 0; r < g.nreps; r++) {
            double tau[3];
            row_mul(g.reps[r], g.B, tau);
            for (int theta = 0; theta <= magnetic; theta++) {
                total++;
                int ok = keeps;
                for (int i = 0; i < g.natoms && ok; i++) {
                    double x[3], xq[3], y[3];
                    row_mul(g.y + 3 * i, g.B, x);
                    for (int a = 0; a < 3; a++) x[a] += m->origin[a];
                    row_mul(x, Q, xq);
                    for (int a = 0; a < 3; a++) xq[a] += t[a] + tau[a] - m->origin[a];
                    row_mul(xq, g.Binv, y);
                    int j = atom_at(&g, g.kind[i], y);
                    if (j < 0) {
                        ok = 0;
                        break;
                    }
                    double dv[3], mv[3];
                    row_mul(disp + 3 * i, Qc, dv);
                    row_mul(mag + 3 * i, Qc, mv);
                    double sm = (det > 0 ? 1 : -1) * (theta ? -1 : 1);
                    for (int a = 0; a < 3 && ok; a++)
                        ok = fabs(dv[a] - disp[3 * j + a]) <= 1e-3 * scale && fabs(sm * mv[a] - mag[3 * j + a]) <= 1e-3 * scale;
                }
                stab += ok;
            }
        }
    }
    verdict = !stab ? IND : total == stab * m->index ? SAT : total < stab * m->index ? IND : REF;
done:
    free(disp), free(mag), free(u);
    sg_free(&g);
    return verdict;
}

/* (g w)(x_j) = R w(g^-1 x_j) on a Bloch field of k (w(x + tau) = exp(-2 pi i k.tau) w(x) for every parent translation
 * tau), reduced to the child atoms: g^-1 x_j = x_i + tau with x_i a child atom. 0, or -1 when a preimage is missing. */
static int bloch_act(const mc_t *m, const sg_t *g, const double *xs, int q, const double *k, int axial, const double *w, double *out)
{
    const double *Q = g->ops[q], *t = g->ops[q] + 9;
    double Qinv[9], Qc[9], T1[9];
    if (inverse3(Q, Qinv)) return -1;
    mat_mul(g->B, Q, T1);
    mat_mul(T1, g->Binv, Qc);
    double sign = 1;
    if (axial) sign = Q[0] * (Q[4] * Q[8] - Q[5] * Q[7]) - Q[1] * (Q[3] * Q[8] - Q[5] * Q[6]) + Q[2] * (Q[3] * Q[7] - Q[4] * Q[6]) > 0 ? 1 : -1;
    for (int j = 0; j < g->natoms; j++) {
        double d[3], x[3];
        for (int a = 0; a < 3; a++) d[a] = xs[3 * j + a] - t[a];
        row_mul(d, Qinv, x);
        int found = -1;
        double tau[3] = {0};
        for (int c = -1; c < g->ncent && found < 0; c++) {
            double y[3], shifted[3];
            for (int a = 0; a < 3; a++) shifted[a] = x[a] - m->origin[a] - (c >= 0 ? g->cent[c][a] : 0);
            row_mul(shifted, g->Binv, y);
            int i = atom_at(g, g->kind[j], y);
            if (i < 0) continue;
            found = i;
            for (int a = 0; a < 3; a++) tau[a] = x[a] - xs[3 * i + a];
        }
        if (found < 0) return -1;
        double phase = -2 * M_PI * (k[0] * tau[0] + k[1] * tau[1] + k[2] * tau[2]), c = cos(phase), s = sin(phase);
        double re[3], im[3];
        double wr[3] = {w[6 * found], w[6 * found + 2], w[6 * found + 4]}, wi[3] = {w[6 * found + 1], w[6 * found + 3], w[6 * found + 5]};
        row_mul(wr, Qc, re);
        row_mul(wi, Qc, im);
        for (int a = 0; a < 3; a++) {
            out[6 * j + 2 * a] = sign * (c * re[a] - s * im[a]);
            out[6 * j + 2 * a + 1] = sign * (s * re[a] + c * im[a]);
        }
    }
    return 0;
}

/* within each k star, the little-group images of the k-Bloch part of every mode of one irrep label are orthogonal to the
 * k-Bloch parts of the modes of every other label (distinct irreps have orthogonal isotypic components under the
 * invariant, Cartesian, inner product) */
int mc_isotypic(void *h)
{
    mc_t *m = h;
    sg_t g;
    int status = structure(m, &g);
    if (status) {
        sg_free(&g);
        return status;
    }
    int N = g.nreps, n6 = 6 * g.natoms, verdict = NA;
    double G[9], L[9] = {0};
    cell_metric(m->cell, G);
    L[0] = sqrt(G[0]);
    L[3] = G[3] / L[0], L[4] = sqrt(G[4] - L[3] * L[3]);
    L[6] = G[6] / L[0], L[7] = (G[7] - L[6] * L[3]) / L[4], L[8] = sqrt(G[8] - L[6] * L[6] - L[7] * L[7]);
    double *xs = malloc(sizeof(double) * 3 * g.natoms), *tp = malloc(sizeof(double) * 3 * N), *u = malloc(sizeof(double) * 3 * g.natoms);
    for (int i = 0; i < g.natoms; i++) {
        row_mul(g.y + 3 * i, g.B, xs + 3 * i);
        for (int a = 0; a < 3; a++) xs[3 * i + a] += m->origin[a];
    }
    for (int r = 0; r < N; r++) row_mul(g.reps[r], g.B, tp + 3 * r);
    int *kind_count = calloc(g.nkinds + 1, sizeof(int)), *kind_start = calloc(g.nkinds + 2, sizeof(int)), *kind_atoms = malloc(sizeof(int) * (g.natoms + 1));
    for (int i = 0; i < g.natoms; i++) kind_count[g.kind[i]]++;
    for (int s = 0; s < g.nkinds; s++) kind_start[s + 1] = kind_start[s] + kind_count[s];
    for (int i = 0, *fill = calloc(g.nkinds + 1, sizeof(int)); i < g.natoms || (free(fill), 0); i++) kind_atoms[kind_start[g.kind[i]] + fill[g.kind[i]]++] = i;
    int *done = calloc(m->ndefs + 1, sizeof(int));
    for (int d0 = 0; d0 < m->ndefs && verdict != REF && verdict != IND; d0++) {
        const def_t *a0 = &m->defs[d0];
        if (done[d0] || !a0->has_k || a0->irrep < 0 || a0->kind > 1) continue;
        /* the definitions of this k and kind */
        int *members = malloc(sizeof(int) * m->ndefs), nm = 0, labels = 0;
        for (int d = d0; d < m->ndefs; d++) {
            const def_t *x = &m->defs[d];
            if (done[d] || !x->has_k || x->irrep < 0 || x->kind != a0->kind) continue;
            if (fabs(x->k[0] - a0->k[0]) > 1e-9 || fabs(x->k[1] - a0->k[1]) > 1e-9 || fabs(x->k[2] - a0->k[2]) > 1e-9) continue;
            done[d] = 1;
            members[nm++] = d;
        }
        for (int i = 0; i < nm; i++) {
            int fresh = 1;
            for (int j = 0; j < i; j++) fresh &= m->defs[members[j]].irrep != m->defs[members[i]].irrep;
            labels += fresh;
        }
        if (labels < 2) {
            free(members);
            continue;
        }
        /* one arm of the star periodic in the child cell, and its little group */
        double k[3];
        int have = 0;
        for (int q = 0; q < g.nrot && !have; q++) {
            double Qinv[9];
            if (inverse3(g.rot[q], Qinv)) continue;
            for (int j = 0; j < 3; j++) k[j] = a0->k[0] * Qinv[3 * j] + a0->k[1] * Qinv[3 * j + 1] + a0->k[2] * Qinv[3 * j + 2];
            have = 1;
            for (int r = 0; r < 3 && have; r++) {
                double v = g.B[3 * r] * k[0] + g.B[3 * r + 1] * k[1] + g.B[3 * r + 2] * k[2];
                have = fabs(v - round(v)) <= 1e-6;
            }
        }
        if (!have) {
            free(members);
            continue;
        }
        int *little = malloc(sizeof(int) * (g.nrot + 1)), nl = 0;
        for (int q = 0; q < g.nrot; q++) {
            double Qinv[9], kk[3];
            if (inverse3(g.rot[q], Qinv)) continue;
            for (int j = 0; j < 3; j++) kk[j] = k[0] * Qinv[3 * j] + k[1] * Qinv[3 * j + 1] + k[2] * Qinv[3 * j + 2] - k[j];
            int in = 1;
            for (int gen = 0; gen < 3 + g.ncent && in; gen++) {
                double v[3] = {gen == 0, gen == 1, gen == 2};
                if (gen >= 3) memcpy(v, g.cent[gen - 3], sizeof v);
                double dot = kk[0] * v[0] + kk[1] * v[1] + kk[2] * v[2];
                in = fabs(dot - round(dot)) <= 1e-6;
            }
            if (in) little[nl++] = q;
        }

        /* k-Bloch parts (interleaved complex, Cartesian for the products) and their little-group images */
        double *w = calloc((size_t)nm * n6, sizeof(double)), *img = malloc(sizeof(double) * (size_t)nm * nl * n6), *tmp = malloc(sizeof(double) * n6);
        double *norm = calloc(nm, sizeof(double)), *unorm = calloc(nm, sizeof(double));
        int failed = 0;
        for (int a = 0; a < nm && !failed; a++) {
            if (def_field(m, &g, &m->defs[members[a]], u)) {
                failed = 1;
                break;
            }
            for (int i = 0; i < g.natoms; i++) {
                double c[3];
                row_mul(u + 3 * i, L, c);
                unorm[a] += c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
            }
            double *wa = w + (size_t)a * n6;
            for (int r = 0; r < N; r++) {
                double ph = -2 * M_PI * (k[0] * tp[3 * r] + k[1] * tp[3 * r + 1] + k[2] * tp[3 * r + 2]), c = cos(ph) / N, s = sin(ph) / N;
                const int *perm = g.perm + (size_t)r * g.natoms;
                for (int i = 0; i < g.natoms; i++)
                    for (int x = 0; x < 3; x++) {
                        wa[6 * perm[i] + 2 * x] += c * u[3 * i + x];
                        wa[6 * perm[i] + 2 * x + 1] += s * u[3 * i + x];
                    }
            }
            for (int l = 0; l < nl && !failed; l++) failed = bloch_act(m, &g, xs, little[l], k, a0->kind == 1, wa, img + ((size_t)a * nl + l) * n6);
        }
        /* to Cartesian: each atom's row vector times L */
        for (size_t v = 0; v < (size_t)nm * (nl + 1) && !failed; v++) {
            double *x = v < (size_t)nm ? w + v * n6 : img + (v - nm) * n6;
            for (int i = 0; i < g.natoms; i++)
                for (int part = 0; part < 2; part++) {
                    double r[3] = {x[6 * i + part], x[6 * i + 2 + part], x[6 * i + 4 + part]}, c[3];
                    row_mul(r, L, c);
                    for (int a = 0; a < 3; a++) x[6 * i + 2 * a + part] = c[a];
                }
        }
        /* a Bloch part below the display precision of its mode is zero; the parent sites it lives on (operations keep
         * each site, so two parts meet only on their common sites) */
        unsigned char *present = calloc((size_t)nm * g.nkinds + 1, 1);
        for (int a = 0; a < nm && !failed; a++) {
            for (int e = 0; e < n6; e++) norm[a] += w[(size_t)a * n6 + e] * w[(size_t)a * n6 + e];
            norm[a] = sqrt(norm[a]);
            if (norm[a] < 1e-3 * sqrt(unorm[a])) norm[a] = 0;
            for (int i = 0; i < g.natoms; i++) {
                const double *x = w + (size_t)a * n6 + 6 * i;
                double s = 0;
                for (int e = 0; e < 6; e++) s += x[e] * x[e];
                if (norm[a] > 0 && sqrt(s) > 1e-4 * norm[a]) present[(size_t)a * g.nkinds + g.kind[i]] = 1;
            }
        }
        double work = 0;
        for (int a = 0; a < nm; a++)
            for (int b = 0; b < nm; b++)
                for (int s = 0; s < g.nkinds; s++)
                    if (present[(size_t)a * g.nkinds + s] && present[(size_t)b * g.nkinds + s]) work += (double)nl * kind_count[s];
        if (work > 4e9) failed = 1;
        if (failed) {
            verdict = IND;
            free(members), free(little), free(w), free(img), free(tmp), free(norm), free(unorm), free(present);
            break;
        }
        for (int a = 0; a < nm && verdict != REF; a++)
            for (int b = 0; b < nm && verdict != REF; b++) {
                if (m->defs[members[a]].irrep == m->defs[members[b]].irrep || norm[a] < 1e-9 || norm[b] < 1e-9) continue;
                const double *wb = w + (size_t)b * n6;
                int shared = 0;
                for (int s = 0; s < g.nkinds; s++) shared |= present[(size_t)a * g.nkinds + s] && present[(size_t)b * g.nkinds + s];
                if (!shared) continue;
                for (int l = 0; l < nl; l++) {
                    const double *v = img + ((size_t)a * nl + l) * n6;
                    double re = 0, im = 0;
                    for (int s = 0; s < g.nkinds; s++) {
                        if (!present[(size_t)a * g.nkinds + s] || !present[(size_t)b * g.nkinds + s]) continue;
                        for (int q = kind_start[s]; q < kind_start[s + 1]; q++) {
                            int e0 = 6 * kind_atoms[q];
                            for (int e = e0; e < e0 + 6; e += 2) {
                                re += v[e] * wb[e] + v[e + 1] * wb[e + 1];
                                im += v[e] * wb[e + 1] - v[e + 1] * wb[e];
                            }
                        }
                    }
                    if (verdict == NA) verdict = SAT;
                    if (sqrt(re * re + im * im) > 1e-2 * norm[a] * norm[b]) {
                        verdict = REF;
                        break;
                    }
                }
            }
        free(members), free(little), free(w), free(img), free(tmp), free(norm), free(unorm), free(present);
    }
    free(done), free(kind_count), free(kind_start), free(kind_atoms), free(xs), free(tp), free(u);
    sg_free(&g);
    return verdict;
}
