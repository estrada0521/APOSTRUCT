/* The three Local theorems (mode basis, group invariance, fixed-space dimension) on one apo mode-details payload,
 * at the display precision of the Python checks in verification/mathematics. Shares no code with core/. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PSCALE 100000

enum { SAT = 0, REF = 1, IND = 2, NA = 3 };

typedef struct {
    const char *p, *end;
} cur;

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
    int kind, label_atom, nf_zero, block, nrows;
    int *vec;
} def_t;

typedef struct {
    int label, multiplicity, rep[3], nexp;
    int *exp;
} orbit_t;

typedef struct {
    char err[256];
    int natoms, nblocks, ndefs, nstrain, norbits, failed;
    int *atom_marker, *atom_pos;
    block_t *blocks;
    def_t *defs;
    int64_t (*strain)[6];
    int *strain_zero;
    orbit_t *orbits;
    strtab markers, label_atoms;
    char *subgroup;
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

static void ws(cur *c)
{
    while (c->p < c->end && (*c->p == ' ' || *c->p == '\n' || *c->p == '\r' || *c->p == '\t')) c->p++;
}

static int skip(cur *c);

static int skip_string(cur *c)
{
    if (c->p >= c->end || *c->p != '"') return -1;
    for (c->p++; c->p < c->end; c->p++) {
        if (*c->p == '\\') c->p++;
        else if (*c->p == '"') {
            c->p++;
            return 0;
        }
    }
    return -1;
}

static int skip(cur *c)
{
    ws(c);
    if (c->p >= c->end) return -1;
    char ch = *c->p;
    if (ch == '"') return skip_string(c);
    if (ch == '{' || ch == '[') {
        char close = ch == '{' ? '}' : ']';
        c->p++;
        ws(c);
        if (c->p < c->end && *c->p == close) {
            c->p++;
            return 0;
        }
        for (;;) {
            if (ch == '{') {
                if (skip_string(c)) return -1;
                ws(c);
                if (c->p >= c->end || *c->p != ':') return -1;
                c->p++;
            }
            if (skip(c)) return -1;
            ws(c);
            if (c->p >= c->end) return -1;
            if (*c->p == ',') {
                c->p++;
                ws(c);
                continue;
            }
            if (*c->p == close) {
                c->p++;
                return 0;
            }
            return -1;
        }
    }
    while (c->p < c->end && !strchr(",]} \n\r\t", *c->p)) c->p++;
    return 0;
}

/* object iteration: key and klen set, c at the value; returns 1 per member, 0 at the end, -1 on error */
static int obj_begin(cur *c)
{
    ws(c);
    if (c->p >= c->end || *c->p != '{') return -1;
    c->p++;
    return 0;
}

static int obj_next(cur *c, const char **key, int *klen, int *first)
{
    ws(c);
    if (c->p < c->end && *c->p == '}') {
        c->p++;
        return 0;
    }
    if (!*first) {
        if (c->p >= c->end || *c->p != ',') return -1;
        c->p++;
        ws(c);
    }
    *first = 0;
    const char *s = c->p + 1;
    if (skip_string(c)) return -1;
    *key = s;
    *klen = (int)(c->p - 1 - s);
    ws(c);
    if (c->p >= c->end || *c->p != ':') return -1;
    c->p++;
    ws(c);
    return 1;
}

static int arr_begin(cur *c)
{
    ws(c);
    if (c->p >= c->end || *c->p != '[') return -1;
    c->p++;
    return 0;
}

static int arr_next(cur *c, int *first)
{
    ws(c);
    if (c->p < c->end && *c->p == ']') {
        c->p++;
        return 0;
    }
    if (!*first) {
        if (c->p >= c->end || *c->p != ',') return -1;
        c->p++;
        ws(c);
    }
    *first = 0;
    return 1;
}

static int keyis(const char *key, int klen, const char *want) { return (int)strlen(want) == klen && !memcmp(key, want, klen); }

static int is_null(cur *c)
{
    ws(c);
    return c->end - c->p >= 4 && !memcmp(c->p, "null", 4);
}

static int number(cur *c, double *out)
{
    ws(c);
    char buf[64];
    int n = 0;
    while (c->p < c->end && n < 63 && strchr("+-0123456789.eE", *c->p)) buf[n++] = *c->p++;
    if (!n) return -1;
    buf[n] = 0;
    char *e;
    *out = strtod(buf, &e);
    return *e || !isfinite(*out) ? -1 : 0;
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

static int string_value(cur *c, const char **s, int *len)
{
    ws(c);
    const char *start = c->p + 1;
    if (skip_string(c)) return -1;
    *s = start;
    *len = (int)(c->p - 1 - start);
    return 0;
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
    m->atom_pos = malloc(sizeof(int) * 3 * cap);
    if (arr_begin(c)) return fail(m, "atoms is not a list");
    for (int r; (r = arr_next(c, &first)) == 1;) {
        if (m->natoms == cap) {
            cap *= 2;
            m->atom_marker = realloc(m->atom_marker, sizeof(int) * cap);
            m->atom_pos = realloc(m->atom_pos, sizeof(int) * 3 * cap);
        }
        int i = m->natoms++, f = 1, seen = 0;
        m->atom_marker[i] = -1;
        const char *key;
        int klen;
        if (obj_begin(c)) return fail(m, "atom row");
        for (int k; (k = obj_next(c, &key, &klen, &f)) == 1;) {
            if (keyis(key, klen, "atom") && !is_null(c)) {
                const char *s;
                int len;
                if (string_value(c, &s, &len) || !len) return fail(m, "atom marker");
                m->atom_marker[i] = intern(&m->markers, s, len);
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
                has_label = 1;
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
    cur root = {text, text + len}, md;
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
                if (skip(&s)) return fail(m, "subgroup_details");
            }
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
    for (strtab *t = &m->markers; t <= &m->label_atoms; t++) {
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
