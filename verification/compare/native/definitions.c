/* The mode-definition comparison of one apo output against one Web Complete Mode Details text: row assignment,
 * signed permutations, mode-span proofs and the mode_definition layers of the verdict. */
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../json.h"

enum { DSP, MAG, OCC, STRAIN, NKIND };
static const char *const KIND[NKIND] = {"dsp", "mag", "occ", "strain"};

enum { S_OK, S_SIGN, S_SIGNED_PERMUTATION, S_FAMILY_BASIS, S_POSITION, S_VECTOR, S_COUNT, S_LABEL, S_NORMFACTOR, S_FAIL, S_TENSOR };
static const char *const STATUS[] = {"ok", "sign", "signed_permutation", "family_basis", "position", "vector", "count", "label", "normfactor", "fail", "tensor"};

typedef struct {
    char *atom;
    double xyz[3], d[3];
    int nx, nd;
} row_t;

typedef struct {
    char *label;
    int has_nf, nrows, pos_id;
    double nf;
    row_t *rows;
} def_t;

typedef struct {
    char *label;
    int has_nf, nc;
    double nf, c[6];
} strain_t;

typedef struct {
    int ok, kind;
    char *exact, *physical, *direction, *base;
} block_t;

typedef struct {
    int block_count, signed_count, family_count, minimax_count, direction_count;
    double max_condition, max_residual, bound, ratio;
} proof_t;

typedef struct {
    int n;
    int *start, *idx;
} cand_t;

typedef struct {
    int64_t key;
    cand_t *val;
} cslot;

typedef struct {
    double *pos;
    int n, id;
    uint64_t hash;
} pslot;

typedef struct vm vm_t;

typedef struct {
    vm_t *v;
    const def_t *W, *L;
    int nw, nl, kind, unamb;
    signed char *fact, *unambiguous;
    block_t *wb, *lb;
    int status, strict, count_match, validated, has_proof;
    int *mst, *mstrict, *label_match, *nf_match;
    int **assign;
    proof_t proof;
} payload_t;

typedef struct {
    int status, strict, count_match, nw, nl;
    int *mst, *mstrict;
} strain_payload_t;

struct vm {
    double tol, vtol;
    char err[512];
    int failed;
    def_t *web[3], *local[3];
    int nweb[3], nlocal[3];
    strain_t *wstrain, *lstrain;
    int nwstrain, nlstrain;
    payload_t P[3];
    strain_payload_t SP;
    pslot *pos;
    int npos, pos_cap;
    cslot *cand;
    int ncand, cand_cap;
    void **owned;
    int nowned, owned_cap;
    char *result;
};

static void *own(vm_t *v, void *p)
{
    if (v->nowned == v->owned_cap) {
        v->owned_cap = v->owned_cap ? 2 * v->owned_cap : 256;
        v->owned = realloc(v->owned, sizeof(void *) * v->owned_cap);
    }
    v->owned[v->nowned++] = p;
    return p;
}

static void *zalloc(vm_t *v, size_t n) { return own(v, calloc(1, n ? n : 1)); }

static char *sdup(vm_t *v, const char *s, size_t n)
{
    char *o = zalloc(v, n + 1);
    memcpy(o, s, n);
    return o;
}

static int fail(vm_t *v, const char *fmt, ...)
{
    if (!v->failed) {
        va_list ap;
        va_start(ap, fmt);
        vsnprintf(v->err, sizeof v->err, fmt, ap);
        va_end(ap);
    }
    v->failed = 1;
    return -1;
}

typedef struct {
    char *s;
    size_t n, cap;
} sb;

static void sb_put(sb *b, const char *s, size_t n)
{
    if (b->n + n + 1 > b->cap) {
        b->cap = 2 * (b->n + n + 1) + 64;
        b->s = realloc(b->s, b->cap);
    }
    memcpy(b->s + b->n, s, n);
    b->n += n;
    b->s[b->n] = 0;
}

static void sb_cat(sb *b, const char *s) { sb_put(b, s, strlen(s)); }

static void sb_fmt(sb *b, const char *fmt, ...)
{
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    sb_put(b, buf, (size_t)n < sizeof buf ? (size_t)n : sizeof buf - 1);
}

static void sb_json(sb *b, const char *s)
{
    if (!s) {
        sb_cat(b, "null");
        return;
    }
    sb_put(b, "\"", 1);
    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
        if (*p == '"' || *p == '\\') {
            sb_put(b, "\\", 1);
            sb_put(b, (const char *)p, 1);
        } else if (*p < 0x20) sb_fmt(b, "\\u%04x", *p);
        else sb_put(b, (const char *)p, 1);
    }
    sb_put(b, "\"", 1);
}

/* Python repr() of an ASCII str */
static void sb_repr(sb *b, const char *s)
{
    char q = strchr(s, '\'') && !strchr(s, '"') ? '"' : '\'';
    sb_put(b, &q, 1);
    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
        if (*p == '\\') sb_cat(b, "\\\\");
        else if (*p == (unsigned char)q) {
            sb_put(b, "\\", 1);
            sb_put(b, &q, 1);
        } else if (*p == '\n') sb_cat(b, "\\n");
        else if (*p == '\r') sb_cat(b, "\\r");
        else if (*p == '\t') sb_cat(b, "\\t");
        else if (*p < 0x20 || *p == 0x7f) sb_fmt(b, "\\x%02x", *p);
        else sb_put(b, (const char *)p, 1);
    }
    sb_put(b, &q, 1);
}

static double pymod1(double x)
{
    double m = fmod(x, 1.0);
    if (m != 0.0) {
        if (m < 0) m += 1.0;
    } else m = 0.0;
    return m;
}

static int close_to(double a, double b, double tol) { return isfinite(a) && isfinite(b) && fabs(a - b) <= tol; }

static int vec_close(const double *a, int na, const double *b, int nb, double tol, double sign)
{
    if (na != nb) return 0;
    for (int i = 0; i < na; i++)
        if (!close_to(a[i], sign * b[i], tol)) return 0;
    return 1;
}

static int frac_close(double a, double b, double tol)
{
    double diff = fabs(pymod1(a - b));
    return fmin(diff, fabs(diff - 1.0)) <= tol;
}

static int frac_vec_close(const row_t *a, const row_t *b, double tol)
{
    if (a->nx != b->nx) return 0;
    for (int i = 0; i < a->nx; i++)
        if (!frac_close(a->xyz[i], b->xyz[i], tol)) return 0;
    return 1;
}

static int isclose_abs(double a, double b, double tol)
{
    if (a == b) return 1;
    if (isinf(a) || isinf(b)) return 0;
    return fabs(b - a) <= tol;
}

static int optional_close(int ha, double a, int hb, double b, double tol)
{
    if (!ha || !hb) return !ha && !hb;
    return isclose_abs(a, b, tol);
}

/* ---- Web text ---- */

/* first occurrence of needle in [h, hend): Horspool over long spans, memchr over short lines */
static const char *find(const char *h, const char *hend, const char *needle)
{
    size_t n = strlen(needle);
    if (!h || hend - h < (ptrdiff_t)n) return 0;
    size_t len = (size_t)(hend - h);
    if (len < 512) {
        for (const char *p = h; (p = memchr(p, needle[0], (size_t)(hend - p))) && p + n <= hend; p++)
            if (!memcmp(p, needle, n)) return p;
        return 0;
    }
    size_t shift[256];
    for (int i = 0; i < 256; i++) shift[i] = n;
    for (size_t i = 0; i + 1 < n; i++) shift[(unsigned char)needle[i]] = n - 1 - i;
    for (size_t i = 0; i + n <= len; i += shift[(unsigned char)h[i + n - 1]])
        if (h[i + n - 1] == needle[n - 1] && !memcmp(h + i, needle, n - 1)) return h + i;
    return 0;
}

typedef struct {
    const char *p, *end;
} span;

static span section(const char *t, const char *tend, const char *start, const char *stop)
{
    span s = {t, t};
    const char *a = find(t, tend, start);
    if (!a) return s;
    s.p = a + strlen(start);
    s.end = tend;
    if (stop) {
        const char *b = find(s.p, tend, stop);
        if (b) s.end = b;
    }
    return s;
}

static void cut(span *s, const char *stop)
{
    const char *q = find(s->p, s->end, stop);
    if (q) s->end = q;
}

static int is_space(char c) { return c == ' ' || (c >= '\t' && c <= '\r') || (c >= 0x1c && c <= 0x1f); }

static int is_digit(char c) { return c >= '0' && c <= '9'; }

/* one match of [-+]?(?:\d+\.\d+|\d+)(?:[eE][-+]?\d+)? at p */
static const char *float_at(const char *p, const char *end)
{
    if (p < end && (*p == '+' || *p == '-')) p++;
    const char *d = p;
    while (p < end && is_digit(*p)) p++;
    if (p == d) return 0;
    if (p + 1 < end && *p == '.' && is_digit(p[1])) {
        p++;
        while (p < end && is_digit(*p)) p++;
    }
    if (p < end && (*p == 'e' || *p == 'E')) {
        const char *q = p + 1;
        if (q < end && (*q == '+' || *q == '-')) q++;
        if (q < end && is_digit(*q)) {
            while (q < end && is_digit(*q)) q++;
            p = q;
        }
    }
    return p;
}

typedef struct {
    double *x;
    int n, cap;
} dvec;

static void dpush(dvec *v, double x)
{
    if (v->n == v->cap) {
        v->cap = v->cap ? 2 * v->cap : 16;
        v->x = realloc(v->x, sizeof(double) * v->cap);
    }
    v->x[v->n++] = x;
}

static void numbers(const char *p, const char *end, dvec *out)
{
    out->n = 0;
    while (p < end) {
        const char *q = float_at(p, end);
        if (!q) {
            p++;
            continue;
        }
        double x = 0;
        decimal(p, (size_t)(q - p), &x);
        dpush(out, x);
        p = q;
    }
}

static void strip(const char **p, const char **end)
{
    while (*p < *end && is_space(**p)) (*p)++;
    while (*end > *p && is_space((*end)[-1])) (*end)--;
}

static int starts(const char *p, const char *end, const char *s)
{
    size_t n = strlen(s);
    return (size_t)(end - p) >= n && !memcmp(p, s, n);
}

typedef struct {
    def_t *m;
    int n, cap;
} mvec;

typedef struct {
    row_t *r;
    int n, cap;
} rvec;

static const char *const DSP_STOPS[] = {"Displacive mode amplitudes", "Magnetic mode definitions", "Magnetic mode amplitudes", "Occupancy mode definitions",
                                        "Occupancy mode amplitudes", "Parent-cell strain mode definitions", "Parent-cell strain mode amplitudes", 0};
static const char *const MAG_STOPS[] = {"Magnetic mode amplitudes", "Occupancy mode definitions", "Occupancy mode amplitudes", "Parent-cell strain mode definitions",
                                        "Parent-cell strain mode amplitudes", 0};
static const char *const OCC_STOPS[] = {"Occupancy mode amplitudes", "Parent-cell strain mode definitions", "Parent-cell strain mode amplitudes", 0};

static void flush_mode(vm_t *v, mvec *ms, char **label, int *has_nf, double *nf, rvec *rows)
{
    if (!*label) return;
    if (ms->n == ms->cap) {
        ms->cap = ms->cap ? 2 * ms->cap : 16;
        ms->m = realloc(ms->m, sizeof(def_t) * ms->cap);
    }
    def_t *m = &ms->m[ms->n++];
    m->label = *label;
    m->has_nf = *has_nf;
    m->nf = *nf;
    m->nrows = rows->n;
    m->rows = zalloc(v, sizeof(row_t) * (rows->n ? rows->n : 1));
    memcpy(m->rows, rows->r, sizeof(row_t) * rows->n);
    m->pos_id = -1;
    *label = 0;
    *has_nf = 0;
    *nf = 0;
    rows->n = 0;
}

static void parse_definitions(vm_t *v, const char *t, const char *tend, int kind)
{
    static const char *const START[] = {"Displacive mode definitions", "Magnetic mode definitions", "Occupancy mode definitions"};
    const char *const *stops = kind == DSP ? DSP_STOPS : kind == MAG ? MAG_STOPS : OCC_STOPS;
    span s = section(t, tend, START[kind], 0);
    for (int i = 0; stops[i]; i++) cut(&s, stops[i]);
    for (int i = 0; DSP_STOPS[i]; i++) cut(&s, DSP_STOPS[i]);
    mvec ms = {0};
    rvec rows = {0};
    dvec nums = {0};
    char *label = 0;
    int has_nf = 0;
    double nf = 0;
    for (const char *p = s.p; p < s.end;) {
        const char *e = memchr(p, '\n', (size_t)(s.end - p));
        if (!e) e = s.end;
        const char *raw = p, *raw_end = e;
        p = e + 1;
        const char *a = raw, *b = raw_end;
        strip(&a, &b);
        if (a == b || starts(a, b, "atom ")) continue;
        const char *nfpos = find(a, b, "normfactor");
        if (nfpos) {
            flush_mode(v, &ms, &label, &has_nf, &nf, &rows);
            const char *la = a, *lb = nfpos;
            strip(&la, &lb);
            int open = 0, closed = 0;
            const char *last_open = 0;
            for (const char *q = la; q < lb; q++) {
                if (*q == '(') open++, last_open = q;
                else if (*q == ')') closed++;
            }
            int fix = open == closed + 1 && last_open && !memchr(last_open, ')', (size_t)(lb - last_open));
            label = zalloc(v, (size_t)(lb - la) + 2);
            memcpy(label, la, (size_t)(lb - la));
            if (fix) label[lb - la] = ')';
            numbers(nfpos + 10, b, &nums);
            has_nf = nums.n > 0;
            nf = has_nf ? nums.x[0] : 0;
            continue;
        }
        numbers(a, b, &nums);
        if (kind == OCC && nums.n >= 4) {
            dpush(&nums, 0.0);
            dpush(&nums, 0.0);
        }
        if (!label || nums.n < 6) continue;
        const char *t0 = a;
        while (t0 < b && !is_space(*t0)) t0++;
        if (rows.n == rows.cap) {
            rows.cap = rows.cap ? 2 * rows.cap : 64;
            rows.r = realloc(rows.r, sizeof(row_t) * rows.cap);
        }
        row_t *r = &rows.r[rows.n++];
        memset(r, 0, sizeof *r);
        if (float_at(a, t0) != t0) r->atom = sdup(v, a, (size_t)(t0 - a));
        const double *x = nums.x + nums.n - 6;
        for (int i = 0; i < 3; i++) r->xyz[i] = x[i], r->d[i] = x[3 + i];
        r->nx = r->nd = 3;
    }
    flush_mode(v, &ms, &label, &has_nf, &nf, &rows);
    v->web[kind] = own(v, ms.m ? ms.m : malloc(1));
    v->nweb[kind] = ms.n;
    free(rows.r);
    free(nums.x);
}

static void parse_strain(vm_t *v, const char *t, const char *tend)
{
    span s = section(t, tend, "Parent-cell strain mode definitions", "Parent-cell strain mode amplitudes");
    dvec nums = {0};
    char *label = 0;
    int has_nf = 0, n = 0, cap = 0;
    double nf = 0;
    strain_t *out = 0;
    for (const char *p = s.p; p < s.end;) {
        const char *e = memchr(p, '\n', (size_t)(s.end - p));
        if (!e) e = s.end;
        const char *a = p, *b = e;
        p = e + 1;
        strip(&a, &b);
        if (a == b || starts(a, b, "e1")) continue;
        const char *nfpos = find(a, b, "normfactor");
        if (nfpos) {
            const char *la = a, *lb = nfpos;
            strip(&la, &lb);
            label = sdup(v, la, (size_t)(lb - la));
            numbers(nfpos + 10, b, &nums);
            has_nf = nums.n > 0;
            nf = has_nf ? nums.x[0] : 0;
            continue;
        }
        numbers(a, b, &nums);
        if (!label || nums.n != 6) continue;
        if (n == cap) {
            cap = cap ? 2 * cap : 8;
            out = realloc(out, sizeof(strain_t) * cap);
        }
        strain_t *m = &out[n++];
        m->label = label;
        m->has_nf = has_nf;
        m->nf = nf;
        m->nc = 6;
        memcpy(m->c, nums.x, sizeof m->c);
        label = 0;
        has_nf = 0;
        nf = 0;
    }
    v->wstrain = own(v, out ? out : malloc(1));
    v->nwstrain = n;
    free(nums.x);
}

/* ---- apo output ---- */

static int read_numbers(cur *c, double *out, int max, int *n)
{
    int first = 1;
    *n = 0;
    if (arr_begin(c)) return -1;
    for (int r; (r = arr_next(c, &first));) {
        if (r < 0) return -1;
        double x;
        if (number(c, &x)) return -1;
        if (*n < max) out[*n] = x;
        (*n)++;
    }
    return 0;
}

static int read_optional_string(vm_t *v, cur *c, char **out)
{
    if (is_null(c)) {
        *out = 0;
        return skip(c);
    }
    if (string_dup(c, out)) return -1;
    own(v, *out);
    return 0;
}

typedef struct {
    char *atom;
    double xyz[3];
    int nx;
} atom_t;

static int read_atoms(vm_t *v, cur *c, atom_t **out, int *n)
{
    int first = 1, cap = 0;
    *n = 0;
    *out = 0;
    if (arr_begin(c)) return -1;
    for (int r; (r = arr_next(c, &first));) {
        if (r < 0) return -1;
        if (*n == cap) {
            cap = cap ? 2 * cap : 64;
            *out = realloc(*out, sizeof(atom_t) * cap);
        }
        atom_t *a = &(*out)[(*n)++];
        memset(a, 0, sizeof *a);
        int f2 = 1;
        const char *key;
        int klen;
        if (obj_begin(c)) return -1;
        for (int q; (q = obj_next(c, &key, &klen, &f2));) {
            if (q < 0) return -1;
            if (keyis(key, klen, "atom")) {
                if (read_optional_string(v, c, &a->atom)) return -1;
            } else if (keyis(key, klen, "xyz")) {
                if (read_numbers(c, a->xyz, 3, &a->nx)) return -1;
            } else if (skip(c)) return -1;
        }
    }
    own(v, *out);
    return 0;
}

typedef struct {
    int *idx;
    int n;
} iblock_t;

static int read_blocks(vm_t *v, cur *c, iblock_t **out, int *n)
{
    int first = 1, cap = 0;
    *n = 0;
    *out = 0;
    if (arr_begin(c)) return -1;
    for (int r; (r = arr_next(c, &first));) {
        if (r < 0) return -1;
        if (*n == cap) {
            cap = cap ? 2 * cap : 16;
            *out = realloc(*out, sizeof(iblock_t) * cap);
        }
        iblock_t *b = &(*out)[(*n)++];
        int f2 = 1, bcap = 0;
        b->n = 0;
        b->idx = 0;
        if (arr_begin(c)) return -1;
        for (int q; (q = arr_next(c, &f2));) {
            if (q < 0) return -1;
            double x;
            if (number(c, &x) || x != floor(x)) return -1;
            if (b->n == bcap) {
                bcap = bcap ? 2 * bcap : 16;
                b->idx = realloc(b->idx, sizeof(int) * bcap);
            }
            b->idx[b->n++] = (int)x;
        }
        own(v, b->idx);
    }
    own(v, *out);
    return 0;
}

static int read_row_objects(vm_t *v, cur *c, def_t *m)
{
    int first = 1, cap = 0;
    m->nrows = 0;
    m->rows = 0;
    if (is_null(c)) return skip(c);
    if (arr_begin(c)) return -1;
    for (int r; (r = arr_next(c, &first));) {
        if (r < 0) return -1;
        if (m->nrows == cap) {
            cap = cap ? 2 * cap : 64;
            m->rows = realloc(m->rows, sizeof(row_t) * cap);
        }
        row_t *row = &m->rows[m->nrows++];
        memset(row, 0, sizeof *row);
        int f2 = 1;
        const char *key;
        int klen;
        if (obj_begin(c)) return -1;
        for (int q; (q = obj_next(c, &key, &klen, &f2));) {
            if (q < 0) return -1;
            if (keyis(key, klen, "atom")) {
                if (read_optional_string(v, c, &row->atom)) return -1;
            } else if (keyis(key, klen, "xyz") && !is_null(c)) {
                if (read_numbers(c, row->xyz, 3, &row->nx)) return -1;
            } else if (keyis(key, klen, "dxyz") && !is_null(c)) {
                if (read_numbers(c, row->d, 3, &row->nd)) return -1;
            } else if (skip(c)) return -1;
        }
    }
    own(v, m->rows);
    return 0;
}

static int read_definitions(vm_t *v, cur *c, int kind, const atom_t *atoms, int natoms, const iblock_t *blocks, int nblocks)
{
    int first = 1, cap = 0, n = 0;
    def_t *out = 0;
    if (arr_begin(c)) return fail(v, "local %s definitions", KIND[kind]);
    for (int r; (r = arr_next(c, &first));) {
        if (r < 0) return fail(v, "local %s definitions", KIND[kind]);
        if (n == cap) {
            cap = cap ? 2 * cap : 16;
            out = realloc(out, sizeof(def_t) * cap);
        }
        def_t *m = &out[n++];
        memset(m, 0, sizeof *m);
        m->pos_id = -1;
        int f2 = 1, block = -1, has_rows = 0;
        cur vectors = {0, 0};
        const char *key;
        int klen;
        if (obj_begin(c)) return fail(v, "local definition");
        for (int q; (q = obj_next(c, &key, &klen, &f2));) {
            if (q < 0) return fail(v, "local definition");
            if (keyis(key, klen, "label")) {
                if (read_optional_string(v, c, &m->label)) return fail(v, "local definition label");
            } else if (keyis(key, klen, "normfactor") && !is_null(c)) {
                if (number(c, &m->nf)) return fail(v, "local normfactor");
                m->has_nf = 1;
            } else if (keyis(key, klen, "atom_block")) {
                double x;
                if (number(c, &x) || x != floor(x)) return fail(v, "local atom_block");
                block = (int)x;
            } else if (keyis(key, klen, "vectors")) {
                vectors = *c;
                if (skip(c)) return fail(v, "local vectors");
            } else if (keyis(key, klen, "rows")) {
                has_rows = 1;
                if (read_row_objects(v, c, m)) return fail(v, "local rows");
            } else if (skip(c)) return fail(v, "local definition");
        }
        if (has_rows) continue;
        if (block < 0 || block >= nblocks) return fail(v, "IndexError: list index out of range");
        const iblock_t *b = &blocks[block];
        m->rows = zalloc(v, sizeof(row_t) * (b->n ? b->n : 1));
        m->nrows = b->n;
        int f3 = 1, k = 0;
        if (!vectors.p || arr_begin(&vectors)) return fail(v, "local vectors");
        for (int q; (q = arr_next(&vectors, &f3));) {
            if (q < 0) return fail(v, "local vectors");
            if (k >= b->n) return fail(v, "ValueError: zip() argument 2 is longer than argument 1");
            int a = b->idx[k];
            if (a < 0 || a >= natoms) return fail(v, "IndexError: list index out of range");
            row_t *row = &m->rows[k++];
            row->atom = atoms[a].atom;
            memcpy(row->xyz, atoms[a].xyz, sizeof row->xyz);
            row->nx = atoms[a].nx;
            ws(&vectors);
            if (vectors.p < vectors.end && *vectors.p == '[') {
                if (read_numbers(&vectors, row->d, 3, &row->nd)) return fail(v, "local vector");
            } else {
                if (number(&vectors, &row->d[0])) return fail(v, "local vector");
                row->nd = 3;
            }
        }
        if (k != b->n) return fail(v, "ValueError: zip() argument 2 is shorter than argument 1");
    }
    v->local[kind] = own(v, out ? out : malloc(1));
    v->nlocal[kind] = n;
    return 0;
}

static int read_strain_definitions(vm_t *v, cur *c)
{
    int first = 1, cap = 0, n = 0;
    strain_t *out = 0;
    if (arr_begin(c)) return fail(v, "local strain definitions");
    for (int r; (r = arr_next(c, &first));) {
        if (r < 0) return fail(v, "local strain definitions");
        if (n == cap) {
            cap = cap ? 2 * cap : 8;
            out = realloc(out, sizeof(strain_t) * cap);
        }
        strain_t *m = &out[n++];
        memset(m, 0, sizeof *m);
        double comps[6], tensor[6];
        int nc = 0, nt = 0, f2 = 1;
        const char *key;
        int klen;
        if (obj_begin(c)) return fail(v, "local strain definition");
        for (int q; (q = obj_next(c, &key, &klen, &f2));) {
            if (q < 0) return fail(v, "local strain definition");
            if (keyis(key, klen, "label")) {
                if (read_optional_string(v, c, &m->label)) return fail(v, "local strain label");
            } else if (keyis(key, klen, "normfactor") && !is_null(c)) {
                if (number(c, &m->nf)) return fail(v, "local strain normfactor");
                m->has_nf = 1;
            } else if (keyis(key, klen, "components") && !is_null(c)) {
                if (read_numbers(c, comps, 6, &nc)) return fail(v, "local strain components");
            } else if (keyis(key, klen, "tensor") && !is_null(c)) {
                if (read_numbers(c, tensor, 6, &nt)) return fail(v, "local strain tensor");
            } else if (skip(c)) return fail(v, "local strain definition");
        }
        if (nc) {
            m->nc = nc;
            memcpy(m->c, comps, sizeof comps);
        } else {
            m->nc = nt;
            memcpy(m->c, tensor, sizeof tensor);
        }
    }
    v->lstrain = own(v, out ? out : malloc(1));
    v->nlstrain = n;
    return 0;
}

static int descend(cur *c, const char *want)
{
    cur probe = *c;
    int first = 1;
    const char *key;
    int klen;
    if (obj_begin(&probe)) return 0;
    for (int r; (r = obj_next(&probe, &key, &klen, &first));) {
        if (r < 0) return 0;
        if (keyis(key, klen, want)) {
            *c = probe;
            return 1;
        }
        if (skip(&probe)) return 0;
    }
    return 0;
}

static int read_local(vm_t *v, const char *text, int64_t len)
{
    cur c = {text, text + len};
    cur root = c;
    if (descend(&c, "preview") && descend(&c, "selected") && descend(&c, "mode_details")) {
    } else if (c = root, descend(&c, "selected") && descend(&c, "mode_details")) {
    } else c = root;
    cur at[6] = {{0, 0}};
    static const char *const keys[6] = {"atoms", "atom_blocks", "displacive_definitions", "magnetic_definitions", "occupancy_definitions", "strain_definitions"};
    int first = 1;
    const char *key;
    int klen;
    if (obj_begin(&c)) return fail(v, "local mode details");
    for (int r; (r = obj_next(&c, &key, &klen, &first));) {
        if (r < 0) return fail(v, "local mode details");
        for (int i = 0; i < 6; i++)
            if (keyis(key, klen, keys[i]) && !is_null(&c)) at[i] = c;
        if (skip(&c)) return fail(v, "local mode details");
    }
    atom_t *atoms = 0;
    iblock_t *blocks = 0;
    int natoms = 0, nblocks = 0;
    if (at[0].p && read_atoms(v, &at[0], &atoms, &natoms)) return fail(v, "local atoms");
    if (at[1].p && read_blocks(v, &at[1], &blocks, &nblocks)) return fail(v, "local atom_blocks");
    for (int k = 0; k < 3; k++) {
        if (!at[2 + k].p) {
            v->local[k] = zalloc(v, sizeof(def_t));
            v->nlocal[k] = 0;
        } else if (read_definitions(v, &at[2 + k], k, atoms, natoms, blocks, nblocks)) return -1;
    }
    if (!at[5].p) {
        v->lstrain = zalloc(v, sizeof(strain_t));
        v->nlstrain = 0;
    } else if (read_strain_definitions(v, &at[5])) return -1;
    return 0;
}

static int cmp_str(const void *a, const void *b) { return strcmp(*(char *const *)a, *(char *const *)b); }

static int cmp_label(const void *a, const void *b) { return strcmp((*(def_t *const *)a)->label, (*(def_t *const *)b)->label); }

/* local modes reordered by unique complete labels when both sides carry the same label set */
static void exact_web_order(vm_t *v, int kind)
{
    int n = v->nweb[kind];
    def_t *W = v->web[kind], *L = v->local[kind];
    if (n != v->nlocal[kind] || !n) return;
    char **ws_ = malloc(sizeof(char *) * n), **ls = malloc(sizeof(char *) * n);
    int ok = 1;
    for (int i = 0; i < n && ok; i++) {
        ws_[i] = W[i].label;
        ls[i] = L[i].label;
        if (!ws_[i] || !*ws_[i] || !ls[i] || !*ls[i]) ok = 0;
    }
    if (ok) {
        qsort(ws_, n, sizeof(char *), cmp_str);
        qsort(ls, n, sizeof(char *), cmp_str);
        for (int i = 0; i < n && ok; i++)
            if (strcmp(ws_[i], ls[i]) || (i && !strcmp(ws_[i], ws_[i - 1])) || (i && !strcmp(ls[i], ls[i - 1]))) ok = 0;
    }
    free(ws_);
    free(ls);
    if (!ok) return;
    def_t *o = zalloc(v, sizeof(def_t) * n);
    def_t **sorted = malloc(sizeof(def_t *) * n);
    for (int i = 0; i < n; i++) sorted[i] = &L[i];
    qsort(sorted, n, sizeof(def_t *), cmp_label);
    for (int i = 0; i < n; i++) {
        int lo = 0, hi = n - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (strcmp(sorted[mid]->label, W[i].label) < 0) lo = mid + 1;
            else hi = mid;
        }
        o[i] = *sorted[lo];
    }
    free(sorted);
    v->local[kind] = o;
}

/* ---- periodic candidates ---- */

static uint64_t mix64(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    return x ^ (x >> 33);
}

static int position_id(vm_t *v, def_t *m)
{
    if (m->pos_id != -1) return m->pos_id;
    double *pos = malloc(sizeof(double) * 3 * (m->nrows ? m->nrows : 1));
    uint64_t h = (uint64_t)m->nrows * 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < m->nrows; i++) {
        if (m->rows[i].nx != 3) {
            free(pos);
            return m->pos_id = -2;
        }
        for (int a = 0; a < 3; a++) {
            double f = pymod1(m->rows[i].xyz[a]);
            if (!isfinite(f)) {
                free(pos);
                return m->pos_id = -2;
            }
            pos[3 * i + a] = f;
            uint64_t bits;
            memcpy(&bits, &f, 8);
            h = mix64(h ^ bits);
        }
    }
    if (2 * (v->npos + 1) > v->pos_cap) {
        int cap = v->pos_cap ? 2 * v->pos_cap : 256;
        pslot *t = calloc(cap, sizeof(pslot));
        for (int i = 0; i < v->pos_cap; i++)
            if (v->pos[i].pos) {
                int j = (int)(v->pos[i].hash & (uint64_t)(cap - 1));
                while (t[j].pos) j = (j + 1) & (cap - 1);
                t[j] = v->pos[i];
            }
        free(v->pos);
        v->pos = t;
        v->pos_cap = cap;
    }
    int j = (int)(h & (uint64_t)(v->pos_cap - 1));
    while (v->pos[j].pos) {
        if (v->pos[j].hash == h && v->pos[j].n == m->nrows && !memcmp(v->pos[j].pos, pos, sizeof(double) * 3 * m->nrows)) {
            free(pos);
            return m->pos_id = v->pos[j].id;
        }
        j = (j + 1) & (v->pos_cap - 1);
    }
    v->pos[j].pos = pos;
    v->pos[j].n = m->nrows;
    v->pos[j].hash = h;
    v->pos[j].id = v->npos++;
    return m->pos_id = v->pos[j].id;
}

static int cmp_int(const void *a, const void *b) { return (*(const int *)a > *(const int *)b) - (*(const int *)a < *(const int *)b); }

static cand_t *compute_candidates(const def_t *w, const def_t *l, double tol)
{
    cand_t *c = calloc(1, sizeof *c);
    c->n = w->nrows;
    c->start = calloc(w->nrows + 1, sizeof(int));
    double bw = 2.0 * tol;
    int count = (int)floor(1.0 / bw);
    if (count < 1) count = 1;
    int nl = l->nrows, hcap = 16;
    while (hcap < 2 * nl + 2) hcap <<= 1;
    int64_t *hkey = malloc(sizeof(int64_t) * hcap);
    int *head = malloc(sizeof(int) * hcap), *tail = malloc(sizeof(int) * hcap), *next = malloc(sizeof(int) * (nl ? nl : 1));
    for (int i = 0; i < hcap; i++) head[i] = -1;
    for (int i = 0; i < nl; i++) {
        next[i] = -1;
        const row_t *r = &l->rows[i];
        if (r->nx != 3) continue;
        int64_t k[3];
        int finite = 1;
        for (int a = 0; a < 3; a++) {
            double f = pymod1(r->xyz[a]);
            if (!isfinite(f)) finite = 0;
            double q = floor(f / bw);
            k[a] = q > count - 1 ? count - 1 : (int64_t)q;
        }
        if (!finite) continue;
        int64_t key = (k[0] * count + k[1]) * count + k[2];
        int j = (int)(mix64((uint64_t)key) & (uint64_t)(hcap - 1));
        while (head[j] >= 0 && hkey[j] != key) j = (j + 1) & (hcap - 1);
        if (head[j] < 0) {
            hkey[j] = key;
            head[j] = tail[j] = i;
        } else {
            next[tail[j]] = i;
            tail[j] = i;
        }
    }
    int cap = 64, n = 0, *buf = malloc(sizeof(int) * 64);
    int *tmp = malloc(sizeof(int) * (nl + 1)), ntmp;
    for (int wi = 0; wi < w->nrows; wi++) {
        c->start[wi] = n;
        const row_t *r = &w->rows[wi];
        int64_t ctr[3];
        for (int a = 0; a < 3; a++) {
            double q = floor(pymod1(r->xyz[a]) / bw);
            ctr[a] = q > count - 1 ? count - 1 : (int64_t)q;
        }
        ntmp = 0;
        for (int d0 = -1; d0 <= 1; d0++)
            for (int d1 = -1; d1 <= 1; d1++)
                for (int d2 = -1; d2 <= 1; d2++) {
                    int64_t k0 = ((ctr[0] + d0) % count + count) % count, k1 = ((ctr[1] + d1) % count + count) % count, k2 = ((ctr[2] + d2) % count + count) % count;
                    int64_t key = (k0 * count + k1) * count + k2;
                    int j = (int)(mix64((uint64_t)key) & (uint64_t)(hcap - 1));
                    while (head[j] >= 0 && hkey[j] != key) j = (j + 1) & (hcap - 1);
                    for (int i = head[j]; i >= 0; i = next[i]) tmp[ntmp++] = i;
                }
        qsort(tmp, ntmp, sizeof(int), cmp_int);
        for (int t = 0; t < ntmp; t++) {
            if (t && tmp[t] == tmp[t - 1]) continue;
            if (!frac_vec_close(r, &l->rows[tmp[t]], tol)) continue;
            if (n == cap) buf = realloc(buf, sizeof(int) * (cap *= 2));
            buf[n++] = tmp[t];
        }
    }
    c->start[w->nrows] = n;
    c->idx = buf;
    free(tmp);
    free(hkey);
    free(head);
    free(tail);
    free(next);
    return c;
}

static void cand_free(cand_t *c)
{
    if (!c) return;
    free(c->start);
    free(c->idx);
    free(c);
}

/* candidates are cached per pair of position contents */
static cand_t *candidates(vm_t *v, def_t *w, def_t *l, int *owned)
{
    int a = position_id(v, w), b = position_id(v, l);
    if (a < 0 || b < 0) {
        *owned = 1;
        return compute_candidates(w, l, v->tol);
    }
    *owned = 0;
    int64_t key = (int64_t)a << 32 | (uint32_t)b;
    if (2 * (v->ncand + 1) > v->cand_cap) {
        int cap = v->cand_cap ? 2 * v->cand_cap : 256;
        cslot *t = calloc(cap, sizeof(cslot));
        for (int i = 0; i < v->cand_cap; i++)
            if (v->cand[i].val) {
                int j = (int)(mix64((uint64_t)v->cand[i].key) & (uint64_t)(cap - 1));
                while (t[j].val) j = (j + 1) & (cap - 1);
                t[j] = v->cand[i];
            }
        free(v->cand);
        v->cand = t;
        v->cand_cap = cap;
    }
    int j = (int)(mix64((uint64_t)key) & (uint64_t)(v->cand_cap - 1));
    while (v->cand[j].val) {
        if (v->cand[j].key == key) return v->cand[j].val;
        j = (j + 1) & (v->cand_cap - 1);
    }
    v->cand[j].key = key;
    v->cand[j].val = compute_candidates(w, l, v->tol);
    v->ncand++;
    return v->cand[j].val;
}

/* ---- row assignment (augmenting paths in Web row order, candidates ascending) ---- */

typedef struct {
    const def_t *w, *l;
    const cand_t *c;
    int pred, stamp;
    double vtol;
    int *l2w, *seen;
} kuhn_t;

static int row_pred(const kuhn_t *k, int wi, int li)
{
    if (k->pred == 2) return 1;
    return vec_close(k->w->rows[wi].d, 3, k->l->rows[li].d, k->l->rows[li].nd, k->vtol, k->pred ? -1.0 : 1.0);
}

static int augment(kuhn_t *k, int wi)
{
    for (int t = k->c->start[wi]; t < k->c->start[wi + 1]; t++) {
        int li = k->c->idx[t];
        if (k->seen[li] == k->stamp || !row_pred(k, wi, li)) continue;
        k->seen[li] = k->stamp;
        int prev = k->l2w[li];
        if (prev < 0 || augment(k, prev)) {
            k->l2w[li] = wi;
            return 1;
        }
    }
    return 0;
}

static int row_matching(const def_t *w, const def_t *l, const cand_t *c, int pred, double vtol, int *w2l)
{
    int n = w->nrows;
    if (n != l->nrows) return 0;
    kuhn_t k = {w, l, c, pred, 0, vtol, malloc(sizeof(int) * (n ? n : 1)), malloc(sizeof(int) * (n ? n : 1))};
    for (int i = 0; i < n; i++) k.l2w[i] = -1, k.seen[i] = -1;
    int ok = 1;
    for (int wi = 0; wi < n && ok; wi++) {
        k.stamp = wi;
        ok = augment(&k, wi);
    }
    if (ok && w2l)
        for (int li = 0; li < n; li++) w2l[k.l2w[li]] = li;
    free(k.l2w);
    free(k.seen);
    return ok;
}

/* S_OK, S_SIGN, S_VECTOR (positions only), S_POSITION, or S_COUNT for unequal row counts */
static int pair_status(payload_t *P, int wi, int li, int **assign)
{
    signed char *f = &P->fact[(size_t)wi * P->nl + li];
    if (*f >= 0 && !assign) return *f;
    def_t *w = (def_t *)&P->W[wi], *l = (def_t *)&P->L[li];
    int code;
    int *a = assign ? malloc(sizeof(int) * (w->nrows ? w->nrows : 1)) : 0;
    if (w->nrows != l->nrows) code = S_COUNT;
    else {
        int owned;
        cand_t *c = candidates(P->v, w, l, &owned);
        if (row_matching(w, l, c, 0, P->v->vtol, a)) code = S_OK;
        else if (row_matching(w, l, c, 1, P->v->vtol, a)) code = S_SIGN;
        else if (row_matching(w, l, c, 2, P->v->vtol, a)) code = S_VECTOR;
        else code = S_POSITION;
        if (owned) cand_free(c);
    }
    *f = (signed char)code;
    if (assign) {
        if (code == S_COUNT || code == S_POSITION) {
            free(a);
            a = 0;
        }
        *assign = a ? own(P->v, a) : 0;
    }
    return code;
}

static int position_vectors_unambiguous(payload_t *P, int wi, int li)
{
    signed char *f = &P->unambiguous[(size_t)wi * P->nl + li];
    if (*f >= 0) return *f;
    def_t *w = (def_t *)&P->W[wi], *l = (def_t *)&P->L[li];
    int owned, ok = 1;
    cand_t *c = candidates(P->v, w, l, &owned);
    for (int r = 0; r < w->nrows && ok; r++) {
        int s = c->start[r], e = c->start[r + 1];
        if (s == e) ok = 0;
        for (int t = s + 1; t < e && ok; t++) {
            const row_t *a = &l->rows[c->idx[s]], *b = &l->rows[c->idx[t]];
            if (!vec_close(a->d, a->nd, b->d, b->nd, P->v->tol, 1.0)) ok = 0;
        }
    }
    if (owned) cand_free(c);
    return *f = (signed char)ok;
}

/* ---- labels ---- */

typedef struct {
    int64_t n, d;
} frac_t;

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

/* Fraction(text) for sign, digits and an optional /digits or .digits[e±digits]; 0 when Python would raise or overflow */
static int parse_fraction(const char *s, size_t n, frac_t *out)
{
    size_t i = 0;
    int neg = 0;
    if (i < n && (s[i] == '+' || s[i] == '-')) neg = s[i++] == '-';
    if (!(i < n && (is_digit(s[i]) || (s[i] == '.' && i + 1 < n && is_digit(s[i + 1]))))) return 0;
    __int128 num = 0, den = 1;
    while (i < n && is_digit(s[i])) {
        num = num * 10 + (s[i++] - '0');
        if (num > INT64_MAX) return 0;
    }
    if (i < n && s[i] == '/') {
        i++;
        if (!(i < n && is_digit(s[i]))) return 0;
        den = 0;
        while (i < n && is_digit(s[i])) {
            den = den * 10 + (s[i++] - '0');
            if (den > INT64_MAX) return 0;
        }
        if (!den) return 0;
    } else {
        if (i < n && s[i] == '.') {
            i++;
            while (i < n && is_digit(s[i])) {
                num = num * 10 + (s[i++] - '0');
                den *= 10;
                if (num > INT64_MAX || den > INT64_MAX) return 0;
            }
        }
        if (i < n && (s[i] == 'e' || s[i] == 'E')) {
            i++;
            int eneg = 0, e = 0;
            if (i < n && (s[i] == '+' || s[i] == '-')) eneg = s[i++] == '-';
            if (!(i < n && is_digit(s[i]))) return 0;
            while (i < n && is_digit(s[i])) {
                e = e * 10 + (s[i++] - '0');
                if (e > 18) return 0;
            }
            for (int k = 0; k < e; k++) {
                if (eneg) den *= 10;
                else num *= 10;
                if (num > INT64_MAX || den > INT64_MAX) return 0;
            }
        }
    }
    if (i != n) return 0;
    int64_t g = gcd64((int64_t)num, (int64_t)den);
    out->n = (int64_t)(num / (g ? g : 1)) * (neg ? -1 : 1);
    out->d = (int64_t)(den / (g ? g : 1));
    return 1;
}

/* physical_k_label: a numeric three-component K prefix folded modulo reciprocal vectors */
static char *physical_k_label(vm_t *v, const char *value)
{
    const char *open = strchr(value, '['), *close = open ? strchr(open + 1, ']') : 0;
    size_t n = strlen(value);
    if (!open || !close || close == open + 1 || !close[1]) return sdup(v, value, n);
    for (const char *q = value; q < open; q++)
        if (*q == ']') return sdup(v, value, n);
    for (const char *q = open + 1; q < close; q++)
        if (*q == '[') return sdup(v, value, n);
    for (const char *q = close + 1; *q; q++)
        if (*q == '[' || *q == ']') return sdup(v, value, n);
    const char *parts[4];
    int np = 0;
    parts[np++] = open + 1;
    for (const char *q = open + 1; q < close; q++)
        if (*q == ',') {
            if (np == 4) return sdup(v, value, n);
            parts[np++] = q + 1;
        }
    if (np != 3) return sdup(v, value, n);
    sb b = {0};
    sb_put(&b, value, (size_t)(open - value));
    sb_put(&b, "[", 1);
    for (int k = 0; k < 3; k++) {
        const char *a = parts[k], *e = k < 2 ? parts[k + 1] - 1 : close;
        strip(&a, &e);
        frac_t f;
        if (!parse_fraction(a, (size_t)(e - a), &f)) {
            free(b.s);
            return sdup(v, value, n);
        }
        int64_t r = f.n % f.d;
        if (r < 0) r += f.d;
        if (k) sb_put(&b, ",", 1);
        if (f.d == 1 || r == 0) sb_fmt(&b, "%lld", (long long)r);
        else {
            int64_t g = gcd64(r, f.d);
            sb_fmt(&b, "%lld/%lld", (long long)(r / g), (long long)(f.d / g));
        }
    }
    sb_put(&b, "]", 1);
    sb_cat(&b, close + 1);
    char *out = sdup(v, b.s, b.n);
    free(b.s);
    return out;
}

/* (?<![A-Za-z0-9_.])(\d*)\.(\d+)(?![\d.]) with redundant zeroes dropped */
static char *canonical_decimals(vm_t *v, const char *s)
{
    sb b = {0};
    size_t n = strlen(s), i = 0;
    while (i < n) {
        int behind_ok = i == 0 || !(is_digit(s[i - 1]) || (s[i - 1] >= 'A' && s[i - 1] <= 'Z') || (s[i - 1] >= 'a' && s[i - 1] <= 'z') || s[i - 1] == '_' || s[i - 1] == '.');
        if (behind_ok) {
            size_t w = i;
            while (w < n && is_digit(s[w])) w++;
            if (w < n && s[w] == '.' && w + 1 < n && is_digit(s[w + 1])) {
                size_t f = w + 1;
                while (f < n && is_digit(s[f])) f++;
                if (!(f < n && (is_digit(s[f]) || s[f] == '.'))) {
                    size_t ws0 = i;
                    while (ws0 < w && s[ws0] == '0') ws0++;
                    size_t fe = f;
                    while (fe > w + 1 && s[fe - 1] == '0') fe--;
                    if (ws0 == w) sb_put(&b, "0", 1);
                    else sb_put(&b, s + ws0, w - ws0);
                    if (fe > w + 1) {
                        sb_put(&b, ".", 1);
                        sb_put(&b, s + w + 1, fe - (w + 1));
                    }
                    i = f;
                    continue;
                }
            }
        }
        sb_put(&b, s + i, 1);
        i++;
    }
    char *out = b.s ? sdup(v, b.s, b.n) : sdup(v, "", 0);
    free(b.s);
    return out;
}

/* (formal)(\([^()]*\))(\[[^\[\]]+:(dsp|mag|occ)\])(tail) with the longest formal, tail = base(_\d+)?\([^()]*\) */
static block_t block_identity(vm_t *v, const char *label)
{
    block_t out = {0};
    if (!label) return out;
    int n = (int)strlen(label);
    for (int fe = n - 1; fe >= 1; fe--) {
        if (label[fe] != '(') continue;
        int q = fe + 1;
        while (q < n && label[q] != '(' && label[q] != ')') q++;
        if (q >= n || label[q] != ')') continue;
        int de = q + 1;
        if (de >= n || label[de] != '[') continue;
        int r = de + 1;
        while (r < n && label[r] != '[' && label[r] != ']') r++;
        if (r >= n || label[r] != ']' || r - (de + 1) < 5) continue;
        int kind = -1;
        for (int k = 0; k < 3; k++)
            if (label[r - 4] == ':' && !memcmp(label + r - 3, KIND[k], 3)) kind = k;
        if (kind < 0 || r - 4 <= de + 1) continue;
        int te = r + 1;
        if (te >= n) continue;
        const char *tail = label + te;
        int tn = n - te;
        if (tail[tn - 1] != ')') return out;
        int po = tn - 2;
        while (po >= 0 && tail[po] != '(' && tail[po] != ')') po--;
        if (po < 0 || tail[po] != '(') return out;
        int be = po, dg = po;
        while (dg > 0 && is_digit(tail[dg - 1])) dg--;
        if (dg < po && dg > 0 && tail[dg - 1] == '_') be = dg - 1;
        char *formal = sdup(v, label, (size_t)fe);
        char *pformal = physical_k_label(v, formal);
        sb e = {0}, p = {0};
        sb_cat(&e, pformal);
        sb_put(&e, label + fe, (size_t)(de - fe));
        sb_put(&e, label + de, (size_t)(te - de));
        sb_cat(&p, pformal);
        sb_put(&p, label + de, (size_t)(te - de));
        out.ok = 1;
        out.kind = kind;
        out.exact = canonical_decimals(v, e.s);
        out.physical = canonical_decimals(v, p.s);
        out.direction = sdup(v, label + fe, (size_t)(de - fe));
        out.base = sdup(v, tail, (size_t)be);
        free(e.s);
        free(p.s);
        return out;
    }
    return out;
}

typedef struct {
    int kind;
    char *irrep, *pirrep, *site, *site_irrep;
    long long ordinal;
} ident_t;

static char *substr(vm_t *v, const char *s, int a, int b) { return sdup(v, s + a, (size_t)(b - a)); }

/* head(\([^()]*\)) with a non-empty head */
static int split_final_direction(const char *s, int n, int *head_end)
{
    if (n < 2 || s[n - 1] != ')') return 0;
    int q = n - 2;
    while (q >= 0 && s[q] != '(' && s[q] != ')') q--;
    if (q < 1 || s[q] != '(') return 0;
    *head_end = q;
    return 1;
}

static int digits_only(const char *s, int a, int b)
{
    if (a >= b) return 0;
    for (int i = a; i < b; i++)
        if (!is_digit(s[i])) return 0;
    return 1;
}

static int ident_fail(sb *err, const char *msg, const char *value)
{
    sb_cat(err, msg);
    sb_repr(err, value);
    return -1;
}

static int parse_identity(vm_t *v, const char *label, int expected, ident_t *out, sb *err)
{
    memset(out, 0, sizeof *out);
    if (!label || !*label) {
        sb_cat(err, "mode label must be non-empty text");
        return -1;
    }
    int n = (int)strlen(label);
    const char *o = strchr(label, '[');
    int first_open = o ? (int)(o - label) : -1;
    const char *c = strchr(label + first_open + 1, ']');
    int first_close = c ? (int)(c - label) : -1;
    if (first_open < 1 || first_close < 0) return ident_fail(err, "mode label has no K vector: ", label);
    sb k = {0};
    for (int i = first_open + 1; i < first_close; i++)
        if (!is_space(label[i])) sb_put(&k, label + i, 1);
    const char *tail = label + first_close + 1;
    int tn = n - first_close - 1;
    sb irrep = {0};
    sb_put(&irrep, "[", 1);
    if (k.s) sb_cat(&irrep, k.s);
    sb_put(&irrep, "]", 1);
    free(k.s);
    if (expected == STRAIN) {
        int found = 0;
        for (int p = 0; p < tn && !found; p++) {
            if (tail[p] != '(') continue;
            int q = p + 1;
            while (q < tn && tail[q] != '(' && tail[q] != ')') q++;
            if (q >= tn || tail[q] != ')') continue;
            int r = q + 1;
            if (tn - r < 6 || memcmp(tail + r, "strain", 6)) continue;
            int s = r + 6;
            long long ordinal = 1;
            int sd = -1;
            if (s < tn && tail[s] == '_') {
                int t = s + 1;
                while (t < tn && is_digit(tail[t])) t++;
                if (t > s + 1) {
                    int u = t;
                    if (u < tn && tail[u] == '(') {
                        int w = u + 1;
                        while (w < tn && tail[w] != '(' && tail[w] != ')') w++;
                        if (w == tn - 1 && tail[w] == ')') sd = u, ordinal = strtoll(tail + s + 1, 0, 10);
                    }
                }
            }
            if (sd < 0 && s < tn && tail[s] == '(') {
                int w = s + 1;
                while (w < tn && tail[w] != '(' && tail[w] != ')') w++;
                if (w == tn - 1 && tail[w] == ')') sd = s;
            }
            if (sd < 0) continue;
            found = 1;
            sb_put(&irrep, tail, (size_t)p);
            out->kind = STRAIN;
            out->ordinal = ordinal;
        }
        if (!found) {
            free(irrep.s);
            return ident_fail(err, "cannot parse strain family label ", label);
        }
        out->irrep = sdup(v, irrep.s, irrep.n);
        out->pirrep = physical_k_label(v, out->irrep);
        free(irrep.s);
        return 0;
    }
    const char *so = strchr(tail, '[');
    if (!so) {
        free(irrep.s);
        return ident_fail(err, "mode label has no site identity: ", label);
    }
    int second_open = (int)(so - tail);
    const char *sc = strchr(tail + second_open + 1, ']');
    if (!sc) {
        free(irrep.s);
        return ident_fail(err, "mode label has an unterminated site identity: ", label);
    }
    int second_close = (int)(sc - tail), head_end;
    if (!split_final_direction(tail, second_open, &head_end)) {
        free(irrep.s);
        sb_cat(err, "cannot parse parent irrep from mode label ");
        char *head = substr(v, tail, 0, second_open);
        sb_repr(err, head);
        return -1;
    }
    sb_put(&irrep, tail, (size_t)head_end);
    int colons[3], nc = 0;
    for (int i = second_open + 1; i < second_close; i++)
        if (tail[i] == ':') {
            if (nc < 3) colons[nc] = i;
            nc++;
        }
    if (nc != 2) {
        free(irrep.s);
        return ident_fail(err, "invalid mode site identity in ", label);
    }
    char *atom = substr(v, tail, second_open + 1, colons[0]), *wyckoff = substr(v, tail, colons[0] + 1, colons[1]);
    char *kind = substr(v, tail, colons[1] + 1, second_close);
    int kk = -1;
    for (int i = 0; i < 3; i++)
        if (!strcmp(kind, KIND[i])) kk = i;
    if (kk != expected || kk < 0) {
        free(irrep.s);
        sb_cat(err, "mode label kind ");
        sb_repr(err, kind);
        sb_cat(err, " does not match section ");
        sb_repr(err, KIND[expected]);
        return -1;
    }
    const char *st = tail + second_close + 1;
    int stn = (int)strlen(st);
    out->ordinal = 1;
    if (stn >= 4 && !memcmp(st, "mode", 4) && (stn == 4 || (st[4] == '_' && digits_only(st, 5, stn)))) {
        out->site_irrep = 0;
        if (stn > 4) out->ordinal = strtoll(st + 5, 0, 10);
    } else {
        int he;
        if (!split_final_direction(st, stn, &he)) {
            free(irrep.s);
            return ident_fail(err, "cannot parse site irrep from mode label ", st);
        }
        int be = he;
        for (int p = 0; p < he; p++)
            if (st[p] == '_' && digits_only(st, p + 1, he)) {
                be = p;
                break;
            }
        if (be == 0) {
            free(irrep.s);
            return ident_fail(err, "invalid site-irrep family ", substr(v, st, 0, he));
        }
        if (be < he) out->ordinal = strtoll(st + be + 1, 0, 10);
        out->site_irrep = substr(v, st, 0, be);
    }
    out->kind = kk;
    out->irrep = sdup(v, irrep.s, irrep.n);
    out->pirrep = physical_k_label(v, out->irrep);
    size_t al = strlen(atom), wl = strlen(wyckoff);
    out->site = zalloc(v, al + wl + 2);
    memcpy(out->site, atom, al);
    out->site[al] = ':';
    memcpy(out->site + al + 1, wyckoff, wl);
    free(irrep.s);
    return 0;
}

/* ---- signed permutations ---- */

static int str_eq(const char *a, const char *b) { return a && b && !strcmp(a, b); }

static int block_key_eq(const block_t *a, const block_t *b) { return str_eq(a->exact, b->exact) && str_eq(a->base, b->base); }

typedef struct {
    int nedge;
    int *start, *to;
    int *match_l, *seen, stamp;
} bip_t;

static int bip_augment(bip_t *g, int w)
{
    for (int t = g->start[w]; t < g->start[w + 1]; t++) {
        int l = g->to[t];
        if (g->seen[l] == g->stamp) continue;
        g->seen[l] = g->stamp;
        if (g->match_l[l] < 0 || bip_augment(g, g->match_l[l])) {
            g->match_l[l] = w;
            return 1;
        }
    }
    return 0;
}

/* a complete within-block Web-to-local mode assignment over the given index lists */
static int signed_permutation(payload_t *P, const int *wi, const int *li, int n, int unamb)
{
    if (!n) return 0;
    for (int i = 0; i < n; i++)
        if (!P->wb[wi[i]].ok || !P->lb[li[i]].ok) return 0;
    int *wg = malloc(sizeof(int) * n), *lg = malloc(sizeof(int) * n), *rep = malloc(sizeof(int) * n), ng = 0, ok = 1;
    for (int i = 0; i < n; i++) {
        int g = 0;
        while (g < ng && !block_key_eq(&P->wb[wi[i]], &P->wb[wi[rep[g]]])) g++;
        if (g == ng) rep[ng++] = i;
        wg[i] = g;
    }
    for (int i = 0; i < n && ok; i++) {
        int g = 0;
        while (g < ng && !block_key_eq(&P->lb[li[i]], &P->wb[wi[rep[g]]])) g++;
        if (g == ng) ok = 0;
        lg[i] = g;
    }
    for (int g = 0; g < ng && ok; g++) {
        int nw = 0, nl = 0;
        for (int i = 0; i < n; i++) nw += wg[i] == g, nl += lg[i] == g;
        if (nw != nl) {
            ok = 0;
            break;
        }
        int *W = malloc(sizeof(int) * nw), *L = malloc(sizeof(int) * nl), a = 0, b = 0;
        for (int i = 0; i < n; i++) {
            if (wg[i] == g) W[a++] = wi[i];
            if (lg[i] == g) L[b++] = li[i];
        }
        bip_t G = {0};
        G.start = malloc(sizeof(int) * (nw + 1));
        G.to = malloc(sizeof(int) * (nw * nl ? nw * nl : 1));
        int e = 0;
        for (int x = 0; x < nw; x++) {
            G.start[x] = e;
            for (int y = 0; y < nl; y++) {
                if (pair_status(P, W[x], L[y], 0) > S_SIGN) continue;
                if (unamb && !position_vectors_unambiguous(P, W[x], L[y])) continue;
                G.to[e++] = y;
            }
        }
        G.start[nw] = e;
        G.match_l = malloc(sizeof(int) * nl);
        G.seen = malloc(sizeof(int) * nl);
        for (int y = 0; y < nl; y++) G.match_l[y] = -1, G.seen[y] = -1;
        for (int x = 0; x < nw && ok; x++) {
            G.stamp = x;
            ok = bip_augment(&G, x);
        }
        free(G.start);
        free(G.to);
        free(G.match_l);
        free(G.seen);
        free(W);
        free(L);
    }
    free(wg);
    free(lg);
    free(rep);
    return ok;
}

/* ---- mode-span proofs ---- */

static void householder_qr(double *A, int m, int d, double *B, int nb)
{
    /* A (m x d) and B (m x nb) column-major; A becomes R in its top d rows, B becomes Q^T B */
    for (int k = 0; k < d; k++) {
        double norm = 0;
        for (int i = k; i < m; i++) norm = hypot(norm, A[k * m + i]);
        if (norm == 0) continue;
        double alpha = A[k * m + k] > 0 ? -norm : norm;
        double *vk = malloc(sizeof(double) * (m - k));
        for (int i = k; i < m; i++) vk[i - k] = A[k * m + i];
        vk[0] -= alpha;
        double vn = 0;
        for (int i = 0; i < m - k; i++) vn += vk[i] * vk[i];
        if (vn > 0) {
            for (int j = k; j < d; j++) {
                double s = 0;
                for (int i = k; i < m; i++) s += vk[i - k] * A[j * m + i];
                s = 2 * s / vn;
                for (int i = k; i < m; i++) A[j * m + i] -= s * vk[i - k];
            }
            for (int j = 0; j < nb; j++) {
                double s = 0;
                for (int i = k; i < m; i++) s += vk[i - k] * B[j * m + i];
                s = 2 * s / vn;
                for (int i = k; i < m; i++) B[j * m + i] -= s * vk[i - k];
            }
        }
        free(vk);
    }
}

/* singular values of a d x d column-major matrix, descending (one-sided Jacobi) */
static void singular_values(const double *M, int d, double *s)
{
    double *U = malloc(sizeof(double) * d * d);
    memcpy(U, M, sizeof(double) * d * d);
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0;
        for (int p = 0; p < d; p++)
            for (int q = p + 1; q < d; q++) {
                double a = 0, b = 0, c = 0;
                for (int i = 0; i < d; i++) {
                    a += U[p * d + i] * U[p * d + i];
                    b += U[q * d + i] * U[q * d + i];
                    c += U[p * d + i] * U[q * d + i];
                }
                if (c == 0 || fabs(c) <= DBL_EPSILON * sqrt(a * b)) continue;
                off = fmax(off, fabs(c) / sqrt(a * b));
                double zeta = (b - a) / (2 * c), t = (zeta >= 0 ? 1 : -1) / (fabs(zeta) + sqrt(1 + zeta * zeta));
                double cs = 1 / sqrt(1 + t * t), sn = cs * t;
                for (int i = 0; i < d; i++) {
                    double x = U[p * d + i], y = U[q * d + i];
                    U[p * d + i] = cs * x - sn * y;
                    U[q * d + i] = sn * x + cs * y;
                }
            }
        if (off <= DBL_EPSILON) break;
    }
    for (int j = 0; j < d; j++) {
        double n = 0;
        for (int i = 0; i < d; i++) n = hypot(n, U[j * d + i]);
        s[j] = n;
    }
    free(U);
    for (int i = 0; i < d; i++)
        for (int j = i + 1; j < d; j++)
            if (s[j] > s[i]) {
                double t = s[i];
                s[i] = s[j];
                s[j] = t;
            }
}

static int solve_square(double *M, double *rhs, int n)
{
    /* M row-major n x n, partial pivoting; rhs overwritten by the solution */
    for (int k = 0; k < n; k++) {
        int p = k;
        for (int i = k + 1; i < n; i++)
            if (fabs(M[i * n + k]) > fabs(M[p * n + k])) p = i;
        if (M[p * n + k] == 0) return -1;
        if (p != k) {
            for (int j = 0; j < n; j++) {
                double t = M[k * n + j];
                M[k * n + j] = M[p * n + j];
                M[p * n + j] = t;
            }
            double t = rhs[k];
            rhs[k] = rhs[p];
            rhs[p] = t;
        }
        for (int i = k + 1; i < n; i++) {
            double f = M[i * n + k] / M[k * n + k];
            if (f == 0) continue;
            for (int j = k; j < n; j++) M[i * n + j] -= f * M[k * n + j];
            rhs[i] -= f * rhs[k];
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        double s = rhs[k];
        for (int j = k + 1; j < n; j++) s -= M[k * n + j] * rhs[j];
        rhs[k] = s / M[k * n + k];
    }
    return 0;
}

/* B x = rhs for the (d+1)-row basis whose columns are col(basis[k]); transposed solves B^T x = rhs */
static int basis_solve(const double *A, int m, int d, const int *basis, const double *rhs, double *x, int transposed)
{
    int rows = d + 1;
    double *M = malloc(sizeof(double) * rows * rows);
    for (int k = 0; k < rows; k++) {
        int j = basis[k];
        for (int r = 0; r < rows; r++) {
            double v = r == d ? 1 : j < m ? A[r * m + j] : j < 2 * m ? -A[r * m + j - m] : 0;
            if (transposed) M[k * rows + r] = v;
            else M[r * rows + k] = v;
        }
    }
    memcpy(x, rhs, sizeof(double) * rows);
    int bad = solve_square(M, x, rows);
    free(M);
    return bad;
}

/* min_q max_i |(A q - b)_i| by the revised simplex method on its dual, max b.(u - v) subject to A^T (u - v) = 0 and
 * sum(u + v) + s = 1, refactoring the basis every step; the dual values of the optimal basis are (q, max residual) */
static int chebyshev(const double *A, int m, int d, const double *b, double *q)
{
    int rows = d + 1, n = 2 * m + 1, ok = 1;
    int *basis = malloc(sizeof(int) * rows);
    double *x = malloc(sizeof(double) * rows), *y = malloc(sizeof(double) * rows), *w = malloc(sizeof(double) * rows);
    double *cb = malloc(sizeof(double) * rows), *rhs = calloc(rows, sizeof(double)), *col = malloc(sizeof(double) * rows);
    double *R = malloc(sizeof(double) * d * m);
    char *in = calloc(n, 1);
    memcpy(R, A, sizeof(double) * d * m);
    /* d rows of A spanning its column space, by elimination with the largest pivot, then the slack */
    for (int k = 0; k < d && ok; k++) {
        int best = -1;
        for (int i = 0; i < m; i++)
            if (!in[i] && (best < 0 || fabs(R[k * m + i]) > fabs(R[k * m + best]))) best = i;
        if (best < 0 || R[k * m + best] == 0) {
            ok = 0;
            break;
        }
        for (int r = k + 1; r < d; r++) {
            double f = R[r * m + best] / R[k * m + best];
            if (f != 0)
                for (int i = 0; i < m; i++) R[r * m + i] -= f * R[k * m + i];
        }
        basis[k] = best;
        in[best] = 1;
    }
    basis[d] = 2 * m;
    in[2 * m] = 1;
    rhs[d] = 1;
    double *g = malloc(sizeof(double) * m);
    for (int iter = 0, degenerate = 0; ok; iter++) {
        if (iter == 100000 || basis_solve(A, m, d, basis, rhs, x, 0)) {
            ok = 0;
            break;
        }
        for (int k = 0; k < rows; k++) {
            int j = basis[k];
            cb[k] = j < m ? b[j] : j < 2 * m ? -b[j - m] : 0;
        }
        if (basis_solve(A, m, d, basis, cb, y, 1)) {
            ok = 0;
            break;
        }
        for (int i = 0; i < m; i++) {
            double acc = 0;
            for (int k = 0; k < d; k++) acc += A[k * m + i] * y[k];
            g[i] = acc;
        }
        int e = -1;
        double de = 1e-12;
        for (int j = 0; j < n; j++) {
            if (in[j]) continue;
            double reduced = j < m ? b[j] - g[j] - y[d] : j < 2 * m ? -b[j - m] + g[j - m] - y[d] : -y[d];
            if (reduced > de || (degenerate >= 50 && reduced > 1e-12 && e < 0)) {
                e = j;
                de = reduced;
                if (degenerate >= 50) break;
            }
        }
        if (e < 0) break;
        for (int r = 0; r < d; r++) col[r] = e < m ? A[r * m + e] : e < 2 * m ? -A[r * m + e - m] : 0;
        col[d] = 1;
        if (basis_solve(A, m, d, basis, col, w, 0)) {
            ok = 0;
            break;
        }
        int r = -1;
        double best = 0;
        for (int k = 0; k < rows; k++)
            if (w[k] > 1e-12) {
                double ratio = fmax(x[k], 0) / w[k];
                if (r < 0 || ratio < best - 1e-15 || (fabs(ratio - best) <= 1e-15 && basis[k] < basis[r])) r = k, best = ratio;
            }
        if (r < 0) {
            ok = 0;
            break;
        }
        degenerate = best > 1e-15 ? 0 : degenerate + 1;
        in[basis[r]] = 0;
        basis[r] = e;
        in[e] = 1;
    }
    if (ok)
        for (int k = 0; k < d; k++) q[k] = y[k];
    free(basis);
    free(x);
    free(y);
    free(w);
    free(cb);
    free(rhs);
    free(col);
    free(R);
    free(in);
    free(g);
    return ok;
}

static double visible_bound(double tol)
{
    double half = nextafter(0.5 * pow(10.0, -4), INFINITY);
    return tol > half ? tol : half;
}

/* the rows of other in the order of ref, when every ref row has exactly one candidate */
static int *reference_order(payload_t *P, def_t *ref, def_t *other)
{
    int owned;
    cand_t *c = candidates(P->v, ref, other, &owned);
    int n = ref->nrows, *a = malloc(sizeof(int) * (n ? n : 1)), ok = 1;
    for (int r = 0; r < n && ok; r++) {
        if (c->start[r + 1] - c->start[r] != 1) ok = 0;
        else a[r] = c->idx[c->start[r]];
    }
    if (ok) {
        char *seen = calloc(other->nrows + 1, 1);
        int distinct = 0;
        for (int r = 0; r < n; r++)
            if (!seen[a[r]]) seen[a[r]] = 1, distinct++;
        free(seen);
        if (distinct != other->nrows) ok = 0;
    }
    if (owned) cand_free(c);
    if (!ok) {
        free(a);
        return 0;
    }
    return a;
}

static int span_block(payload_t *P, const int *wi, const int *li, int dim, proof_t *pr)
{
    def_t *ref = (def_t *)&P->W[wi[0]];
    int nref = ref->nrows, m = 3 * nref;
    if (!nref) return 0;
    double *web = malloc(sizeof(double) * m * dim), *loc = malloc(sizeof(double) * m * dim);
    int ok = 1;
    for (int side = 0; side < 2 && ok; side++)
        for (int j = 0; j < dim && ok; j++) {
            def_t *mm = side ? (def_t *)&P->L[li[j]] : (def_t *)&P->W[wi[j]];
            int *ord = reference_order(P, ref, mm);
            if (!ord) {
                ok = 0;
                break;
            }
            double *col = (side ? loc : web) + (size_t)j * m;
            for (int r = 0; r < nref && ok; r++) {
                const row_t *row = &mm->rows[ord[r]];
                if (row->nd != 3) ok = 0;
                for (int a = 0; a < 3; a++) {
                    col[3 * r + a] = row->d[a];
                    if (!isfinite(row->d[a])) ok = 0;
                }
            }
            free(ord);
        }
    double tol = P->v->tol, bound = visible_bound(tol);
    double *R = 0, *QB = 0, *X = 0, *res = 0, *s = 0;
    if (ok) {
        R = malloc(sizeof(double) * m * dim);
        QB = malloc(sizeof(double) * m * dim);
        memcpy(R, loc, sizeof(double) * m * dim);
        memcpy(QB, web, sizeof(double) * m * dim);
        householder_qr(R, m, dim, QB, dim);
        double *Rs = malloc(sizeof(double) * dim * dim);
        for (int j = 0; j < dim; j++)
            for (int i = 0; i < dim; i++) Rs[j * dim + i] = i <= j && i < m ? R[j * m + i] : 0;
        s = malloc(sizeof(double) * dim);
        singular_values(Rs, dim, s);
        double cutoff = DBL_EPSILON * (m > dim ? m : dim) * s[0];
        int rank = 0;
        for (int i = 0; i < dim; i++) rank += s[i] > cutoff;
        if (rank != dim || !(s[dim - 1] > 0)) ok = 0;
        double lc = ok ? s[0] / s[dim - 1] : 0;
        if (ok && (!isfinite(lc) || lc > 1.0 / tol)) ok = 0;
        X = malloc(sizeof(double) * dim * dim);
        for (int c = 0; c < dim && ok; c++)
            for (int i = dim - 1; i >= 0; i--) {
                double acc = QB[c * m + i];
                for (int j = i + 1; j < dim; j++) acc -= Rs[j * dim + i] * X[c * dim + j];
                X[c * dim + i] = acc / Rs[i * dim + i];
            }
        free(Rs);
        res = malloc(sizeof(double) * m * dim);
        double residual = 0;
        for (int c = 0; c < dim && ok; c++)
            for (int i = 0; i < m; i++) {
                double acc = 0;
                for (int j = 0; j < dim; j++) acc += loc[j * m + i] * X[c * dim + j];
                res[c * m + i] = acc - web[c * m + i];
                if (!isfinite(res[c * m + i])) ok = 0;
                residual = fmax(residual, fabs(res[c * m + i]));
            }
        if (ok && residual > bound) {
            double box = nextafter(sqrt((double)m) * bound, INFINITY);
            for (int c = 0; c < dim && ok; c++) {
                double nn = 0;
                for (int i = 0; i < m; i++) nn = hypot(nn, res[c * m + i]);
                if (nn > box) ok = 0;
            }
            for (int c = 0; c < dim && ok; c++) {
                double cm = 0;
                for (int i = 0; i < m; i++) cm = fmax(cm, fabs(res[c * m + i]));
                if (cm > bound && !chebyshev(loc, m, dim, web + (size_t)c * m, X + (size_t)c * dim)) ok = 0;
            }
            residual = 0;
            for (int c = 0; c < dim && ok; c++)
                for (int i = 0; i < m; i++) {
                    double acc = 0;
                    for (int j = 0; j < dim; j++) acc += loc[j * m + i] * X[c * dim + j];
                    double r = acc - web[c * m + i];
                    if (!isfinite(r)) ok = 0;
                    residual = fmax(residual, fabs(r));
                }
            if (ok && residual > bound) ok = 0;
            if (ok) pr->minimax_count++;
        }
        if (ok) {
            double *sq = malloc(sizeof(double) * dim);
            singular_values(X, dim, sq);
            double qc = sq[dim - 1] > 0 ? sq[0] / sq[dim - 1] : INFINITY;
            free(sq);
            double cond = fmax(qc, lc);
            if (!isfinite(cond) || cond > 1.0 / tol) ok = 0;
            else {
                pr->block_count++;
                pr->family_count++;
                pr->max_condition = fmax(pr->max_condition, cond);
                pr->max_residual = fmax(pr->max_residual, residual);
                pr->ratio = fmax(pr->ratio, residual / bound);
            }
        }
    }
    free(web);
    free(loc);
    free(R);
    free(QB);
    free(X);
    free(res);
    free(s);
    return ok;
}

typedef struct {
    char **key;
    int n, cap;
    int **members, *count;
} groups_t;

static int group_index(groups_t *g, const char *key, int add)
{
    for (int i = 0; i < g->n; i++)
        if (!strcmp(g->key[i], key)) return i;
    if (!add) return -1;
    if (g->n == g->cap) {
        g->cap = g->cap ? 2 * g->cap : 16;
        g->key = realloc(g->key, sizeof(char *) * g->cap);
        g->members = realloc(g->members, sizeof(int *) * g->cap);
        g->count = realloc(g->count, sizeof(int) * g->cap);
    }
    g->key[g->n] = strdup(key);
    g->members[g->n] = 0;
    g->count[g->n] = 0;
    return g->n++;
}

static void group_add(groups_t *g, int i, int member)
{
    g->members[i] = realloc(g->members[i], sizeof(int) * (g->count[i] + 1));
    g->members[i][g->count[i]++] = member;
}

static void groups_free(groups_t *g)
{
    for (int i = 0; i < g->n; i++) {
        free(g->key[i]);
        free(g->members[i]);
    }
    free(g->key);
    free(g->members);
    free(g->count);
}

static char *pair_key(const char *a, const char *b)
{
    size_t na = strlen(a), nb = b ? strlen(b) : 0;
    char *k = malloc(na + nb + 2);
    memcpy(k, a, na);
    k[na] = b ? '\x01' : 0;
    if (b) {
        memcpy(k + na + 1, b, nb);
        k[na + 1 + nb] = 0;
    }
    return k;
}

static int same_direction_sets(groups_t *a, int ia, groups_t *b, int ib, const block_t *ab, const block_t *bb)
{
    /* the direction sets of one physical key on both sides */
    for (int s = 0; s < 2; s++) {
        groups_t *x = s ? b : a, *y = s ? a : b;
        int ix = s ? ib : ia, iy = s ? ia : ib;
        const block_t *xb = s ? bb : ab, *yb = s ? ab : bb;
        for (int t = 0; t < x->count[ix]; t++) {
            const char *dir = xb[x->members[ix][t]].direction;
            int found = 0;
            for (int u = 0; u < y->count[iy] && !found; u++) found = !strcmp(yb[y->members[iy][u]].direction, dir);
            if (!found) return 0;
        }
    }
    return 1;
}

static int family_span(payload_t *P, int direction_gauge, int site_basis, proof_t *pr)
{
    int n = P->nw;
    memset(pr, 0, sizeof *pr);
    if (P->kind == STRAIN || P->nw != P->nl || !n || !isfinite(P->v->tol) || P->v->tol <= 0) return 0;
    for (int i = 0; i < n; i++)
        if (!P->wb[i].ok || P->wb[i].kind != P->kind || !P->lb[i].ok || P->lb[i].kind != P->kind) return 0;
    groups_t wd = {0}, ld = {0};
    for (int i = 0; i < n; i++) {
        char *k = pair_key(P->wb[i].physical, P->wb[i].base);
        group_add(&wd, group_index(&wd, k, 1), i);
        free(k);
        k = pair_key(P->lb[i].physical, P->lb[i].base);
        group_add(&ld, group_index(&ld, k, 1), i);
        free(k);
    }
    int ok = wd.n == ld.n;
    for (int i = 0; i < wd.n && ok; i++) ok = group_index(&ld, wd.key[i], 0) >= 0;
    char *gauge = calloc(wd.n + 1, 1);
    int ngauge = 0;
    for (int i = 0; i < wd.n && ok; i++)
        if (direction_gauge && !same_direction_sets(&wd, i, &ld, group_index(&ld, wd.key[i], 0), P->wb, P->lb)) gauge[i] = 1, ngauge++;
    groups_t wg = {0}, lg = {0};
    for (int side = 0; side < 2 && ok; side++) {
        const block_t *bl = side ? P->lb : P->wb;
        groups_t *g = side ? &lg : &wg;
        for (int i = 0; i < n; i++) {
            char *pk = pair_key(bl[i].physical, bl[i].base);
            int gi = group_index(&wd, pk, 0);
            free(pk);
            int phys = gi >= 0 && gauge[gi];
            char *k = pair_key(phys ? bl[i].physical : bl[i].exact, site_basis ? 0 : bl[i].base);
            group_add(g, group_index(g, k, 1), i);
            free(k);
        }
    }
    if (ok) {
        ok = wg.n == lg.n;
        for (int i = 0; i < wg.n && ok; i++) ok = group_index(&lg, wg.key[i], 0) >= 0;
    }
    if (ok && direction_gauge && !ngauge) ok = 0;
    for (int i = 0; i < wg.n && ok; i++) {
        int j = group_index(&lg, wg.key[i], 0), dim = wg.count[i];
        if (!dim || lg.count[j] != dim) {
            ok = 0;
            break;
        }
        if (signed_permutation(P, wg.members[i], lg.members[j], dim, 0)) {
            pr->block_count++;
            pr->signed_count++;
            continue;
        }
        if (!span_block(P, wg.members[i], lg.members[j], dim, pr)) ok = 0;
    }
    pr->direction_count = ngauge;
    pr->bound = visible_bound(P->v->tol);
    free(gauge);
    groups_free(&wd);
    groups_free(&ld);
    groups_free(&wg);
    groups_free(&lg);
    return ok;
}

static int site_scope_broader(const block_t *b, int n, int kind)
{
    for (int i = 0; i < n; i++)
        if (!b[i].ok || b[i].kind != kind) return 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
            if (!strcmp(b[i].exact, b[j].exact) && strcmp(b[i].base, b[j].base)) return 1;
    return 0;
}

/* ---- per-section payloads ---- */

static void vector_payload(vm_t *v, payload_t *P, const def_t *W, int nw, const def_t *L, int nl, int kind, int unamb)
{
    memset(P, 0, sizeof *P);
    P->v = v;
    P->W = W;
    P->L = L;
    P->nw = nw;
    P->nl = nl;
    P->kind = kind;
    P->unamb = unamb;
    P->count_match = nw == nl;
    P->fact = zalloc(v, (size_t)nw * nl + 1);
    P->unambiguous = zalloc(v, (size_t)nw * nl + 1);
    memset(P->fact, -1, (size_t)nw * nl);
    memset(P->unambiguous, -1, (size_t)nw * nl);
    P->wb = zalloc(v, sizeof(block_t) * (nw + 1));
    P->lb = zalloc(v, sizeof(block_t) * (nl + 1));
    for (int i = 0; i < nw; i++) P->wb[i] = block_identity(v, W[i].label);
    for (int i = 0; i < nl; i++) P->lb[i] = block_identity(v, L[i].label);
    int nz = nw < nl ? nw : nl;
    P->mst = zalloc(v, sizeof(int) * (nz + 1));
    P->mstrict = zalloc(v, sizeof(int) * (nz + 1));
    P->label_match = zalloc(v, sizeof(int) * (nz + 1));
    P->nf_match = zalloc(v, sizeof(int) * (nz + 1));
    P->assign = zalloc(v, sizeof(int *) * (nz + 1));
    unsigned hard = 0, strict = 0;
    for (int i = 0; i < nz; i++) {
        int code = pair_status(P, i, i, &P->assign[i]);
        P->mst[i] = code;
        P->label_match[i] = str_eq(W[i].label, L[i].label);
        P->nf_match[i] = optional_close(W[i].has_nf, W[i].nf, L[i].has_nf, L[i].nf, v->tol);
        P->mstrict[i] = code != S_OK ? code : !P->label_match[i] ? S_LABEL : !P->nf_match[i] ? S_NORMFACTOR : S_OK;
        hard |= 1u << code;
        strict |= 1u << P->mstrict[i];
    }
    int signed_perm = 0, family = 0;
    if (P->count_match && (hard & ~((1u << S_OK) | (1u << S_SIGN)))) {
        int *idx = malloc(sizeof(int) * (nw + 1));
        for (int i = 0; i < nw; i++) idx[i] = i;
        signed_perm = signed_permutation(P, idx, idx, nw, unamb);
        free(idx);
        if (!signed_perm) {
            int scope = site_scope_broader(P->wb, nw, kind) && site_scope_broader(P->lb, nl, kind);
            static const int variants[4][2] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
            for (int t = 0; t < 4 && !family; t++) {
                if (variants[t][1] && !scope) continue;
                family = family_span(P, variants[t][0], variants[t][1], &P->proof);
            }
            P->has_proof = family;
        }
    }
    if (!P->count_match) P->status = S_COUNT;
    else if (!(hard & ~(1u << S_OK))) P->status = S_OK;
    else if (!(hard & ~((1u << S_OK) | (1u << S_SIGN)))) P->status = S_SIGN;
    else if (signed_perm) P->status = S_SIGNED_PERMUTATION;
    else if (family) P->status = S_FAMILY_BASIS;
    else if (hard & (1u << S_POSITION)) P->status = S_POSITION;
    else P->status = S_VECTOR;
    if (!P->count_match) P->strict = S_COUNT;
    else if (P->status == S_SIGNED_PERMUTATION) P->strict = S_SIGNED_PERMUTATION;
    else if (P->status == S_FAMILY_BASIS) P->strict = S_FAMILY_BASIS;
    else if (!(strict & ~(1u << S_OK))) P->strict = S_OK;
    else {
        static const int order[] = {S_POSITION, S_VECTOR, S_SIGN, S_LABEL, S_NORMFACTOR};
        P->strict = S_FAIL;
        for (int i = 0; i < 5; i++)
            if (strict & (1u << order[i])) {
                P->strict = order[i];
                break;
            }
    }
}

static void strain_payload(vm_t *v)
{
    strain_payload_t *S = &v->SP;
    int nw = v->nwstrain, nl = v->nlstrain, nz = nw < nl ? nw : nl;
    S->nw = nw;
    S->nl = nl;
    S->count_match = nw == nl;
    S->mst = zalloc(v, sizeof(int) * (nz + 1));
    S->mstrict = zalloc(v, sizeof(int) * (nz + 1));
    unsigned hard = 0, strict = 0;
    for (int i = 0; i < nz; i++) {
        const strain_t *w = &v->wstrain[i], *l = &v->lstrain[i];
        int same = vec_close(w->c, 6, l->c, l->nc, v->vtol, 1.0), sign = vec_close(w->c, 6, l->c, l->nc, v->vtol, -1.0);
        int lm = str_eq(w->label, l->label), nm = optional_close(w->has_nf, w->nf, l->has_nf, l->nf, v->tol);
        S->mst[i] = same ? S_OK : sign ? S_SIGN : S_TENSOR;
        S->mstrict[i] = S->mst[i] != S_OK ? S->mst[i] : !lm ? S_LABEL : !nm ? S_NORMFACTOR : S_OK;
        hard |= 1u << S->mst[i];
        strict |= 1u << S->mstrict[i];
    }
    if (!S->count_match) S->status = S->strict = S_COUNT;
    else {
        S->status = !(hard & ~(1u << S_OK)) ? S_OK : !(hard & ~((1u << S_OK) | (1u << S_SIGN))) ? S_SIGN : S_TENSOR;
        if (!(strict & ~(1u << S_OK))) S->strict = S_OK;
        else {
            static const int order[] = {S_TENSOR, S_SIGN, S_LABEL, S_NORMFACTOR};
            S->strict = S_FAIL;
            for (int i = 0; i < 4; i++)
                if (strict & (1u << order[i])) {
                    S->strict = order[i];
                    break;
                }
        }
    }
}

static int accepted(int status) { return status == S_OK || status == S_SIGN || status == S_SIGNED_PERMUTATION || status == S_FAMILY_BASIS; }

/* the definitions expressed in the request child frame (x U + shift, d U, unit max per column), adopted only when they close */
static int transport_try(vm_t *v, int kind, const int *U, const double *shift)
{
    int n = v->nlocal[kind];
    const def_t *L = v->local[kind];
    if (kind != DSP && kind != MAG) return 0;
    if (n != v->nweb[kind] || !n) return 0;
    def_t *T = zalloc(v, sizeof(def_t) * n);
    for (int i = 0; i < n; i++) {
        const def_t *s = &L[i];
        if (!s->label || !*s->label || !s->nrows) return 0;
        def_t *t = &T[i];
        *t = *s;
        t->pos_id = -1;
        t->rows = zalloc(v, sizeof(row_t) * s->nrows);
        double scale = 0;
        for (int r = 0; r < s->nrows; r++) {
            const row_t *a = &s->rows[r];
            row_t *b = &t->rows[r];
            if (a->nx != 3 || a->nd != 3) return 0;
            for (int k = 0; k < 3; k++)
                if (!isfinite(a->xyz[k]) || !isfinite(a->d[k])) return 0;
            *b = *a;
            for (int col = 0; col < 3; col++) {
                double x = 0, d = 0;
                for (int row = 0; row < 3; row++) {
                    x += a->xyz[row] * (double)U[3 * row + col];
                    d += a->d[row] * (double)U[3 * row + col];
                }
                b->xyz[col] = pymod1(x + shift[col]);
                b->d[col] = d;
                if (!isfinite(b->xyz[col]) || !isfinite(b->d[col])) return 0;
                scale = fmax(scale, fabs(d));
            }
        }
        if (!isfinite(scale) || scale <= 0) return 0;
        if (t->has_nf) {
            t->nf = s->nf * scale;
            if (!isfinite(t->nf)) return 0;
        }
        for (int r = 0; r < s->nrows; r++)
            for (int k = 0; k < 3; k++) t->rows[r].d[k] /= scale;
    }
    payload_t Q;
    vector_payload(v, &Q, v->web[kind], n, T, n, kind, 1);
    if (!Q.count_match) return 0;
    if (Q.status == S_OK || Q.status == S_SIGN) {
        for (int i = 0; i < n; i++)
            if (!position_vectors_unambiguous(&Q, i, i)) return 0;
    } else if (Q.status != S_SIGNED_PERMUTATION && !(Q.status == S_FAMILY_BASIS && Q.has_proof)) return 0;
    Q.validated = 1;
    v->P[kind] = Q;
    return 1;
}

/* a failed try leaves nothing behind (a domain search tries many frames) */
static int transport(vm_t *v, int kind, const int *U, const double *shift)
{
    int mark = v->nowned, ok = transport_try(v, kind, U, shift);
    if (!ok) {
        for (int i = mark; i < v->nowned; i++) free(v->owned[i]);
        v->nowned = mark;
    }
    return ok;
}

/* ---- the mode_definition layers ---- */

typedef struct {
    char **s;
    int n, cap;
} slist;

static void sl_add(slist *l, char *s)
{
    if (l->n == l->cap) {
        l->cap = l->cap ? 2 * l->cap : 32;
        l->s = realloc(l->s, sizeof(char *) * l->cap);
    }
    l->s[l->n++] = s;
}

static void sl_sort_unique(slist *l)
{
    qsort(l->s, l->n, sizeof(char *), cmp_str);
    int k = 0;
    for (int i = 0; i < l->n; i++)
        if (!k || strcmp(l->s[i], l->s[k - 1])) l->s[k++] = l->s[i];
    l->n = k;
}

static int sl_has(const slist *l, const char *s)
{
    int lo = 0, hi = l->n - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2, c = strcmp(l->s[mid], s);
        if (!c) return 1;
        if (c < 0) lo = mid + 1;
        else hi = mid - 1;
    }
    return 0;
}

static void sb_list(sb *b, const slist *l)
{
    sb_put(b, "[", 1);
    for (int i = 0; i < l->n; i++) {
        if (i) sb_cat(b, ", ");
        sb_json(b, l->s[i]);
    }
    sb_put(b, "]", 1);
}

/* {"pass": ..., "missing": [...], "extra": [...]} of two token sets; returns pass */
static int diff_layer(slist *w, slist *l, slist *missing, slist *extra)
{
    sl_sort_unique(w);
    sl_sort_unique(l);
    for (int i = 0; i < w->n; i++)
        if (!sl_has(l, w->s[i])) sl_add(missing, w->s[i]);
    for (int i = 0; i < l->n; i++)
        if (!sl_has(w, l->s[i])) sl_add(extra, l->s[i]);
    return !missing->n && !extra->n;
}

static void sb_diff(sb *b, int pass, const slist *missing, const slist *extra)
{
    sb_fmt(b, "{\"pass\": %s", pass ? "true" : "false");
    if (missing->n) {
        sb_cat(b, ", \"missing\": ");
        sb_list(b, missing);
    }
    if (extra->n) {
        sb_cat(b, ", \"extra\": ");
        sb_list(b, extra);
    }
}

static char *token(vm_t *v, const ident_t *id, int level, int physical)
{
    sb b = {0};
    sb_cat(&b, KIND[id->kind]);
    sb_put(&b, ":", 1);
    sb_cat(&b, physical ? id->pirrep : id->irrep);
    if (level >= 1) {
        sb_put(&b, ":", 1);
        sb_cat(&b, id->site ? id->site : "-");
    }
    if (level >= 2) {
        sb_put(&b, ":", 1);
        sb_cat(&b, id->site_irrep ? id->site_irrep : "-");
    }
    char *o = sdup(v, b.s, b.n);
    free(b.s);
    return o;
}

typedef struct {
    ident_t *id;
    int n, cap;
} ilist;

static void il_add(ilist *l, ident_t id)
{
    if (l->n == l->cap) {
        l->cap = l->cap ? 2 * l->cap : 32;
        l->id = realloc(l->id, sizeof(ident_t) * l->cap);
    }
    l->id[l->n++] = id;
}

static const char *web_label(vm_t *v, int kind, int i) { return kind == STRAIN ? v->wstrain[i].label : v->P[kind].W[i].label; }

static const char *local_label(vm_t *v, int kind, int i) { return kind == STRAIN ? v->lstrain[i].label : v->P[kind].L[i].label; }

static int web_count(vm_t *v, int kind) { return kind == STRAIN ? v->SP.nw : v->P[kind].nw; }

static int local_count(vm_t *v, int kind) { return kind == STRAIN ? v->SP.nl : v->P[kind].nl; }

static int kind_status(vm_t *v, int kind) { return kind == STRAIN ? v->SP.status : v->P[kind].status; }

static int kind_count_match(vm_t *v, int kind) { return kind == STRAIN ? v->SP.count_match : v->P[kind].count_match; }

/* the labels of the paired modes (missing ones skipped) and of the unpaired tails; strain has no tails */
static int identities(vm_t *v, int local, ilist *out, sb *err)
{
    for (int k = 0; k < NKIND; k++) {
        int nw = web_count(v, k), nl = local_count(v, k), nz = nw < nl ? nw : nl, n = k == STRAIN ? nz : local ? nl : nw;
        for (int i = 0; i < n; i++) {
            const char *label = local ? local_label(v, k, i) : web_label(v, k, i);
            if (!label && local && i >= nz) label = "";
            if (!label) continue;
            ident_t id;
            if (parse_identity(v, label, k, &id, err)) return -1;
            il_add(out, id);
        }
    }
    return 0;
}

static int presentation_proof(const payload_t *P)
{
    if (!P->count_match || !P->nw) return 0;
    if (P->status == S_SIGNED_PERMUTATION) return P->strict == S_SIGNED_PERMUTATION;
    if (P->status != S_FAMILY_BASIS || P->strict != S_FAMILY_BASIS || !P->has_proof) return 0;
    const proof_t *p = &P->proof;
    if (p->block_count <= 0 || p->family_count <= 0 || p->signed_count + p->family_count != p->block_count || p->direction_count > p->block_count) return 0;
    double vals[4] = {p->max_condition, p->max_residual, p->bound, p->ratio};
    for (int i = 0; i < 4; i++)
        if (!isfinite(vals[i]) || vals[i] < 0) return 0;
    if (p->bound <= 0) return 0;
    double expect = p->max_residual / p->bound;
    return p->max_residual <= p->bound && p->ratio <= 1.0 && fabs(p->ratio - expect) <= fmax(1e-12 * fmax(fabs(p->ratio), fabs(expect)), 1e-15);
}

/* the length of an atom label without its _\d+ orbit ordinal */
static size_t without_ordinal(const char *s)
{
    size_t len = strlen(s), e = len;
    while (e > 0 && is_digit(s[e - 1])) e--;
    return e < len && e > 0 && s[e - 1] == '_' ? e - 1 : len;
}

/* the parent atoms of two row labels, a missing label standing for the mode's site atom */
static int parents_equal(const char *a, const char *site_a, const char *b, const char *site_b)
{
    const char *x = a ? a : site_a, *y = b ? b : site_b;
    if (!x || !y) return 0;
    size_t nx = a ? without_ordinal(a) : strlen(site_a), ny = b ? without_ordinal(b) : strlen(site_b);
    return nx == ny && !memcmp(x, y, nx);
}

static char *site_parent(vm_t *v, const ident_t *id)
{
    if (!id->site) return 0;
    const char *c = strchr(id->site, ':');
    return sdup(v, id->site, c ? (size_t)(c - id->site) : strlen(id->site));
}

static int row_atom_layer(vm_t *v, sb *b, sb *err)
{
    long exact = 0, gauge = 0, real = 0, positions = 0, presentation = 0;
    for (int k = 0; k < 3; k++) {
        payload_t *P = &v->P[k];
        int proof = presentation_proof(P), nz = P->nw < P->nl ? P->nw : P->nl;
        for (int i = 0; i < nz; i++) {
            ident_t wid, lid;
            if (parse_identity(v, P->W[i].label, k, &wid, err) || parse_identity(v, P->L[i].label, k, &lid, err)) return -1;
            if (P->W[i].nrows != P->L[i].nrows || P->mst[i] == S_POSITION) {
                if (proof) presentation++;
                else positions++;
                continue;
            }
            char *wp = site_parent(v, &wid), *lp = site_parent(v, &lid);
            const def_t *W = &P->W[i], *L = &P->L[i];
            const char *wcur = 0, *lcur = 0;
            for (int r = 0; r < W->nrows && !wcur; r++)
                if (W->rows[r].atom && *W->rows[r].atom) wcur = W->rows[r].atom;
            for (int r = 0; r < L->nrows && !lcur; r++)
                if (L->rows[r].atom && *L->rows[r].atom) lcur = L->rows[r].atom;
            const char **lfill = malloc(sizeof(char *) * (L->nrows + 1));
            for (int r = 0; r < L->nrows; r++) {
                if (L->rows[r].atom && *L->rows[r].atom) lcur = L->rows[r].atom;
                lfill[r] = lcur;
            }
            for (int r = 0; r < W->nrows; r++) {
                if (W->rows[r].atom && *W->rows[r].atom) wcur = W->rows[r].atom;
                int li = P->assign[i] ? P->assign[i][r] : r;
                const char *la = li < L->nrows ? lfill[li] : 0;
                if (wcur && la && !strcmp(wcur, la)) continue;
                exact++;
                if (parents_equal(wcur, wp, la, lp)) gauge++;
                else real++;
            }
            free(lfill);
        }
    }
    int strict_pass = !positions && !exact && !presentation;
    sb_fmt(b, "{\"pass\": %s", !positions && !real ? "true" : "false");
    if (!strict_pass) sb_cat(b, ", \"strict_pass\": false");
    if (gauge) sb_fmt(b, ", \"label_gauge_mismatches\": %ld", gauge);
    if (presentation) sb_fmt(b, ", \"presentation_assignment_modes\": %ld", presentation);
    if (gauge || presentation) {
        sb_cat(b, ", \"gauge\": [");
        if (gauge) sb_cat(b, "\"orbit_ordinal\"");
        if (presentation) sb_cat(b, gauge ? ", \"mode_family_assignment\"" : "\"mode_family_assignment\"");
        sb_cat(b, "]");
    }
    if (real) sb_fmt(b, ", \"label_real_mismatches\": %ld", real);
    if (positions) sb_fmt(b, ", \"position_failures\": %ld", positions);
    sb_cat(b, "}");
    return (!positions && !real) | (strict_pass << 1);
}

static void mode_definition(vm_t *v, sb *b)
{
    sb err = {0}, L = {0};
    ilist wid = {0}, lid = {0};
    const char *divergence = 0;
    int physical = 0, strict = 0, nlayers = 0;
    static const char *const names[8] = {"kind_set", "irrep_set", "site_set", "site_irrep", "mode_count", "normfactor", "row_atom", "vector"};
    if (identities(v, 0, &wid, &err) || identities(v, 1, &lid, &err)) {
        sb_cat(b, "{\"error\": ");
        sb_json(b, err.s);
        sb_cat(b, "}");
        free(err.s);
        free(wid.id);
        free(lid.id);
        return;
    }
    slist w = {0}, l = {0}, mi = {0}, ex = {0};
    for (int k = 0; k < NKIND; k++) {
        if (web_count(v, k) > 0) sl_add(&w, (char *)KIND[k]);
        if (local_count(v, k) > 0) sl_add(&l, (char *)KIND[k]);
    }
    int pass = diff_layer(&w, &l, &mi, &ex);
    sb_diff(&L, pass, &mi, &ex);
    sb_cat(&L, "}");
    nlayers = 1;
    sb layer_text[8] = {{0}};
    layer_text[0] = L;
    int irreps_strict = 1, norms_pass = 1, rows_strict = 1, vector_gauge = 0;
    if (!pass) divergence = "mode.mode_definition.kind_set";
    for (int level = 0; level < 3 && !divergence; level++) {
        slist a = {0}, c = {0}, m2 = {0}, e2 = {0};
        for (int i = 0; i < wid.n; i++) sl_add(&a, token(v, &wid.id[i], level, 1));
        for (int i = 0; i < lid.n; i++) sl_add(&c, token(v, &lid.id[i], level, 1));
        int p = diff_layer(&a, &c, &m2, &e2);
        sb t = {0};
        sb_diff(&t, p, &m2, &e2);
        if (level == 0) {
            slist pa = {0}, pc = {0}, pm = {0}, pe = {0};
            for (int i = 0; i < wid.n; i++) sl_add(&pa, token(v, &wid.id[i], 0, 0));
            for (int i = 0; i < lid.n; i++) sl_add(&pc, token(v, &lid.id[i], 0, 0));
            if (!diff_layer(&pa, &pc, &pm, &pe)) {
                irreps_strict = 0;
                sb_cat(&t, ", \"strict_pass\": false, \"gauge\": [\"reciprocal_coordinate\"]");
                if (pm.n) {
                    sb_cat(&t, ", \"presentation_missing\": ");
                    sb_list(&t, &pm);
                }
                if (pe.n) {
                    sb_cat(&t, ", \"presentation_extra\": ");
                    sb_list(&t, &pe);
                }
            }
            free(pa.s);
            free(pc.s);
            free(pm.s);
            free(pe.s);
        }
        sb_cat(&t, "}");
        layer_text[nlayers++] = t;
        free(a.s);
        free(c.s);
        free(m2.s);
        free(e2.s);
        if (!p) divergence = level == 0 ? "mode.mode_definition.irrep_set" : level == 1 ? "mode.mode_definition.site_set" : "mode.mode_definition.site_irrep";
    }
    if (!divergence) {
        slist a = {0}, c = {0};
        for (int i = 0; i < wid.n; i++) sl_add(&a, token(v, &wid.id[i], 2, 1));
        for (int i = 0; i < lid.n; i++) sl_add(&c, token(v, &lid.id[i], 2, 1));
        qsort(a.s, a.n, sizeof(char *), cmp_str);
        qsort(c.s, c.n, sizeof(char *), cmp_str);
        sb t = {0}, mm = {0};
        int nmis = 0;
        for (int i = 0, j = 0; i < a.n || j < c.n;) {
            const char *key = i < a.n && (j >= c.n || strcmp(a.s[i], c.s[j]) <= 0) ? a.s[i] : c.s[j];
            int x = 0, y = 0;
            while (i < a.n && !strcmp(a.s[i], key)) i++, x++;
            while (j < c.n && !strcmp(c.s[j], key)) j++, y++;
            if (x == y) continue;
            sb_cat(&mm, nmis++ ? ", " : "");
            sb_json(&mm, key);
            sb_fmt(&mm, ": {\"web\": %d, \"local\": %d}", x, y);
        }
        int section = 0;
        sb sm = {0};
        for (int k = 0; k < NKIND; k++)
            if (!kind_count_match(v, k)) {
                sb_cat(&sm, section++ ? ", " : "");
                sb_fmt(&sm, "\"%s\": {\"web\": %d, \"local\": %d}", KIND[k], web_count(v, k), local_count(v, k));
            }
        int p = !nmis && !section;
        sb_fmt(&t, "{\"pass\": %s", p ? "true" : "false");
        if (nmis) {
            sb_cat(&t, ", \"mismatches\": {");
            sb_cat(&t, mm.s);
            sb_cat(&t, "}");
        }
        if (section) {
            sb_cat(&t, ", \"section_mismatches\": {");
            sb_cat(&t, sm.s);
            sb_cat(&t, "}");
        }
        sb_cat(&t, "}");
        layer_text[nlayers++] = t;
        free(mm.s);
        free(sm.s);
        free(a.s);
        free(c.s);
        if (!p) divergence = "mode.mode_definition.mode_count";
    }
    if (!divergence) {
        sb t = {0}, mm = {0};
        int nmis = 0;
        for (int k = 0; k < NKIND; k++) {
            int nz = web_count(v, k) < local_count(v, k) ? web_count(v, k) : local_count(v, k);
            for (int i = 0; i < nz; i++) {
                int eq;
                if (k == STRAIN) eq = optional_close(v->wstrain[i].has_nf, v->wstrain[i].nf, v->lstrain[i].has_nf, v->lstrain[i].nf, v->tol);
                else eq = optional_close(v->P[k].W[i].has_nf, v->P[k].W[i].nf, v->P[k].L[i].has_nf, v->P[k].L[i].nf, v->tol);
                if (eq) continue;
                sb_cat(&mm, nmis++ ? ", " : "");
                sb_fmt(&mm, "\"%s:%d\"", KIND[k], i + 1);
            }
        }
        norms_pass = !nmis;
        sb_fmt(&t, "{\"pass\": %s", norms_pass ? "true" : "false");
        if (nmis) {
            sb_cat(&t, ", \"mismatches\": [");
            sb_cat(&t, mm.s);
            sb_cat(&t, "]");
        }
        sb_cat(&t, "}");
        layer_text[nlayers++] = t;
        free(mm.s);
        sb r = {0};
        int rr = row_atom_layer(v, &r, &err);
        if (rr < 0) {
            sb_cat(b, "{\"error\": ");
            sb_json(b, err.s);
            sb_cat(b, "}");
            for (int i = 0; i < nlayers; i++) free(layer_text[i].s);
            free(r.s);
            free(err.s);
            free(wid.id);
            free(lid.id);
            free(w.s);
            free(l.s);
            free(mi.s);
            free(ex.s);
            return;
        }
        layer_text[nlayers++] = r;
        rows_strict = rr >> 1 & 1;
        if (!(rr & 1)) divergence = "mode.mode_definition.row_atom";
    }
    if (!divergence) {
        sb t = {0}, fl = {0};
        int gauge[4] = {0}, nfail = 0;
        for (int k = 0; k < NKIND; k++) {
            int s = kind_status(v, k);
            if (s == S_SIGN) gauge[0] = 1;
            else if (s == S_SIGNED_PERMUTATION) gauge[3] = 1;
            else if (s == S_FAMILY_BASIS) gauge[2] = 1;
            else if (s != S_OK && !(web_count(v, k) == 0 && local_count(v, k) == 0)) {
                sb_cat(&fl, nfail++ ? ", " : "");
                sb_fmt(&fl, "\"%s\": \"%s\"", KIND[k], STATUS[s]);
            }
            if (k != STRAIN && v->P[k].validated) gauge[1] = 1;
        }
        static const char *const gnames[4] = {"basis_frame", "family_basis", "sign", "signed_permutation"};
        int g2[4] = {gauge[1], gauge[2], gauge[0], gauge[3]}, ng = 0;
        sb_fmt(&t, "{\"pass\": %s", nfail ? "false" : "true");
        for (int i = 0; i < 4; i++)
            if (g2[i]) {
                sb_cat(&t, ng++ ? ", " : ", \"gauge\": [");
                sb_json(&t, gnames[i]);
            }
        if (ng) sb_cat(&t, "]");
        if (nfail) {
            sb_cat(&t, ", \"failures\": {");
            sb_cat(&t, fl.s);
            sb_cat(&t, "}");
        }
        sb_cat(&t, "}");
        layer_text[nlayers++] = t;
        free(fl.s);
        vector_gauge = ng > 0;
        physical = !nfail;
        strict = physical && irreps_strict && !vector_gauge && rows_strict && norms_pass;
        if (nfail) divergence = "mode.mode_definition.vector";
        else if (!irreps_strict) divergence = "mode.mode_definition.irrep_set.gauge";
        else if (!norms_pass) divergence = "mode.mode_definition.normfactor";
        else if (!rows_strict) divergence = "mode.mode_definition.row_atom.gauge";
        else if (vector_gauge) divergence = "mode.mode_definition.vector.gauge";
    }
    sb_cat(b, "{\"layers\": {");
    for (int i = 0; i < 8; i++) {
        if (i) sb_cat(b, ", ");
        sb_json(b, names[i]);
        sb_cat(b, ": ");
        sb_cat(b, i < nlayers ? layer_text[i].s : "null");
    }
    sb_cat(b, "}, \"divergence\": ");
    sb_json(b, divergence);
    sb_fmt(b, ", \"physical\": %s, \"strict\": %s}", physical ? "true" : "false", strict ? "true" : "false");
    for (int i = 0; i < nlayers; i++) free(layer_text[i].s);
    free(err.s);
    free(wid.id);
    free(lid.id);
    free(w.s);
    free(l.s);
    free(mi.s);
    free(ex.s);
}

/* ---- API ---- */

void *vd_open(const char *web, int64_t wlen, const char *local, int64_t llen, double tol)
{
    vm_t *v = calloc(1, sizeof *v);
    v->tol = tol;
    v->vtol = visible_bound(tol);
    for (int k = 0; k < 3; k++) parse_definitions(v, web, web + wlen, k);
    parse_strain(v, web, web + wlen);
    if (read_local(v, local, llen)) return v;
    for (int k = 0; k < 3; k++) {
        exact_web_order(v, k);
        vector_payload(v, &v->P[k], v->web[k], v->nweb[k], v->local[k], v->nlocal[k], k, 0);
    }
    strain_payload(v);
    return v;
}

const char *vd_error(void *h) { return ((vm_t *)h)->failed ? ((vm_t *)h)->err : 0; }

int vd_web_count(void *h, int kind)
{
    vm_t *v = h;
    return kind == STRAIN ? v->nwstrain : v->nweb[kind];
}

const char *vd_status(void *h, int kind)
{
    vm_t *v = h;
    return STATUS[kind_status(v, kind)];
}

const char *vd_strict_status(void *h, int kind)
{
    vm_t *v = h;
    return STATUS[kind == STRAIN ? v->SP.strict : v->P[kind].strict];
}

int vd_accepted(void *h, int kind) { return accepted(kind_status(h, kind)); }

int vd_transport(void *h, int kind, const int64_t *basis, const double *shift)
{
    int U[9];
    for (int i = 0; i < 9; i++) U[i] = (int)basis[i];
    return transport(h, kind, U, shift);
}

const char *vd_result(void *h)
{
    vm_t *v = h;
    if (!v->result) {
        sb b = {0};
        mode_definition(v, &b);
        v->result = b.s;
    }
    return v->result;
}

void vd_close(void *h)
{
    vm_t *v = h;
    for (int i = 0; i < v->nowned; i++) free(v->owned[i]);
    free(v->owned);
    for (int i = 0; i < v->pos_cap; i++) free(v->pos[i].pos);
    free(v->pos);
    for (int i = 0; i < v->cand_cap; i++) cand_free(v->cand[i].val);
    free(v->cand);
    free(v->result);
    free(v);
}
