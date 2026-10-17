/* SPEC.md at the Gamma point for displacements, from the declared symmetry
 * operations and the sites alone (no Source tables): the group P = G/T with
 * its closure checked, the characters of P (Dixon-Schneider modulo a prime,
 * values recovered as sums of roots of unity and checked by the
 * orthogonality relations over Q(zeta_m)), the physically irreducible
 * characters, their multiplicities in each site's displacements (Frobenius
 * reciprocity), one irreducible copy W of each, and every G-orbit of nonzero
 * isotropy subspaces of W with its stabilizer. */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cyclo.h"
#include "json.h"
#include "md.h"
#include "modp.h"
#include "pyconv.h"
#include "rational.h"

#define SP_TOL 1e-3 /* stage 0: matching declared coordinates */
#define SP_MAX_DIM 8

typedef struct {
    int rot[3][3];
    int tran[3]; /* 1/24 units, 0..23 */
} sp_op;

typedef struct {
    int n, natoms, nclasses, ncent;
    sp_op *ops;      /* the elements of P */
    int *mult, *inv; /* n x n, n */
    int *perm;       /* n x nclasses: atom class images */
    int *class_site; /* site index of each atom class */
    int *shift;      /* n x nclasses x 3: g x_a - x_{g a} in 1/24 units (a lattice vector) */
    int cent[8][3];
    int ncls, *cls, *cls_size, *cls_rep, *cls_inv, *order;
} sp_group;

static int fail_status;

static int spfail(const char *what)
{
    if (!fail_status) md_fail(MDE_VALUE, "%s", what);
    fail_status = 1;
    return 1;
}

static int mod24(int x) { return ((x % 24) + 24) % 24; }

static void compose(const sp_op *a, const sp_op *b, sp_op *out)
{
    sp_op r;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            r.rot[i][j] = 0;
            for (int k = 0; k < 3; k++) r.rot[i][j] += a->rot[i][k] * b->rot[k][j];
        }
        int t = a->tran[i];
        for (int k = 0; k < 3; k++) t += a->rot[i][k] * b->tran[k];
        r.tran[i] = mod24(t);
    }
    *out = r;
}

/* the translation coset key: the least of t + c over the centering vectors */
static void canonical(sp_op *op, const int (*cent)[3], int ncent)
{
    int best[3] = {24, 24, 24};
    for (int c = 0; c < ncent; c++) {
        int t[3];
        for (int i = 0; i < 3; i++) t[i] = mod24(op->tran[i] + cent[c][i]);
        if (memcmp(t, best, sizeof t) < 0) memcpy(best, t, sizeof best);
    }
    memcpy(op->tran, best, sizeof best);
}

static int same_op(const sp_op *a, const sp_op *b) { return !memcmp(a, b, sizeof *a); }

static void apply_op(const sp_op *op, const double *x, double *out)
{
    for (int i = 0; i < 3; i++) {
        double v = op->tran[i] / 24.0;
        for (int k = 0; k < 3; k++) v += op->rot[i][k] * x[k];
        out[i] = v - floor(v);
    }
}

static double frac_dist(const double *a, const double *b)
{
    double worst = 0;
    for (int i = 0; i < 3; i++) {
        double d = fabs(fmod(a[i] - b[i] + 10.5, 1.0) - 0.5);
        if (d > worst) worst = d;
    }
    return worst;
}

static int find_atom(const double (*atoms)[3], int natoms, const double *y)
{
    for (int i = 0; i < natoms; i++)
        if (frac_dist(atoms[i], y) < SP_TOL) return i;
    return -1;
}

static int read_text_double(const jvalue *v, double *out)
{
    if (!jis(v, J_STRING)) return py_float_of(v, out);
    char buf[64];
    int n = v->count < 63 ? (int)v->count : 63;
    memcpy(buf, v->u.s, (size_t)n);
    buf[n] = 0;
    char *paren = strchr(buf, '(');
    if (paren) *paren = 0;
    char *slash = strchr(buf, '/');
    if (slash) {
        *slash = 0;
        *out = strtod(buf, 0) / strtod(slash + 1, 0);
        return 0;
    }
    char *end;
    *out = strtod(buf, &end);
    return end == buf;
}

/* identity, inverses, conjugacy classes and element orders from the multiplication table */
static int finish_group(sp_group *g)
{
    int n = g->n, identity = -1;
    for (int e = 0; e < n && identity < 0; e++) {
        int ok = 1;
        for (int b = 0; b < n; b++) ok &= g->mult[e * n + b] == b;
        if (ok) identity = e;
    }
    if (identity < 0) return spfail("no identity");
    if (!g->inv) g->inv = malloc(sizeof(int) * (size_t)n);
    if (!g->inv) return spfail("out of memory");
    for (int a = 0; a < n; a++) {
        g->inv[a] = -1;
        for (int b = 0; b < n; b++)
            if (g->mult[a * n + b] == identity) g->inv[a] = b;
        if (g->inv[a] < 0) return spfail("not a group");
    }
    /* conjugacy classes and element orders */
    g->cls = malloc(sizeof(int) * (size_t)n);
    g->cls_size = calloc((size_t)n, sizeof(int));
    g->cls_rep = malloc(sizeof(int) * (size_t)n);
    g->cls_inv = malloc(sizeof(int) * (size_t)n);
    g->order = malloc(sizeof(int) * (size_t)n);
    if (!g->cls || !g->cls_size || !g->cls_rep || !g->cls_inv || !g->order) return spfail("out of memory");
    for (int a = 0; a < n; a++) g->cls[a] = -1;
    g->ncls = 0;
    for (int a = 0; a < n; a++) {
        if (g->cls[a] >= 0) continue;
        for (int h = 0; h < n; h++) {
            int c = g->mult[g->mult[h * n + a] * n + g->inv[h]];
            if (g->cls[c] < 0) {
                g->cls[c] = g->ncls;
                g->cls_size[g->ncls]++;
            }
        }
        g->cls_rep[g->ncls++] = a;
    }
    for (int k = 0; k < g->ncls; k++) g->cls_inv[k] = g->cls[g->inv[g->cls_rep[k]]];
    for (int a = 0; a < n; a++) {
        int x = a, o = 1;
        while (x != identity) {
            x = g->mult[x * n + a];
            o++;
        }
        g->order[a] = o;
    }
    return 0;
}

/* ---- stage 0: the group, its closure, the atoms ---- */

static int build_group(const jvalue *cif_info, const int *sites, int nsites, sp_group *g)
{
    const jvalue *ops = jget(cif_info, "symmetry_operations"), *atom_sites = jget(cif_info, "atom_sites");
    if (!jis(ops, J_ARRAY) || !ops->count) return spfail("no symmetry operations");
    int nops = (int)ops->count;
    sp_op *all = calloc((size_t)nops, sizeof *all);
    int (*cent)[3] = calloc((size_t)nops, sizeof *cent), ncent = 0;
    if (!all || !cent) return spfail("out of memory");
    for (int o = 0; o < nops; o++) {
        const jvalue *t = jat(ops, o);
        if (!jis(t, J_STRING)) return spfail("symmetry operation is not a string");
        int rot[3][3], tran[3];
        char notation;
        if (gemmi_parse_triplet(t->u.s, t->count, rot, tran, &notation)) return spfail("unreadable symmetry operation");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (rot[i][j] % 24) return spfail("symmetry operation with a fractional rotation");
                all[o].rot[i][j] = rot[i][j] / 24;
            }
            all[o].tran[i] = mod24(tran[i]);
        }
        int ident = 1;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) ident &= all[o].rot[i][j] == (i == j);
        if (ident) memcpy(cent[ncent++], all[o].tran, sizeof cent[0]);
    }
    if (!ncent) return spfail("symmetry operations lack the identity");
    if (ncent > 8) return spfail("too many centering translations");
    memcpy(g->cent, cent, sizeof(int[3]) * (size_t)ncent);
    /* P: the operations modulo the centering translations, closed under products */
    g->ops = calloc((size_t)nops, sizeof *g->ops);
    if (!g->ops) return spfail("out of memory");
    g->n = 0;
    for (int o = 0; o < nops; o++) {
        sp_op op = all[o];
        canonical(&op, (const int (*)[3])cent, ncent);
        int seen = 0;
        for (int e = 0; e < g->n && !seen; e++) seen = same_op(&g->ops[e], &op);
        if (!seen) g->ops[g->n++] = op;
    }
    int n = g->n;
    g->ncent = ncent;
    g->mult = malloc(sizeof(int) * (size_t)n * n);
    g->inv = malloc(sizeof(int) * (size_t)n);
    if (!g->mult || !g->inv) return spfail("out of memory");
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++) {
            sp_op r;
            compose(&g->ops[a], &g->ops[b], &r);
            canonical(&r, (const int (*)[3])cent, ncent);
            int found = -1;
            for (int e = 0; e < n && found < 0; e++)
                if (same_op(&g->ops[e], &r)) found = e;
            if (found < 0) return spfail("symmetry operations are not closed under composition");
            g->mult[a * n + b] = found;
        }
    /* atoms of the selected sites, every image, and their classes modulo centering */
    int cap = nsites * nops + 1, natoms = 0;
    double (*atoms)[3] = malloc(sizeof *atoms * (size_t)cap);
    int *atom_site = malloc(sizeof(int) * (size_t)cap);
    if (!atoms || !atom_site) return spfail("out of memory");
    for (int s = 0; s < nsites; s++) {
        const jvalue *site = jat(atom_sites, sites[s] - 1), *fract = jget(site, "fract");
        double x[3];
        if (!jis(fract, J_ARRAY) || fract->count != 3) return spfail("site without coordinates");
        for (int i = 0; i < 3; i++)
            if (read_text_double(jat(fract, i), &x[i])) return spfail("unreadable site coordinate");
        for (int o = 0; o < nops; o++) {
            double y[3];
            apply_op(&all[o], x, y);
            int at = find_atom((const double (*)[3])atoms, natoms, y);
            if (at < 0) {
                memcpy(atoms[natoms], y, sizeof y);
                atom_site[natoms++] = s;
            } else if (atom_site[at] != s)
                return spfail("two selected sites overlap");
        }
    }
    int *class_of = malloc(sizeof(int) * (size_t)natoms);
    g->class_site = malloc(sizeof(int) * (size_t)natoms);
    if (!class_of || !g->class_site) return spfail("out of memory");
    int nclasses = 0;
    for (int a = 0; a < natoms; a++) class_of[a] = -1;
    for (int a = 0; a < natoms; a++) {
        if (class_of[a] >= 0) continue;
        for (int c = 0; c < ncent; c++) {
            double y[3];
            for (int i = 0; i < 3; i++) {
                double v = atoms[a][i] + cent[c][i] / 24.0;
                y[i] = v - floor(v);
            }
            int b = find_atom((const double (*)[3])atoms, natoms, y);
            if (b < 0) return spfail("a centering translation does not map the structure onto itself");
            class_of[b] = nclasses;
        }
        g->class_site[nclasses++] = atom_site[a];
    }
    g->natoms = natoms;
    g->nclasses = nclasses;
    g->perm = malloc(sizeof(int) * (size_t)n * nclasses);
    g->shift = malloc(sizeof(int) * 3 * (size_t)n * nclasses);
    int *rep_atom = malloc(sizeof(int) * (size_t)nclasses);
    if (!g->perm || !g->shift || !rep_atom) return spfail("out of memory");
    for (int a = natoms - 1; a >= 0; a--) rep_atom[class_of[a]] = a;
    for (int e = 0; e < n; e++)
        for (int a = 0; a < natoms; a++) {
            double y[3];
            apply_op(&g->ops[e], atoms[a], y);
            int b = find_atom((const double (*)[3])atoms, natoms, y);
            if (b < 0) return spfail("a symmetry operation does not map the structure onto itself");
            g->perm[e * nclasses + class_of[a]] = class_of[b];
            if (a != rep_atom[class_of[a]]) continue;
            const double *x = atoms[a], *xb = atoms[rep_atom[class_of[b]]];
            for (int i = 0; i < 3; i++) {
                double v = g->ops[e].tran[i] / 24.0 - xb[i];
                for (int k2 = 0; k2 < 3; k2++) v += g->ops[e].rot[i][k2] * x[k2];
                double q = nearbyint(v * 24.0);
                if (fabs(v * 24.0 - q) > 24.0 * SP_TOL) return spfail("an atom image is off the lattice");
                g->shift[(e * nclasses + class_of[a]) * 3 + i] = (int)q;
            }
        }
    free(rep_atom);
    if (finish_group(g)) return 1;
    free(all);
    free(cent);
    free(atoms);
    free(atom_site);
    free(class_of);
    return 0;
}

/* ---- characters modulo p ---- */

static uint64_t PR;

static uint64_t pmul(uint64_t a, uint64_t b) { return mulmod(a, b, PR); }
static uint64_t padd(uint64_t a, uint64_t b) { return (a + b) % PR; }
static uint64_t psub(uint64_t a, uint64_t b) { return (a + PR - b) % PR; }
static uint64_t pinv(uint64_t a) { return powmod(a, PR - 2, PR); }
static uint64_t pint(int64_t v) { return (uint64_t)(((v % (int64_t)PR) + (int64_t)PR) % (int64_t)PR); }

/* reduced row echelon form in place; returns the rank, pivots[] the pivot columns */
static int prref(uint64_t *m, int rows, int cols, int *pivots)
{
    int r = 0;
    for (int c = 0; c < cols && r < rows; c++) {
        int p = -1;
        for (int i = r; i < rows && p < 0; i++)
            if (m[i * cols + c]) p = i;
        if (p < 0) continue;
        for (int j = 0; j < cols; j++) {
            uint64_t t = m[p * cols + j];
            m[p * cols + j] = m[r * cols + j];
            m[r * cols + j] = t;
        }
        uint64_t s = pinv(m[r * cols + c]);
        for (int j = 0; j < cols; j++) m[r * cols + j] = pmul(m[r * cols + j], s);
        for (int i = 0; i < rows; i++) {
            if (i == r || !m[i * cols + c]) continue;
            uint64_t f = m[i * cols + c];
            for (int j = 0; j < cols; j++) m[i * cols + j] = psub(m[i * cols + j], pmul(f, m[r * cols + j]));
        }
        pivots[r++] = c;
    }
    return r;
}

/* the nullspace of rows x cols m as row vectors; returns the count */
static int pnull(const uint64_t *m, int rows, int cols, uint64_t *out)
{
    uint64_t *w = malloc(sizeof(uint64_t) * (size_t)(rows * cols + 1));
    int *piv = malloc(sizeof(int) * (size_t)(cols + 1)), count = 0;
    memcpy(w, m, sizeof(uint64_t) * (size_t)rows * cols);
    int rank = prref(w, rows, cols, piv);
    for (int f = 0; f < cols; f++) {
        int is_pivot = 0;
        for (int i = 0; i < rank; i++) is_pivot |= piv[i] == f;
        if (is_pivot) continue;
        uint64_t *v = out + (size_t)count * cols;
        memset(v, 0, sizeof(uint64_t) * (size_t)cols);
        v[f] = 1;
        for (int i = 0; i < rank; i++) v[piv[i]] = psub(0, w[i * cols + f]);
        count++;
    }
    free(w);
    free(piv);
    return count;
}

typedef struct {
    int deg;
    uint64_t c[512];
} poly;

static void ptrim(poly *a)
{
    while (a->deg > 0 && !a->c[a->deg]) a->deg--;
}

static void pmod_poly(poly *a, const poly *b)
{
    uint64_t lead = pinv(b->c[b->deg]);
    while (a->deg >= b->deg && (a->deg > 0 || a->c[0])) {
        uint64_t f = pmul(a->c[a->deg], lead);
        int shift = a->deg - b->deg;
        for (int i = 0; i <= b->deg; i++) a->c[i + shift] = psub(a->c[i + shift], pmul(f, b->c[i]));
        if (a->deg == 0) break;
        a->deg--;
        ptrim(a);
        if (a->deg < b->deg) break;
    }
    ptrim(a);
}

static void pmulmod(const poly *a, const poly *b, const poly *m, poly *out)
{
    poly r = {a->deg + b->deg, {0}};
    for (int i = 0; i <= a->deg; i++)
        for (int j = 0; j <= b->deg; j++) r.c[i + j] = padd(r.c[i + j], pmul(a->c[i], b->c[j]));
    ptrim(&r);
    pmod_poly(&r, m);
    *out = r;
}

static void pgcd(poly a, poly b, poly *out)
{
    ptrim(&a);
    ptrim(&b);
    while (b.deg > 0 || b.c[0]) {
        poly r = a;
        pmod_poly(&r, &b);
        a = b;
        b = r;
    }
    *out = a;
}

/* every root in F_p of a polynomial that splits into linear factors */
static int proots(poly f, uint64_t *roots, uint64_t seed)
{
    ptrim(&f);
    if (f.deg == 0) return 0;
    if (f.deg == 1) {
        roots[0] = psub(0, pmul(f.c[0], pinv(f.c[1])));
        return 1;
    }
    for (int attempt = 0; attempt < 64; attempt++) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        poly base = {1, {seed % PR, 1}}, acc = {0, {1}};
        uint64_t e = (PR - 1) / 2;
        while (e) {
            if (e & 1) pmulmod(&acc, &base, &f, &acc);
            pmulmod(&base, &base, &f, &base);
            e >>= 1;
        }
        acc.c[0] = psub(acc.c[0], 1);
        ptrim(&acc);
        poly d;
        pgcd(f, acc, &d);
        if (d.deg > 0 && d.deg < f.deg) {
            poly rest = {0, {0}}, num = f;
            int dq = f.deg - d.deg;
            rest.deg = dq;
            uint64_t lead = pinv(d.c[d.deg]);
            for (int k = dq; k >= 0; k--) {
                uint64_t c = pmul(num.c[k + d.deg], lead);
                rest.c[k] = c;
                for (int i = 0; i <= d.deg; i++) num.c[k + i] = psub(num.c[k + i], pmul(c, d.c[i]));
            }
            int n1 = proots(d, roots, seed + 1);
            return n1 + proots(rest, roots + n1, seed + 2);
        }
    }
    return -1;
}

/* the characteristic polynomial of an s x s matrix (Faddeev-LeVerrier) */
static void charpoly(const uint64_t *b, int s, poly *out)
{
    uint64_t *m = calloc((size_t)s * s + 1, sizeof(uint64_t)), *t = calloc((size_t)s * s + 1, sizeof(uint64_t));
    out->deg = s;
    out->c[s] = 1;
    for (int k = 1; k <= s; k++) {
        /* m = b (m_prev + c_{s-k+1} I) */
        for (int i = 0; i < s; i++) m[i * s + i] = padd(m[i * s + i], out->c[s - k + 1]);
        for (int i = 0; i < s; i++)
            for (int j = 0; j < s; j++) {
                uint64_t v = 0;
                for (int l = 0; l < s; l++) v = padd(v, pmul(b[i * s + l], m[l * s + j]));
                t[i * s + j] = v;
            }
        memcpy(m, t, sizeof(uint64_t) * (size_t)s * s);
        uint64_t tr = 0;
        for (int i = 0; i < s; i++) tr = padd(tr, m[i * s + i]);
        out->c[s - k] = psub(0, pmul(tr, pinv(pint(k))));
    }
    free(m);
    free(t);
}

typedef struct {
    int count;
    cy *values; /* count x ncls: the exact characters by class */
    int *conj;  /* index of the complex conjugate character */
} sp_chars;

static int characters(const sp_group *g, const cy_field *f, uint64_t w, sp_chars *out)
{
    int n = g->n, k = g->ncls, m = f->m;
    if (k > 255) return spfail("too many conjugacy classes");
    /* a[(i * k + j) * k + l]: pairs (x, y) in K_i x K_j with x y = rep of K_l */
    int *a = calloc((size_t)k * k * k, sizeof(int));
    uint64_t *space = malloc(sizeof(uint64_t) * (size_t)k * k * 2), *b = malloc(sizeof(uint64_t) * (size_t)k * k);
    uint64_t *ev = malloc(sizeof(uint64_t) * (size_t)k * k * 2), *null = malloc(sizeof(uint64_t) * (size_t)k * k * 2);
    int *dims = malloc(sizeof(int) * (size_t)(k + 1)), *piv = malloc(sizeof(int) * (size_t)(k + 1));
    if (!a || !space || !b || !ev || !null || !dims || !piv) return spfail("out of memory");
    for (int l = 0; l < k; l++) {
        int z = g->cls_rep[l];
        for (int x = 0; x < n; x++) {
            int y = g->mult[g->inv[x] * n + z];
            a[(g->cls[x] * k + g->cls[y]) * k + l]++;
        }
    }
    /* subspaces as consecutive blocks of row vectors in space; start with the whole of F_p^k */
    int nsub = 1, total = k;
    memset(space, 0, sizeof(uint64_t) * (size_t)k * k);
    for (int i = 0; i < k; i++) space[i * k + i] = 1;
    dims[0] = k;
    for (int j = 1; j < k; j++) {
        uint64_t *next = ev;
        int nnext = 0, rows_next = 0, *nd = malloc(sizeof(int) * (size_t)(k + 1));
        int offset = 0;
        for (int s = 0; s < nsub; s++) {
            int dim = dims[s];
            uint64_t *basis = space + (size_t)offset * k;
            offset += dim;
            if (dim == 1) {
                memcpy(next + (size_t)rows_next * k, basis, sizeof(uint64_t) * (size_t)k);
                rows_next++;
                nd[nnext++] = 1;
                continue;
            }
            /* A_j on the subspace: image of each basis vector in basis coordinates */
            uint64_t *ech = malloc(sizeof(uint64_t) * (size_t)dim * k);
            memcpy(ech, basis, sizeof(uint64_t) * (size_t)dim * k);
            int rank = prref(ech, dim, k, piv);
            if (rank != dim) return spfail("character subspace lost rank");
            /* re-express the basis as the echelon rows so coordinates read off at the pivots */
            for (int r = 0; r < dim; r++) {
                for (int i = 0; i < k; i++) {
                    uint64_t v = 0;
                    for (int l = 0; l < k; l++)
                        if (a[(i * k + j) * k + l]) v = padd(v, pmul(pint(a[(i * k + j) * k + l]), ech[r * k + l]));
                    null[i] = v;
                }
                for (int c = 0; c < dim; c++) b[c * dim + r] = null[piv[c]];
            }
            poly cp;
            charpoly(b, dim, &cp);
            poly sq, deriv = {cp.deg - 1, {0}};
            for (int i = 1; i <= cp.deg; i++) deriv.c[i - 1] = pmul(cp.c[i], pint(i));
            pgcd(cp, deriv, &sq);
            /* the squarefree part */
            poly num = cp, part = {cp.deg - sq.deg, {0}};
            uint64_t lead = pinv(sq.c[sq.deg]);
            for (int q = part.deg; q >= 0; q--) {
                uint64_t c = pmul(num.c[q + sq.deg], lead);
                part.c[q] = c;
                for (int i = 0; i <= sq.deg; i++) num.c[q + i] = psub(num.c[q + i], pmul(c, sq.c[i]));
            }
            uint64_t roots[256];
            int nr = proots(part, roots, 0x9e3779b97f4a7c15ULL + (uint64_t)j);
            if (nr != part.deg) return spfail("class matrix eigenvalues outside the field");
            for (int r = 0; r < nr; r++) {
                uint64_t *shifted = malloc(sizeof(uint64_t) * (size_t)dim * dim);
                for (int x = 0; x < dim * dim; x++) shifted[x] = b[x];
                for (int x = 0; x < dim; x++) shifted[x * dim + x] = psub(shifted[x * dim + x], roots[r]);
                uint64_t *coords = malloc(sizeof(uint64_t) * (size_t)dim * dim);
                int cnt = pnull(shifted, dim, dim, coords);
                for (int v = 0; v < cnt; v++) {
                    uint64_t *row = next + (size_t)rows_next * k;
                    memset(row, 0, sizeof(uint64_t) * (size_t)k);
                    for (int c = 0; c < dim; c++)
                        for (int l = 0; l < k; l++) row[l] = padd(row[l], pmul(coords[v * dim + c], ech[c * k + l]));
                    rows_next++;
                }
                nd[nnext++] = cnt;
                free(shifted);
                free(coords);
            }
            free(ech);
        }
        memcpy(space, next, sizeof(uint64_t) * (size_t)rows_next * k);
        memcpy(dims, nd, sizeof(int) * (size_t)nnext);
        nsub = nnext;
        total = rows_next;
        free(nd);
    }
    if (nsub != k || total != k) return spfail("characters not separated by the class matrices");
    int ident = 0;
    for (int e = 0; e < n; e++)
        if (g->order[e] == 1) ident = g->cls[e];
    out->count = k;
    out->values = calloc((size_t)k * k, sizeof(cy));
    out->conj = malloc(sizeof(int) * (size_t)k);
    if (!out->values || !out->conj) return spfail("out of memory");
    uint64_t *chi = malloc(sizeof(uint64_t) * (size_t)k * k);
    for (int c = 0; c < k; c++) {
        uint64_t *omega = space + (size_t)c * k;
        if (!omega[ident]) return spfail("character vector vanishes at the identity");
        uint64_t s = pinv(omega[ident]);
        for (int l = 0; l < k; l++) omega[l] = pmul(omega[l], s);
        uint64_t sum = 0;
        for (int l = 0; l < k; l++) sum = padd(sum, pmul(pmul(omega[l], omega[g->cls_inv[l]]), pinv(pint(g->cls_size[l]))));
        uint64_t d2 = pmul(pint(n), pinv(sum));
        int d = 0;
        for (int t = 1; t * t <= n && !d; t++)
            if (pint(t * t) == d2) d = t;
        if (!d) return spfail("character degree not found");
        for (int l = 0; l < k; l++) chi[c * k + l] = pmul(pmul(pint(d), omega[l]), pinv(pint(g->cls_size[l])));
    }
    /* exact values: chi(g) = sum_j mu_j zeta_o^j, mu_j = (1/o) sum_t chi(g^t) zeta_o^(-jt) */
    for (int c = 0; c < k; c++)
        for (int l = 0; l < k; l++) {
            int rep = g->cls_rep[l], o = g->order[rep];
            uint64_t zo = powmod(w, (uint64_t)(m / o), PR), zinv = pinv(zo);
            cy value, term;
            cy_zero(f, &value);
            for (int jj = 0; jj < o; jj++) {
                uint64_t mu = 0, step = powmod(zinv, (uint64_t)jj, PR), z = 1;
                int power = -1;
                for (int t = 0; t < o; t++) {
                    /* g^t */
                    if (t == 0) {
                        for (int e = 0; e < n; e++)
                            if (g->order[e] == 1) power = e;
                    } else
                        power = g->mult[power * n + rep];
                    mu = padd(mu, pmul(chi[c * k + g->cls[power]], z));
                    z = pmul(z, step);
                }
                mu = pmul(mu, pinv(pint(o)));
                if (mu > (uint64_t)n) return spfail("character value is not a sum of roots of unity");
                if (!mu) continue;
                cy root;
                cy_root(f, (int64_t)jj * (m / o), &root);
                cy_rat(f, (i128)mu, 1, &term);
                cy_mul(f, &term, &root, &term);
                cy_add(f, &value, &term, &value);
            }
            out->values[c * k + l] = value;
        }
    /* exact orthogonality: sum_l |K_l| chi(l) chi'(l*) = |G| delta */
    for (int c = 0; c < k; c++)
        for (int d = 0; d <= c; d++) {
            cy acc, term, size;
            cy_zero(f, &acc);
            for (int l = 0; l < k; l++) {
                cy_mul(f, &out->values[c * k + l], &out->values[d * k + g->cls_inv[l]], &term);
                cy_rat(f, g->cls_size[l], 1, &size);
                cy_mul(f, &term, &size, &term);
                cy_add(f, &acc, &term, &acc);
            }
            cy want;
            cy_rat(f, c == d ? n : 0, 1, &want);
            cy_sub(f, &acc, &want, &acc);
            if (!cy_is_zero(f, &acc) || cy_overflow) return spfail("characters fail the orthogonality relations");
        }
    for (int c = 0; c < k; c++) {
        out->conj[c] = -1;
        for (int d = 0; d < k && out->conj[c] < 0; d++) {
            int same = 1;
            for (int l = 0; l < k && same; l++) {
                cy diff;
                cy_sub(f, &out->values[d * k + l], &out->values[c * k + g->cls_inv[l]], &diff);
                same = cy_is_zero(f, &diff);
            }
            if (same) out->conj[c] = d;
        }
        if (out->conj[c] < 0) return spfail("character without a conjugate");
    }
    free(a);
    free(space);
    free(b);
    free(ev);
    free(null);
    free(dims);
    free(piv);
    free(chi);
    return 0;
}

/* ---- rationals ---- */

static int rat_of_cy(const cy_field *f, const cy *x, rat *out)
{
    for (int i = 1; i < f->phi; i++)
        if (x->c[i]) return 1;
    *out = rat_make(x->c[0], x->den);
    return 0;
}

static void put_int_list(jbuf *b, const int *v, int n)
{
    jb_lit(b, "[");
    for (int i = 0; i < n; i++) {
        char t[24];
        snprintf(t, sizeof t, i ? ", %d" : "%d", v[i]);
        jb_lit(b, t);
    }
    jb_lit(b, "]");
}

/* ---- k != 0: the little group modulo the kernel of the arm character ---- */

/* k . t for t in 1/24 units */
static rat kdot(const rat *k, const int *t24)
{
    rat s = rat_int(0);
    for (int i = 0; i < 3; i++)
        if (t24[i]) s = rat_add(s, rat_mul(k[i], rat_make(t24[i], 24)));
    return s;
}

/* a == b as characters of the lattice: (a - b) . t integral on Z^3 and the centering vectors */
static int k_same(const sp_group *g, const rat *a, const rat *b)
{
    rat d[3];
    for (int i = 0; i < 3; i++) {
        d[i] = rat_sub(a[i], b[i]);
        if (d[i].d != 1) return 0;
    }
    for (int c = 0; c < g->ncent; c++)
        if (kdot(d, g->cent[c]).d != 1) return 0;
    return 1;
}

/* (k R)_j = sum_i k_i R_ij: the arm moved by the inverse of R */
static void k_times(const rat *k, const int (*rot)[3], rat *out)
{
    for (int j = 0; j < 3; j++) {
        out[j] = rat_int(0);
        for (int i = 0; i < 3; i++)
            if (rot[i][j]) out[j] = rat_add(out[j], rat_mul(k[i], rat_int(rot[i][j])));
    }
}

typedef struct {
    sp_group h;      /* H = G_k / ker chi: elements (little index, j) at little * N + j */
    int nlittle, N, star;
    int *little;     /* P indices of the little co-group elements */
    int *little_of;  /* P index -> little index or -1 */
} sp_little;

/* c(g1, g2) = N k . (tau1 + R1 tau2 - tau12) */
static int cocycle(const sp_group *g, const rat *k, int N, int a, int b, int *out)
{
    int ab = g->mult[a * g->n + b], t[3];
    for (int i = 0; i < 3; i++) {
        t[i] = g->ops[a].tran[i] - g->ops[ab].tran[i];
        for (int q = 0; q < 3; q++) t[i] += g->ops[a].rot[i][q] * g->ops[b].tran[q];
    }
    rat c = rat_mul(rat_int(N), kdot(k, t));
    if (c.d != 1) return spfail("translation part outside the lattice");
    *out = (int)(((c.n % N) + N) % N);
    return 0;
}

static int build_little(const sp_group *g, const rat *k, sp_little *L)
{
    int n = g->n;
    memset(L, 0, sizeof *L);
    L->little = malloc(sizeof(int) * (size_t)n);
    L->little_of = malloc(sizeof(int) * (size_t)n);
    if (!L->little || !L->little_of) return spfail("out of memory");
    for (int e = 0; e < n; e++) {
        rat kr[3];
        k_times(k, (const int (*)[3])g->ops[e].rot, kr);
        L->little_of[e] = -1;
        if (k_same(g, kr, k)) {
            L->little_of[e] = L->nlittle;
            L->little[L->nlittle++] = e;
        }
    }
    L->star = n / L->nlittle;
    /* N: the order of chi on the lattice */
    int64_t N = 1;
    int axes[3][3] = {{24, 0, 0}, {0, 24, 0}, {0, 0, 24}};
    for (int v = 0; v < 3 + g->ncent; v++) {
        rat d = kdot(k, v < 3 ? axes[v] : g->cent[v - 3]);
        int64_t den = (int64_t)d.d, a = N, b = den;
        while (b) {
            int64_t t = a % b;
            a = b;
            b = t;
        }
        N = N / a * den;
    }
    if (N > 1024 || (int64_t)L->nlittle * N > 4096) return spfail("k denominator above the reference limit");
    L->N = (int)N;
    int nh = L->nlittle * L->N;
    L->h.n = nh;
    L->h.mult = malloc(sizeof(int) * (size_t)nh * nh);
    if (!L->h.mult) return spfail("out of memory");
    for (int a = 0; a < L->nlittle; a++)
        for (int b = 0; b < L->nlittle; b++) {
            int c;
            if (cocycle(g, k, L->N, L->little[a], L->little[b], &c)) return 1;
            int ab = L->little_of[g->mult[L->little[a] * n + L->little[b]]];
            for (int ja = 0; ja < L->N; ja++)
                for (int jb = 0; jb < L->N; jb++)
                    L->h.mult[(a * L->N + ja) * nh + b * L->N + jb] = ab * L->N + (ja + jb + c) % L->N;
        }
    /* associativity: H must be a group */
    for (int x = 0; x < nh; x++)
        for (int y = 0; y < nh; y++)
            for (int z = 0; z < nh; z += nh > 64 ? 7 : 1)
                if (L->h.mult[L->h.mult[x * nh + y] * nh + z] != L->h.mult[x * nh + L->h.mult[y * nh + z]])
                    return spfail("the little group law is not associative");
    return finish_group(&L->h);
}

/* the H element y x y^-1 for y in P taking k to -k and x = (e, j) in H, by affine composition */
static int conjugate_in_h(const sp_group *g, const sp_little *L, const rat *k, int y, int x, int *out)
{
    int N = L->N, e = L->little[x / N], j = x % N;
    int t1[3], have = 0;
    {
        for (int c = 0; c < g->ncent && !have; c++)
            for (int a = 0; a < N && !have; a++)
                for (int b2 = 0; b2 < N && !have; b2++)
                    for (int d = 0; d < N && !have; d++) {
                        int t[3] = {24 * a + g->cent[c][0], 24 * b2 + g->cent[c][1], 24 * d + g->cent[c][2]};
                        rat v = rat_mul(rat_int(N), kdot(k, t));
                        if (v.d == 1 && ((v.n % N) + N) % N == 1 % N) {
                            memcpy(t1, t, sizeof t1);
                            have = 1;
                        }
                    }
        if (!have) return spfail("no lattice vector generates the arm character");
    }
    /* x as an affine map: (R_e, tau_e + j t1) */
    int rx[3][3], tx[3], ry[3][3], ty[3], ri[3][3], ti[3];
    memcpy(rx, g->ops[e].rot, sizeof rx);
    for (int i = 0; i < 3; i++) tx[i] = g->ops[e].tran[i] + j * t1[i];
    memcpy(ry, g->ops[y].rot, sizeof ry);
    memcpy(ty, g->ops[y].tran, sizeof ty);
    int yi = g->inv[y];
    memcpy(ri, g->ops[yi].rot, sizeof ri);
    for (int i = 0; i < 3; i++) {
        ti[i] = 0;
        for (int q = 0; q < 3; q++) ti[i] -= ri[i][q] * ty[q];
    }
    /* z = y x y^-1 */
    int r1[3][3], t1x[3], rz[3][3], tz[3];
    for (int i = 0; i < 3; i++) {
        t1x[i] = ty[i];
        for (int q = 0; q < 3; q++) t1x[i] += ry[i][q] * tx[q];
        for (int c = 0; c < 3; c++) {
            r1[i][c] = 0;
            for (int q = 0; q < 3; q++) r1[i][c] += ry[i][q] * rx[q][c];
        }
    }
    for (int i = 0; i < 3; i++) {
        tz[i] = t1x[i];
        for (int q = 0; q < 3; q++) tz[i] += r1[i][q] * ti[q];
        for (int c = 0; c < 3; c++) {
            rz[i][c] = 0;
            for (int q = 0; q < 3; q++) rz[i][c] += r1[i][q] * ri[q][c];
        }
    }
    for (int q = 0; q < L->nlittle; q++) {
        int p = L->little[q];
        if (memcmp(g->ops[p].rot, rz, sizeof rz)) continue;
        int d[3];
        for (int i = 0; i < 3; i++) d[i] = tz[i] - g->ops[p].tran[i];
        int in_lattice = 0;
        for (int c = 0; c < g->ncent && !in_lattice; c++) {
            in_lattice = 1;
            for (int i = 0; i < 3; i++) in_lattice &= ((d[i] - g->cent[c][i]) % 24 + 24) % 24 == 0;
        }
        if (!in_lattice) continue;
        rat v = rat_mul(rat_int(N), kdot(k, d));
        if (v.d != 1) return spfail("conjugate translation off the lattice");
        *out = q * N + (int)(((v.n % N) + N) % N);
        return 0;
    }
    return spfail("conjugate outside the little group");
}

/* ---- OPDs at k != 0, modulo p in one embedding ---- */

typedef struct {
    int r[3][3];
    int t[3]; /* 1/24 units, not reduced */
} sp_aff;

static void aff_mul(const sp_aff *a, const sp_aff *b, sp_aff *out)
{
    sp_aff z;
    for (int i = 0; i < 3; i++) {
        z.t[i] = a->t[i];
        for (int q = 0; q < 3; q++) z.t[i] += a->r[i][q] * b->t[q];
        for (int c = 0; c < 3; c++) {
            z.r[i][c] = 0;
            for (int q = 0; q < 3; q++) z.r[i][c] += a->r[i][q] * b->r[q][c];
        }
    }
    *out = z;
}

static int p_of_rot(const sp_group *g, const int (*rot)[3])
{
    for (int e = 0; e < g->n; e++)
        if (!memcmp(g->ops[e].rot, rot, sizeof(int[3][3]))) return e;
    return -1;
}

static void aff_inv(const sp_group *g, const sp_aff *a, sp_aff *out)
{
    int e = p_of_rot(g, (const int (*)[3])a->r), ei = g->inv[e];
    memcpy(out->r, g->ops[ei].rot, sizeof out->r);
    for (int i = 0; i < 3; i++) {
        out->t[i] = 0;
        for (int q = 0; q < 3; q++) out->t[i] -= out->r[i][q] * a->t[q];
    }
}

static int in_lattice(const sp_group *g, const int *d)
{
    for (int c = 0; c < g->ncent; c++) {
        int ok = 1;
        for (int i = 0; i < 3; i++) ok &= ((d[i] - g->cent[c][i]) % 24 + 24) % 24 == 0;
        if (ok) return 1;
    }
    return 0;
}

/* the H element of an affine map in G_k */
static int aff_to_h(const sp_group *g, const sp_little *L, const rat *k, const sp_aff *a, int *out)
{
    int e = p_of_rot(g, (const int (*)[3])a->r);
    if (e < 0 || L->little_of[e] < 0) return spfail("element outside the little group");
    int d[3];
    for (int i = 0; i < 3; i++) d[i] = a->t[i] - g->ops[e].tran[i];
    if (!in_lattice(g, d)) return spfail("translation outside the lattice");
    rat v = rat_mul(rat_int(L->N), kdot(k, d));
    if (v.d != 1) return spfail("translation phase outside Z_N");
    *out = L->little_of[e] * L->N + (int)(((v.n % L->N) + L->N) % L->N);
    return 0;
}

static uint64_t cy_mod(const cy_field *f, const cy *x, uint64_t w)
{
    uint64_t v = 0, z = 1;
    for (int i = 0; i < f->phi; i++) {
        i128 c = x->c[i] % (i128)PR;
        if (c < 0) c += PR;
        v = padd(v, pmul((uint64_t)c, z));
        z = pmul(z, w);
    }
    i128 dd = x->den % (i128)PR;
    if (dd < 0) dd += PR;
    return pmul(v, pinv((uint64_t)dd));
}

/* add row v (len) to a reduced echelon basis mod p; returns 1 if independent */
static int padd_vector(uint64_t *basis, int *piv, int *rank, int len, const uint64_t *v0)
{
    uint64_t *v = malloc(sizeof(uint64_t) * (size_t)len);
    memcpy(v, v0, sizeof(uint64_t) * (size_t)len);
    for (int r = 0; r < *rank; r++) {
        uint64_t c = v[piv[r]];
        if (!c) continue;
        for (int j = 0; j < len; j++) v[j] = psub(v[j], pmul(c, basis[(size_t)r * len + j]));
    }
    int p = -1;
    for (int j = 0; j < len && p < 0; j++)
        if (v[j]) p = j;
    if (p < 0) {
        free(v);
        return 0;
    }
    uint64_t sc = pinv(v[p]);
    for (int j = 0; j < len; j++) v[j] = pmul(v[j], sc);
    for (int r = 0; r < *rank; r++) {
        uint64_t c = basis[(size_t)r * len + p];
        if (!c) continue;
        for (int j = 0; j < len; j++) basis[(size_t)r * len + j] = psub(basis[(size_t)r * len + j], pmul(c, v[j]));
    }
    memcpy(basis + (size_t)*rank * len, v, sizeof(uint64_t) * (size_t)len);
    piv[(*rank)++] = p;
    free(v);
    return 1;
}

/* canonical rows of a subspace: reduced echelon sorted by pivot */
static int pcanon(uint64_t *rows, int r, int d)
{
    uint64_t *basis = malloc(sizeof(uint64_t) * (size_t)(r + 1) * d);
    int *piv = malloc(sizeof(int) * (size_t)(r + 1)), rank = 0;
    for (int i = 0; i < r; i++) padd_vector(basis, piv, &rank, d, rows + (size_t)i * d);
    for (int i = 0; i < rank; i++)
        for (int j = i + 1; j < rank; j++)
            if (piv[j] < piv[i]) {
                int t = piv[i];
                piv[i] = piv[j];
                piv[j] = t;
                for (int x = 0; x < d; x++) {
                    uint64_t u = basis[(size_t)i * d + x];
                    basis[(size_t)i * d + x] = basis[(size_t)j * d + x];
                    basis[(size_t)j * d + x] = u;
                }
            }
    memcpy(rows, basis, sizeof(uint64_t) * (size_t)rank * d);
    free(basis);
    free(piv);
    return rank;
}

typedef struct {
    int width; /* words per item */
    size_t count, cap, slots;
    uint64_t *items;
    int64_t *table;
} pset;

static uint64_t phash(const uint64_t *v, int n)
{
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < n; i++) h = (h ^ v[i]) * 1099511628211ULL;
    return h ^ (h >> 31);
}

/* index of v in the set, inserting it when absent (*added = 1) */
static int64_t pset_put(pset *s, const uint64_t *v, int *added)
{
    if (!s->slots || s->count * 2 >= s->slots) {
        size_t slots = s->slots ? s->slots * 2 : 1024;
        int64_t *table = malloc(sizeof(int64_t) * slots);
        if (!table) return -1;
        for (size_t i = 0; i < slots; i++) table[i] = -1;
        for (size_t i = 0; i < s->count; i++) {
            size_t slot = phash(s->items + i * s->width, s->width) & (slots - 1);
            while (table[slot] >= 0) slot = (slot + 1) & (slots - 1);
            table[slot] = (int64_t)i;
        }
        free(s->table);
        s->table = table;
        s->slots = slots;
    }
    size_t slot = phash(v, s->width) & (s->slots - 1);
    for (; s->table[slot] >= 0; slot = (slot + 1) & (s->slots - 1))
        if (!memcmp(s->items + (size_t)s->table[slot] * s->width, v, sizeof(uint64_t) * (size_t)s->width)) {
            *added = 0;
            return s->table[slot];
        }
    if (s->count == s->cap) {
        size_t cap = s->cap ? s->cap * 2 : 256;
        uint64_t *items = realloc(s->items, sizeof(uint64_t) * cap * s->width);
        if (!items) return -1;
        s->items = items;
        s->cap = cap;
    }
    memcpy(s->items + s->count * s->width, v, sizeof(uint64_t) * (size_t)s->width);
    s->table[slot] = (int64_t)s->count;
    *added = 1;
    return (int64_t)s->count++;
}

static void pset_free(pset *s)
{
    free(s->items);
    free(s->table);
}

/* rows x -> x M^T for d x d M */
static void prow_image(const uint64_t *x, int r, const uint64_t *mt, int d, uint64_t *out)
{
    for (int i = 0; i < r; i++)
        for (int a = 0; a < d; a++) {
            uint64_t t = 0;
            for (int q = 0; q < d; q++)
                if (x[(size_t)i * d + q] && mt[(size_t)a * d + q]) t = padd(t, pmul(mt[(size_t)a * d + q], x[(size_t)i * d + q]));
            out[(size_t)i * d + a] = t;
        }
}

/* the OPDs of Gamma = Ind sigma: one copy of sigma in the site's Bloch displacements, the induced matrices on the
 * image of G, and the closure of fixed-space intersections (complex dimensions) */
static int det3(const int (*r)[3])
{
    return r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1]) - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0]) +
           r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
}

/* the strain action eps -> R eps R^T on (11, 22, 33, 23, 13, 12) */
static void strain_matrix(const int (*r)[3], int out[6][6])
{
    static const int ij[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
    for (int c = 0; c < 6; c++) {
        int e[3][3] = {{0}}, t[3][3];
        e[ij[c][0]][ij[c][1]] = e[ij[c][1]][ij[c][0]] = 1;
        for (int a = 0; a < 3; a++)
            for (int b2 = 0; b2 < 3; b2++) {
                t[a][b2] = 0;
                for (int p2 = 0; p2 < 3; p2++)
                    for (int q = 0; q < 3; q++) t[a][b2] += r[a][p2] * e[p2][q] * r[b2][q];
            }
        for (int row = 0; row < 6; row++) out[row][c] = t[ij[row][0]][ij[row][1]];
    }
}

/* site >= 0: the Bloch displacements (axial: moments, times det R) of one site; site < 0: the strain (k = 0) */
static int k_opds(const sp_group *g, const sp_little *L, const rat *k, const cy_field *f, uint64_t w, const cy *sigma, int dim, int site,
                  int axial, int magnetic, int nu, jbuf *b)
{
    const sp_group *h = &L->h;
    int nh = h->n, N = L->N, m = f->m, n = g->n;
    uint64_t zN = powmod(w, (uint64_t)(m / N), PR);
    uint64_t *sig = malloc(sizeof(uint64_t) * (size_t)nh);
    for (int x = 0; x < nh; x++) sig[x] = cy_mod(f, &sigma[h->cls[x]], w);
    int nloc = 0, *local = malloc(sizeof(int) * (size_t)g->nclasses), *global = malloc(sizeof(int) * (size_t)g->nclasses);
    for (int cl = 0; cl < g->nclasses; cl++) {
        local[cl] = -1;
        if (g->class_site[cl] == site) {
            local[cl] = nloc;
            global[nloc++] = cl;
        }
    }
    int len = site < 0 ? 6 : 3 * nloc;
    if (site < 0) nloc = 2;
    int *sdense = site < 0 ? malloc(sizeof(int) * 36 * (size_t)nh) : 0;
    for (int x = 0; sdense && x < nh; x++) strain_matrix((const int (*)[3])g->ops[L->little[x / N]].rot, (int (*)[6])(sdense + 36 * x));
    /* rho(x) on the Bloch displacements: block (gx a, a) = zeta_N^(-j - N k.L) R */
    uint64_t *phase = malloc(sizeof(uint64_t) * (size_t)nh * nloc);
    int *target = malloc(sizeof(int) * (size_t)nh * nloc);
    for (int x = 0; x < nh && site >= 0; x++) {
        int e = L->little[x / N], j = x % N, sign = axial ? det3((const int (*)[3])g->ops[e].rot) : 1;
        for (int c = 0; c < nloc; c++) {
            int cl = global[c];
            rat ph = rat_mul(rat_int(N), kdot(k, g->shift + (e * g->nclasses + cl) * 3));
            int ex = (int)((((-(int64_t)ph.n - j) % N) + N) % N);
            phase[x * nloc + c] = sign < 0 ? psub(0, powmod(zN, (uint64_t)ex, PR)) : powmod(zN, (uint64_t)ex, PR);
            target[x * nloc + c] = local[g->perm[e * g->nclasses + cl]];
        }
    }
#define ACT(x, v, u)                                                                                           \
    do {                                                                                                       \
        memset(u, 0, sizeof(uint64_t) * (size_t)len);                                                          \
        if (sdense) {                                                                                          \
            for (int a_ = 0; a_ < 6; a_++)                                                                     \
                for (int b_ = 0; b_ < 6; b_++)                                                                 \
                    if (sdense[36 * (x) + 6 * a_ + b_]) u[a_] = padd(u[a_], pmul(pint(sdense[36 * (x) + 6 * a_ + b_]), (v)[b_])); \
            break;                                                                                             \
        }                                                                                                      \
        int e_ = L->little[(x) / N];                                                                           \
        for (int c_ = 0; c_ < nloc; c_++) {                                                                    \
            int t_ = target[(x) * nloc + c_];                                                                  \
            for (int a_ = 0; a_ < 3; a_++) {                                                                   \
                uint64_t s_ = 0;                                                                               \
                for (int b_ = 0; b_ < 3; b_++)                                                                 \
                    if (g->ops[e_].rot[a_][b_]) s_ = padd(s_, pmul(pint(g->ops[e_].rot[a_][b_]), (v)[3 * c_ + b_])); \
                u[3 * t_ + a_] = pmul(s_, phase[(x) * nloc + c_]);                                             \
            }                                                                                                  \
        }                                                                                                      \
    } while (0)
    uint64_t *v0 = malloc(sizeof(uint64_t) * (size_t)len), *v = malloc(sizeof(uint64_t) * (size_t)len), *u = malloc(sizeof(uint64_t) * (size_t)len);
    uint64_t *e0 = malloc(sizeof(uint64_t) * (size_t)len), *basis = malloc(sizeof(uint64_t) * (size_t)(dim + 1) * len);
    int *piv = malloc(sizeof(int) * (size_t)(dim + 2)), rank = 0, found = 0;
    /* candidate subgroups K: trivial, cyclic, and pairs of lifted little elements; keep those with dim sigma^K = 1 */
    /* first: an element x and an eigenvalue zeta_o^j of sigma(x) of multiplicity one; the eigenvectors of the
     * isotypic part for it are pure tensors */
    for (int x = 0; x < nh && !found; x++) {
        int o = h->order[x];
        if (m % o) continue;
        uint64_t zo = powmod(w, (uint64_t)(m / o), PR);
        int *pw = malloc(sizeof(int) * (size_t)o);
        int identity = 0;
        for (int y = 0; y < nh; y++)
            if (h->order[y] == 1) identity = y;
        pw[0] = identity;
        for (int t2 = 1; t2 < o; t2++) pw[t2] = h->mult[pw[t2 - 1] * nh + x];
        for (int j = 0; j < o && !found; j++) {
            uint64_t mu = 0, zinv = pinv(powmod(zo, (uint64_t)j, PR)), z = 1;
            for (int t2 = 0; t2 < o; t2++) {
                mu = padd(mu, pmul(sig[pw[t2]], z));
                z = pmul(z, zinv);
            }
            mu = pmul(mu, pinv(pint(o)));
            if (mu != 1) continue;
            for (int start = 0; start < len && !found; start++) {
                memset(e0, 0, sizeof(uint64_t) * (size_t)len);
                e0[start] = 1;
                memset(v0, 0, sizeof(uint64_t) * (size_t)len);
                for (int y = 0; y < nh; y++) {
                    ACT(y, e0, u);
                    uint64_t c = sig[h->inv[y]];
                    for (int i = 0; i < len; i++) v0[i] = padd(v0[i], pmul(c, u[i]));
                }
                memset(v, 0, sizeof(uint64_t) * (size_t)len);
                z = 1;
                for (int t2 = 0; t2 < o; t2++) {
                    ACT(pw[t2], v0, u);
                    for (int i = 0; i < len; i++) v[i] = padd(v[i], pmul(z, u[i]));
                    z = pmul(z, zinv);
                }
                int nonzero = 0;
                for (int i = 0; i < len; i++) nonzero |= v[i] != 0;
                if (!nonzero) continue;
                rank = 0;
                int ok = 1;
                for (int y = 0; y < nh && ok; y++) {
                    ACT(y, v, u);
                    padd_vector(basis, piv, &rank, len, u);
                    ok = rank <= dim;
                }
                found = ok && rank == dim;
            }
        }
        free(pw);
    }
    /* pairs over all of H when it is small, else over the lifts with every phase */
    int pool = nh <= 192 ? nh : L->nlittle, pool_step = nh <= 192 ? 1 : N;
    int64_t ncand = 1 + nh + (int64_t)pool * pool;
    for (int64_t cand = 0; cand < ncand && !found; cand++) {
        unsigned char *in = calloc((size_t)nh, 1);
        int gens[2] = {-1, -1};
        if (cand > 0 && cand <= nh) gens[0] = (int)cand - 1;
        else if (cand > nh) {
            int a = (int)((cand - nh - 1) / pool), c2 = (int)((cand - nh - 1) % pool);
            if (c2 <= a) {
                free(in);
                continue;
            }
            gens[0] = a * pool_step;
            gens[1] = c2 * pool_step;
        }
        int identity = 0;
        for (int x = 0; x < nh; x++)
            if (h->order[x] == 1) identity = x;
        in[identity] = 1;
        int grown = 1;
        while (grown) {
            grown = 0;
            for (int x = 0; x < nh; x++)
                if (in[x])
                    for (int q = 0; q < 2; q++)
                        if (gens[q] >= 0 && !in[h->mult[x * nh + gens[q]]]) {
                            in[h->mult[x * nh + gens[q]]] = 1;
                            grown = 1;
                        }
        }
        uint64_t fix = 0;
        int size = 0;
        for (int x = 0; x < nh; x++)
            if (in[x]) {
                fix = padd(fix, sig[x]);
                size++;
            }
        fix = pmul(fix, pinv(pint(size)));
        if (fix != 1) {
            free(in);
            continue;
        }
        for (int start = 0; start < len && !found; start++) {
            memset(e0, 0, sizeof(uint64_t) * (size_t)len);
            e0[start] = 1;
            memset(v0, 0, sizeof(uint64_t) * (size_t)len);
            for (int x = 0; x < nh; x++) {
                ACT(x, e0, u);
                uint64_t c = sig[h->inv[x]];
                for (int i = 0; i < len; i++) v0[i] = padd(v0[i], pmul(c, u[i]));
            }
            memset(v, 0, sizeof(uint64_t) * (size_t)len);
            for (int x = 0; x < nh; x++)
                if (in[x]) {
                    ACT(x, v0, u);
                    for (int i = 0; i < len; i++) v[i] = padd(v[i], u[i]);
                }
            int nonzero = 0;
            for (int i = 0; i < len; i++) nonzero |= v[i] != 0;
            if (!nonzero) continue;
            rank = 0;
            int ok = 1;
            for (int x = 0; x < nh && ok; x++) {
                ACT(x, v, u);
                padd_vector(basis, piv, &rank, len, u);
                ok = rank <= dim;
            }
            found = ok && rank == dim;
        }
        free(in);
    }
    if (!found) {
        jb_lit(b, ", \"opds\": null");
        free(sig), free(local), free(global), free(phase), free(target), free(v0), free(v), free(u), free(e0), free(basis), free(piv);
        return 0;
    }
    /* sigma's matrices on the copy: coordinates at the pivots */
    uint64_t *smat = malloc(sizeof(uint64_t) * (size_t)nh * dim * dim);
    for (int x = 0; x < nh; x++)
        for (int c = 0; c < dim; c++) {
            ACT(x, basis + (size_t)c * len, u);
            for (int r = 0; r < dim; r++) smat[((size_t)x * dim + r) * dim + c] = u[piv[r]];
        }
#undef ACT
    /* arms and coset representatives */
    int s = L->star, *rep = malloc(sizeof(int) * (size_t)s), narms = 0;
    rat (*arms)[3] = malloc(sizeof(rat[3]) * (size_t)s);
    for (int e = 0; e < n && narms < s; e++) {
        rat a[3];
        k_times(k, (const int (*)[3])g->ops[g->inv[e]].rot, a);
        int known = 0;
        for (int q = 0; q < narms && !known; q++) known = k_same(g, arms[q], a);
        if (!known) {
            memcpy(arms[narms], a, sizeof a);
            rep[narms++] = e;
        }
    }
    int D = dim * s;
    /* Gamma(x) for an affine x: block (pi(i), i) = sigma(r_pi(i)^-1 x r_i) */
    sp_aff *raff = malloc(sizeof(sp_aff) * (size_t)s), *rinv = malloc(sizeof(sp_aff) * (size_t)s);
    for (int i = 0; i < s; i++) {
        memcpy(raff[i].r, g->ops[rep[i]].rot, sizeof raff[i].r);
        memcpy(raff[i].t, g->ops[rep[i]].tran, sizeof raff[i].t);
        aff_inv(g, &raff[i], &rinv[i]);
    }
    int ngens = 0;
    sp_aff gens[64];
    /* generators: a generating set of P (greedy) and the lattice translations */
    {
        unsigned char *in = calloc((size_t)n, 1);
        int identity = 0;
        for (int e = 0; e < n; e++)
            if (g->order[e] == 1) identity = e;
        in[identity] = 1;
        for (int e = 0; e < n && ngens < 48; e++) {
            if (in[e]) continue;
            memcpy(gens[ngens].r, g->ops[e].rot, sizeof gens[ngens].r);
            memcpy(gens[ngens].t, g->ops[e].tran, sizeof gens[ngens].t);
            ngens++;
            int grown = 1;
            in[e] = 1;
            while (grown) {
                grown = 0;
                for (int x = 0; x < n; x++)
                    if (in[x])
                        for (int y = 0; y < n; y++)
                            if (in[y] && !in[g->mult[x * n + y]]) {
                                in[g->mult[x * n + y]] = 1;
                                grown = 1;
                            }
            }
        }
        free(in);
        for (int a = 0; a < 3 + g->ncent && ngens < 64; a++) {
            memset(&gens[ngens], 0, sizeof gens[ngens]);
            for (int i = 0; i < 3; i++) gens[ngens].r[i][i] = 1;
            if (a < 3) gens[ngens].t[a] = 24;
            else memcpy(gens[ngens].t, g->cent[a - 3], sizeof gens[ngens].t);
            ngens++;
        }
    }
    uint64_t *gmat = calloc((size_t)(ngens + 1) * D * D, sizeof(uint64_t));
    for (int q = 0; q < ngens; q++)
        for (int i = 0; i < s; i++) {
            sp_aff y, hh;
            aff_mul(&gens[q], &raff[i], &y);
            rat a[3];
            int ey = p_of_rot(g, (const int (*)[3])y.r);
            k_times(k, (const int (*)[3])g->ops[g->inv[ey]].rot, a);
            int pi = -1;
            for (int q2 = 0; q2 < s && pi < 0; q2++)
                if (k_same(g, arms[q2], a)) pi = q2;
            if (pi < 0) return spfail("arm image outside the star"), MDE_VALUE;
            aff_mul(&rinv[pi], &y, &hh);
            int x;
            if (aff_to_h(g, L, k, &hh, &x)) return MDE_VALUE;
            for (int r = 0; r < dim; r++)
                for (int c = 0; c < dim; c++) gmat[((size_t)q * D + pi * dim + r) * D + i * dim + c] = smat[((size_t)x * dim + r) * dim + c];
        }
    if (magnetic) {
        for (int i = 0; i < D; i++) gmat[((size_t)ngens * D + i) * D + i] = psub(0, 1);
        ngens++;
    }
    /* the image group by closure */
    pset image = {D * D, 0, 0, 0, 0, 0};
    uint64_t *ident = calloc((size_t)D * D, sizeof(uint64_t)), *prod = malloc(sizeof(uint64_t) * (size_t)D * D);
    for (int i = 0; i < D; i++) ident[i * D + i] = 1;
    int added;
    pset_put(&image, ident, &added);
    for (size_t at = 0; at < image.count; at++) {
        if (image.count * (size_t)image.width > ((size_t)40 << 20)) return spfail("image group above the reference limit"), MDE_VALUE;
        for (int q = 0; q < ngens; q++) {
            const uint64_t *a = gmat + (size_t)q * D * D;
            for (int i = 0; i < D; i++)
                for (int c = 0; c < D; c++) {
                    uint64_t t2 = 0;
                    for (int l = 0; l < D; l++) {
                        uint64_t x1 = a[i * D + l], x2 = image.items[at * image.width + (size_t)l * D + c];
                        if (x1 && x2) t2 = padd(t2, pmul(x1, x2));
                    }
                    prod[i * D + c] = t2;
                }
            if (pset_put(&image, prod, &added) < 0) return spfail("out of memory"), MDE_MEMORY;
        }
    }
    size_t order = image.count;
    /* closure of isotropy subspaces; subspaces as canonical rows, key = rank then rows */
    pset seen = {1 + D * D, 0, 0, 0, 0, 0};
    uint64_t *key = calloc((size_t)1 + D * D, sizeof(uint64_t)), *rows = malloc(sizeof(uint64_t) * (size_t)D * D);
    uint64_t *img = malloc(sizeof(uint64_t) * (size_t)D * D), *diff = malloc(sizeof(uint64_t) * (size_t)D * D);
    size_t rcap = 1024, nreps = 0;
    uint64_t *reps = malloc(sizeof(uint64_t) * rcap * (1 + D * D));
    memcpy(rows, ident, sizeof(uint64_t) * (size_t)D * D);
    key[0] = (uint64_t)D;
    memcpy(key + 1, rows, sizeof(uint64_t) * (size_t)D * D);
    memcpy(reps, key, sizeof(uint64_t) * (size_t)(1 + D * D));
    nreps = 1;
    pset_put(&seen, key, &added);
    for (size_t head = 0; head < nreps; head++) {
        int r = (int)reps[head * (1 + D * D)];
        const uint64_t *x = reps + head * (1 + D * D) + 1;
        for (size_t q = 0; q < order; q++) {
            const uint64_t *mat = image.items + q * image.width;
            /* D' = X M^T - X; coefficients c with c D' = 0 */
            prow_image(x, r, mat, D, img);
            for (int i = 0; i < r * D; i++) diff[i] = psub(img[i], x[i]);
            uint64_t *tr = malloc(sizeof(uint64_t) * (size_t)D * r), *cs = malloc(sizeof(uint64_t) * (size_t)r * r);
            for (int i = 0; i < r; i++)
                for (int a = 0; a < D; a++) tr[a * r + i] = diff[i * D + a];
            int cnt = pnull(tr, D, r, cs);
            free(tr);
            if (!cnt || cnt == r) {
                free(cs);
                continue;
            }
            memset(rows, 0, sizeof(uint64_t) * (size_t)cnt * D);
            for (int v2 = 0; v2 < cnt; v2++)
                for (int i = 0; i < r; i++)
                    if (cs[v2 * r + i])
                        for (int a = 0; a < D; a++) rows[v2 * D + a] = padd(rows[v2 * D + a], pmul(cs[v2 * r + i], x[i * D + a]));
            free(cs);
            int rr = pcanon(rows, cnt, D);
            memset(key, 0, sizeof(uint64_t) * (size_t)(1 + D * D));
            key[0] = (uint64_t)rr;
            memcpy(key + 1, rows, sizeof(uint64_t) * (size_t)rr * D);
            pset_put(&seen, key, &added);
            if (!added) continue;
            if (seen.count * (size_t)seen.width > ((size_t)40 << 20)) return spfail("isotropy orbits above the reference limit"), MDE_VALUE;
            if (nreps == rcap) {
                rcap *= 2;
                reps = realloc(reps, sizeof(uint64_t) * rcap * (1 + D * D));
            }
            if (nreps > 20000) return spfail("too many isotropy subspaces"), MDE_VALUE;
            memcpy(reps + nreps * (1 + D * D), key, sizeof(uint64_t) * (size_t)(1 + D * D));
            nreps++;
            /* the whole orbit is seen */
            uint64_t *orbit = malloc(sizeof(uint64_t) * (size_t)(1 + D * D));
            for (size_t q2 = 0; q2 < order; q2++) {
                prow_image(rows, rr, image.items + q2 * image.width, D, img);
                int r2 = pcanon(img, rr, D);
                memset(orbit, 0, sizeof(uint64_t) * (size_t)(1 + D * D));
                orbit[0] = (uint64_t)r2;
                memcpy(orbit + 1, img, sizeof(uint64_t) * (size_t)r2 * D);
                pset_put(&seen, orbit, &added);
            }
            free(orbit);
        }
    }
    jb_lit(b, ", \"image\": ");
    char t[96];
    snprintf(t, sizeof t, "%zu, \"opds\": [", order);
    jb_lit(b, t);
    for (size_t q = 0; q < nreps; q++) {
        int r = (int)reps[q * (1 + D * D)];
        const uint64_t *x = reps + q * (1 + D * D) + 1;
        size_t stab = 0;
        for (size_t q2 = 0; q2 < order; q2++) {
            prow_image(x, r, image.items + q2 * image.width, D, img);
            stab += !memcmp(img, x, sizeof(uint64_t) * (size_t)r * D);
        }
        snprintf(t, sizeof t, q ? ", {\"dim\": %d, \"index\": %zu}" : "{\"dim\": %d, \"index\": %zu}", r * (nu == 1 ? 1 : 2), order / stab);
        jb_lit(b, t);
    }
    jb_lit(b, "]");
    pset_free(&image);
    pset_free(&seen);
    free(key), free(rows), free(img), free(diff), free(reps), free(ident), free(prod), free(gmat), free(raff), free(rinv), free(rep), free(arms);
    free(smat), free(sig), free(local), free(global), free(phase), free(target), free(v0), free(v), free(u), free(e0), free(basis), free(piv), free(sdense);
    return 0;
}

/* irreps at k != 0: per small irrep of H appearing in the displacements, its Herring type, the
 * full and real dimensions and the multiplicity per site */
static int spec_k(const sp_group *g, const rat *k, int nsites, const int *kinds, int strain, jbuf *b)
{
    sp_little L;
    if (build_little(g, k, &L)) return MDE_VALUE;
    const sp_group *h = &L.h;
    int nh = h->n, N = L.N, exponent = 1;
    for (int e = 0; e < nh; e++) {
        int a = exponent, c = h->order[e];
        while (c) {
            int t = a % c;
            a = c;
            c = t;
        }
        exponent = exponent / a * h->order[e];
    }
    int m = exponent % 4 ? exponent * (exponent % 2 ? 4 : 2) : exponent;
    const cy_field *f = cy_field_for(m);
    uint64_t w;
    const int *units;
    const uint64_t *vinv;
    if (!f || modp_field(m, f->phi, 0, &PR, &w, &units, &vinv)) return spfail("no prime for the character field"), MDE_VALUE;
    w = powmod(w, (uint64_t)units[0], PR);
    sp_chars ch;
    if (characters(h, f, w, &ch)) return MDE_VALUE;
    int kc = h->ncls, identity = 0;
    for (int e = 0; e < nh; e++)
        if (h->order[e] == 1) identity = e;
    /* the displacement character per site on H: (g, j) -> zeta_N^(-j) sum_{a fixed} chi(-L_{g,a}) tr R */
    /* rows: kind (0 displacements, 1 moments: tr R times det R) x site, then the strain (k = 0) */
    cy *vchar = calloc((size_t)(2 * nsites + 1) * nh, sizeof(cy));
    if (!vchar) return spfail("out of memory"), MDE_MEMORY;
    for (int row = 0; row < 2 * nsites; row++)
        for (int x = 0; x < nh; x++) {
            int s = row % nsites, kind = row / nsites;
            int e = L.little[x / N], j = x % N;
            cy acc, term, root;
            cy_zero(f, &acc);
            for (int cl = 0; cl < g->nclasses; cl++) {
                if (g->class_site[cl] != s || g->perm[e * g->nclasses + cl] != cl) continue;
                rat ph = rat_mul(rat_int(N), kdot(k, g->shift + (e * g->nclasses + cl) * 3));
                if (ph.d != 1) return spfail("atom shift outside the lattice"), MDE_VALUE;
                int tr = g->ops[e].rot[0][0] + g->ops[e].rot[1][1] + g->ops[e].rot[2][2];
                if (kind) tr *= det3((const int (*)[3])g->ops[e].rot);
                cy_root(f, -(int64_t)(ph.n + j) * (m / N), &root);
                cy_rat(f, tr, 1, &term);
                cy_mul(f, &term, &root, &term);
                cy_add(f, &acc, &term, &acc);
            }
            vchar[(size_t)row * nh + x] = acc;
        }
    if (strain && N == 1)
        for (int x = 0; x < nh; x++) {
            int sm[6][6], tr = 0;
            strain_matrix((const int (*)[3])g->ops[L.little[x]].rot, sm);
            for (int i = 0; i < 6; i++) tr += sm[i][i];
            cy_rat(f, tr, 1, &vchar[(size_t)2 * nsites * nh + x]);
        }
    /* elements of P taking k to -k, for Herring's sum */
    rat minus[3];
    for (int i = 0; i < 3; i++) minus[i] = rat_neg(k[i]);
    int *flip = malloc(sizeof(int) * (size_t)g->n), nflip = 0;
    for (int e = 0; e < g->n; e++) {
        rat kr[3];
        k_times(k, (const int (*)[3])g->ops[e].rot, kr);
        if (k_same(g, kr, minus)) flip[nflip++] = e;
    }
    jb_reset(b);
    char t[200];
    snprintf(t, sizeof t, "{\"schema\": \"APOSTRUCT.spec.k\", \"order\": %d, \"little\": %d, \"N\": %d, \"star\": %d, \"minus_k_in_star\": %s, \"irreps\": [",
             g->n, L.nlittle, N, L.star, nflip ? "true" : "false");
    jb_lit(b, t);
    int written = 0, *skip = calloc((size_t)kc, sizeof(int));
    for (int c = 0; c < kc; c++) {
        if (skip[c]) continue;
        const cy *sigma = ch.values + (size_t)c * kc;
        int dim = 0;
        {
            rat r;
            if (rat_of_cy(f, &sigma[h->cls[identity]], &r)) return spfail("irrational degree"), MDE_VALUE;
            dim = (int)r.n;
        }
        /* multiplicity per site: (1/|H|) sum chi_V(x) sigma(x^-1) */
        /* multiplicities: (1/|H|) sum chi_V(x) sigma(x^-1), per kind and site, and in the strain */
        int *mult = calloc((size_t)(2 * nsites + 1) + 1, sizeof(int)), any = 0;
        for (int row = 0; row <= 2 * nsites; row++) {
            int s2 = row % (nsites ? nsites : 1), kind = row == 2 * nsites ? 2 : row / nsites;
            if (kind < 2 && !(kinds[s2] & (1 << kind))) continue;
            if (kind == 2 && !(strain && N == 1)) continue;
            cy acc, term;
            cy_zero(f, &acc);
            for (int x = 0; x < nh; x++) {
                cy_mul(f, &vchar[(size_t)row * nh + x], &sigma[h->cls[h->inv[x]]], &term);
                cy_add(f, &acc, &term, &acc);
            }
            rat r;
            if (rat_of_cy(f, &acc, &r)) return spfail("irrational multiplicity"), MDE_VALUE;
            r = rat_div(r, rat_int(nh));
            if (r.d != 1 || r.n < 0) return spfail("non-integral multiplicity"), MDE_VALUE;
            mult[row] = (int)r.n;
            any |= mult[row] > 0;
        }
        if (!any) {
            free(mult);
            continue;
        }
        /* Herring: (1/|P_k|) sum over g with g k = -k of sigma((g^2, c(g, g))) */
        int nu = 0, partner = -1;
        if (nflip) {
            cy acc;
            cy_zero(f, &acc);
            for (int q = 0; q < nflip; q++) {
                int e = flip[q], e2 = g->mult[e * g->n + e], c2;
                if (cocycle(g, k, N, e, e, &c2)) return MDE_VALUE;
                int x = L.little_of[e2] * N + c2;
                cy_add(f, &acc, &sigma[h->cls[x]], &acc);
            }
            rat r;
            if (rat_of_cy(f, &acc, &r)) return spfail("irrational Herring sum"), MDE_VALUE;
            r = rat_div(r, rat_int(L.nlittle));
            if (r.d != 1 || r.n < -1 || r.n > 1) return spfail("Herring sum outside {-1, 0, 1}"), MDE_VALUE;
            nu = (int)r.n;
            if (!nu) {
                /* the partner sigma'(x) = conj(sigma(y x y^-1)) for a y taking k to -k pairs with sigma */
                int *conjcls = malloc(sizeof(int) * (size_t)kc);
                for (int l = 0; l < kc; l++) {
                    int z;
                    if (conjugate_in_h(g, &L, k, flip[0], h->cls_rep[l], &z)) return MDE_VALUE;
                    conjcls[l] = h->cls_inv[h->cls[z]];
                }
                for (int d = 0; d < kc && partner < 0; d++) {
                    int same = 1;
                    for (int l = 0; l < kc && same; l++) {
                        cy diff;
                        cy_sub(f, &ch.values[(size_t)d * kc + l], &sigma[conjcls[l]], &diff);
                        same = cy_is_zero(f, &diff);
                    }
                    if (same) partner = d;
                }
                free(conjcls);
                if (partner < 0) return spfail("no partner for a complex small irrep"), MDE_VALUE;
                if (partner != c) skip[partner] = 1;
            }
        }
        int full = dim * L.star, real = nu == 1 ? full : 2 * full;
        for (int kind = 0; kind < 2; kind++) {
            int site = -1, total = 0;
            for (int s2 = 0; s2 < nsites; s2++) {
                total += mult[kind * nsites + s2];
                if (site < 0 && mult[kind * nsites + s2]) site = s2;
            }
            int strained = kind == 0 ? mult[2 * nsites] : 0;
            if (!total && !strained) continue;
            jb_lit(b, written++ ? ", {" : "{");
            snprintf(t, sizeof t, "\"kind\": \"%s\", \"small_dim\": %d, \"full_dim\": %d, \"real_dim\": %d, \"herring\": %d, \"partner\": %s, \"strain\": %d, \"multiplicity\": ",
                     kind ? "magnetic" : "displacive", dim, full, real, nu, partner >= 0 && partner != c ? "true" : "false", strained);
            jb_lit(b, t);
            put_int_list(b, mult + kind * nsites, nsites);
            if (k_opds(g, &L, k, f, w, sigma, dim, site, kind, kind, nu, b)) return MDE_VALUE;
            jb_lit(b, "}");
        }
        free(mult);
    }
    jb_lit(b, "]}");
    free(skip);
    free(flip);
    free(vchar);
    return b->failed ? md_fail(MDE_MEMORY, "out of memory") : 0;
}

int spec_gamma(const jvalue *cif_info, const int *sites, const int *kinds, int nsites, int strain, const rat *kvec, jbuf *b)
{
    fail_status = 0;
    sp_group g;
    memset(&g, 0, sizeof g);
    if (build_group(cif_info, sites, nsites, &g)) return MDE_VALUE;
    return spec_k(&g, kvec, nsites, kinds, strain, b);
}
