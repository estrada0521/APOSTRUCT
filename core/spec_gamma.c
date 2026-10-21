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
#include "geometry.h"
#include "irreps.h"
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

typedef struct {
    const sp_group *g;
    const sp_little *L;
    const rat *k;
    int s, dim, D;
    const rat (*arms)[3];
    const sp_aff *raff, *rinv;
    const uint64_t *smat;
} sp_induced;

/* Gamma(x) = Ind sigma (x): block (pi(i), i) = sigma(r_pi(i)^-1 x r_i) */
static int gamma_of(const sp_induced *c, const sp_aff *x, uint64_t *out)
{
    int D = c->D, dim = c->dim;
    memset(out, 0, sizeof(uint64_t) * (size_t)D * D);
    for (int i = 0; i < c->s; i++) {
        sp_aff y, hh;
        aff_mul(x, &c->raff[i], &y);
        rat a[3];
        int ey = p_of_rot(c->g, (const int (*)[3])y.r);
        if (ey < 0) return spfail("rotation outside the group");
        k_times(c->k, (const int (*)[3])c->g->ops[c->g->inv[ey]].rot, a);
        int pi = -1;
        for (int q2 = 0; q2 < c->s && pi < 0; q2++)
            if (k_same(c->g, c->arms[q2], a)) pi = q2;
        if (pi < 0) return spfail("arm image outside the star");
        aff_mul(&c->rinv[pi], &y, &hh);
        int h;
        if (aff_to_h(c->g, c->L, c->k, &hh, &h)) return 1;
        for (int r = 0; r < dim; r++)
            for (int col = 0; col < dim; col++) out[((size_t)pi * dim + r) * D + i * dim + col] = c->smat[((size_t)h * dim + r) * dim + col];
    }
    return 0;
}

/* ---- naming bridge: H as Source records (pml frame) for id_subgroup ---- */

static int SP_SG, SP_SETTING;
int id_subgroup_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t *out);
int isotropy_row_subgroup(int64_t row_id);

static int64_t igcd(int64_t a, int64_t b)
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

/* rows: a basis of the lattice spanned by n integer vectors (Euclid on columns); returns the rank */
static int lattice_basis(int64_t (*v)[3], int n, int64_t (*out)[3])
{
    int r = 0;
    for (int col = 0; col < 3 && r < n; col++) {
        for (;;) {
            int best = -1;
            for (int i = r; i < n; i++)
                if (v[i][col] && (best < 0 || llabs(v[i][col]) < llabs(v[best][col]))) best = i;
            if (best < 0) break;
            int64_t t[3];
            memcpy(t, v[best], sizeof t);
            memcpy(v[best], v[r], sizeof t);
            memcpy(v[r], t, sizeof t);
            int done = 1;
            for (int i = r + 1; i < n; i++) {
                int64_t q = v[i][col] / v[r][col];
                for (int j = 0; j < 3; j++) v[i][j] -= q * v[r][j];
                if (v[i][col]) done = 0;
            }
            if (done) break;
        }
        if (v[r][col]) r++;
    }
    memcpy(out, v, sizeof(int64_t[3]) * (size_t)r);
    return r;
}

/* the sublattice of basis b (rows, 1/24 units) where k . t is integral */
static int lattice_kernel(int64_t (*b)[3], const rat *k)
{
    int64_t vals[3], d = 1;
    rat f[3];
    for (int j = 0; j < 3; j++) {
        int t[3] = {(int)b[j][0], (int)b[j][1], (int)b[j][2]};
        f[j] = rat_mod1(kdot(k, t));
        d = d / igcd(d, (int64_t)f[j].d) * (int64_t)f[j].d;
    }
    if (d == 1) return 0;
    for (int j = 0; j < 3; j++) vals[j] = (int64_t)(f[j].n * (d / f[j].d)) % d;
    /* unimodular u with vals . u = (gcd, 0, 0) */
    int64_t u[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int j = 1; j < 3; j++)
        while (vals[j]) {
            int64_t q = vals[0] / vals[j];
            vals[0] -= q * vals[j];
            for (int i = 0; i < 3; i++) u[i][0] -= q * u[i][j];
            int64_t t = vals[0];
            vals[0] = vals[j];
            vals[j] = t;
            for (int i = 0; i < 3; i++) {
                t = u[i][0];
                u[i][0] = u[i][j];
                u[i][j] = t;
            }
        }
    int64_t gg = igcd(vals[0], d), scale[3] = {d / gg, 1, 1}, nb[3][3];
    for (int c = 0; c < 3; c++)
        for (int j = 0; j < 3; j++) {
            nb[c][j] = 0;
            for (int i = 0; i < 3; i++) nb[c][j] += scale[c] * u[i][c] * b[i][j];
        }
    memcpy(b, nb, sizeof nb);
    return 0;
}

/* the ordinary subgroup of the OPD with canonical rows x (r x D): records of one element per rotation in H, the
 * lattice of H, both in the pml frame; out as id_subgroup_identify (subgroup, basis, origin); returns found */
int mag_nonmagnetic_op(int op);
int mag_time_reversed(int op);
int mag_tables(void);
int64_t mag_point_product_count(void);
int id_subgroup_magnetic_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t capacity, int64_t *out);

/* magnetic: elements (R, t) with or without time reversal (Gamma negated) that fix U; out[0] the magnetic group */
static int identify_opd(const sp_induced *c, const uint64_t *x, int r, int magnetic, int64_t *out)
{
    const sp_group *g = c->g;
    int D = c->D, dim = c->dim, N = c->L->N;
    /* arms carrying U */
    int carried[64] = {0};
    for (int i = 0; i < c->s; i++)
        for (int q = 0; q < r && !carried[i]; q++)
            for (int a = 0; a < dim && !carried[i]; a++) carried[i] = x[(size_t)q * D + i * dim + a] != 0;
    /* T_H: the lattice where k_i . t is integral on every carried arm */
    int64_t gens[16][3], basis[3][3];
    int ng = 0;
    for (int a = 0; a < 3; a++) {
        memset(gens[ng], 0, sizeof gens[ng]);
        gens[ng++][a] = 24;
    }
    for (int q = 0; q < g->ncent; q++) {
        for (int a = 0; a < 3; a++) gens[ng][a] = g->cent[q][a];
        ng++;
    }
    if (lattice_basis(gens, ng, basis) != 3) return spfail("lattice of rank below 3");
    for (int i = 0; i < c->s; i++)
        if (carried[i]) lattice_kernel(basis, c->arms[i]);
    /* translation classes T / T_H by phase on the carried arms */
    int64_t (*reps)[3] = malloc(sizeof(int64_t[3]) * 4096);
    int nrep = 0;
    int *keys = malloc(sizeof(int) * 4096 * 64);
    for (int q = 0; q < g->ncent && nrep < 4096; q++)
        for (int a = 0; a < N && nrep < 4096; a++)
            for (int b2 = 0; b2 < N && nrep < 4096; b2++)
                for (int c2 = 0; c2 < N && nrep < 4096; c2++) {
                    int t[3] = {24 * a + g->cent[q][0], 24 * b2 + g->cent[q][1], 24 * c2 + g->cent[q][2]}, key[64];
                    for (int i = 0; i < c->s; i++) {
                        rat v = rat_mul(rat_int(N), kdot(c->arms[i], t));
                        key[i] = carried[i] ? (int)(((v.n % N) + N) % N) : 0;
                    }
                    int known = 0;
                    for (int z = 0; z < nrep && !known; z++) known = !memcmp(keys + z * 64, key, sizeof(int) * (size_t)c->s);
                    if (known) continue;
                    memcpy(keys + nrep * 64, key, sizeof(int) * (size_t)c->s);
                    for (int z = 0; z < 3; z++) reps[nrep][z] = t[z];
                    nrep++;
                }
    free(keys);
    /* one element of H per rotation */
    rat m[3][3], origin[3], minv[3][3];
    if (geo_setting_to_cinter(SP_SG, 1, SP_SETTING, m, origin) || rat_inverse3((const rat (*)[3])m, minv)) {
        free(reps);
        return spfail("parent setting tables are inconsistent");
    }
    int64_t records[2 * GEO_MAX_OPS][5], src[GEO_MAX_OPS][5];
    int nsrc = geo_sg_records(SP_SG, src), nkept = 0;
    uint64_t *gm = malloc(sizeof(uint64_t) * (size_t)D * D), *img = malloc(sizeof(uint64_t) * (size_t)r * D);
    if (magnetic && mag_tables()) return spfail("magnetic tables are missing");
    for (int pass = 0; pass < (magnetic ? 2 : 1) * g->n; pass++) {
        int e = pass % g->n, reversed = pass >= g->n, got = -1;
        for (int q = 0; q < nrep && got < 0; q++) {
            sp_aff y;
            memcpy(y.r, g->ops[e].rot, sizeof y.r);
            for (int z = 0; z < 3; z++) y.t[z] = g->ops[e].tran[z] + (int)reps[q][z];
            if (gamma_of(c, &y, gm)) return 1;
            if (reversed)
                for (int i = 0; i < D * D; i++) gm[i] = psub(0, gm[i]);
            prow_image(x, r, gm, D, img);
            if (!memcmp(img, x, sizeof(uint64_t) * (size_t)r * D)) got = q;
        }
        if (got < 0) continue;
        /* pml: point = m R^T m^-1 (rows act), t_pml = (t_C - origin + origin r1) m^-1 with r1 = R^T */
        rat rt[3][3], tmp[3][3], point[3][3], tc[3], o1[3], tp[3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) rt[i][j] = rat_int(g->ops[e].rot[j][i]);
        for (int i = 0; i < 3; i++) rat_row_mul(m[i], (const rat (*)[3])rt, tmp[i]);
        for (int i = 0; i < 3; i++) rat_row_mul(tmp[i], (const rat (*)[3])minv, point[i]);
        rat_row_mul(origin, (const rat (*)[3])rt, o1);
        for (int z = 0; z < 3; z++) tc[z] = rat_add(rat_sub(rat_make(g->ops[e].tran[z] + reps[got][z], 24), origin[z]), o1[z]);
        rat_row_mul(tc, (const rat (*)[3])minv, tp);
        int op = -1;
        for (int q = 0; q < nsrc && op < 0; q++) {
            int same = 1;
            for (int i = 0; i < 3 && same; i++) {
                rat unit[3] = {rat_int(i == 0), rat_int(i == 1), rat_int(i == 2)}, row[3];
                geo_vrot(SP_SG, (int)src[q][4], unit, row);
                for (int j = 0; j < 3; j++) same &= rat_eq(row[j], point[i][j]);
            }
            if (same) op = (int)src[q][4];
        }
        if (op < 0) {
            free(reps), free(gm), free(img);
            return spfail("rotation without a Source point operation");
        }
        int64_t den = 1;
        for (int z = 0; z < 3; z++) den = den / igcd(den, (int64_t)tp[z].d) * (int64_t)tp[z].d;
        for (int z = 0; z < 3; z++) records[nkept][z] = (int64_t)(tp[z].n * (den / tp[z].d));
        records[nkept][3] = den;
        records[nkept][4] = op;
        if (magnetic) {
            int code = 0;
            for (int64_t mc = 1; mc <= mag_point_product_count() && mc <= 1000 && !code; mc++)
                if (mag_nonmagnetic_op((int)mc) == op && mag_time_reversed((int)mc) == reversed) code = (int)mc;
            if (!code) {
                free(reps), free(gm), free(img);
                return spfail("no magnetic point operation");
            }
            records[nkept][4] = code;
        }
        nkept++;
    }
    free(reps), free(gm), free(img);
    int64_t bpml[9];
    for (int i = 0; i < 3; i++) {
        rat row[3], v[3];
        for (int z = 0; z < 3; z++) row[z] = rat_make(basis[i][z], 24);
        rat_row_mul(row, (const rat (*)[3])minv, v);
        for (int z = 0; z < 3; z++) {
            if (v[z].d != 1) return spfail("lattice of H is not integral in the pml frame");
            bpml[i * 3 + z] = (int64_t)v[z].n;
        }
    }
    if (rat_overflow) return spfail("rational overflow in the naming bridge");
    if (magnetic) {
        static int64_t flat[8192];
        int64_t length = id_subgroup_magnetic_identify(SP_SG, bpml, nkept, &records[0][0], 0, 8192, flat);
        if (length == -2 || length < -7) return spfail("magnetic subgroup identification failed");
        if (length <= 0) return 0;
        out[0] = flat[0];
        return 1;
    }
    int found = id_subgroup_identify(SP_SG, bpml, nkept, &records[0][0], 0, out);
    if (found < 0) return spfail("subgroup identification failed");
    return found;
}

/* ---- naming bridge: irrep symbols from characters ---- */


typedef struct {
    int count;
    int gid[64];
    char symbol[64][24];
    int64_t kparam[4];
} sp_names;

static const sp_names *SP_NAMES;

int exact_record_matrix(int gid, int dim, int arms, const int64_t *k, int subtract_generated, int element_count, const int *elements,
                        int position, const int *record, int want_m, const cy_field **field_out, cy *w);
int bw_find_position(int element_count, const int *elements, const int *record);

/* x in Q(zeta_m) as an element of Q(zeta_M), m | M */
static void cy_lift(const cy_field *from, const cy *x, const cy_field *to, cy *out)
{
    cy acc, term, root;
    cy_zero(to, &acc);
    for (int i = 0; i < from->phi; i++) {
        if (!x->c[i]) continue;
        cy_root(to, (int64_t)i * (to->m / from->m), &root);
        cy_rat(to, x->c[i], 1, &term);
        cy_mul(to, &term, &root, &term);
        cy_add(to, &acc, &term, &acc);
    }
    cy scale;
    cy_rat(to, 1, x->den, &scale);
    cy_mul(to, &acc, &scale, out);
}

static void cy_conj(const cy_field *f, const cy *x, cy *out)
{
    cy acc, term, root;
    cy_zero(f, &acc);
    for (int i = 0; i < f->phi; i++) {
        if (!x->c[i]) continue;
        cy_root(f, -(int64_t)i, &root);
        cy_rat(f, x->c[i], 1, &term);
        cy_mul(f, &term, &root, &term);
        cy_add(f, &acc, &term, &acc);
    }
    cy scale;
    cy_rat(f, 1, x->den, &scale);
    cy_mul(f, &acc, &scale, out);
}

/* an affine map of G (CIF frame) as a Source record (pml frame) */
static int aff_to_record(const sp_group *g, const sp_aff *a, int64_t *rec)
{
    rat m[3][3], origin[3], minv[3][3];
    if (geo_setting_to_cinter(SP_SG, 1, SP_SETTING, m, origin) || rat_inverse3((const rat (*)[3])m, minv))
        return spfail("parent setting tables are inconsistent");
    rat rt[3][3], tmp[3][3], point[3][3], tc[3], o1[3], tp[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) rt[i][j] = rat_int(a->r[j][i]);
    for (int i = 0; i < 3; i++) rat_row_mul(m[i], (const rat (*)[3])rt, tmp[i]);
    for (int i = 0; i < 3; i++) rat_row_mul(tmp[i], (const rat (*)[3])minv, point[i]);
    rat_row_mul(origin, (const rat (*)[3])rt, o1);
    for (int z = 0; z < 3; z++) tc[z] = rat_add(rat_sub(rat_make(a->t[z], 24), origin[z]), o1[z]);
    rat_row_mul(tc, (const rat (*)[3])minv, tp);
    int64_t src[GEO_MAX_OPS][5];
    int nsrc = geo_sg_records(SP_SG, src), op = -1;
    for (int q = 0; q < nsrc && op < 0; q++) {
        int same = 1;
        for (int i = 0; i < 3 && same; i++) {
            rat unit[3] = {rat_int(i == 0), rat_int(i == 1), rat_int(i == 2)}, row[3];
            geo_vrot(SP_SG, (int)src[q][4], unit, row);
            for (int j = 0; j < 3; j++) same &= rat_eq(row[j], point[i][j]);
        }
        if (same) op = (int)src[q][4];
    }
    if (op < 0) return spfail("rotation without a Source point operation");
    int64_t den = 1;
    for (int z = 0; z < 3; z++) den = den / igcd(den, (int64_t)tp[z].d) * (int64_t)tp[z].d;
    for (int z = 0; z < 3; z++) rec[z] = (int64_t)(tp[z].n * (den / tp[z].d));
    rec[3] = den;
    rec[4] = op;
    (void)g;
    return rat_overflow ? spfail("rational overflow in the naming bridge") : 0;
}

/* the exact character of Gamma = Ind sigma at x: sum over fixed arms of sigma(r_i^-1 x r_i) */
static int chi_exact(const sp_induced *c, const cy_field *f, const cy *sigma, const sp_aff *x, cy *out)
{
    cy_zero(f, out);
    for (int i = 0; i < c->s; i++) {
        sp_aff y, hh;
        aff_mul(x, &c->raff[i], &y);
        rat a[3];
        int ey = p_of_rot(c->g, (const int (*)[3])y.r);
        k_times(c->k, (const int (*)[3])c->g->ops[c->g->inv[ey]].rot, a);
        if (!k_same(c->g, c->arms[i], a)) continue;
        aff_mul(&c->rinv[i], &y, &hh);
        int h;
        if (aff_to_h(c->g, c->L, c->k, &hh, &h)) return 1;
        cy_add(f, out, &sigma[c->L->h.cls[h]], out);
    }
    return 0;
}

/* the candidate whose Source matrices have the physically irreducible character of Gamma, or -1 */
static int name_irrep(const sp_induced *c, const cy_field *f, const cy *sigma, int nu, int real_dim, const sp_names *names, char *why)
{
    const sp_group *g = c->g;
    int found = -1, matches = 0, alive[64];
    for (int q = 0; q < names->count; q++) alive[q] = 1;
    /* translations from {0, 1}^3 first; a wider box only when several candidates remain */
    for (int pass = 0; pass < 2; pass++) {
    int side = pass ? (c->L->N < 6 ? c->L->N : 6) : 2;
    if (side < 2) side = 2;
    found = -1;
    matches = 0;
    for (int q = 0; q < names->count; q++) {
        if (!alive[q]) continue;
        irrep_context ctx;
        if (irrep_context_for(names->gid[q], 4, names->kparam, &ctx) || ctx.dim != real_dim) {
            alive[q] = 0;
            continue;
        }
        int ok = 1;
        cy *w = malloc(sizeof(cy) * (size_t)ctx.dim * ctx.dim);
        for (int e = 0; e < g->n && ok; e++)
            for (int tcode = 0; tcode < side * side * side && ok; tcode++)
                for (int cc = 0; cc < g->ncent && ok; cc++) {
                    sp_aff x;
                    int tv[3] = {tcode % side, (tcode / side) % side, tcode / (side * side)};
                    memcpy(x.r, g->ops[e].rot, sizeof x.r);
                    for (int z = 0; z < 3; z++) x.t[z] = g->ops[e].tran[z] + 24 * tv[z] + g->cent[cc][z];
                    int64_t rec64[5];
                    int rec[5];
                    if (aff_to_record(g, &x, rec64)) {
                        free(w);
                        return -1;
                    }
                    for (int z = 0; z < 5; z++) rec[z] = (int)rec64[z];
                    int position = bw_find_position(ctx.element_count, ctx.elements, rec);
                    const cy_field *fs;
                    int need = -exact_record_matrix(names->gid[q], ctx.dim, ctx.arms, ctx.k, ctx.subtract, ctx.element_count, ctx.elements,
                                                    position, rec, -1, &fs, 0);
                    if (!position || need <= 0) {
                        ok = 0;
                        break;
                    }
                    int M = need;
                    {
                        int a = M, b2 = f->m;
                        while (b2) {
                            int t = a % b2;
                            a = b2;
                            b2 = t;
                        }
                        M = M / a * f->m;
                    }
                    const cy_field *fm = cy_field_for(M);
                    if (!fm || exact_record_matrix(names->gid[q], ctx.dim, ctx.arms, ctx.k, ctx.subtract, ctx.element_count, ctx.elements,
                                                   position, rec, M, &fs, w)) {
                        snprintf(why, 96, "no common field for %s", names->symbol[q]);
                        free(w);
                        return -1;
                    }
                    cy trace, ours, mine, cj;
                    cy_zero(fm, &trace);
                    for (int i = 0; i < ctx.dim; i++) cy_add(fm, &trace, &w[i * ctx.dim + i], &trace);
                    if (chi_exact(c, f, sigma, &x, &ours)) {
                        free(w);
                        return -1;
                    }
                    if (nu != 1) {
                        cy_conj(f, &ours, &cj);
                        cy_add(f, &ours, &cj, &ours);
                    }
                    cy_lift(f, &ours, fm, &mine);
                    cy_sub(fm, &mine, &trace, &mine);
                    ok = cy_is_zero(fm, &mine);
                }
        free(w);
        if (ok) {
            matches++;
            found = q;
        } else
            alive[q] = 0;
    }
    if (matches <= 1) break;
    }
    if (matches != 1) {
        snprintf(why, 96, "%d candidates match", matches);
        return -1;
    }
    return found;
}

/* site >= 0: the Bloch displacements (axial: moments, times det R) of one site; site < 0: the strain (k = 0) */
#include "tables.h"

int mag_isotropy_row_dim(int64_t row_id);
int mag_isotropy_row_free(int64_t row_id);
int parent_operation_records(int parent_sg, int child_sg, const int64_t *basis, const int64_t *origin, int64_t (*out)[5]);

/* a Source record (pml frame) as an affine map of G in the CIF frame */
static int record_to_aff(const int64_t *rec, sp_aff *a)
{
    rat m[3][3], origin[3], minv[3][3];
    if (geo_setting_to_cinter(SP_SG, 1, SP_SETTING, m, origin) || rat_inverse3((const rat (*)[3])m, minv))
        return spfail("parent setting tables are inconsistent");
    rat point[3][3], tmp[3][3], rt[3][3], tp[3], tc[3], o1[3];
    for (int i = 0; i < 3; i++) {
        rat unit[3] = {rat_int(i == 0), rat_int(i == 1), rat_int(i == 2)};
        geo_vrot(SP_SG, (int)rec[4], unit, point[i]);
    }
    for (int i = 0; i < 3; i++) rat_row_mul(minv[i], (const rat (*)[3])point, tmp[i]);
    for (int i = 0; i < 3; i++) rat_row_mul(tmp[i], (const rat (*)[3])m, rt[i]);
    for (int z = 0; z < 3; z++) tp[z] = rat_make(rec[z], rec[3]);
    rat_row_mul(tp, (const rat (*)[3])m, tc);
    rat_row_mul(origin, (const rat (*)[3])rt, o1);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            if (rt[j][i].d != 1) return spfail("record rotation is not integral in the CIF frame");
            a->r[i][j] = (int)rt[j][i].n;
        }
    for (int z = 0; z < 3; z++) {
        rat v = rat_mul(rat_sub(rat_add(tc[z], origin[z]), o1[z]), rat_int(24));
        if (v.d != 1) return spfail("record translation off the 1/24 grid");
        a->t[z] = (int)v.n;
    }
    return rat_overflow ? spfail("rational overflow in the naming bridge") : 0;
}

/* intersect canonical rows x (r of them) with the fixed space of mat; returns the new rank */
static int meet_fixed(uint64_t *x, int r, const uint64_t *mat, int D)
{
    uint64_t *img = malloc(sizeof(uint64_t) * (size_t)r * D + 1), *tr = malloc(sizeof(uint64_t) * (size_t)D * r + 1);
    uint64_t *cs = malloc(sizeof(uint64_t) * (size_t)r * r + 1), *rows = calloc((size_t)r * D + 1, sizeof(uint64_t));
    prow_image(x, r, mat, D, img);
    for (int i = 0; i < r; i++)
        for (int a = 0; a < D; a++) tr[a * r + i] = psub(img[i * D + a], x[i * D + a]);
    int cnt = pnull(tr, D, r, cs);
    for (int v2 = 0; v2 < cnt; v2++)
        for (int i = 0; i < r; i++)
            if (cs[v2 * r + i])
                for (int a = 0; a < D; a++) rows[v2 * D + a] = padd(rows[v2 * D + a], pmul(cs[v2 * r + i], x[i * D + a]));
    int rr = cnt ? pcanon(rows, cnt, D) : 0;
    memcpy(x, rows, sizeof(uint64_t) * (size_t)rr * D);
    free(img), free(tr), free(cs), free(rows);
    return rr;
}

/* the number of data_isotropy rows of a tabulated irrep (0 when untabulated) */
static int table_rows(int gid, int magnetic)
{
    section_view old, pointer;
    if (section_find("little.little_irr_old", &old) ||
        section_find(magnetic ? "magnetic.mag_iso_irrep_pointer" : "isotropy.isotropy_irrep_pointer", &pointer))
        return 0;
    int64_t old_id = section_int(&old, gid - 1);
    if (old_id <= 0) return 0;
    return (int)(section_int(&pointer, old_id) - section_int(&pointer, old_id - 1));
}

/* labels from the data_isotropy rows of a tabulated irrep: each row's subgroup (number, basis, origin) acts on
 * Gamma, its fixed space is an isotropy subspace, and its orbit is one of ours */
static int table_labels(const sp_induced *c, int gid, const pset *seen, const int *orbit_of, char (*labels)[16], size_t nreps)
{
    static section_view old, pointer, basis, origin, label;
    static int ready;
    if (!ready) {
        if (section_find("little.little_irr_old", &old) || section_find("isotropy.isotropy_irrep_pointer", &pointer) ||
            section_find("isotropy.isotropy_basis", &basis) ||
            section_find("isotropy.isotropy_origin", &origin) || section_find("isotropy.isotropy_orderparam_label", &label))
            return spfail("isotropy tables are missing");
        ready = 1;
    }
    int64_t old_id = section_int(&old, gid - 1);
    if (old_id <= 0) return 0;
    int64_t start = section_int(&pointer, old_id - 1), end = section_int(&pointer, old_id);
    int D = c->D;
    uint64_t *x = malloc(sizeof(uint64_t) * ((size_t)D * D + 1)), *gm = malloc(sizeof(uint64_t) * (size_t)D * D), *key = malloc(sizeof(uint64_t) * (1 + (size_t)D * D));
    for (int64_t row = start; row < end; row++) {
        int64_t idx = row - 1, b9[9], o4[4], records[GEO_MAX_OPS][5];
        for (int i = 0; i < 9; i++) b9[i] = section_int(&basis, idx * 9 + i);
        for (int i = 0; i < 4; i++) o4[i] = section_int(&origin, idx * 4 + i);
        int n = parent_operation_records(SP_SG, isotropy_row_subgroup(idx + 1), b9, o4, records);
        if (n <= 0) return spfail("table subgroup does not map into the parent");
        memset(x, 0, sizeof(uint64_t) * (size_t)D * D);
        for (int i = 0; i < D; i++) x[i * D + i] = 1;
        int r = D;
        for (int q = 0; q < n + 3 && r; q++) {
            sp_aff a;
            if (q < n) {
                if (record_to_aff(records[q], &a)) return 1;
            } else {
                /* a lattice vector of the subgroup: a pml basis row */
                int64_t rec[5] = {b9[3 * (q - n)], b9[3 * (q - n) + 1], b9[3 * (q - n) + 2], 1, 0};
                int64_t src[GEO_MAX_OPS][5];
                int nsrc = geo_sg_records(SP_SG, src);
                for (int z = 0; z < nsrc; z++) {
                    int ident = 1;
                    rat unit[3][3];
                    for (int i = 0; i < 3; i++) {
                        rat e[3] = {rat_int(i == 0), rat_int(i == 1), rat_int(i == 2)};
                        geo_vrot(SP_SG, (int)src[z][4], e, unit[i]);
                        for (int j = 0; j < 3; j++) ident &= rat_eq(unit[i][j], rat_int(i == j));
                    }
                    if (ident) {
                        rec[4] = src[z][4];
                        break;
                    }
                }
                if (record_to_aff(rec, &a)) return 1;
            }
            if (gamma_of(c, &a, gm)) return 1;
            r = meet_fixed(x, r, gm, D);
        }
        if (!r) continue;
        memset(key, 0, sizeof(uint64_t) * (1 + (size_t)D * D));
        key[0] = (uint64_t)r;
        memcpy(key + 1, x, sizeof(uint64_t) * (size_t)r * D);
        size_t slot = phash(key, seen->width) & (seen->slots - 1);
        int64_t found = -1;
        for (; seen->table[slot] >= 0; slot = (slot + 1) & (seen->slots - 1))
            if (!memcmp(seen->items + (size_t)seen->table[slot] * seen->width, key, sizeof(uint64_t) * (size_t)seen->width)) {
                found = seen->table[slot];
                break;
            }
        if (found < 0) continue;
        int rep = orbit_of[found];
        const char *text;
        int64_t length;
        if (rep < 0 || (size_t)rep >= nreps || section_string(&label, idx, &text, &length)) continue;
        while (length && (*text == ' ')) text++, length--;
        while (length && text[length - 1] == ' ') length--;
        if (length > 15) length = 15;
        memcpy(labels[rep], text, (size_t)length);
        labels[rep][length] = 0;
    }
    free(x), free(gm), free(key);
    return 0;
}

#include "kernels.h"

/* the walk's candidate order (the order rule of the dynamic labels): Source records in table order, each with the
 * translation classes modulo the irrep's kernel lattice in op_new_fractionals order; distinct Gamma matrices kept
 * in first-seen order */
int magnetic_records_iso(int group, int capacity, int64_t *records_out);
int mag_point_action(int mode, int op, int *sign);
int mag_orderparam_group(int sg);
int mag_tables(void);

/* the walk's elements in candidate order: Source records in table order, each with the translation classes modulo the
 * irrep's kernel lattice in op_new_fractionals order; records keep the spatial point operation, signs the time reversal */
static int walk_records(int gid, const int64_t *kparam, int magnetic, int64_t (**out)[5], int **signs, int *count)
{
    irrep_context ctx;
    if (irrep_context_for(gid, 4, kparam, &ctx)) return spfail("irrep context for the walk order");
    if (magnetic && mag_tables()) return spfail("magnetic tables are missing");
    int64_t records[192][5];
    int nrec = magnetic ? magnetic_records_iso(mag_orderparam_group(SP_SG), 192, &records[0][0]) : geo_sg_records(SP_SG, records);
    int cops[192], csigns[192];
    if (nrec <= 0) return spfail("records for the walk order");
    for (int i = 0; i < nrec; i++) cops[i] = mag_point_action(magnetic, (int)records[i][4], &csigns[i]);
    double *full = calloc((size_t)ctx.dim * 48 + 1, sizeof(double));
    int64_t basis[9], kept[5 * 192 + 5];
    int kept_count = 0;
    for (int row = 0; row < ctx.dim; row++) full[row * 48 + row] = 1.0;
    int status = subgroup_rows_in_context(ctx.dim, ctx.arms, ctx.k, ctx.subtract, ctx.element_count, ctx.elements, ctx.sparse, ctx.available,
                                          ctx.identity_op, full, (int64_t)ctx.dim * 48, ctx.dim, nrec, &records[0][0], cops, csigns, 256, basis,
                                          kept, &kept_count);
    free(full);
    if (status) return spfail("kernel lattice for the walk order");
    i128 det = (i128)basis[0] * (basis[4] * basis[8] - basis[5] * basis[7]) - (i128)basis[1] * (basis[3] * basis[8] - basis[5] * basis[6]) +
               (i128)basis[2] * (basis[3] * basis[7] - basis[4] * basis[6]);
    if (det < 0) det = -det;
    if (!det || det > (1 << 20)) return spfail("kernel lattice out of range");
    int64_t *fractions = malloc(sizeof(int64_t) * 4 * (size_t)(det + 1));
    int nf = 0;
    if (op_new_fractionals(basis, fractions, (int)det + 1, &nf)) {
        free(fractions);
        return spfail("translation classes for the walk order");
    }
    *out = malloc(sizeof(int64_t[5]) * (size_t)nrec * nf + 1);
    *signs = malloc(sizeof(int) * (size_t)nrec * nf + 1);
    *count = 0;
    for (int sp = 0; sp < nrec; sp++)
        for (int f = 0; f < nf; f++) {
            int64_t *op = (*out)[*count];
            if (op_vadd(fractions + 4 * f, records[sp], op)) return spfail("walk-order element");
            op[4] = cops[sp];
            (*signs)[(*count)++] = csigns[sp];
        }
    free(fractions);
    return 0;
}

/* the walk's candidate order (the order rule of the dynamic labels): distinct Gamma matrices in first-seen order */
static int walk_order(const sp_induced *c, int gid, const int64_t *kparam, int magnetic, pset *out)
{
    int64_t (*records)[5];
    int *signs, count;
    if (walk_records(gid, kparam, magnetic, &records, &signs, &count)) return 1;
    uint64_t *gm = malloc(sizeof(uint64_t) * (size_t)c->D * c->D);
    int added;
    for (int q = 0; q < count; q++) {
        sp_aff a;
        if (record_to_aff(records[q], &a) || gamma_of(c, &a, gm)) {
            free(records), free(signs), free(gm);
            return spfail("walk-order element");
        }
        if (signs[q] < 0)
            for (int i = 0; i < c->D * c->D; i++) gm[i] = psub(0, gm[i]);
        if (pset_put(out, gm, &added) < 0) return spfail("out of memory");
    }
    free(records), free(signs), free(gm);
    return 0;
}

/* invert a d x d matrix mod p; nonzero when singular */
static int pinvert(const uint64_t *a, int d, uint64_t *out)
{
    uint64_t *m = malloc(sizeof(uint64_t) * (size_t)d * 2 * d);
    int *piv = malloc(sizeof(int) * (size_t)(d + 1));
    for (int i = 0; i < d; i++)
        for (int j = 0; j < 2 * d; j++) m[i * 2 * d + j] = j < d ? a[i * d + j] : (uint64_t)(j - d == i);
    int rank = prref(m, d, 2 * d, piv);
    int ok = rank == d;
    for (int i = 0; i < d && ok; i++) ok = piv[i] == i;
    for (int i = 0; i < d && ok; i++)
        for (int j = 0; j < d; j++) out[i * d + j] = m[i * 2 * d + d + j];
    free(m), free(piv);
    return !ok;
}

static void pmatmul(const uint64_t *a, const uint64_t *b2, int n, int k, int m2, uint64_t *out)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m2; j++) {
            uint64_t t = 0;
            for (int l = 0; l < k; l++)
                if (a[i * k + l] && b2[l * m2 + j]) t = padd(t, pmul(a[i * k + l], b2[l * m2 + j]));
            out[i * m2 + j] = t;
        }
}

/* the BYU gauge from Gamma: Source's W on the generators fixes an intertwiner T (W T = T Gamma_big, with Gamma_big =
 * Gamma, or Gamma + its dual for a complex or quaternionic type); then T Gamma_big(x) T^-1 is compared with Source's W(x)
 * on every walk element. counts: [checked, differ, status] */
static int gauge_check(const sp_induced *c, const cy_field *f, uint64_t w, int gid, int nu, int magnetic, const sp_aff *gens, int ngens, int *counts)
{
    counts[0] = counts[1] = counts[2] = 0;
    irrep_context ctx;
    if (irrep_context_for(gid, 4, SP_NAMES->kparam, &ctx)) return (counts[2] = 1), 0;
    int D = c->D, Db = nu == 1 ? D : 2 * D, rd = ctx.dim;
    if (rd != Db) return (counts[2] = 2), 0;
    if ((size_t)rd * Db > 2304) return (counts[2] = 3), 0;
    cy *wc = malloc(sizeof(cy) * (size_t)rd * rd);
    uint64_t *gm = malloc(sizeof(uint64_t) * (size_t)D * D), *gi = malloc(sizeof(uint64_t) * (size_t)D * D);
    uint64_t *big = calloc((size_t)Db * Db, sizeof(uint64_t)), *wm = malloc(sizeof(uint64_t) * (size_t)rd * rd);
    /* Source W(x) mod p in this field; transposed when Source's matrices act on rows */
    int transposed = -1;
#define SOURCE_W(aff, out)                                                                                                       \
    do {                                                                                                                       \
        int64_t r64[5];                                                                                                        \
        int rec[5];                                                                                                            \
        if (aff_to_record(c->g, (aff), r64)) return 1;                                                                         \
        for (int z_ = 0; z_ < 5; z_++) rec[z_] = (int)r64[z_];                                                                 \
        int pos_ = bw_find_position(ctx.element_count, ctx.elements, rec);                                                     \
        const cy_field *fs_;                                                                                                   \
        if (!pos_ || exact_record_matrix(gid, ctx.dim, ctx.arms, ctx.k, ctx.subtract, ctx.element_count, ctx.elements, pos_, rec, \
                                         f->m, &fs_, wc))                                                                      \
            return (counts[2] = 4), 0;                                                                                         \
        for (int i_ = 0; i_ < rd; i_++)                                                                                        \
            for (int j_ = 0; j_ < rd; j_++) (out)[i_ * rd + j_] = cy_mod(f, &wc[transposed == 1 ? j_ * rd + i_ : i_ * rd + j_], w); \
    } while (0)
#define BIG(aff)                                                                                                                 \
    do {                                                                                                                       \
        if (gamma_of(c, (aff), gm)) return 1;                                                                                  \
        memset(big, 0, sizeof(uint64_t) * (size_t)Db * Db);                                                                    \
        for (int i_ = 0; i_ < D; i_++)                                                                                         \
            for (int j_ = 0; j_ < D; j_++) big[i_ * Db + j_] = gm[i_ * D + j_];                                                \
        if (Db > D) {                                                                                                          \
            sp_aff inv_;                                                                                                       \
            aff_inv(c->g, (aff), &inv_);                                                                                       \
            if (gamma_of(c, &inv_, gi)) return 1;                                                                              \
            for (int i_ = 0; i_ < D; i_++)                                                                                     \
                for (int j_ = 0; j_ < D; j_++) big[(D + i_) * Db + D + j_] = gi[j_ * D + i_];                                  \
        }                                                                                                                      \
    } while (0)
    /* orientation: W(a) W(b) = W(ab) for two generators */
    {
        uint64_t *w1 = malloc(sizeof(uint64_t) * (size_t)rd * rd), *w2 = malloc(sizeof(uint64_t) * (size_t)rd * rd), *w12 = malloc(sizeof(uint64_t) * (size_t)rd * rd);
        uint64_t *pr = malloc(sizeof(uint64_t) * (size_t)rd * rd);
        sp_aff ab;
        int a0 = 0, b0 = ngens > 1 ? 1 : 0;
        aff_mul(&gens[a0], &gens[b0], &ab);
        transposed = 0;
        SOURCE_W(&gens[a0], w1);
        SOURCE_W(&gens[b0], w2);
        SOURCE_W(&ab, w12);
        pmatmul(w1, w2, rd, rd, rd, pr);
        if (memcmp(pr, w12, sizeof(uint64_t) * (size_t)rd * rd)) transposed = 1;
        free(w1), free(w2), free(w12), free(pr);
    }
    /* T: rd x Db unknowns, W(g) T - T Big(g) = 0 for each generator */
    int unknowns = rd * Db, rows = 0;
    uint64_t *sys = calloc((size_t)ngens * unknowns * unknowns + 1, sizeof(uint64_t));
    for (int q = 0; q < ngens; q++) {
        SOURCE_W(&gens[q], wm);
        BIG(&gens[q]);
        for (int i = 0; i < rd; i++)
            for (int j = 0; j < Db; j++, rows++) {
                uint64_t *row = sys + (size_t)rows * unknowns;
                for (int l = 0; l < rd; l++) row[l * Db + j] = padd(row[l * Db + j], wm[i * rd + l]);
                for (int l = 0; l < Db; l++) row[i * Db + l] = psub(row[i * Db + l], big[l * Db + j]);
            }
    }
    uint64_t *null = malloc(sizeof(uint64_t) * (size_t)unknowns * unknowns + 1);
    int cnt = pnull(sys, rows, unknowns, null);
    free(sys);
    if (!cnt) return (counts[2] = 5), 0;
    uint64_t *t = calloc((size_t)unknowns, sizeof(uint64_t)), *tinv = malloc(sizeof(uint64_t) * (size_t)unknowns), seed = 0x9e3779b97f4a7c15ULL;
    for (int v = 0; v < cnt; v++) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        uint64_t coef = (seed >> 11) % PR;
        for (int i = 0; i < unknowns; i++) t[i] = padd(t[i], pmul(coef, null[(size_t)v * unknowns + i]));
    }
    free(null);
    if (pinvert(t, rd, tinv)) return (counts[2] = 6), 0;
    /* every walk element */
    int64_t (*records)[5];
    int *signs, count;
    if (walk_records(gid, SP_NAMES->kparam, magnetic, &records, &signs, &count)) return 1;
    uint64_t *tmp = malloc(sizeof(uint64_t) * (size_t)rd * Db), *ours = malloc(sizeof(uint64_t) * (size_t)rd * rd);
    for (int q = 0; q < count; q++) {
        sp_aff a;
        if (record_to_aff(records[q], &a)) return 1;
        SOURCE_W(&a, wm);
        BIG(&a);
        pmatmul(t, big, rd, Db, Db, tmp);
        pmatmul(tmp, tinv, rd, Db, rd, ours);
        counts[0]++;
        if (memcmp(ours, wm, sizeof(uint64_t) * (size_t)rd * rd)) counts[1]++;
    }
#undef SOURCE_W
#undef BIG
    free(records), free(signs), free(tmp), free(ours), free(t), free(tinv), free(wc), free(gm), free(gi), free(big), free(wm);
    return 0;
}

int orderparam_to_subgroup_rows(int gid, int param_count, const int64_t *kparam, const double *values, int64_t value_count, int row_count,
                                int candidate_count, const int64_t *candidates, const int *candidate_matrix_ops, const int *candidate_signs,
                                int max_iter, int64_t *basis_out, int64_t *kept_out, int *kept_count);

/* labels from the magnetic isotropy rows of a tabulated irrep: each row's direction gives its magnetic subgroup (records
 * with time reversal, and a lattice) as the table path derives it; its fixed space on Gamma names one of our orbits */
static int mag_table_labels(const sp_induced *c, int gid, const pset *seen, const int *orbit_of, char (*labels)[16], size_t nreps)
{
    static section_view old, pointer, opptr, ops, label, consts;
    static int ready;
    if (!ready) {
        if (mag_tables() || section_find("little.little_irr_old", &old) || section_find("magnetic.mag_iso_irrep_pointer", &pointer) ||
            section_find("magnetic.mag_iso_orderparam_pointer", &opptr) || section_find("magnetic.mag_iso_orderparam", &ops) ||
            section_find("magnetic.mag_iso_orderparam_label", &label) || section_find("const.values", &consts))
            return spfail("magnetic isotropy tables are missing");
        ready = 1;
    }
    int64_t old_id = section_int(&old, gid - 1);
    if (old_id <= 0) return 0;
    int64_t start = section_int(&pointer, old_id - 1), end = section_int(&pointer, old_id);
    const double *constants = (const double *)consts.data;
    int64_t cand[192][5];
    int ncand = magnetic_records_iso(mag_orderparam_group(SP_SG), 192, &cand[0][0]), cops[192], csigns[192];
    if (ncand <= 0) return spfail("magnetic records");
    for (int i = 0; i < ncand; i++) cops[i] = mag_point_action(1, (int)cand[i][4], &csigns[i]);
    int D = c->D;
    uint64_t *x = malloc(sizeof(uint64_t) * ((size_t)D * D + 1)), *gm = malloc(sizeof(uint64_t) * (size_t)D * D), *key = malloc(sizeof(uint64_t) * (1 + (size_t)D * D));
    for (int64_t row = start; row < end; row++) {
        int64_t idx = row - 1, dim = mag_isotropy_row_dim(row), fn = mag_isotropy_row_free(row), ptr = section_int(&opptr, idx) - 1;
        if (fn <= 0 || dim <= 0 || dim > 48) continue;
        double *values = calloc((size_t)fn * 48 + 1, sizeof(double));
        for (int64_t p = 0; p < fn; p++)
            for (int64_t q = 0; q < dim; q++) values[p * 48 + q] = constants[section_int(&ops, ptr + p * dim + q)];
        int64_t sub_basis[9], kept[192][5], none = 0;
        int kept_count = 0;
        int code = orderparam_to_subgroup_rows(gid, 0, &none, values, fn * 48, (int)fn, ncand, &cand[0][0], cops, csigns, 256, sub_basis,
                                               &kept[0][0], &kept_count);
        free(values);
        if (code) continue;
        memset(x, 0, sizeof(uint64_t) * (size_t)D * D);
        for (int i = 0; i < D; i++) x[i * D + i] = 1;
        int r = D;
        for (int q = 0; q < kept_count + 3 && r; q++) {
            sp_aff a;
            int sign = 1;
            int64_t rec[5];
            if (q < kept_count) {
                memcpy(rec, kept[q], sizeof rec);
                rec[4] = mag_point_action(1, (int)kept[q][4], &sign);
            } else {
                int64_t src[GEO_MAX_OPS][5];
                int nsrc = geo_sg_records(SP_SG, src);
                rec[0] = sub_basis[3 * (q - kept_count)];
                rec[1] = sub_basis[3 * (q - kept_count) + 1];
                rec[2] = sub_basis[3 * (q - kept_count) + 2];
                rec[3] = 1;
                rec[4] = 0;
                for (int z = 0; z < nsrc && !rec[4]; z++) {
                    int ident = 1;
                    for (int i = 0; i < 3; i++) {
                        rat e[3] = {rat_int(i == 0), rat_int(i == 1), rat_int(i == 2)}, unit[3];
                        geo_vrot(SP_SG, (int)src[z][4], e, unit);
                        for (int j = 0; j < 3; j++) ident &= rat_eq(unit[j], rat_int(i == j));
                    }
                    if (ident) rec[4] = src[z][4];
                }
            }
            if (record_to_aff(rec, &a) || gamma_of(c, &a, gm)) return 1;
            if (sign < 0)
                for (int i = 0; i < D * D; i++) gm[i] = psub(0, gm[i]);
            r = meet_fixed(x, r, gm, D);
        }
        if (!r) continue;
        memset(key, 0, sizeof(uint64_t) * (1 + (size_t)D * D));
        key[0] = (uint64_t)r;
        memcpy(key + 1, x, sizeof(uint64_t) * (size_t)r * D);
        size_t slot = phash(key, seen->width) & (seen->slots - 1);
        int64_t found = -1;
        for (; seen->table[slot] >= 0; slot = (slot + 1) & (seen->slots - 1))
            if (!memcmp(seen->items + (size_t)seen->table[slot] * seen->width, key, sizeof(uint64_t) * (size_t)seen->width)) {
                found = seen->table[slot];
                break;
            }
        if (found < 0) continue;
        int rep = orbit_of[found];
        const char *text;
        int64_t length;
        if (rep < 0 || (size_t)rep >= nreps || section_string(&label, idx, &text, &length)) continue;
        while (length && (*text == ' ')) text++, length--;
        while (length && text[length - 1] == ' ') length--;
        if (length > 15) length = 15;
        memcpy(labels[rep], text, (size_t)length);
        labels[rep][length] = 0;
    }
    free(x), free(gm), free(key);
    return 0;
}


static int SP_NAMED = -1;
static char SP_WHY[96];

static int k_opds_in(const sp_group *g, const sp_little *L, const rat *k, const cy_field *f, uint64_t w, const cy *sigma, int dim, int site,
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
    sp_induced ind = {g, L, k, s, dim, D, (const rat (*)[3])arms, raff, rinv, smat};
    int named_gid = 0;
    if (SP_NAMES && SP_NAMES->count) {
        char why[96];
        int q = SP_NAMED;
        snprintf(why, sizeof why, "%s", SP_WHY);
        jb_lit(b, ", \"symbol\": ");
        if (q >= 0) named_gid = SP_NAMES->gid[q];
        if (q >= 0) {
            char sym[32];
            snprintf(sym, sizeof sym, "%s%s", magnetic ? "m" : "", SP_NAMES->symbol[q]);
            jb_cstring(b, sym);
        } else {
            jb_lit(b, "null, \"symbol_error\": ");
            jb_cstring(b, why);
        }
    }
    if (named_gid > 0 && getenv("APO_SPEC_GAUGE")) {
        int counts[3];
        if (gauge_check(&ind, f, w, named_gid, nu, magnetic, gens, ngens, counts)) return MDE_VALUE;
        char gt[96];
        snprintf(gt, sizeof gt, ", \"gauge\": [%d, %d, %d]", counts[0], counts[1], counts[2]);
        jb_lit(b, gt);
    }
    uint64_t *gmat = calloc((size_t)(ngens + 1) * D * D, sizeof(uint64_t));
    for (int q = 0; q < ngens; q++)
        if (gamma_of(&ind, &gens[q], gmat + (size_t)q * D * D)) return MDE_VALUE;
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
    size_t ocap = 1024;
    int *orbit_of = malloc(sizeof(int) * ocap);
    int64_t at0 = pset_put(&seen, key, &added);
    orbit_of[at0] = 0;
#define MARK(index, rep)                                                   \
    do {                                                                   \
        if ((size_t)(index) >= ocap) {                                     \
            while ((size_t)(index) >= ocap) ocap *= 2;                     \
            orbit_of = realloc(orbit_of, sizeof(int) * ocap);              \
        }                                                                  \
        orbit_of[index] = (rep);                                           \
    } while (0)
    /* candidate order: the walk's for an untabulated ordinary irrep (its labels count in that order), else the image's */
    pset cand = {D * D, 0, 0, 0, 0, 0};
    const pset *candidates = &image;
    if (named_gid > 0 && !table_rows(named_gid, magnetic)) {
        if (walk_order(&ind, named_gid, SP_NAMES->kparam, magnetic, &cand)) return MDE_VALUE;
        if (cand.count != order) return spfail("walk order and image differ"), MDE_VALUE;
        candidates = &cand;
    }
    for (size_t head = 0; head < nreps; head++) {
        int r = (int)reps[head * (1 + D * D)];
        const uint64_t *x = reps + head * (1 + D * D) + 1;
        for (size_t q = 0; q < candidates->count; q++) {
            const uint64_t *mat = candidates->items + q * candidates->width;
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
            int64_t at = pset_put(&seen, key, &added);
            if (!added) continue;
            MARK(at, (int)nreps);
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
                int64_t at2 = pset_put(&seen, orbit, &added);
                if (added) MARK(at2, (int)nreps - 1);
            }
            free(orbit);
        }
    }
#undef MARK
    /* names of the OPDs: the table rows of the named irrep, each located by the fixed space of its subgroup */
    char (*labels)[16] = calloc(nreps, sizeof *labels);
    if (named_gid > 0 && (magnetic ? mag_table_labels(&ind, named_gid, &seen, orbit_of, labels, nreps)
                                   : table_labels(&ind, named_gid, &seen, orbit_of, labels, nreps)))
        return MDE_VALUE;
    {
        /* untabulated irreps: P, C, S or nD by real dimension, numbered in discovery order */
        int any = 0, counts[100] = {0};
        for (size_t q = 0; q < nreps; q++) any |= labels[q][0] != 0;
        for (size_t q = 0; q < nreps && !any; q++) {
            int rd = (int)reps[q * (1 + D * D)] * (nu == 1 ? 1 : 2);
            if (rd < 1 || rd > 99) continue;
            int i = ++counts[rd];
            if (rd == 1) snprintf(labels[q], 16, "P%d", i);
            else if (rd == 2) snprintf(labels[q], 16, "C%d", i);
            else if (rd == 3) snprintf(labels[q], 16, "S%d", i);
            else snprintf(labels[q], 16, "%dD%d", rd, i);
        }
    }
    jb_lit(b, ", \"image\": ");
    char t[160];
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
        int64_t ident14[14] = {0};
        int found = identify_opd(&ind, x, r, magnetic, ident14);
        if (found < 0 || fail_status) return MDE_VALUE;
        snprintf(t, sizeof t, "%s{\"dim\": %d, \"index\": %zu, \"subgroup\": %lld, \"label\": ", q ? ", " : "", r * (nu == 1 ? 1 : 2), order / stab,
                 found > 0 ? (long long)ident14[0] : -1LL);
        jb_lit(b, t);
        if (labels[q][0]) jb_cstring(b, labels[q]);
        else jb_lit(b, "null");
        jb_lit(b, "}");
    }
    jb_lit(b, "]");
    free(labels);
    free(orbit_of);
    pset_free(&cand);
    pset_free(&image);
    pset_free(&seen);
    free(key), free(rows), free(img), free(diff), free(reps), free(ident), free(prod), free(gmat), free(raff), free(rinv), free(rep), free(arms);
    free(smat), free(sig), free(local), free(global), free(phase), free(target), free(v0), free(v), free(u), free(e0), free(basis), free(piv), free(sdense);
    return 0;
}

/* names the irrep first (exactly), then runs the OPDs modulo a prime of a field that also holds Source's matrices */
static int k_opds(const sp_group *g, const sp_little *L, const rat *k, const cy_field *f, uint64_t w, const cy *sigma, int dim, int site,
                  int axial, int magnetic, int nu, jbuf *b)
{
    SP_NAMED = -1;
    SP_WHY[0] = 0;
    if (!SP_NAMES || !SP_NAMES->count) return k_opds_in(g, L, k, f, w, sigma, dim, site, axial, magnetic, nu, b);
    int s = L->star, *rep = malloc(sizeof(int) * (size_t)s), narms = 0;
    rat (*arms)[3] = malloc(sizeof(rat[3]) * (size_t)s);
    for (int e = 0; e < g->n && narms < s; e++) {
        rat a[3];
        k_times(k, (const int (*)[3])g->ops[g->inv[e]].rot, a);
        int known = 0;
        for (int q = 0; q < narms && !known; q++) known = k_same(g, arms[q], a);
        if (!known) {
            memcpy(arms[narms], a, sizeof a);
            rep[narms++] = e;
        }
    }
    sp_aff *raff = malloc(sizeof(sp_aff) * (size_t)s), *rinv = malloc(sizeof(sp_aff) * (size_t)s);
    for (int i = 0; i < s; i++) {
        memcpy(raff[i].r, g->ops[rep[i]].rot, sizeof raff[i].r);
        memcpy(raff[i].t, g->ops[rep[i]].tran, sizeof raff[i].t);
        aff_inv(g, &raff[i], &rinv[i]);
    }
    int D = dim * s;
    sp_induced pre = {g, L, k, s, dim, D, (const rat (*)[3])arms, raff, rinv, 0};
    SP_NAMED = name_irrep(&pre, f, sigma, nu, nu == 1 ? D : 2 * D, SP_NAMES, SP_WHY);
    free(rep), free(arms), free(raff), free(rinv);
    if (fail_status) return MDE_VALUE;
    if (SP_NAMED < 0) return k_opds_in(g, L, k, f, w, sigma, dim, site, axial, magnetic, nu, b);
    /* the field of Source's matrices on the walk's elements */
    int gid = SP_NAMES->gid[SP_NAMED], M = f->m;
    irrep_context ctx;
    int64_t (*records)[5];
    int *signs, count;
    if (irrep_context_for(gid, 4, SP_NAMES->kparam, &ctx) || walk_records(gid, SP_NAMES->kparam, 0, &records, &signs, &count)) return MDE_VALUE;
    for (int q = 0; q < count; q++) {
        int rec[5];
        for (int z = 0; z < 5; z++) rec[z] = (int)records[q][z];
        int pos = bw_find_position(ctx.element_count, ctx.elements, rec);
        const cy_field *fs;
        int need = pos ? -exact_record_matrix(gid, ctx.dim, ctx.arms, ctx.k, ctx.subtract, ctx.element_count, ctx.elements, pos, rec, -1, &fs, 0) : 0;
        if (need <= 0) continue;
        int a = M, b2 = need;
        while (b2) {
            int t = a % b2;
            a = b2;
            b2 = t;
        }
        M = M / a * need;
    }
    free(records), free(signs);
    if (M == f->m) return k_opds_in(g, L, k, f, w, sigma, dim, site, axial, magnetic, nu, b);
    const cy_field *fm = cy_field_for(M);
    uint64_t saved = PR, wm;
    const int *units;
    const uint64_t *vinv;
    if (!fm || modp_field(M, fm->phi, 0, &PR, &wm, &units, &vinv)) {
        PR = saved;
        return spfail("no prime for the gauge field"), MDE_VALUE;
    }
    wm = powmod(wm, (uint64_t)units[0], PR);
    int kc = L->h.ncls;
    cy *lifted = malloc(sizeof(cy) * (size_t)kc);
    for (int l = 0; l < kc; l++) cy_lift(f, &sigma[l], fm, &lifted[l]);
    int status = k_opds_in(g, L, k, fm, wm, lifted, dim, site, axial, magnetic, nu, b);
    free(lifted);
    PR = saved;
    return status;
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

int spec_gamma(const jvalue *cif_info, const int *sites, const int *kinds, int nsites, int strain, const rat *kvec, int sg, int setting,
               const void *names, jbuf *b)
{
    fail_status = 0;
    SP_NAMES = names;
    SP_SG = sg;
    SP_SETTING = setting;
    if (geo_load()) return spfail("Source tables are inconsistent"), MDE_VALUE;
    sp_group g;
    memset(&g, 0, sizeof g);
    if (build_group(cif_info, sites, nsites, &g)) return MDE_VALUE;
    return spec_k(&g, kvec, nsites, kinds, strain, b);
}

/* the naming candidates of the selected k point: its catalog irreps and the k parameters */
void *spec_names_new(int count, const int *gids, const char *const *symbols, const int64_t *kparam)
{
    sp_names *n = calloc(1, sizeof *n);
    if (!n) return 0;
    n->count = count < 64 ? count : 64;
    for (int i = 0; i < n->count; i++) {
        n->gid[i] = gids[i];
        snprintf(n->symbol[i], sizeof n->symbol[i], "%s", symbols[i]);
    }
    memcpy(n->kparam, kparam, sizeof n->kparam);
    return n;
}
