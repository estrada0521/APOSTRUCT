#include "cyclo.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

int cy_overflow;

static i128 ck_mul(i128 a, i128 b)
{
    i128 out;
    if (__builtin_mul_overflow(a, b, &out)) {
        cy_overflow = 1;
        return 0;
    }
    return out;
}

static i128 ck_add(i128 a, i128 b)
{
    i128 out;
    if (__builtin_add_overflow(a, b, &out)) {
        cy_overflow = 1;
        return 0;
    }
    return out;
}

static i128 ck_gcd(i128 a, i128 b)
{
    i128 g = gcd_i128(a, b);
    if (g < 0) {
        cy_overflow = 1;
        return 1;
    }
    return g;
}

static int totient(int m)
{
    int out = m, n = m;
    for (int p = 2; p * p <= n; p++)
        if (n % p == 0) {
            while (n % p == 0) n /= p;
            out -= out / p;
        }
    if (n > 1) out -= out / n;
    return out;
}

static int mobius(int n)
{
    int out = 1;
    for (int p = 2; p * p <= n; p++)
        if (n % p == 0) {
            n /= p;
            if (n % p == 0) return 0;
            out = -out;
        }
    if (n > 1) out = -out;
    return out;
}

/* Phi_m = prod_{d | m} (x^d - 1)^{mu(m/d)}, multiplying the positive factors
 * and then dividing exactly by the negative ones. */
static int cyclotomic(int m, int64_t *poly, int phi)
{
    int size = 1;
    for (int d = 1; d <= m; d++)
        if (m % d == 0) size += d;
    int64_t *num = calloc((size_t)size, sizeof *num), *work = calloc((size_t)size, sizeof *work);
    if (!num || !work) {
        free(num);
        free(work);
        return 1;
    }
    num[0] = 1;
    int deg = 0;
    for (int d = 1; d <= m; d++) {
        if (m % d || mobius(m / d) != 1) continue;
        memset(work, 0, sizeof *work * (size_t)size);
        for (int i = 0; i <= deg; i++) {
            work[i] -= num[i];
            work[i + d] += num[i];
        }
        deg += d;
        memcpy(num, work, sizeof *work * (size_t)size);
    }
    for (int d = 1; d <= m; d++) {
        if (m % d || mobius(m / d) != -1) continue;
        /* divide by x^d - 1: q_i = num_{i+d} + q_{i+d} from the top */
        memset(work, 0, sizeof *work * (size_t)size);
        for (int i = deg - d; i >= 0; i--) work[i] = num[i + d] + (i + d <= deg - d ? work[i + d] : 0);
        deg -= d;
        memcpy(num, work, sizeof *work * (size_t)size);
    }
    int ok = deg == phi && num[phi] == 1;
    for (int i = 0; i <= phi && ok; i++) poly[i] = num[i];
    free(num);
    free(work);
    return !ok;
}

const cy_field *cy_field_for(int m)
{
    static cy_field *cache[64];
    static int count;
    for (int i = 0; i < count; i++)
        if (cache[i]->m == m) return cache[i];
    if (m < 1) return 0;
    int phi = totient(m);
    if (phi > CY_MAX_PHI) return 0;
    cy_field *f = calloc(1, sizeof *f);
    if (!f) return 0;
    f->m = m;
    f->phi = phi;
    if (cyclotomic(m, f->poly, phi)) {
        free(f);
        return 0;
    }
    if (count < 64) cache[count++] = f;
    return f;
}

void cy_copy(const cy_field *f, cy *out, const cy *a)
{
    if (out == a) return;
    memcpy(out->c, a->c, sizeof(i128) * (size_t)f->phi);
    out->den = a->den;
}

static void normalize(const cy_field *f, cy *a)
{
    if (a->den < 0) {
        a->den = -a->den;
        for (int i = 0; i < f->phi; i++) a->c[i] = -a->c[i];
    }
    i128 g = a->den;
    for (int i = 0; i < f->phi && g != 1; i++)
        if (a->c[i]) g = ck_gcd(g, a->c[i]);
    if (g > 1) {
        a->den /= g;
        for (int i = 0; i < f->phi; i++) a->c[i] /= g;
    }
    int zero = 1;
    for (int i = 0; i < f->phi && zero; i++) zero = !a->c[i];
    if (zero) a->den = 1;
}

void cy_zero(const cy_field *f, cy *out)
{
    memset(out->c, 0, sizeof(i128) * (size_t)f->phi);
    out->den = 1;
}

void cy_rat(const cy_field *f, i128 n, i128 d, cy *out)
{
    cy_zero(f, out);
    out->c[0] = n;
    out->den = d;
    normalize(f, out);
}

/* coefficients r[0..len) of a polynomial reduced in place modulo Phi_m */
static void reduce(const cy_field *f, i128 *r, int len)
{
    for (int k = len - 1; k >= f->phi; k--) {
        i128 t = r[k];
        if (!t) continue;
        r[k] = 0;
        for (int j = 0; j < f->phi; j++)
            if (f->poly[j]) r[k - f->phi + j] = ck_add(r[k - f->phi + j], -ck_mul(t, f->poly[j]));
    }
}

void cy_root(const cy_field *f, int64_t e, cy *out)
{
    int64_t m = f->m;
    e %= m;
    if (e < 0) e += m;
    cy_zero(f, out);
    if (!f->powers) {
        /* z^(e+1) = z * z^e: shift and reduce the one overflowing term */
        int64_t *table = calloc((size_t)m * f->phi, sizeof *table);
        if (!table) {
            cy_overflow = 1;
            return;
        }
        table[0] = 1;
        for (int64_t k = 1; k < m; k++) {
            const int64_t *prev = table + (k - 1) * f->phi;
            int64_t *next = table + k * f->phi, top = prev[f->phi - 1];
            for (int i = f->phi - 1; i > 0; i--) next[i] = prev[i - 1];
            next[0] = 0;
            for (int i = 0; i < f->phi; i++) next[i] -= top * f->poly[i];
        }
        ((cy_field *)f)->powers = table;
    }
    const int64_t *row = f->powers + e * f->phi;
    for (int i = 0; i < f->phi; i++) out->c[i] = row[i];
}

/* every coefficient and the denominator below 2^40: products and sums of a
 * few hundred of them stay far inside int128, so no overflow checks */
static int small(const cy_field *f, const cy *a)
{
    const i128 limit = (i128)1 << 40;
    if (a->den >= limit) return 0;
    for (int i = 0; i < f->phi; i++)
        if (a->c[i] >= limit || a->c[i] <= -limit) return 0;
    return 1;
}

void cy_add(const cy_field *f, const cy *a, const cy *b, cy *out)
{
    if (small(f, a) && small(f, b)) {
        cy r;
        if (a->den == b->den) {
            r.den = a->den;
            for (int i = 0; i < f->phi; i++) r.c[i] = a->c[i] + b->c[i];
        } else {
            i128 g = gcd_i128(a->den, b->den), la = b->den / g, lb = a->den / g;
            r.den = a->den * la;
            for (int i = 0; i < f->phi; i++) r.c[i] = a->c[i] * la + b->c[i] * lb;
        }
        normalize(f, &r);
        cy_copy(f, out, &r);
        return;
    }
    i128 g = ck_gcd(a->den, b->den), la = b->den / g, lb = a->den / g;
    cy r;
    r.den = ck_mul(a->den, la);
    for (int i = 0; i < f->phi; i++) r.c[i] = ck_add(ck_mul(a->c[i], la), ck_mul(b->c[i], lb));
    normalize(f, &r);
    cy_copy(f, out, &r);
}

void cy_sub(const cy_field *f, const cy *a, const cy *b, cy *out)
{
    cy nb;
    for (int i = 0; i < f->phi; i++) nb.c[i] = -b->c[i];
    nb.den = b->den;
    cy_add(f, a, &nb, out);
}

void cy_mul(const cy_field *f, const cy *a, const cy *b, cy *out)
{
    i128 r[2 * CY_MAX_PHI];
    memset(r, 0, sizeof(i128) * (size_t)(2 * f->phi));
    if (small(f, a) && small(f, b)) {
        for (int i = 0; i < f->phi; i++) {
            if (!a->c[i]) continue;
            for (int j = 0; j < f->phi; j++) r[i + j] += a->c[i] * b->c[j];
        }
        for (int k = 2 * f->phi - 1; k >= f->phi; k--) {
            i128 t = r[k];
            if (!t) continue;
            r[k] = 0;
            for (int j = 0; j < f->phi; j++)
                if (f->poly[j]) r[k - f->phi + j] -= t * f->poly[j];
        }
        cy o;
        memcpy(o.c, r, sizeof(i128) * (size_t)f->phi);
        o.den = a->den * b->den;
        normalize(f, &o);
        cy_copy(f, out, &o);
        return;
    }
    for (int i = 0; i < f->phi; i++) {
        if (!a->c[i]) continue;
        for (int j = 0; j < f->phi; j++)
            if (b->c[j]) r[i + j] = ck_add(r[i + j], ck_mul(a->c[i], b->c[j]));
    }
    reduce(f, r, 2 * f->phi);
    cy o;
    memcpy(o.c, r, sizeof(i128) * (size_t)f->phi);
    o.den = ck_mul(a->den, b->den);
    normalize(f, &o);
    cy_copy(f, out, &o);
}

int cy_is_zero(const cy_field *f, const cy *a)
{
    for (int i = 0; i < f->phi; i++)
        if (a->c[i]) return 0;
    return 1;
}

void cy_value(const cy_field *f, const cy *a, double *re, double *im)
{
    long double x = 0, y = 0;
    for (int i = 0; i < f->phi; i++) {
        if (!a->c[i]) continue;
        long double angle = 2.0L * 3.14159265358979323846264338327950288L * i / f->m;
        x += (long double)a->c[i] * cosl(angle);
        y += (long double)a->c[i] * sinl(angle);
    }
    *re = (double)(x / (long double)a->den);
    *im = (double)(y / (long double)a->den);
}

static void half_pair(const cy_field *f, int64_t e1, int64_t e2, int sign, cy *out)
{
    cy a, b, half;
    cy_root(f, e1, &a);
    cy_root(f, e2, &b);
    if (sign < 0)
        cy_sub(f, &a, &b, &a);
    else
        cy_add(f, &a, &b, &a);
    cy_rat(f, 1, 2, &half);
    cy_mul(f, &a, &half, out);
}

void cy_from_sqrt23(const cy_field *f, const int64_t *words, cy *out)
{
    int m = f->m;
    cy_zero(f, out);
    if (m % 24) {
        cy_overflow = 1;
        return;
    }
    cy basis[4], two, term;
    cy_rat(f, 1, 1, &basis[0]);
    half_pair(f, m / 8, -m / 8, 1, &basis[1]);
    half_pair(f, m / 12, -m / 12, 1, &basis[2]);
    cy_rat(f, 2, 1, &two);
    cy_mul(f, &basis[1], &two, &basis[1]);
    cy_mul(f, &basis[2], &two, &basis[2]);
    cy_mul(f, &basis[1], &basis[2], &basis[3]);
    for (int i = 0; i < 4; i++) {
        if (!words[2 * i]) continue;
        cy q;
        cy_rat(f, words[2 * i], words[2 * i + 1], &q);
        cy_mul(f, &q, &basis[i], &term);
        cy_add(f, out, &term, out);
    }
}

int cy_cos_sin(const cy_field *f, int64_t n, int64_t d, cy *c, cy *s)
{
    int64_t m = f->m;
    if (d <= 0 || m % d || m % 4) return 1;
    int64_t e = (n % d) * (m / d);
    half_pair(f, e, -e, 1, c);
    half_pair(f, e + 3 * m / 4, -e + 3 * m / 4, -1, s);
    return 0;
}

/* a^{-1}: solve mult(a) x = e_0 over Q by Gauss-Jordan on (numerator, denominator)
 * pairs; returns nonzero for a = 0 or on overflow. */
int cy_inv(const cy_field *f, const cy *a, cy *out)
{
    int n = f->phi;
    if (cy_is_zero(f, a)) return 1;
    i128 (*mn)[CY_MAX_PHI + 1] = calloc((size_t)n, sizeof *mn), (*md)[CY_MAX_PHI + 1] = calloc((size_t)n, sizeof *md);
    if (!mn || !md) {
        free(mn);
        free(md);
        return 1;
    }
    cy col, zj;
    for (int j = 0; j < n; j++) {
        cy_root(f, j, &zj);
        cy_mul(f, a, &zj, &col);
        for (int i = 0; i < n; i++) {
            mn[i][j] = col.c[i];
            md[i][j] = col.den;
        }
    }
    for (int i = 0; i < n; i++) {
        mn[i][n] = i == 0;
        md[i][n] = 1;
    }
    int status = 0;
    for (int c = 0; c < n && !status; c++) {
        int p = -1;
        for (int r = c; r < n && p < 0; r++)
            if (mn[r][c]) p = r;
        if (p < 0) {
            status = 1;
            break;
        }
        if (p != c)
            for (int k = 0; k <= n; k++) {
                i128 t = mn[p][k], u = md[p][k];
                mn[p][k] = mn[c][k];
                md[p][k] = md[c][k];
                mn[c][k] = t;
                md[c][k] = u;
            }
        i128 pn = mn[c][c], pd = md[c][c];
        for (int k = 0; k <= n; k++) {
            /* row_c[k] / pivot */
            i128 nn = ck_mul(mn[c][k], pd), dd = ck_mul(md[c][k], pn);
            if (dd < 0) {
                nn = -nn;
                dd = -dd;
            }
            i128 g = ck_gcd(nn, dd);
            if (g > 1) {
                nn /= g;
                dd /= g;
            }
            mn[c][k] = nn;
            md[c][k] = nn ? dd : 1;
        }
        for (int r = 0; r < n; r++) {
            if (r == c || !mn[r][c]) continue;
            i128 fn = mn[r][c], fd = md[r][c];
            for (int k = 0; k <= n; k++) {
                if (!mn[c][k]) continue;
                /* row_r[k] -= f * row_c[k] */
                i128 tn = ck_mul(fn, mn[c][k]), td = ck_mul(fd, md[c][k]);
                i128 g = ck_gcd(td, md[r][k]);
                i128 nn = ck_add(ck_mul(mn[r][k], td / g), -ck_mul(tn, md[r][k] / g)), dd = ck_mul(md[r][k], td / g);
                i128 h = ck_gcd(nn, dd);
                if (h > 1) {
                    nn /= h;
                    dd /= h;
                }
                mn[r][k] = nn;
                md[r][k] = nn ? dd : 1;
            }
        }
        if (cy_overflow) status = 1;
    }
    if (!status) {
        i128 l = 1;
        for (int i = 0; i < n; i++) l = ck_mul(l / ck_gcd(l, md[i][n]), md[i][n]);
        out->den = l;
        for (int i = 0; i < n; i++) out->c[i] = ck_mul(mn[i][n], l / md[i][n]);
        normalize(f, out);
        if (cy_overflow) status = 1;
    }
    free(mn);
    free(md);
    return status;
}

/* double-double arithmetic for correctly rounded values */
typedef struct {
    double hi, lo;
} dd;

static dd dd_two_sum(double a, double b)
{
    double s = a + b, v = s - a;
    return (dd){s, (a - (s - v)) + (b - v)};
}

static dd dd_quick(double a, double b)
{
    double s = a + b;
    return (dd){s, b - (s - a)};
}

static dd dd_add(dd a, dd b)
{
    dd s = dd_two_sum(a.hi, b.hi), t = dd_two_sum(a.lo, b.lo);
    s.lo += t.hi;
    s = dd_quick(s.hi, s.lo);
    s.lo += t.lo;
    return dd_quick(s.hi, s.lo);
}

static dd dd_mul(dd a, dd b)
{
    double p = a.hi * b.hi, e = fma(a.hi, b.hi, -p);
    e += a.hi * b.lo + a.lo * b.hi;
    return dd_quick(p, e);
}

static dd dd_div(dd a, dd b)
{
    double q1 = a.hi / b.hi;
    dd r = dd_add(a, dd_mul((dd){-q1, 0.0}, b));
    double q2 = r.hi / b.hi;
    r = dd_add(r, dd_mul((dd){-q2, 0.0}, b));
    double q3 = r.hi / b.hi;
    return dd_add(dd_quick(q1, q2), (dd){q3, 0.0});
}

static dd dd_from_i128(i128 v)
{
    double hi = (double)v;
    return dd_quick(hi, (double)(v - (i128)hi));
}

static const dd DD_2PI = {6.283185307179586232, 2.449293598294706414e-16};

/* sin and cos of 2 pi n / d for 0 <= n / d <= 1 / 8 by Taylor series */
static void dd_sincos_small(int64_t n, int64_t d, dd *s, dd *c)
{
    dd x = dd_div(dd_mul(DD_2PI, dd_from_i128(n)), dd_from_i128(d)), x2 = dd_mul(x, x);
    dd term = x, sum = x;
    for (int k = 1; k < 30; k++) {
        term = dd_div(dd_mul(term, x2), dd_from_i128(-(i128)(2 * k) * (2 * k + 1)));
        sum = dd_add(sum, term);
    }
    *s = sum;
    term = (dd){1.0, 0.0};
    sum = term;
    for (int k = 1; k < 30; k++) {
        term = dd_div(dd_mul(term, x2), dd_from_i128(-(i128)(2 * k - 1) * (2 * k)));
        sum = dd_add(sum, term);
    }
    *c = sum;
}

/* cos and sin of 2 pi n / d */
static void dd_sincos(int64_t n, int64_t d, dd *s, dd *c)
{
    n %= d;
    if (n < 0) n += d;
    int64_t q = 4 * n / d, un = 4 * n - q * d, ud = 4 * d; /* u = un / ud in [0, 1/4) */
    dd su, cu;
    if (2 * un <= d)
        dd_sincos_small(un, ud, &su, &cu);
    else
        dd_sincos_small(d - un, ud, &cu, &su);
    dd ns = {-su.hi, -su.lo}, nc = {-cu.hi, -cu.lo};
    *c = q == 0 ? cu : q == 1 ? ns : q == 2 ? nc : su;
    *s = q == 0 ? su : q == 1 ? cu : q == 2 ? ns : nc;
}

static dd field_cos(const cy_field *f, int i)
{
    if (!f->cosines) {
        double *table = malloc(sizeof(double) * 2 * (size_t)f->phi);
        if (!table) {
            dd s, c;
            dd_sincos(i, f->m, &s, &c);
            return c;
        }
        for (int k = 0; k < f->phi; k++) {
            dd s, c;
            dd_sincos(k, f->m, &s, &c);
            table[2 * k] = c.hi;
            table[2 * k + 1] = c.lo;
        }
        ((cy_field *)f)->cosines = table;
    }
    return (dd){f->cosines[2 * i], f->cosines[2 * i + 1]};
}

/* The correctly rounded double of a real element; nonzero when it is not
 * real or lies too close to a rounding boundary to decide. */
int cy_to_double(const cy_field *f, const cy *a, double *out)
{
    if (cy_is_zero(f, a)) {
        *out = 0.0;
        return 0;
    }
    /* real exactly: a equals its complex conjugate z -> z^-1 */
    cy conj, term;
    cy_zero(f, &conj);
    for (int i = 0; i < f->phi; i++) {
        if (!a->c[i]) continue;
        cy_root(f, -i, &term);
        cy q;
        cy_rat(f, a->c[i], a->den, &q);
        cy_mul(f, &q, &term, &term);
        cy_add(f, &conj, &term, &conj);
    }
    cy diff;
    cy_sub(f, a, &conj, &diff);
    if (!cy_is_zero(f, &diff) || cy_overflow) return 1;
    dd sum = {0.0, 0.0};
    double magnitude = 0.0;
    for (int i = 0; i < f->phi; i++) {
        if (!a->c[i]) continue;
        sum = dd_add(sum, dd_mul(dd_from_i128(a->c[i]), field_cos(f, i)));
        magnitude += fabs((double)a->c[i]);
    }
    sum = dd_div(sum, dd_from_i128(a->den));
    magnitude /= fabs((double)a->den);
    double hi = sum.hi, lo = sum.lo, err = magnitude * 0x1p-96;
    double ulp = fabs(nextafter(hi, (lo > 0) == (hi > 0) ? hi * 2 : 0.0) - hi);
    if (fabs(fabs(lo) - ulp / 2) <= err || fabs(lo) > ulp / 2 + err) return 2;
    if (hi == 0.0 && fabs(lo) > err) return 2;
    *out = hi;
    return 0;
}

/* The sign of a real element: -1, 0 or 1; nonzero return when it cannot be
 * decided at double-double precision. */
int cy_sign(const cy_field *f, const cy *a, int *sign)
{
    if (cy_is_zero(f, a)) {
        *sign = 0;
        return 0;
    }
    dd sum = {0.0, 0.0};
    double magnitude = 0.0;
    for (int i = 0; i < f->phi; i++) {
        if (!a->c[i]) continue;
        sum = dd_add(sum, dd_mul(dd_from_i128(a->c[i]), field_cos(f, i)));
        magnitude += fabs((double)a->c[i]);
    }
    if (fabs(sum.hi) <= magnitude * 0x1p-96) return 1;
    *sign = sum.hi > 0 ? 1 : -1;
    return 0;
}
