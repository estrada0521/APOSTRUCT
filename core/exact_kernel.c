/* The kernel of a matrix over Q(zeta_m) in reduced row echelon form, computed
 * modulo a prime p = 1 mod m in every embedding, lifted coefficientwise by
 * rational reconstruction and then verified exactly by multiplication. The
 * verified vectors certify the kernel: a minor nonzero mod p is nonzero over
 * Q(zeta_m), so the true rank is at least the rank mod p. */

#include <stdlib.h>
#include <string.h>

#include "cyclo.h"


typedef unsigned __int128 u128;

static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t p) { return (uint64_t)((u128)a * b % p); }

static uint64_t powmod(uint64_t a, uint64_t e, uint64_t p)
{
    uint64_t r = 1;
    a %= p;
    while (e) {
        if (e & 1) r = mulmod(r, a, p);
        a = mulmod(a, a, p);
        e >>= 1;
    }
    return r;
}

static uint64_t invmod(uint64_t a, uint64_t p) { return powmod(a, p - 2, p); }

static int is_prime(uint64_t n)
{
    if (n < 2) return 0;
    static const uint64_t bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++)
        if (n % bases[i] == 0) return n == bases[i];
    uint64_t d = n - 1;
    int s = 0;
    while (!(d & 1)) {
        d >>= 1;
        s++;
    }
    for (int i = 0; i < 12; i++) {
        uint64_t x = powmod(bases[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (int r = 1; r < s && composite; r++) {
            x = mulmod(x, x, n);
            composite = x != n - 1;
        }
        if (composite) return 0;
    }
    return 1;
}

/* The attempt-th prime p = 1 mod m below 2^62 and a primitive m-th root of unity. */
static int prime_for(int m, int attempt, uint64_t *p_out, uint64_t *root_out)
{
    uint64_t t = ((uint64_t)1 << 62) / (uint64_t)m;
    int found = 0;
    for (; t > 0; t--) {
        uint64_t p = t * (uint64_t)m + 1;
        if (!is_prime(p)) continue;
        if (found++ < attempt) continue;
        /* a generator's (p-1)/m power has exact order m when it is not a
         * smaller root for any prime divisor of m */
        for (uint64_t g = 2; g < 1000; g++) {
            uint64_t w = powmod(g, (p - 1) / (uint64_t)m, p);
            int ok = 1;
            for (int q = 2; q <= m && ok; q++) {
                if (m % q) continue;
                int prime = 1;
                for (int r = 2; r * r <= q && prime; r++) prime = q % r != 0;
                if (prime && powmod(w, (uint64_t)(m / q), p) == 1) ok = 0;
            }
            if (ok) {
                *p_out = p;
                *root_out = w;
                return 0;
            }
        }
    }
    return 1;
}

static uint64_t to_mod(i128 a, uint64_t p)
{
    i128 r = a % (i128)p;
    return (uint64_t)(r < 0 ? r + (i128)p : r);
}

/* x / y = a mod p with |x|, y < sqrt(p / 2); nonzero when none exists. */
static int reconstruct(uint64_t a, uint64_t p, i128 *num, i128 *den)
{
    i128 r0 = p, r1 = a, t0 = 0, t1 = 1;
    u128 bound = 1;
    while (bound * bound < p / 2) bound <<= 1;
    bound >>= 1;
    while ((u128)(r1 < 0 ? -r1 : r1) >= bound) {
        i128 q = r0 / r1, r2 = r0 - q * r1, t2 = t0 - q * t1;
        r0 = r1;
        r1 = r2;
        t0 = t1;
        t1 = t2;
    }
    if (!t1 || (u128)(t1 < 0 ? -t1 : t1) >= bound) return 1;
    if (t1 < 0) {
        t1 = -t1;
        r1 = -r1;
    }
    *num = r1;
    *den = t1;
    return 0;
}

static int units_of(int m, int *units)
{
    int n = 0;
    for (int j = 1; j < m; j++) {
        int a = j, b = m;
        while (b) {
            int t = a % b;
            a = b;
            b = t;
        }
        if (a == 1) units[n++] = j;
    }
    if (m == 1) units[n++] = 0;
    return n;
}

/* rows x cols matrix a (row-major) over f. On success *count kernel vectors
 * (count x cols, row-major) in reduced echelon form: vector t has a 1 at the
 * t-th free column, zeros at the other free columns. pivots_out (cols) marks
 * pivot columns. Returns 0, or 1 when no prime certified the kernel. */
int exact_kernel(const cy_field *f, int rows, int cols, const cy *a, int *count, cy *kernel, unsigned char *pivots_out)
{
    int m = f->m, phi = f->phi;
    int *units = malloc(sizeof(int) * (size_t)(m + 1));
    uint64_t *mat = malloc(sizeof(uint64_t) * (size_t)rows * cols + 1);
    uint64_t *values = malloc(sizeof(uint64_t) * (size_t)phi * cols * cols + 1);
    uint64_t *vinv = malloc(sizeof(uint64_t) * (size_t)phi * phi + 1), *work = malloc(sizeof(uint64_t) * (size_t)phi * 2 * phi + 1);
    int *pivot_col = malloc(sizeof(int) * (size_t)cols + 1);
    unsigned char *is_pivot = malloc((size_t)cols + 1);
    int status = 1;
    if (!units || !mat || !values || !vinv || !work || !pivot_col || !is_pivot) goto done;
    if (units_of(m, units) != phi) goto done;
    for (int attempt = 0; attempt < 4 && status; attempt++) {
        uint64_t p, w;
        if (prime_for(m, attempt, &p, &w)) break;
        int rank = -1, good = 1;
        memset(values, 0, sizeof(uint64_t) * (size_t)phi * cols * cols);
        for (int e = 0; e < phi && good; e++) {
            uint64_t x = powmod(w, (uint64_t)units[e], p), xpow[CY_MAX_PHI];
            xpow[0] = 1;
            for (int i = 1; i < phi; i++) xpow[i] = mulmod(xpow[i - 1], x, p);
            for (int r = 0; r < rows && good; r++)
                for (int c = 0; c < cols; c++) {
                    const cy *v = &a[(size_t)r * cols + c];
                    uint64_t d = to_mod(v->den, p), s = 0;
                    if (!d) good = 0;
                    for (int i = 0; i < phi; i++)
                        if (v->c[i]) s = (s + mulmod(to_mod(v->c[i], p), xpow[i], p)) % p;
                    mat[(size_t)r * cols + c] = good ? mulmod(s, invmod(d, p), p) : 0;
                }
            /* Gauss-Jordan mod p */
            int prow = 0;
            memset(is_pivot, 0, (size_t)cols);
            for (int c = 0; c < cols && prow < rows; c++) {
                int pr = -1;
                for (int r = prow; r < rows && pr < 0; r++)
                    if (mat[(size_t)r * cols + c]) pr = r;
                if (pr < 0) continue;
                if (pr != prow)
                    for (int k = 0; k < cols; k++) {
                        uint64_t t = mat[(size_t)pr * cols + k];
                        mat[(size_t)pr * cols + k] = mat[(size_t)prow * cols + k];
                        mat[(size_t)prow * cols + k] = t;
                    }
                uint64_t inv = invmod(mat[(size_t)prow * cols + c], p);
                for (int k = 0; k < cols; k++) mat[(size_t)prow * cols + k] = mulmod(mat[(size_t)prow * cols + k], inv, p);
                for (int r = 0; r < rows; r++) {
                    uint64_t factor = mat[(size_t)r * cols + c];
                    if (r == prow || !factor) continue;
                    for (int k = 0; k < cols; k++)
                        mat[(size_t)r * cols + k] = (mat[(size_t)r * cols + k] + p - mulmod(factor, mat[(size_t)prow * cols + k], p)) % p;
                }
                pivot_col[prow] = c;
                is_pivot[c] = 1;
                prow++;
            }
            if (e == 0) {
                rank = prow;
                memcpy(pivots_out, is_pivot, (size_t)cols);
            } else if (prow != rank || memcmp(pivots_out, is_pivot, (size_t)cols)) {
                good = 0;
            }
            /* kernel vector t for the t-th free column: 1 there, -row entries at pivots */
            int t = 0;
            for (int fc = 0; fc < cols && good; fc++) {
                if (is_pivot[fc]) continue;
                uint64_t *vec = values + ((size_t)t * cols) * phi;
                for (int c = 0; c < cols; c++) vec[(size_t)c * phi + e] = c == fc;
                for (int r = 0; r < prow; r++) vec[(size_t)pivot_col[r] * phi + e] = (p - mat[(size_t)r * cols + fc]) % p;
                t++;
            }
        }
        if (!good) continue;
        int kcount = cols - rank;
        /* inverse of the Vandermonde matrix V[e][i] = x_e^i */
        for (int e = 0; e < phi; e++) {
            uint64_t x = powmod(w, (uint64_t)units[e], p), xp = 1;
            for (int i = 0; i < phi; i++) {
                work[(size_t)e * 2 * phi + i] = xp;
                work[(size_t)e * 2 * phi + phi + i] = e == i;
                xp = mulmod(xp, x, p);
            }
        }
        for (int c = 0; c < phi && good; c++) {
            int pr = -1;
            for (int r = c; r < phi && pr < 0; r++)
                if (work[(size_t)r * 2 * phi + c]) pr = r;
            if (pr < 0) {
                good = 0;
                break;
            }
            for (int k = 0; k < 2 * phi; k++) {
                uint64_t tmp = work[(size_t)pr * 2 * phi + k];
                work[(size_t)pr * 2 * phi + k] = work[(size_t)c * 2 * phi + k];
                work[(size_t)c * 2 * phi + k] = tmp;
            }
            uint64_t inv = invmod(work[(size_t)c * 2 * phi + c], p);
            for (int k = 0; k < 2 * phi; k++) work[(size_t)c * 2 * phi + k] = mulmod(work[(size_t)c * 2 * phi + k], inv, p);
            for (int r = 0; r < phi; r++) {
                uint64_t factor = work[(size_t)r * 2 * phi + c];
                if (r == c || !factor) continue;
                for (int k = 0; k < 2 * phi; k++)
                    work[(size_t)r * 2 * phi + k] = (work[(size_t)r * 2 * phi + k] + p - mulmod(factor, work[(size_t)c * 2 * phi + k], p)) % p;
            }
        }
        if (!good) continue;
        for (int e = 0; e < phi; e++)
            for (int i = 0; i < phi; i++) vinv[(size_t)i * phi + e] = work[(size_t)i * 2 * phi + phi + e];
        /* lift: coefficients c = V^{-1} values, then rational reconstruction */
        cy_overflow = 0;
        for (int t = 0; t < kcount && good; t++)
            for (int c = 0; c < cols && good; c++) {
                const uint64_t *vals = values + ((size_t)t * cols + c) * phi;
                i128 nums[CY_MAX_PHI], dens[CY_MAX_PHI], l = 1;
                for (int i = 0; i < phi && good; i++) {
                    uint64_t s = 0;
                    for (int e = 0; e < phi; e++) s = (s + mulmod(vinv[(size_t)i * phi + e], vals[e], p)) % p;
                    if (!s) {
                        nums[i] = 0;
                        dens[i] = 1;
                        continue;
                    }
                    if (reconstruct(s, p, &nums[i], &dens[i])) good = 0;
                    i128 g = gcd_i128(l, dens[i]);
                    l = l / g * dens[i];
                }
                if (!good) break;
                cy *out = &kernel[(size_t)t * cols + c];
                out->den = l;
                for (int i = 0; i < phi; i++) out->c[i] = nums[i] * (l / dens[i]);
                cy one;
                cy_rat(f, 1, 1, &one);
                cy_mul(f, out, &one, out);
            }
        if (!good) continue;
        /* exact verification: a v = 0 for every vector */
        for (int t = 0; t < kcount && good; t++)
            for (int r = 0; r < rows && good; r++) {
                cy acc, term;
                cy_zero(f, &acc);
                for (int c = 0; c < cols; c++) {
                    const cy *x = &a[(size_t)r * cols + c], *y = &kernel[(size_t)t * cols + c];
                    if (cy_is_zero(f, x) || cy_is_zero(f, y)) continue;
                    cy_mul(f, x, y, &term);
                    cy_add(f, &acc, &term, &acc);
                }
                good = cy_is_zero(f, &acc) && !cy_overflow;
            }
        if (good) {
            *count = kcount;
            status = 0;
        }
    }
done:
    free(units);
    free(mat);
    free(values);
    free(vinv);
    free(work);
    free(pivot_col);
    free(is_pivot);
    return status;
}
