/* Exact rationals over int128 with overflow tracking, one status per
 * translation unit (reset by each kernel entry). */

#ifndef APOSTRUCT_RATIONAL_H
#define APOSTRUCT_RATIONAL_H

#include <stdint.h>

#include "intmath.h"

typedef struct {
    i128 n, d;
} rat;

static int rat_overflow;

static inline i128 rat_gcd(i128 a, i128 b)
{
    i128 g = gcd_i128(a, b);
    if (g < 0) {
        rat_overflow = 1;
        return 1;
    }
    return g;
}

static inline i128 rat_neg_i(i128 a)
{
    if (a == I128_MIN) {
        rat_overflow = 1;
        return 0;
    }
    return -a;
}

static inline i128 rat_mul_i(i128 a, i128 b)
{
    i128 out;
    if (__builtin_mul_overflow(a, b, &out)) {
        rat_overflow = 1;
        return 0;
    }
    return out;
}

static inline i128 rat_add_i(i128 a, i128 b)
{
    i128 out;
    if (__builtin_add_overflow(a, b, &out)) {
        rat_overflow = 1;
        return 0;
    }
    return out;
}

/* a / b for b > 0, in 64-bit arithmetic whenever both fit. */
static inline i128 rat_quot(i128 a, i128 b)
{
    if (a == (int64_t)a && b == (int64_t)b) return (int64_t)a / (int64_t)b;
    return a / b;
}

static inline rat rat_make(i128 n, i128 d)
{
    if (!d) {
        rat_overflow = 1;
        return (rat){0, 1};
    }
    if (d < 0) {
        n = rat_neg_i(n);
        d = rat_neg_i(d);
    }
    i128 g = rat_gcd(n, d);
    if (g > 1) {
        n = rat_quot(n, g);
        d = rat_quot(d, g);
    }
    return (rat){n, d};
}

static inline rat rat_int(i128 v) { return (rat){v, 1}; }
static inline rat rat_neg(rat a) { return (rat){rat_neg_i(a.n), a.d}; }
static inline int rat_zero(rat a) { return a.n == 0; }
static inline int rat_eq(rat a, rat b) { return a.n == b.n && a.d == b.d; }

static inline rat rat_add(rat a, rat b)
{
    i128 g = rat_gcd(a.d, b.d);
    return rat_make(rat_add_i(rat_mul_i(a.n, rat_quot(b.d, g)), rat_mul_i(b.n, rat_quot(a.d, g))),
                    rat_mul_i(rat_quot(a.d, g), b.d));
}

static inline rat rat_sub(rat a, rat b) { return rat_add(a, rat_neg(b)); }
static inline rat rat_mul(rat a, rat b) { return rat_make(rat_mul_i(a.n, b.n), rat_mul_i(a.d, b.d)); }
static inline rat rat_div(rat a, rat b) { return rat_make(rat_mul_i(a.n, b.d), rat_mul_i(a.d, b.n)); }

static inline rat rat_mod1(rat a)
{
    i128 r = a.n % a.d;
    if (r < 0) r += a.d;
    return (rat){r, a.d};
}

static inline int rat_cmp(rat a, rat b)
{
    i128 left = rat_mul_i(a.n, b.d), right = rat_mul_i(b.n, a.d);
    return (left > right) - (left < right);
}

static inline int rat_fits64(rat a)
{
    return a.n <= INT64_MAX && a.n >= INT64_MIN && a.d <= INT64_MAX;
}

/* An int128 as (high, low) int64 words. */
static inline void write_i128(int64_t *out, i128 value)
{
    out[0] = (int64_t)(uint64_t)((u128)value >> 64);
    out[1] = (int64_t)(uint64_t)value;
}

static inline i128 read_i128(const int64_t *in) { return (i128)(((u128)(uint64_t)in[0] << 64) | (uint64_t)in[1]); }

/* Three rationals as (numerator, denominator) int128 words. */
static inline void write_rat3(int64_t *out, const rat *v)
{
    for (int axis = 0; axis < 3; axis++) {
        write_i128(out + 4 * axis, v[axis].n);
        write_i128(out + 4 * axis + 2, v[axis].d);
    }
}

static inline void read_rat3(const int64_t *in, rat *out)
{
    for (int axis = 0; axis < 3; axis++) out[axis] = rat_make(read_i128(in + 4 * axis), read_i128(in + 4 * axis + 2));
}

/* row . matrix for a 3-vector and a row-major 3x3. */
static inline void rat_row_mul(const rat *v, const rat (*m)[3], rat *out)
{
    rat t[3];
    for (int col = 0; col < 3; col++) {
        rat acc = rat_int(0);
        for (int row = 0; row < 3; row++) acc = rat_add(acc, rat_mul(v[row], m[row][col]));
        t[col] = acc;
    }
    for (int col = 0; col < 3; col++) out[col] = t[col];
}

static inline int rat_inverse3(const rat (*m)[3], rat (*out)[3])
{
    rat det = rat_int(0);
    for (int col = 0; col < 3; col++)
        det = rat_add(det, rat_mul(m[0][col], rat_sub(rat_mul(m[1][(col + 1) % 3], m[2][(col + 2) % 3]),
                                                       rat_mul(m[1][(col + 2) % 3], m[2][(col + 1) % 3]))));
    if (rat_zero(det)) return 1;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            int r1 = (col + 1) % 3, r2 = (col + 2) % 3, c1 = (row + 1) % 3, c2 = (row + 2) % 3;
            out[row][col] = rat_div(rat_sub(rat_mul(m[r1][c1], m[r2][c2]), rat_mul(m[r1][c2], m[r2][c1])), det);
        }
    return 0;
}

#endif
