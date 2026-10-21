/* CPython number semantics on doubles and exact rationals: float(Fraction),
 * float %, Fraction(str(x)), limit_denominator and _input_fraction. */

#ifndef APOSTRUCT_PYNUM_H
#define APOSTRUCT_PYNUM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rational.h"

/* int / int true division, correctly rounded as CPython does. */
static inline double rat_to_double(rat a)
{
    if (!a.n) return 0.0;
    int negative = a.n < 0;
    u128 n = abs_u128(a.n), d = (u128)a.d;
    if (n < ((u128)1 << 53) && d < ((u128)1 << 53)) {
        double value = (double)(uint64_t)n / (double)(uint64_t)d;
        return negative ? -value : value;
    }
    int n_bits = 128 - (int)(n >> 64 ? __builtin_clzll((uint64_t)(n >> 64)) : 64 + __builtin_clzll((uint64_t)n));
    int d_bits = 128 - (int)(d >> 64 ? __builtin_clzll((uint64_t)(d >> 64)) : 64 + __builtin_clzll((uint64_t)d));
    int exponent = n_bits - d_bits - 55;
    u128 quotient = 0, remainder = 0;
    if (exponent >= 0) {
        u128 divisor = d << exponent;
        quotient = n / divisor;
        remainder = n % divisor;
    } else {
        int shift = -exponent;
        u128 high = shift >= 128 ? 0 : n >> (128 - shift), low = shift >= 128 ? 0 : n << shift;
        for (int bit = 255; bit >= 0; bit--) {
            int set = bit >= 128 ? (int)((high >> (bit - 128)) & 1) : (int)((low >> bit) & 1);
            remainder = (remainder << 1) | (u128)set;
            quotient <<= 1;
            if (remainder >= d) {
                remainder -= d;
                quotient |= 1;
            }
        }
    }
    while (quotient >= ((u128)1 << 55)) {
        remainder |= quotient & 1;
        quotient >>= 1;
        exponent++;
    }
    int guard = (int)(quotient & 3);
    quotient >>= 2;
    exponent += 2;
    if (guard > 2 || (guard == 2 && (remainder || (quotient & 1)))) quotient += 1;
    double value = ldexp((double)(uint64_t)quotient, exponent);
    return negative ? -value : value;
}

static inline double pymod1(double x)
{
    double mod = fmod(x, 1.0);
    if (mod != 0.0) {
        if (mod < 0.0) mod += 1.0;
    } else {
        mod = copysign(0.0, 1.0);
    }
    return mod;
}

/* Fraction(str(value)) of a decimal string. */
static inline rat decimal_rat(const char *text)
{
    int negative = 0;
    if (*text == '-' || *text == '+') negative = *text++ == '-';
    i128 digits = 0;
    int scale = 0, seen_point = 0;
    for (; *text && *text != 'e' && *text != 'E'; text++) {
        if (*text == '.') {
            seen_point = 1;
            continue;
        }
        digits = rat_add_i(rat_mul_i(digits, 10), *text - '0');
        if (seen_point) scale++;
    }
    if (*text) scale -= atoi(text + 1);
    i128 power = 1;
    for (int i = 0; i < (scale < 0 ? -scale : scale); i++) power = rat_mul_i(power, 10);
    rat out = scale >= 0 ? rat_make(digits, power) : rat_int(rat_mul_i(digits, power));
    return negative ? rat_neg(out) : out;
}

static inline rat limit_denominator(rat value, i128 max_denominator)
{
    if (value.d <= max_denominator) return value;
    i128 p0 = 0, q0 = 1, p1 = 1, q1 = 0, n = value.n, d = value.d;
    for (;;) {
        i128 a = n / d;
        if ((n % d) && ((n < 0) != (d < 0))) a -= 1;
        i128 q2 = q0 + a * q1;
        if (q2 > max_denominator) break;
        i128 p2 = p0 + a * p1;
        p0 = p1;
        q0 = q1;
        p1 = p2;
        q1 = q2;
        i128 next = n - a * d;
        n = d;
        d = next;
    }
    i128 k = (max_denominator - q0) / q1;
    if (rat_mul_i(rat_mul_i(2, d), q0 + k * q1) <= value.d) return rat_make(p1, q1);
    return rat_make(p0 + k * p1, q0 + k * q1);
}

/* Fraction(str(value)): the decimal value of the shortest repr of a double. */
static inline rat repr_rat(double value)
{
    char buffer[64];
    for (int precision = 1; precision <= 17; precision++) {
        snprintf(buffer, sizeof buffer, "%.*g", precision, value);
        if (strtod(buffer, 0) == value) break;
    }
    return decimal_rat(buffer);
}

/* ordinary_presentation._input_fraction of a presentation float, mod 1. */
static inline rat input_fraction_mod1(double value)
{
    char buffer[64];
    snprintf(buffer, sizeof buffer, "%.9f", value);
    rat decimal = repr_rat(strtod(buffer, 0));
    rat crystallographic = limit_denominator(decimal, 384);
    rat gap = rat_sub(crystallographic, decimal);
    if (gap.n < 0) gap.n = -gap.n;
    return rat_mod1(rat_cmp(gap, rat_make(1, 200000000)) <= 0 ? crystallographic : decimal);
}

/* round(value, digits), digits in {10, 12}: the correctly rounded decimal,
 * then the nearest double. fma recovers the exact product; only values
 * within 1e-9 of a decimal tie go through the C library conversion. */
static inline double python_round(double value, int digits)
{
    if (value == trunc(value) && fabs(value) < 4.0e15) return value;
    double scale = digits == 12 ? 1e12 : 1e10;
    double scaled = value * scale;
    if (isfinite(scaled) && fabs(scaled) < 4.0e15) {
        double error = fma(value, scale, -scaled);
        double lower = floor(scaled);
        double fraction = (scaled - lower) + error;
        if (fraction < 0.0) { lower -= 1.0; fraction += 1.0; }
        if (fraction >= 1.0) { lower += 1.0; fraction -= 1.0; }
        if (fabs(fraction - 0.5) > 1e-9) {
            double integer = fraction > 0.5 ? lower + 1.0 : lower;
            if (integer == 0.0) return copysign(0.0, value);
            return integer / scale;
        }
    }
    char text[80];
    snprintf(text, sizeof text, "%.*f", digits, value);
    return strtod(text, NULL);
}

/* Unsigned big integers for the exact sums that leave int128. */
#define BIG_LIMBS 96

typedef struct {
    uint32_t limb[BIG_LIMBS];
} big;

static inline void big_set(big *a, uint64_t v)
{
    memset(a, 0, sizeof *a);
    a->limb[0] = (uint32_t)v;
    a->limb[1] = (uint32_t)(v >> 32);
}

static inline int big_mul_small(big *a, uint32_t m)
{
    uint64_t carry = 0;
    for (int i = 0; i < BIG_LIMBS; i++) {
        uint64_t t = (uint64_t)a->limb[i] * m + carry;
        a->limb[i] = (uint32_t)t;
        carry = t >> 32;
    }
    return carry != 0;
}

static inline int big_mul_u128(big *a, unsigned __int128 m)
{
    big copy = *a, part;
    int over = 0;
    memset(a, 0, sizeof *a);
    for (int word = 0; word < 4 && m; word++, m >>= 32) {
        part = copy;
        over |= big_mul_small(&part, (uint32_t)m);
        for (int shift = 0; shift < word; shift++) {
            over |= part.limb[BIG_LIMBS - 1] != 0;
            memmove(part.limb + 1, part.limb, sizeof(uint32_t) * (BIG_LIMBS - 1));
            part.limb[0] = 0;
        }
        uint64_t carry = 0;
        for (int i = 0; i < BIG_LIMBS; i++) {
            uint64_t t = (uint64_t)a->limb[i] + part.limb[i] + carry;
            a->limb[i] = (uint32_t)t;
            carry = t >> 32;
        }
        over |= carry != 0;
    }
    return over;
}

static inline int big_cmp(const big *a, const big *b)
{
    for (int i = BIG_LIMBS - 1; i >= 0; i--)
        if (a->limb[i] != b->limb[i]) return a->limb[i] < b->limb[i] ? -1 : 1;
    return 0;
}

static inline void big_add(big *a, const big *b)
{
    uint64_t carry = 0;
    for (int i = 0; i < BIG_LIMBS; i++) {
        uint64_t t = (uint64_t)a->limb[i] + b->limb[i] + carry;
        a->limb[i] = (uint32_t)t;
        carry = t >> 32;
    }
}

/* a -= b for a >= b */
static inline void big_sub(big *a, const big *b)
{
    int64_t borrow = 0;
    for (int i = 0; i < BIG_LIMBS; i++) {
        int64_t t = (int64_t)a->limb[i] - b->limb[i] - borrow;
        borrow = t < 0;
        a->limb[i] = (uint32_t)(t + (borrow ? ((int64_t)1 << 32) : 0));
    }
}

static inline int big_bits(const big *a)
{
    for (int i = BIG_LIMBS - 1; i >= 0; i--)
        if (a->limb[i]) return 32 * i + 32 - __builtin_clz(a->limb[i]);
    return 0;
}

static inline void big_shl(big *a, int bits)
{
    int words = bits / 32, rest = bits % 32;
    for (int i = BIG_LIMBS - 1; i >= 0; i--) {
        uint64_t v = i - words >= 0 ? (uint64_t)a->limb[i - words] << rest : 0;
        if (rest && i - words - 1 >= 0) v |= (uint64_t)a->limb[i - words - 1] >> (32 - rest);
        a->limb[i] = (uint32_t)v;
    }
}

/* num / den correctly rounded (both nonzero magnitudes). */
static inline double big_ratio(const big *num, const big *den)
{
    big a = *num, b = *den;
    int shift = 55 - (big_bits(&a) - big_bits(&b));
    if (shift > 0) big_shl(&a, shift);
    else if (shift < 0) big_shl(&b, -shift);
    uint64_t quotient = 0;
    for (int i = 57; i >= 0; i--) {
        big t = b;
        big_shl(&t, i);
        if (big_bits(&b) + i <= BIG_LIMBS * 32 && big_cmp(&a, &t) >= 0) {
            big_sub(&a, &t);
            quotient |= (uint64_t)1 << i;
        }
    }
    int sticky = big_bits(&a) != 0, exponent = -shift;
    while (quotient >= ((uint64_t)1 << 55)) {
        sticky |= (int)(quotient & 1);
        quotient >>= 1;
        exponent++;
    }
    int guard = (int)(quotient & 3);
    quotient >>= 2;
    exponent += 2;
    if (guard > 2 || (guard == 2 && (sticky || (quotient & 1)))) quotient += 1;
    return ldexp((double)quotient, exponent);
}

/* The shortest repr of a double as sign * mantissa * 10^exponent. */
static inline void repr_decimal(double value, int *negative, uint64_t *mantissa, int *exponent)
{
    char buffer[64];
    for (int precision = 1; precision <= 17; precision++) {
        snprintf(buffer, sizeof buffer, "%.*e", precision - 1, value);
        if (strtod(buffer, 0) == value) break;
    }
    const char *text = buffer;
    *negative = *text == '-';
    if (*text == '-' || *text == '+') text++;
    uint64_t digits = 0;
    int scale = 0, seen_point = 0;
    for (; *text && *text != 'e'; text++) {
        if (*text == '.') {
            seen_point = 1;
            continue;
        }
        digits = digits * 10 + (uint64_t)(*text - '0');
        if (seen_point) scale++;
    }
    *mantissa = digits;
    *exponent = (*text ? atoi(text + 1) : 0) - scale;
}

/* float(sum((Fraction(str(values[c])) - origin[c]) * matrix[c][axis])) for
 * each axis, exactly at any decimal exponent; returns 1 on a range error. */
static inline int affine_repr_double(const double *values, const rat *origin, const rat (*matrix)[3], double *out)
{
    int negative[3], exponent[3], low = 0;
    uint64_t mantissa[3];
    for (int c = 0; c < 3; c++) {
        repr_decimal(values[c], &negative[c], &mantissa[c], &exponent[c]);
        if (mantissa[c] && exponent[c] < low) low = exponent[c];
    }
    for (int axis = 0; axis < 3; axis++) {
        rat constant = rat_int(0);
        for (int c = 0; c < 3; c++) constant = rat_sub(constant, rat_mul(origin[c], matrix[c][axis]));
        i128 lcm = constant.d;
        for (int c = 0; c < 3; c++) lcm = rat_mul_i(lcm / rat_gcd(lcm, matrix[c][axis].d), matrix[c][axis].d);
        if (rat_overflow) return 1;
        big positive, negative_sum, term, den;
        big_set(&positive, 0);
        big_set(&negative_sum, 0);
        int over = 0;
        for (int c = 0; c <= 3; c++) {
            i128 factor;
            int sign, power;
            if (c < 3) {
                rat m = matrix[c][axis];
                if (!mantissa[c] || !m.n) continue;
                factor = m.n * (lcm / m.d);
                big_set(&term, mantissa[c]);
                sign = negative[c] ^ (factor < 0);
                power = exponent[c] - low;
            } else {
                if (!constant.n) continue;
                factor = constant.n * (lcm / constant.d);
                big_set(&term, 1);
                sign = factor < 0;
                power = -low;
            }
            over |= big_mul_u128(&term, (unsigned __int128)(factor < 0 ? -factor : factor));
            for (int i = 0; i < power; i++) over |= big_mul_small(&term, 10);
            big_add(sign ? &negative_sum : &positive, &term);
        }
        big_set(&den, 1);
        over |= big_mul_u128(&den, (unsigned __int128)lcm);
        for (int i = 0; i < -low; i++) over |= big_mul_small(&den, 10);
        if (over || big_bits(&positive) >= BIG_LIMBS * 32 - 64 || big_bits(&negative_sum) >= BIG_LIMBS * 32 - 64) return 1;
        int sign = big_cmp(&positive, &negative_sum) < 0;
        if (sign) {
            big_sub(&negative_sum, &positive);
            positive = negative_sum;
        } else
            big_sub(&positive, &negative_sum);
        double value = big_bits(&positive) ? big_ratio(&positive, &den) : 0.0;
        out[axis] = sign ? -value : value;
    }
    return 0;
}

/* sum() over floats: int 0 start, Neumaier compensation (CPython >= 3.12). */
static inline double py_sum(const double *values, int count)
{
    if (count <= 0) return 0.0;
    double total = 0.0 + values[0];
    double compensation = 0.0;
    for (int index = 1; index < count; index++) {
        double value = values[index];
        double step = total + value;
        if (fabs(total) >= fabs(value)) compensation += (total - step) + value;
        else compensation += (value - step) + total;
        total = step;
    }
    if (compensation != 0.0 && isfinite(compensation)) total += compensation;
    return total;
}

#endif
