/* Integer helpers shared by the exact kernels. */

#ifndef APOSTRUCT_INTMATH_H
#define APOSTRUCT_INTMATH_H

#include <stdint.h>

typedef __int128 i128;
typedef unsigned __int128 u128;

#define I128_MIN (-(i128)(((u128)1 << 127) - 1) - 1)

static inline u128 abs_u128(i128 a) { return a < 0 ? -(u128)a : (u128)a; }

static inline uint64_t gcd_u64(uint64_t a, uint64_t b)
{
    if (!a) return b;
    if (!b) return a;
    int shift = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b) {
            uint64_t t = a;
            a = b;
            b = t;
        }
        b -= a;
    } while (b);
    return a << shift;
}

/* gcd(|a|, |b|), negative when it is 2^127; the software int128 division
 * only runs while a value exceeds 64 bits. */
static inline i128 gcd_i128(i128 signed_a, i128 signed_b)
{
    u128 a = abs_u128(signed_a), b = abs_u128(signed_b);
    while (b) {
        if (!(a >> 64) && !(b >> 64)) return (i128)gcd_u64((uint64_t)a, (uint64_t)b);
        u128 t = a % b;
        a = b;
        b = t;
    }
    return a >> 127 ? I128_MIN : (i128)a;
}

#endif
