/* Arithmetic modulo a prime below 2^62, with Montgomery multiplication;
 * each translation unit keeps its own current prime (mont_setup). */
#ifndef APOSTRUCT_MODP_H
#define APOSTRUCT_MODP_H

#include <stdint.h>

typedef unsigned __int128 u128;

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t p) { return (uint64_t)((u128)a * b % p); }

static inline uint64_t powmod(uint64_t a, uint64_t e, uint64_t p)
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


/* Montgomery arithmetic modulo the current prime: values are kept as x R mod p, R = 2^64 */
static uint64_t MP, MPINV, MR2, MONE;

static inline uint64_t mm(uint64_t a, uint64_t b)
{
    u128 t = (u128)a * b;
    uint64_t m = (uint64_t)t * MPINV;
    u128 u = (t + (u128)m * MP) >> 64;
    return u >= MP ? (uint64_t)(u - MP) : (uint64_t)u;
}

static inline void mont_setup(uint64_t p)
{
    uint64_t inv = p; /* Newton for p^-1 mod 2^64 */
    for (int i = 0; i < 6; i++) inv *= 2 - p * inv;
    MP = p;
    MPINV = -inv;
    uint64_t r = (uint64_t)(((u128)1 << 64) % p);
    MR2 = mulmod(r, r, p);
    MONE = r;
}

static inline uint64_t to_mont(uint64_t x) { return mm(x, MR2); }
static inline uint64_t from_mont(uint64_t x) { return mm(x, 1); }
static inline uint64_t madd(uint64_t a, uint64_t b) { uint64_t s = a + b; return s >= MP ? s - MP : s; }
static inline uint64_t msub(uint64_t a, uint64_t b) { return a >= b ? a - b : a + MP - b; }

static inline uint64_t minv(uint64_t a)
{
    uint64_t r = MONE, e = MP - 2;
    while (e) {
        if (e & 1) r = mm(r, a);
        a = mm(a, a);
        e >>= 1;
    }
    return r;
}


/* (math/exact_kernel.c) the attempt-th prime p = 1 mod m for f with a primitive
 * m-th root w, the field's units (exponents of the embeddings), and the
 * Montgomery-form inverse Vandermonde matrix vinv[i][e] of x_e = w^units[e];
 * cached. Nonzero when none is available. */
int modp_field(int m, int phi, int attempt, uint64_t *p, uint64_t *w, const int **units, const uint64_t **vinv);
/* x / y = a mod p with |x|, y < sqrt(p / 2) */
int modp_reconstruct(uint64_t a, uint64_t p, __int128 *num, __int128 *den);

#endif
