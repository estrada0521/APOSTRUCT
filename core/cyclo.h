/* Exact arithmetic in the cyclotomic field Q(zeta_M): an element is
 * (c_0 + c_1 z + ... + c_{phi-1} z^{phi-1}) / den over int128, reduced
 * modulo the M-th cyclotomic polynomial. Overflow sets cy_overflow. */
#ifndef APOSTRUCT_CYCLO_H
#define APOSTRUCT_CYCLO_H

#include <stdint.h>

#include "intmath.h"

#define CY_MAX_PHI 192

typedef struct {
    int m, phi;
    int64_t poly[CY_MAX_PHI + 1]; /* monic Phi_M, poly[phi] = 1 */
    int64_t *powers;              /* m x phi: z^e reduced, built on first use */
    double *cosines;              /* phi x (hi, lo): cos(2 pi i / m) in double-double, built on first use */
} cy_field;

typedef struct {
    i128 c[CY_MAX_PHI];
    i128 den;
} cy;

extern int cy_overflow;

/* The field of the m-th roots of unity; 0 when phi(m) exceeds CY_MAX_PHI. */
const cy_field *cy_field_for(int m);

void cy_zero(const cy_field *f, cy *out);
void cy_copy(const cy_field *f, cy *out, const cy *a);
void cy_rat(const cy_field *f, i128 n, i128 d, cy *out);
/* z^e for any integer e. */
void cy_root(const cy_field *f, int64_t e, cy *out);
void cy_add(const cy_field *f, const cy *a, const cy *b, cy *out);
void cy_sub(const cy_field *f, const cy *a, const cy *b, cy *out);
void cy_mul(const cy_field *f, const cy *a, const cy *b, cy *out);
int cy_is_zero(const cy_field *f, const cy *a);
/* a^{-1}; nonzero for a = 0 or overflow. */
int cy_inv(const cy_field *f, const cy *a, cy *out);
/* The correctly rounded double of a real element; nonzero when not real or
 * undecidable at double-double precision. */
int cy_to_double(const cy_field *f, const cy *a, double *out);
/* The sign of a real element; nonzero when undecidable. */
int cy_sign(const cy_field *f, const cy *a, int *sign);
/* a^{-1} certified by exact multiplication (exact_kernel.c); nonzero for 0 or failure. */
int exact_inverse(const cy_field *f, const cy *a, cy *out);
/* The value under z -> exp(2 pi i / m). */
void cy_value(const cy_field *f, const cy *a, double *re, double *im);

/* a + b sqrt2 + c sqrt3 + d sqrt6 from four rationals (n, d pairs). */
void cy_from_sqrt23(const cy_field *f, const int64_t *words, cy *out);
/* cos(2 pi q) and sin(2 pi q) for a rational q = n / d with d | m. */
int cy_cos_sin(const cy_field *f, int64_t n, int64_t d, cy *c, cy *s);

#endif
