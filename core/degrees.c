/* Lattice-angle trigonometry without libm, so every platform gets the same
 * bits. The kernels and acos are FreeBSD msun's (k_sin.c, k_cos.c, e_acos.c):
 *
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this software is freely
 * granted, provided that this notice is preserved. */

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "geometry.h"

static double kernel_sin(double x, double y)
{
    static const double S1 = -1.66666666666666324348e-01, S2 = 8.33333333332248946124e-03,
                        S3 = -1.98412698298579493134e-04, S4 = 2.75573137070700676789e-06,
                        S5 = -2.50507602534068634195e-08, S6 = 1.58969099521155010221e-10;
    double z = x * x, w = z * z, r = S2 + z * (S3 + z * S4) + z * w * (S5 + z * S6), v = z * x;
    return x - ((z * (0.5 * y - v * r) - y) - v * S1);
}

static double kernel_cos(double x, double y)
{
    static const double C1 = 4.16666666666666019037e-02, C2 = -1.38888888888741095749e-03,
                        C3 = 2.48015872894767294178e-05, C4 = -2.75573143513906633035e-07,
                        C5 = 2.08757232129817482790e-09, C6 = -1.13596475577881948265e-11;
    double z = x * x, w = z * z, r = z * (C1 + z * (C2 + z * C3)) + w * w * (C4 + z * (C5 + z * C6));
    double hz = 0.5 * z;
    w = 1.0 - hz;
    return w + (((1.0 - w) - hz) + (z * r - x * y));
}

/* degrees - 90 q (exact) as radians in two parts, and q mod 4 */
static int reduce(double degrees, double *hi, double *lo)
{
    static const double C_HI = 0x1.1df46a2529d39p-6, C_LO = 0x1.5c1d8becdd291p-62;
    double q = nearbyint(degrees / 90.0), t = degrees - 90.0 * q;
    double p = t * C_HI, e = fma(t, C_HI, -p) + t * C_LO;
    *hi = p + e;
    *lo = e - (*hi - p);
    int64_t n = (int64_t)q % 4;
    return (int)(n < 0 ? n + 4 : n);
}

double geo_cos_deg(double degrees)
{
    if (degrees == 0.0) return 1.0;
    if (degrees == 60.0) return 0.5;
    if (degrees == 90.0) return 0.0;
    if (degrees == 120.0) return -0.5;
    if (degrees == 180.0) return -1.0;
    double hi, lo;
    switch (reduce(degrees, &hi, &lo)) {
    case 0: return kernel_cos(hi, lo);
    case 1: return -kernel_sin(hi, lo);
    case 2: return -kernel_cos(hi, lo);
    default: return kernel_sin(hi, lo);
    }
}

double geo_sin_deg(double degrees)
{
    if (degrees == 0.0 || degrees == 180.0) return 0.0;
    if (degrees == 60.0 || degrees == 120.0) return sqrt(3.0) * 0.5;
    if (degrees == 90.0) return 1.0;
    double hi, lo;
    switch (reduce(degrees, &hi, &lo)) {
    case 0: return kernel_sin(hi, lo);
    case 1: return kernel_cos(hi, lo);
    case 2: return -kernel_sin(hi, lo);
    default: return -kernel_cos(hi, lo);
    }
}

static double acos_radians(double x)
{
    static const double pi = 3.14159265358979311600e+00, pio2_hi = 1.57079632679489655800e+00,
                        pio2_lo = 6.12323399573676603587e-17, pS0 = 1.66666666666666657415e-01,
                        pS1 = -3.25565818622400915405e-01, pS2 = 2.01212532134862925881e-01,
                        pS3 = -4.00555345006794114027e-02, pS4 = 7.91534994289814532176e-04,
                        pS5 = 3.47933107596021167570e-05, qS1 = -2.40339491173441421878e+00,
                        qS2 = 2.02094576023350569471e+00, qS3 = -6.88283971605453293030e-01,
                        qS4 = 7.70381505559019352791e-02;
    if (!(fabs(x) <= 1.0)) return (x - x) / (x - x);
    if (x == 1.0) return 0.0;
    if (x == -1.0) return pi + 2.0 * pio2_lo;
    if (fabs(x) < 0.5) {
        if (fabs(x) <= 0x1p-57) return pio2_hi + pio2_lo;
        double z = x * x;
        double p = z * (pS0 + z * (pS1 + z * (pS2 + z * (pS3 + z * (pS4 + z * pS5)))));
        double q = 1.0 + z * (qS1 + z * (qS2 + z * (qS3 + z * qS4)));
        return pio2_hi - (x - (pio2_lo - x * (p / q)));
    }
    double z = (x < 0 ? 1.0 + x : 1.0 - x) * 0.5;
    double p = z * (pS0 + z * (pS1 + z * (pS2 + z * (pS3 + z * (pS4 + z * pS5)))));
    double q = 1.0 + z * (qS1 + z * (qS2 + z * (qS3 + z * qS4)));
    double s = sqrt(z), r = p / q;
    if (x < 0) return pi - 2.0 * (s + (r * s - pio2_lo));
    uint64_t bits;
    double df = s;
    memcpy(&bits, &df, sizeof bits);
    bits &= 0xffffffff00000000ULL;
    memcpy(&df, &bits, sizeof df);
    double c = (z - df * df) / (s + df);
    return 2.0 * (df + (r * s + c));
}

double geo_acos_deg(double cosine)
{
    if (cosine == 1.0) return 0.0;
    if (cosine == 0.5) return 60.0;
    if (cosine == 0.0) return 90.0;
    if (cosine == -0.5) return 120.0;
    if (cosine == -1.0) return 180.0;
    static const double R_HI = 0x1.ca5dc1a63c1f8p+5, R_LO = -0x1.1e7ab456405f9p-49;
    double r = acos_radians(cosine);
    return fma(r, R_HI, r * R_LO);
}
