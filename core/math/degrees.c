/* Lattice-angle trigonometry without libm, so every platform gets the same
 * bits. sin, cos and acos are evaluated in double-double and rounded once
 * (correctly rounded unless the value lies within about 2^-100 of a rounding
 * boundary); cosd and sind follow gfortran 11's (libgfortran trigd.inc), the
 * functions upstream findsym and iso call. */

#include <math.h>

#include "geometry.h"

typedef struct {
    double hi, lo;
} dd;

static const dd SIN_C[15] = {{0x1.0000000000000p+0, 0x0.0p+0}, {-0x1.5555555555555p-3, -0x1.5555555555555p-57}, {0x1.1111111111111p-7, 0x1.1111111111111p-63}, {-0x1.a01a01a01a01ap-13, -0x1.a01a01a01a01ap-73}, {0x1.71de3a556c734p-19, -0x1.c154f8ddc6c00p-73}, {-0x1.ae64567f544e4p-26, 0x1.c062e06d1f209p-80}, {0x1.6124613a86d09p-33, 0x1.f28e0cc748ebep-87}, {-0x1.ae7f3e733b81fp-41, -0x1.1d8656b0ee8cbp-97}, {0x1.952c77030ad4ap-49, 0x1.ac981465ddc6cp-103}, {-0x1.2f49b46814157p-57, -0x1.2650f61dbdcb4p-112}, {0x1.71b8ef6dcf572p-66, -0x1.d043ae40c4647p-120}, {-0x1.761b41316381ap-75, 0x1.3423c7d91404fp-130}, {0x1.3f3ccdd165fa9p-84, -0x1.58ddadf344487p-139}, {-0x1.d1ab1c2dccea3p-94, -0x1.054d0c78aea14p-149}, {0x1.259f98b4358adp-103, 0x1.eaf8c39dd9bc5p-157}};
static const dd COS_C[15] = {{0x1.0000000000000p+0, 0x0.0p+0}, {-0x1.0000000000000p-1, 0x0.0p+0}, {0x1.5555555555555p-5, 0x1.5555555555555p-59}, {-0x1.6c16c16c16c17p-10, 0x1.f49f49f49f49fp-65}, {0x1.a01a01a01a01ap-16, 0x1.a01a01a01a01ap-76}, {-0x1.27e4fb7789f5cp-22, -0x1.cbbc05b4fa99ap-76}, {0x1.1eed8eff8d898p-29, -0x1.2aec959e14c06p-83}, {-0x1.93974a8c07c9dp-37, -0x1.05d6f8a2efd1fp-92}, {0x1.ae7f3e733b81fp-45, 0x1.1d8656b0ee8cbp-101}, {-0x1.6827863b97d97p-53, -0x1.eec01221a8b0bp-107}, {0x1.e542ba4020225p-62, 0x1.ea72b4afe3c2fp-120}, {-0x1.0ce396db7f853p-70, 0x1.aebcdbd20331cp-124}, {0x1.f2cf01972f578p-80, -0x1.9ada5fcc1ab14p-135}, {-0x1.88e85fc6a4e5ap-89, 0x1.71c37ebd16540p-143}, {0x1.0a18a2635085dp-98, 0x1.b9e2e28e1aa54p-153}};
static const dd ASIN_C[56] = {{0x1.0000000000000p+0, 0x0.0p+0}, {0x1.5555555555555p-3, 0x1.5555555555555p-57}, {0x1.3333333333333p-4, 0x1.999999999999ap-59}, {0x1.6db6db6db6db7p-5, -0x1.2492492492492p-60}, {0x1.f1c71c71c71c7p-6, 0x1.c71c71c71c71cp-62}, {0x1.6e8ba2e8ba2e9p-6, -0x1.1745d1745d174p-60}, {0x1.1c4ec4ec4ec4fp-6, -0x1.d89d89d89d89ep-61}, {0x1.c99999999999ap-7, -0x1.999999999999ap-61}, {0x1.7a87878787878p-7, 0x1.e1e1e1e1e1e1ep-61}, {0x1.3fde50d79435ep-7, 0x1.435e50d79435ep-61}, {0x1.12ef3cf3cf3cfp-7, 0x1.e79e79e79e79ep-62}, {0x1.df3bd37a6f4dfp-8, -0x1.90b21642c8591p-62}, {0x1.a6863d70a3d71p-8, -0x1.70a3d70a3d70ap-62}, {0x1.782dda12f684cp-8, -0x1.2f684bda12f68p-63}, {0x1.51ba308d3dcb1p-8, -0x1.cb08d3dcb08d4p-62}, {0x1.31683bdef7bdfp-8, -0x1.0842108421084p-65}, {0x1.15ee9d45d1746p-8, -0x1.745d1745d1746p-63}, {0x1.fcaf8fb6db6dbp-9, 0x1.b6db6db6db6dbp-63}, {0x1.d3d2a8e0dd67dp-9, -0x1.d67c8a60dd67dp-63}, {0x1.b026f57b13b14p-9, -0x1.3b13b13b13b14p-63}, {0x1.90cb77f60c7cep-9, 0x1.8f9c18f9c18fap-66}, {0x1.750de64d7d05fp-9, 0x1.05f417d05f418p-63}, {0x1.5c5f56efaaaabp-9, -0x1.5555555555555p-63}, {0x1.464c0950f7d47p-9, -0x1.882b931057262p-64}, {0x1.3275586c5f2f0p-9, 0x1.4e5e0a72f0539p-63}, {0x1.208d3570ae5a6p-9, -0x1.6969696969697p-63}, {0x1.1052bc5fa960ap-9, -0x1.5bc609a90e7d9p-63}, {0x1.018f963c229bfp-9, -0x1.4f2094f2094f2p-64}, {0x1.e82be60d9127ep-10, -0x1.f7047dc11f704p-65}, {0x1.cf7dea5b6e830p-10, -0x1.15b1e5f75270dp-64}, {0x1.b8d2e5667ce6cp-10, 0x1.0c9714fbcda3bp-66}, {0x1.a3f1ef82137eep-10, 0x1.c71c71c71c71cp-65}, {0x1.90a9f747db95dp-10, -0x1.6a56a56a56a57p-65}, {0x1.7ed079ed4c037p-10, -0x1.03d226357e16fp-66}, {0x1.6e40790442038p-10, 0x1.a6f4de9bd37a7p-67}, {0x1.5ed9a0bd901b6p-10, -0x1.040e6c2b4481dp-64}, {0x1.507f94c2470bdp-10, 0x1.56070381c0e07p-65}, {0x1.43195bf54e5d7p-10, -0x1.b9e4b17e4b17ep-64}, {0x1.3690e51f04536p-10, 0x1.9d9745d1745d1p-64}, {0x1.2ad29fcd49d54p-10, 0x1.82fbf309b8b57p-64}, {0x1.1fcd25ae4a26ep-10, 0x1.ca730fcd6e9e0p-65}, {0x1.1570f16ece10ap-10, 0x1.a9e4bf3a9a378p-64}, {0x1.0bb020bc1eaa0p-10, 0x1.8c4b078787878p-64}, {0x1.027e3f7fcd8bfp-10, 0x1.67f632c234f73p-64}, {0x1.f3a03591c2915p-11, 0x1.22b2be05c0b81p-65}, {0x1.e337343f5c1f5p-11, 0x1.ef13d589d89d9p-65}, {0x1.d3af3c78ce2e4p-11, 0x1.8eec64a5294a5p-68}, {0x1.c4f7c88f08b5ep-11, 0x1.4d45526bca1afp-66}, {0x1.b701d9e1f038ep-11, -0x1.2bb23e09fab8cp-65}, {0x1.a9bfcd93a26fdp-11, 0x1.437ecd445d174p-67}, {0x1.9d2536c99619ap-11, -0x1.fd25b093f5dc8p-70}, {0x1.9126bdad380c3p-11, -0x1.b9945c5ec813ep-65}, {0x1.85ba027ef48d7p-11, -0x1.dd7be76e2e2bep-65}, {0x1.7ad5842a041b9p-11, 0x1.2880cf6e197edp-65}, {0x1.707089ddf7dc4p-11, 0x1.f99269e145605p-65}, {0x1.66830f45c9d3fp-11, -0x1.ad72816baa298p-68}};
static const dd PI = {0x1.921fb54442d18p+1, 0x1.1a62633145c07p-53}, PI_2 = {0x1.921fb54442d18p+0, 0x1.1a62633145c07p-54};

static dd two_sum(double a, double b)
{
    double s = a + b, v = s - a;
    return (dd){s, (a - (s - v)) + (b - v)};
}

static dd dd_add(dd a, dd b)
{
    dd s = two_sum(a.hi, b.hi);
    double lo = s.lo + a.lo + b.lo;
    double hi = s.hi + lo;
    return (dd){hi, lo - (hi - s.hi)};
}

static dd dd_mul(dd a, dd b)
{
    double p = a.hi * b.hi, e = fma(a.hi, b.hi, -p) + (a.hi * b.lo + a.lo * b.hi);
    double hi = p + e;
    return (dd){hi, e - (hi - p)};
}

static dd dd_series(const dd *c, int n, dd z)
{
    dd p = c[n - 1];
    for (int k = n - 2; k >= 0; k--) p = dd_add(dd_mul(p, z), c[k]);
    return p;
}

static double dd_round(dd a) { return a.hi + a.lo; }

/* 0 <= x <= pi/4 */
static double small_sin(double x)
{
    dd z = dd_mul((dd){x, 0.0}, (dd){x, 0.0});
    return dd_round(dd_mul(dd_series(SIN_C, 15, z), (dd){x, 0.0}));
}

static double small_cos(double x) { return dd_round(dd_series(COS_C, 15, dd_mul((dd){x, 0.0}, (dd){x, 0.0}))); }

static dd dd_asin(dd y)
{
    return dd_mul(dd_series(ASIN_C, 56, dd_mul(y, y)), y);
}

static dd dd_sqrt(double t)
{
    double s = sqrt(t), e = fma(-s, s, t);
    double lo = e / (2.0 * s), hi = s + lo;
    return (dd){hi, lo - (hi - s)};
}

static double exact_acos(double x)
{
    if (x == 1.0) return 0.0;
    if (x == -1.0) return dd_round(PI);
    if (!(fabs(x) <= 1.0)) return (x - x) / (x - x);
    if (fabs(x) <= 0.5) {
        dd a = dd_asin((dd){x, 0.0});
        return dd_round(dd_add(PI_2, (dd){-a.hi, -a.lo}));
    }
    dd a = dd_asin(dd_sqrt((x > 0 ? 1.0 - x : 1.0 + x) * 0.5));
    a = (dd){2.0 * a.hi, 2.0 * a.lo};
    return x > 0 ? dd_round(a) : dd_round(dd_add(PI, (dd){-a.hi, -a.lo}));
}

static double d2r(double x) { return fma(x, 1.7453283071517944e-02, x * 9.4484253514332993e-09); }

double geo_cos_deg(double x)
{
    if (!isfinite(x)) return x - x;
    x = fabs(x);
    if (x <= 0x1p-21) return 1.0;
    x = fmod(x, 360.0);
    int n = (int)x;
    if (n == x && n % 30 == 0) {
        if (n % 180 == 0) return n == 180 ? -1.0 : 1.0;
        if (n % 90 == 0) return 0.0;
        if (n % 60 == 0) return n == 60 || n == 300 ? 0.5 : -0.5;
        return n == 30 || n == 330 ? 8.6602540378443860e-01 : -8.6602540378443860e-01;
    }
    int neg = 0, fold_sin = 0;
    if (x <= 180) {
        if (x <= 90) {
            if (x > 45) x = 90 - x, fold_sin = 1;
        } else {
            if (x <= 135) x = x - 90, fold_sin = 1;
            else x = 180 - x;
            neg = 1;
        }
    } else if (x <= 270) {
        if (x <= 225) x = x - 180;
        else x = 270 - x, fold_sin = 1;
        neg = 1;
    } else {
        if (x <= 315) x = x - 270, fold_sin = 1;
        else x = 360 - x;
    }
    x = d2r(x);
    x = fold_sin ? small_sin(x) : small_cos(x);
    return neg ? -x : x;
}

double geo_sin_deg(double x)
{
    if (!isfinite(x)) return x - x;
    double s = copysign(1.0, x);
    if (fabs(x) < 0x1p-19) return d2r(x);
    x = fmod(fabs(x), 360.0);
    int n = (int)x;
    if (n == x && n % 30 == 0) {
        if (n % 180 == 0) {
            x = 0.0;
            if (n == 180) s = -s;
        } else if (n % 90 == 0) {
            x = n == 90 ? 1.0 : -1.0;
        } else if (n % 60 == 0) {
            x = n >= 180 ? -8.6602540378443860e-01 : 8.6602540378443860e-01;
        } else {
            x = n < 180 ? 0.5 : -0.5;
        }
        return x * s;
    }
    int fold_cos = 0;
    if (x <= 180) {
        if (x <= 90) {
            if (x > 45) x = 90 - x, fold_cos = 1;
        } else {
            if (x <= 135) x = x - 90, fold_cos = 1;
            else x = 180 - x;
        }
    } else if (x <= 270) {
        if (x <= 225) x = x - 180;
        else x = 270 - x, fold_cos = 1;
        s = -s;
    } else {
        if (x <= 315) x = x - 270, fold_cos = 1;
        else x = 360 - x;
        s = -s;
    }
    x = d2r(x);
    x = fold_cos ? small_cos(x) : small_sin(x);
    return x * s;
}

double geo_acos_upstream(double cosine) { return exact_acos(cosine) * 57.29577951308232; }

double geo_acos_deg(double cosine)
{
    if (cosine == 1.0) return 0.0;
    if (cosine == 0.5) return 60.0;
    if (cosine == 0.0) return 90.0;
    if (cosine == -0.5) return 120.0;
    if (cosine == -1.0) return 180.0;
    return geo_acos_upstream(cosine);
}
