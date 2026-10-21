/* Small dense linear algebra in place of BLAS/LAPACK. Every sum runs in index
 * order with fma from +0.0, so a result depends only on its inputs. */

#include "linalg.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef CMPLX
#define CMPLX(x, y) __builtin_complex((double)(x), (double)(y))
#endif

void la_dgemm(int trans_a, int trans_b, int m, int n, int k, const double *a, int lda, const double *b, int ldb, double *c, int ldc)
{
    for (int i = 0; i < m; i++) {
        double *row = c + (size_t)i * ldc;
        for (int j = 0; j < n; j++) row[j] = 0.0;
        for (int p = 0; p < k; p++) {
            double x = trans_a ? a[p * lda + i] : a[i * lda + p];
            if (trans_b)
                for (int j = 0; j < n; j++) row[j] = fma(x, b[j * ldb + p], row[j]);
            else
                for (int j = 0; j < n; j++) row[j] = fma(x, b[p * ldb + j], row[j]);
        }
    }
}

void la_zgemm(int m, int n, int k, const double complex *a, int lda, const double complex *b, int ldb, double complex *c, int ldc)
{
    double *re = malloc(sizeof(double) * (size_t)(2 * n + 1)), *im = re + n;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) re[j] = im[j] = 0.0;
        for (int p = 0; p < k; p++) {
            double ar = creal(a[i * lda + p]), ai = cimag(a[i * lda + p]);
            const double complex *brow = b + (size_t)p * ldb;
            for (int j = 0; j < n; j++) {
                double br = creal(brow[j]), bi = cimag(brow[j]);
                re[j] = fma(-ai, bi, fma(ar, br, re[j]));
                im[j] = fma(ai, br, fma(ar, bi, im[j]));
            }
        }
        for (int j = 0; j < n; j++) c[i * ldc + j] = CMPLX(re[j], im[j]);
    }
    free(re);
}

double la_ddot(int n, const double *x, const double *y)
{
    double acc = 0.0;
    for (int i = 0; i < n; i++) acc = fma(x[i], y[i], acc);
    return acc;
}

int la_inverse(int n, const double *a, double *out)
{
    double *lu = malloc(sizeof(double) * (size_t)(n * n + 1));
    int *pivot = malloc(sizeof(int) * (size_t)(n + 1));
    if (!lu || !pivot) {
        free(lu);
        free(pivot);
        return 2;
    }
    memcpy(lu, a, sizeof(double) * (size_t)(n * n));
    int singular = 0;
    for (int col = 0; col < n && !singular; col++) {
        int best = col;
        for (int r = col + 1; r < n; r++)
            if (fabs(lu[r * n + col]) > fabs(lu[best * n + col])) best = r;
        pivot[col] = best;
        if (lu[best * n + col] == 0.0) {
            singular = 1;
            break;
        }
        if (best != col)
            for (int j = 0; j < n; j++) {
                double t = lu[col * n + j];
                lu[col * n + j] = lu[best * n + j];
                lu[best * n + j] = t;
            }
        double inv = 1.0 / lu[col * n + col];
        for (int r = col + 1; r < n; r++) {
            double f = lu[r * n + col] * inv;
            lu[r * n + col] = f;
            for (int j = col + 1; j < n; j++) lu[r * n + j] = fma(-f, lu[col * n + j], lu[r * n + j]);
        }
    }
    if (!singular)
        for (int e = 0; e < n; e++) {
            double *x = out;
            for (int i = 0; i < n; i++) x[i * n + e] = i == e ? 1.0 : 0.0;
            for (int i = 0; i < n; i++)
                if (pivot[i] != i) {
                    double t = x[i * n + e];
                    x[i * n + e] = x[pivot[i] * n + e];
                    x[pivot[i] * n + e] = t;
                }
            for (int i = 0; i < n; i++) {
                double acc = x[i * n + e];
                for (int j = 0; j < i; j++) acc = fma(-lu[i * n + j], x[j * n + e], acc);
                x[i * n + e] = acc;
            }
            for (int i = n - 1; i >= 0; i--) {
                double acc = x[i * n + e];
                for (int j = i + 1; j < n; j++) acc = fma(-lu[i * n + j], x[j * n + e], acc);
                x[i * n + e] = acc / lu[i * n + i];
            }
        }
    free(lu);
    free(pivot);
    return singular;
}

/* One-sided Jacobi on the n columns of A (m x n): w = A V with orthogonal columns,
 * V accumulated when v is not NULL. Column-major work arrays. */
static void jacobi(int m, int n, double *w, double *v)
{
    for (int sweep = 0; sweep < 80; sweep++) {
        int rotated = 0;
        for (int p = 0; p < n - 1; p++)
            for (int q = p + 1; q < n; q++) {
                double *wp = w + (size_t)p * m, *wq = w + (size_t)q * m;
                double alpha = la_ddot(m, wp, wp), beta = la_ddot(m, wq, wq), gamma = la_ddot(m, wp, wq);
                if (gamma == 0.0 || fabs(gamma) <= 1e-15 * sqrt(alpha * beta)) continue;
                rotated = 1;
                double zeta = (beta - alpha) / (2.0 * gamma);
                double t = (zeta >= 0 ? 1.0 : -1.0) / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
                double c = 1.0 / sqrt(1.0 + t * t), s = c * t;
                for (int i = 0; i < m; i++) {
                    double x = wp[i], y = wq[i];
                    wp[i] = c * x - s * y;
                    wq[i] = s * x + c * y;
                }
                if (v) {
                    double *vp = v + (size_t)p * n, *vq = v + (size_t)q * n;
                    for (int i = 0; i < n; i++) {
                        double x = vp[i], y = vq[i];
                        vp[i] = c * x - s * y;
                        vq[i] = s * x + c * y;
                    }
                }
            }
        if (!rotated) return;
    }
}

static int descending(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a < b) - (a > b);
}

int la_singular_values(int m, int n, const double *a, double *s)
{
    int transpose = m < n, rows = transpose ? n : m, cols = transpose ? m : n;
    double *w = malloc(sizeof(double) * (size_t)rows * cols + 1), *norms = malloc(sizeof(double) * (size_t)cols + 1);
    if (!w || !norms) {
        free(w);
        free(norms);
        return 1;
    }
    for (int c = 0; c < cols; c++)
        for (int r = 0; r < rows; r++) w[(size_t)c * rows + r] = transpose ? a[c * n + r] : a[r * n + c];
    jacobi(rows, cols, w, NULL);
    for (int c = 0; c < cols; c++) norms[c] = sqrt(la_ddot(rows, w + (size_t)c * rows, w + (size_t)c * rows));
    qsort(norms, (size_t)cols, sizeof *norms, descending);
    memcpy(s, norms, sizeof(double) * (size_t)cols);
    free(w);
    free(norms);
    return 0;
}

int la_lstsq(int m, int n, int nrhs, const double *a, const double *b, double rcond, double *x, int *rank)
{
    double *w = malloc(sizeof(double) * (size_t)m * n + 1), *v = calloc((size_t)n * n + 1, sizeof(double));
    double *sigma = malloc(sizeof(double) * (size_t)n + 1);
    if (!w || !v || !sigma) {
        free(w);
        free(v);
        free(sigma);
        return 1;
    }
    for (int c = 0; c < n; c++)
        for (int r = 0; r < m; r++) w[(size_t)c * m + r] = a[r * n + c];
    for (int i = 0; i < n; i++) v[(size_t)i * n + i] = 1.0;
    jacobi(m, n, w, v);
    double largest = 0.0;
    for (int c = 0; c < n; c++) {
        sigma[c] = sqrt(la_ddot(m, w + (size_t)c * m, w + (size_t)c * m));
        if (sigma[c] > largest) largest = sigma[c];
    }
    double cut = rcond * largest;
    memset(x, 0, sizeof(double) * (size_t)n * nrhs);
    *rank = 0;
    for (int c = 0; c < n; c++) {
        if (!(sigma[c] > cut)) continue;
        (*rank)++;
        const double *wc = w + (size_t)c * m, *vc = v + (size_t)c * n;
        double s2 = sigma[c] * sigma[c];
        for (int j = 0; j < nrhs; j++) {
            double acc = 0.0;
            for (int r = 0; r < m; r++) acc = fma(wc[r], b[r * nrhs + j], acc);
            double f = acc / s2;
            for (int i = 0; i < n; i++) x[i * nrhs + j] = fma(vc[i], f, x[i * nrhs + j]);
        }
    }
    free(w);
    free(v);
    free(sigma);
    return 0;
}
