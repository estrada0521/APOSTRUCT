/* Small dense linear algebra in place of BLAS/LAPACK (math/linalg.c): the same
 * results on every platform. Matrices are row-major. */

#ifndef APOSTRUCT_LINALG_H
#define APOSTRUCT_LINALG_H

#include <complex.h>

/* C = op(A) op(B); op transposes when the flag is set; lda/ldb are row strides of the stored arrays. */
void la_dgemm(int trans_a, int trans_b, int m, int n, int k, const double *a, int lda, const double *b, int ldb, double *c, int ldc);
void la_zgemm(int m, int n, int k, const double complex *a, int lda, const double complex *b, int ldb, double complex *c, int ldc);
double la_ddot(int n, const double *x, const double *y);
/* inverse of an n x n matrix by LU with partial pivoting; 1 when singular */
int la_inverse(int n, const double *a, double *out);
/* singular values (descending) of an m x n matrix into s (min(m, n) values); 1 when out of memory */
int la_singular_values(int m, int n, const double *a, double *s);
/* minimum-norm least squares A X = B (A m x n, B m x nrhs, X n x nrhs), singular values at
 * most rcond * the largest treated as zero; *rank is the effective rank; 1 when out of memory */
int la_lstsq(int m, int n, int nrhs, const double *a, const double *b, double rcond, double *x, int *rank);

#endif
