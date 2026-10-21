/* transform_irrep2complex_: the type-2/3 real-to-complex transform of a
 * little irrep and its conjugate transpose, with numpy's complex products
 * (zgemm on the active block) and its small-entry clearing. */

#include "linalg.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "tables.h"

enum { LD = 48 };

static struct {
    section_view irr_sg, irr_k, star_count, minusk, full_dim, type, real2_pointer, real2, const_values;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (section_find("irreps.little.sg", &T.irr_sg) || section_find("irreps.little.k", &T.irr_k) ||
        section_find("k.star.arm_count", &T.star_count) || section_find("k.star.minus_k", &T.minusk) ||
        section_find("irreps.little.dim", &T.full_dim) || section_find("irreps.little.type", &T.type) ||
        section_find("irreps.little.real_form_start", &T.real2_pointer) || section_find("irreps.little.real_form", &T.real2) ||
        section_find("const.values", &T.const_values))
        return 1;
    T.ready = 1;
    return 0;
}

static double const_value(int64_t code);

static void zmatmlt(const double complex *left, double complex *right, int active)
{
    double complex *out = calloc((size_t)active * active + 1, sizeof(double complex));
    la_zgemm(active, active, active, left, LD, right, LD, out, active);
    for (int r = 0; r < active; r++)
        for (int c = 0; c < active; c++) right[r * LD + c] = out[r * active + c];
    free(out);
    for (int r = 0; r < LD; r++)
        for (int c = 0; c < LD; c++) {
            if (r >= active || c >= active) right[r * LD + c] = 0;
            if (cabs(right[r * LD + c]) < 1e-12) right[r * LD + c] = 0;
        }
}

/* The inverse transform's active dim x dim block (row-major) of gid:
 * 1 for type 1 (no transform), 0 with inverse, -1 for an invalid block
 * count, -2 for missing tables. */
int irrep_transform_inverse(int gid, double complex *inverse)
{
    if (load()) return -2;
    int type = (int)section_int(&T.type, gid - 1);
    if (type != 2 && type != 3) return 1;
    int active = (int)section_int(&T.full_dim, gid - 1);
    int sg = (int)section_int(&T.irr_sg, gid - 1), kslot = (int)section_int(&T.irr_k, gid - 1);
    int block_count = (int)section_int(&T.star_count, (int64_t)(sg - 1) * 27 + kslot - 1);
    if (block_count <= 0) return -1;
    int permutation = (int)section_int(&T.minusk, (int64_t)(sg - 1) * 27 + kslot - 1) != 0;
    int64_t real2 = section_int(&T.real2_pointer, gid - 1);
    int block_dim = active / block_count, half = block_dim / 2;
    double root2 = sqrt(2.0);
    double complex *matrix = calloc(LD * LD, sizeof(double complex));
    if (!matrix) return -2;
    for (int i = 0; i < active; i++) matrix[i * LD + i] = 1;
    if (!real2 && permutation && block_dim == 4 && block_count % 2 == 1) {
        memset(matrix, 0, sizeof(double complex) * LD * LD);
        int pairs = active / 2;
        for (int pair = 0; pair < pairs; pair++) {
            int top = pair, bottom = pair + pairs, col0 = (pair / 2) * 4 + (pair % 2), col1 = col0 + 2;
            if (pair % 2 == 0) {
                matrix[top * LD + col0] = 1.0 / root2;
                matrix[bottom * LD + col0] = 1.0 / root2;
                matrix[top * LD + col1] = (-1.0 * I) / root2;
                matrix[bottom * LD + col1] = (1.0 * I) / root2;
            } else {
                matrix[top * LD + col0] = (-1.0 * I) / root2;
                matrix[bottom * LD + col0] = (1.0 * I) / root2;
                matrix[top * LD + col1] = 1.0 / root2;
                matrix[bottom * LD + col1] = 1.0 / root2;
            }
        }
    } else {
        if (real2 && permutation) {
            static const int perm[4][4] = {{1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 1, 0, 0}};
            memset(matrix, 0, sizeof(double complex) * LD * LD);
            if (block_dim > 3) {
                int repeat = block_dim / 4;
                for (int block = 0; block < block_count; block++) {
                    int base = block * block_dim;
                    for (int outer = 0; outer < 4; outer++) {
                        int row_start = base + (outer * block_dim) / 4;
                        for (int inner = 0; inner < 4; inner++) {
                            int col_start = base + (inner * block_dim) / 4;
                            for (int offset = 0; offset < repeat; offset++)
                                matrix[(row_start + offset) * LD + col_start + offset] = perm[outer][inner];
                        }
                    }
                }
            }
        }
        double complex *first = calloc(LD * LD, sizeof(double complex));
        if (!first) return -2;
        const double complex values[4] = {1.0, -1.0 * I, 1.0, 1.0 * I};
        int base = 0;
        for (int block = 0; block < block_count; block++) {
            for (int outer = 0; outer < 2; outer++) {
                int row_start = base + (outer * block_dim) / 2;
                for (int inner = 0; inner < 2; inner++) {
                    int col_start = base + (inner * block_dim) / 2;
                    double complex value = values[outer * 2 + inner] / root2;
                    for (int offset = 0; offset < half; offset++) first[(row_start + offset) * LD + col_start + offset] = value;
                }
            }
            base += block_dim;
        }
        zmatmlt(first, matrix, active);
        memset(first, 0, sizeof(double complex) * LD * LD);
        int wrapped = 1;
        base = 0;
        for (int step = 0; step < block_count * 2; step++) {
            int row_start = base / 2, col_start = ((wrapped - 1) * block_dim) / 2;
            for (int offset = 0; offset < half; offset++) first[(row_start + offset) * LD + col_start + offset] = 1.0;
            wrapped += 2;
            if (block_count * 2 < wrapped) wrapped = 2;
            base += block_dim;
        }
        zmatmlt(first, matrix, active);
        if (type == 2 && real2) {
            memset(first, 0, sizeof(double complex) * LD * LD);
            int final_half = active / 2;
            for (int offset = 0; offset < final_half; offset++) first[offset * LD + offset] = 1.0;
            for (int column = 0; column < active / 2; column++) {
                int64_t at = (real2 - 1) * 3 + (int64_t)column * 3;
                int row = final_half + (int)section_int(&T.real2, at) - 1, col = final_half + column;
                first[row * LD + col] = const_value(section_int(&T.real2, at + 1)) + const_value(section_int(&T.real2, at + 2)) * I;
            }
            zmatmlt(first, matrix, active);
        }
        free(first);
    }
    for (int r = 0; r < active; r++)
        for (int c = 0; c < active; c++) inverse[r * active + c] = conj(matrix[c * LD + r]);
    free(matrix);
    return 0;
}

static double const_value(int64_t code) { return ((const double *)T.const_values.data)[code]; }
