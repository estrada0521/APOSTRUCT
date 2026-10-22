/* Complete-mode Cartesian Gram-Schmidt within print families
 * Every arithmetic step
 * keeps the operation order of the former numpy implementation: Cartesian rows
 * use the CPython sum() replica, row dots add rows upward from +0.0. */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "pysum.h"

enum {
    GS_OK = 0,
    GS_SHAPE = 1,
    GS_MEMORY = 2,
};

static void gs_cartesian(const double *d, int64_t rows, const double *basis, double *out)
{
    for (int64_t r = 0; r < rows; r++)
        for (int col = 0; col < 3; col++)
            out[3 * r + col] = apostruct_builtin_sum3(d[3 * r] * basis[col], d[3 * r + 1] * basis[3 + col],
                                                      d[3 * r + 2] * basis[6 + col]);
}

static double gs_row_dot(const double *left, const double *right, int64_t rows)
{
    double total = 0.0;
    for (int64_t r = 0; r < rows; r++)
        total = total + apostruct_builtin_sum3(left[3 * r] * right[3 * r], left[3 * r + 1] * right[3 * r + 1],
                                               left[3 * r + 2] * right[3 * r + 2]);
    return total;
}

static void gs_support(const double *d, int64_t rows, unsigned char *support)
{
    for (int64_t r = 0; r < rows; r++)
        support[r] = d[3 * r] != 0.0 || d[3 * r + 1] != 0.0 || d[3 * r + 2] != 0.0;
}

/* d -= coefficient * ref; returns whether any bit of d changed */
static int subtract(double *d, double coefficient, const double *ref, int64_t rows)
{
    int changed = 0;
    for (int64_t e = 0; e < 3 * rows; e++) {
        double value = d[e] - coefficient * ref[e];
        changed |= memcmp(&value, &d[e], sizeof value) != 0;
        d[e] = value;
    }
    return changed;
}

/* groups[i] < 0: only the normfactor; otherwise orthogonalize against the
 * earlier modes of the same group.  dxyz rows offsets[i]..offsets[i+1]. */
int gram_schmidt_modes(int mode_count, const int *groups, const int64_t *offsets, double *dxyz,
                       const double *basis, int has_centering, double centering,
                       double *normfactor, int *has_normfactor)
{
    int64_t total_rows = offsets[mode_count];
    double *cartesian = malloc(sizeof(double) * 3 * (size_t)(total_rows ? total_rows : 1));
    double *denominators = malloc(sizeof(double) * (size_t)(mode_count ? mode_count : 1));
    unsigned char *support = malloc((size_t)(total_rows ? total_rows : 1));
    double *scratch = malloc(sizeof(double) * 3 * (size_t)(total_rows ? total_rows : 1));
    int *previous = malloc(sizeof(int) * (size_t)(mode_count ? mode_count : 1));
    int *last_in_group = NULL;
    int group_limit = 0;
    for (int i = 0; i < mode_count; i++)
        if (groups[i] + 1 > group_limit) group_limit = groups[i] + 1;
    last_in_group = malloc(sizeof(int) * (size_t)(group_limit ? group_limit : 1));
    int status = GS_OK;
    if (!cartesian || !denominators || !support || !scratch || !previous || !last_in_group) {
        status = GS_MEMORY;
        goto done;
    }
    for (int g = 0; g < group_limit; g++) last_in_group[g] = -1;
    for (int i = 0; i < mode_count; i++) {
        int64_t start = offsets[i], rows = offsets[i + 1] - start;
        double *d = dxyz + 3 * start;
        int group = groups[i];
        if (group >= 0) {
            unsigned char *own = support + start;
            gs_support(d, rows, own);
            int chain[4096], length = 0;
            for (int p = last_in_group[group]; p >= 0; p = previous[p]) {
                if (length == 4096) {
                    status = GS_SHAPE;
                    goto done;
                }
                chain[length++] = p;
            }
            int dirty = 1;
            for (int k = length - 1; k >= 0; k--) {
                int p = chain[k];
                int64_t pstart = offsets[p];
                if (offsets[p + 1] - pstart != rows) {
                    status = GS_SHAPE;
                    goto done;
                }
                const double *ref = dxyz + 3 * pstart;
                const unsigned char *ref_support = support + pstart;
                int overlap = 0;
                for (int64_t r = 0; r < rows && !overlap; r++) overlap = own[r] && ref_support[r];
                if (!overlap) {
                    dirty |= subtract(d, 0.0, ref, rows);
                    continue;
                }
                if (dirty) gs_cartesian(d, rows, basis, scratch);
                dirty = 0;
                double numerator = gs_row_dot(scratch, cartesian + 3 * pstart, rows);
                if (denominators[p] <= 1e-15) continue;
                double coefficient = numerator / denominators[p];
                dirty |= subtract(d, coefficient, ref, rows);
                if (coefficient != 0.0) gs_support(d, rows, own);
            }
            double scale = 0.0;
            int nan = 0;
            for (int64_t e = 0; e < 3 * rows; e++) {
                double value = fabs(d[e]);
                if (value != value) nan = 1;
                else if (value > scale) scale = value;
            }
            if (!nan && rows > 0 && scale > 1e-15)
                for (int64_t e = 0; e < 3 * rows; e++) d[e] = d[e] / scale;
            gs_support(d, rows, own);
            previous[i] = last_in_group[group];
            last_in_group[group] = i;
        }
        gs_cartesian(d, rows, basis, cartesian + 3 * start);
        double denominator = gs_row_dot(cartesian + 3 * start, cartesian + 3 * start, rows);
        denominators[i] = denominator;
        if (denominator <= 0) {
            has_normfactor[i] = 0;
            normfactor[i] = 0.0;
        } else {
            double value = 1.0 / sqrt(denominator);
            if (has_centering) value *= centering;
            has_normfactor[i] = 1;
            normfactor[i] = value;
        }
    }
done:
    free(cartesian);
    free(denominators);
    free(support);
    free(scratch);
    free(previous);
    free(last_in_group);
    return status;
}

int family_presentation_order(int count, const int64_t *modes, int *perm);

/* One (spec, parent site) block of mode definitions: Gram-Schmidt within
 * the print families keyed by (site pg, pg irrep, direction component),
 * numbered by first appearance, then the Source family presentation order
 * (k/family_order.c rows).  dxyz: mode_count x rows x 3, orthogonalized in
 * place; perm[target] = the mode the target presents. */
int present_block(int mode_count, int64_t rows, double *dxyz, const int64_t *keys, const int64_t *family_rows,
                  const double *basis, int has_centering, double centering, double *normfactor,
                  int *has_normfactor, int *perm)
{
    if (mode_count <= 0) return GS_OK;
    int *groups = malloc(sizeof(int) * (size_t)mode_count);
    int64_t *offsets = malloc(sizeof(int64_t) * ((size_t)mode_count + 1));
    if (!groups || !offsets) {
        free(groups);
        free(offsets);
        return GS_MEMORY;
    }
    int group_count = 0;
    for (int i = 0; i < mode_count; i++) {
        offsets[i] = (int64_t)i * rows;
        groups[i] = -1;
        for (int j = 0; j < i && groups[i] < 0; j++)
            if (keys[3 * j] == keys[3 * i] && keys[3 * j + 1] == keys[3 * i + 1] && keys[3 * j + 2] == keys[3 * i + 2])
                groups[i] = groups[j];
        if (groups[i] < 0) groups[i] = group_count++;
    }
    offsets[mode_count] = (int64_t)mode_count * rows;
    int status = gram_schmidt_modes(mode_count, groups, offsets, dxyz, basis, has_centering, centering, normfactor,
                                    has_normfactor);
    if (!status && family_presentation_order(mode_count, family_rows, perm)) status = GS_MEMORY;
    free(groups);
    free(offsets);
    return status;
}
