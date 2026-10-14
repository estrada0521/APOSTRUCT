/* subgroup_to_orderparam_ for a little irrep over a list of operations: the
 * D - I constraint rows (columns reversed), their xrowop2_ reduction and
 * the invariant directions from the free columns (rowops.py). */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "irreps.h"
#include "kernels.h"
#include "magnetic.h"

enum {
    SD_MEMORY = 1,
    SD_STACK = 2,
};

/* mode_case: 1 for get_irreps_matrix_for_case (params: three rationals as
 * (numerator, denominator)), 0 for the isotropy get_irreps_matrix_for_record
 * (params: the k parameter record).  point_mode: the records' point action
 * (mag_point_action).  directions_out: up to dim directions of dim
 * coordinates.  Returns the direction count or a negative status. */
int subduction_directions(int gid, int mode_case, int point_mode, int param_count, const int64_t *params,
                          int record_count, const int *records, double *directions_out)
{
    irrep_context c;
    const int64_t none[1] = {0};
    if (irrep_context_for(gid, 0, none, &c)) return -SD_STACK;
    int dim = c.dim, rows = record_count * dim;
    double *stack = malloc(sizeof(double) * (size_t)(rows + 1) * dim);
    double *m = malloc(sizeof(double) * (size_t)(rows + 1) * dim);
    int *mapped = malloc(sizeof(int) * (5 * (size_t)record_count + 1)), *signs = malloc(sizeof(int) * ((size_t)record_count + 1));
    if (!stack || !m || !mapped || !signs) {
        free(stack);
        free(m);
        free(mapped);
        free(signs);
        return -SD_MEMORY;
    }
    if (point_mode && mag_tables()) record_count = -1;
    for (int r = 0; r < record_count; r++) {
        memcpy(mapped + 5 * r, records + 5 * r, sizeof(int) * 4);
        mapped[5 * r + 4] = mag_point_action(point_mode, records[5 * r + 4], signs + r);
    }
    int status = record_count < 0 || (point_mode && mag_take_error());
    if (!status)
        status = mode_case ? case_irrep_stack(gid, param_count, params, record_count, mapped,
                                              (int64_t)record_count * dim * dim, stack)
                           : irrep_matrix_stack(gid, param_count, params, record_count, mapped,
                                                record_count * dim * dim, stack);
    if (status) {
        free(stack);
        free(m);
        free(mapped);
        free(signs);
        return -SD_STACK;
    }
    for (int r = 0; r < record_count; r++)
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                double entry = stack[((size_t)r * dim + i) * dim + (dim - 1 - j)];
                double value = (signs[r] < 0 ? -entry : entry) - ((dim - 1 - j) == i ? 1.0 : 0.0);
                m[((size_t)r * dim + i) * dim + j] = fabs(value) < 1e-12 ? 0.0 : value;
            }
    const double tol = 1e-10;
    int pivot_row = 0;
    for (int col = 0; col < dim && pivot_row < rows; col++) {
        int pivot = -1;
        for (int row = pivot_row; row < rows && pivot < 0; row++)
            if (fabs(m[(size_t)row * dim + col]) > tol) pivot = row;
        if (pivot < 0) continue;
        if (pivot != pivot_row)
            for (int k = 0; k < dim; k++) {
                double t = m[(size_t)pivot_row * dim + k];
                m[(size_t)pivot_row * dim + k] = m[(size_t)pivot * dim + k];
                m[(size_t)pivot * dim + k] = t;
            }
        double divisor = m[(size_t)pivot_row * dim + col];
        for (int k = 0; k < dim; k++) m[(size_t)pivot_row * dim + k] = m[(size_t)pivot_row * dim + k] / divisor;
        for (int row = 0; row < rows; row++) {
            if (row == pivot_row) continue;
            double factor = m[(size_t)row * dim + col];
            if (fabs(factor) > tol)
                for (int k = 0; k < dim; k++) {
                    double product = factor * m[(size_t)pivot_row * dim + k];
                    m[(size_t)row * dim + k] = m[(size_t)row * dim + k] - product;
                }
        }
        pivot_row++;
    }
    for (int i = 0; i < rows * dim; i++)
        if (fabs(m[i]) < tol) m[i] = 0.0;
    int nonzero_rows[48], rank = 0, pivots[48], pivot_count = 0;
    for (int row = 0; row < rows; row++) {
        double scale = 0.0;
        for (int k = 0; k < dim; k++)
            if (fabs(m[(size_t)row * dim + k]) > scale) scale = fabs(m[(size_t)row * dim + k]);
        if (!(scale > tol)) continue;
        if (rank < 48) nonzero_rows[rank] = row;
        rank++;
    }
    int count = dim - rank;
    for (int i = 0; i < rank && i < 48; i++) {
        const double *row = m + (size_t)nonzero_rows[i] * dim;
        for (int k = 0; k < dim; k++)
            if (fabs(row[k]) > tol) {
                pivots[pivot_count++] = k;
                break;
            }
    }
    if (count > 0) memset(directions_out, 0, sizeof(double) * (size_t)count * dim);
    int free_index = 0;
    for (int write_col = 0, matrix_col = dim - 1; matrix_col >= 0; write_col++, matrix_col--) {
        int is_pivot = 0;
        for (int p = 0; p < pivot_count && !is_pivot; p++) is_pivot = pivots[p] == matrix_col;
        if (is_pivot) continue;
        if (free_index >= count) break;
        double *direction = directions_out + (size_t)free_index * dim;
        direction[write_col] = 1.0;
        for (int p = 0; p < pivot_count; p++) {
            const double *row = m + (size_t)nonzero_rows[p] * dim;
            double pivot_value = row[pivots[p]];
            if (fabs(pivot_value) <= tol) continue;
            direction[dim - pivots[p] - 1] = -row[matrix_col] / pivot_value;
        }
        free_index++;
    }
    for (int i = 0; i < (count > 0 ? count : 0) * dim; i++)
        if (fabs(directions_out[i]) < tol) directions_out[i] = 0.0;
    free(stack);
    free(m);
    free(mapped);
    free(signs);
    return count;
}
