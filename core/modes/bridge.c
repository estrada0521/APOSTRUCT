/* Batched Source bridge projection: coefficients, basis function,
 * project_vector_ and the final setting transform, one lane per displayed
 * atom, in the elementwise order of the Python batch these replaced; build
 * with -ffp-contract=off. */

#include <math.h>
#include <stdlib.h>
#include "pysum.h"

enum {
    BM_OK = 0,
    BM_MEMORY = 1,
    BM_CAPACITY = 2,
};

#define COEFFICIENT_LENGTH 181
#define BASIS_LENGTH 360
#define OUTPUT_LENGTH 720
#define FAMILY_STRIDE 144

/* project_vector_ over lanes: basis is lanes x basis_length, output lanes x
 * output_length (pre-zeroed). Returns the emitted count. */
static int project_vector_lanes(int lanes, const double *basis, int basis_length, int atom_count,
                                const int *reps, int target_rep, int vector_dim,
                                const double *constants, const double *points, int project_count,
                                double *output, int output_length)
{
    int emitted = 0;
    int step = vector_dim > 1 ? vector_dim : 1;
    for (int start = 1; start < 4; start += step) {
        int next = emitted;
        if (target_rep == reps[start - 1]) {
            if (project_count == 0) {
                next = emitted + 1;
            } else if (project_count > 0) {
                const double *constant = constants + (size_t)(start - 1) * vector_dim * 3;
                int family_offset = emitted * FAMILY_STRIDE;
                for (int family = 0; family < project_count; family++, family_offset += FAMILY_STRIDE) {
                    next++;
                    for (int atom = 0; atom < atom_count; atom++)
                        for (int lane = 0; lane < lanes; lane++) {
                            double local[3] = {0.0, 0.0, 0.0};
                            for (int column = 0; column < vector_dim; column++) {
                                int input = family * FAMILY_STRIDE + atom * 3 + column;
                                double value = input < basis_length ? basis[(size_t)lane * basis_length + input] : 0.0;
                                for (int axis = 0; axis < 3; axis++)
                                    local[axis] = local[axis] + constant[column * 3 + axis] * value;
                            }
                            const double *point = points + lane * 9;
                            for (int axis = 0; axis < 3; axis++) {
                                int index = family_offset + axis + atom * 3;
                                if (index >= 0 && index < output_length)
                                    output[(size_t)lane * output_length + index] = apostruct_builtin_sum3(
                                        point[axis] * local[0], point[3 + axis] * local[1],
                                        point[6 + axis] * local[2]);
                            }
                        }
                }
            }
        }
        emitted = next;
    }
    return emitted;
}

/* weights: lanes x dim x dim; out: lanes x 181. */
void bridge_coefficients(int lanes, int dim, const double *weights, int active_length,
                         const double *active, int vector_loop_count, int row_count,
                         double *out)
{
    for (int i = 0; i < lanes * COEFFICIENT_LENGTH; i++) out[i] = 0.0;
    for (int lane = 0; lane < lanes; lane++) {
        const double *w = weights + (size_t)lane * dim * dim;
        for (int row = 0; row < row_count; row++)
            for (int vector = 0; vector < vector_loop_count; vector++) {
                double value = 0.0;
                for (int column = 0; column < row_count; column++) {
                    int source = column * 3 + vector;
                    double source_value = source < active_length ? active[source] : 0.0;
                    double weight = row < dim && column < dim ? w[row * dim + column] : 0.0;
                    value = value + source_value * weight;
                }
                int output = row * 3 + vector;
                if (output < COEFFICIENT_LENGTH) out[lane * COEFFICIENT_LENGTH + output] = value;
            }
    }
}

/* One OPD group: rows are row_count x row_stride with per-row lengths;
 * constants[c] holds vector_dim x 3 basis constants for component start
 * c + 1 and reps[c] its vector representation; points: lanes x 9 (already
 * signed); transforms: transform_count 3x3 matrices applied in turn, then an
 * optional negation. modes: capacity x lanes x 3. Returns the mode count
 * through mode_count. */
int bridge_source_modes(int lanes, const double *coefficients, int row_count,
                        int row_stride, const int *row_lengths, const double *rows,
                        int vector_loop_count, int bridge_project_count,
                        const int *reps, int target_rep, int vector_dim,
                        const double *constants, const double *points, int project_count,
                        int split_row_slots, int transform_count, const double *transforms,
                        int negate, int capacity, double *modes, int *mode_count)
{
    double *basis = calloc((size_t)lanes * BASIS_LENGTH, sizeof(double));
    double *output = calloc((size_t)lanes * OUTPUT_LENGTH, sizeof(double));
    if (!basis || !output) {
        free(basis);
        free(output);
        return BM_MEMORY;
    }
    int term_count = bridge_project_count + 1;
    for (int lane = 0; lane < lanes; lane++) {
        const double *c = coefficients + (size_t)lane * COEFFICIENT_LENGTH;
        for (int row = 0; row < row_count; row++)
            for (int vector = 0; vector < vector_loop_count; vector++) {
                double value = 0.0;
                for (int term = 0; term < term_count; term++) {
                    double source = term < row_lengths[row] ? rows[row * row_stride + term] : 0.0;
                    int index = term * 3 + vector;
                    double coefficient = index < COEFFICIENT_LENGTH ? c[index] : 0.0;
                    value = value + source * coefficient;
                }
                int output_index = row * 3 + vector;
                if (output_index < BASIS_LENGTH) basis[lane * BASIS_LENGTH + output_index] = value;
            }
    }

    int emitted = project_vector_lanes(lanes, basis, BASIS_LENGTH, row_count, reps, target_rep,
                                       vector_dim, constants, points, project_count, output,
                                       OUTPUT_LENGTH);

    int count = 0;
    int per_emitted = split_row_slots ? row_count : 1;
    for (int e = 0; e < emitted; e++)
        for (int slot = 0; slot < per_emitted; slot++) {
            if (count >= capacity) {
                free(basis);
                free(output);
                return BM_CAPACITY;
            }
            int offset = e * FAMILY_STRIDE + (split_row_slots ? 3 * slot : 0);
            double *mode = modes + (size_t)count * lanes * 3;
            double scale = 0.0;
            int has_nan = 0;
            for (int lane = 0; lane < lanes; lane++) {
                double v[3];
                for (int axis = 0; axis < 3; axis++) {
                    int index = offset + axis;
                    v[axis] = index < OUTPUT_LENGTH ? output[lane * OUTPUT_LENGTH + index] : 0.0;
                }
                for (int t = 0; t < transform_count; t++) {
                    const double *m = transforms + 9 * t;
                    double r[3];
                    for (int column = 0; column < 3; column++)
                        r[column] = apostruct_builtin_sum3(v[0] * m[column], v[1] * m[3 + column], v[2] * m[6 + column]);
                    v[0] = r[0];
                    v[1] = r[1];
                    v[2] = r[2];
                }
                for (int axis = 0; axis < 3; axis++) {
                    double value = negate ? -v[axis] : v[axis];
                    mode[lane * 3 + axis] = value;
                    double magnitude = fabs(value);
                    if (isnan(magnitude)) has_nan = 1;
                    else if (magnitude > scale) scale = magnitude;
                }
            }
            if (!has_nan && scale > 1e-12 && fabs(scale - 1.0) > 1e-12)
                for (int i = 0; i < lanes * 3; i++) mode[i] = mode[i] / scale;
            count++;
        }
    *mode_count = count;
    free(basis);
    free(output);
    return BM_OK;
}

/* Source DISPLAY DISTORTION mode_b0 per occurrence lane: irreps lanes x dim
 * x dim (real), direction rows x dim, basis project_count x 144 flat, out
 * lanes x (project_count * 144). Both products sum upward from +0.0 without
 * fused products. */
void display_mode_b0(int lanes, int dim, const double *irreps, int direction_rows,
                     const double *direction, int project_count, const double *basis,
                     double *out)
{
    int length = project_count * FAMILY_STRIDE;
    double source[48], transformed[48];
    for (int i = 0; i < lanes * length; i++) out[i] = 0.0;
    for (int lane = 0; lane < lanes; lane++) {
        const double *irrep = irreps + (size_t)lane * dim * dim;
        for (int family = 0; family < project_count; family++)
            for (int component = 0; component < 3; component++) {
                for (int row = 0; row < dim; row++) source[row] = basis[family * FAMILY_STRIDE + row * 3 + component];
                for (int row = 0; row < dim; row++) {
                    double acc = 0.0;
                    for (int k = 0; k < dim; k++) acc = acc + irrep[row * dim + k] * source[k];
                    transformed[row] = acc;
                }
                for (int row = 0; row < direction_rows; row++) {
                    double acc = 0.0;
                    for (int k = 0; k < dim; k++) acc = acc + direction[row * dim + k] * transformed[k];
                    out[(size_t)lane * length + family * FAMILY_STRIDE + row * 3 + component] = acc;
                }
            }
    }
}

/* project_vector_ for one boundary state per lane (see project_vector_lanes). */
int project_vectors(int lanes, const double *basis, int basis_length, int atom_count,
                    const int *reps, int target_rep, int vector_dim, const double *constants,
                    const double *points, int project_count, double *output, int output_length)
{
    for (int i = 0; i < lanes * output_length; i++) output[i] = 0.0;
    return project_vector_lanes(lanes, basis, basis_length, atom_count, reps, target_rep,
                                vector_dim, constants, points, project_count, output, output_length);
}
