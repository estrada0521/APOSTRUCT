/* Canonical child-atom vectors for every mode of one row plan:
 * optional presentation inverse (CPython-sum replica per axis) with max-abs
 * scaling, gather onto canonical atoms, Source-equivalence and finiteness
 * checks, and per-atom negation. */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

double apostruct_builtin_sum3(double a, double b, double c);

enum {
    MR_OK = 0,
    MR_EQUIVALENT = 1,
    MR_NONFINITE = 2,
};

/* pre (when has_pre): the presentation inverse applied to the raw vectors
 * first, unscaled; inverse (when has_inverse): applied on the plan with the
 * per-mode scale.  Only one of them applies to a block. */
int mode_rows_from_plan(int mode_count, int raw_count, const double *raw, int has_pre, const double *pre,
                        int has_inverse, const double *inverse, int atom_count, const int64_t *raw_indices,
                        const int64_t *equivalent_offsets, const int64_t *equivalents,
                        const unsigned char *negated, double *scratch, double *out)
{
    if (has_pre && has_inverse) return MR_EQUIVALENT;
    for (int m = 0; m < mode_count; m++) {
        const double *source = raw + (int64_t)m * raw_count * 3;
        if (has_pre) {
            for (int r = 0; r < raw_count; r++)
                for (int axis = 0; axis < 3; axis++)
                    scratch[3 * r + axis] = apostruct_builtin_sum3(source[3 * r] * pre[axis], source[3 * r + 1] * pre[3 + axis],
                                                                   source[3 * r + 2] * pre[6 + axis]);
            source = scratch;
        }
        double *dxyz = out + (int64_t)m * atom_count * 3;
        if (has_inverse) {
            double scale = 0.0;
            int nan = 0;
            for (int r = 0; r < raw_count; r++)
                for (int axis = 0; axis < 3; axis++) {
                    double value = apostruct_builtin_sum3(source[3 * r] * inverse[axis],
                                                          source[3 * r + 1] * inverse[3 + axis],
                                                          source[3 * r + 2] * inverse[6 + axis]);
                    scratch[3 * r + axis] = value;
                    double size = fabs(value);
                    if (size != size) nan = 1;
                    else if (size > scale) scale = size;
                }
            if (!nan && scale > 0.0)
                for (int e = 0; e < 3 * raw_count; e++) scratch[e] = scratch[e] / scale;
            for (int a = 0; a < atom_count; a++)
                for (int axis = 0; axis < 3; axis++) dxyz[3 * a + axis] = scratch[3 * raw_indices[a] + axis];
            for (int a = 0; a < atom_count; a++)
                for (int64_t k = equivalent_offsets[a]; k < equivalent_offsets[a + 1]; k++)
                    for (int axis = 0; axis < 3; axis++)
                        if (fabs(dxyz[3 * a + axis] - scratch[3 * equivalents[k] + axis]) > 1e-8)
                            return MR_EQUIVALENT;
        } else {
            for (int a = 0; a < atom_count; a++)
                for (int axis = 0; axis < 3; axis++) dxyz[3 * a + axis] = source[3 * raw_indices[a] + axis];
        }
        for (int e = 0; e < 3 * atom_count; e++)
            if (!isfinite(dxyz[e])) return MR_NONFINITE;
        if (negated)
            for (int a = 0; a < atom_count; a++)
                if (negated[a])
                    for (int axis = 0; axis < 3; axis++) dxyz[3 * a + axis] = -dxyz[3 * a + axis];
    }
    return MR_OK;
}

/* int(Decimal(format(value, ".4f")) * 10000) for every value: the
 * correctly rounded four-decimal value, as an integer count of 1e-4. */
int quantize4(int count, const double *values, int64_t *out)
{
    char text[64];
    for (int i = 0; i < count; i++) {
        int n = snprintf(text, sizeof text, "%.4f", values[i]);
        if (n <= 0 || n >= (int)sizeof text) return 1;
        int64_t integer = 0;
        int negative = 0;
        for (const char *c = text; *c; c++) {
            if (*c == '-') negative = 1;
            else if (*c >= '0' && *c <= '9') integer = integer * 10 + (*c - '0');
        }
        out[i] = negative ? -integer : integer;
    }
    return 0;
}
