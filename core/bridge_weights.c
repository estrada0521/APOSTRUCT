/* Phase-adjusted little-irrep matrices for many operation records at once.
 *
 * W(record) = P(phases) @ S(position), matching irreps.py: phases are exact
 * rationals k_arm . delta, P is the real star-arm phase operator and S the
 * decoded Source matrix at the record's generated position. The product sums
 * k upward from +0.0 without fused products (representation.phase_product).
 * Build with -ffp-contract=off. */

#include <math.h>
#include <stdint.h>

#include "cyclo.h"
#include "intmath.h"
#include <string.h>

enum {
    BW_OK = 0,
    BW_NO_POSITION = 1,
    BW_ARM_COUNT = 2,
    BW_SCALAR_PHASE = 3,
    BW_ODD_BLOCK_PHASE = 4,
    BW_OVERFLOW = 5,
};


static i128 gcd128(i128 a, i128 b) { return gcd_i128(a, b); }

static int64_t mod_floor(int64_t a, int64_t m)
{
    int64_t r = a % m;
    return r < 0 ? r + m : r;
}

/* Reduced (num, den) of value/den mod 1 for each axis, for position matching. */
static void tau_key(const int *record, int64_t *key)
{
    for (int axis = 0; axis < 3; axis++) {
        int64_t den = record[3];
        int64_t num = record[axis];
        if (den < 0) {
            den = -den;
            num = -num;
        }
        num = mod_floor(num, den);
        int64_t g = (int64_t)gcd128(num, den);
        if (g == 0) g = 1;
        key[2 * axis] = num / g;
        key[2 * axis + 1] = den / g;
    }
}

static int find_position(int element_count, const int *elements, const int *record)
{
    int64_t want[6], have[6];
    tau_key(record, want);
    for (int index = 0; index < element_count; index++) {
        const int *element = elements + 5 * index;
        if (element[4] != record[4]) continue;
        tau_key(element, have);
        if (!memcmp(want, have, sizeof want)) return index + 1;
    }
    return 0;
}

static int first_with_op(int element_count, const int *elements, int op)
{
    for (int index = 0; index < element_count; index++)
        if (elements[5 * index + 4] == op) return index;
    return -1;
}

/* phase = sum_i (kn_i / kd_i) * (dn_i / dd_i) exactly, reduced mod 1 to
 * (num, den) (both kept below 2^53). */
static int arm_phase(const int64_t *k, const i128 *dn, const i128 *dd, int64_t *out)
{
    i128 num = 0, den = 1;
    for (int axis = 0; axis < 3; axis++) {
        i128 tn = (i128)k[2 * axis] * dn[axis];
        i128 td = (i128)k[2 * axis + 1] * dd[axis];
        i128 g = gcd128(tn, td);
        if (g > 1) {
            tn /= g;
            td /= g;
        }
        g = gcd128(den, td);
        i128 scale_left = td / g, scale_right = den / g;
        num = num * scale_left + tn * scale_right;
        den = den * scale_left;
        g = gcd128(num, den);
        if (g > 1) {
            num /= g;
            den /= g;
        }
        if (den > ((i128)1 << 100) || num > ((i128)1 << 100) || num < -((i128)1 << 100))
            return BW_OVERFLOW;
    }
    const i128 limit = (i128)1 << 53;
    if (den > limit || num > limit || num < -limit) return BW_OVERFLOW;
    num %= den;
    if (num < 0) num += den;
    out[0] = (int64_t)num;
    out[1] = (int64_t)den;
    return BW_OK;
}

/* cos and sin of 2 pi n / d, correctly rounded from Q(zeta_d); libm only when
 * that field is out of reach */
static void phase_cos_sin(const int64_t *phase, double *c, double *s)
{
    typedef struct {
        int64_t n, d;
        double c, s;
    } entry;
    static entry memo[512];
    static int used;
    for (int i = 0; i < used; i++)
        if (memo[i].n == phase[0] && memo[i].d == phase[1]) {
            *c = memo[i].c;
            *s = memo[i].s;
            return;
        }
    const cy_field *f = phase[1] <= 100000 ? cy_field_for((int)phase[1]) : 0;
    cy exact_c, exact_s;
    cy_overflow = 0;
    if (!f || cy_cos_sin(f, phase[0], phase[1], &exact_c, &exact_s) || cy_overflow || cy_to_double(f, &exact_c, c) ||
        cy_to_double(f, &exact_s, s)) {
        double angle = 2.0 * M_PI * ((double)phase[0] / (double)phase[1]);
        *c = cos(angle);
        *s = sin(angle);
        if (fabs(*c) < 1e-12) *c = 0.0;
        if (fabs(*s) < 1e-12) *s = 0.0;
    }
    if (*c == 0.0) *c = 0.0;
    if (*s == 0.0) *s = 0.0;
    if (used < 512) memo[used++] = (entry){phase[0], phase[1], *c, *s};
}

static int phase_operator(int dim, int arms, const int64_t *phases, double *p)
{
    for (int i = 0; i < dim * dim; i++) p[i] = 0.0;
    if (arms == 0) {
        for (int i = 0; i < dim; i++) p[i * dim + i] = 1.0;
        return BW_OK;
    }
    if (dim % arms) return BW_ARM_COUNT;
    int block = dim / arms;
    for (int arm = 0; arm < arms; arm++) {
        int start = arm * block;
        double c, s;
        phase_cos_sin(phases + 2 * arm, &c, &s);
        if (block == 1) {
            if (s != 0.0) return BW_SCALAR_PHASE;
            p[start * dim + start] = c;
            continue;
        }
        if (block % 2) {
            if (s != 0.0) return BW_ODD_BLOCK_PHASE;
            for (int a = 0; a < block; a++)
                for (int b = 0; b < block; b++)
                    p[(start + a) * dim + start + b] = a == b ? c : c * 0.0;
            continue;
        }
        int half = block / 2;
        for (int a = 0; a < half; a++)
            for (int b = 0; b < half; b++) {
                int same = a == b;
                p[(start + a) * dim + start + b] = same ? c : c * 0.0;
                p[(start + a) * dim + start + half + b] = same ? s : s * 0.0;
                p[(start + half + a) * dim + start + b] = same ? -s : -s * 0.0;
                p[(start + half + a) * dim + start + half + b] = same ? c : c * 0.0;
            }
    }
    return BW_OK;
}

/* W(record) at a known generated position into w (dim x dim). */
int bw_record_matrix(int dim, int arms, const int64_t *k, int subtract_generated,
                     int element_count, const int *elements, int position,
                     const double *sparse, const int *record, double *w)
{
    int64_t phases[96];
    double p[48 * 48];
    if (dim <= 0 || dim > 48 || arms < 0 || arms > 48) return BW_ARM_COUNT;
    i128 dn[3], dd[3];
    int generated = -1;
    if (subtract_generated) {
        generated = first_with_op(element_count, elements, record[4]);
        if (generated < 0) return BW_NO_POSITION;
    }
    for (int axis = 0; axis < 3; axis++) {
        if (generated < 0) {
            dn[axis] = record[axis];
            dd[axis] = record[3];
        } else {
            const int *element = elements + 5 * generated;
            dn[axis] = (i128)record[axis] * element[3] - (i128)element[axis] * record[3];
            dd[axis] = (i128)record[3] * element[3];
        }
        if (dd[axis] < 0) {
            dd[axis] = -dd[axis];
            dn[axis] = -dn[axis];
        }
    }
    for (int arm = 0; arm < arms; arm++) {
        int status = arm_phase(k + 6 * arm, dn, dd, phases + 2 * arm);
        if (status) return status;
    }
    int status = phase_operator(dim, arms, phases, p);
    if (status) return status;
    /* Zero entries of p contribute only signed zeros, which never change a
     * sum started from +0.0, so the ascending sum over the nonzero ones is
     * the dense product bit for bit. */
    const double *s = sparse + (size_t)(position - 1) * dim * dim;
    for (int i = 0; i < dim; i++) {
        int columns[48], count = 0;
        for (int q = 0; q < dim; q++)
            if (p[i * dim + q] != 0.0) columns[count++] = q;
        for (int j = 0; j < dim; j++) {
            double acc = 0.0;
            for (int t = 0; t < count; t++) acc = acc + p[i * dim + columns[t]] * s[columns[t] * dim + j];
            w[i * dim + j] = acc;
        }
    }
    return BW_OK;
}

/* bw_record_matrix in two steps: the per-arm phase coefficients (c, s) after
 * its zero snapping, with the status it would return, and then any one row
 * of W from them, bit for bit as the dense product. */
int bw_record_phase_coefficients(int dim, int arms, const int64_t *k, int subtract_generated, int element_count,
                                 const int *elements, const int *record, double *coefficients)
{
    int64_t phases[96];
    if (dim <= 0 || dim > 48 || arms < 0 || arms > 48) return BW_ARM_COUNT;
    i128 dn[3], dd[3];
    int generated = -1;
    if (subtract_generated) {
        generated = first_with_op(element_count, elements, record[4]);
        if (generated < 0) return BW_NO_POSITION;
    }
    for (int axis = 0; axis < 3; axis++) {
        if (generated < 0) {
            dn[axis] = record[axis];
            dd[axis] = record[3];
        } else {
            const int *element = elements + 5 * generated;
            dn[axis] = (i128)record[axis] * element[3] - (i128)element[axis] * record[3];
            dd[axis] = (i128)record[3] * element[3];
        }
        if (dd[axis] < 0) {
            dd[axis] = -dd[axis];
            dn[axis] = -dn[axis];
        }
    }
    for (int arm = 0; arm < arms; arm++) {
        int status = arm_phase(k + 6 * arm, dn, dd, phases + 2 * arm);
        if (status) return status;
    }
    if (arms == 0) return BW_OK;
    if (dim % arms) return BW_ARM_COUNT;
    int block = dim / arms;
    for (int arm = 0; arm < arms; arm++) {
        double c, sn;
        phase_cos_sin(phases + 2 * arm, &c, &sn);
        if (block == 1 && sn != 0.0) return BW_SCALAR_PHASE;
        if (block > 1 && block % 2 && sn != 0.0) return BW_ODD_BLOCK_PHASE;
        coefficients[2 * arm] = c;
        coefficients[2 * arm + 1] = sn;
    }
    return BW_OK;
}

void bw_record_row(int dim, int arms, const double *coefficients, const double *sparse, int i, double *out)
{
    int columns[2], count = 0;
    double values[2];
    if (arms == 0) {
        columns[count] = i;
        values[count++] = 1.0;
    } else {
        int block = dim / arms, arm = i / block, a = i % block, start = arm * block;
        double c = coefficients[2 * arm], sn = coefficients[2 * arm + 1];
        if (block == 1 || block % 2) {
            if (c != 0.0) {
                columns[count] = i;
                values[count++] = c;
            }
        } else {
            int half = block / 2, base = a < half ? a : a - half;
            double first = a < half ? c : -sn, second = a < half ? sn : c;
            if (first != 0.0) {
                columns[count] = start + base;
                values[count++] = first;
            }
            if (second != 0.0) {
                columns[count] = start + half + base;
                values[count++] = second;
            }
        }
    }
    for (int j = 0; j < dim; j++) {
        double acc = 0.0;
        for (int t = 0; t < count; t++) acc = acc + values[t] * sparse[columns[t] * dim + j];
        out[j] = acc;
    }
}

int bw_find_position(int element_count, const int *elements, const int *record)
{
    return find_position(element_count, elements, record);
}

/* k: arms x 3 (num, den); elements: element_count x 5 generated records;
 * canonical_position[op] is the fallback position for a record whose
 * translation is not generated (0: none); sparse: element_count x dim x dim,
 * usable where available[position - 1]; records: record_count x 5;
 * out: record_count x dim x dim. */
int bridge_weight_matrices(int dim, int arms, const int64_t *k, int subtract_generated,
                           int element_count, const int *elements,
                           int op_limit, const int *canonical_position,
                           const double *sparse, const int *available,
                           int record_count, const int *records, double *out)
{
    if (dim <= 0 || dim > 48 || arms < 0 || arms > 48) return BW_ARM_COUNT;
    for (int r = 0; r < record_count; r++) {
        const int *record = records + 5 * r;
        int position = find_position(element_count, elements, record);
        if (!position && record[4] >= 0 && record[4] <= op_limit)
            position = canonical_position[record[4]];
        if (!position || !available[position - 1]) return BW_NO_POSITION;
        int status = bw_record_matrix(dim, arms, k, subtract_generated, element_count, elements,
                                      position, sparse, record, out + (size_t)r * dim * dim);
        if (status) return status;
    }
    return BW_OK;
}
