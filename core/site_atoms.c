/* The atom table of one parent Wyckoff site in a supercell: DISPLAY
 * DISTORTION's Point column, the operation record of every atom, and the
 * atom order rebased on the parent representative. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"

enum {
    SA_TABLES = 1,
    SA_RANGE = 2,
    SA_SINGULAR = 3,
    SA_TRANSPORT = 4,
    SA_OVERLAP = 5,
    SA_UNIQUE = 6,
    SA_PERMUTATION = 7,
    SA_MISSING = 8,
    SA_CAPACITY = 9,
};

static void record_values(const int64_t *record, rat *out)
{
    for (int axis = 0; axis < 3; axis++) out[axis] = rat_make(record[axis], record[3]);
}

/* A rational vector as (x, y, z, den) over its least common denominator. */
static void vector_record(const rat *v, int64_t *out)
{
    i128 den = 1;
    for (int axis = 0; axis < 3; axis++) den = rat_mul_i(den / rat_gcd(den, v[axis].d), v[axis].d);
    for (int axis = 0; axis < 3; axis++) {
        i128 n = rat_mul_i(v[axis].n, den / v[axis].d);
        if (n > INT64_MAX || n < INT64_MIN) rat_overflow = 1;
        out[axis] = (int64_t)n;
    }
    if (den > INT64_MAX) rat_overflow = 1;
    out[3] = (int64_t)den;
}

/* round() of a rational: nearest integer, ties to even. */
static i128 rat_round(rat a)
{
    i128 q = a.n / a.d, r = a.n % a.d;
    if (r < 0) {
        q -= 1;
        r += a.d;
    }
    i128 twice = 2 * r;
    if (twice > a.d || (twice == a.d && (q & 1))) q += 1;
    return q;
}


/* a * b as a 256-bit (high, low) pair. */
static void mul256(u128 a, u128 b, u128 *high, u128 *low)
{
    const u128 mask = UINT64_MAX;
    u128 a0 = a & mask, a1 = a >> 64, b0 = b & mask, b1 = b >> 64;
    u128 p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    u128 middle = (p00 >> 64) + (p01 & mask) + (p10 & mask);
    *low = (p00 & mask) | (middle << 64);
    *high = p11 + (p01 >> 64) + (p10 >> 64) + (middle >> 64);
}

/* |a| <= 1e-5 compared exactly against the double 1e-5 = m / 2^69. */
static int within_tolerance(rat a)
{
    const u128 mantissa = 0x14f8b588e368f1;
    u128 n = (u128)(a.n < 0 ? -a.n : a.n);
    u128 left_high = n >> 59, left_low = n << 69, right_high, right_low;
    mul256(mantissa, (u128)a.d, &right_high, &right_low);
    return left_high < right_high || (left_high == right_high && left_low <= right_low);
}

typedef struct {
    int op;
    rat t[3];
} keyed;

static int same_key(const keyed *a, const keyed *b)
{
    return a->op == b->op && rat_eq(a->t[0], b->t[0]) && rat_eq(a->t[1], b->t[1]) && rat_eq(a->t[2], b->t[2]);
}

static void key_of(int op, const rat *translation, const rat (*inverse)[3], keyed *out)
{
    out->op = op;
    rat_row_mul(translation, inverse, out->t);
    for (int axis = 0; axis < 3; axis++) out->t[axis] = rat_mod1(out->t[axis]);
}

static int rebased_order(int sg, int row_id, const rat *params, const int64_t *basis, const rat *representative,
                         const int64_t (*site)[5], int site_count, const int64_t (*cosets)[5], int coset_count,
                         const int64_t *kernel, int64_t kernel_count, const int64_t (*records)[5], int64_t count,
                         int64_t *order)
{
    int64_t generated[GEO_MAX_OPS][5];
    int generated_count = geo_sg_records(sg, generated);
    rat vectors[4][3], point[3];
    geo_wyckoff_vectors(row_id, vectors);
    for (int axis = 0; axis < 3; axis++) {
        point[axis] = vectors[0][axis];
        for (int j = 0; j < 3; j++) point[axis] = rat_add(point[axis], rat_mul(params[j], vectors[j + 1][axis]));
    }
    int64_t transport[5];
    int found = 0;
    for (int g = 0; g < generated_count && !found; g++) {
        rat rotated[3], t[3], shifted[3];
        geo_vrot(sg, (int)generated[g][4], point, rotated);
        record_values(generated[g], t);
        i128 integers[3];
        int close = 1;
        for (int axis = 0; axis < 3; axis++) {
            rat delta = rat_sub(rat_sub(representative[axis], rotated[axis]), t[axis]);
            integers[axis] = rat_round(delta);
            close &= within_tolerance(rat_sub(delta, rat_int(integers[axis])));
        }
        if (!close) continue;
        for (int axis = 0; axis < 3; axis++) shifted[axis] = rat_add(t[axis], rat_int(integers[axis]));
        vector_record(shifted, transport);
        transport[4] = generated[g][4];
        found = 1;
    }
    if (!found) return SA_TRANSPORT;
    if (!transport[0] && !transport[1] && !transport[2] && transport[3] == 1 && transport[4] == generated[0][4]) {
        for (int64_t i = 0; i < count; i++) order[i] = i;
        return 0;
    }
    int reduction_op[73] = {0};
    rat correction[73][3];
    for (int c = 0; c < coset_count; c++)
        for (int s = 0; s < site_count; s++) {
            int64_t point_op = geo_point_product((int)site[s][4], (int)cosets[c][4]);
            if (point_op < 1 || point_op > 72) return SA_TABLES;
            if (reduction_op[point_op]) return SA_OVERLAP;
            rat t[3];
            record_values(site[s], t);
            reduction_op[point_op] = (int)cosets[c][4];
            geo_vrot(sg, (int)cosets[c][4], t, correction[point_op]);
        }
    rat rows[3][3], inverse[3][3];
    for (int i = 0; i < 9; i++) rows[i / 3][i % 3] = rat_int(basis[i]);
    if (rat_inverse3((const rat (*)[3])rows, inverse)) return SA_SINGULAR;
    keyed *keys = malloc(sizeof(keyed) * (size_t)(count + 1));
    if (!keys) return SA_RANGE;
    for (int64_t i = 0; i < count; i++) {
        rat t[3];
        record_values(records[i], t);
        key_of((int)records[i][4], t, (const rat (*)[3])inverse, &keys[i]);
        for (int64_t j = 0; j < i; j++)
            if (same_key(&keys[i], &keys[j])) {
                free(keys);
                return SA_UNIQUE;
            }
    }
    int seen[73] = {0}, status = 0;
    int64_t ordered = 0;
    rat transported[3];
    record_values(transport, transported);
    for (int g = 0; g < generated_count && !status; g++) {
        int64_t product = geo_point_product((int)transport[4], (int)generated[g][4]);
        if (product < 1 || product > 72 || !reduction_op[product]) {
            status = SA_MISSING;
            break;
        }
        int point_op = reduction_op[product];
        if (seen[point_op]) continue;
        seen[point_op] = 1;
        rat rotated[3], t[3], translation[3];
        geo_vrot(sg, (int)generated[g][4], transported, rotated);
        record_values(generated[g], t);
        for (int axis = 0; axis < 3; axis++)
            translation[axis] = rat_sub(rat_add(rotated[axis], t[axis]), correction[product][axis]);
        for (int64_t e = 0; e < kernel_count && !status; e++) {
            rat shifted[3];
            for (int axis = 0; axis < 3; axis++)
                shifted[axis] = rat_add(translation[axis], rat_make(kernel[4 * e + axis], kernel[4 * e + 3]));
            keyed wanted;
            key_of(point_op, shifted, (const rat (*)[3])inverse, &wanted);
            int64_t hit = -1;
            for (int64_t i = 0; i < count && hit < 0; i++)
                if (same_key(&keys[i], &wanted)) hit = i;
            if (hit < 0 || ordered >= count) status = hit < 0 ? SA_MISSING : SA_PERMUTATION;
            else order[ordered++] = hit;
        }
    }
    free(keys);
    if (status) return status;
    if (ordered != count) return SA_PERMUTATION;
    for (int64_t i = 0; i < count; i++)
        for (int64_t j = 0; j < i; j++)
            if (order[i] == order[j]) return SA_PERMUTATION;
    return 0;
}

/* params and representative: three rationals, each (numerator, denominator)
 * as int128 words; basis: the supercell basis in pml, column convention.
 * points: count x 3 rationals in the same words; records: count x 5; order
 * (when representative is given): the rebased display permutation.  Returns
 * the atom count or a negative status. */
int64_t site_atom_table(int sg, int row_id, int site_pg, const int64_t *param_words, const int64_t *basis,
                        const int64_t *representative_words, int64_t capacity, int64_t *points, int64_t *records,
                        int64_t *order)
{
    rat_overflow = 0;
    if (geo_load()) return -SA_TABLES;
    geo_take_error();
    rat params[3];
    read_rat3(param_words, params);
    int64_t site[GEO_MAX_OPS][5], cosets[GEO_MAX_OPS][5];
    int site_count = geo_site_records(sg, row_id, site_pg, site);
    int coset_count = geo_cosets(sg, (const int64_t (*)[5])site, site_count, cosets);
    if (coset_count <= 0) return -SA_TABLES;
    int64_t *kernel = 0;
    int64_t kernel_count = geo_kernel_fractions(basis, &kernel);
    if (kernel_count <= 0) return kernel_count == -GEO_RANGE ? -SA_SINGULAR : -SA_TABLES;
    int64_t count = coset_count * kernel_count;
    if (count > capacity) {
        free(kernel);
        return -SA_CAPACITY;
    }
    rat vectors[4][3], to_cinter[3][3], cml_to_cinter[3][3], origin[3], cml_origin[3], to_cml[3][3];
    geo_wyckoff_vectors(row_id, vectors);
    geo_setting_to_cinter(sg, 1, 0, to_cinter, origin);
    geo_setting_to_cinter(sg, 0, 0, cml_to_cinter, cml_origin);
    int64_t record[36];
    int64_t ml_den = geo_ml_setting(geo_sg_lattice(sg), record);
    if (!ml_den) {
        free(kernel);
        return -SA_TABLES;
    }
    for (int i = 0; i < 9; i++) to_cml[i / 3][i % 3] = rat_make(record[9 + i], ml_den);
    int64_t index = 0;
    for (int c = 0; c < coset_count; c++) {
        int op = (int)cosets[c][4];
        rat branch[4][3], internal[3], cml[3], display[3], offset[3], base[3], t[3];
        for (int v = 0; v < 4; v++) geo_vrot(sg, op, vectors[v], branch[v]);
        for (int axis = 0; axis < 3; axis++) {
            internal[axis] = rat_int(0);
            for (int i = 0; i < 3; i++) internal[axis] = rat_add(internal[axis], rat_mul(params[i], branch[i + 1][axis]));
        }
        rat_row_mul(internal, (const rat (*)[3])to_cml, cml);
        rat_row_mul(cml, (const rat (*)[3])cml_to_cinter, display);
        for (int axis = 0; axis < 3; axis++) offset[axis] = rat_add(origin[axis], display[axis]);
        record_values(cosets[c], t);
        for (int axis = 0; axis < 3; axis++) base[axis] = rat_add(branch[0][axis], t[axis]);
        for (int64_t e = 0; e < kernel_count; e++, index++) {
            rat atom[3], shifted[3], point[3];
            for (int axis = 0; axis < 3; axis++) {
                rat lattice = rat_make(kernel[4 * e + axis], kernel[4 * e + 3]);
                atom[axis] = rat_add(lattice, t[axis]);
                shifted[axis] = rat_add(lattice, base[axis]);
            }
            vector_record(atom, records + 5 * index);
            records[5 * index + 4] = op;
            rat_row_mul(shifted, (const rat (*)[3])to_cinter, point);
            for (int axis = 0; axis < 3; axis++) {
                rat value = rat_add(point[axis], offset[axis]);
                write_i128(points + 12 * index + 4 * axis, value.n);
                write_i128(points + 12 * index + 4 * axis + 2, value.d);
            }
        }
    }
    int status = 0;
    if (representative_words) {
        rat representative[3];
        read_rat3(representative_words, representative);
        status = rebased_order(sg, row_id, params, basis, representative, (const int64_t (*)[5])site, site_count,
                               (const int64_t (*)[5])cosets, coset_count, kernel, kernel_count,
                               (const int64_t (*)[5])records, count, order);
    }
    free(kernel);
    if (!status && geo_take_error()) status = SA_TABLES;
    if (!status && rat_overflow) status = SA_RANGE;
    return status ? -status : count;
}
