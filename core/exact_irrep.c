/* Exact little-irrep matrices W(record) = P(phases) S(position) over
 * Q(zeta_m): the Source entries from const.exact, the star-arm phases as
 * rationals k_arm . delta. */

#include <stdlib.h>
#include <string.h>

#include "cyclo.h"
#include "source_blob.h"
#include "tables.h"

int bw_find_position(int element_count, const int *elements, const int *record);

/* The Source matrix of a little irrep at a generated position as const.exact
 * words (8 per entry, a zero denominator where the entry has no exact form). */
static int little_exact_words(int gid, int dim, int position, int64_t *words)
{
    static section_view irr_pointer, pointer, matrices, exact;
    static int ready;
    if (!ready) {
        if (section_find("little.little_irr_full_matrices_irr_pointer", &irr_pointer)
            || section_find("little.little_irr_full_matrices_pointer", &pointer)
            || section_find("little.little_irr_full_matrices", &matrices) || section_find("const.exact", &exact)
            || !section_is_int(&exact))
            return 1;
        ready = 1;
    }
    int64_t irr_length = (int64_t)irr_pointer.entry->count, pointer_length = (int64_t)pointer.entry->count;
    if (gid < 1 || gid > irr_length || position < 1) return 1;
    int64_t index = section_int(&irr_pointer, gid - 1) + position - 1;
    if (index < 1 || index >= pointer_length) return 1;
    int64_t start = section_int(&pointer, index - 1), end = section_int(&pointer, index);
    memset(words, 0, sizeof(int64_t) * 8 * (size_t)dim * dim);
    for (int i = 0; i < dim * dim; i++)
        for (int j = 0; j < 4; j++) words[8 * i + 2 * j + 1] = 1;
    int64_t codes = (int64_t)exact.entry->count / 8;
    for (int64_t offset = (start - 1) * 3; offset < (end - 1) * 3; offset += 3) {
        int64_t row = section_int(&matrices, offset), col = section_int(&matrices, offset + 1);
        int64_t code = section_int(&matrices, offset + 2);
        if (row < 1 || row > dim || col < 1 || col > dim || code < 0 || code >= codes) return 1;
        int64_t *w = words + 8 * ((row - 1) * dim + col - 1);
        for (int j = 0; j < 8; j++) w[j] = section_int(&exact, 8 * code + j);
        if (!w[1] && code > 1) return 2;
        if (!w[1]) w[1] = w[3] = w[5] = w[7] = 1;
    }
    return 0;
}

static int64_t gcd64(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        int64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/* 0 ok, 1 tables, 2 no exact constant, 3 field too large or overflow, 4 shape */
/* The arm phases k_arm . delta of a record (pn / pd per arm, 0 <= pn < pd) and
 * the least m = lcm(24, pd); 0, or 3 (range) or 4 (shape). */
static int record_phases(int dim, int arms, const int64_t *k, int subtract_generated, int element_count, const int *elements,
                         const int *record, int64_t *pn, int64_t *pd, int64_t *m_out)
{
    if (dim <= 0 || dim > 48 || arms < 0 || (arms && dim % arms)) return 4;
    i128 dn[3], dd[3];
    int generated = -1;
    if (subtract_generated)
        for (int i = 0; i < element_count && generated < 0; i++)
            if (elements[5 * i + 4] == record[4]) generated = i;
    if (subtract_generated && generated < 0) return 4;
    for (int axis = 0; axis < 3; axis++) {
        if (generated < 0) {
            dn[axis] = record[axis];
            dd[axis] = record[3];
        } else {
            const int *e = elements + 5 * generated;
            dn[axis] = (i128)record[axis] * e[3] - (i128)e[axis] * record[3];
            dd[axis] = (i128)record[3] * e[3];
        }
    }
    int64_t m = 24;
    for (int arm = 0; arm < arms; arm++) {
        i128 num = 0, den = 1;
        for (int axis = 0; axis < 3; axis++) {
            i128 tn = (i128)k[6 * arm + 2 * axis] * dn[axis], td = (i128)k[6 * arm + 2 * axis + 1] * dd[axis];
            num = num * td + tn * den;
            den = den * td;
            i128 g = gcd_i128(num, den);
            if (g > 1) {
                num /= g;
                den /= g;
            }
        }
        if (den < 0) {
            num = -num;
            den = -den;
        }
        if (den > 1000000) return 3;
        pd[arm] = (int64_t)den;
        pn[arm] = (int64_t)(((num % den) + den) % den);
        m = m / gcd64(m, pd[arm]) * pd[arm];
        if (m > 100000) return 3;
    }
    *m_out = m;
    return 0;
}

int exact_record_phases(int dim, int arms, const int64_t *k, int subtract_generated, int element_count, const int *elements,
                        const int *record, int64_t *pn, int64_t *pd, int64_t *m_out)
{
    return record_phases(dim, arms, k, subtract_generated, element_count, elements, record, pn, pd, m_out);
}

int exact_source_words(int gid, int dim, int position, int64_t *words)
{
    int status = little_exact_words(gid, dim, position, words);
    return status == 1 ? 1 : status ? 2 : 0;
}

/* The pieces of W(record) = P S: the arm phases, the Source entries as
 * const.exact words (8 per entry), and the least field m. */
int exact_record_parts(int gid, int dim, int arms, const int64_t *k, int subtract_generated, int element_count, const int *elements,
                       int position, const int *record, int64_t *pn, int64_t *pd, int64_t *words, int64_t *m_out)
{
    int status = record_phases(dim, arms, k, subtract_generated, element_count, elements, record, pn, pd, m_out);
    if (status) return status;
    status = little_exact_words(gid, dim, position, words);
    return status == 1 ? 1 : status ? 2 : 0;
}

int exact_record_matrix(int gid, int dim, int arms, const int64_t *k, int subtract_generated, int element_count,
                        const int *elements, int position, const int *record, int want_m, const cy_field **field_out, cy *w)
{
    if (dim <= 0 || dim > 48 || arms < 0 || (arms && dim % arms)) return 4;
    i128 dn[3], dd[3];
    int generated = -1;
    if (subtract_generated)
        for (int i = 0; i < element_count && generated < 0; i++)
            if (elements[5 * i + 4] == record[4]) generated = i;
    if (subtract_generated && generated < 0) return 4;
    for (int axis = 0; axis < 3; axis++) {
        if (generated < 0) {
            dn[axis] = record[axis];
            dd[axis] = record[3];
        } else {
            const int *e = elements + 5 * generated;
            dn[axis] = (i128)record[axis] * e[3] - (i128)e[axis] * record[3];
            dd[axis] = (i128)record[3] * e[3];
        }
    }
    int64_t pn[48], pd[48], m = 24;
    for (int arm = 0; arm < arms; arm++) {
        i128 num = 0, den = 1;
        for (int axis = 0; axis < 3; axis++) {
            i128 tn = (i128)k[6 * arm + 2 * axis] * dn[axis], td = (i128)k[6 * arm + 2 * axis + 1] * dd[axis];
            num = num * td + tn * den;
            den = den * td;
            i128 g = gcd_i128(num, den);
            if (g > 1) {
                num /= g;
                den /= g;
            }
        }
        if (den < 0) {
            num = -num;
            den = -den;
        }
        if (den > 1000000) return 3;
        pd[arm] = (int64_t)den;
        pn[arm] = (int64_t)(num % den);
        m = m / gcd64(m, pd[arm]) * pd[arm];
        if (m > 100000) return 3;
    }
    if (want_m < 0) {
        *field_out = 0;
        return -(int)m;
    }
    if (want_m) {
        if (want_m % m) return 3;
        m = want_m;
    }
    const cy_field *f = cy_field_for((int)m);
    if (!f) return 3;
    *field_out = f;
    if (!w) return 0;
    int64_t *words = malloc(sizeof(int64_t) * 8 * (size_t)dim * dim);
    cy *s = malloc(sizeof(cy) * (size_t)dim * dim);
    if (!words || !s) {
        free(words);
        free(s);
        return 3;
    }
    int status = little_exact_words(gid, dim, position, words);
    cy_overflow = 0;
    for (int i = 0; i < dim * dim && !status; i++) cy_from_sqrt23(f, words + 8 * i, &s[i]);
    int block = arms ? dim / arms : 0;
    for (int i = 0; i < dim && !status; i++) {
        /* row i of P has at most two nonzero entries */
        int cols[2], n = 0;
        cy vals[2];
        if (!arms) {
            cols[n] = i;
            cy_rat(f, 1, 1, &vals[n++]);
        } else {
            int arm = i / block, a = i % block, start = arm * block;
            cy c, sn;
            cy_cos_sin(f, pn[arm], pd[arm], &c, &sn);
            if (block == 1 || block % 2) {
                if (!cy_is_zero(f, &sn)) status = 4;
                cols[n] = i;
                cy_copy(f, &vals[n++], &c);
            } else {
                int half = block / 2, base = a < half ? a : a - half;
                cols[n] = start + base;
                if (a < half)
                    cy_copy(f, &vals[n++], &c);
                else
                    cy_sub(f, &(cy){{0}, 1}, &sn, &vals[n++]);
                cols[n] = start + half + base;
                cy_copy(f, &vals[n++], a < half ? &sn : &c);
            }
        }
        for (int j = 0; j < dim && !status; j++) {
            cy acc, t;
            cy_zero(f, &acc);
            for (int q = 0; q < n; q++) {
                const cy *entry = &s[cols[q] * dim + j];
                if (cy_is_zero(f, entry) || cy_is_zero(f, &vals[q])) continue;
                cy_mul(f, &vals[q], entry, &t);
                if (cy_is_zero(f, &acc))
                    cy_copy(f, &acc, &t);
                else
                    cy_add(f, &acc, &t, &acc);
            }
            cy_copy(f, &w[i * dim + j], &acc);
        }
    }
    if (!status && cy_overflow) status = 3;
    free(words);
    free(s);
    return status;
}

#include "irreps.h"

int exact_kernel(const cy_field *f, int rows, int cols, const cy *a, int *count, cy *kernel, unsigned char *pivots_out);

/* subgroup_to_orderparam_ over Q(zeta_m): the invariant directions of a
 * little irrep under the records (mapped point actions, signs for the
 * antiunitary ones) from the reduced echelon form of the stacked D - I rows
 * with columns reversed, as count x dim correctly rounded doubles.  Returns
 * the count, or -1 (tables), -2 (field or arithmetic range), -3 (no
 * certified kernel), -4 (a value without a correctly rounded double). */
int exact_subduction_directions(int gid, int param_count, const int64_t *params, int record_count, const int *mapped,
                                const int *signs, double *directions_out)
{
    irrep_context c;
    if (irrep_context_for(gid, param_count, params, &c)) return -1;
    int dim = c.dim, status = 0, kcount = 0;
    cy *w = malloc(sizeof(cy) * (size_t)dim * dim);
    cy *mat = malloc(sizeof(cy) * ((size_t)record_count * dim * dim + 1));
    cy *kernel = malloc(sizeof(cy) * ((size_t)dim * dim + 1));
    int *positions = malloc(sizeof(int) * ((size_t)record_count + 1));
    unsigned char pivots[48];
    if (!w || !mat || !kernel || !positions) {
        status = -2;
        goto out;
    }
    int m = 24;
    for (int r = 0; r < record_count && !status; r++) {
        const cy_field *g;
        positions[r] = bw_find_position(c.element_count, c.elements, mapped + 5 * r);
        int need = positions[r] ? exact_record_matrix(gid, dim, c.arms, c.k, c.subtract, c.element_count, c.elements, positions[r],
                                                      mapped + 5 * r, -1, &g, 0)
                                : 0;
        if (!positions[r]) status = -1;
        else if (need >= 0) status = -2;
        else {
            int a = m, b = -need;
            while (b) {
                int t = a % b;
                a = b;
                b = t;
            }
            m = m / a * -need;
        }
    }
    const cy_field *f = status ? 0 : cy_field_for(m);
    if (!status && !f) status = -2;
    cy_overflow = 0;
    for (int r = 0; r < record_count && !status; r++) {
        const cy_field *g;
        if (exact_record_matrix(gid, dim, c.arms, c.k, c.subtract, c.element_count, c.elements, positions[r], mapped + 5 * r, m, &g, w)) {
            status = -2;
            break;
        }
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                cy *dst = &mat[((size_t)r * dim + i) * dim + j];
                cy_copy(f, dst, &w[i * dim + (dim - 1 - j)]);
                if (signs[r] < 0 && !cy_is_zero(f, dst)) cy_sub(f, &(cy){{0}, 1}, dst, dst);
                if (dim - 1 - j == i) {
                    cy one;
                    cy_rat(f, 1, 1, &one);
                    cy_sub(f, dst, &one, dst);
                }
            }
    }
    if (!status && cy_overflow) status = -2;
    if (!status && exact_kernel(f, record_count * dim, dim, mat, &kcount, kernel, pivots)) status = -3;
    if (!status && kcount > 0) memset(directions_out, 0, sizeof(double) * (size_t)kcount * dim);
    /* the free columns from the last one down, as the float reduction wrote them */
    for (int t = 0; t < kcount && !status; t++) {
        const cy *vec = &kernel[(size_t)(kcount - 1 - t) * dim];
        for (int col = 0; col < dim && !status; col++)
            if (cy_to_double(f, &vec[col], &directions_out[(size_t)t * dim + (dim - 1 - col)])) status = -4;
    }
out:
    free(w);
    free(mat);
    free(kernel);
    free(positions);
    return status ? status : kcount;
}
