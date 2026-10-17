/* Exact little-irrep matrices W(record) = P(phases) S(position) over
 * Q(zeta_m): the Source entries from const.exact, the star-arm phases as
 * rationals k_arm . delta. */

#include <stdlib.h>
#include <string.h>

#include "cyclo.h"
#include "geometry.h"
#include "rational.h"
#include "source_blob.h"
#include "tables.h"

int bw_find_position(int element_count, const int *elements, const int *record);

/* The table's constant codes of a little irrep at a generated position, dense (code 1 is zero); nonzero when the table
 * does not hold the position. */
int little_table_codes(int gid, int dim, int position, int *codes)
{
    static section_view irr_pointer, pointer, matrices;
    static int ready;
    if (!ready) {
        if (section_find("little.little_irr_full_matrices_irr_pointer", &irr_pointer)
            || section_find("little.little_irr_full_matrices_pointer", &pointer)
            || section_find("little.little_irr_full_matrices", &matrices))
            return 1;
        ready = 1;
    }
    int64_t irr_length = (int64_t)irr_pointer.entry->count, pointer_length = (int64_t)pointer.entry->count;
    if (gid < 1 || gid > irr_length || position < 1 || dim <= 0) return 1;
    int64_t index = section_int(&irr_pointer, gid - 1) + position - 1;
    if (index < 1 || index >= pointer_length) return 1;
    int64_t start = section_int(&pointer, index - 1), end = section_int(&pointer, index);
    if (start < 1 || end < start || (end - 1) * 3 > (int64_t)matrices.entry->count) return 1;
    for (int i = 0; i < dim * dim; i++) codes[i] = 1;
    for (int64_t offset = (start - 1) * 3; offset < (end - 1) * 3; offset += 3) {
        int64_t row = section_int(&matrices, offset), col = section_int(&matrices, offset + 1);
        if (row < 1 || row > dim || col < 1 || col > dim) return 1;
        codes[(row - 1) * dim + col - 1] = (int)section_int(&matrices, offset + 2);
    }
    return 0;
}

int little_codes(int gid, int dim, int position, int *codes);

/* The Source matrix of a little irrep at a generated position as const.exact words (8 per entry), from its codes. */
static int little_exact_words(int gid, int dim, int position, int64_t *words)
{
    static section_view exact;
    static int ready_exact;
    if (!ready_exact) {
        if (section_find("const.exact", &exact) || !section_is_int(&exact)) return 1;
        ready_exact = 1;
    }
    int *codes = malloc(sizeof(int) * (size_t)dim * dim + 1);
    if (!codes || little_codes(gid, dim, position, codes)) {
        free(codes);
        return 1;
    }
    int64_t ncodes = (int64_t)exact.entry->count / 8;
    for (int i = 0; i < dim * dim; i++) {
        int64_t *w = words + 8 * i, code = codes[i];
        if (code < 0 || code >= ncodes) {
            free(codes);
            return 1;
        }
        for (int j = 0; j < 8; j++) w[j] = code == 1 ? (j % 2) : section_int(&exact, 8 * code + j);
        if (!w[1] && code > 1) {
            free(codes);
            return 2;
        }
        if (!w[1]) w[1] = w[3] = w[5] = w[7] = 1;
    }
    free(codes);
    return 0;
}

/* the table's words, for the generators of the derivation */
static int little_table_words(int gid, int dim, int position, int64_t *words)
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

/* the phase part P of W(record) = P S (dim x dim, zero outside its blocks) */
static int record_phase(int dim, int arms, const int64_t *k, int subtract_generated, int element_count, const int *elements, const int *record,
                        int want_m, const cy_field **field_out, cy *p)
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
    if (!p) return 0;
    int status = 0;
    for (int i = 0; i < dim * dim; i++) cy_zero(f, &p[i]);
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
        for (int q = 0; q < n && !status; q++) cy_copy(f, &p[i * dim + cols[q]], &vals[q]);
    }
    if (!status && cy_overflow) status = 3;
    return status;
}


int exact_record_matrix(int gid, int dim, int arms, const int64_t *k, int subtract_generated, int element_count,
                        const int *elements, int position, const int *record, int want_m, const cy_field **field_out, cy *w)
{
    int status = record_phase(dim, arms, k, subtract_generated, element_count, elements, record, want_m, field_out, 0);
    if (status || !w) return status;
    const cy_field *f = *field_out;
    cy *p = malloc(sizeof(cy) * (size_t)dim * dim), *s = malloc(sizeof(cy) * (size_t)dim * dim);
    int64_t *words = malloc(sizeof(int64_t) * 8 * (size_t)dim * dim);
    if (!p || !s || !words) {
        free(p), free(s), free(words);
        return 3;
    }
    status = record_phase(dim, arms, k, subtract_generated, element_count, elements, record, f->m, field_out, p);
    if (!status) status = little_exact_words(gid, dim, position, words);
    cy_overflow = 0;
    for (int i = 0; i < dim * dim && !status; i++) cy_from_sqrt23(f, words + 8 * i, &s[i]);
    for (int i = 0; i < dim && !status; i++)
        for (int j = 0; j < dim; j++) {
            cy acc, t;
            cy_zero(f, &acc);
            for (int q = 0; q < dim; q++) {
                if (cy_is_zero(f, &p[i * dim + q]) || cy_is_zero(f, &s[q * dim + j])) continue;
                cy_mul(f, &p[i * dim + q], &s[q * dim + j], &t);
                cy_add(f, &acc, &t, &acc);
            }
            cy_copy(f, &w[i * dim + j], &acc);
        }
    if (!status && cy_overflow) status = 3;
    free(p), free(s), free(words);
    return status;
}

#include "irreps.h"
int irrep_context_current(int gid, irrep_context *out);

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

/* const.exact words (8) of a Source constant code; nonzero without an exact form */
int exact_const_words(int64_t code, int64_t *words)
{
    static section_view exact;
    static int ready;
    if (!ready) {
        if (section_find("const.exact", &exact) || !section_is_int(&exact)) return 1;
        ready = 1;
    }
    if (code < 0 || code >= (int64_t)exact.entry->count / 8) return 1;
    for (int j = 0; j < 8; j++) words[j] = section_int(&exact, 8 * code + j);
    if (!words[1]) {
        if (code > 1) return 1;
        for (int j = 0; j < 8; j++) words[j] = j % 2;
    }
    return 0;
}

/* ---- the little irrep matrices from those at a generating set ---- */

/* x -> R x + t on pml columns, R the transpose of geo_vrot's rows */
static void position_rotation(int sg, int op, rat (*r)[3])
{
    for (int i = 0; i < 3; i++) {
        rat unit[3] = {rat_int(i == 0), rat_int(i == 1), rat_int(i == 2)}, row[3];
        geo_vrot(sg, op, unit, row);
        for (int j = 0; j < 3; j++) r[j][i] = row[j];
    }
}

/* S at every generated position of gid from S at a generating set: for r = r_a r_g, W(r) = W(r_a) W(r_g) with W = P S,
 * so S(p) = P(r)^T W(r_a) W(r_g) (P is orthogonal).  counts: [positions, generators, derived, matching the table,
 * differing]; status nonzero when the context or arithmetic fails. */
typedef struct {
    int gid, n, dim, ng;
    const cy_field *f;
    cy *S;               /* n x dim x dim */
    unsigned char *seed; /* positions read from the table */
} derived_t;

static int derive_positions(int gid, const int64_t *kparam, derived_t *out, int *counts)
{
    memset(counts, 0, sizeof(int) * 5);
    memset(out, 0, sizeof *out);
    irrep_context c;
    if (irrep_context_current(gid, &c) && irrep_context_for(gid, 4, kparam, &c)) return 1;
    int sg = irrep_space_group(gid), n = c.element_count, dim = c.dim;
    rat (*rot)[3][3] = malloc(sizeof(rat[3][3]) * (size_t)n);
    for (int p = 0; p < n; p++) position_rotation(sg, c.elements[5 * p + 4], rot[p]);
    /* products of positions by rotation */
    int *product = malloc(sizeof(int) * (size_t)n * n);
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++) {
            rat m[3][3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) {
                    m[i][j] = rat_int(0);
                    for (int q = 0; q < 3; q++) m[i][j] = rat_add(m[i][j], rat_mul(rot[a][i][q], rot[b][q][j]));
                }
            product[a * n + b] = -1;
            for (int q = 0; q < n && product[a * n + b] < 0; q++) {
                int same = 1;
                for (int i = 0; i < 9 && same; i++) same = rat_eq(m[i / 3][i % 3], rot[q][i / 3][i % 3]);
                if (same) product[a * n + b] = q;
            }
            if (product[a * n + b] < 0) return free(rot), free(product), 2;
        }
    /* a generating set, greedily in position order */
    unsigned char *in = calloc((size_t)n, 1);
    int gens[48], ng = 0, identity = -1;
    for (int p = 0; p < n && identity < 0; p++) {
        int ok = 1;
        for (int i = 0; i < 9; i++) ok &= rat_eq(rot[p][i / 3][i % 3], rat_int(i / 3 == i % 3));
        if (ok) identity = p;
    }
    in[identity] = 1;
    for (int p = 0; p < n; p++) {
        if (in[p]) continue;
        gens[ng++] = p;
        in[p] = 1;
        for (int grown = 1; grown;) {
            grown = 0;
            for (int a = 0; a < n; a++)
                if (in[a])
                    for (int b = 0; b < n; b++)
                        if (in[b] && !in[product[a * n + b]]) in[product[a * n + b]] = grown = 1;
        }
    }
    counts[0] = n;
    counts[1] = ng;
    /* the field: every product record's phases */
    const cy_field *f = 0;
    int m = 24;
    for (int a = 0; a < n; a++) {
        const cy_field *fs;
        int need = -record_phase(dim, c.arms, c.k, c.subtract, n, c.elements, c.elements + 5 * a, -1, &fs, 0);
        if (need <= 0) return free(rot), free(product), free(in), 3;
        m = m / (int)gcd64(m, need) * need;
    }
    size_t sq = (size_t)dim * dim;
    cy *S = malloc(sizeof(cy) * sq * n), *W = malloc(sizeof(cy) * sq * n), *P = malloc(sizeof(cy) * sq), *prod = malloc(sizeof(cy) * sq);
    cy *table = malloc(sizeof(cy) * sq);
    int64_t *words = malloc(sizeof(int64_t) * 8 * sq);
    unsigned char *have = calloc((size_t)n, 1);
    int status = 0;
    /* W at the base records, S from the table at the generators (and the identity) */
    for (int q = -1; q < ng && !status; q++) {
        int p = q < 0 ? identity : gens[q];
        int ph = record_phase(dim, c.arms, c.k, c.subtract, n, c.elements, c.elements + 5 * p, m, &f, P);
        if (ph || little_table_words(gid, dim, p + 1, words)) {
            status = ph ? 40 + ph : 7;
            break;
        }
        for (size_t i = 0; i < sq; i++) cy_from_sqrt23(f, words + 8 * i, &S[(size_t)p * sq + i]);
        have[p] = 1;
    }
    /* W(base record) = P S for every known position; BFS by right multiplication with generators */
    int *queue = malloc(sizeof(int) * (size_t)n), head = 0, tail = 0;
    for (int p = 0; p < n; p++)
        if (have[p]) queue[tail++] = p;
#define WOF(p)                                                                                                            \
    do {                                                                                                                  \
        record_phase(dim, c.arms, c.k, c.subtract, n, c.elements, c.elements + 5 * (p), m, &f, P);                       \
        for (size_t i_ = 0; i_ < (size_t)dim; i_++)                                                                       \
            for (size_t j_ = 0; j_ < (size_t)dim; j_++) {                                                                 \
                cy acc_, t_;                                                                                              \
                cy_zero(f, &acc_);                                                                                        \
                for (size_t q_ = 0; q_ < (size_t)dim; q_++) {                                                             \
                    if (cy_is_zero(f, &P[i_ * dim + q_])) continue;                                                       \
                    cy_mul(f, &P[i_ * dim + q_], &S[(size_t)(p) * sq + q_ * dim + j_], &t_);                             \
                    cy_add(f, &acc_, &t_, &acc_);                                                                         \
                }                                                                                                         \
                W[(size_t)(p) * sq + i_ * dim + j_] = acc_;                                                              \
            }                                                                                                             \
    } while (0)
    for (int i = 0; i < tail; i++) WOF(queue[i]);
    while (head < tail && !status) {
        int a = queue[head++];
        for (int q = 0; q < ng && !status; q++) {
            int g = gens[q], p = product[a * n + g];
            if (have[p]) continue;
            /* r = r_a r_g: rotation R_a R_g, translation R_a t_g + t_a */
            const int *ra = c.elements + 5 * a, *rg = c.elements + 5 * g;
            rat t[3];
            int64_t den = 1;
            for (int i = 0; i < 3; i++) {
                t[i] = rat_make(ra[i], ra[3]);
                for (int j = 0; j < 3; j++) t[i] = rat_add(t[i], rat_mul(rot[a][i][j], rat_make(rg[j], rg[3])));
                den = den / gcd64(den, (int64_t)t[i].d) * (int64_t)t[i].d;
            }
            int rec[5] = {(int)(t[0].n * (den / t[0].d)), (int)(t[1].n * (den / t[1].d)), (int)(t[2].n * (den / t[2].d)), (int)den,
                          c.elements[5 * p + 4]};
            if (record_phase(dim, c.arms, c.k, c.subtract, n, c.elements, rec, m, &f, P)) {
                status = 5;
                break;
            }
            /* W(r_a) W(r_g), then S = P^T that */
            for (size_t i = 0; i < (size_t)dim; i++)
                for (size_t j = 0; j < (size_t)dim; j++) {
                    cy acc, tt;
                    cy_zero(f, &acc);
                    for (size_t l = 0; l < (size_t)dim; l++) {
                        if (cy_is_zero(f, &W[(size_t)a * sq + i * dim + l])) continue;
                        cy_mul(f, &W[(size_t)a * sq + i * dim + l], &W[(size_t)g * sq + l * dim + j], &tt);
                        cy_add(f, &acc, &tt, &acc);
                    }
                    prod[i * dim + j] = acc;
                }
            for (size_t i = 0; i < (size_t)dim; i++)
                for (size_t j = 0; j < (size_t)dim; j++) {
                    cy acc, tt;
                    cy_zero(f, &acc);
                    for (size_t l = 0; l < (size_t)dim; l++) {
                        if (cy_is_zero(f, &P[l * dim + i])) continue;
                        cy_mul(f, &P[l * dim + i], &prod[l * dim + j], &tt);
                        cy_add(f, &acc, &tt, &acc);
                    }
                    S[(size_t)p * sq + i * dim + j] = acc;
                }
            have[p] = 1;
            WOF(p);
            queue[tail++] = p;
            counts[2]++;
        }
    }
#undef WOF
    for (int p = 0; p < n; p++)
        if (!have[p]) status = status ? status : 8;
    if (cy_overflow) status = 6;
    out->gid = gid;
    out->n = n;
    out->dim = dim;
    out->ng = ng;
    out->f = f;
    out->S = S;
    out->seed = calloc((size_t)n, 1);
    if (out->seed) {
        out->seed[identity] = 1;
        for (int q = 0; q < ng; q++) out->seed[gens[q]] = 1;
    }
    free(rot), free(product), free(in), free(W), free(P), free(prod), free(table), free(words), free(have), free(queue);
    if (status) {
        free(out->S), free(out->seed);
        memset(out, 0, sizeof *out);
    }
    return status;
}

static const int64_t DERIVE_KPARAM[4] = {1, 2, 3, 5};

int little_matrices_from_generators(int gid, const int64_t *kparam, int *counts)
{
    derived_t d;
    int status = derive_positions(gid, kparam, &d, counts);
    if (status) return status;
    int64_t *words = malloc(sizeof(int64_t) * 8 * (size_t)d.dim * d.dim);
    cy table;
    size_t sq = (size_t)d.dim * d.dim;
    for (int p = 0; p < d.n; p++) {
        if (little_table_words(gid, d.dim, p + 1, words)) continue;
        int same = 1;
        for (size_t i = 0; i < sq && same; i++) {
            cy diff;
            cy_from_sqrt23(d.f, words + 8 * i, &table);
            cy_sub(d.f, &table, &d.S[(size_t)p * sq + i], &diff);
            same = cy_is_zero(d.f, &diff);
        }
        counts[same ? 3 : 4]++;
    }
    free(words), free(d.S), free(d.seed);
    return 0;
}

/* codes of the little irrep matrices: the table's at the generators and the identity, derived elsewhere (each entry
 * matched exactly against the constant table); APO_TABLE_LITTLE reads the table everywhere */
int little_codes(int gid, int dim, int position, int *codes)
{
    static int table_only = -1, busy;
    if (table_only < 0) table_only = getenv("APO_TABLE_LITTLE") != 0;
    if (table_only || busy) return little_table_codes(gid, dim, position, codes);
    enum { SLOTS = 8 };
    static derived_t cache[SLOTS];
    static int *cache_codes[SLOTS], next;
    int slot = -1;
    for (int i = 0; i < SLOTS && slot < 0; i++)
        if (cache[i].gid == gid && cache[i].dim == dim) slot = i;
    if (slot < 0) {
        slot = next++ % SLOTS;
        free(cache[slot].S), free(cache[slot].seed), free(cache_codes[slot]);
        cache_codes[slot] = 0;
        int counts[5];
        busy = 1;
        int derived = derive_positions(gid, DERIVE_KPARAM, &cache[slot], counts);
        busy = 0;
        if (derived || cache[slot].dim != dim) {
            free(cache[slot].S), free(cache[slot].seed);
            memset(&cache[slot], 0, sizeof cache[slot]);
            return little_table_codes(gid, dim, position, codes);
        }
        /* each derived entry as a constant code */
        static section_view exact;
        static int ready;
        if (!ready && (section_find("const.exact", &exact) || !section_is_int(&exact))) return 1;
        ready = 1;
        const cy_field *f = cache[slot].f;
        int64_t ncodes = (int64_t)exact.entry->count / 8;
        cy *value = malloc(sizeof(cy) * (size_t)ncodes);
        unsigned char *exactly = calloc((size_t)ncodes, 1);
        for (int64_t q = 1; q < ncodes; q++) {
            int64_t w[8];
            for (int j = 0; j < 8; j++) w[j] = section_int(&exact, 8 * q + j);
            if (!w[1] && q > 1) continue;
            if (!w[1]) w[1] = w[3] = w[5] = w[7] = 1;
            cy_from_sqrt23(f, w, &value[q]);
            exactly[q] = 1;
        }
        size_t sq = (size_t)dim * dim, n = (size_t)cache[slot].n;
        cache_codes[slot] = malloc(sizeof(int) * sq * n);
        for (size_t p = 0; p < n; p++)
            for (size_t i = 0; i < sq; i++) {
                const cy *x = &cache[slot].S[p * sq + i];
                int code = 0;
                if (cy_is_zero(f, x)) code = 1;
                for (int64_t q = 2; q < ncodes && !code; q++) {
                    if (!exactly[q]) continue;
                    cy diff;
                    cy_sub(f, x, &value[q], &diff);
                    if (cy_is_zero(f, &diff)) code = (int)q;
                }
                cache_codes[slot][p * sq + i] = code;
            }
        free(value), free(exactly);
        free(cache[slot].S);
        cache[slot].S = 0;
    }
    derived_t *d = &cache[slot];
    if (position < 1 || position > d->n) return 1;
    if (d->seed[position - 1]) return little_table_codes(gid, dim, position, codes);
    const int *src = cache_codes[slot] + (size_t)(position - 1) * dim * dim;
    for (int i = 0; i < dim * dim; i++) {
        if (!src[i]) return little_table_codes(gid, dim, position, codes);
        codes[i] = src[i];
    }
    return 0;
}
