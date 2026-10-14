/* The (space operation, kernel coset) multiplication tables fed to
 * isotropy_traverse: space_product[l][r] is the pair of op_l * op_r and
 * space_fraction[s][f] the coset of op_s applied to kernel fraction f.
 * Translations are exact rationals; a coset is identified by its residue
 * (v K^-1 mod 1) * modulus against the kernel fractions' residues. */

#include <stdint.h>
#include <stdlib.h>

typedef __int128 i128;

enum {
    TT_OK = 0,
    TT_RANGE = 1,
    TT_COSET = 2,
    TT_PAIR = 3,
    TT_POINT = 4,
    TT_DISTINCT = 5,
};

typedef struct {
    int64_t r[3];
    int index;
} tt_residue;

typedef struct {
    const int64_t *adj;
    i128 det, modulus;
    const tt_residue *sorted;
    int count;
} tt_cosets;

static int tt_cmp(const void *left, const void *right)
{
    const tt_residue *a = left, *b = right;
    for (int i = 0; i < 3; i++)
        if (a->r[i] != b->r[i]) return a->r[i] < b->r[i] ? -1 : 1;
    return 0;
}

static int tt_mul(i128 a, i128 b, i128 *out) { return __builtin_mul_overflow(a, b, out); }
static int tt_add(i128 a, i128 b, i128 *out) { return __builtin_add_overflow(a, b, out); }

/* Coset index (1-based) of num/den, 0 when not a kernel coset, -1 on overflow. */
static int tt_coset(const tt_cosets *c, const i128 *num, i128 den)
{
    tt_residue key;
    i128 scale;
    if (tt_mul(den, c->det, &scale)) return -1;
    for (int col = 0; col < 3; col++) {
        i128 total = 0, product;
        for (int row = 0; row < 3; row++)
            if (tt_mul(num[row], c->adj[row * 3 + col], &product) || tt_add(total, product, &total)) return -1;
        if (tt_mul(total, c->modulus, &total)) return -1;
        if (total % scale) return 0;
        i128 value = (total / scale) % c->modulus;
        if (value < 0) value += c->modulus;
        key.r[col] = (int64_t)value;
    }
    const tt_residue *hit = bsearch(&key, c->sorted, (size_t)c->count, sizeof key, tt_cmp);
    return hit ? hit->index + 1 : 0;
}

/* f R / (f_den * rot_den) for the composite rotation R of one space op. */
static int tt_rotate(const int64_t *rot, const int64_t *fraction, i128 rot_den, i128 *num, i128 *den)
{
    for (int col = 0; col < 3; col++) {
        i128 total = 0, product;
        for (int row = 0; row < 3; row++)
            if (tt_mul(fraction[row], rot[row * 3 + col], &product) || tt_add(total, product, &total)) return 1;
        num[col] = total;
    }
    return tt_mul(fraction[3], rot_den, den);
}

int traversal_tables(int space_count, const int64_t *translations, const int *ops,
                     const int64_t *rotations, int64_t rotation_den, int magnetic,
                     const int64_t *point_matrices, const int *mlt, int fraction_count,
                     const int64_t *fractions, const int64_t *kernel_adjugate, int64_t kernel_det,
                     int64_t *residues, int64_t *modulus_out, int *space_product, int *space_fraction)
{
    if (space_count <= 0 || fraction_count <= 0 || !kernel_det) return TT_RANGE;
    i128 modulus = 1;
    for (int f = 0; f < fraction_count; f++)
        for (int col = 0; col < 3; col++) {
            i128 num = 0, den = (i128)fractions[4 * f + 3] * kernel_det;
            for (int row = 0; row < 3; row++) num += (i128)fractions[4 * f + row] * kernel_adjugate[row * 3 + col];
            if (!den) return TT_RANGE;
            i128 a = num < 0 ? -num : num, b = den < 0 ? -den : den;
            while (a) {
                i128 t = b % a;
                b = a;
                a = t;
            }
            i128 reduced = (den < 0 ? -den : den) / b, g = modulus, h = reduced;
            while (h) {
                i128 t = g % h;
                g = h;
                h = t;
            }
            if (tt_mul(modulus / g, reduced, &modulus) || modulus > ((i128)1 << 40)) return TT_RANGE;
        }
    tt_residue *sorted = malloc(sizeof *sorted * (size_t)fraction_count);
    if (!sorted) return TT_RANGE;
    for (int f = 0; f < fraction_count; f++) {
        for (int col = 0; col < 3; col++) {
            i128 num = 0, den = (i128)fractions[4 * f + 3] * kernel_det;
            for (int row = 0; row < 3; row++) num += (i128)fractions[4 * f + row] * kernel_adjugate[row * 3 + col];
            i128 value = (num * modulus / den) % modulus;
            if (value < 0) value += modulus;
            residues[3 * f + col] = sorted[f].r[col] = (int64_t)value;
        }
        sorted[f].index = f;
    }
    qsort(sorted, (size_t)fraction_count, sizeof *sorted, tt_cmp);
    for (int f = 1; f < fraction_count; f++)
        if (!tt_cmp(&sorted[f - 1], &sorted[f])) {
            free(sorted);
            return TT_DISTINCT;
        }
    *modulus_out = (int64_t)modulus;
    tt_cosets cosets = {kernel_adjugate, kernel_det, modulus, sorted, fraction_count};
    int status = TT_OK;
    for (int s = 0; s < space_count && !status; s++)
        for (int f = 0; f < fraction_count; f++) {
            i128 num[3], den;
            if (tt_rotate(rotations + 9 * s, fractions + 4 * f, rotation_den, num, &den)) {
                status = TT_RANGE;
                break;
            }
            int index = tt_coset(&cosets, num, den);
            if (index <= 0) {
                status = index < 0 ? TT_RANGE : TT_COSET;
                break;
            }
            space_fraction[s * fraction_count + f] = index;
        }
    for (int l = 0; l < space_count && !status; l++)
        for (int r = 0; r < space_count; r++) {
            const int64_t *left = translations + 4 * l, *right = translations + 4 * r;
            i128 rnum[3], rden;
            if (tt_rotate(rotations + 9 * l, right, rotation_den, rnum, &rden)) {
                status = TT_RANGE;
                break;
            }
            int point = 0;
            if (magnetic) {
                point = mlt[(ops[r] - 1) * 144 + (ops[l] - 1)];
            } else {
                int64_t product[9];
                const int64_t *lm = point_matrices + 9 * l, *rm = point_matrices + 9 * r;
                for (int row = 0; row < 3; row++)
                    for (int col = 0; col < 3; col++) {
                        int64_t total = 0;
                        for (int k = 0; k < 3; k++) total += rm[row * 3 + k] * lm[k * 3 + col];
                        product[row * 3 + col] = total;
                    }
                for (int s = space_count - 1; s >= 0 && !point; s--) {
                    int same = 1;
                    for (int i = 0; i < 9 && same; i++) same = point_matrices[9 * s + i] == product[i];
                    if (same) point = ops[s];
                }
                if (!point) {
                    status = TT_POINT;
                    break;
                }
            }
            int found = 0;
            for (int s = 0; s < space_count && !found; s++) {
                if (ops[s] != point) continue;
                const int64_t *base = translations + 4 * s;
                i128 den, num[3], a, b, c;
                if (tt_mul(left[3], rden, &den) || tt_mul(den, base[3], &den)) {
                    status = TT_RANGE;
                    break;
                }
                for (int axis = 0; axis < 3; axis++) {
                    if (tt_mul((i128)left[axis] * rden, base[3], &a) || tt_mul(rnum[axis] * left[3], base[3], &b)
                        || tt_mul((i128)base[axis] * left[3], rden, &c)) {
                        status = TT_RANGE;
                        break;
                    }
                    num[axis] = a + b - c;
                }
                if (status) break;
                int index = tt_coset(&cosets, num, den);
                if (index < 0) {
                    status = TT_RANGE;
                    break;
                }
                if (index) {
                    space_product[2 * (l * space_count + r)] = s + 1;
                    space_product[2 * (l * space_count + r) + 1] = index;
                    found = 1;
                }
            }
            if (status) break;
            if (!found) {
                status = TT_PAIR;
                break;
            }
        }
    free(sorted);
    return status;
}
