/* Parameter solutions of one k vector on the arms of a little-k star.
 *
 * Each arm is origin + params . directions; its system is pre-eliminated into
 * integer left-inverse rows and consistency rows over one scale
 * (k/families.c, solver_for). For every reciprocal shift in
 * (-1, 0, 1)^3, x-major, the rhs target + shift - origin either fails the
 * consistency rows or yields params = left . rhs / scale. All arithmetic is
 * exact; fractions are returned reduced with positive denominators. */

#include <stdint.h>

#include "intmath.h"


enum {
    KS_FIRST = 0,
    KS_ALL_MOD_ONE = 1,
};

enum {
    KS_OVERFLOW = -1,
    KS_CAPACITY = -2,
};

static i128 gcd_abs(i128 a, i128 b) { return gcd_i128(a, b); }

static int reduce(i128 *num, i128 *den)
{
    if (*den < 0) {
        *den = -*den;
        *num = -*num;
    }
    i128 g = gcd_abs(*num, *den);
    if (g > 1) {
        *num /= g;
        *den /= g;
    }
    const i128 limit = (i128)1 << 62;
    return *den > limit || *num > limit || *num < -limit;
}

/* Returns the number of solutions written (KS_FIRST: 0 or 1), or a negative
 * error. out: capacity x kdim (num, den) pairs. */
int kstar_solutions(int mode, int kdim, int arm_count, const int64_t *origin,
                    const int *valid, const int64_t *left, const int64_t *consistency,
                    const int64_t *scale, const int64_t *target, int capacity, int64_t *out)
{
    int count = 0;
    int consistency_rows = 3 - kdim;
    for (int arm = 0; arm < arm_count; arm++) {
        if (!valid[arm]) continue;
        const int64_t *o = origin + 6 * arm;
        const int64_t *l = left + 9 * arm;
        const int64_t *c = consistency + 9 * arm;
        for (int sx = -1; sx <= 1; sx++)
            for (int sy = -1; sy <= 1; sy++)
                for (int sz = -1; sz <= 1; sz++) {
                    int shift[3] = {sx, sy, sz};
                    i128 num[3], den[3];
                    i128 common = 1;
                    for (int axis = 0; axis < 3; axis++) {
                        i128 tn = target[2 * axis], td = target[2 * axis + 1];
                        i128 on = o[2 * axis], od = o[2 * axis + 1];
                        num[axis] = (tn + shift[axis] * td) * od - on * td;
                        den[axis] = td * od;
                        if (reduce(&num[axis], &den[axis])) return KS_OVERFLOW;
                        common = common / gcd_abs(common, den[axis]) * den[axis];
                        if (common > ((i128)1 << 62)) return KS_OVERFLOW;
                    }
                    i128 x[3];
                    for (int axis = 0; axis < 3; axis++) x[axis] = num[axis] * (common / den[axis]);
                    int consistent = 1;
                    for (int row = 0; row < consistency_rows && consistent; row++) {
                        const int64_t *r = c + 3 * row;
                        if ((i128)r[0] * x[0] + (i128)r[1] * x[1] + (i128)r[2] * x[2]) consistent = 0;
                    }
                    if (!consistent) continue;
                    i128 solved_num[3], solved_den[3];
                    for (int row = 0; row < kdim; row++) {
                        const int64_t *r = l + 3 * row;
                        solved_num[row] = (i128)r[0] * x[0] + (i128)r[1] * x[1] + (i128)r[2] * x[2];
                        solved_den[row] = (i128)scale[arm] * common;
                        if (reduce(&solved_num[row], &solved_den[row])) return KS_OVERFLOW;
                    }
                    if (mode == KS_FIRST) {
                        for (int row = 0; row < kdim; row++) {
                            out[2 * row] = (int64_t)solved_num[row];
                            out[2 * row + 1] = (int64_t)solved_den[row];
                        }
                        return 1;
                    }
                    for (int row = 0; row < kdim; row++) {
                        i128 n = solved_num[row] % solved_den[row];
                        if (n < 0) n += solved_den[row];
                        i128 d = solved_den[row];
                        reduce(&n, &d);
                        solved_num[row] = n;
                        solved_den[row] = d;
                    }
                    int seen = 0;
                    for (int previous = 0; previous < count && !seen; previous++) {
                        int same = 1;
                        for (int row = 0; row < kdim && same; row++)
                            same = out[(previous * kdim + row) * 2] == (int64_t)solved_num[row]
                                   && out[(previous * kdim + row) * 2 + 1] == (int64_t)solved_den[row];
                        seen = same;
                    }
                    if (seen) continue;
                    if (count >= capacity) return KS_CAPACITY;
                    for (int row = 0; row < kdim; row++) {
                        out[(count * kdim + row) * 2] = (int64_t)solved_num[row];
                        out[(count * kdim + row) * 2 + 1] = (int64_t)solved_den[row];
                    }
                    count++;
                }
    }
    return count;
}
