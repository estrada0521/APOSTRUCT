/* project_: the normalized real candidate blocks of one little irrep at one
 * site irrep, the selection of project_count independent ones (by little
 * irrep type), and their basis values as (index, value) pairs.
 *
 * Candidate block i is sum_r P_r[i][j] * S_r[0][k] over the parent records,
 * scaled by its largest magnitude; it is reproduced bit for bit with the
 * complex arithmetic of the original (numpy) evaluation.  Independence is
 * decided by complex Gauss-Jordan elimination against the tolerance. */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double re, im;
} cplx;

enum {
    PS_OK = 0,
    PS_SHORT = 1,
    PS_MEMORY = 2,
    PS_SHAPE = 3,
    PS_CAPACITY = 4,
};

static cplx cmul(cplx a, cplx b) { return (cplx){fma(a.re, b.re, -(a.im * b.im)), fma(a.re, b.im, a.im * b.re)}; }
static cplx csub(cplx a, cplx b) { return (cplx){a.re - b.re, a.im - b.im}; }
static cplx cadd(cplx a, cplx b) { return (cplx){a.re + b.re, a.im + b.im}; }
static double cnorm1(cplx a) { return fabs(a.re) + fabs(a.im); }

static cplx cdiv(cplx a, cplx b)
{
    if (fabs(b.re) >= fabs(b.im)) {
        double rat = b.im / b.re, scl = 1.0 / fma(b.im, rat, b.re);
        return (cplx){fma(a.im, rat, a.re) * scl, fma(-a.re, rat, a.im) * scl};
    }
    double rat = b.re / b.im, scl = 1.0 / fma(b.re, rat, b.im);
    return (cplx){fma(a.re, rat, a.im) * scl, fma(a.im, rat, -a.re) * scl};
}

/* zrowop2_ in place on rows x cols. */
static void zrowop2(cplx *m, int rows, int cols, double tol)
{
    int pivot_row = 0;
    for (int col = 0; col < cols && pivot_row < rows; col++) {
        int pivot = -1;
        for (int row = pivot_row; row < rows && pivot < 0; row++)
            if (cnorm1(m[row * cols + col]) > tol) pivot = row;
        if (pivot < 0) continue;
        if (pivot != pivot_row)
            for (int j = 0; j < cols; j++) {
                cplx t = m[pivot_row * cols + j];
                m[pivot_row * cols + j] = m[pivot * cols + j];
                m[pivot * cols + j] = t;
            }
        cplx divisor = m[pivot_row * cols + col];
        for (int j = 0; j < cols; j++) m[pivot_row * cols + j] = cdiv(m[pivot_row * cols + j], divisor);
        for (int row = 0; row < rows; row++) {
            if (row == pivot_row) continue;
            cplx factor = m[row * cols + col];
            if (cnorm1(factor) <= tol) continue;
            for (int j = 0; j < cols; j++) m[row * cols + j] = csub(m[row * cols + j], cmul(factor, m[pivot_row * cols + j]));
        }
        pivot_row++;
    }
    for (int i = 0; i < rows * cols; i++)
        if (hypot(m[i].re, m[i].im) < tol) m[i] = (cplx){0.0, 0.0};
}

static int row_nonzero(const cplx *row, int cols, double tol)
{
    double best = 0.0;
    for (int j = 0; j < cols; j++) {
        double v = cnorm1(row[j]);
        if (v > best) best = v;
    }
    return best > tol;
}

typedef struct {
    int pd, sd;
    const cplx *inverse; /* pd x pd */
    cplx *blocks;        /* pd blocks of pd x sd */
    double tol;
    cplx *kept;
    int kept_rows, kept_cols;
} state;

static const cplx *block_of(const state *s, int index) { return s->blocks + (size_t)index * s->pd * s->sd; }

static double block_max_abs(const state *s, int index)
{
    const cplx *b = block_of(s, index);
    double best = 0.0;
    for (int i = 0; i < s->pd * s->sd; i++) {
        double v = hypot(b[i].re, b[i].im);
        if (v > best) best = v;
    }
    return best;
}

/* block.T @ inverse[:pd, start:start+width] flattened (sd x width). */
static void transformed(const state *s, int index, int start, int width, cplx *out)
{
    const cplx *b = block_of(s, index);
    for (int k = 0; k < s->sd; k++)
        for (int h = 0; h < width; h++) {
            cplx acc = {0.0, 0.0};
            for (int j = 0; j < s->pd; j++) acc = cadd(acc, cmul(b[j * s->sd + k], s->inverse[j * s->pd + start + h]));
            out[k * width + h] = acc;
        }
}

/* Append rows to the kept matrix and reduce; returns the reduced trial in
 * *trial (caller owns) with its row count. */
static int reduce_with(const state *s, const cplx *rows, int count, int cols, cplx **trial, int *trial_rows)
{
    if (s->kept && s->kept_cols != cols) return PS_SHAPE;
    int total = s->kept_rows + count;
    cplx *t = malloc(sizeof(cplx) * (size_t)total * cols);
    if (!t) return PS_MEMORY;
    if (s->kept) memcpy(t, s->kept, sizeof(cplx) * (size_t)s->kept_rows * cols);
    memcpy(t + (size_t)s->kept_rows * cols, rows, sizeof(cplx) * (size_t)count * cols);
    zrowop2(t, total, cols, s->tol);
    *trial = t;
    *trial_rows = total;
    return PS_OK;
}

static void keep(state *s, cplx *trial, int rows, int cols)
{
    free(s->kept);
    s->kept = trial;
    s->kept_rows = rows;
    s->kept_cols = cols;
}

static int nonzero_rows(const cplx *m, int rows, int cols, double tol)
{
    int count = 0;
    for (int r = 0; r < rows; r++) count += row_nonzero(m + (size_t)r * cols, cols, tol);
    return count;
}

static int support_of(const cplx *flat, int length, int block_count, int width, double tol, int *support)
{
    int count = 0;
    for (int b = 0; b < block_count; b++) {
        double best = 0.0;
        for (int i = b * width; i < (b + 1) * width && i < length; i++) {
            double v = cnorm1(flat[i]);
            if (v > best) best = v;
        }
        if (best > tol) support[count++] = b;
    }
    return count;
}

static int same_support(const int *a, int na, const int *b, int nb)
{
    return na == nb && !memcmp(a, b, sizeof(int) * (size_t)na);
}

static int select_type1(state *s, int project_count, int *selected)
{
    int chosen = 0, cols = s->pd * s->sd;
    for (int index = 0; index < s->pd && chosen < project_count; index++) {
        if (block_max_abs(s, index) <= s->tol) continue;
        cplx *trial;
        int rows, status = reduce_with(s, block_of(s, index), 1, cols, &trial, &rows);
        if (status) return status;
        if (!row_nonzero(trial + (size_t)(rows - 1) * cols, cols, s->tol)) {
            free(trial);
            continue;
        }
        keep(s, trial, rows, cols);
        selected[chosen++] = index;
    }
    return chosen == project_count ? PS_OK : PS_SHORT;
}

static int select_type2(state *s, int project_count, int *selected)
{
    int half = s->pd / 2, width = s->sd * half, chosen = 0;
    cplx *rows2 = malloc(sizeof(cplx) * 2 * (size_t)width);
    if (!rows2) return PS_MEMORY;
    for (int index = 0; index < s->pd && chosen < project_count; index++) {
        if (block_max_abs(s, index) <= s->tol) continue;
        transformed(s, index, 0, half, rows2);
        transformed(s, index, half, s->pd - half, rows2 + width);
        cplx *trial;
        int rows, status = reduce_with(s, rows2, 2, width, &trial, &rows);
        if (status) {
            free(rows2);
            return status;
        }
        if (!row_nonzero(trial + (size_t)(rows - 1) * width, width, s->tol)) {
            free(trial);
            continue;
        }
        keep(s, trial, rows, width);
        selected[chosen++] = index;
    }
    free(rows2);
    return chosen == project_count ? PS_OK : PS_SHORT;
}

static int select_type3(state *s, int block_count, int project_count, int *selected)
{
    int pd = s->pd, sd = s->sd, half0 = pd / 2, length = sd * half0;
    double tol = s->tol;
    int status = PS_OK;
    cplx *first = malloc(sizeof(cplx) * (size_t)(length + 1));
    int *supports = calloc((size_t)pd * (block_count + 1) + 1, sizeof(int));
    int *support_sizes = calloc((size_t)pd + 1, sizeof(int));
    int *singleton = calloc((size_t)pd + 1, sizeof(int));
    int *exact = calloc((size_t)pd * (length + 1) + 1, sizeof(int));
    int *exact_sizes = calloc((size_t)pd + 1, sizeof(int));
    int *exact_index = calloc((size_t)pd + 1, sizeof(int));
    if (!first || !supports || !support_sizes || !singleton || !exact || !exact_sizes || !exact_index) {
        status = PS_MEMORY;
        goto done;
    }
    int support_count = 0, use_block_support = 0;
    if (pd > 0 && block_count > 1) {
        int ok = 1;
        for (int index = 0; index < pd; index++) {
            transformed(s, index, 0, half0, first);
            if (length % block_count) {
                ok = 0;
                break;
            }
            support_sizes[index] = support_of(first, length, block_count, length / block_count, tol,
                                              supports + (size_t)index * (block_count + 1));
            support_count++;
        }
        int distinct = 0;
        for (int i = 0; i < support_count; i++) {
            if (!support_sizes[i]) continue;
            int seen = 0;
            for (int j = 0; j < i && !seen; j++)
                seen = support_sizes[j] && same_support(supports + (size_t)i * (block_count + 1), support_sizes[i],
                                                        supports + (size_t)j * (block_count + 1), support_sizes[j]);
            distinct += !seen;
        }
        use_block_support = ok && distinct >= project_count;
    }
    int use_singleton = 0;
    if (pd > 0 && sd == 1 && block_count > 1 && half0 % block_count == 0) {
        int ok = 1, width = half0 / block_count;
        for (int index = 0; index < pd && ok; index++) {
            transformed(s, index, 0, half0, first);
            int found[64];
            int n = block_count <= 64 ? support_of(first, length, block_count, width, tol, found) : 0;
            if (n != 1)
                ok = 0;
            else
                singleton[index] = found[0];
        }
        int distinct = 0;
        for (int i = 0; ok && i < pd; i++) {
            int seen = 0;
            for (int j = 0; j < i && !seen; j++) seen = singleton[j] == singleton[i];
            distinct += !seen;
        }
        use_singleton = ok && distinct >= project_count;
    }
    if (support_count) {
        int exact_count = 0;
        for (int index = 0; index < pd; index++) {
            if (block_max_abs(s, index) <= tol) continue;
            transformed(s, index, 0, half0, first);
            int *e = exact + (size_t)exact_count * (length + 1), n = 0;
            for (int i = 0; i < length; i++)
                if (cnorm1(first[i]) > tol) e[n++] = i;
            if (!n) continue;
            int seen = 0;
            for (int j = 0; j < exact_count && !seen; j++)
                seen = same_support(exact + (size_t)j * (length + 1), exact_sizes[j], e, n);
            if (seen) continue;
            exact_sizes[exact_count] = n;
            exact_index[exact_count++] = index;
        }
        if (exact_count >= project_count) {
            int chosen = 0;
            for (int k = 0; k < exact_count && chosen < project_count && !status; k++) {
                int index = exact_index[k];
                transformed(s, index, 0, half0, first);
                cplx *trial;
                int rows;
                status = reduce_with(s, first, 1, length, &trial, &rows);
                if (status) break;
                if (nonzero_rows(trial, rows, length, tol) <= chosen) {
                    free(trial);
                    continue;
                }
                keep(s, trial, rows, length);
                selected[chosen++] = index;
            }
            if (status || chosen == project_count) goto done;
            free(s->kept);
            s->kept = 0;
            s->kept_rows = 0;
        }
        int groups = 0, group_first[64], group_size[64], group_members[64][64];
        for (int index = 0; index < support_count && groups >= 0; index++) {
            if (!support_sizes[index]) continue;
            int g = -1;
            for (int j = 0; j < groups && g < 0; j++)
                if (same_support(supports + (size_t)group_first[j] * (block_count + 1), support_sizes[group_first[j]],
                                 supports + (size_t)index * (block_count + 1), support_sizes[index]))
                    g = j;
            if (g < 0) {
                if (groups >= 64) {
                    groups = -1;
                    break;
                }
                g = groups++;
                group_first[g] = index;
                group_size[g] = 0;
            }
            if (group_size[g] < 64) group_members[g][group_size[g]++] = index;
        }
        if (groups > 1) {
            int distributed = 0, list[256];
            for (int g = 0; g < groups; g++) {
                int take = sd > 1 ? 1 : support_sizes[group_first[g]];
                for (int i = 0; i < take && i < group_size[g] && distributed < 256; i++) list[distributed++] = group_members[g][i];
            }
            if (distributed == project_count) {
                memcpy(selected, list, sizeof(int) * (size_t)project_count);
                goto done;
            }
        }
    }
    {
        int chosen = 0;
        cplx *row = malloc(sizeof(cplx) * (size_t)(length + block_count + 1));
        if (!row) {
            status = PS_MEMORY;
            goto done;
        }
        for (int index = 0; index < pd && chosen < project_count && !status; index++) {
            if (block_max_abs(s, index) <= tol) continue;
            transformed(s, index, 0, half0, first);
            int cols;
            if (use_block_support) {
                cols = block_count;
                for (int b = 0; b < cols; b++) row[b] = (cplx){0.0, 0.0};
                const int *sup = supports + (size_t)index * (block_count + 1);
                for (int i = 0; i < support_sizes[index]; i++) row[sup[i]] = (cplx){1.0, 0.0};
            } else if (use_singleton && sd == 1 && block_count > 1 && half0 % block_count == 0) {
                cols = block_count;
                int found[64], n = support_of(first, length, block_count, half0 / block_count, tol, found);
                for (int b = 0; b < cols; b++) row[b] = (cplx){0.0, 0.0};
                for (int i = 0; i < n; i++) row[found[i]] = (cplx){1.0, 0.0};
            } else if (sd > 1 && block_count == 2) {
                cols = (length + 1) / 2;
                for (int i = 0; i < cols; i++) row[i] = first[2 * i];
            } else {
                cols = length;
                memcpy(row, first, sizeof(cplx) * (size_t)length);
            }
            if (!row_nonzero(row, cols, tol)) continue;
            int before = s->kept_rows;
            cplx *trial;
            int rows;
            status = reduce_with(s, row, 1, cols, &trial, &rows);
            if (status) break;
            if (nonzero_rows(trial, rows, cols, tol) <= before) {
                free(trial);
                continue;
            }
            keep(s, trial, rows, cols);
            selected[chosen++] = index;
        }
        free(row);
        if (!status && chosen != project_count) status = PS_SHORT;
    }
done:
    free(first);
    free(supports);
    free(support_sizes);
    free(singleton);
    free(exact);
    free(exact_sizes);
    free(exact_index);
    return status;
}

/* parent: record_count x pd x pd; site_rows: record_count x sd (first row of
 * each site matrix); inverse: pd x pd complex (re, im) for types 2 and 3.
 * selected[project_count] gets the chosen block indices (-1 where none was
 * found); active_index /
 * active_value the basis values within basis_length (family stride 144). */
int project_selection(int type, int pd, int sd, int record_count, const double *parent, const double *site_rows,
                      const double *inverse, int block_count, int project_count, int basis_length, int capacity,
                      int *selected, int64_t *active_index, double *active_value, int *active_count)
{
    if (pd <= 0 || sd <= 0 || project_count <= 0) return PS_SHAPE;
    for (int i = 0; i < project_count; i++) selected[i] = -1;
    state s = {pd, sd, (const cplx *)inverse, 0, 1e-10, 0, 0, 0};
    s.blocks = malloc(sizeof(cplx) * (size_t)pd * pd * sd);
    if (!s.blocks) return PS_MEMORY;
    for (int i = 0; i < pd; i++)
        for (int j = 0; j < pd; j++)
            for (int k = 0; k < sd; k++) {
                double acc = 0.0;
                for (int r = 0; r < record_count; r++)
                    acc = acc + parent[((size_t)r * pd + i) * pd + j] * site_rows[(size_t)r * sd + k];
                s.blocks[((size_t)i * pd + j) * sd + k] = (cplx){acc, 0.0};
            }
    for (int i = 0; i < pd; i++) {
        double scale = 0.0;
        cplx *b = s.blocks + (size_t)i * pd * sd;
        for (int e = 0; e < pd * sd; e++)
            if (fabs(b[e].re) > scale) scale = fabs(b[e].re);
        if (scale > 1e-12)
            for (int e = 0; e < pd * sd; e++) b[e] = cdiv(b[e], (cplx){scale, 0.0});
    }
    int status = type == 1 ? select_type1(&s, project_count, selected)
               : type == 2 ? select_type2(&s, project_count, selected)
                           : select_type3(&s, block_count, project_count, selected);
    int count = 0;
    if (!status) {
        for (int family = 0; family < project_count && !status; family++) {
            int start = family * 144;
            if (start >= basis_length) break;
            const cplx *b = block_of(&s, selected[family]);
            for (int row = 0; row < pd && !status; row++)
                for (int col = 0; col < sd; col++) {
                    int index = start + row * 3 + col;
                    if (index >= basis_length) continue;
                    if (count >= capacity) {
                        status = PS_CAPACITY;
                        break;
                    }
                    active_index[count] = index;
                    active_value[count++] = b[row * sd + col].re;
                }
        }
    }
    *active_count = count;
    free(s.blocks);
    free(s.kept);
    return status;
}
