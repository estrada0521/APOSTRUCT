/* The first Wyckoff row of an inter setting holding a fractional point: the modular integer solve of the
 * formula through its Smith normal form, else the least-squares solve
 * (minimum-norm least squares, math/linalg.c), then the Source-branch parameters of
 * the original coordinates.  Float arithmetic follows the original order. */

#include "linalg.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "pynum.h"
#include "pysum.h"

enum {
    WM_FOUND = 1,
    WM_NONE = 0,
    WM_TABLES = -1,
    WM_MEMORY = -2,
    WM_LAPACK = -3,
};

typedef struct {
    double base[3], cols[3][3];
    int nonzero[3];
} formula;

typedef struct {
    int64_t multiplicity;
    formula cinter;
    int branch_count;
    formula *branches;
} row_data;

typedef struct {
    int ready, sg, row_count;
    row_data *rows;
} setting_data;

static setting_data *settings;
static int setting_capacity;

static double to_double(rat r) { return (double)(int64_t)r.n / (double)(int64_t)r.d; }

static void formula_from(rat (*v)[3], formula *f)
{
    for (int a = 0; a < 3; a++) f->base[a] = to_double(v[0][a]);
    for (int i = 0; i < 3; i++) {
        f->nonzero[i] = 0;
        for (int a = 0; a < 3; a++) {
            f->cols[i][a] = to_double(v[i + 1][a]);
            f->nonzero[i] |= f->cols[i][a] != 0.0;
        }
    }
}

static setting_data *setting_for(int sg, int setting_id)
{
    if (geo_load()) return 0;
    if (!settings) {
        setting_capacity = 1024;
        settings = calloc((size_t)setting_capacity, sizeof *settings);
        if (!settings) return 0;
    }
    if (setting_id < 1 || setting_id >= setting_capacity) return 0;
    setting_data *s = &settings[setting_id];
    if (s->ready) return s->sg == sg ? s : 0;
    int start, count;
    geo_wyckoff_row_range(sg, &start, &count);
    rat matrix[3][3], origin[3];
    if (geo_setting_to_cinter(sg, 1, setting_id, matrix, origin)) return 0;
    int64_t records[GEO_MAX_OPS][5];
    int record_count = geo_sg_records(sg, records);
    s->rows = calloc((size_t)count + 1, sizeof *s->rows);
    if (!s->rows) return 0;
    for (int r = 0; r < count; r++) {
        row_data *row = &s->rows[r];
        row->multiplicity = geo_multiplicity(sg, geo_wyckoff_site_pg(start + r));
        rat pml[4][3], cinter[4][3];
        geo_inter_wyckoff_vectors(setting_id, r, pml);
        for (int i = 0; i < 4; i++) rat_row_mul(pml[i], (const rat (*)[3])matrix, cinter[i]);
        for (int a = 0; a < 3; a++) cinter[0][a] = rat_add(cinter[0][a], origin[a]);
        formula_from(cinter, &row->cinter);
        row->branches = malloc(sizeof(formula) * (size_t)(record_count + 1));
        if (!row->branches) return 0;
        row->branch_count = record_count;
        for (int b = 0; b < record_count; b++) {
            rat v[4][3];
            int op = (int)records[b][4];
            geo_vrot(sg, op, pml[0], v[0]);
            for (int a = 0; a < 3; a++) v[0][a] = rat_add(v[0][a], rat_make(records[b][a], records[b][3]));
            for (int i = 1; i < 4; i++) geo_vrot(sg, op, pml[i], v[i]);
            formula_from(v, &row->branches[b]);
        }
    }
    if (geo_take_error()) return 0;
    s->sg = sg;
    s->row_count = count;
    s->ready = 1;
    return s;
}

static double mod_distance(double x) { return fabs(x - rint(x)); }

static double dot3(const double *a, const double *b) { return (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]; }

static int64_t trunc_div(int64_t n, int64_t d)
{
    int64_t q = llabs(n) / llabs(d);
    return ((n < 0) ^ (d < 0)) ? -q : q;
}

/* row_echelon_form_ on an n x n integer matrix with an n x n transform. */
static void row_echelon(int64_t m[3][3], int64_t t[3][3])
{
    int64_t aug[3][6];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            aug[r][c] = m[r][c];
            aug[r][3 + c] = t[r][c];
        }
    for (int p = 0; p < 3; p++) {
        for (;;) {
            int best = -1, nonzero = 0;
            for (int r = p; r < 3; r++)
                if (aug[r][p] != 0) {
                    nonzero++;
                    if (best < 0 || llabs(aug[r][p]) < llabs(aug[best][p])) best = r;
                }
            if (!nonzero) break;
            if (best != p)
                for (int c = 0; c < 6; c++) {
                    int64_t x = aug[p][c];
                    aug[p][c] = aug[best][c];
                    aug[best][c] = x;
                }
            if (aug[p][p] < 0)
                for (int c = 0; c < 6; c++) aug[p][c] = -aug[p][c];
            if (nonzero == 1) break;
            int64_t pivot = aug[p][p];
            for (int r = p + 1; r < 3; r++) {
                if (!aug[r][p]) continue;
                int64_t q = trunc_div(aug[r][p], pivot);
                for (int c = p; c < 6; c++) aug[r][c] -= q * aug[p][c];
            }
        }
    }
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            m[r][c] = aug[r][c];
            t[r][c] = aug[r][3 + c];
        }
}

static int offdiagonal(int64_t m[3][3])
{
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            if (r != c && m[r][c]) return 1;
    return 0;
}

/* _smith_normal_decomp_iso: diagonal, left, right (right returned transposed). */
static int smith(const int64_t rows[3][3], double diag[3][3], double left[3][3], double right[3][3])
{
    int64_t m[3][3], l[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, r[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    memcpy(m, rows, sizeof m);
    int done = 0;
    for (int it = 0; it < 100 && !done; it++) {
        row_echelon(m, l);
        if (!offdiagonal(m)) {
            done = 1;
            break;
        }
        int64_t t[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) t[i][j] = m[j][i];
        row_echelon(t, r);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) m[i][j] = t[j][i];
        if (!offdiagonal(m)) done = 1;
    }
    if (!done) return 1;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            diag[i][j] = (double)m[i][j];
            left[i][j] = (double)l[i][j];
            right[i][j] = (double)r[j][i];
        }
    return 0;
}

/* _solve_integer_modular_system for a 3x3 matrix[axis][param]. */
static int modular_solve(const double m[3][3], const double *rhs, double tol, double *solution)
{
    int64_t rounded[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double nearest = rint(m[i][j]);
            if (!(fabs(m[i][j] - nearest) <= tol + 1e-05 * fabs(nearest))) return 0;
            rounded[i][j] = (int64_t)nearest;
        }
    double diag[3][3], left[3][3], right[3][3];
    if (smith((const int64_t (*)[3])rounded, diag, left, right)) return 0;
    double transformed[3], y[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++) transformed[i] = dot3(left[i], rhs);
    for (int i = 0; i < 3; i++) {
        double d = diag[i][i];
        double value = pymod1(transformed[i]);
        if (fabs(value - 1.0) <= tol) value = 0.0;
        if (fabs(d) <= tol) {
            if (fabs(value) > tol) return 0;
            continue;
        }
        y[i] = value / d;
    }
    for (int i = 0; i < 3; i++) solution[i] = dot3(right[i], y);
    for (int i = 0; i < 3; i++)
        if (!(mod_distance(dot3(m[i], solution) - rhs[i]) <= tol)) return 0;
    return 1;
}

static int least_squares(const double m[3][3], const double *rhs, double *solution)
{
    int rank;
    return la_lstsq(3, 3, 1, &m[0][0], rhs, DBL_EPSILON * 3, solution, &rank) ? WM_MEMORY : 0;
}

/* _solve_wyckoff_params_from_vectors: 1 with params (mask of named
 * columns) when the formula holds target, 0 when not, negative on error. */
static int solve(const formula *f, const double *target, double tol, int modulo, double *params)
{
    if (!f->nonzero[0] && !f->nonzero[1] && !f->nonzero[2]) {
        for (int a = 0; a < 3; a++)
            if (!(mod_distance(target[a] - f->base[a]) <= tol)) return 0;
        return 1;
    }
    double m[3][3], rhs[3], solution[3];
    for (int a = 0; a < 3; a++) {
        for (int i = 0; i < 3; i++) m[a][i] = f->cols[i][a];
        rhs[a] = target[a] - f->base[a];
    }
    if (!modular_solve((const double (*)[3])m, rhs, tol, solution)) {
        double wrapped[3];
        for (int a = 0; a < 3; a++) {
            double w = pymod1(rhs[a]);
            wrapped[a] = fabs(w - 1.0) <= tol + 1e-05 ? 0.0 : w;
        }
        int status = least_squares((const double (*)[3])m, wrapped, solution);
        if (status) return status;
        double worst = 0.0;
        int unordered = 0;
        for (int a = 0; a < 3; a++) {
            double r = fabs(dot3(m[a], solution) - wrapped[a]);
            unordered |= isnan(r);
            if (r > worst) worst = r;
        }
        if (!unordered && worst > tol) return 0;
    }
    for (int i = 0; i < 3; i++) params[i] = f->nonzero[i] ? (modulo ? pymod1(solution[i]) : solution[i]) : 0.0;
    return 1;
}

static double source_param_value(double value, double tol)
{
    double n = value;
    while (fabs(n - 1.0) < fabs(n)) n -= 1.0;
    while (fabs(n + 1.0) < fabs(n)) n += 1.0;
    return fabs(n) <= tol ? 0.0 : n;
}


/* out: row offset (0-based); params[3] with mask[3] of named columns;
 * representative[3].  Returns 1 found, 0 none, negative on error. */
/* source_values: the Source-setting coordinates, brought to pml exactly
 * through Fraction(str(value)) before branch matching. */
int match_inter_wyckoff(int sg, int setting_id, const double *fract, const double *source_values, int wanted_mult,
                        double tol, int *row_out, double *params, int *mask, double *representative)
{
    setting_data *s = setting_for(sg, setting_id);
    if (!s) return WM_TABLES;
    rat to_cinter[3][3], inverse[3][3], origin[3], shifted[3], pml[3];
    if (geo_setting_to_cinter(sg, 1, setting_id, to_cinter, origin) || rat_inverse3((const rat (*)[3])to_cinter, inverse))
        return WM_TABLES;
    rat_overflow = 0;
    for (int axis = 0; axis < 3; axis++) shifted[axis] = rat_sub(repr_rat(source_values[axis]), origin[axis]);
    rat_row_mul(shifted, (const rat (*)[3])inverse, pml);
    double source_target[3];
    for (int axis = 0; axis < 3; axis++) source_target[axis] = rat_to_double(pml[axis]);
    if (rat_overflow) {
        rat_overflow = 0;
        if (affine_repr_double(source_values, origin, (const rat (*)[3])inverse, source_target)) return WM_TABLES;
    }
    for (int r = 0; r < s->row_count; r++) {
        const row_data *row = &s->rows[r];
        if (wanted_mult > 0 && row->multiplicity != wanted_mult) continue;
        double values[3];
        int found = solve(&row->cinter, fract, tol, 1, values);
        if (found < 0) return found;
        if (!found) continue;
        for (int b = 0; b < row->branch_count; b++) {
            double branch_values[3];
            int hit = solve(&row->branches[b], source_target, tol, 0, branch_values);
            if (hit < 0) return hit;
            if (!hit) continue;
            for (int i = 0; i < 3; i++)
                values[i] = row->branches[b].nonzero[i] ? source_param_value(branch_values[i], tol) : values[i];
            break;
        }
        for (int i = 0; i < 3; i++) {
            mask[i] = row->cinter.nonzero[i];
            params[i] = values[i];
        }
        for (int a = 0; a < 3; a++) {
            double p[3];
            for (int i = 0; i < 3; i++) p[i] = row->cinter.cols[i][a] * (mask[i] ? values[i] : 0.0);
            representative[a] = row->cinter.base[a] + apostruct_builtin_sum3(p[0], p[1], p[2]);
        }
        *row_out = r;
        return 1;
    }
    return 0;
}

/* _solve_wyckoff_params_from_vectors for one formula: base[3], cols[3][3]
 * (parameter column i at cols[3 * i]); params/mask out. */
int solve_wyckoff_formula(const double *base, const double *cols, const double *target, double tol, int modulo,
                          double *params, int *mask)
{
    formula f;
    for (int a = 0; a < 3; a++) f.base[a] = base[a];
    for (int i = 0; i < 3; i++) {
        f.nonzero[i] = 0;
        for (int a = 0; a < 3; a++) {
            f.cols[i][a] = cols[3 * i + a];
            f.nonzero[i] |= f.cols[i][a] != 0.0;
        }
        mask[i] = f.nonzero[i];
    }
    return solve(&f, target, tol, modulo, params);
}

/* match_inter_wyckoff for every image of one site: rows_out[i] the row
 * offset or -1 when image i matches no row.  Returns 0 or a status. */
int match_inter_wyckoff_images(int sg, int setting_id, int image_count, const double *fracts,
                               const double *source_values, int wanted_mult, double tol, int *rows_out,
                               double *params, int *masks, double *representatives)
{
    for (int i = 0; i < image_count; i++) {
        int found = match_inter_wyckoff(sg, setting_id, fracts + 3 * i, source_values, wanted_mult, tol, rows_out + i,
                                        params + 3 * i, masks + 3 * i, representatives + 3 * i);
        if (found < 0) return found;
        if (!found) rows_out[i] = -1;
    }
    return 0;
}
