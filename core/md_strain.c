/* definition_presentation._strain_mode_definitions with modes.engine.tensor:
 * the rank-2 tensor copies of a GM irrep and the strain definitions they
 * give, in numpy's operation order (linalg.c in place of BLAS / LAPACK, with the
 * transpose flags numpy passes for each operand layout). */


#include "linalg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "md.h"
#include "tables.h"

const char *k_catalog_json(int sg);
const char *direction_matrix_text(int rows, const int64_t *row_lengths, const double *values, double tol, int group_width);
int sp_matrix_rank(const double *m, int rows, int cols, double tol, int *rank);
int little_sparse_matrix(int gid, int dim, int position, double *out);

/* A 2-D array: rows x cols values in C order, or column-major (a numpy
 * transposed view) when fortran. */
typedef struct {
    int rows, cols, fortran;
    double *v;
} nd;

static double at(const nd *a, int i, int j) { return a->fortran ? a->v[j * a->rows + i] : a->v[i * a->cols + j]; }

static int nd_new(nd *a, int rows, int cols)
{
    a->rows = rows;
    a->cols = cols;
    a->fortran = 0;
    a->v = calloc((size_t)(rows * cols + 1), sizeof(double));
    return a->v ? 0 : md_fail(MDE_MEMORY, "out of memory");
}

static nd transposed(const nd *a)
{
    nd out = {a->cols, a->rows, !a->fortran, a->v};
    return out;
}

/* numpy matmul of two 2-D float64 arrays */
static int matmul(const nd *a, const nd *b, nd *out)
{
    MD_TRY(nd_new(out, a->rows, b->cols));
    if (!a->rows || !b->cols) return 0;
    if (!a->cols) return 0;
    int lda = a->fortran ? a->rows : a->cols, ldb = b->fortran ? b->rows : b->cols;
    la_dgemm(a->fortran, b->fortran, a->rows, b->cols, a->cols, a->v, lda, b->v, ldb, out->v, b->cols);
    return 0;
}

/* numpy.linalg.inv */
static int inverse(const nd *a, nd *out, int *singular)
{
    int n = a->rows;
    *singular = 0;
    double *m = malloc(sizeof(double) * (size_t)(n * n + 1));
    if (!m) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) m[i * n + j] = at(a, i, j);
    int status = nd_new(out, n, n);
    if (!status) status = la_inverse(n, m, out->v);
    free(m);
    if (status == 2) return md_fail(MDE_MEMORY, "out of memory");
    if (status) {
        *singular = 1;
        return md_fail(MDE_LINALG, "Singular matrix");
    }
    return 0;
}

static int rank_of(const nd *a, double tol, int *rank)
{
    double *c = malloc(sizeof(double) * (size_t)(a->rows * a->cols + 1));
    if (!c) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++) c[i * a->cols + j] = at(a, i, j);
    int status = sp_matrix_rank(c, a->rows, a->cols, tol, rank);
    free(c);
    return status ? md_fail(MDE_LINALG, "SVD did not converge") : 0;
}

/* ------------------------------------------------------------ tables */

static struct {
    section_view lattice, psettings, full_dim, irr_type, irr_sg, irr_k, star_count, star_ml_pointer, star_special, irr_pointer, matrices_pointer,
        matrices, constants, elements_pointer, elements, point_group, group_order;
    int ready;
} T;

static int tables(void)
{
    if (T.ready) return 0;
    if (section_find("space.ispace_lattice", &T.lattice) || section_find("space.ipoint_op_psettings", &T.psettings) ||
        section_find("little.little_irr_full_dim", &T.full_dim) || section_find("little.little_irr_type", &T.irr_type) ||
        section_find("little.little_irr_space_group", &T.irr_sg) || section_find("little.little_irr_k", &T.irr_k) ||
        section_find("little.little_k_star_count", &T.star_count) || section_find("little.little_k_star_ml_pointer", &T.star_ml_pointer) ||
        section_find("little.little_k_star_special", &T.star_special) ||
        section_find("space.ispace_elements_pointer", &T.elements_pointer) || section_find("space.ispace_elements", &T.elements) ||
        section_find("space.ispace_point_group", &T.point_group) || section_find("space.ipoint_group_order", &T.group_order))
        return 1;
    T.ready = 1;
    return 0;
}

/* ------------------------------------------------------------ tensors */

static int cartesian_basis(int sg, nd *out)
{
    MD_TRY(nd_new(out, 3, 3));
    if (sg >= 143 && sg <= 194) {
        double values[9] = {1.0, -0.5, 0.0, 0.0, sqrt(3.0) / 2.0, 0.0, 0.0, 0.0, 1.0};
        memcpy(out->v, values, sizeof values);
    } else {
        for (int i = 0; i < 3; i++) out->v[4 * i] = 1.0;
    }
    return 0;
}

/* tensor._cartesian_point_matrix */
static int cartesian_point_matrix(int sg, int op, nd *out)
{
    int64_t lattice = section_int(&T.lattice, sg - 1);
    int64_t code = section_int(&T.psettings, (lattice - 1) * 72 + op - 1);
    nd pml, p2c, p2c_inverse, left, conjugated, basis, basis_inverse, product;
    MD_TRY(nd_new(&pml, 3, 3));
    for (int i = 0; i < 9; i++) {
        int64_t r = code % 3;
        if (r < 0) r += 3;
        pml.v[i] = (double)(r - 1);
        code = (code - r) / 3;
    }
    rat matrix[3][3], origin[3];
    if (geo_setting_to_cinter(sg, 1, 0, matrix, origin)) return md_fail(MDE_VALUE, "parent setting Source tables are inconsistent");
    MD_TRY(nd_new(&p2c, 3, 3));
    for (int i = 0; i < 9; i++)
        if (py_float_of_fraction(matrix[i / 3][i % 3], &p2c.v[i])) return md_fail(MDE_ARITH, "setting value out of range");
    int singular;
    MD_TRY(inverse(&p2c, &p2c_inverse, &singular));
    MD_TRY(matmul(&p2c_inverse, &pml, &left));
    MD_TRY(matmul(&left, &p2c, &conjugated));
    nd fractional = transposed(&conjugated);
    MD_TRY(cartesian_basis(sg, &basis));
    MD_TRY(inverse(&basis, &basis_inverse, &singular));
    MD_TRY(matmul(&basis, &fractional, &product));
    MD_TRY(matmul(&product, &basis_inverse, out));
    free(pml.v);
    free(p2c.v);
    free(p2c_inverse.v);
    free(left.v);
    free(conjugated.v);
    free(basis.v);
    free(basis_inverse.v);
    free(product.v);
    return 0;
}

static const int INDICES[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};

/* tensor._rank2_representation: columns[k] = transformed unit tensor k, returned as columns.T */
static int rank2_representation(const nd *rotation, nd *out)
{
    nd columns;
    MD_TRY(nd_new(&columns, 6, 6));
    nd rt = transposed(rotation);
    for (int k = 0; k < 6; k++) {
        nd tensor, left, transformed;
        MD_TRY(nd_new(&tensor, 3, 3));
        tensor.v[INDICES[k][0] * 3 + INDICES[k][1]] = 1.0;
        tensor.v[INDICES[k][1] * 3 + INDICES[k][0]] = 1.0;
        MD_TRY(matmul(rotation, &tensor, &left));
        MD_TRY(matmul(&left, &rt, &transformed));
        for (int c = 0; c < 6; c++) columns.v[k * 6 + c] = transformed.v[INDICES[c][0] * 3 + INDICES[c][1]];
        free(tensor.v);
        free(left.v);
        free(transformed.v);
    }
    *out = transposed(&columns);
    return 0;
}

/* The phase-adjusted real irrep matrix of a generated record (Gamma phases
 * are zero; source.representation.real_phase_operator then phase_product). */
static int irrep_matrix(int64_t gid, int position, int dim, nd *out)
{
    nd sparse;
    MD_TRY(nd_new(&sparse, dim, dim));
    if (little_sparse_matrix((int)gid, dim, position, sparse.v)) return md_fail(MDE_VALUE, "the little irrep has no native irrep tables");
    int64_t sg = section_int(&T.irr_sg, gid - 1), kslot = section_int(&T.irr_k, gid - 1);
    int64_t nmod = section_int(&T.star_count, (sg - 1) * 27 + kslot - 1);
    nd phase;
    MD_TRY(nd_new(&phase, dim, dim));
    if (!nmod) {
        for (int i = 0; i < dim; i++) phase.v[i * dim + i] = 1.0;
    } else {
        if (dim % nmod) return md_fail(MDE_VALUE, "phase count %lld does not divide dim %d for gid=%lld", (long long)nmod, dim, (long long)gid);
        int block = dim / (int)nmod;
        for (int arm = 0; arm < nmod; arm++) {
            int s0 = arm * block;
            double c = 1.0, s = 0.0;
            if (block == 1) {
                phase.v[s0 * dim + s0] = c;
                continue;
            }
            if (block % 2) {
                for (int i = 0; i < block; i++)
                    for (int j = 0; j < block; j++) phase.v[(s0 + i) * dim + s0 + j] = c * (i == j ? 1.0 : 0.0);
                continue;
            }
            int half = block / 2;
            for (int i = 0; i < half; i++)
                for (int j = 0; j < half; j++) {
                    double e = i == j ? 1.0 : 0.0;
                    phase.v[(s0 + i) * dim + s0 + j] = c * e;
                    phase.v[(s0 + i) * dim + s0 + half + j] = s * e;
                    phase.v[(s0 + half + i) * dim + s0 + j] = -s * e;
                    phase.v[(s0 + half + i) * dim + s0 + half + j] = c * e;
                }
        }
    }
    MD_TRY(nd_new(out, dim, dim));
    for (int inner = 0; inner < dim; inner++)
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) out->v[i * dim + j] = out->v[i * dim + j] + phase.v[i * dim + inner] * sparse.v[inner * dim + j];
    free(sparse.v);
    free(phase.v);
    return 0;
}

/* tensor._canonical_copy (in place, layout kept) */
static void canonical_copy(nd *a, double tol)
{
    int n = a->rows * a->cols;
    for (int i = 0; i < n; i++)
        if (fabs(a->v[i]) <= tol) a->v[i] = 0.0;
    int have = 0;
    double low = 0.0;
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++) {
            double value = fabs(at(a, i, j));
            if (value > tol && (!have || value < low)) {
                low = value;
                have = 1;
            }
        }
    if (!have) return;
    for (int i = 0; i < n; i++) a->v[i] /= low;
    for (int i = 0; i < n; i++)
        if (fabs(a->v[i]) <= tol) a->v[i] = 0.0;
    double first = 0.0;
    int found = 0;
    for (int i = 0; i < a->rows && !found; i++)
        for (int j = 0; j < a->cols && !found; j++)
            if (fabs(at(a, i, j)) > tol) {
                first = at(a, i, j);
                found = 1;
            }
    if (first < 0)
        for (int i = 0; i < n; i++) a->v[i] *= -1.0;
}

static int vstack(const nd *items, int count, const nd *extra, nd *out)
{
    int rows = 0, cols = extra ? extra->cols : items[0].cols;
    for (int i = 0; i < count; i++) rows += items[i].rows;
    if (extra) rows += extra->rows;
    MD_TRY(nd_new(out, rows, cols));
    int at_row = 0;
    for (int i = 0; i <= count; i++) {
        const nd *a = i < count ? items + i : extra;
        if (!a) break;
        for (int r = 0; r < a->rows; r++, at_row++)
            for (int c = 0; c < cols; c++) out->v[at_row * cols + c] = at(a, r, c);
    }
    return 0;
}

/* tensor._hexagonal_source_carrier_gauge */
static int hexagonal_gauge(const nd *anchor, nd *out, double tol)
{
    int rows = anchor->rows;
    int *pivots = malloc(sizeof(int) * (size_t)(anchor->cols + 1));
    if (!pivots) return md_fail(MDE_MEMORY, "out of memory");
    int npivots = 0, current_rank = 0;
    for (int col = 0; col < anchor->cols; col++) {
        nd candidate;
        MD_TRY(nd_new(&candidate, rows, npivots + 1));
        for (int r = 0; r < rows; r++) {
            for (int k = 0; k < npivots; k++) candidate.v[r * (npivots + 1) + k] = at(anchor, r, pivots[k]);
            candidate.v[r * (npivots + 1) + npivots] = at(anchor, r, col);
        }
        int rank;
        MD_TRY(rank_of(&candidate, tol, &rank));
        free(candidate.v);
        if (rank <= current_rank) continue;
        pivots[npivots++] = col;
        current_rank = rank;
        if (current_rank == rows) break;
    }
    if (current_rank != rows) {
        free(pivots);
        MD_TRY(nd_new(out, rows, rows));
        for (int i = 0; i < rows; i++) out->v[i * rows + i] = 1.0;
        return 0;
    }
    nd pivot;
    MD_TRY(nd_new(&pivot, rows, npivots));
    for (int r = 0; r < rows; r++)
        for (int k = 0; k < npivots; k++) pivot.v[r * npivots + k] = at(anchor, r, pivots[k]);
    free(pivots);
    int columns_one = 1, rows_one = 1;
    for (int k = 0; k < npivots; k++) {
        int n = 0;
        for (int r = 0; r < rows; r++) n += fabs(pivot.v[r * npivots + k]) > tol;
        columns_one &= n == 1;
    }
    for (int r = 0; r < rows; r++) {
        int n = 0;
        for (int k = 0; k < npivots; k++) n += fabs(pivot.v[r * npivots + k]) > tol;
        rows_one &= n == 1;
    }
    if (columns_one && rows_one) {
        free(pivot.v);
        MD_TRY(nd_new(out, rows, rows));
        for (int i = 0; i < rows; i++) out->v[i * rows + i] = 1.0;
        return 0;
    }
    int singular;
    int status = inverse(&pivot, out, &singular);
    free(pivot.v);
    return status;
}

typedef struct {
    nd *items;
    int count;
} copies;

/* tensor._rank2_tensor_copies */
static int rank2_tensor_copies(int sg, int64_t gid, copies *out)
{
    out->items = 0;
    out->count = 0;
    int dim = (int)section_int(&T.full_dim, gid - 1);
    int64_t pg = section_int(&T.point_group, sg - 1), order_count = section_int(&T.group_order, pg - 1);
    int64_t pointer = section_int(&T.elements_pointer, sg - 1);
    int nops = (int)order_count;
    nd *tensors = calloc((size_t)nops + 1, sizeof(nd)), *irreps = calloc((size_t)nops + 1, sizeof(nd));
    if (!tensors || !irreps) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < nops; i++) {
        int op = (int)section_int(&T.elements, 5 * (pointer - 1 + i) + 4);
        nd rotation;
        MD_TRY(cartesian_point_matrix(sg, op, &rotation));
        MD_TRY(rank2_representation(&rotation, &tensors[i]));
        free(rotation.v);
        MD_TRY(irrep_matrix(gid, i + 1, dim, &irreps[i]));
    }
    double order = (double)nops;
    out->items = calloc(36 * (size_t)dim + 1, sizeof(nd));
    if (!out->items) return md_fail(MDE_MEMORY, "out of memory");
    int accepted_rank = 0;
    for (int tc = 0; tc < 6; tc++)
        for (int sc = 0; sc < dim; sc++) {
            nd stacked;
            MD_TRY(nd_new(&stacked, 6, dim));
            for (int target = 0; target < dim; target++) {
                double column[6];
                for (int e = 0; e < 6; e++) column[e] = 0.0;
                for (int o = 0; o < nops; o++) {
                    double scale = irreps[o].v[target * dim + sc];
                    for (int e = 0; e < 6; e++) {
                        double term = scale * at(&tensors[o], e, tc);
                        column[e] = o == 0 ? 0 + term : column[e] + term;
                    }
                }
                for (int e = 0; e < 6; e++) stacked.v[e * dim + target] = column[e] * (double)dim / order;
            }
            nd copy = transposed(&stacked);
            canonical_copy(&copy, 1e-8);
            if (sg >= 143 && sg <= 194)
                for (int r = 0; r < copy.rows; r++) {
                    if (copy.fortran) copy.v[1 * copy.rows + r] *= 2.0;
                    else copy.v[r * copy.cols + 1] *= 2.0;
                }
            nd combined;
            MD_TRY(vstack(out->items, out->count, &copy, &combined));
            int rank;
            MD_TRY(rank_of(&combined, 1e-8, &rank));
            free(combined.v);
            if (rank > accepted_rank) {
                out->items[out->count++] = copy;
                accepted_rank = rank;
            } else {
                free(copy.v);
            }
        }
    if (out->count && section_int(&T.irr_type, gid - 1) == 1 && sg >= 143 && sg <= 194) {
        nd gauge;
        MD_TRY(hexagonal_gauge(&out->items[0], &gauge, 1e-8));
        for (int k = 0; k < out->count; k++) {
            nd product;
            MD_TRY(matmul(&gauge, &out->items[k], &product));
            for (int i = 0; i < product.rows * product.cols; i++)
                if (fabs(product.v[i]) <= 1e-8) product.v[i] = 0.0;
            free(out->items[k].v);
            out->items[k] = product;
        }
        free(gauge.v);
    }
    for (int i = 0; i < nops; i++) {
        free(tensors[i].v);
        free(irreps[i].v);
    }
    free(tensors);
    free(irreps);
    return 0;
}

static void copies_free(copies *c)
{
    for (int i = 0; i < c->count; i++) free(c->items[i].v);
    free(c->items);
}

/* A group of strain rows (six values each). */
typedef struct {
    int count;
    double (*rows)[6];
} row_group;

typedef struct {
    int count;
    row_group *groups;
} row_groups;

static int push_row(row_groups *g, int new_group, const double *values)
{
    if (new_group) {
        row_group *grown = realloc(g->groups, sizeof(row_group) * (size_t)(g->count + 1));
        if (!grown) return md_fail(MDE_MEMORY, "out of memory");
        g->groups = grown;
        g->groups[g->count].count = 0;
        g->groups[g->count].rows = 0;
        g->count++;
    }
    row_group *last = g->groups + g->count - 1;
    double (*rows)[6] = realloc(last->rows, sizeof(double[6]) * (size_t)(last->count + 1));
    if (!rows) return md_fail(MDE_MEMORY, "out of memory");
    last->rows = rows;
    memcpy(last->rows[last->count++], values, sizeof(double[6]));
    return 0;
}

static void groups_free(row_groups *g)
{
    for (int i = 0; i < g->count; i++) free(g->groups[i].rows);
    free(g->groups);
    memset(g, 0, sizeof *g);
}

/* tensor.selected_rank2_tensor_row_groups */
static int selected_row_groups(int sg, int64_t gid, const md_matrix *direction_matrix, row_groups *out)
{
    memset(out, 0, sizeof *out);
    copies c;
    MD_TRY(rank2_tensor_copies(sg, gid, &c));
    nd direction = {direction_matrix->rows, direction_matrix->cols, 0, direction_matrix->values};
    nd dt = transposed(&direction);
    nd accepted = {0, 6, 0, 0};
    for (int k = 0; k < c.count; k++) {
        nd selected;
        MD_TRY(matmul(&dt, &c.items[k], &selected));
        int new_group = 1;
        for (int r = 0; r < selected.rows; r++) {
            double values[6];
            for (int e = 0; e < 6; e++) {
                values[e] = selected.v[r * 6 + e];
                if (fabs(values[e]) <= 1e-8) values[e] = 0.0;
            }
            int independent;
            if (!accepted.rows) {
                independent = sqrt(la_ddot(6, values, values)) > 1e-8;
            } else {
                int before, after;
                MD_TRY(rank_of(&accepted, 1e-8, &before));
                nd row = {1, 6, 0, values};
                nd stacked;
                MD_TRY(vstack(&accepted, 1, &row, &stacked));
                MD_TRY(rank_of(&stacked, 1e-8, &after));
                free(stacked.v);
                independent = after > before;
            }
            if (!independent) continue;
            double *grown = realloc(accepted.v, sizeof(double) * (size_t)(6 * (accepted.rows + 1)));
            if (!grown) return md_fail(MDE_MEMORY, "out of memory");
            accepted.v = grown;
            memcpy(accepted.v + 6 * accepted.rows, values, sizeof values);
            accepted.rows++;
            MD_TRY(push_row(out, new_group, values));
            new_group = 0;
        }
        free(selected.v);
    }
    free(accepted.v);
    copies_free(&c);
    return 0;
}

/* tensor.totally_symmetric_rank2_rows, each row its own group */
static int totally_symmetric_rows(int sg, int64_t gid, row_groups *out)
{
    memset(out, 0, sizeof *out);
    copies c;
    MD_TRY(rank2_tensor_copies(sg, gid, &c));
    int dim = (int)section_int(&T.full_dim, gid - 1);
    nd identity;
    MD_TRY(nd_new(&identity, dim, dim));
    for (int i = 0; i < dim; i++) identity.v[i * dim + i] = 1.0;
    nd dt = transposed(&identity);
    for (int k = 0; k < c.count; k++) {
        nd selected;
        MD_TRY(matmul(&dt, &c.items[k], &selected));
        for (int r = 0; r < selected.rows; r++) {
            double values[6];
            int any = 0;
            for (int e = 0; e < 6; e++) {
                values[e] = selected.v[r * 6 + e];
                if (fabs(values[e]) <= 1e-8) values[e] = 0.0;
                any |= fabs(values[e]) > 1e-8;
            }
            if (any) MD_TRY(push_row(out, 1, values));
        }
        free(selected.v);
    }
    free(identity.v);
    copies_free(&c);
    return 0;
}

/* ------------------------------------------------------------ definitions */

typedef struct {
    char *symbol, *direction;
    row_groups groups;
    int has_order;
    int64_t order;
} strain_spec;

static int strain_definitions_of(md_request *r, const strain_spec *specs, int nspecs, jbuf *b, int *first)
{
    static const int WEB_ORDER[6] = {0, 3, 5, 4, 2, 1};
    for (int si = 0; si < nspecs; si++) {
        const strain_spec *s = specs + si;
        int npresented = 0, presented_cap = 16;
        double (*presented)[6] = malloc(sizeof(double[6]) * (size_t)presented_cap);
        row_groups shown = {0};
        if (!presented) return md_fail(MDE_MEMORY, "out of memory");
        for (int g = 0; g < s->groups.count; g++) {
            int new_group = 1;
            for (int q = 0; q < s->groups.groups[g].count; q++) {
                const double *source = s->groups.groups[g].rows[q];
                double scale = 0.0;
                for (int e = 0; e < 6; e++) {
                    double value = fabs(source[e]);
                    if (e == 0 || value > scale) scale = value;
                }
                if (scale <= 1e-15) continue;
                double normalized[6];
                for (int e = 0; e < 6; e++) {
                    double value = source[e] / scale;
                    double rounded = nearbyint(value);
                    normalized[e] = fabs(value - rounded) <= 1e-12 ? rounded : value;
                }
                double pivot = 1.0;
                for (int e = 0; e < 6; e++)
                    if (fabs(fabs(normalized[e]) - 1.0) <= 1e-10) {
                        pivot = normalized[e];
                        break;
                    }
                if (pivot < 0.0)
                    for (int e = 0; e < 6; e++) normalized[e] = -normalized[e];
                int seen = 0;
                for (int p = 0; p < npresented && !seen; p++) {
                    int same = 1;
                    for (int e = 0; e < 6 && same; e++) same = fabs(normalized[e] - presented[p][e]) <= 1e-10;
                    seen = same;
                }
                if (seen) continue;
                if (npresented == presented_cap) {
                    presented_cap *= 2;
                    double (*grown)[6] = realloc(presented, sizeof(double[6]) * (size_t)presented_cap);
                    if (!grown) return md_fail(MDE_MEMORY, "out of memory");
                    presented = grown;
                }
                memcpy(presented[npresented++], normalized, sizeof normalized);
                MD_TRY(push_row(&shown, new_group, normalized));
                new_group = 0;
            }
        }
        free(presented);
        char parameters[32];
        int nparameters = 0;
        const char *d = s->direction;
        for (int i = 0; d[i]; i++) {
            char c = d[i];
            if (c < 'a' || c > 'z') continue;
            int before = i > 0 && ((d[i - 1] >= 'a' && d[i - 1] <= 'z') || (d[i - 1] >= 'A' && d[i - 1] <= 'Z'));
            int after = (d[i + 1] >= 'a' && d[i + 1] <= 'z') || (d[i + 1] >= 'A' && d[i + 1] <= 'Z');
            if (before || after) continue;
            int known = 0;
            for (int k = 0; k < nparameters; k++) known |= parameters[k] == c;
            if (!known && nparameters < 32) parameters[nparameters++] = c;
        }
        for (int g = 0; g < shown.count; g++) {
            const row_group *group = shown.groups + g;
            char suffix[32] = "";
            if (shown.count > 1) snprintf(suffix, sizeof suffix, "_%d", g + 1);
            int parameterized = nparameters == group->count;
            for (int q = 0; q < group->count; q++) {
                const double *n = group->rows[q];
                char parameter = parameterized ? parameters[q] : 'a';
                if (!*first) jb_lit(b, ", ");
                *first = 0;
                char *label = md_format("%s[0,0,0]%s%sstrain%s(%c)", r->parent_symbol, s->symbol, s->direction, suffix, parameter);
                if (!label) return md_fail(MDE_MEMORY, "out of memory");
                jb_lit(b, "{\"label\": ");
                jb_cstring(b, label);
                free(label);
                double xx = n[0], xy = n[1], xz = n[2], yy = n[3], yz = n[4], zz = n[5];
                double norm = sqrt(xx * xx + yy * yy + zz * zz + 0.5 * (xy * xy + xz * xz + yz * yz));
                jb_lit(b, ", \"normfactor\": ");
                jb_float(b, norm <= 1e-15 ? 0.0 : 1.0 / norm);
                jb_lit(b, ", \"components\": [");
                for (int e = 0; e < 6; e++) {
                    if (e) jb_lit(b, ", ");
                    jb_float(b, n[WEB_ORDER[e]]);
                }
                jb_lit(b, "], \"tensor_components\": {\"xx\": ");
                jb_float(b, xx);
                jb_lit(b, ", \"xy\": ");
                jb_float(b, xy / 2.0);
                jb_lit(b, ", \"xz\": ");
                jb_float(b, xz / 2.0);
                jb_lit(b, ", \"yy\": ");
                jb_float(b, yy);
                jb_lit(b, ", \"yz\": ");
                jb_float(b, yz / 2.0);
                jb_lit(b, ", \"zz\": ");
                jb_float(b, zz);
                jb_lit(b, "}, \"_spec_order\": ");
                if (s->has_order) jb_int(b, s->order);
                else jb_lit(b, "null");
                MD_TRY(md_factor_metadata(r, s->has_order ? (int)s->order : -1, b));
                jb_lit(b, "}");
            }
        }
        groups_free(&shown);
    }
    return 0;
}

/* The strain specs of build_mode_details and _strain_mode_definitions over them. */
int md_strain_definitions(md_request *r, const md_spec *const *specs, const int *orders, int nspecs, jbuf *b)
{
    if (tables()) return md_fail(MDE_VALUE, "mode Source tables are inconsistent");
    strain_spec *selected = calloc((size_t)nspecs + 1, sizeof(strain_spec));
    if (!selected) return md_fail(MDE_MEMORY, "out of memory");
    int count = 0;
    for (int i = 0; i < nspecs; i++) {
        const md_spec *s = specs[i];
        char *k = s->k_label && s->k_label[0] ? (char *)md_strdup(s->k_label, (int64_t)strlen(s->k_label)) : md_k_label_from_irrep_label(s->label);
        int gm = k && !strcmp(k, "GM");
        free(k);
        if (!gm || s->gid_state != MD_SET || s->dir_state != MD_SET || !s->dir.rows) continue;
        strain_spec *t = selected + count++;
        t->symbol = (char *)md_strdup(s->label, (int64_t)strlen(s->label));
        int64_t *lengths = malloc(sizeof(int64_t) * (size_t)s->dir.rows);
        if (!lengths || !t->symbol) return md_fail(MDE_MEMORY, "out of memory");
        for (int q = 0; q < s->dir.rows; q++) lengths[q] = s->dir.cols;
        const char *text = direction_matrix_text(s->dir.rows, lengths, s->dir.values, 1e-7, 0);
        free(lengths);
        if (!text) return md_fail(MDE_MEMORY, "out of memory");
        t->direction = (char *)md_strdup(text, (int64_t)strlen(text));
        MD_TRY(selected_row_groups(r->sg, s->gid, &s->dir, &t->groups));
        t->has_order = orders[i] >= 0;
        t->order = orders[i];
    }
    if (!count) {
        const char *catalog = k_catalog_json(r->sg);
        jdoc *doc = catalog ? jdoc_parse(catalog, strlen(catalog)) : 0;
        if (!doc) return md_fail(MDE_VALUE, "k-point catalog Source tables are inconsistent");
        const jvalue *points = jget(jdoc_root(doc), "kpoints"), *gamma = 0;
        for (int64_t i = 0; jtruthy(points) && i < points->count && !gamma; i++)
            if (jstr_is(jget(jat(points, i), "label"), "GM")) gamma = jat(points, i);
        const jvalue *irrep = 0;
        const jvalue *irreps = gamma ? jget(gamma, "irreps") : 0;
        for (int64_t i = 0; jtruthy(irreps) && i < irreps->count && !irrep; i++) {
            const jvalue *symbol = jget(jat(irreps, i), "symbol");
            if (jstr_is(symbol, "GM1+") || jstr_is(symbol, "GM1")) irrep = jat(irreps, i);
        }
        if (irrep) {
            int64_t gid;
            if (py_int_of(jget(irrep, "gid"), &gid)) {
                jdoc_free(doc);
                return md_fail(MDE_VALUE, "invalid literal for int()");
            }
            strain_spec *t = selected + count++;
            t->symbol = md_py_str(jget(irrep, "symbol"));
            t->direction = (char *)md_strdup("(a)", 3);
            MD_TRY(totally_symmetric_rows(r->sg, gid, &t->groups));
            t->has_order = 0;
        }
        jdoc_free(doc);
    }
    jb_lit(b, "[");
    int first = 1;
    int status = strain_definitions_of(r, selected, count, b, &first);
    jb_lit(b, "]");
    for (int i = 0; i < count; i++) {
        free(selected[i].symbol);
        free(selected[i].direction);
        groups_free(&selected[i].groups);
    }
    free(selected);
    return status;
}
