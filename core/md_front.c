/* build_mode_details before its site loop: the selected texts, the child
 * lattice and presentation frame, the viewer placement and the selected
 * Source frame flags. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "md.h"
#include "tables.h"

const char *evaluated_kvector_text(const char *text, int64_t length, int nparams, const char *names, const int64_t *name_lengths,
                                   const char *values, const int64_t *value_lengths);
int mode_text_status(void);
const char *opd_direction_text(const char *text, int64_t length, int little_dim);
int little_star_arm_width(int64_t gid);

static const double DEG_TO_RAD = M_PI / 180.0, RAD_TO_DEG = 180.0 / M_PI;

static struct {
    section_view inter_choice, label_short, freeparam;
    int ready;
} T;

static int tables(void)
{
    if (T.ready) return 0;
    if (section_find("space.ispace_inter_choice", &T.inter_choice) || section_find("space.space_inter_label_short", &T.label_short) ||
        section_find("isotropy.isotropy_orderparam_freeparam", &T.freeparam))
        return 1;
    T.ready = 1;
    return 0;
}

char *md_default_symbol(int sg)
{
    if (tables()) return 0;
    const char *s;
    int64_t n;
    if (section_string(&T.label_short, section_int(&T.inter_choice, sg - 1) - 1, &s, &n)) return 0;
    while (n && (*s == ' ' || *s == '\t')) {
        s++;
        n--;
    }
    while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    return (char *)md_strdup(s, n);
}

/* ------------------------------------------------------------ texts */

static int kvector_text(md_request *r)
{
    const jvalue *k = r->selected_k;
    const char *keys[2] = {"display_kvector", "kvector"};
    const jvalue *text = 0;
    for (int i = 0; i < 2 && !text; i++)
        if (jtruthy(jget(k, keys[i]))) text = jget(k, keys[i]);
    char *vector = text ? md_py_str(text) : (char *)md_strdup("?", 1);
    const jvalue *base = jget(k, "display_parameters");
    if (!jtruthy(base)) base = jget(k, "parameters");
    jbuf names = {0}, values = {0};
    jb_reset(&names);
    jb_reset(&values);
    int64_t *name_lengths = 0, *value_lengths = 0;
    int count = 0;
    const jvalue *sources[2] = {jtruthy(base) ? base : 0, jtruthy(r->k_params) ? r->k_params : 0};
    int64_t total = (sources[0] ? sources[0]->count : 0) + (sources[1] ? sources[1]->count : 0);
    name_lengths = malloc(sizeof(int64_t) * (size_t)(total + 1));
    value_lengths = malloc(sizeof(int64_t) * (size_t)(total + 1));
    int status = !vector || !name_lengths || !value_lengths;
    for (int s = 0; s < 2 && !status; s++) {
        if (!jis(sources[s], J_OBJECT)) continue;
        for (int64_t i = 0; i < sources[s]->count; i++) {
            const jvalue *key = sources[s]->u.members + 2 * i, *value = key + 1;
            size_t before = values.used;
            jb_raw(&names, key->u.s, (size_t)key->count);
            py_str_of(&values, value);
            name_lengths[count] = key->count;
            value_lengths[count] = (int64_t)(values.used - before);
            count++;
        }
    }
    if (status || names.failed || values.failed) {
        free(vector);
        free(name_lengths);
        free(value_lengths);
        jb_free(&names);
        jb_free(&values);
        return md_fail(MDE_MEMORY, "out of memory");
    }
    const char *out = evaluated_kvector_text(vector, (int64_t)strlen(vector), count, names.data, name_lengths, values.data, value_lengths);
    free(vector);
    free(name_lengths);
    free(value_lengths);
    jb_free(&names);
    jb_free(&values);
    if (!out) return mode_text_status() == 2 ? md_fail(MDE_ARITH, "k-vector value out of range") : md_fail(MDE_MEMORY, "out of memory");
    r->display_kvector = (char *)md_strdup(out, (int64_t)strlen(out));
    return r->display_kvector ? 0 : md_fail(MDE_MEMORY, "out of memory");
}

static int direction_text(md_request *r)
{
    int group_width = 0;
    if (r->has_selected_gid) {
        group_width = little_star_arm_width(r->selected_gid);
        if (group_width < 0) return md_fail(group_width == -4 ? MDE_INDEX : MDE_VALUE, group_width == -4 ? "index out of range" : "mode label Source tables are inconsistent");
    }
    if (!jis(r->opd, J_OBJECT)) {
        r->opd_direction = (char *)md_strdup("", 0);
        return r->opd_direction ? 0 : md_fail(MDE_MEMORY, "out of memory");
    }
    const jvalue *direction = jget(r->opd, "direction"), *isotropy = jget(r->opd, "isotropy");
    const jvalue *text = 0;
    if (jtruthy(isotropy) && jtruthy(jget(isotropy, "display_opd"))) text = jget(isotropy, "display_opd");
    else if (jtruthy(direction) && jtruthy(jget(direction, "opd"))) text = jget(direction, "opd");
    char *value = text ? md_py_str(text) : (char *)md_strdup("", 0);
    int64_t little_dim = 0;
    const jvalue *dim = jget(r->irrep, "little_dim");
    if (jtruthy(dim)) {
        if (py_int_of(dim, &little_dim)) {
            free(value);
            return md_fail(MDE_VALUE, "invalid literal for int()");
        }
    } else {
        little_dim = group_width;
    }
    const char *out = value ? opd_direction_text(value, (int64_t)strlen(value), (int)little_dim) : 0;
    r->opd_direction = out ? (char *)md_strdup(out, (int64_t)strlen(out)) : 0;
    free(value);
    return r->opd_direction ? 0 : md_fail(MDE_MEMORY, "out of memory");
}

/* -------------------------------------------------------- float frames */

static double dot_metric(const double *u, const double *v, const double (*m)[3])
{
    double terms[9];
    int n = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) terms[n++] = u[i] * m[i][j] * v[j];
    return py_sum(terms, 9);
}

/* common._cell_from_basis */
static void cell_from_basis(const double *params, int has_basis, const double (*basis)[3], double *out)
{
    if (!has_basis) {
        memcpy(out, params, 6 * sizeof(double));
        return;
    }
    double a = params[0], b = params[1], c = params[2];
    double ca = cos(params[3] * DEG_TO_RAD), cb = cos(params[4] * DEG_TO_RAD), cg = cos(params[5] * DEG_TO_RAD);
    double metric[3][3] = {{a * a, a * b * cg, a * c * cb}, {a * b * cg, b * b, b * c * ca}, {a * c * cb, b * c * ca, c * c}};
    double lengths[3];
    for (int i = 0; i < 3; i++) {
        double d = dot_metric(basis[i], basis[i], (const double (*)[3])metric);
        lengths[i] = sqrt(0.0 > d ? 0.0 : d);
    }
    int pairs[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    out[0] = lengths[0];
    out[1] = lengths[1];
    out[2] = lengths[2];
    for (int k = 0; k < 3; k++) {
        int i = pairs[k][0], j = pairs[k][1];
        double denom = lengths[i] * lengths[j];
        if (denom == 0) {
            out[3 + k] = 0.0;
            continue;
        }
        double value = dot_metric(basis[i], basis[j], (const double (*)[3])metric) / denom;
        if (!(value < 1.0)) value = 1.0;
        if (!(value > -1.0)) value = -1.0;
        out[3 + k] = acos(value) * RAD_TO_DEG;
    }
}

static void row_multiply(const double *row, const double (*matrix)[3], double *out)
{
    for (int axis = 0; axis < 3; axis++) {
        double terms[3] = {row[0] * matrix[0][axis], row[1] * matrix[1][axis], row[2] * matrix[2][axis]};
        out[axis] = py_sum(terms, 3);
    }
}

/* presentation.child_lattice_cartesian with the identity basis */
static void child_lattice_cartesian(const double *lattice, double (*out)[3])
{
    double a = lattice[0], b = lattice[1], c = lattice[2];
    double ar = lattice[3] * DEG_TO_RAD, br = lattice[4] * DEG_TO_RAD, gr = lattice[5] * DEG_TO_RAD;
    double t = c * (cos(ar) - cos(br) * cos(gr)) / sin(gr);
    double inner = 1.0 - py_pow(cos(br), 2.0) - py_pow((cos(ar) - cos(br) * cos(gr)) / sin(gr), 2.0);
    double parent[3][3] = {{a, 0.0, 0.0}, {b * cos(gr), b * sin(gr), 0.0}, {c * cos(br), t, c * sqrt(inner > 0.0 ? inner : 0.0)}};
    double identity[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int i = 0; i < 3; i++) row_multiply(identity[i], (const double (*)[3])parent, out[i]);
}

/* presentation.centering_translations */
static int centering_translations(char letter, double (*out)[3])
{
    out[0][0] = out[0][1] = out[0][2] = 0.0;
    int n = 1;
    switch (letter) {
    case 'A': out[n][0] = 0.0, out[n][1] = 0.5, out[n++][2] = 0.5; break;
    case 'B': out[n][0] = 0.5, out[n][1] = 0.0, out[n++][2] = 0.5; break;
    case 'C': out[n][0] = 0.5, out[n][1] = 0.5, out[n++][2] = 0.0; break;
    case 'I': out[n][0] = 0.5, out[n][1] = 0.5, out[n++][2] = 0.5; break;
    case 'F':
        out[n][0] = 0.0, out[n][1] = 0.5, out[n++][2] = 0.5;
        out[n][0] = 0.5, out[n][1] = 0.0, out[n++][2] = 0.5;
        out[n][0] = 0.5, out[n][1] = 0.5, out[n++][2] = 0.0;
        break;
    case 'R':
        out[n][0] = 2.0 / 3.0, out[n][1] = 1.0 / 3.0, out[n++][2] = 1.0 / 3.0;
        out[n][0] = 1.0 / 3.0, out[n][1] = 2.0 / 3.0, out[n++][2] = 2.0 / 3.0;
        break;
    default: break;
    }
    return n;
}

static char centering_letter(int sg)
{
    char *symbol = md_default_symbol(sg);
    char letter = symbol && *symbol ? (char)(*symbol & ~0x20) : 0;
    free(symbol);
    return letter;
}

/* common._float_matrix_inverse_3 */
static int float_inverse(const double (*m)[3], double (*out)[3])
{
    double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                 m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if (fabs(det) < 1e-12) return 1;
    out[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
    out[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
    out[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
    out[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
    out[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
    out[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
    out[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
    out[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
    out[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
    return 0;
}

/* presentation._inverse3 */
static int presentation_inverse3(const double (*a)[3], double (*out)[3])
{
    double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                 a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if (fabs(det) <= 1e-15) return 1;
    out[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    out[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
    out[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    out[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    out[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
    out[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
    out[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    out[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
    out[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
    return 0;
}

static double round12(double value)
{
    char text[400];
    snprintf(text, sizeof text, "%.12f", value);
    return strtod(text, 0);
}

/* viewer.parent_cell_placement */
static int parent_cell_placement(md_request *r, const double (*basis)[3], const double *origin, char letter)
{
    double parent[3][3];
    if (presentation_inverse3(basis, parent)) return md_fail(MDE_VALUE, "viewer basis is singular");
    double projected[3], base[3], offset[3];
    row_multiply(origin, (const double (*)[3])parent, projected);
    for (int axis = 0; axis < 3; axis++) {
        base[axis] = -projected[axis];
        double terms[3] = {parent[0][axis], parent[1][axis], parent[2][axis]};
        offset[axis] = 0.5 * py_sum(terms, 3);
    }
    double centerings[4][3];
    int ncentering = centering_translations(letter, centerings);
    int have = 0;
    double best_score = 0, best_center[3] = {0}, best_origin[3] = {0};
    int best_index = 0, best_shift[3] = {0};
    for (int k = 0; k < ncentering; k++) {
        double untranslated[3];
        int low[3], high[3];
        for (int axis = 0; axis < 3; axis++) {
            untranslated[axis] = base[axis] + centerings[k][axis] + offset[axis];
            double value = 0.5 - untranslated[axis];
            low[axis] = (int)floor(value);
            high[axis] = (int)ceil(value);
        }
        for (int x = low[0]; x <= high[0]; x++)
            for (int y = low[1]; y <= high[1]; y++)
                for (int z = low[2]; z <= high[2]; z++) {
                    int shift[3] = {x, y, z};
                    double candidate[3], center[3], squares[3];
                    for (int axis = 0; axis < 3; axis++) {
                        candidate[axis] = base[axis] + centerings[k][axis] + (double)shift[axis];
                        center[axis] = candidate[axis] + offset[axis];
                        squares[axis] = py_pow(center[axis] - 0.5, 2.0);
                    }
                    double score = round12(py_sum(squares, 3));
                    int better = !have;
                    if (have) {
                        int c = (score > best_score) - (score < best_score);
                        for (int axis = 0; axis < 3 && !c; axis++) c = (center[axis] > best_center[axis]) - (center[axis] < best_center[axis]);
                        if (!c) c = (k > best_index) - (k < best_index);
                        for (int axis = 0; axis < 3 && !c; axis++) c = (shift[axis] > best_shift[axis]) - (shift[axis] < best_shift[axis]);
                        better = c < 0;
                    }
                    if (better) {
                        have = 1;
                        best_score = score;
                        memcpy(best_center, center, sizeof center);
                        best_index = k;
                        memcpy(best_shift, shift, sizeof shift);
                        memcpy(best_origin, candidate, sizeof candidate);
                    }
                }
    }
    memcpy(r->viewer_basis, parent, sizeof parent);
    memcpy(r->viewer_origin, best_origin, sizeof best_origin);
    return 0;
}

/* common._origin_vector of the selected origin (as structure_runtime reads it) */
void md_origin_vector(const md_request *r, double *out);

static int64_t source_freeparam(md_request *r, int *has)
{
    *has = 0;
    if (r->iso) {
        const jvalue *value = jget(r->iso, "free");
        if (value && !jis(value, J_NULL)) {
            int64_t out;
            if (py_int_of(value, &out)) return 0;
            *has = 1;
            return out;
        }
    }
    int64_t row_id;
    if (!r->iso || py_int_of(jget(r->iso, "row_id"), &row_id)) return 0;
    int64_t count = section_count(&T.freeparam), index = row_id - 1;
    if (index < 0) index += count;
    if (index < 0 || index >= count) return 0;
    *has = 1;
    return section_int(&T.freeparam, index);
}

int md_front(md_request *r)
{
    if (tables()) return md_fail(MDE_VALUE, "mode Source tables are inconsistent");
    r->parent_symbol = md_default_symbol(r->sg);
    if (!r->parent_symbol) return md_fail(MDE_MEMORY, "out of memory");
    MD_TRY(kvector_text(r));
    MD_TRY(md_dynamic_gid(r, &r->has_selected_gid, &r->selected_gid));
    MD_TRY(direction_text(r));
    double identity[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    cell_from_basis(r->cell, r->has_display_basis, (const double (*)[3])r->display_basis, r->child_lattice);
    child_lattice_cartesian(r->child_lattice, r->presentation_basis);
    r->has_centering = 0;
    if (r->has_child_sg && r->child_sg >= 1 && r->child_sg <= 230) {
        double translations[4][3];
        r->has_centering = 1;
        r->centering_factor = sqrt((double)centering_translations(centering_letter(r->child_sg), translations));
    }
    r->has_presentation_inverse =
        r->has_display_basis && !float_inverse((const double (*)[3])r->display_basis, r->presentation_inverse);
    r->freeparam = source_freeparam(r, &r->has_freeparam);
    r->full_parameter_opd = 0;
    if (jis(r->opd, J_OBJECT) && r->has_freeparam) {
        const jvalue *direction = jget(r->opd, "direction");
        int64_t dimension, subduction;
        if (jis(direction, J_OBJECT) && !py_int_of(jget(direction, "dimension"), &dimension) &&
            !py_int_of(jget(direction, "subduction"), &subduction))
            r->full_parameter_opd = r->freeparam == dimension && subduction == dimension;
    }
    double origin[3];
    md_origin_vector(r, origin);
    char letter = r->has_child_sg ? centering_letter(r->child_sg) : 'P';
    MD_TRY(parent_cell_placement(r, r->has_display_basis ? (const double (*)[3])r->display_basis : (const double (*)[3])identity, origin,
                                 letter));
    return 0;
}
