/* The child structure: the selected presentation, the parent setting bridge and every parent site's canonical child layout. */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"
#include "md.h"
#include "tables.h"

int child_presentation_setting(int parent_sg, int child_sg, int parent_setting, const int64_t *basis, const int64_t *origin,
                               const int64_t *presentation_pml, const int64_t *origin_words, int capacity, int *setting_ids,
                               int64_t *matrix_words, int64_t *origin_out_words);
int child_site_layout(int parent_sg, int child_sg, int parent_row_id, const int64_t *split_basis, const int64_t *split_origin,
                      const int64_t *parent_param_words, int setting_count, const int *setting_ids, const int64_t *matrix_words,
                      const int64_t *origin_words, int64_t *label_memo, int memo_capacity, int source_row_id, int source_site_pg,
                      const int64_t *site_param_words, const int64_t *presentation_pml, const int64_t *representative_words,
                      const int64_t *bridge_words, const int64_t *presentation_words, const double *presentation_basis,
                      const double *presentation_origin,
                      int centering_letter, int point_operation_count, int magnetic_group, int capacity, int64_t *sites_out,
                      int64_t *atoms_out, double *positions_out, int64_t *rows_out, int *counts);
int solve_wyckoff_formula(const double *base, const double *cols, const double *target, double tol, int modulo, double *params,
                          int *mask);

md_error md_err;

int md_fail(int kind, const char *format, ...)
{
    md_err.kind = kind;
    md_err.status = 0;
    va_list args;
    va_start(args, format);
    vsnprintf(md_err.message, sizeof md_err.message, format, args);
    va_end(args);
    return kind;
}

static int md_kernel(int kind, int status)
{
    md_err.kind = kind;
    md_err.status = status;
    md_err.message[0] = 0;
    return kind;
}

static struct {
    section_view inter_number, inter_choice, label_short, wyckoff_label, wyckoff_pointer, wyckoff_count,
        inter_wyckoff_pointer, inter_wyckoff_fract, mag_type4_pointer, mag_type4_supergroup, mag_nlabel;
    int ready;
} T;

int md_tables(void)
{
    if (T.ready) return 0;
    if (geo_load() || section_find("settings.inter_settings.number", &T.inter_number) ||
        section_find("group.space_groups.inter_choice", &T.inter_choice) || section_find("settings.inter_settings.short", &T.label_short) ||
        section_find("settings.wyckoff.label", &T.wyckoff_label) || section_find("settings.wyckoff.by_sg_start", &T.wyckoff_pointer) ||
        section_find("settings.wyckoff.by_sg_count", &T.wyckoff_count) ||
        section_find("settings.inter_settings.wyckoff_start", &T.inter_wyckoff_pointer) ||
        section_find("settings.inter_settings.wyckoff", &T.inter_wyckoff_fract) || mag_tables() ||
        section_find("group.magnetic_groups.type4", &T.mag_type4_pointer) ||
        section_find("group.type4_supergroups.supergroup", &T.mag_type4_supergroup) || section_find("group.magnetic_groups.og", &T.mag_nlabel))
        return 1;
    T.ready = 1;
    return 0;
}

const char *md_strdup(const char *text, int64_t length)
{
    char *out = malloc((size_t)length + 1);
    if (out) {
        memcpy(out, text, (size_t)length);
        out[length] = 0;
    }
    return out;
}

char *md_format(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    int n = vsnprintf(0, 0, format, args);
    va_end(args);
    char *out = malloc((size_t)n + 1);
    if (!out) return 0;
    va_start(args, format);
    vsnprintf(out, (size_t)n + 1, format, args);
    va_end(args);
    return out;
}

char *md_py_str(const jvalue *value)
{
    jbuf b = {0};
    jb_reset(&b);
    py_str_of(&b, value);
    if (b.failed) {
        jb_free(&b);
        return 0;
    }
    return b.data;
}

static int strip_equal(const section_view *view, int64_t index, const char *text, int64_t length)
{
    const char *s;
    int64_t n;
    if (section_string(view, index, &s, &n)) return 0;
    while (n && (*s == ' ' || *s == '\t')) {
        s++;
        n--;
    }
    while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    return n == length && !memcmp(s, text, (size_t)length);
}

static md_str section_text(const section_view *view, int64_t index)
{
    md_str out = {"", 0};
    const char *s;
    int64_t n;
    if (section_string(view, index, &s, &n)) return out;
    while (n && (*s == ' ' || *s == '\t')) {
        s++;
        n--;
    }
    while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    out.data = s;
    out.length = n;
    return out;
}

const jvalue *md_isotropy(const jvalue *opd)
{
    if (!jis(opd, J_OBJECT)) return 0;
    const jvalue *iso = jget(opd, "isotropy");
    if (!jtruthy(iso)) iso = opd;
    return jis(iso, J_OBJECT) ? iso : 0;
}

static int exact_int(const jvalue *v, int64_t *out)
{
    if (!jis(v, J_INT)) return 1;
    *out = v->u.i;
    return 0;
}
static int exact_integer_tuple(const jvalue *v, int size, const char *name, int64_t *out)
{
    if (!jis(v, J_ARRAY) || v->count != size) return md_fail(MDE_VALUE, "%s requires %d entries", name, size);
    for (int i = 0; i < size; i++)
        if (exact_int(jat(v, i), out + i)) return md_fail(MDE_TYPE, "%s must contain exact integers", name);
    return 0;
}

/* site_transport._parent_setting_bridge */
int md_parent_setting_bridge(int sg, int has_setting, int64_t setting_id, int *has_bridge, rat (*matrix)[3], rat *shift)
{
    *has_bridge = 0;
    if (!has_setting) return 0;
    if (setting_id < 1 || setting_id > section_count(&T.inter_number))
        return md_fail(MDE_INDEX, "inter setting %lld out of range", (long long)setting_id);
    if (section_int(&T.inter_number, setting_id - 1) != sg) return md_fail(MDE_VALUE, "parent setting does not belong to the parent space group");
    int default_id = (int)section_int(&T.inter_choice, sg - 1);
    if (setting_id == default_id) return 0;
    rat default_matrix[3][3], default_origin[3], selected_matrix[3][3], selected_origin[3], inverse[3][3];
    if (geo_setting_to_cinter(sg, 0, default_id, default_matrix, default_origin) ||
        geo_setting_to_cinter(sg, 0, (int)setting_id, selected_matrix, selected_origin) ||
        rat_inverse3((const rat (*)[3])default_matrix, inverse))
        return md_fail(MDE_VALUE, "parent setting Source tables are inconsistent");
    for (int row = 0; row < 3; row++) rat_row_mul(inverse[row], (const rat (*)[3])selected_matrix, matrix[row]);
    for (int col = 0; col < 3; col++) {
        rat sum = rat_int(0);
        for (int row = 0; row < 3; row++) sum = rat_add(sum, rat_mul(default_origin[row], matrix[row][col]));
        shift[col] = rat_sub(selected_origin[col], sum);
    }
    *has_bridge = 1;
    return 0;
}
int md_cell_params(md_request *r)
{
    const jvalue *lattice = jget(r->cif, "lattice");
    static const char *keys[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
    for (int i = 0; i < 6; i++)
        if (float_cif_number(jtruthy(lattice) ? jget(lattice, keys[i]) : 0, &r->cell[i])) {
            jbuf b = {0};
            jb_reset(&b);
            if (jtruthy(lattice)) py_repr_of(&b, lattice);
            else jb_lit(&b, "{}");
            md_fail(MDE_VALUE, "cannot build lattice parameters from CIF summary: %s", b.data ? b.data : "");
            jb_free(&b);
            return MDE_VALUE;
        }
    return 0;
}
static void selected_subgroup_number(md_request *r)
{
    r->has_child_sg = 0;
    if (!r->iso) return;
    const jvalue *subgroup = jget(r->iso, "subgroup");
    if (jis(subgroup, J_OBJECT)) {
        const jvalue *number = jget(subgroup, "ordinary_number");
        subgroup = jtruthy(number) ? number : jget(subgroup, "number");
    }
    int64_t value;
    if (!py_int_of(subgroup, &value)) {
        r->has_child_sg = 1;
        r->child_sg = (int)value;
    }
}

/* "(a,b,c),(d,e,f),(g,h,i)" of rationals; nonzero when it does not read */
static int exact_basis_text(const jvalue *text, rat (*out)[3])
{
    if (!jis(text, J_STRING)) return 1;
    const char *c = text->u.s, *e = c + text->count;
    for (int row = 0; row < 3; row++) {
        if (row && (c >= e || *c++ != ',')) return 1;
        if (c >= e || *c++ != '(') return 1;
        for (int col = 0; col < 3; col++) {
            const char *start = c;
            while (c < e && *c != ',' && *c != ')') c++;
            if (c >= e || *c != (col < 2 ? ',' : ')') || py_fraction_of_text(start, c - start, &out[row][col])) return 1;
            c++;
        }
    }
    return c != e;
}
static void display_basis(md_request *r)
{
    r->has_display_basis = 0;
    const jvalue *basis = r->iso ? jget(r->iso, "basis") : 0;
    if (!jis(basis, J_ARRAY) || basis->count != 3) return;
    for (int row = 0; row < 3; row++)
        if (!jis(jat(basis, row), J_ARRAY) || jat(basis, row)->count != 3) return;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
            if (py_float_of(jat(jat(basis, row), col), &r->display_basis[row][col])) return;
    r->has_display_basis = 1;
    r->has_exact_basis = !exact_basis_text(r->iso ? jget(r->iso, "basis_text") : 0, r->exact_basis);
    for (int row = 0; row < 3 && r->has_exact_basis; row++)
        for (int col = 0; col < 3; col++)
            if (py_float_of_fraction(r->exact_basis[row][col], &r->display_basis[row][col])) r->has_exact_basis = 0;
}
static void selected_origin(md_request *r)
{
    r->origin_kind = 0;
    const jvalue *origin = r->iso ? jget(r->iso, "origin") : 0;
    if (jis(origin, J_STRING)) {
        const char *s = origin->u.s, *e = s + origin->count;
        while (s < e && (*s == ' ' || (*s >= '\t' && *s <= '\r'))) s++;
        if (s < e) {
            r->origin_kind = 1;
            r->origin_text.data = origin->u.s;
            r->origin_text.length = origin->count;
        }
        return;
    }
    if (!jis(origin, J_ARRAY) || origin->count < 4) return;
    for (int i = 0; i < 4; i++)
        if (py_int_of(jat(origin, i), &r->origin_record[i])) return;
    r->origin_kind = 2;
}
static int basis_cinter_to_pml(md_request *r)
{
    r->has_presentation_pml = 0;
    rat rows[3][3];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            rat raw, snapped;
            if (r->has_exact_basis) {
                rows[row][col] = r->exact_basis[row][col];
                continue;
            }
            if (py_fraction_of_float(r->display_basis[row][col], &raw))
                return md_fail(MDE_ARITH, "presentation basis value out of range");
            snapped = py_limit_denominator(raw, 384);
            double gap;
            rat_overflow = 0;
            rat difference = rat_sub(raw, snapped);
            if (rat_overflow || py_float_of_fraction(difference, &gap)) return md_fail(MDE_ARITH, "presentation basis value out of range");
            rows[row][col] = fabs(gap) <= 5e-9 ? snapped : raw;
        }
    if (r->has_bridge) {
        rat inverse[3][3], product[3][3];
        if (rat_inverse3((const rat (*)[3])r->bridge_matrix, inverse)) return md_fail(MDE_VALUE, "singular parent setting bridge");
        for (int row = 0; row < 3; row++) rat_row_mul(rows[row], (const rat (*)[3])inverse, product[row]);
        memcpy(rows, product, sizeof rows);
    }
    rat pml[3][3], origin[3], inverse[3][3];
    if (geo_setting_to_cinter(r->sg, 1, 0, pml, origin) || rat_inverse3((const rat (*)[3])pml, inverse))
        return md_fail(MDE_VALUE, "parent setting Source tables are inconsistent");
    for (int row = 0; row < 3; row++) {
        rat out[3];
        rat_row_mul(rows[row], (const rat (*)[3])inverse, out);
        for (int col = 0; col < 3; col++) {
            if (out[col].d != 1) return 0;
            r->presentation_pml[3 * row + col] = (int64_t)out[col].n;
        }
    }
    if (rat_overflow) return md_fail(MDE_ARITH, "presentation basis value out of range");
    r->has_presentation_pml = 1;
    return 0;
}
static void selected_magnetic_group(md_request *r)
{
    r->magnetic_group = 0;
    if (!r->iso) return;
    const jvalue *subgroup = jget(r->iso, "subgroup");
    if (!jis(subgroup, J_OBJECT)) return;
    const jvalue *ordinary = jget(subgroup, "ordinary_number");
    if (!ordinary || jis(ordinary, J_NULL)) return;
    int64_t group;
    const jvalue *number = jget(subgroup, "number");
    if (!number || py_int_of(number, &group)) return;
    int64_t total = section_count(&T.mag_nlabel);
    if (group < 1 || group > total) return;
    int64_t ordinary_sg = mag_ordinary_group((int)group);
    int64_t type = mag_group_type((int)group);
    if (type == 4 && section_int(&T.mag_type4_pointer, group - 1) <= 0) return;
    if (!r->has_child_sg || ordinary_sg != r->child_sg) return;
    r->magnetic_group = (int)group;
}

static int reference_space_group(int group)
{
    int64_t type = mag_group_type(group);
    int64_t sg = mag_ordinary_group(group);
    if (type == 4) sg = section_int(&T.mag_type4_supergroup, section_int(&T.mag_type4_pointer, group - 1) - 1);
    return (int)sg;
}
static int exact_presentation_origin(const md_request *r, rat *out)
{
    for (int i = 0; i < 3; i++) out[i] = rat_int(0);
    if (r->origin_kind == 0) return 0;
    if (r->origin_kind == 1) {
        const char *s = r->origin_text.data, *e = s + r->origin_text.length;
        while (s < e && (*s == ' ' || (*s >= '\t' && *s <= '\r'))) s++;
        while (e > s && (e[-1] == ' ' || (e[-1] >= '\t' && e[-1] <= '\r'))) e--;
        jbuf b = {0};
        jb_reset(&b);
        py_repr_text(&b, r->origin_text.data, r->origin_text.length);
        if (!(e - s >= 2 && *s == '(' && e[-1] == ')')) {
            md_fail(MDE_VALUE, "invalid presentation origin: %s", b.data);
            jb_free(&b);
            return MDE_VALUE;
        }
        const char *bounds[4], *c = s + 1;
        int parts = 0;
        bounds[0] = c;
        for (; c < e - 1; c++)
            if (*c == ',') {
                if (++parts < 3) bounds[parts] = c + 1;
            }
        if (parts != 2) {
            md_fail(MDE_VALUE, "invalid presentation origin: %s", b.data);
            jb_free(&b);
            return MDE_VALUE;
        }
        jb_free(&b);
        bounds[3] = e;
        for (int i = 0; i < 3; i++) {
            const char *a = bounds[i], *z = bounds[i + 1] - 1;
            if (py_fraction_of_text(a, z - a, &out[i])) {
                jbuf t = {0};
                jb_reset(&t);
                const char *x = a, *y = z;
                while (x < y && (*x == ' ' || (*x >= '\t' && *x <= '\r'))) x++;
                while (y > x && (y[-1] == ' ' || (y[-1] >= '\t' && y[-1] <= '\r'))) y--;
                py_repr_text(&t, x, y - x);
                md_fail(MDE_VALUE, "Invalid literal for Fraction: %s", t.data);
                jb_free(&t);
                return MDE_VALUE;
            }
        }
        return 0;
    }
    if (r->origin_record[3] <= 0) return md_fail(MDE_VALUE, "presentation origin denominator must be positive");
    for (int i = 0; i < 3; i++) out[i] = rat_make(r->origin_record[i], r->origin_record[3]);
    return 0;
}

/* the selected origin, exactly */
int md_origin_rational(const md_request *r, rat *out)
{
    for (int i = 0; i < 3; i++) out[i] = rat_int(0);
    if (r->origin_kind == 2) {
        if (!r->origin_record[3]) return 0;
        rat_overflow = 0;
        for (int i = 0; i < 3; i++) out[i] = rat_make(r->origin_record[i], r->origin_record[3]);
        return !rat_overflow;
    }
    if (r->origin_kind != 1) return 1;
    const char *s = r->origin_text.data, *e = s + r->origin_text.length;
    while (s < e && (*s == ' ' || (*s >= '\t' && *s <= '\r'))) s++;
    while (e > s && (e[-1] == ' ' || (e[-1] >= '\t' && e[-1] <= '\r'))) e--;
    if (!(e - s >= 2 && *s == '(' && e[-1] == ')')) return 1;
    const char *c = s + 1;
    for (int i = 0; i < 3; i++) {
        const char *start = c;
        while (c < e - 1 && *c != ',') c++;
        if ((i < 2) != (c < e - 1) || py_fraction_of_text(start, c - start, &out[i])) return 0;
        c++;
    }
    return 1;
}

void md_origin_vector(const md_request *r, double *out)
{
    out[0] = out[1] = out[2] = 0.0;
    if (r->origin_kind == 2) {
        double den = (double)r->origin_record[3];
        if (den != 0.0)
            for (int i = 0; i < 3; i++) out[i] = (double)r->origin_record[i] / den;
        return;
    }
    if (r->origin_kind != 1) return;
    const char *s = r->origin_text.data, *e = s + r->origin_text.length;
    while (s < e && (*s == ' ' || (*s >= '\t' && *s <= '\r'))) s++;
    while (e > s && (e[-1] == ' ' || (e[-1] >= '\t' && e[-1] <= '\r'))) e--;
    if (!(e - s >= 2 && *s == '(' && e[-1] == ')')) return;
    const char *bounds[4], *c = s + 1;
    int parts = 1;
    bounds[0] = c;
    for (; c < e - 1; c++)
        if (*c == ',') {
            if (parts < 3) bounds[parts] = c + 1;
            parts++;
        }
    if (parts < 3) return;
    const char *ends[3];
    for (int i = 0; i < 3; i++) {
        const char *z = bounds[i];
        while (z < e - 1 && *z != ',') z++;
        ends[i] = z;
    }
    double values[3];
    for (int i = 0; i < 3; i++) {
        rat value;
        if (py_fraction_of_text(bounds[i], ends[i] - bounds[i], &value) || py_float_of_fraction(value, &values[i])) return;
    }
    memcpy(out, values, sizeof values);
}

static void site_params(md_site *s)
{
    s->site_param_count = -1;
    if (!s->has_parent_params || !s->nparams) return;
    static const char *keys[3] = {"x", "y", "z"};
    double values[3];
    int n = 0;
    for (int index = 0; index < 3; index++) {
        int found = -1;
        for (int p = 0; p < s->nparams; p++)
            if (!strcmp(s->param_names[p], keys[index])) found = p;
        if (found >= 0) {
            values[n++] = s->param_values[found];
            continue;
        }
        int later = 0;
        for (int next = index + 1; next < 3; next++)
            for (int p = 0; p < s->nparams; p++) later |= !strcmp(s->param_names[p], keys[next]);
        if (later) values[n++] = 0.0;
    }
    while (n && fabs(values[n - 1]) < 1e-15) n--;
    if (!n) return;
    s->site_param_count = n;
    memcpy(s->site_params, values, sizeof(double) * (size_t)n);
}

static int wyckoff_row_of(int sg, const char *label, int64_t label_length, int *matches, int *row_id, int *site_pg, int *offset0)
{
    int64_t start = section_int(&T.wyckoff_pointer, sg - 1), count = section_int(&T.wyckoff_count, sg - 1);
    *matches = 0;
    for (int64_t k = 0; k < count; k++)
        if (strip_equal(&T.wyckoff_label, start + k - 1, label, label_length)) {
            (*matches)++;
            *row_id = (int)(start + k);
            *site_pg = geo_wyckoff_site_pg(start + k);
            *offset0 = (int)k;
        }
    return 0;
}

/* The parent site's Wyckoff parameters as parsed from its JSON dict. */
static int raw_params(md_site *s)
{
    const jvalue *params = jget(s->site, "wyckoff_params");
    s->has_parent_params = jis(params, J_OBJECT) && params->count;
    s->nparams = 0;
    if (!s->has_parent_params) return 0;
    for (int64_t i = 0; i < params->count && s->nparams < 3; i++) {
        const jvalue *key = params->u.members + 2 * i, *value = key + 1;
        if (key->count >= 8) continue;
        memcpy(s->param_names[s->nparams], key->u.s, (size_t)key->count);
        s->param_names[s->nparams][key->count] = 0;
        if (py_float_of(value, &s->param_values[s->nparams])) return md_fail(MDE_TYPE, "Wyckoff parameter is not a number");
        s->param_raw[s->nparams] = value;
        s->nparams++;
    }
    return 0;
}

/* site_transport._source_default_wyckoff_params */
static int default_wyckoff_params(md_request *r, md_site *s)
{
    MD_TRY(raw_params(s));
    if (!r->has_bridge) return 0;
    int64_t row_id;
    if (py_int_of(jget(s->site, "wyckoff_row_id"), &row_id)) return 0;
    const jvalue *fract = jget(s->site, "fract");
    if (!jtruthy(fract) || !jis(fract, J_ARRAY) || fract->count != 3) return 0;
    double point[3];
    for (int i = 0; i < 3; i++)
        if (float_cif_number(jat(fract, i), &point[i])) return 0;
    int default_id = (int)section_int(&T.inter_choice, r->sg - 1);
    int64_t start = section_int(&T.wyckoff_pointer, r->sg - 1), count = section_int(&T.wyckoff_count, r->sg - 1);
    int64_t offset0 = row_id - start;
    if (offset0 < 0 || offset0 >= count) return 0;
    int64_t pointer = section_int(&T.inter_wyckoff_pointer, default_id - 1);
    rat vectors[4][3], pml[3][3], origin[3];
    for (int v = 0; v < 4; v++) {
        int64_t base = (pointer - 1 + offset0) * 16 + 4 * v, den = section_int(&T.inter_wyckoff_fract, base + 3);
        for (int axis = 0; axis < 3; axis++)
            vectors[v][axis] = den ? rat_make(section_int(&T.inter_wyckoff_fract, base + axis), den) : rat_int(0);
    }
    if (geo_setting_to_cinter(r->sg, 1, default_id, pml, origin)) return md_fail(MDE_VALUE, "parent setting Source tables are inconsistent");
    rat cml[3][3], cml_origin[3];
    if (geo_setting_to_cinter(r->sg, 0, default_id, cml, cml_origin)) return md_fail(MDE_VALUE, "parent setting Source tables are inconsistent");
    double base[3], columns[9];
    for (int v = 0; v < 4; v++) {
        rat out[3];
        rat_row_mul(vectors[v], (const rat (*)[3])pml, out);
        if (!v)
            for (int axis = 0; axis < 3; axis++) out[axis] = rat_add(out[axis], cml_origin[axis]);
        for (int axis = 0; axis < 3; axis++) {
            double value;
            if (py_float_of_fraction(out[axis], &value)) return md_fail(MDE_ARITH, "Wyckoff vector out of range");
            if (!v) base[axis] = value;
            else columns[3 * (v - 1) + axis] = value;
        }
    }
    rat inverse[3][3];
    if (rat_inverse3((const rat (*)[3])r->bridge_matrix, inverse)) return md_fail(MDE_VALUE, "singular parent setting bridge");
    const jvalue *operations = jget(r->cif, "symmetry_operations");
    int64_t total = jtruthy(operations) && jis(operations, J_ARRAY) ? operations->count : 0;
    for (int64_t candidate = -1; candidate < total; candidate++) {
        double xyz[3];
        if (candidate < 0) {
            memcpy(xyz, point, sizeof xyz);
        } else {
            char *text = md_py_str(jat(operations, candidate));
            int rot[3][3], tran[3];
            char notation;
            int bad = !text || gemmi_parse_triplet(text, (int64_t)strlen(text), rot, tran, &notation) || notation == ('h' & ~3);
            free(text);
            if (bad) continue;
            gemmi_apply_to_xyz((const int (*)[3])rot, tran, point, xyz);
        }
        rat values[3], shifted[3], local[3];
        int out_of_range = 0;
        for (int axis = 0; axis < 3; axis++) out_of_range |= py_fraction_of_float(xyz[axis], &values[axis]);
        if (out_of_range) return md_fail(MDE_ARITH, "parent point out of range");
        rat_overflow = 0;
        for (int axis = 0; axis < 3; axis++) shifted[axis] = rat_sub(values[axis], r->bridge_shift[axis]);
        rat_row_mul(shifted, (const rat (*)[3])inverse, local);
        double target[3];
        for (int axis = 0; axis < 3; axis++)
            if (rat_overflow || py_float_of_fraction(rat_mod1(local[axis]), &target[axis])) return md_fail(MDE_ARITH, "parent point out of range");
        double params[3];
        int mask[3];
        int found = solve_wyckoff_formula(base, columns, target, 1e-5, 1, params, mask);
        if (found < 0) return md_fail(MDE_VALUE, "Wyckoff formula solve failed (%d)", found);
        if (!found) continue;
        static const char *names[3] = {"x", "y", "z"};
        s->has_parent_params = 1;
        s->nparams = 0;
        for (int i = 0; i < 3; i++)
            if (mask[i]) {
                strcpy(s->param_names[s->nparams], names[i]);
                s->param_raw[s->nparams] = 0;
                s->param_values[s->nparams++] = params[i];
            }
        return 0;
    }
    return 0;
}

static void write_rat_words(int64_t *out, const rat *values, int count)
{
    for (int i = 0; i < count; i++) {
        write_i128(out + 4 * i, values[i].n);
        write_i128(out + 4 * i + 2, values[i].d);
    }
}

static int64_t determinant3(const int64_t *b)
{
    int64_t d = b[0] * (b[4] * b[8] - b[5] * b[7]) - b[1] * (b[3] * b[8] - b[5] * b[6]) + b[2] * (b[3] * b[7] - b[4] * b[6]);
    return d < 0 ? -d : d;
}
static int child_layout(md_request *r, md_site *s, const rat *representative)
{
    const jvalue *row_value = jget(s->site, "wyckoff_row_id");
    int64_t row_id;
    if (exact_int(row_value, &row_id)) return md_fail(MDE_TYPE, "parent Wyckoff row id must be an exact integer");
    char *wyckoff = md_py_str(jget(s->site, "wyckoff"));
    if (!wyckoff) return md_fail(MDE_MEMORY, "out of memory");
    int matches, source_row = 0, source_pg = 0, offset0 = 0;
    wyckoff_row_of(r->sg, wyckoff, (int64_t)strlen(wyckoff), &matches, &source_row, &source_pg, &offset0);
    if (matches != 1) {
        md_fail(MDE_VALUE, "expected one Source Wyckoff row for SG%d %s, got %d", r->sg, wyckoff, matches);
        free(wyckoff);
        return MDE_VALUE;
    }
    free(wyckoff);
    rat parent_params[3], params3[3];
    static const char *names[3] = {"x", "y", "z"};
    for (int i = 0; i < 3; i++) {
        parent_params[i] = rat_int(0);
        for (int p = 0; p < (s->has_parent_params ? s->nparams : 0); p++) {
            if (strcmp(s->param_names[p], names[i])) continue;
            /* Fraction(str(value)): a JSON string or int is read exactly, a float through its repr */
            const jvalue *raw = s->param_raw[p];
            int bad;
            if (jis(raw, J_STRING)) bad = py_fraction_of_text(raw->u.s, raw->count, &parent_params[i]);
            else if (jis(raw, J_INT)) bad = (parent_params[i] = rat_int(raw->u.i), 0);
            else bad = py_fraction_of_float(s->param_values[p], &parent_params[i]);
            if (bad) return md_fail(bad == 2 ? MDE_ARITH : MDE_VALUE, "Wyckoff parameter out of range");
        }
    }
    for (int i = 0; i < 3; i++) {
        params3[i] = rat_int(0);
        if (i < s->site_param_count && py_fraction_of_float(s->site_params[i], &params3[i]))
            return md_fail(MDE_ARITH, "Wyckoff parameter out of range");
    }
    int64_t capacity = determinant3(r->source_basis), other = determinant3(r->presentation_pml);
    if (other > capacity) capacity = other;
    if (capacity < 1) capacity = 1;
    capacity *= 4 * 48;
    int64_t *sites = calloc((size_t)capacity * 5, sizeof(int64_t)), *atoms = calloc((size_t)capacity * 38, sizeof(int64_t));
    double *positions = calloc((size_t)capacity * 3, sizeof(double));
    int64_t *rows = calloc((size_t)capacity * 4 * 8, sizeof(int64_t));
    if (!sites || !atoms || !positions || !rows) {
        free(sites);
        free(atoms);
        free(positions);
        free(rows);
        return md_fail(MDE_MEMORY, "child layout ran out of memory");
    }
    int64_t parent_words[12], matrix_words[36], origin_words[12], site_words[12], representative_words[12], bridge_words[48];
    write_rat_words(parent_words, parent_params, 3);
    for (int row = 0; row < 3; row++) write_rat_words(matrix_words + 12 * row, r->child_matrix[row], 3);
    write_rat_words(origin_words, r->child_origin, 3);
    write_rat_words(site_words, params3, 3);
    write_rat_words(representative_words, representative, 3);
    if (r->has_bridge) {
        for (int row = 0; row < 3; row++) write_rat_words(bridge_words + 12 * row, r->bridge_matrix[row], 3);
        write_rat_words(bridge_words + 36, r->bridge_shift, 3);
    }
    double presentation[9], origin[3];
    for (int i = 0; i < 9; i++) presentation[i] = r->display_basis[i / 3][i % 3];
    md_origin_vector(r, origin);
    int64_t exact_words[48];
    rat exact_origin[3];
    int exact = r->has_exact_basis && md_origin_rational(r, exact_origin);
    if (exact) {
        for (int row = 0; row < 3; row++) write_rat_words(exact_words + 12 * row, r->exact_basis[row], 3);
        write_rat_words(exact_words + 36, exact_origin, 3);
    }
    int counts[3] = {0, 0, 0};
    int status = child_site_layout(r->sg, r->child_sg, (int)row_id, r->source_basis, r->source_origin, parent_words, r->setting_count,
                                   r->setting_ids, matrix_words, origin_words, r->label_memo, 256, source_row, source_pg, site_words,
                                   r->presentation_pml, representative_words, r->has_bridge ? bridge_words : 0, exact ? exact_words : 0, presentation, origin,
                                   r->centering_letter, r->point_operation_count, r->magnetic_group, (int)capacity, sites, atoms,
                                   positions, rows, counts);
    free(rows);
    if (status) {
        free(sites);
        free(atoms);
        free(positions);
        return md_kernel(MDE_CHILD_LAYOUT, status);
    }
    int site_count = counts[0], atom_count = counts[1];
    s->presentation_rows = counts[2];
    s->nchild = site_count;
    s->natoms = atom_count;
    s->children = calloc((size_t)site_count + 1, sizeof(md_child_site));
    s->atom_ids = calloc((size_t)atom_count + 1, sizeof(char *));
    s->atom_site = calloc((size_t)atom_count + 1, sizeof(int));
    s->member_order = calloc((size_t)atom_count + 1, sizeof(int64_t));
    s->topology = calloc((size_t)atom_count * 10 + 1, sizeof(int64_t));
    s->positions = positions;
    int failed = !s->children || !s->atom_ids || !s->atom_site || !s->member_order || !s->topology;
    int reference = r->magnetic_group ? reference_space_group(r->magnetic_group) : r->child_sg;
    int64_t child_start = section_int(&T.wyckoff_pointer, reference - 1);
    for (int a = 0; a < atom_count && !failed; a++) {
        const int64_t *t = atoms + 38 * (size_t)a;
        s->atom_ids[a] = md_format("%s:%lld:%lld", s->label_prefix, (long long)t[0], (long long)t[1]);
        s->member_order[a] = t[12];
        memcpy(s->topology + 10 * (size_t)a, t + 3, 9 * sizeof(int64_t));
        s->topology[10 * (size_t)a + 9] = t[1];
        failed |= !s->atom_ids[a];
    }
    for (int k = 0; k < site_count && !failed; k++) {
        const int64_t *row = sites + 5 * (size_t)k;
        md_child_site *c = s->children + k;
        int64_t offset = row[0], label_index = r->magnetic_group ? offset - 1 : offset;
        md_str label = section_text(&T.wyckoff_label, child_start + label_index - 1);
        c->id = md_format("%s_%d", s->label_prefix, k + 1);
        c->wyckoff_site = md_format("%lld%.*s", (long long)row[1], (int)label.length, label.data);
        const int64_t *words = atoms + 38 * (size_t)row[2] + 25;
        for (int axis = 0; axis < 3; axis++) c->representative[axis] = (rat){read_i128(words + 4 * axis), read_i128(words + 4 * axis + 2)};
        c->first = (int)row[3];
        c->count = (int)row[4];
        for (int a = c->first; a < c->first + c->count; a++) s->atom_site[a] = k;
        failed |= !c->id || !c->wyckoff_site;
    }
    free(sites);
    free(atoms);
    if (failed) return md_fail(MDE_MEMORY, "out of memory");
    return 0;
}

void md_selected(md_request *r)
{
    r->point_operation_count = geo_point_op_count();
    selected_subgroup_number(r);
    display_basis(r);
    selected_origin(r);
    selected_magnetic_group(r);
}

int md_child_structure(md_request *r)
{
    MD_TRY(md_tables());
    if (!r->has_child_sg) return md_fail(MDE_VALUE, "selected subgroup has no ordinary space-group identity");
    if (!r->has_display_basis) return md_fail(MDE_VALUE, "selected subgroup has no presentation basis");
    MD_TRY(md_parent_setting_bridge(r->sg, r->has_setting, r->setting_id, &r->has_bridge, r->bridge_matrix, r->bridge_shift));
    MD_TRY(basis_cinter_to_pml(r));
    if (!r->has_presentation_pml) return md_fail(MDE_VALUE, "selected subgroup has no exact presentation basis");
    if (!r->iso) return md_fail(MDE_VALUE, "selected OPD has no Source subgroup transform");
    MD_TRY(exact_integer_tuple(jget(r->iso, "source_basis_values"), 9, "Source subgroup basis", r->source_basis));
    MD_TRY(exact_integer_tuple(jget(r->iso, "source_origin_values"), 4, "Source subgroup origin", r->source_origin));
    if (r->source_origin[3] <= 0) return md_fail(MDE_VALUE, "Source subgroup origin denominator must be positive");
    if (r->child_sg < 1 || r->child_sg > 230) return md_fail(MDE_VALUE, "child space group must be between 1 and 230");
    md_str symbol = section_text(&T.label_short, section_int(&T.inter_choice, r->child_sg - 1) - 1);
    r->centering_letter = symbol.length ? (char)(symbol.data[0] & ~0x20) : 0;
    rat pml[3][3];
    rat unused[3][3];
    if (geo_setting_to_cinter(r->sg, 1, r->has_setting ? (int)r->setting_id : 0, pml, r->parent_pml_origin) ||
        rat_inverse3((const rat (*)[3])pml, r->parent_pml_inverse) ||
        geo_setting_to_cinter(r->sg, 0, r->has_setting ? (int)r->setting_id : 0, unused, r->parent_pml_origin))
        return md_fail(MDE_VALUE, "parent setting Source tables are inconsistent");
    rat presentation_origin[3];
    MD_TRY(exact_presentation_origin(r, presentation_origin));
    int64_t presentation_words[12], matrix_words[36], origin_words[12];
    write_rat_words(presentation_words, presentation_origin, 3);
    int count = child_presentation_setting(r->sg, r->child_sg, r->has_setting ? (int)r->setting_id : 0, r->source_basis, r->source_origin,
                                           r->presentation_pml, presentation_words, 32, r->setting_ids, matrix_words, origin_words);
    if (count < 0) return md_kernel(MDE_CHILD_SETTING, -count);
    r->setting_count = count;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
            r->child_matrix[row][col] = (rat){read_i128(matrix_words + 12 * row + 4 * col), read_i128(matrix_words + 12 * row + 4 * col + 2)};
    for (int axis = 0; axis < 3; axis++) r->child_origin[axis] = (rat){read_i128(origin_words + 4 * axis), read_i128(origin_words + 4 * axis + 2)};
    r->label_memo = calloc(1 + 3 * 256, sizeof(int64_t));
    if (!r->label_memo) return md_fail(MDE_MEMORY, "out of memory");
    for (int index = 0; index < r->nsites; index++) {
        md_site *s = r->sites + index;
        MD_TRY(default_wyckoff_params(r, s));
        site_params(s);
        const jvalue *matched = jget(s->site, "representative_matched_fract");
        if (!matched) return md_fail(MDE_KEY, "'representative_matched_fract'");
        rat shifted[3], representative[3];
        for (int axis = 0; axis < 3; axis++) {
            const jvalue *value = jat(matched, axis);
            if (!value) return md_fail(MDE_INDEX, "list index out of range");
            char *text = md_py_str(value);
            rat parsed;
            int bad = !text || py_fraction_of_text(text, (int64_t)strlen(text), &parsed);
            free(text);
            if (bad) return md_fail(MDE_VALUE, "invalid representative coordinate");
            shifted[axis] = rat_sub(parsed, r->parent_pml_origin[axis]);
        }
        rat_row_mul(shifted, (const rat (*)[3])r->parent_pml_inverse, representative);
        MD_TRY(child_layout(r, s, representative));
        for (int other = 0; other < index; other++)
            for (int a = 0; a < s->natoms; a++)
                for (int b = 0; b < r->sites[other].natoms; b++)
                    if (!strcmp(s->atom_ids[a], r->sites[other].atom_ids[b]))
                        return md_fail(MDE_VALUE, "canonical atom identities cross parent sites: ['%s']", s->atom_ids[a]);
    }
    return 0;
}

void md_release(md_request *r)
{
    for (int i = 0; i < r->nspecs; i++) md_spec_free(r->specs + i);
    free(r->specs);
    for (int i = 0; i < r->nsites; i++) {
        md_site *s = r->sites + i;
        free(s->atom_label);
        free(s->label_prefix);
        for (int k = 0; k < s->nchild; k++) {
            free(s->children[k].id);
            free(s->children[k].wyckoff_site);
        }
        free(s->children);
        for (int a = 0; a < s->natoms; a++) free(s->atom_ids[a]);
        free(s->atom_ids);
        free(s->atom_site);
        free(s->member_order);
        free(s->topology);
        free(s->positions);
    }
    free(r->sites);
    free(r->label_memo);
    free(r->parent_symbol);
    free(r->display_kvector);
    free(r->opd_direction);
    memset(r, 0, sizeof *r);
}

static void put_rat_float(jbuf *b, rat value)
{
    double out;
    if (py_float_of_fraction(value, &out)) b->failed = 1;
    else jb_float(b, out);
}

/* The child structure as the dev comparison reads it. */
void md_child_atoms_json(const md_request *r, jbuf *b)
{
    jb_lit(b, "[");
    int first = 1;
    for (int i = 0; i < r->nsites; i++) {
        const md_site *s = r->sites + i;
        for (int k = 0; k < s->nchild; k++) {
            const md_child_site *c = s->children + k;
            if (!first) jb_lit(b, ", ");
            first = 0;
            jb_lit(b, "{\"label\": ");
            jb_cstring(b, c->id);
            jb_lit(b, ", \"child_site\": ");
            jb_cstring(b, c->id);
            jb_lit(b, ", \"site\": ");
            jb_cstring(b, c->wyckoff_site);
            jb_lit(b, ", \"multiplicity\": ");
            jb_int(b, c->count);
            jb_lit(b, ", \"atom_ids\": [");
            for (int a = c->first; a < c->first + c->count; a++) {
                if (a > c->first) jb_lit(b, ", ");
                jb_cstring(b, s->atom_ids[a]);
            }
            jb_lit(b, "], \"xyz\": [");
            for (int axis = 0; axis < 3; axis++) {
                if (axis) jb_lit(b, ", ");
                put_rat_float(b, c->representative[axis]);
            }
            jb_lit(b, "], \"_presentation_orbit_points\": [");
            for (int a = c->first; a < c->first + c->count; a++) {
                if (a > c->first) jb_lit(b, ", ");
                jb_lit(b, "[");
                for (int axis = 0; axis < 3; axis++) {
                    if (axis) jb_lit(b, ", ");
                    jb_float(b, s->positions[3 * (size_t)a + axis]);
                }
                jb_lit(b, "]");
            }
            jb_lit(b, "]}");
        }
    }
    jb_lit(b, "]");
}
