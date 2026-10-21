/* isotropy.catalog._dynamic_opd_rows: the get_isotropy_ rows of a
 * parametric gid in their displayed settings. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"
#include "sel.h"
#include "tables.h"

int isotropy_rows_for(int gid, int param_count, const int64_t *kparam, int magnetic, int occurrence);
int64_t isotropy_rows_ints(const int64_t **data);
int64_t isotropy_rows_doubles(const double **data);
const char *k_display_kvector(int sg, int kslot);
int k_source_kparam_of_visible(int sg, int kslot, const int64_t *visible_words, int64_t *out);

static struct {
    section_view irr_sg, irr_k, irr_full_dim, k_dim;
    int ready;
} T;

static int tables(void)
{
    if (T.ready) return 0;
    if (sel_tables() || sel_mag_tables() || section_find("little.little_irr_space_group", &T.irr_sg) || section_find("little.little_irr_k", &T.irr_k) ||
        section_find("little.little_irr_full_dim", &T.irr_full_dim) || section_find("little.little_k_dim", &T.k_dim))
        return 1;
    T.ready = 1;
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

static int strip_bounds(const char *raw, const char **a, const char **e)
{
    *a = raw;
    *e = raw + strlen(raw);
    while (*a < *e && sel_is_space(**a)) (*a)++;
    while (*e > *a && sel_is_space((*e)[-1])) (*e)--;
    return *a < *e;
}

/* catalog._dynamic_kparam: 1 with out, 0 for None */
static int dynamic_kparam(int64_t gid, const sel_params *p, int64_t *out)
{
    if (!p->count) return 0;
    static const char *keys[3] = {"_source_kparam", "source_kparam", "kparam"};
    for (int k = 0; k < 3; k++) {
        const char *raw = sel_param(p, keys[k]);
        if (!raw) continue;
        int64_t values[4];
        int n = 0;
        for (const char *c = raw; *c && n < 4;) {
            int negative = *c == '-' && c[1] >= '0' && c[1] <= '9';
            if (negative || (*c >= '0' && *c <= '9')) {
                const char *d = c + negative;
                int64_t v = 0;
                while (*d >= '0' && *d <= '9') v = v * 10 + (*d++ - '0');
                values[n++] = negative ? -v : v;
                c = d;
            } else {
                c++;
            }
        }
        if (n >= 4 && values[3]) {
            memcpy(out, values, sizeof values);
            return 1;
        }
    }
    int sg = (int)section_int(&T.irr_sg, gid - 1), kslot = (int)section_int(&T.irr_k, gid - 1);
    int64_t dim = section_int(&T.k_dim, (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1);
    if (dim > 0) {
        int64_t visible[6];
        int has = sel_visible_k_vector(sg, kslot, p, visible);
        if (has < 0) {
            md_fail(-has, "k-vector Source tables are inconsistent");
            return -1;
        }
        if (has) {
            int found = k_source_kparam_of_visible(sg, kslot, visible, out);
            if (found < 0) {
            md_fail(found == -2 ? MDE_ARITH : MDE_VALUE, "k-vector Source tables are inconsistent");
            return -1;
        }
            if (found) return 1;
        }
    }
    const char *display = k_display_kvector(sg, kslot);
    if (!display) {
            md_fail(MDE_VALUE, "k-vector Source tables are inconsistent");
            return -1;
        }
    char display_keys[3];
    int nkeys = 0;
    size_t length = strlen(display);
    for (size_t i = 0; i < length; i++) {
        char c = display[i];
        if (c != 'a' && c != 'b' && c != 'g') continue;
        int letter_before = i && ((display[i - 1] >= 'a' && display[i - 1] <= 'z') || (display[i - 1] >= 'A' && display[i - 1] <= 'Z'));
        int letter_after = i + 1 < length && ((display[i + 1] >= 'a' && display[i + 1] <= 'z') || (display[i + 1] >= 'A' && display[i + 1] <= 'Z'));
        if (letter_before || letter_after) continue;
        int seen = 0;
        for (int k = 0; k < nkeys; k++) seen |= display_keys[k] == c;
        if (!seen && nkeys < 3) display_keys[nkeys++] = c;
    }
    static const char names[3] = {'a', 'b', 'g'};
    rat values[3] = {rat_int(0), rat_int(0), rat_int(0)};
    for (int i = 0; i < (dim < 3 ? dim : 3); i++) {
        char key[2] = {i < nkeys ? display_keys[i] : names[i], 0}, fallback[2] = {names[i], 0};
        const char *raw = sel_param(p, key);
        if (!raw) raw = sel_param(p, fallback);
        const char *a, *e;
        if (!raw || !strip_bounds(raw, &a, &e)) continue;
        if (py_fraction_of_text(a, e - a, &values[i])) {
            md_fail(MDE_VALUE, "Invalid literal for Fraction: '%.*s'", (int)(e - a), a);
            return -1;
        }
    }
    int64_t denominator = 1;
    for (int i = 0; i < 3; i++) {
        int64_t d = (int64_t)values[i].d;
        denominator = denominator / gcd64(denominator, d) * d;
    }
    for (int i = 0; i < 3; i++) out[i] = (int64_t)(values[i].n * (denominator / (int64_t)values[i].d));
    out[3] = denominator;
    return 1;
}

int sel_dynamic_kparam(int64_t gid, const sel_params *p, int64_t *out)
{
    if (tables()) return -md_fail(MDE_VALUE, "isotropy Source tables are inconsistent");
    int found = dynamic_kparam(gid, p, out);
    return found < 0 ? -md_err.kind : found;
}

/* get_isotropy._fmt_matrix_value */
static char *matrix_value(double value)
{
    char text[64];
    if (fabs(value) < 5e-7) value = 0.0;
    double rounded = nearbyint(value);
    if (fabs(value - rounded) < 5e-7) {
        snprintf(text, sizeof text, "%lld", (long long)rounded);
    } else {
        snprintf(text, sizeof text, "%.6f", value);
        size_t n = strlen(text);
        while (n && text[n - 1] == '0') n--;
        if (n && text[n - 1] == '.') n--;
        text[n] = 0;
    }
    return (char *)md_strdup(text, (int64_t)strlen(text));
}

typedef struct {
    int64_t free, subgroup, basis[9], origin[4], index;
    int64_t *flat, *kept, nkept;
    char ***matrix;
    int *widths;
    char label[16];
} dyn_row;

static void label_of(int64_t free, int64_t index, char *out)
{
    if (free == 1) snprintf(out, 16, "P%lld", (long long)index);
    else if (free == 2) snprintf(out, 16, "C%lld", (long long)index);
    else if (free == 3) snprintf(out, 16, "S%lld", (long long)index);
    else snprintf(out, 16, "%lldD%lld", (long long)free, (long long)index);
}

static int generate(int64_t gid, const int64_t *kparam, int magnetic, int occurrence, dyn_row **rows_out, int *count_out)
{
    *rows_out = 0;
    *count_out = 0;
    int status = isotropy_rows_for((int)gid, 4, kparam, magnetic, occurrence);
    if (status) {
        md_err.kind = MDE_ISOTROPY_ROWS;
        md_err.status = status;
        md_err.message[0] = 0;
        return MDE_ISOTROPY_ROWS;
    }
    const int64_t *ints;
    const double *doubles;
    int64_t nints = isotropy_rows_ints(&ints);
    isotropy_rows_doubles(&doubles);
    int64_t dim = section_int(&T.irr_full_dim, gid - 1);
    int cap = 0, n = 0;
    dyn_row *rows = 0;
    int64_t counts[100] = {0};
    int64_t cursor = 0, mcursor = 0;
    while (cursor < nints) {
        int64_t free = ints[cursor], length = ints[cursor + 1];
        const int64_t *payload = ints + cursor + 2;
        cursor += 2 + length;
        if (free < 1 || free > 99) return md_fail(MDE_VALUE, "invalid order-parameter dimension: %lld", (long long)free);
        if (n == cap) {
            cap = cap ? 2 * cap : 16;
            dyn_row *grown = realloc(rows, (size_t)cap * sizeof *rows);
            if (!grown) return md_fail(MDE_MEMORY, "out of memory");
            rows = grown;
        }
        dyn_row *r = rows + n++;
        memset(r, 0, sizeof *r);
        *rows_out = rows;
        *count_out = n;
        r->free = free;
        r->index = ++counts[free];
        label_of(free, r->index, r->label);
        if (magnetic) {
            r->subgroup = payload[0];
            memcpy(r->basis, payload + 2, sizeof r->basis);
            memcpy(r->origin, payload + 11, sizeof r->origin);
            r->flat = malloc(sizeof(int64_t) * (size_t)(length ? length : 1));
            r->nkept = ints[cursor];
            r->kept = malloc(sizeof(int64_t) * 5 * (size_t)(r->nkept ? r->nkept : 1));
            if (!r->flat || !r->kept) return md_fail(MDE_MEMORY, "out of memory");
            memcpy(r->flat, payload, sizeof(int64_t) * (size_t)length);
            memcpy(r->kept, ints + cursor + 1, sizeof(int64_t) * 5 * (size_t)r->nkept);
            cursor += 1 + 5 * r->nkept;
        } else {
            r->subgroup = payload[0];
            memcpy(r->basis, payload + 1, sizeof r->basis);
            memcpy(r->origin, payload + 10, sizeof r->origin);
        }
        r->matrix = calloc((size_t)free, sizeof(char **));
        r->widths = calloc((size_t)free, sizeof(int));
        if (!r->matrix || !r->widths) return md_fail(MDE_MEMORY, "out of memory");
        for (int64_t row = 0; row < free; row++) {
            r->matrix[row] = calloc((size_t)dim + 1, sizeof(char *));
            if (!r->matrix[row]) return md_fail(MDE_MEMORY, "out of memory");
            r->widths[row] = (int)dim;
            for (int64_t col = 0; col < dim; col++) r->matrix[row][col] = matrix_value(doubles[mcursor + row * dim + col]);
        }
        mcursor += free * dim;
    }
    *rows_out = rows;
    *count_out = n;
    return 0;
}

static void rows_free(dyn_row *rows, int n)
{
    for (int i = 0; i < n; i++) {
        for (int64_t r = 0; r < rows[i].free && rows[i].matrix; r++) {
            for (int c = 0; rows[i].matrix[r] && c < rows[i].widths[r]; c++) free(rows[i].matrix[r][c]);
            free(rows[i].matrix[r]);
        }
        free(rows[i].matrix);
        free(rows[i].widths);
        free(rows[i].flat);
        free(rows[i].kept);
    }
    free(rows);
}

static int by_file_order(const void *a, const void *b)
{
    const dyn_row *x = a, *y = b;
    if (x->free != y->free) return x->free < y->free ? -1 : 1;
    return (x->index > y->index) - (x->index < y->index);
}

static void put_matrix(jbuf *b, const dyn_row *r)
{
    jb_lit(b, "[");
    for (int64_t row = 0; row < r->free; row++) {
        if (row) jb_lit(b, ", ");
        jb_lit(b, "[");
        for (int c = 0; c < r->widths[row]; c++) {
            if (c) jb_lit(b, ", ");
            jb_cstring(b, r->matrix[row][c]);
        }
        jb_lit(b, "]");
    }
    jb_lit(b, "]");
}

static int dynamic_row(jbuf *b, int magnetic, int parent_sg, int64_t gid, int64_t index, const dyn_row *r, const int64_t *kparam,
                       const sel_params *active_params, int has_setting, int64_t setting_id, int has_cell, const double *cell)
{
    int group = (int)r->subgroup;
    int subgroup = magnetic ? sel_mag_ordinary(group) : group;
    int parent_group = magnetic ? mag_orderparam_group((int)section_int(&T.irr_sg, gid - 1)) : 0;
    sel_setting displayed;
    rat rows[3][3];
    MD_TRY(sel_display(parent_sg, subgroup, r->basis, r->origin, has_setting, setting_id, has_cell, cell, 1, 0, magnetic ? group : 0,
                       &displayed, rows));
    const int64_t *m = r->basis;
    int64_t det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    if (det < 0) det = -det;
    jbuf opd = {0}, display = {0};
    jb_reset(&opd);
    jb_reset(&display);
    if (r->free) sel_opd_from_rows(&opd, r->matrix, (int)r->free, r->widths);
    else jb_lit(&opd, "()");
    sel_display_opd(&display, opd.data);
    strings active = {0};
    int64_t full_dim = section_int(&T.irr_full_dim, gid - 1);
    int status = sel_active_k_vectors(parent_sg, gid, full_dim, r->matrix, (int)r->free, r->widths, active_params, &active);
    if (status) {
        jb_free(&opd);
        jb_free(&display);
        return status;
    }
    jb_lit(b, "{\"direction\": {\"id\": 0, \"label\": ");
    jb_cstring(b, r->label);
    jb_lit(b, ", \"dimension\": ");
    jb_int(b, r->free);
    jb_lit(b, ", \"subduction\": 0, \"opd\": ");
    jb_string(b, display.data, display.used);
    jb_lit(b, ", \"vectors\": ");
    put_matrix(b, r);
    jb_lit(b, magnetic ? ", \"dynamic\": true, \"magnetic\": true}, \"isotropy\": {\"row_id\": " : ", \"dynamic\": true}, \"isotropy\": {\"row_id\": ");
    jb_int(b, index);
    jb_lit(b, ", \"dynamic\": true, \"dynamic_gid\": ");
    jb_int(b, gid);
    if (magnetic) jb_lit(b, ", \"magnetic\": true");
    jb_lit(b, ", \"opd_label\": ");
    jb_cstring(b, r->label);
    jb_lit(b, ", \"direction_label\": ");
    jb_cstring(b, r->label);
    const char *opd_keys[3] = {", \"direction_opd\": ", ", \"stored_opd\": ", ", \"display_opd\": "};
    for (int k = 0; k < 3; k++) {
        jb_lit(b, opd_keys[k]);
        if (k == 2) jb_string(b, display.data, display.used);
        else jb_string(b, opd.data, opd.used);
    }
    jb_lit(b, ", \"source_opd\": ");
    jb_string(b, opd.data, opd.used);
    jb_lit(b, ", \"source_matrix\": ");
    put_matrix(b, r);
    jb_lit(b, ", \"source_display_rows\": ");
    put_matrix(b, r);
    jb_lit(b, ", \"source_kparam\": [");
    for (int i = 0; i < 4; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, kparam[i]);
    }
    jb_lit(b, "], \"op_rows\": ");
    jb_int(b, r->free);
    jb_lit(b, ", \"free\": ");
    jb_int(b, r->free);
    jb_lit(b, ", \"arms\": ");
    jb_int(b, det);
    jb_lit(b, ", \"s\": ");
    jb_int(b, det);
    jb_lit(b, ", \"i\": ");
    int has_index;
    int64_t subgroup_index = magnetic ? sel_mag_index(parent_group, group, det, &has_index) : sel_subgroup_index(parent_sg, subgroup, det, &has_index);
    if (has_index) jb_int(b, subgroup_index);
    else jb_lit(b, "null");
    jb_lit(b, ", \"k_active_vectors\": ");
    sel_put_strings(b, &active);
    jb_lit(b, ", \"k_active\": ");
    sel_put_joined(b, &active);
    jbuf subgroup_json = {0}, extra = {0};
    if (magnetic) {
        jb_reset(&subgroup_json);
        jb_reset(&extra);
        sel_mag_subgroup_json(&subgroup_json, group);
        if (r->nkept) {
            jb_lit(&extra, ", \"source_operation_records\": ");
            sel_put_records(&extra, r->kept, r->nkept);
        }
    }
    sel_put_subgroup(b, parent_sg, subgroup, r->basis, r->origin, &displayed, (const rat (*)[3])rows, det, magnetic ? &subgroup_json : 0,
                     magnetic ? &extra : 0);
    jb_free(&subgroup_json);
    jb_free(&extra);
    jb_lit(b, ", \"origin_values\": [");
    for (int i = 0; i < 4; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, r->origin[i]);
    }
    if (magnetic) {
        jb_lit(b, "], \"ferroic\": null, \"ferroic_properties\": ");
        status = sel_mag_ferroic(b, group);
        jb_lit(b, "}}");
    } else {
        jb_lit(b, "], \"ferroic\": null, \"ferroic_properties\": []}}");
    }
    strings_free(&active);
    jb_free(&opd);
    jb_free(&display);
    return status;
}

/* catalog._dynamic_point_occurrence */
static int point_occurrence(const jvalue *selector)
{
    if (!jis(selector, J_STRING)) return 0;
    const char *a = selector->u.s, *e = a + selector->count;
    while (a < e && sel_is_space(*a)) a++;
    while (e > a && sel_is_space(e[-1])) e--;
    if (e - a < 2 || *a != 'P' || a[1] < '1' || a[1] > '9') return 0;
    int64_t value = 0;
    for (const char *c = a + 1; c < e; c++) {
        if (*c < '0' || *c > '9') return 0;
        value = value * 10 + (*c - '0');
        if (value > 1000000) return 0;
    }
    return (int)value;
}

int sel_dynamic_rows(int magnetic, int parent_sg, int64_t gid, const sel_params *k_params, const sel_params *display, const jvalue *selected_orderparam,
                     int has_setting, int64_t setting_id, int has_cell, const double *cell, jbuf *out)
{
    if (tables()) return md_fail(MDE_VALUE, "isotropy Source tables are inconsistent");
    int64_t kparam[4];
    int found = dynamic_kparam(gid, k_params, kparam);
    if (found < 0) return md_err.kind;
    if (!found) {
        jb_lit(out, "[]");
        return 0;
    }
    const sel_params *active_params = display ? display : k_params;
    int occurrence = point_occurrence(selected_orderparam);
    dyn_row *rows = 0;
    int n = 0, status = 0, selected = -1;
    if (occurrence) {
        status = generate(gid, kparam, magnetic, occurrence, &rows, &n);
        char label[16];
        label_of(1, occurrence, label);
        for (int i = 0; !status && i < n && selected < 0; i++)
            if (rows[i].free == 1 && !strcmp(rows[i].label, label)) selected = i;
        if (!status && selected < 0) {
            rows_free(rows, n);
            rows = 0;
            n = 0;
        }
    }
    if (!status && selected < 0) {
        status = generate(gid, kparam, magnetic, 0, &rows, &n);
        if (!status) qsort(rows, (size_t)n, sizeof *rows, by_file_order);
    }
    jb_lit(out, "[");
    int first = 1;
    for (int i = 0; !status && i < n; i++) {
        if (selected >= 0 && i != selected) continue;
        if (!first) jb_lit(out, ", ");
        first = 0;
        status = dynamic_row(out, magnetic, parent_sg, gid, selected >= 0 ? occurrence : i + 1, rows + i, kparam, active_params, has_setting,
                             setting_id, has_cell, cell);
    }
    jb_lit(out, "]");
    rows_free(rows, n);
    return status;
}
