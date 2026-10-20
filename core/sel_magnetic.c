/* isotropy.catalog.magnetic_opd_rows for a Source irrep: the data_magnetic
 * isotropy rows, their fixed BNS selection and ferromagnetism. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"
#include "sel.h"
#include "tables.h"

int mag_isotropy_row_dim(int64_t row_id);
int mag_isotropy_row_free(int64_t row_id);

void exact_register_codes(int rows, int cols, const double *doubles, const int64_t *codes);

int magnetic_records(int group, int setting, int capacity, int64_t *records_out);
int magnetic_records_iso(int group, int capacity, int64_t *records_out);
int orderparam_to_subgroup_rows(int gid, int param_count, const int64_t *kparam, const double *values, int64_t value_count, int row_count,
                                int candidate_count, const int64_t *candidates, const int *candidate_matrix_ops, const int *candidate_signs,
                                int max_iter, int64_t *basis_out, int64_t *kept_out, int *kept_count);
int id_subgroup_magnetic_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t capacity,
                                  int64_t *out);

static struct {
    section_view point_group, point_group_order, op_r, mag2nonmag, bns, nlabel, space_group, iso_pointer, iso_parent, iso_subgroup,
        iso_subgroup_old, iso_basis, iso_origin, iso_label, iso_op_pointer, iso_op, point_op, const_values, irr_sg;
    int ready;
} T;

int sel_mag_tables(void)
{
    if (T.ready) return 0;
    struct {
        const char *name;
        section_view *view;
    } list[] = {{"magnetic.mag_point_group", &T.point_group},
                {"magnetic.mag_point_group_order", &T.point_group_order},
                {"magnetic.mag_point_op_r", &T.op_r},
                {"magnetic.mag_point_op_mag2nonmag", &T.mag2nonmag},
                {"magnetic.mag_bns_label", &T.bns},
                {"magnetic.mag_nlabel", &T.nlabel},
                {"magnetic.mag_space_group", &T.space_group},
                {"magnetic.mag_iso_irrep_pointer", &T.iso_pointer},
                {"magnetic.mag_iso_parent", &T.iso_parent},
                {"magnetic.mag_iso_subgroup", &T.iso_subgroup},
                {"magnetic.mag_iso_subgroup_old", &T.iso_subgroup_old},
                {"magnetic.mag_iso_basis", &T.iso_basis},
                {"magnetic.mag_iso_origin", &T.iso_origin},
                {"magnetic.mag_iso_orderparam_label", &T.iso_label},
                {"magnetic.mag_iso_orderparam_pointer", &T.iso_op_pointer},
                {"magnetic.mag_iso_orderparam", &T.iso_op},
                {"space.ipoint_op", &T.point_op},
                {"const.values", &T.const_values},
                {"little.little_irr_space_group", &T.irr_sg}};
    if (sel_tables() || mag_tables()) return 1;
    for (size_t i = 0; i < sizeof list / sizeof list[0]; i++)
        if (section_find(list[i].name, list[i].view)) return 1;
    T.ready = 1;
    return 0;
}

static void put_label(jbuf *b, const section_view *view, int64_t index)
{
    const char *s;
    int64_t n;
    if (section_string(view, index, &s, &n)) {
        jb_lit(b, "\"\"");
        return;
    }
    while (n && (*s == ' ' || *s == '\t')) {
        s++;
        n--;
    }
    while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    jb_string(b, s, (size_t)n);
}

int sel_mag_ordinary(int group) { return (int)section_int(&T.space_group, group - 1); }

/* The magnetic subgroup dict {number, symbol, ordinary_number, display_label}. */
void sel_mag_subgroup_json(jbuf *b, int group)
{
    jb_lit(b, "{\"number\": ");
    jb_int(b, group);
    jb_lit(b, ", \"symbol\": ");
    put_label(b, &T.bns, group - 1);
    jb_lit(b, ", \"ordinary_number\": ");
    jb_int(b, sel_mag_ordinary(group));
    jb_lit(b, ", \"display_label\": ");
    put_label(b, &T.nlabel, group - 1);
    jb_lit(b, "}");
}

/* catalog._magnetic_subgroup_index_for_groups */
int64_t sel_mag_index(int parent_group, int subgroup_group, int64_t size, int *has)
{
    int64_t parent_pg = section_int(&T.point_group, parent_group - 1), sub_pg = section_int(&T.point_group, subgroup_group - 1);
    int64_t parent_order = section_int(&T.point_group_order, parent_pg - 1), sub_order = section_int(&T.point_group_order, sub_pg - 1);
    int64_t value = size * parent_order;
    *has = sub_order > 0 && value % sub_order == 0;
    return *has ? value / sub_order : 0;
}

static int rank3(int64_t (*m)[3], int rows)
{
    int rank = 0;
    for (int col = 0; col < 3 && rank < rows; col++) {
        int pivot = -1;
        for (int r = rank; r < rows; r++)
            if (m[r][col]) {
                pivot = r;
                break;
            }
        if (pivot < 0) continue;
        for (int c = 0; c < 3; c++) {
            int64_t t = m[rank][c];
            m[rank][c] = m[pivot][c];
            m[pivot][c] = t;
        }
        for (int r = 0; r < rows; r++) {
            if (r == rank || !m[r][col]) continue;
            int64_t a = m[rank][col], f = m[r][col];
            int64_t g = 0;
            for (int c = 0; c < 3; c++) {
                m[r][c] = m[r][c] * a - f * m[rank][c];
                int64_t v = m[r][c] < 0 ? -m[r][c] : m[r][c];
                while (v) {
                    int64_t t = g % v;
                    g = v;
                    v = t;
                }
            }
            if (g > 1)
                for (int c = 0; c < 3; c++) m[r][c] /= g;
        }
        rank++;
    }
    return rank;
}

/* catalog.magnetic_ferroic_properties into b */
int sel_mag_ferroic(jbuf *b, int group)
{
    static signed char cache[256];
    int64_t total = section_count(&T.point_group);
    if (group < 1 || group > total) return md_fail(MDE_VALUE, "unknown magnetic group: %d", group);
    int64_t pg = section_int(&T.point_group, group - 1);
    if (pg < 1 || pg > 255) return md_fail(MDE_VALUE, "unknown magnetic point group: %lld", (long long)pg);
    if (!cache[pg]) {
        int representative = 0;
        for (int64_t g = 0; g < total && !representative; g++)
            if (section_int(&T.point_group, g) == pg) representative = (int)g + 1;
        if (!representative) return md_fail(MDE_VALUE, "unknown magnetic point group: %lld", (long long)pg);
        int64_t records[192][5];
        int count = magnetic_records(representative, 1, 192, &records[0][0]);
        if (count < 0) {
            md_err.kind = MDE_MAGNETIC;
            md_err.status = -count;
            md_err.message[0] = 0;
            return MDE_MAGNETIC;
        }
        int64_t expected = section_int(&T.point_group_order, pg - 1);
        if (count != expected)
            return md_fail(MDE_VALUE, "magnetic point-group closure is incomplete: group=%d expected=%lld actual=%d", representative,
                           (long long)expected, count);
        int64_t (*constraints)[3] = malloc(sizeof(int64_t) * 3 * 3 * (size_t)(count ? count : 1));
        if (!constraints) return md_fail(MDE_MEMORY, "out of memory");
        for (int i = 0; i < count; i++) {
            int op = (int)records[i][4];
            int64_t nonmag = section_int(&T.mag2nonmag, op - 1);
            int64_t r[9];
            for (int k = 0; k < 9; k++) r[k] = section_int(&T.point_op, (nonmag - 1) * 9 + k);
            int64_t det = r[0] * (r[4] * r[8] - r[5] * r[7]) - r[1] * (r[3] * r[8] - r[5] * r[6]) + r[2] * (r[3] * r[7] - r[4] * r[6]);
            int64_t axial = (section_int(&T.op_r, op - 1) ? -1 : 1) * det;
            for (int row = 0; row < 3; row++)
                for (int col = 0; col < 3; col++) constraints[3 * i + row][col] = axial * r[3 * row + col] - (row == col);
        }
        cache[pg] = rank3(constraints, 3 * count) < 3 ? 1 : 2;
        free(constraints);
    }
    jb_lit(b, cache[pg] == 1 ? "[\"ferromagnetic\"]" : "[]");
    return 0;
}

static void put_ints(jbuf *b, const int64_t *values, int64_t n)
{
    jb_lit(b, "[");
    for (int64_t i = 0; i < n; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, values[i]);
    }
    jb_lit(b, "]");
}

void sel_put_records(jbuf *b, const int64_t *records, int64_t count)
{
    jb_lit(b, "[");
    for (int64_t i = 0; i < count; i++) {
        if (i) jb_lit(b, ", ");
        put_ints(b, records + 5 * i, 5);
    }
    jb_lit(b, "]");
}

/* catalog._fixed_magnetic_subgroup_selection: *found 1 with the flat result
 * in flat and the kept records in kept; errors Python catches leave it 0. */
static int fixed_selection(int parent_sg, int64_t gid, const double *const *numeric, int free_count, int64_t full_dim, int subgroup_group,
                           const int64_t *basis, const int64_t *origin, int *found, int64_t *flat, int64_t (*kept)[5], int *kept_count)
{
    *found = 0;
    *kept_count = 0;
    double *stride = calloc((size_t)(free_count > 0 ? free_count : 1) * 48, sizeof(double));
    int64_t (*candidates)[5] = malloc(sizeof(int64_t) * 5 * 192);
    int *ops = malloc(sizeof(int) * 192), *signs = malloc(sizeof(int) * 192);
    int status = 0;
    if (!stride || !candidates || !ops || !signs) status = md_fail(MDE_MEMORY, "out of memory");
    for (int r = 0; !status && r < free_count; r++)
        for (int64_t c = 0; c < full_dim && c < 48; c++) stride[r * 48 + c] = numeric[r][c];
    int count = 0;
    if (!status) {
        int sg = (int)section_int(&T.irr_sg, gid - 1);
        count = magnetic_records_iso(mag_orderparam_group(sg), 192, &candidates[0][0]);
        if (count < 0) {
            int code = -count;
            if (code == 3 || code > 15) {
                md_err.kind = MDE_MAGNETIC;
                md_err.status = code;
                md_err.message[0] = 0;
                status = MDE_MAGNETIC;
            }
            count = -1;
        }
    }
    int64_t sub_basis[9];
    if (!status && count >= 0) {
        for (int i = 0; i < count; i++) ops[i] = mag_point_action(1, (int)candidates[i][4], &signs[i]);
        int64_t none = 0;
        int code = orderparam_to_subgroup_rows((int)gid, 0, &none, stride, (int64_t)free_count * 48, free_count, count, &candidates[0][0], ops,
                                               signs, 256, sub_basis, &kept[0][0], kept_count);
        if (code == 103 || (code && (code < 101 || code > 106))) {
            md_err.kind = MDE_ORDERPARAM;
            md_err.status = code;
            md_err.message[0] = 0;
            status = MDE_ORDERPARAM;
        } else if (!code) {
            int64_t length = id_subgroup_magnetic_identify(parent_sg, sub_basis, *kept_count, &kept[0][0], 0, 8192, flat);
            if (length == -2 || length < -7) {
                md_err.kind = MDE_ID_SUBGROUP;
                md_err.status = (int)length;
                md_err.message[0] = 0;
                status = MDE_ID_SUBGROUP;
            } else if (length > 0) {
                int same = flat[0] == subgroup_group;
                for (int i = 0; i < 9; i++) same &= flat[2 + i] == basis[i];
                for (int i = 0; i < 4; i++) same &= flat[11 + i] == origin[i];
                *found = same;
            }
        }
    }
    if (!*found) *kept_count = 0;
    free(stride);
    free(candidates);
    free(ops);
    free(signs);
    return status;
}

static int magnetic_row(jbuf *b, int64_t old_id, int64_t row_id, int canonical, int has_parent, int parent_sg_in, int has_gid, int64_t gid,
                        const sel_params *display, int has_setting, int64_t setting_id, int has_cell, const double *cell, int64_t *flat,
                        int64_t (*kept)[5])
{
    int64_t idx = row_id - 1;
    int64_t basis[9], origin[4];
    for (int i = 0; i < 9; i++) basis[i] = section_int(&T.iso_basis, idx * 9 + i);
    for (int i = 0; i < 4; i++) origin[i] = section_int(&T.iso_origin, idx * 4 + i);
    int parent_group = (int)section_int(&T.iso_parent, idx), sub_group = (int)section_int(&T.iso_subgroup, idx);
    int ordinary_parent = sel_mag_ordinary(parent_group), ordinary_sub = sel_mag_ordinary(sub_group);
    int64_t s0 = basis[0] * (basis[4] * basis[8] - basis[5] * basis[7]) - basis[1] * (basis[3] * basis[8] - basis[5] * basis[6]) +
                 basis[2] * (basis[3] * basis[7] - basis[4] * basis[6]);
    int64_t native_size = s0 < 0 ? -s0 : s0;
    sel_setting displayed;
    rat rows[3][3];
    MD_TRY(sel_display(ordinary_parent, ordinary_sub, basis, origin, has_setting, setting_id, has_cell, cell, 0, 0, sub_group, &displayed, rows));
    int64_t dim = mag_isotropy_row_dim(idx + 1), free_count = mag_isotropy_row_free(idx + 1);
    int64_t pointer = section_int(&T.iso_op_pointer, idx);
    int64_t end = idx + 1 < section_count(&T.iso_op_pointer) ? section_int(&T.iso_op_pointer, idx + 1) : section_count(&T.iso_op) + 1;
    int64_t raw_start = pointer - 1, raw_end = end - 1, total = section_count(&T.iso_op);
    if (raw_end > total) raw_end = total;
    const double *constants = (const double *)T.const_values.data;
    int fn = free_count > 0 ? (int)free_count : 0;
    double **numeric = calloc((size_t)fn + 1, sizeof(double *));
    char ***matrix = calloc((size_t)fn + 1, sizeof(char **));
    int *widths = calloc((size_t)fn + 1, sizeof(int));
    if (!numeric || !matrix || !widths) return md_fail(MDE_MEMORY, "out of memory");
    for (int p = 0; p < fn; p++) {
        int64_t a = raw_start + p * dim, e = a + dim;
        if (e > raw_end) e = raw_end;
        int w = e > a ? (int)(e - a) : 0;
        numeric[p] = calloc((size_t)w + 1, sizeof(double));
        matrix[p] = calloc((size_t)w + 1, sizeof(char *));
        widths[p] = w;
        for (int c = 0; c < w; c++) {
            numeric[p][c] = constants[section_int(&T.iso_op, a + c)];
            jbuf t = {0};
            jb_reset(&t);
            sel_fmt_display_value(&t, numeric[p][c]);
            matrix[p][c] = t.data;
        }
    }
    {
        /* the exact numeric rows, for the coupled orbit */
        int complete = fn > 0;
        for (int p = 0; p < fn && complete; p++) complete = widths[p] == dim;
        double *flat = complete ? malloc(sizeof(double) * ((size_t)fn * dim + 1)) : 0;
        int64_t *codes = complete ? malloc(sizeof(int64_t) * ((size_t)fn * dim + 1)) : 0;
        if (flat && codes) {
            for (int p = 0; p < fn; p++)
                for (int c = 0; c < dim; c++) {
                    flat[(size_t)p * dim + c] = numeric[p][c];
                    codes[(size_t)p * dim + c] = section_int(&T.iso_op, raw_start + p * dim + c);
                }
            exact_register_codes(fn, (int)dim, flat, codes);
        }
        free(flat);
        free(codes);
    }
    jbuf opd = {0};
    jb_reset(&opd);
    sel_opd_from_rows(&opd, matrix, fn, widths);
    int parent_sg = has_parent && parent_sg_in ? parent_sg_in : ordinary_parent;
    int found = 0, kept_count = 0;
    strings active = {0};
    int status = 0;
    if (has_gid) {
        status = fixed_selection(parent_sg, gid, (const double *const *)numeric, fn, dim, sub_group, basis, origin, &found, flat, kept,
                                 &kept_count);
        int64_t little = status ? 0 : sel_first_little_gid(old_id, parent_sg);
        sel_params none = {0};
        if (!status && little) status = sel_active_k_vectors(parent_sg, little, dim, matrix, fn, widths, display ? display : &none, &active);
    }
    const char *label;
    int64_t label_length;
    if (section_string(&T.iso_label, idx, &label, &label_length)) {
        label = "";
        label_length = 0;
    }
    while (label_length && (*label == ' ' || *label == '\t')) {
        label++;
        label_length--;
    }
    while (label_length && (label[label_length - 1] == ' ' || label[label_length - 1] == '\t')) label_length--;
    if (!status) {
        jb_lit(b, "{\"direction\": {\"id\": ");
        jb_int(b, row_id);
        jb_lit(b, ", \"label\": ");
        jb_string(b, label, (size_t)label_length);
        jb_lit(b, ", \"opd\": ");
        jb_string(b, opd.data, opd.used);
        jb_lit(b, ", \"dimension\": ");
        jb_int(b, dim);
        jb_lit(b, ", \"subduction\": ");
        jb_int(b, free_count);
        jb_lit(b, ", \"vectors\": ");
        jbuf vectors = {0};
        jb_reset(&vectors);
        jb_lit(&vectors, "[");
        for (int p = 0; p < fn; p++) {
            if (p) jb_lit(&vectors, ", ");
            jb_lit(&vectors, "[");
            for (int c = 0; c < widths[p]; c++) {
                if (c) jb_lit(&vectors, ", ");
                jb_cstring(&vectors, matrix[p][c]);
            }
            jb_lit(&vectors, "]");
        }
        jb_lit(&vectors, "]");
        jb_raw(b, vectors.data, vectors.used);
        jb_lit(b, "}, \"isotropy\": {\"row_id\": ");
        jb_int(b, row_id);
        jb_lit(b, ", \"magnetic\": true, \"opd_label\": ");
        jb_string(b, label, (size_t)label_length);
        jb_lit(b, ", \"direction_label\": ");
        jb_string(b, label, (size_t)label_length);
        static const char *opd_keys[4] = {", \"direction_opd\": ", ", \"stored_opd\": ", ", \"display_opd\": ", ", \"source_opd\": "};
        for (int k = 0; k < 4; k++) {
            jb_lit(b, opd_keys[k]);
            jb_string(b, opd.data, opd.used);
        }
        jb_lit(b, ", \"source_matrix\": ");
        jb_raw(b, vectors.data, vectors.used);
        jb_lit(b, ", \"source_display_rows\": ");
        jb_raw(b, vectors.data, vectors.used);
        jb_free(&vectors);
        jb_lit(b, ", \"source_numeric_rows\": [");
        for (int p = 0; p < fn; p++) {
            if (p) jb_lit(b, ", ");
            jb_lit(b, "[");
            for (int c = 0; c < widths[p]; c++) {
                if (c) jb_lit(b, ", ");
                jb_float(b, numeric[p][c]);
            }
            jb_lit(b, "]");
        }
        jb_lit(b, "], \"op_rows\": ");
        jb_int(b, dim);
        jb_lit(b, ", \"free\": ");
        jb_int(b, free_count);
        jb_lit(b, ", \"arms\": ");
        jb_int(b, displayed.size);
        jb_lit(b, ", \"s\": ");
        jb_int(b, native_size);
        jb_lit(b, ", \"i\": ");
        int has_index;
        int64_t index = sel_mag_index(parent_group, sub_group, native_size, &has_index);
        if (has_index) jb_int(b, index);
        else jb_lit(b, "null");
        jb_lit(b, ", \"k_active_vectors\": ");
        sel_put_strings(b, &active);
        jb_lit(b, ", \"k_active\": ");
        sel_put_joined(b, &active);
        jbuf subgroup = {0}, extra = {0};
        jb_reset(&subgroup);
        jb_reset(&extra);
        sel_mag_subgroup_json(&subgroup, sub_group);
        if (found) {
            if (kept_count) {
                jb_lit(&extra, ", \"source_operation_records\": ");
                sel_put_records(&extra, &kept[0][0], kept_count);
            }
        }
        sel_put_subgroup(b, ordinary_parent, ordinary_sub, basis, origin, &displayed, (const rat (*)[3])rows, displayed.size, &subgroup, &extra);
        jb_free(&subgroup);
        jb_free(&extra);
        jb_lit(b, ", \"ferroic_properties\": ");
        status = sel_mag_ferroic(b, sub_group);
        jb_lit(b, "}}");
    }
    strings_free(&active);
    jb_free(&opd);
    for (int p = 0; p < fn; p++) {
        for (int c = 0; c < widths[p]; c++) free(matrix[p][c]);
        free(matrix[p]);
        free(numeric[p]);
    }
    free(matrix);
    free(numeric);
    free(widths);
    return status;
}

int sel_magnetic_static_rows(int old_id, int has_parent, int parent_sg, int has_gid, int64_t gid, const sel_params *display,
                             const jvalue *selected_orderparam, int has_setting, int64_t setting_id, int has_cell, const double *cell, jbuf *out)
{
    if (sel_mag_tables()) return md_fail(MDE_VALUE, "magnetic Source tables are inconsistent");
    int64_t start = 0, end = 0;
    if (old_id >= 1 && old_id < section_count(&T.iso_pointer)) {
        start = section_int(&T.iso_pointer, old_id - 1);
        end = section_int(&T.iso_pointer, old_id);
    }
    int selected_any = 0;
    const char *sa = 0, *se = 0;
    if (jis(selected_orderparam, J_STRING)) {
        sa = selected_orderparam->u.s;
        se = sa + selected_orderparam->count;
        while (sa < se && sel_is_space(*sa)) sa++;
        while (se > sa && sel_is_space(se[-1])) se--;
        for (int64_t row_id = start; row_id < end && !selected_any; row_id++) {
            const char *label;
            int64_t n;
            if (section_string(&T.iso_label, row_id - 1, &label, &n)) continue;
            while (n && (*label == ' ' || *label == '\t')) {
                label++;
                n--;
            }
            while (n && (label[n - 1] == ' ' || label[n - 1] == '\t')) n--;
            selected_any = n == se - sa && !memcmp(label, sa, (size_t)n);
        }
    }
    int64_t *flat = malloc(sizeof(int64_t) * 8192);
    int64_t (*kept)[5] = malloc(sizeof(int64_t) * 5 * 193);
    if (!flat || !kept) {
        free(flat);
        free(kept);
        return md_fail(MDE_MEMORY, "out of memory");
    }
    int status = 0, first = 1;
    jb_lit(out, "[");
    for (int64_t row_id = start; row_id < end && !status; row_id++) {
        if (selected_any) {
            const char *label;
            int64_t n;
            if (section_string(&T.iso_label, row_id - 1, &label, &n)) continue;
            while (n && (*label == ' ' || *label == '\t')) {
                label++;
                n--;
            }
            while (n && (label[n - 1] == ' ' || label[n - 1] == '\t')) n--;
            if (!(n == se - sa && !memcmp(label, sa, (size_t)n))) continue;
        }
        if (!first) jb_lit(out, ", ");
        first = 0;
        status = magnetic_row(out, old_id, row_id, row_id == start, has_parent, parent_sg, has_gid, gid, display, has_setting, setting_id,
                              has_cell, cell, flat, kept);
    }
    jb_lit(out, "]");
    free(flat);
    free(kept);
    return status;
}
