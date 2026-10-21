/* The render specs of one request from the spec stage (spec/build.c), decoded into md_spec. */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "magnetic.h"
#include "md.h"
#include "tables.h"

int64_t spec_result_ints(const int64_t **ints, const double **doubles);
int64_t spec_result_occurrences(int64_t *out, int capacity);
int spec_result_occurrence_label(int index, const char **text, int64_t *length);
const char *spec_error(int *kind);
int spec_ordinary_dynamic(int sg, int child_sg, const int64_t *basis, const int64_t *origin, int selected_gid, int has_kp,
                          const int64_t *selected_kp, int has_fallback_direction);
int spec_static(int sg, int selected_old_id, int selected_row_id, int has_frame, int child_sg, const int64_t *basis,
                const int64_t *origin);
int spec_magnetic_embedding_ordinary(int sg, const int64_t *basis, int nops, const int64_t *magnetic_ops);
int spec_magnetic(int sg, const int64_t *basis, int nops, const int64_t *ops, int selected_gid, const int64_t *selected_kp,
                  const char *k_label, int64_t k_label_length, int sel_dim, int sel_count, const double *selected_direction);
int spec_coupled(int sg, const int64_t *basis, int nops, const int64_t *ops, int magnetic_embedding, int include_displacive,
                 int include_magnetic, int nslots, const int64_t *slot_words, const double *directions, const char *const *labels,
                 const int64_t *label_lengths, int npreferred, const int64_t *preferred_words, int nmag_cases,
                 const int64_t *mag_case_words);
int magnetic_records_iso(int group, int capacity, int64_t *records_out);
int orderparam_to_subgroup_rows(int gid, int param_count, const int64_t *kparam, const double *values, int64_t value_count,
                                int row_count, int candidate_count, const int64_t *candidates, const int *candidate_matrix_ops,
                                const int *candidate_signs, int max_iter, int64_t *basis_out, int64_t *kept_out, int *kept_count);

enum { K_DYN = 1, K_STATIC, K_MEMB, K_MAG, K_CPL, K_FALLBACK };

static struct {
    section_view old, full_dim;
    int ready;
} L;

static int little_tables(void)
{
    if (L.ready) return 0;
    if (section_find("irreps.little.old", &L.old) || section_find("irreps.little.dim", &L.full_dim)) return 1;
    L.ready = 1;
    return 0;
}

/* ------------------------------------------------------------- matrices */

static void matrix_free(md_matrix *m)
{
    free(m->values);
    memset(m, 0, sizeof *m);
}

static int matrix_copy(md_matrix *out, const md_matrix *in)
{
    out->rows = in->rows;
    out->cols = in->cols;
    out->values = malloc(sizeof(double) * (size_t)(in->rows * in->cols + 1));
    if (!out->values) return md_fail(MDE_MEMORY, "out of memory");
    if (in->rows * in->cols) memcpy(out->values, in->values, sizeof(double) * (size_t)(in->rows * in->cols));
    return 0;
}

static int matrix_transpose(md_matrix *out, const md_matrix *in)
{
    out->rows = in->cols;
    out->cols = in->rows;
    out->values = malloc(sizeof(double) * (size_t)(in->rows * in->cols + 1));
    if (!out->values) return md_fail(MDE_MEMORY, "out of memory");
    for (int r = 0; r < in->rows; r++)
        for (int c = 0; c < in->cols; c++) out->values[c * in->rows + r] = in->values[r * in->cols + c];
    return 0;
}

/* np.asarray(rows, dtype=float) of a JSON list: *ndim 1 for a flat or
 * empty list, 2 for a rectangular list of lists. */
int exact_refine_rows(int rows, int cols, double *values);

int md_float_rows(const jvalue *rows, md_matrix *out, int *ndim)
{
    memset(out, 0, sizeof *out);
    *ndim = 1;
    if (!jis(rows, J_ARRAY) || !rows->count) return 0;
    const jvalue *first = jat(rows, 0);
    if (!jis(first, J_ARRAY)) {
        out->rows = (int)rows->count;
        out->cols = 1;
        out->values = malloc(sizeof(double) * (size_t)rows->count);
        if (!out->values) return md_fail(MDE_MEMORY, "out of memory");
        for (int64_t i = 0; i < rows->count; i++)
            if (py_float_of(jat(rows, i), &out->values[i])) return md_fail(MDE_VALUE, "could not convert Source rows to float");
        return 0;
    }
    int64_t width = first->count;
    for (int64_t i = 0; i < rows->count; i++)
        if (!jis(jat(rows, i), J_ARRAY) || jat(rows, i)->count != width)
            return md_fail(MDE_VALUE, "setting an array element with a sequence. The requested array has an inhomogeneous shape");
    out->rows = (int)rows->count;
    out->cols = (int)width;
    out->values = malloc(sizeof(double) * (size_t)(rows->count * width + 1));
    if (!out->values) return md_fail(MDE_MEMORY, "out of memory");
    for (int64_t i = 0; i < rows->count; i++)
        for (int64_t j = 0; j < width; j++)
            if (py_float_of(jat(jat(rows, i), j), &out->values[i * width + j]))
                return md_fail(MDE_VALUE, "could not convert Source rows to float");
    exact_refine_rows(out->rows, out->cols, out->values);
    *ndim = 2;
    return 0;
}

/* ------------------------------------------------------- request helpers */

static int selected_dynamic_gid(const md_request *r, int64_t *gid)
{
    const char *keys[2] = {"gid", "dynamic_gid"};
    for (int i = 0; i < 2; i++) {
        const jvalue *raw = jget(r->irrep, keys[i]);
        int64_t value;
        if (raw && !jis(raw, J_NULL)) {
            if (py_int_of(raw, &value)) return md_fail(MDE_VALUE, "invalid literal for int()");
            if (value > 0) {
                *gid = value;
                return 1 << 8;
            }
        }
    }
    const jvalue *raw = r->iso ? jget(r->iso, "dynamic_gid") : 0;
    int64_t value;
    if (raw && !jis(raw, J_NULL)) {
        if (py_int_of(raw, &value)) return md_fail(MDE_VALUE, "invalid literal for int()");
        if (value > 0) {
            *gid = value;
            return 1 << 8;
        }
    }
    return 0;
}

int md_dynamic_gid(const md_request *r, int *has, int64_t *gid)
{
    int status = selected_dynamic_gid(r, gid);
    if (status && status != 1 << 8) return status;
    *has = status == 1 << 8;
    return 0;
}

/* the Source k parameter record (four entries) */
static int source_kparam_record(const md_request *r, int *has, int64_t *kp)
{
    *has = 0;
    const jvalue *raw = r->iso ? jget(r->iso, "source_kparam") : 0;
    if (!raw || jis(raw, J_NULL)) return 0;
    if (!jis(raw, J_ARRAY) || raw->count != 4) return md_fail(MDE_TYPE, "source_kparam is not a four-entry record");
    for (int i = 0; i < 4; i++)
        if (py_int_of(jat(raw, i), &kp[i])) return md_fail(MDE_VALUE, "invalid literal for int()");
    *has = 1;
    return 0;
}

static void isotropy_row_id(const md_request *r, int *state, int64_t *row_id)
{
    *state = MD_NONE;
    if (r->iso && !py_int_of(jget(r->iso, "row_id"), row_id)) *state = MD_SET;
}

char *md_k_label_from_irrep_label(const char *label)
{
    const char *s = label, *e = label + strlen(label);
    while (s < e && isspace((unsigned char)*s)) s++;
    while (e > s && isspace((unsigned char)e[-1])) e--;
    if (e - s > 1 && s[0] == 'm' && isupper((unsigned char)s[1])) s++;
    const char *z = s;
    while (z < e && isalpha((unsigned char)*z)) z++;
    if (z == s) z = e;
    return (char *)md_strdup(s, z - s);
}

/* the integer Source basis of the selected OPD row */
int md_integer_source_basis(const md_request *r, int64_t *out)
{
    const jvalue *basis = r->iso ? jget(r->iso, "source_basis_values") : 0;
    if (!jis(basis, J_ARRAY) || basis->count != 9) return 0;
    double values[9];
    for (int i = 0; i < 9; i++)
        if (py_float_of(jat(basis, i), &values[i])) return 0;
    for (int i = 0; i < 9; i++) {
        rat raw;
        if (py_fraction_of_float(values[i], &raw)) return 0;
        rat snapped = py_limit_denominator(raw, 1000000);
        if (snapped.d != 1) return 0;
        out[i] = (int64_t)snapped.n;
    }
    return 1;
}

static int source_split_origin(const md_request *r, int64_t *out)
{
    const jvalue *origin = r->iso ? jget(r->iso, "source_origin_values") : 0;
    if (!jis(origin, J_ARRAY) || origin->count != 4) return 0;
    for (int i = 0; i < 4; i++)
        if (py_int_of(jat(origin, i), &out[i])) return 0;
    return 1;
}

/* iso.source_numeric_rows or iso.source_matrix or [] */
static const jvalue *selected_source_rows(const md_request *r)
{
    if (!r->iso) return 0;
    const jvalue *rows = jget(r->iso, "source_numeric_rows");
    if (jtruthy(rows)) return rows;
    rows = jget(r->iso, "source_matrix");
    return jtruthy(rows) ? rows : 0;
}

/* ---------------------------------------------------------------- specs */

void md_spec_free(md_spec *s)
{
    free(s->label);
    free(s->display_label);
    free(s->k_label);
    free(s->heading_label);
    free(s->irrep_name);
    for (int i = 0; i < s->nocc; i++) free(s->occ[i].k_label);
    free(s->occ);
    matrix_free(&s->dir);
    matrix_free(&s->rows);
    matrix_free(&s->recomputed);
    memset(s, 0, sizeof *s);
}

static md_spec *push_spec(md_request *r)
{
    if (r->nspecs == r->spec_cap) {
        int cap = r->spec_cap ? 2 * r->spec_cap : 16;
        md_spec *grown = realloc(r->specs, sizeof(md_spec) * (size_t)cap);
        if (!grown) {
            md_fail(MDE_MEMORY, "out of memory");
            return 0;
        }
        r->specs = grown;
        r->spec_cap = cap;
    }
    md_spec *s = r->specs + r->nspecs++;
    memset(s, 0, sizeof *s);
    return s;
}

typedef struct {
    const int64_t *ints;
    const double *doubles;
    int64_t at, used;
} stream;

static int64_t take(stream *s)
{
    return s->at < s->used ? s->ints[s->at++] : 0;
}

static char *take_string(stream *s)
{
    int64_t size = take(s);
    char *out = malloc((size_t)size + 1);
    if (!out) return 0;
    for (int64_t i = 0; i < size; i += 8) {
        int64_t word = take(s);
        memcpy(out + i, &word, (size_t)(size - i < 8 ? size - i : 8));
    }
    out[size] = 0;
    return out;
}

static int take_rats(stream *s, rat *out, int capacity)
{
    int64_t n = take(s);
    for (int64_t i = 0; i < n; i++) {
        int64_t num = take(s), den = take(s);
        if (i < capacity) out[i] = rat_make(num, den);
    }
    return (int)n;
}

static int take_matrix(stream *s, md_matrix *out)
{
    int64_t rows = take(s), cols = take(s), offset = take(s);
    out->rows = (int)rows;
    out->cols = (int)cols;
    out->values = malloc(sizeof(double) * (size_t)(rows * cols + 1));
    if (!out->values) return md_fail(MDE_MEMORY, "out of memory");
    if (rows * cols) memcpy(out->values, s->doubles + offset, sizeof(double) * (size_t)(rows * cols));
    return 0;
}

enum {
    H_KIND, H_MODE_KIND, H_OLD_ID, H_GID, H_LITTLE_TYPE, H_MAGNETIC, H_PRIMARY, H_HAS_KP, H_HAS_CARRIER, H_HAS_PML, H_HAS_REQ,
    H_REQ_NONE, H_HAS_OCC, H_HAS_DIR, H_HAS_FREE, H_FREE_COUNT, H_OBS, H_HAS_ANCHOR, H_ANCHOR_GID, H_HAS_CAND, H_HAS_HEADING,
    H_HEADING_GID, H_STOKES, H_CASE_NONE, H_DCASE_NONE, H_SLOT_ORDER, H_REPLACED, H_ROW_ID,
    H_USE_SELECTED_ROWS, H_STATIC_UPDATED, H_COUNT
};

typedef struct {
    const md_spec *fallback;
    int64_t selected_free;
    int has_selected_rows;
    const md_matrix *selected_rows; /* floats of the selected Source rows */
    int selected_dir_state;
    const md_matrix *selected_dir;
} spec_context;

static int spec_error_of(int status)
{
    int kind = 0;
    const char *message = spec_error(&kind);
    return md_fail(kind ? kind : (status ? status : MDE_RUNTIME), "%s", message);
}

static int copy_spec(md_spec *out, const md_spec *in)
{
    *out = *in;
    out->label = in->label ? (char *)md_strdup(in->label, (int64_t)strlen(in->label)) : 0;
    out->display_label = in->display_label ? (char *)md_strdup(in->display_label, (int64_t)strlen(in->display_label)) : 0;
    out->k_label = in->k_label ? (char *)md_strdup(in->k_label, (int64_t)strlen(in->k_label)) : 0;
    out->heading_label = in->heading_label ? (char *)md_strdup(in->heading_label, (int64_t)strlen(in->heading_label)) : 0;
    out->irrep_name = in->irrep_name ? (char *)md_strdup(in->irrep_name, (int64_t)strlen(in->irrep_name)) : 0;
    out->occ = 0;
    if (in->nocc) {
        out->occ = malloc(sizeof(md_occurrence) * (size_t)in->nocc);
        if (!out->occ) return md_fail(MDE_MEMORY, "out of memory");
        for (int i = 0; i < in->nocc; i++) {
            out->occ[i] = in->occ[i];
            out->occ[i].k_label = (char *)md_strdup(in->occ[i].k_label, (int64_t)strlen(in->occ[i].k_label));
        }
    }
    memset(&out->dir, 0, sizeof out->dir);
    memset(&out->rows, 0, sizeof out->rows);
    memset(&out->recomputed, 0, sizeof out->recomputed);
    if (in->dir.values) MD_TRY(matrix_copy(&out->dir, &in->dir));
    if (in->rows.values) MD_TRY(matrix_copy(&out->rows, &in->rows));
    if (in->recomputed.values) MD_TRY(matrix_copy(&out->recomputed, &in->recomputed));
    return 0;
}

/* _specs: one md_spec per record of the last spec-stage result. */
static int decode_specs(md_request *r, const spec_context *ctx)
{
    stream s = {0};
    s.used = spec_result_ints(&s.ints, &s.doubles);
    int64_t occurrence_count = spec_result_occurrences(0, 0);
    int64_t *pool = malloc(sizeof(int64_t) * 18 * (size_t)(occurrence_count + 1));
    if (!pool) return md_fail(MDE_MEMORY, "out of memory");
    spec_result_occurrences(pool, (int)occurrence_count);
    int64_t count = take(&s);
    int status = 0;
    for (int64_t index = 0; index < count && !status; index++) {
        int64_t width = take(&s), h[H_COUNT] = {0};
        for (int64_t k = 0; k < width; k++) {
            int64_t value = take(&s);
            if (k < H_COUNT) h[k] = value;
        }
        char *label = take_string(&s), *k_label = take_string(&s), *heading = take_string(&s);
        int64_t kp[4], carrier[4], anchor_kp[4], cand_kp[4], heading_kp[4];
        for (int k = 0; k < 4; k++) kp[k] = take(&s);
        for (int k = 0; k < 4; k++) carrier[k] = take(&s);
        for (int k = 0; k < 4; k++) anchor_kp[k] = take(&s);
        for (int k = 0; k < 4; k++) cand_kp[k] = take(&s);
        for (int k = 0; k < 4; k++) heading_kp[k] = take(&s);
        rat pml[3], req[3], cases[3], dcases[3];
        int pml_n = take_rats(&s, pml, 3), req_n = take_rats(&s, req, 3), case_n = take_rats(&s, cases, 3),
            dcase_n = take_rats(&s, dcases, 3);
        int64_t nocc = take(&s);
        int64_t *occ = malloc(sizeof(int64_t) * (size_t)(nocc + 1));
        for (int64_t k = 0; k < nocc; k++) {
            int64_t value = take(&s);
            if (occ) occ[k] = value;
        }
        md_matrix dir = {0}, recomputed = {0};
        if (!label || !k_label || !heading || !occ || take_matrix(&s, &dir) || take_matrix(&s, &recomputed)) {
            free(label);
            free(k_label);
            free(heading);
            free(occ);
            matrix_free(&dir);
            matrix_free(&recomputed);
            status = md_fail(MDE_MEMORY, "out of memory");
            break;
        }
        md_spec *spec = push_spec(r);
        if (!spec) {
            status = MDE_MEMORY;
            break;
        }
        int kind = (int)h[H_KIND];
        if (kind == K_FALLBACK) {
            status = copy_spec(spec, ctx->fallback);
            free(label);
            free(k_label);
            free(heading);
            free(occ);
            matrix_free(&dir);
            matrix_free(&recomputed);
            continue;
        }
        spec->kind = kind;
        spec->old_id = h[H_OLD_ID];
        spec->gid_state = MD_SET;
        spec->gid = h[H_GID];
        spec->label = label;
        spec->k_label = k_label;
        spec->heading_label = heading;
        spec->primary = h[H_PRIMARY] != 0;
        spec->row_id_state = MD_NONE;
        spec->occ_state = MD_SET;
        spec->nocc = (int)nocc;
        spec->occ = calloc((size_t)nocc + 1, sizeof(md_occurrence));
        for (int64_t k = 0; k < nocc && spec->occ; k++) {
            const int64_t *w = pool + 18 * occ[k];
            md_occurrence *o = spec->occ + k;
            o->mesh = w[0];
            memcpy(o->kp, w + 1, sizeof o->kp);
            for (int axis = 0; axis < 3; axis++) o->pml[axis] = rat_make(w[5 + 2 * axis], w[6 + 2 * axis]);
            memcpy(o->bucket, w + 11, sizeof o->bucket);
            const char *text;
            int64_t length;
            if (spec_result_occurrence_label((int)occ[k], &text, &length)) length = 0;
            o->k_label = (char *)md_strdup(length ? text : "", length);
        }
        free(occ);
        int has_dir_rows = 1;
        spec->kp_state = MD_SET;
        memcpy(spec->kp, kp, sizeof kp);
        spec->carrier_state = h[H_HAS_CARRIER] ? MD_SET : MD_NONE;
        memcpy(spec->carrier, carrier, sizeof carrier);
        spec->pml_state = MD_SET;
        spec->pml_n = pml_n;
        memcpy(spec->pml, pml, sizeof pml);
        spec->req_state = h[H_REQ_NONE] ? MD_NONE : MD_SET;
        spec->req_n = req_n;
        memcpy(spec->req, req, sizeof req);
        switch (kind) {
        case K_DYN:
            spec->magnetic_state = 1;
            if (h[H_HAS_FREE]) {
                spec->free_state = MD_SET;
                spec->free_count = ctx->selected_free;
            }
            if (h[H_USE_SELECTED_ROWS]) {
                has_dir_rows = 0;
                spec->rows_state = MD_SET;
                if (ctx->has_selected_rows) MD_TRY(matrix_copy(&spec->rows, ctx->selected_rows));
                spec->dir_state = ctx->selected_dir_state;
                if (ctx->selected_dir_state == MD_SET) MD_TRY(matrix_copy(&spec->dir, ctx->selected_dir));
            }
            break;
        case K_STATIC:
            spec->carrier_state = MD_ABSENT;
            spec->req_state = MD_ABSENT;
            spec->occ_state = MD_ABSENT;
            spec->kp_state = MD_ABSENT;
            spec->pml_state = MD_ABSENT;
            spec->gid_state = h[H_GID] ? MD_SET : MD_NONE;
            spec->row_id_state = h[H_ROW_ID] > 0 ? MD_SET : MD_NONE;
            spec->row_id = h[H_ROW_ID];
            has_dir_rows = 0;
            if (h[H_USE_SELECTED_ROWS]) {
                spec->rows_state = MD_SET;
                if (ctx->has_selected_rows) MD_TRY(matrix_copy(&spec->rows, ctx->selected_rows));
                spec->dir_state = ctx->selected_dir_state;
                if (ctx->selected_dir_state == MD_SET) MD_TRY(matrix_copy(&spec->dir, ctx->selected_dir));
            }
            spec->free_state = MD_SET;
            spec->free_count = h[H_FREE_COUNT];
            if (h[H_STATIC_UPDATED]) {
                spec->gid_state = MD_SET;
                spec->kp_state = MD_SET;
                spec->pml_state = MD_SET;
                matrix_free(&spec->dir);
                matrix_free(&spec->rows);
                has_dir_rows = 1;
            }
            break;
        case K_MEMB:
            spec->carrier_state = MD_ABSENT;
            spec->primary = 0;
            break;
        case K_MAG:
            spec->display_label = md_format("m%s", label);
            spec->magnetic_state = 2;
            spec->free_state = MD_SET;
            spec->free_count = h[H_FREE_COUNT];
            if (h[H_STOKES]) {
                spec->stokes = 1;
                spec->recomputed = recomputed;
                memset(&recomputed, 0, sizeof recomputed);
            }
            break;
        default:
            spec->little_type_state = MD_SET;
            spec->little_type = h[H_LITTLE_TYPE];
            spec->display_label = md_format("%s%s", h[H_MAGNETIC] ? "m" : "", label);
            spec->magnetic_state = h[H_MAGNETIC] ? 2 : 1;
            spec->free_state = MD_SET;
            spec->free_count = h[H_FREE_COUNT];
            spec->coupled = 1;
            spec->case_state = h[H_CASE_NONE] ? MD_NONE : MD_SET;
            spec->case_n = case_n;
            memcpy(spec->case_k, cases, sizeof cases);
            spec->dcase_state = h[H_DCASE_NONE] ? MD_NONE : MD_SET;
            spec->dcase_n = dcase_n;
            memcpy(spec->dcase_k, dcases, sizeof dcases);
            spec->slot_order_state = h[H_SLOT_ORDER] < 0 ? MD_NONE : MD_SET;
            spec->slot_order = h[H_SLOT_ORDER];
            spec->replaced_state = MD_SET;
            spec->replaced = h[H_REPLACED] != 0;
            spec->req_state = MD_ABSENT;
            break;
        }
        if (has_dir_rows) {
            matrix_free(&spec->dir);
            matrix_free(&spec->rows);
            spec->dir_state = MD_SET;
            spec->dir = dir;
            memset(&dir, 0, sizeof dir);
            spec->rows_state = MD_SET;
            MD_TRY(matrix_transpose(&spec->rows, &spec->dir));
        }
        matrix_free(&dir);
        matrix_free(&recomputed);
        if (h[H_OBS]) spec->obs = 1;
        if (h[H_HAS_ANCHOR]) {
            spec->has_anchor = 1;
            spec->anchor_gid = h[H_ANCHOR_GID];
            memcpy(spec->anchor_kp, anchor_kp, sizeof anchor_kp);
        }
        if (h[H_HAS_CAND]) {
            spec->has_cand = 1;
            memcpy(spec->cand_kp, cand_kp, sizeof cand_kp);
        }
        if (h[H_HAS_HEADING]) {
            spec->has_heading = 1;
            spec->heading_gid = h[H_HEADING_GID];
            memcpy(spec->heading_kp, heading_kp, sizeof heading_kp);
        }
        if ((spec->display_label == 0 && (kind == K_MAG || kind == K_CPL)) || !spec->occ) status = md_fail(MDE_MEMORY, "out of memory");
    }
    free(pool);
    return status;
}

static int set_render(md_spec *s, int magnetic, const char *label)
{
    s->magnetic_kind = magnetic;
    s->vector_setting = magnetic ? 2 : 1;
    free(s->irrep_name);
    s->irrep_name = (char *)md_strdup(label, (int64_t)strlen(label));
    return s->irrep_name ? 0 : md_fail(MDE_MEMORY, "out of memory");
}

/* ------------------------------------------------- magnetic frame */

typedef struct {
    int found;
    int64_t gid, kp[4];
    md_matrix rows;
    int64_t basis[9];
    int64_t (*ops)[5];
    int nops;
} magnetic_frame;

static struct {
    section_view irr_sg;
    int ready;
} G;

static int build_magnetic_frame(md_request *r, int require_gid_first, magnetic_frame *f)
{
    memset(f, 0, sizeof *f);
    int has_gid = 0, has_kp = 0;
    MD_TRY(md_dynamic_gid(r, &has_gid, &f->gid));
    MD_TRY(source_kparam_record(r, &has_kp, f->kp));
    const jvalue *raw_rows = selected_source_rows(r);
    if (little_tables()) return md_fail(MDE_VALUE, "mode Source tables are inconsistent");
    int zero_kp = 0;
    if (require_gid_first) {
        if (!has_gid) return 0;
        zero_kp = !has_kp && section_int(&L.old, f->gid - 1) > 0;
        if (zero_kp) has_kp = 1;
        if (!has_kp || !raw_rows) return 0;
    } else {
        zero_kp = has_gid && !has_kp && section_int(&L.old, f->gid - 1) > 0;
        if (zero_kp) has_kp = 1;
        if (!has_gid || !has_kp || !raw_rows) return 0;
    }
    if (zero_kp) {
        f->kp[0] = f->kp[1] = f->kp[2] = 0;
        f->kp[3] = 1;
    }
    int full_dim = (int)section_int(&L.full_dim, f->gid - 1);
    int64_t row_count = raw_rows->count;
    double *orderparam = calloc((size_t)(48 * row_count + 1), sizeof(double));
    if (!orderparam) return md_fail(MDE_MEMORY, "out of memory");
    for (int64_t i = 0; i < row_count; i++) {
        const jvalue *row = jat(raw_rows, i);
        for (int64_t j = 0; jis(row, J_ARRAY) && j < row->count; j++) {
            double value;
            if (py_float_of(jat(row, j), &value)) {
                free(orderparam);
                return md_fail(MDE_VALUE, "could not convert Source rows to float");
            }
            if (j < full_dim && j < 48) orderparam[48 * i + j] = value;
        }
    }
    {
        /* six-decimal dynamic rows: their exact values */
        double *compact = full_dim <= 48 ? malloc(sizeof(double) * ((size_t)row_count * full_dim + 1)) : 0;
        if (compact) {
            for (int64_t i = 0; i < row_count; i++)
                for (int j = 0; j < full_dim; j++) compact[i * full_dim + j] = orderparam[48 * i + j];
            if (exact_refine_rows((int)row_count, full_dim, compact))
                for (int64_t i = 0; i < row_count; i++)
                    for (int j = 0; j < full_dim; j++) orderparam[48 * i + j] = compact[i * full_dim + j];
            free(compact);
        }
    }
    if (!G.ready) {
        if (section_find("irreps.little.sg", &G.irr_sg)) {
            free(orderparam);
            return md_fail(MDE_VALUE, "mode Source tables are inconsistent");
        }
        G.ready = 1;
    }
    int sg = (int)section_int(&G.irr_sg, f->gid - 1);
    int64_t (*candidates)[5] = malloc(sizeof(int64_t) * 5 * 192);
    int *matrix_ops = malloc(sizeof(int) * 192), *signs = malloc(sizeof(int) * 192);
    f->ops = malloc(sizeof(int64_t) * 5 * 193);
    if (!candidates || !matrix_ops || !signs || !f->ops) {
        free(orderparam);
        free(candidates);
        free(matrix_ops);
        free(signs);
        return md_fail(MDE_MEMORY, "out of memory");
    }
    int status = 0;
    if (mag_tables()) status = md_fail(MDE_VALUE, "magnetic Source tables are inconsistent");
    int count = status ? 0 : magnetic_records_iso(mag_orderparam_group(sg), 192, &candidates[0][0]);
    if (!status && count < 0) {
        md_err.kind = MDE_MAGNETIC;
        md_err.status = -count;
        md_err.message[0] = 0;
        status = MDE_MAGNETIC;
    }
    for (int i = 0; i < count && !status; i++) matrix_ops[i] = mag_point_action(1, (int)candidates[i][4], &signs[i]);
    if (!status) {
        int kept = 0;
        int code = orderparam_to_subgroup_rows((int)f->gid, 4, f->kp, orderparam, 48 * row_count, (int)row_count, count,
                                               &candidates[0][0], matrix_ops, signs, 256, f->basis, &f->ops[0][0], &kept);
        if (code) {
            md_err.kind = MDE_ORDERPARAM;
            md_err.status = code;
            md_err.message[0] = 0;
            status = MDE_ORDERPARAM;
        }
        f->nops = kept;
    }
    free(orderparam);
    free(candidates);
    free(matrix_ops);
    free(signs);
    if (status) return status;
    int ndim;
    MD_TRY(md_float_rows(raw_rows, &f->rows, &ndim));
    f->found = 1;
    return 0;
}

static void frame_free(magnetic_frame *f)
{
    free(f->ops);
    matrix_free(&f->rows);
}

/* _ordinary_subduced_mode_specs_from_magnetic_embedding */
static int magnetic_embedding_specs(md_request *r)
{
    magnetic_frame f;
    int status = build_magnetic_frame(r, 1, &f);
    if (!status && f.found) {
        int code = spec_magnetic_embedding_ordinary(r->sg, f.basis, f.nops, &f.ops[0][0]);
        spec_context ctx = {0};
        status = code ? spec_error_of(code) : decode_specs(r, &ctx);
    }
    frame_free(&f);
    return status;
}

/* _magnetic_subduced_mode_specs */
static int magnetic_specs(md_request *r)
{
    magnetic_frame f;
    int status = build_magnetic_frame(r, 0, &f);
    if (!status && f.found) {
        char *label;
        const jvalue *k = jget(r->irrep, "k_label");
        if (jtruthy(k)) {
            label = md_py_str(k);
        } else {
            const jvalue *symbol = jget(r->irrep, "ordinary_symbol");
            char *text = jtruthy(symbol) ? md_py_str(symbol) : (char *)md_strdup("", 0);
            label = text ? md_k_label_from_irrep_label(text) : 0;
            free(text);
        }
        int ndim;
        md_matrix selected = {0}, transposed = {0};
        const jvalue *raw_rows = selected_source_rows(r);
        if (!label) status = md_fail(MDE_MEMORY, "out of memory");
        if (!status) status = md_float_rows(raw_rows, &selected, &ndim);
        if (!status) {
            int finite = 1;
            for (int i = 0; i < selected.rows * selected.cols; i++) finite &= isfinite(selected.values[i]) != 0;
            if (ndim != 2 || selected.rows < 1 || selected.cols < 1 || !finite)
                status = md_fail(MDE_VALUE, "selected Stokes direction is malformed");
        }
        if (!status) status = matrix_transpose(&transposed, &selected);
        if (!status) {
            int code = spec_magnetic(r->sg, f.basis, f.nops, &f.ops[0][0], (int)f.gid, f.kp, label, (int64_t)strlen(label),
                                     transposed.rows, transposed.cols, transposed.values);
            spec_context ctx = {0};
            status = code ? spec_error_of(code) : decode_specs(r, &ctx);
        }
        free(label);
        matrix_free(&selected);
        matrix_free(&transposed);
    }
    frame_free(&f);
    return status;
}

/* ------------------------------------------------------------ coupled */

static int k_params_of(const jvalue *params, rat *out, int *count)
{
    *count = 0;
    if (!jtruthy(params) || !jis(params, J_OBJECT)) return 0;
    static const char *keys[6] = {"a", "b", "g", "alpha", "beta", "gamma"};
    for (int i = 0; i < 6; i++) {
        const jvalue *value = jget(params, keys[i]);
        if (!value || jis(value, J_NULL)) continue;
        char *text = md_py_str(value);
        if (!text) return md_fail(MDE_MEMORY, "out of memory");
        const char *a = text, *e = text + strlen(text);
        while (a < e && isspace((unsigned char)*a)) a++;
        while (e > a && isspace((unsigned char)e[-1])) e--;
        if (a == e) {
            free(text);
            continue;
        }
        rat value_rat;
        int status = py_fraction_of_text(a, e - a, &value_rat);
        if (status) {
            jbuf b = {0};
            jb_reset(&b);
            py_repr_text(&b, a, e - a);
            if (status == 2) md_fail(MDE_ARITH, "k parameter out of range");
            else md_fail(MDE_VALUE, "Invalid literal for Fraction: %s", b.data);
            jb_free(&b);
            free(text);
            return md_err.kind;
        }
        free(text);
        if (*count < 3) out[*count] = value_rat;
        (*count)++;
    }
    if (*count > 3) return md_fail(MDE_VALUE, "too many k parameters");
    return 0;
}

static int strict_kp(const jvalue *raw, int64_t *out)
{
    if (!jis(raw, J_ARRAY) || raw->count != 4) return 0;
    for (int i = 0; i < 4; i++) {
        const jvalue *v = jat(raw, i);
        if (!jis(v, J_INT)) return 0;
        out[i] = v->u.i;
    }
    return out[3] != 0;
}

typedef struct {
    char **items;
    int count, cap;
} label_list;

static int label_index(label_list *l, const char *text)
{
    for (int i = 0; i < l->count; i++)
        if (!strcmp(l->items[i], text)) return i;
    if (l->count == l->cap) {
        int cap = l->cap ? 2 * l->cap : 8;
        char **grown = realloc(l->items, sizeof(char *) * (size_t)cap);
        if (!grown) return -1;
        l->items = grown;
        l->cap = cap;
    }
    l->items[l->count] = (char *)md_strdup(text, (int64_t)strlen(text));
    return l->items[l->count] ? l->count++ : -1;
}

typedef struct {
    int64_t *data;
    size_t used, cap;
} words;

static int put_word(words *w, int64_t value)
{
    if (w->used == w->cap) {
        size_t cap = w->cap ? 2 * w->cap : 64;
        int64_t *grown = realloc(w->data, sizeof(int64_t) * cap);
        if (!grown) return md_fail(MDE_MEMORY, "out of memory");
        w->data = grown;
        w->cap = cap;
    }
    w->data[w->used++] = value;
    return 0;
}

static int put_fraction_words(words *w, const rat *values, int count)
{
    for (int i = 0; i < 3; i++) {
        rat v = i < count ? values[i] : rat_int(0);
        if (!rat_fits64(v)) return md_fail(MDE_ARITH, "k parameter out of range");
        MD_TRY(put_word(w, (int64_t)v.n));
        MD_TRY(put_word(w, i < count ? (int64_t)v.d : 1));
    }
    return 0;
}

int md_coupled_specs(md_request *r, int include_displacive, int include_magnetic)
{
    int64_t basis[9];
    if (!md_integer_source_basis(r, basis)) return 0;
    const jvalue *raw_ops = jget(r->iso, "source_operation_records");
    int64_t nraw = jtruthy(raw_ops) && jis(raw_ops, J_ARRAY) ? raw_ops->count : 0;
    int64_t *ops = malloc(sizeof(int64_t) * 5 * (size_t)(nraw + 1));
    if (!ops) return md_fail(MDE_MEMORY, "out of memory");
    int nops = 0;
    for (int64_t i = 0; i < nraw; i++) {
        const jvalue *record = jat(raw_ops, i);
        if (!jis(record, J_ARRAY) || record->count != 5) continue;
        for (int k = 0; k < 5; k++)
            if (py_int_of(jat(record, k), &ops[5 * nops + k])) {
                free(ops);
                return md_fail(MDE_VALUE, "invalid literal for int()");
            }
        nops++;
    }
    if (!nops) {
        free(ops);
        return 0;
    }
    const jvalue *slot_rows = jget(r->iso, "_slot_source_numeric_rows");
    int64_t nslots = jis(r->slots, J_ARRAY) ? r->slots->count : 0;
    int64_t nrows = jis(slot_rows, J_ARRAY) ? slot_rows->count : 0;
    if (nrows != nslots) {
        free(ops);
        return md_fail(MDE_VALUE, "coupled OPD lost its per-slot Source rows");
    }
    const jvalue *subgroup = jget(r->iso, "subgroup");
    int magnetic_embedding = jis(subgroup, J_OBJECT) && jtruthy(jget(subgroup, "ordinary_number"));
    label_list labels = {0};
    words slot_words = {0}, preferred = {0}, cases = {0};
    double *directions = 0;
    size_t ndirections = 0;
    int nslot_used = 0, npreferred = 0, ncases = 0, status = 0;
    for (int64_t index = 0; index < nslots && !status; index++) {
        const jvalue *slot = jat(r->slots, index);
        const jvalue *irrep = jget(slot, "irrep");
        int magnetic = jis(irrep, J_OBJECT) && jtruthy(jget(irrep, "magnetic"));
        int64_t gid = 0, kp[4] = {0, 0, 0, 1};
        rat case_k[3];
        int case_n = 0;
        if (magnetic) {
            const jvalue *raw_gid = jget(irrep, "gid");
            const jvalue *raw_kp = jget(slot, "source_kparam");
            int64_t mkp[4] = {0, 0, 0, 1};
            int ok_gid = jis(raw_gid, J_INT) && raw_gid->u.i > 0;
            int ok_kp = (!raw_kp || jis(raw_kp, J_NULL)) ? 1 : strict_kp(raw_kp, mkp);
            if (!ok_gid || !ok_kp) continue;
            gid = raw_gid->u.i;
            memcpy(kp, mkp, sizeof kp);
            status = k_params_of(jget(slot, "k_params"), case_k, &case_n);
            if (status) break;
            status = put_word(&cases, gid);
            for (int k = 0; k < 4 && !status; k++) status = put_word(&cases, kp[k]);
            if (!status) status = put_word(&cases, case_n);
            if (!status) status = put_fraction_words(&cases, case_k, case_n);
            ncases++;
        } else {
            const jvalue *raw_gid = jis(irrep, J_OBJECT) ? jget(irrep, "gid") : 0;
            if (jtruthy(raw_gid) && py_int_of(raw_gid, &gid)) status = md_fail(MDE_VALUE, "invalid literal for int()");
            const jvalue *raw_kp = jget(slot, "source_kparam");
            if (!status && jtruthy(raw_kp)) {
                for (int k = 0; k < 4 && !status; k++)
                    if (py_int_of(jat(raw_kp, k), &kp[k])) status = md_fail(MDE_VALUE, "invalid literal for int()");
            }
        }
        if (status) break;
        const jvalue *rows = jat(slot_rows, index);
        md_matrix m = {0};
        int ndim;
        status = md_float_rows(rows, &m, &ndim);
        if (status) break;
        int drows = 0, dcols = 0;
        if (ndim == 2) {
            drows = m.cols;
            dcols = m.rows;
        }
        double *grown = realloc(directions, sizeof(double) * (ndirections + (size_t)(drows * dcols) + 1));
        if (!grown) {
            matrix_free(&m);
            status = md_fail(MDE_MEMORY, "out of memory");
            break;
        }
        directions = grown;
        size_t offset = ndirections;
        for (int i = 0; i < drows; i++)
            for (int j = 0; j < dcols; j++) directions[ndirections++] = m.values[j * m.cols + i];
        matrix_free(&m);
        const jvalue *kpoint = jget(slot, "kpoint");
        const jvalue *klabel = jis(kpoint, J_OBJECT) ? jget(kpoint, "label") : 0;
        char *k_text = jtruthy(klabel) ? md_py_str(klabel) : (char *)md_strdup("", 0);
        int li = k_text ? label_index(&labels, k_text) : -1;
        free(k_text);
        if (li < 0) {
            status = md_fail(MDE_MEMORY, "out of memory");
            break;
        }
        rat slot_case[3], slot_dcase[3];
        int slot_case_n, slot_dcase_n;
        status = k_params_of(jget(slot, "k_params"), slot_case, &slot_case_n);
        if (status) break;
        const jvalue *display = jget(slot, "display_k_params");
        status = k_params_of(jtruthy(display) ? display : jget(slot, "k_params"), slot_dcase, &slot_dcase_n);
        if (status) break;
        int64_t header[] = {index, magnetic, gid, kp[0], kp[1], kp[2], kp[3], li, slot_case_n};
        for (size_t k = 0; k < sizeof header / sizeof header[0] && !status; k++) status = put_word(&slot_words, header[k]);
        if (!status) status = put_fraction_words(&slot_words, slot_case, slot_case_n);
        if (!status) status = put_word(&slot_words, slot_dcase_n);
        if (!status) status = put_fraction_words(&slot_words, slot_dcase, slot_dcase_n);
        if (!status) status = put_word(&slot_words, drows);
        if (!status) status = put_word(&slot_words, dcols);
        if (!status) status = put_word(&slot_words, (int64_t)offset);
        nslot_used++;
    }
    /* preferred: per k label in first-appearance order, every ordinary slot's kparam */
    label_list order = {0};
    for (int64_t index = 0; index < nslots && !status; index++) {
        const jvalue *slot = jat(r->slots, index), *irrep = jget(slot, "irrep");
        if (jis(irrep, J_OBJECT) && jtruthy(jget(irrep, "magnetic"))) continue;
        const jvalue *kpoint = jget(slot, "kpoint"), *klabel = jis(kpoint, J_OBJECT) ? jget(kpoint, "label") : 0;
        char *k_text = jtruthy(klabel) ? md_py_str(klabel) : (char *)md_strdup("", 0);
        if (!k_text || label_index(&order, k_text) < 0) status = md_fail(MDE_MEMORY, "out of memory");
        free(k_text);
    }
    for (int g = 0; g < order.count && !status; g++)
        for (int64_t index = 0; index < nslots && !status; index++) {
            const jvalue *slot = jat(r->slots, index), *irrep = jget(slot, "irrep");
            if (jis(irrep, J_OBJECT) && jtruthy(jget(irrep, "magnetic"))) continue;
            const jvalue *kpoint = jget(slot, "kpoint"), *klabel = jis(kpoint, J_OBJECT) ? jget(kpoint, "label") : 0;
            char *k_text = jtruthy(klabel) ? md_py_str(klabel) : (char *)md_strdup("", 0);
            int same = k_text && !strcmp(k_text, order.items[g]);
            free(k_text);
            if (!same) continue;
            int64_t kp[4] = {0, 0, 0, 1};
            const jvalue *raw_kp = jget(slot, "source_kparam");
            if (jtruthy(raw_kp))
                for (int k = 0; k < 4 && !status; k++)
                    if (py_int_of(jat(raw_kp, k), &kp[k])) status = md_fail(MDE_VALUE, "invalid literal for int()");
            int li = label_index(&labels, order.items[g]);
            if (!status) status = li < 0 ? md_fail(MDE_MEMORY, "out of memory") : put_word(&preferred, li);
            for (int k = 0; k < 4 && !status; k++) status = put_word(&preferred, kp[k]);
            npreferred++;
        }
    if (!status) {
        int64_t *lengths = malloc(sizeof(int64_t) * (size_t)(labels.count + 1));
        if (!lengths) status = md_fail(MDE_MEMORY, "out of memory");
        for (int i = 0; i < labels.count && !status; i++) lengths[i] = (int64_t)strlen(labels.items[i]);
        if (!status) {
            int64_t zero = 0;
            double no_direction = 0.0;
            int code = spec_coupled(r->sg, basis, nops, ops, magnetic_embedding, include_displacive, include_magnetic, nslot_used,
                                    slot_words.used ? slot_words.data : &zero, directions ? directions : &no_direction,
                                    (const char *const *)labels.items, lengths, npreferred, preferred.used ? preferred.data : &zero,
                                    ncases, cases.used ? cases.data : &zero);
            spec_context ctx = {0};
            int first = r->nspecs;
            status = code ? spec_error_of(code) : decode_specs(r, &ctx);
            for (int i = first; i < r->nspecs && !status; i++)
                status = set_render(r->specs + i, r->specs[i].magnetic_state == 2, r->specs[i].display_label);
        }
        free(lengths);
    }
    for (int i = 0; i < labels.count; i++) free(labels.items[i]);
    free(labels.items);
    for (int i = 0; i < order.count; i++) free(order.items[i]);
    free(order.items);
    free(slot_words.data);
    free(preferred.data);
    free(cases.data);
    free(directions);
    free(ops);
    return status;
}

/* ---------------------------------------------------------- entry */

static int fallback_spec(md_request *r, md_spec *fallback, md_matrix *selected_rows, int *dir_state, md_matrix *selected_dir)
{
    memset(fallback, 0, sizeof *fallback);
    const jvalue *old = jget(r->irrep, "old_id");
    int64_t old_id = 0;
    if (jtruthy(old) && py_int_of(old, &old_id)) return md_fail(MDE_VALUE, "invalid literal for int()");
    const char *keys[3] = {"ordinary_symbol", "symbol", "label"};
    const jvalue *label = 0;
    for (int i = 0; i < 3 && !label; i++)
        if (jtruthy(jget(r->irrep, keys[i]))) label = jget(r->irrep, keys[i]);
    fallback->kind = K_FALLBACK;
    fallback->old_id = old_id;
    fallback->label = label ? md_py_str(label) : (char *)md_strdup("", 0);
    if (!fallback->label) return md_fail(MDE_MEMORY, "out of memory");
    fallback->k_label = md_k_label_from_irrep_label(fallback->label);
    int has_gid;
    MD_TRY(md_dynamic_gid(r, &has_gid, &fallback->gid));
    fallback->gid_state = has_gid ? MD_SET : MD_NONE;
    isotropy_row_id(r, &fallback->row_id_state, &fallback->row_id);
    fallback->primary = 1;
    int has_kp;
    MD_TRY(source_kparam_record(r, &has_kp, fallback->kp));
    fallback->kp_state = has_kp ? MD_SET : MD_NONE;
    const jvalue *rows = selected_source_rows(r);
    int ndim = 1;
    MD_TRY(md_float_rows(rows, selected_rows, &ndim));
    fallback->rows_state = MD_SET;
    MD_TRY(matrix_copy(&fallback->rows, selected_rows));
    const jvalue *free_value = r->iso ? jget(r->iso, "free") : 0;
    fallback->free_state = MD_SET;
    fallback->free_count = 0;
    if (jtruthy(free_value) && py_int_of(free_value, &fallback->free_count)) return md_fail(MDE_VALUE, "invalid literal for int()");
    *dir_state = ndim == 2 && selected_rows->rows * selected_rows->cols ? MD_SET : MD_NONE;
    if (*dir_state == MD_SET) MD_TRY(matrix_transpose(selected_dir, selected_rows));
    fallback->dir_state = *dir_state;
    if (*dir_state == MD_SET) MD_TRY(matrix_copy(&fallback->dir, selected_dir));
    fallback->magnetic_state = MD_ABSENT;
    fallback->carrier_state = MD_ABSENT;
    fallback->pml_state = MD_ABSENT;
    fallback->req_state = MD_ABSENT;
    fallback->occ_state = MD_ABSENT;
    return 0;
}

/* _subduced_mode_specs */
static int subduced_mode_specs(md_request *r)
{
    md_spec fallback;
    md_matrix selected_rows = {0}, selected_dir = {0};
    int dir_state = MD_NONE, status = fallback_spec(r, &fallback, &selected_rows, &dir_state, &selected_dir);
    if (!status && jtruthy(jget(r->irrep, "magnetic"))) {
        status = magnetic_embedding_specs(r);
    } else if (!status) {
        int64_t basis[9], origin[4];
        int has_basis = md_integer_source_basis(r, basis), has_origin = source_split_origin(r, origin);
        int has_child = r->has_child_sg;
        spec_context ctx = {&fallback, fallback.free_count, 1, &selected_rows, dir_state, &selected_dir};
        if (fallback.old_id <= 0) {
            if (!has_child || !has_basis || !has_origin) {
                md_spec *spec = push_spec(r);
                status = spec ? copy_spec(spec, &fallback) : MDE_MEMORY;
            } else {
                int64_t kp[4] = {0, 0, 0, 1};
                int has_kp = fallback.kp_state == MD_SET;
                if (has_kp) memcpy(kp, fallback.kp, sizeof kp);
                int code = spec_ordinary_dynamic(r->sg, r->child_sg, basis, origin, fallback.gid_state == MD_SET ? (int)fallback.gid : 0,
                                                 has_kp, kp, dir_state == MD_SET);
                status = code ? spec_error_of(code) : decode_specs(r, &ctx);
            }
        } else {
            if (fallback.row_id_state != MD_SET) {
                status = md_fail(MDE_TYPE, "int() argument must be a string, a bytes-like object or a real number, not 'NoneType'");
            } else {
                int has_frame = has_child && has_basis && has_origin;
                int64_t zero_basis[9] = {0}, zero_origin[4] = {0, 0, 0, 1};
                int code = spec_static(r->sg, (int)fallback.old_id, (int)fallback.row_id, has_frame, has_frame ? r->child_sg : 0,
                                       has_frame ? basis : zero_basis, has_frame ? origin : zero_origin);
                status = code ? spec_error_of(code) : decode_specs(r, &ctx);
            }
        }
    }
    md_spec_free(&fallback);
    matrix_free(&selected_rows);
    matrix_free(&selected_dir);
    return status;
}

int md_render_specs(md_request *r)
{
    const jvalue *coupled = r->iso ? jget(r->iso, "coupled") : 0;
    r->coupled = jtruthy(coupled) && jtruthy(r->slots);
    const jvalue *include_magnetic = jget(r->request, "include_magnetic");
    if (r->coupled) {
        const jvalue *displacive_rows = jget(r->request, "displacive_row_ids"), *magnetic_rows = jget(r->request, "magnetic_row_ids");
        int include_displacive = !displacive_rows || jis(displacive_rows, J_NULL) || jtruthy(displacive_rows);
        int include_mag = jtruthy(include_magnetic) && (!magnetic_rows || jis(magnetic_rows, J_NULL) || jtruthy(magnetic_rows));
        return md_coupled_specs(r, include_displacive, include_mag);
    }
    int first = r->nspecs;
    MD_TRY(subduced_mode_specs(r));
    for (int i = first; i < r->nspecs; i++) MD_TRY(set_render(r->specs + i, 0, r->specs[i].label));
    if (jtruthy(include_magnetic) && jtruthy(jget(r->irrep, "magnetic"))) {
        first = r->nspecs;
        MD_TRY(magnetic_specs(r));
        for (int i = first; i < r->nspecs; i++)
            MD_TRY(set_render(r->specs + i, 1, r->specs[i].display_label ? r->specs[i].display_label : r->specs[i].label));
    }
    return 0;
}
