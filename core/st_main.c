/* The pipeline entry: pipeline.build_state_from_cif_info and
 * pipeline.attach_mode_details in one call. The parent state, the selection
 * stage and the mode details run in C; mode-detail vectors stay in the
 * mode_details_vectors side buffer. A mode-detail failure is returned in
 * place as {"status": "error", "reason": "Type: text"}. */

#include <string.h>

#include "md.h"
#include "par.h"
#include "sel.h"

int sel_stage(const jvalue *state, const jvalue *opts, jbuf *b);
int md_render_text(const jvalue *state, const double *vectors, int64_t nvectors, jbuf *b);
int64_t mode_details_vectors(const double **data);
int md_details(const jvalue *request, jbuf *b);
void exact_rows_reset(void);

static jbuf out;

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

static void put_or_null(jbuf *b, const jvalue *v)
{
    if (v) jb_value(b, v);
    else jb_lit(b, "null");
}

void md_reason(int kind, int status, const char *message, jbuf *out);

static void put_error_object(jbuf *b)
{
    jbuf reason = {0};
    jb_reset(&reason);
    md_reason(md_err.kind, md_err.status, md_err.message, &reason);
    jb_lit(b, "{\"status\": \"error\", \"reason\": ");
    jb_string(b, reason.data, reason.used);
    jb_lit(b, "}");
    jb_free(&reason);
}

/* pipeline.attach_mode_details: state with selected.mode_details */
static int attach(const jvalue *state, const jvalue *k_params, jbuf *b)
{
    const jvalue *selected = member(state, "selected"), *internal = member(state, "_internal"),
                 *distortion = member(state, "distortion_selection");
    jbuf request = {0}, details = {0};
    jb_reset(&request);
    jb_reset(&details);
    jb_lit(&request, "{\"cif_info\": ");
    put_or_null(&request, member(state, "input"));
    static const char *selected_keys[][2] = {{"selected_k", "kpoint"},
                                             {"selected_irrep", "irrep"},
                                             {"selected_opd", "orderparam"},
                                             {"selected_slots", "selection_slots"},
                                             {"selected_factors", "opd_factors"}};
    for (size_t i = 0; i < sizeof selected_keys / sizeof selected_keys[0]; i++) {
        jb_lit(&request, ", \"");
        jb_lit(&request, selected_keys[i][0]);
        jb_lit(&request, "\": ");
        put_or_null(&request, member(selected, selected_keys[i][1]));
    }
    jb_lit(&request, ", \"k_params\": ");
    put_or_null(&request, k_params);
    jb_lit(&request, jtruthy(member(distortion, "strain")) ? ", \"include_strain\": true" : ", \"include_strain\": false");
    jb_lit(&request, jtruthy(member(distortion, "magnetic_sites")) ? ", \"include_magnetic\": true" : ", \"include_magnetic\": false");
    static const char *internal_keys[] = {"displacive_row_ids", "magnetic_row_ids", "displacive_site_indexes", "magnetic_site_indexes",
                                          "parent_inter_setting_id"};
    for (size_t i = 0; i < sizeof internal_keys / sizeof internal_keys[0]; i++) {
        jb_lit(&request, ", \"");
        jb_lit(&request, internal_keys[i]);
        jb_lit(&request, "\": ");
        put_or_null(&request, member(internal, internal_keys[i]));
    }
    jb_lit(&request, "}");
    if (request.failed) {
        jb_free(&request);
        return md_fail(MDE_MEMORY, "out of memory");
    }
    jdoc *doc = jdoc_parse(request.data, request.used);
    jb_free(&request);
    if (!doc) return md_fail(MDE_MEMORY, "out of memory");
    if (md_details(jdoc_root(doc), &details)) {
        jb_reset(&details);
        put_error_object(&details);
    }
    jdoc_free(doc);
    jb_lit(b, "{");
    for (int64_t i = 0; i < state->count; i++) {
        const jvalue *key = state->u.members + 2 * i, *value = key + 1;
        if (i) jb_lit(b, ", ");
        jb_value(b, key);
        jb_lit(b, ": ");
        if (!jstr_is(key, "selected") || !jis(value, J_OBJECT)) {
            jb_value(b, value);
            continue;
        }
        jb_lit(b, "{");
        int first = 1;
        for (int64_t k = 0; k < value->count; k++) {
            const jvalue *sk = value->u.members + 2 * k;
            if (jstr_is(sk, "mode_details")) continue;
            if (!first) jb_lit(b, ", ");
            first = 0;
            jb_value(b, sk);
            jb_lit(b, ": ");
            jb_value(b, sk + 1);
        }
        if (!first) jb_lit(b, ", ");
        jb_lit(b, "\"mode_details\": ");
        jb_raw(b, details.data, details.used);
        jb_lit(b, "}");
    }
    jb_lit(b, "}");
    jb_free(&details);
    return b->failed ? md_fail(MDE_MEMORY, "out of memory") : 0;
}

/* Appends "$mode_detail_text" to a state whose mode details are ok. */
static int add_mode_text(jbuf *b)
{
    jdoc *doc = jdoc_parse(b->data, b->used);
    if (!doc) return md_fail(MDE_MEMORY, "out of memory");
    const jvalue *root = jdoc_root(doc), *details = member(member(root, "selected"), "mode_details");
    int status = 0;
    if (jstr_is(member(details, "status"), "ok")) {
        jbuf text = {0};
        jb_reset(&text);
        const double *vectors;
        int64_t nvectors = mode_details_vectors(&vectors);
        status = md_render_text(root, vectors, nvectors, &text);
        if (!status && b->used && b->data[b->used - 1] == '}') {
            b->used--;
            jb_lit(b, ", \"$mode_detail_text\": ");
            jb_string(b, text.data, text.used);
            jb_lit(b, "}");
        }
        jb_free(&text);
    }
    jdoc_free(doc);
    return status ? status : b->failed ? md_fail(MDE_MEMORY, "out of memory") : 0;
}

int st_state(const jvalue *request, jbuf *b)
{
    exact_rows_reset();
    const jvalue *options = jget(request, "options");
    jbuf parent = {0}, selection = {0};
    jb_reset(&parent);
    jb_reset(&selection);
    int status = par_parent_state(jget(request, "cif_info"), &parent);
    jdoc *pdoc = status ? 0 : jdoc_parse(parent.data, parent.used);
    if (!status && !pdoc) status = md_fail(MDE_MEMORY, "out of memory");
    int attach_modes = jtruthy(member(options, "include_opd")) && jtruthy(jget(request, "include_mode_details"));
    if (!status) status = sel_stage(jdoc_root(pdoc), options, attach_modes ? &selection : b);
    if (!status && attach_modes) {
        jdoc *sdoc = jdoc_parse(selection.data, selection.used);
        if (!sdoc) status = md_fail(MDE_MEMORY, "out of memory");
        else status = attach(jdoc_root(sdoc), member(options, "k_params"), b);
        jdoc_free(sdoc);
    }
    if (!status && attach_modes && jtruthy(jget(request, "include_mode_text"))) status = add_mode_text(b);
    jdoc_free(pdoc);
    jb_free(&parent);
    jb_free(&selection);
    return status;
}

/* {"call": "state", "cif_info", "options", "include_mode_details"} or
 * {"call": "attach", "state", "k_params"}. */
const char *state_json(const char *text, int64_t length)
{
    if (out.cap > ((size_t)16 << 20)) jb_free(&out);
    jb_reset(&out);
    jdoc *doc = jdoc_parse(text, (size_t)length);
    int status;
    if (!doc) {
        status = md_fail(MDE_VALUE, "state request is not JSON");
    } else {
        const jvalue *request = jdoc_root(doc), *call = jget(request, "call");
        if (jstr_is(call, "state")) status = st_state(request, &out);
        else if (jstr_is(call, "attach")) status = attach(jget(request, "state"), jget(request, "k_params"), &out);
        else status = md_fail(MDE_VALUE, "unknown state call");
    }
    if (status) {
        jb_reset(&out);
        jb_lit(&out, "{\"error\": {\"kind\": ");
        jb_int(&out, md_err.kind);
        jb_lit(&out, ", \"status\": ");
        jb_int(&out, md_err.status);
        jb_lit(&out, ", \"message\": ");
        jb_cstring(&out, md_err.message);
        jb_lit(&out, "}}");
    }
    jdoc_free(doc);
    return out.failed ? 0 : out.data;
}
