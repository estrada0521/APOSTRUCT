/* apo serve and apo show: the local GUI server and its launch client. */

#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include "net.h"
#include "platform.h"
#include "json.h"
#include "md.h"
#include "pyconv.h"

const char *cif_summary_json(const char *data, int64_t length, int from_file);
int st_state(const jvalue *request, jbuf *b);
int64_t mode_details_vectors(const double **data);
const char *apo_error(int kind, int status, const char *message, jbuf *text);

/* ------------------------------------------------------------ errors */

typedef struct {
    char message[32768];
    const char *type;
} api_error;

static int api_fail(api_error *e, const char *format, ...) __attribute__((format(printf, 2, 3)));

static int api_fail(api_error *e, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    vsnprintf(e->message, sizeof e->message, format, args);
    va_end(args);
    e->type = "ValueError";
    return 1;
}

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

static int text_order(const void *a, const void *b) { return strcmp(*(char *const *)a, *(char *const *)b); }

static int is_space(char c) { return c == ' ' || (c >= '\t' && c <= '\r') || (c >= 0x1c && c <= 0x1f); }

/* ------------------------------------------------------------ state request */

static void put_repr_list_keys(jbuf *b, const jvalue **keys, int n)
{
    jb_lit(b, "[");
    for (int i = 0; i < n; i++) {
        if (i) jb_lit(b, ", ");
        py_repr_text(b, keys[i]->u.s, keys[i]->count);
    }
    jb_lit(b, "]");
}

static int key_order(const void *a, const void *b)
{
    const jvalue *x = *(const jvalue *const *)a, *y = *(const jvalue *const *)b;
    int64_t n = x->count < y->count ? x->count : y->count;
    int c = memcmp(x->u.s, y->u.s, (size_t)n);
    return c ? c : (x->count > y->count) - (x->count < y->count);
}

/* _reject_unknown_fields */
static int reject_unknown(const jvalue *value, const char *const *allowed, int nallowed, const char *field, api_error *e)
{
    const jvalue *unknown[64];
    int n = 0;
    for (int64_t i = 0; i < value->count && n < 64; i++) {
        const jvalue *key = value->u.members + 2 * i;
        int ok = 0;
        for (int a = 0; a < nallowed && !ok; a++) ok = jstr_is(key, allowed[a]);
        if (!ok) unknown[n++] = key;
    }
    if (!n) return 0;
    qsort(unknown, (size_t)n, sizeof *unknown, key_order);
    jbuf list = {0};
    jb_reset(&list);
    put_repr_list_keys(&list, unknown, n);
    api_fail(e, "unknown %s fields: %.*s", field, (int)list.used, list.data);
    jb_free(&list);
    return 1;
}

/* _exact_text */
static int exact_text(const jvalue *v, const char *field, api_error *e)
{
    int ok = jis(v, J_STRING) && v->count && !is_space(v->u.s[0]) && !is_space(v->u.s[v->count - 1]);
    return ok ? 0 : api_fail(e, "%s must be a non-empty string without outer whitespace", field);
}

static int positive_index(const jvalue *v, const char *field, api_error *e)
{
    return jis(v, J_INT) && v->u.i >= 1 ? 0 : api_fail(e, "%s must be a positive integer", field);
}

static int bool_field(const jvalue *v, const char *field, api_error *e)
{
    return jis(v, J_TRUE) || jis(v, J_FALSE) ? 0 : api_fail(e, "%s must be boolean", field);
}

/* _rational_params, written as an object */
static int rational_params(const jvalue *v, const char *field, jbuf *b, api_error *e)
{
    jb_lit(b, "{");
    if (!v || jis(v, J_NULL)) {
        jb_lit(b, "}");
        return 0;
    }
    if (!jis(v, J_OBJECT)) return api_fail(e, "%s must be an object", field);
    int first = 1;
    for (int64_t i = 0; i < v->count; i++) {
        const jvalue *key = v->u.members + 2 * i, *value = key + 1;
        char key_field[256];
        snprintf(key_field, sizeof key_field, "%s key", field);
        if (exact_text(key, key_field, e)) return 1;
        jbuf text = {0};
        jb_reset(&text);
        if (jis(value, J_STRING)) {
            const char *s = value->u.s, *end = value->u.s + value->count;
            while (s < end && is_space(*s)) s++;
            while (end > s && is_space(end[-1])) end--;
            if (s == end) {
                jb_free(&text);
                continue;
            }
            jb_raw(&text, s, (size_t)(end - s));
        } else if (jis(value, J_INT) || (jis(value, J_FLOAT) && value->u.f == value->u.f && value->u.f - value->u.f == 0.0)) {
            py_str_of(&text, value);
        } else {
            jb_free(&text);
            return api_fail(e, "%s.%.*s must be a finite rational value", field, (int)key->count, key->u.s);
        }
        rat r;
        if (py_fraction_of_text(text.data, (int64_t)text.used, &r) == 1) {
            jb_free(&text);
            return api_fail(e, "%s.%.*s must be a finite rational value", field, (int)key->count, key->u.s);
        }
        if (!first) jb_lit(b, ", ");
        first = 0;
        jb_value(b, key);
        jb_lit(b, ": ");
        jb_string(b, text.data, text.used);
        jb_free(&text);
    }
    jb_lit(b, "}");
    return 0;
}

static int site_indexes(const jvalue *v, const char *field, jbuf *b, api_error *e)
{
    if (jis(v, J_STRING)) {
        if (!jstr_is(v, "all") && !jstr_is(v, "none")) return api_fail(e, "%s must be 'all', 'none', or a list of indexes", field);
        jb_value(b, v);
        return 0;
    }
    if (!jis(v, J_ARRAY)) return api_fail(e, "%s must be 'all', 'none', or a list of indexes", field);
    char item_field[128];
    snprintf(item_field, sizeof item_field, "%s item", field);
    for (int64_t i = 0; i < v->count; i++) {
        if (positive_index(jat(v, i), item_field, e)) return 1;
        for (int64_t k = 0; k < i; k++)
            if (jat(v, k)->u.i == jat(v, i)->u.i) return api_fail(e, "%s must not contain duplicate indexes", field);
    }
    jb_value(b, v);
    return 0;
}

static int distortion_selection(const jvalue *v, jbuf *b, api_error *e)
{
    if (!v || jis(v, J_NULL)) {
        jb_lit(b, "null");
        return 0;
    }
    if (!jis(v, J_OBJECT)) return api_fail(e, "distortion_selection must be an object");
    static const char *const allowed[] = {"strain", "displacive_sites", "magnetic_sites"};
    if (reject_unknown(v, allowed, 3, "distortion_selection", e)) return 1;
    jb_lit(b, "{");
    int first = 1;
    const jvalue *strain = jget(v, "strain");
    if (strain) {
        if (bool_field(strain, "distortion_selection.strain", e)) return 1;
        jb_lit(b, "\"strain\": ");
        jb_value(b, strain);
        first = 0;
    }
    static const char *const names[] = {"displacive_sites", "magnetic_sites"};
    for (int i = 0; i < 2; i++) {
        const jvalue *sites = jget(v, names[i]);
        if (!sites) continue;
        char field[64];
        snprintf(field, sizeof field, "distortion_selection.%s", names[i]);
        if (!first) jb_lit(b, ", ");
        first = 0;
        jb_lit(b, "\"");
        jb_lit(b, names[i]);
        jb_lit(b, "\": ");
        if (site_indexes(sites, field, b, e)) return 1;
    }
    jb_lit(b, "}");
    return 0;
}

static int selection(const jvalue *v, int ordinal, jbuf *b, api_error *e)
{
    char field[64], sub[96];
    snprintf(field, sizeof field, "selections[%d]", ordinal);
    if (!jis(v, J_OBJECT)) return api_fail(e, "%s must be an object", field);
    static const char *const allowed[] = {"k_label", "k_index", "k_params", "display_k_params", "irrep", "irrep_index"};
    if (reject_unknown(v, allowed, 6, field, e)) return 1;
    if (!jget(v, "k_label") && !jget(v, "k_index")) return api_fail(e, "%s requires k_label or k_index", field);
    jb_lit(b, "{");
    int first = 1;
    for (int i = 0; i < 6; i++) {
        const jvalue *value = jget(v, allowed[i]);
        if (!value) continue;
        snprintf(sub, sizeof sub, "%s.%s", field, allowed[i]);
        if (!first) jb_lit(b, ", ");
        first = 0;
        jb_lit(b, "\"");
        jb_lit(b, allowed[i]);
        jb_lit(b, "\": ");
        if (i == 0) {
            if (exact_text(value, sub, e)) return 1;
            jb_value(b, value);
        } else if (i == 1 || i == 5) {
            if (positive_index(value, sub, e)) return 1;
            jb_value(b, value);
        } else if (i == 2 || i == 3) {
            if (rational_params(value, sub, b, e)) return 1;
        } else {
            if (!jis(value, J_NULL) && exact_text(value, sub, e)) return 1;
            jb_value(b, value);
        }
    }
    jb_lit(b, "}");
    return 0;
}

/* APOSTRUCT._state_request -> the native pipeline request */
static int state_request(const jvalue *payload, int allow_empty, jbuf *request, jbuf *cif, int *with_text, api_error *e)
{
    if (!jis(payload, J_OBJECT)) return api_fail(e, "request body must be a JSON object");
    static const char *const allowed[] = {"cif", "mode_set", "distortion_selection", "k_label", "k_index", "k_params", "irrep",
                                          "irrep_index", "orderparam", "orderparam_index", "selections", "include_opd", "include_mode_details"};
    if (reject_unknown(payload, allowed, 13, "state", e)) return 1;
    const jvalue *text = jget(payload, "cif");
    int blank = 1;
    for (int64_t i = 0; jis(text, J_STRING) && i < text->count && blank; i++) blank = is_space(text->u.s[i]);
    if (!jis(text, J_STRING) || blank) return api_fail(e, "cif must be non-empty text");
    jb_raw(cif, text->u.s, (size_t)text->count);
    const jvalue *mode_set = jget(payload, "mode_set");
    static const jvalue default_mode_set = {J_STRING, 17, {.s = "strain_displacive"}};
    if (!mode_set) mode_set = &default_mode_set;
    if (exact_text(mode_set, "mode_set", e)) return 1;
    if (!jstr_is(mode_set, "displacive") && !jstr_is(mode_set, "strain_displacive") && !jstr_is(mode_set, "magnetic") && !jstr_is(mode_set, "all"))
        return api_fail(e, "unknown mode_set: %.*s", (int)mode_set->count, mode_set->u.s);
    jbuf selections = {0};
    jb_reset(&selections);
    const jvalue *raw = jget(payload, "selections");
    if (raw && !jis(raw, J_NULL)) {
        if (!jis(raw, J_ARRAY) || raw->count < 1 || raw->count > 9) return api_fail(e, "selections must contain one to nine ordered objects");
        jb_lit(&selections, "[");
        for (int64_t i = 0; i < raw->count; i++) {
            if (i) jb_lit(&selections, ", ");
            if (selection(jat(raw, i), (int)i + 1, &selections, e)) return 1;
        }
        jb_lit(&selections, "]");
    } else {
        jb_lit(&selections, "null");
    }
    const jvalue *orderparam = jget(payload, "orderparam");
    if (orderparam && !jis(orderparam, J_NULL)) {
        if (jis(orderparam, J_INT)) {
            if (positive_index(orderparam, "orderparam", e)) return 1;
        } else if (exact_text(orderparam, "orderparam", e)) {
            return 1;
        }
    }
    const jvalue *include_opd = jget(payload, "include_opd");
    static const jvalue truth = {J_TRUE, 0, {0}};
    if (!include_opd) include_opd = &truth;
    if (bool_field(include_opd, "include_opd", e)) return 1;
    const jvalue *include_modes = jget(payload, "include_mode_details");
    if (!include_modes) include_modes = include_opd;
    if (bool_field(include_modes, "include_mode_details", e)) return 1;
    if (jis(include_modes, J_TRUE) && !jis(include_opd, J_TRUE)) return api_fail(e, "include_mode_details requires include_opd");
    jbuf distortion = {0}, k_params = {0};
    jb_reset(&distortion);
    jb_reset(&k_params);
    if (distortion_selection(jget(payload, "distortion_selection"), &distortion, e)) return 1;
    const jvalue *k_label = jget(payload, "k_label");
    static const jvalue gamma = {J_STRING, 2, {.s = "GM"}};
    if (!k_label) k_label = &gamma;
    if (exact_text(k_label, "k_label", e)) return 1;
    static const jvalue one = {J_INT, 0, {.i = 1}};
    const jvalue *k_index = jget(payload, "k_index") ? jget(payload, "k_index") : &one;
    if (positive_index(k_index, "k_index", e)) return 1;
    if (rational_params(jget(payload, "k_params"), "k_params", &k_params, e)) return 1;
    const jvalue *irrep = jget(payload, "irrep");
    if (irrep && !jis(irrep, J_NULL) && exact_text(irrep, "irrep", e)) return 1;
    const jvalue *irrep_index = jget(payload, "irrep_index") ? jget(payload, "irrep_index") : &one;
    if (positive_index(irrep_index, "irrep_index", e)) return 1;
    const jvalue *orderparam_index = jget(payload, "orderparam_index") ? jget(payload, "orderparam_index") : &one;
    if (positive_index(orderparam_index, "orderparam_index", e)) return 1;
    *with_text = jis(include_modes, J_TRUE);
    jb_lit(request, "\"options\": {\"mode_set\": ");
    jb_value(request, mode_set);
    jb_lit(request, ", \"distortion_selection\": ");
    jb_raw(request, distortion.data, distortion.used);
    jb_lit(request, allow_empty ? ", \"allow_empty_distortion\": true" : ", \"allow_empty_distortion\": false");
    jb_lit(request, ", \"k_label\": ");
    jb_value(request, k_label);
    jb_lit(request, ", \"k_index\": ");
    jb_value(request, k_index);
    jb_lit(request, ", \"k_params\": ");
    jb_raw(request, k_params.data, k_params.used);
    jb_lit(request, ", \"display_k_params\": ");
    jb_raw(request, k_params.data, k_params.used);
    jb_lit(request, ", \"irrep\": ");
    if (irrep) jb_value(request, irrep);
    else jb_lit(request, "null");
    jb_lit(request, ", \"irrep_index\": ");
    jb_value(request, irrep_index);
    jb_lit(request, ", \"orderparam\": ");
    if (orderparam) jb_value(request, orderparam);
    else jb_lit(request, "null");
    jb_lit(request, ", \"orderparam_index\": ");
    jb_value(request, orderparam_index);
    jb_lit(request, ", \"selections\": ");
    jb_raw(request, selections.data, selections.used);
    jb_lit(request, jis(include_opd, J_TRUE) ? ", \"include_opd\": true" : ", \"include_opd\": false");
    jb_lit(request, ", \"opd_row_limit\": null, \"selected_opd_only\": false}");
    jb_lit(request, jis(include_modes, J_TRUE) ? ", \"include_mode_details\": true" : ", \"include_mode_details\": false");
    jb_lit(request, jis(include_modes, J_TRUE) ? ", \"include_mode_text\": true" : ", \"include_mode_text\": false");
    jb_free(&selections);
    jb_free(&distortion);
    jb_free(&k_params);
    return 0;
}

static void put_vectors(jbuf *b, const jvalue *placeholder, const double *vectors, int64_t nvectors)
{
    const jvalue *pair = jget(placeholder, "$f64");
    int64_t offset = jat(pair, 0)->u.i, atoms = jat(pair, 1)->u.i;
    jb_lit(b, "[");
    for (int64_t a = 0; a < atoms; a++) {
        if (a) jb_lit(b, ", ");
        jb_lit(b, "[");
        for (int k = 0; k < 3; k++) {
            if (k) jb_lit(b, ", ");
            int64_t at = offset + 3 * a + k;
            jb_float(b, at < nvectors ? vectors[at] : 0.0);
        }
        jb_lit(b, "]");
    }
    jb_lit(b, "]");
}

static void put_expanded(jbuf *b, const jvalue *v, const double *vectors, int64_t nvectors)
{
    if (jis(v, J_OBJECT) && v->count == 1 && jstr_is(v->u.members, "$f64")) {
        put_vectors(b, v, vectors, nvectors);
    } else if (jis(v, J_ARRAY)) {
        jb_lit(b, "[");
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_lit(b, ", ");
            put_expanded(b, v->u.items + i, vectors, nvectors);
        }
        jb_lit(b, "]");
    } else if (jis(v, J_OBJECT)) {
        jb_lit(b, "{");
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_lit(b, ", ");
            jb_value(b, v->u.members + 2 * i);
            jb_lit(b, ": ");
            put_expanded(b, v->u.members + 2 * i + 1, vectors, nvectors);
        }
        jb_lit(b, "}");
    } else {
        jb_value(b, v);
    }
}

static void put_mode_text(jbuf *out, const jvalue *details, const jvalue *mode_text)
{
    jb_lit(out, "\"mode_detail_text\": ");
    const jvalue *reason = member(details, "reason");
    if (jstr_is(member(details, "status"), "ok") && mode_text) {
        jb_value(out, mode_text);
    } else if (jtruthy(reason)) {
        jbuf t = {0};
        jb_reset(&t);
        py_str_of(&t, reason);
        jb_string(out, t.data, t.used);
        jb_free(&t);
    } else {
        jb_lit(out, "\"Mode details are unavailable\"");
    }
}

/* APOSTRUCT._state_result as JSON */
static int gui_state(const jvalue *payload, int allow_empty, jbuf *out, api_error *e)
{
    jbuf request = {0}, cif = {0};
    jb_reset(&request);
    jb_reset(&cif);
    int with_text = 0;
    jb_lit(&request, "{\"call\": \"state\", ");
    if (state_request(payload, allow_empty, &request, &cif, &with_text, e)) {
        jb_free(&request);
        jb_free(&cif);
        return 1;
    }
    const char *summary = cif_summary_json(cif.data, (int64_t)cif.used, 0);
    jb_free(&cif);
    jdoc *sdoc = summary ? jdoc_parse(summary, strlen(summary)) : 0;
    const jvalue *error = sdoc ? jget(jdoc_root(sdoc), "error") : 0;
    if (!sdoc || error) {
        const jvalue *message = member(error, "message");
        int64_t kind = 0;
        if (error) py_int_of(jget(error, "kind"), &kind);
        if (!sdoc) api_fail(e, "CIF reader ran out of memory");
        else api_fail(e, "%.*s", jis(message, J_STRING) ? (int)message->count : 0, jis(message, J_STRING) ? message->u.s : "");
        e->type = !sdoc ? "MemoryError" : kind == 1 ? "ValueError" : kind == 3 ? "IndexError" : "RuntimeError";
        jdoc_free(sdoc);
        jb_free(&request);
        return 1;
    }
    jb_lit(&request, ", \"cif_info\": ");
    jb_value(&request, jdoc_root(sdoc));
    jb_lit(&request, "}");
    jdoc_free(sdoc);
    jdoc *rdoc = jdoc_parse(request.data, request.used);
    jb_free(&request);
    jbuf state = {0};
    jb_reset(&state);
    int status = rdoc ? st_state(jdoc_root(rdoc), &state) : md_fail(MDE_MEMORY, "out of memory");
    jdoc_free(rdoc);
    if (status) {
        jbuf text = {0};
        jb_reset(&text);
        e->type = apo_error(md_err.kind, md_err.status, md_err.message, &text);
        snprintf(e->message, sizeof e->message, "%.*s", (int)text.used, text.data ? text.data : "");
        jb_free(&text);
        jb_free(&state);
        return 1;
    }
    jdoc *doc = jdoc_parse(state.data, state.used);
    jb_free(&state);
    if (!doc) return api_fail(e, "state is not JSON");
    const jvalue *root = jdoc_root(doc), *mode_text = jget(root, "$mode_detail_text");
    const double *vectors;
    int64_t nvectors = mode_details_vectors(&vectors);
    jb_lit(out, "{");
    int first = 1, placed = 0;
    for (int64_t i = 0; i < root->count; i++) {
        const jvalue *key = root->u.members + 2 * i, *value = key + 1;
        if (jstr_is(key, "$mode_detail_text")) continue;
        if (!first) jb_lit(out, ", ");
        first = 0;
        jb_value(out, key);
        jb_lit(out, ": ");
        if (!jstr_is(key, "selected") || !with_text) {
            put_expanded(out, value, vectors, nvectors);
            continue;
        }
        placed = 1;
        jb_lit(out, "{");
        int n = jis(value, J_OBJECT) ? (int)value->count : 0;
        for (int k = 0; k < n; k++) {
            if (k) jb_lit(out, ", ");
            jb_value(out, value->u.members + 2 * k);
            jb_lit(out, ": ");
            put_expanded(out, value->u.members + 2 * k + 1, vectors, nvectors);
        }
        if (n) jb_lit(out, ", ");
        put_mode_text(out, member(value, "mode_details"), mode_text);
        jb_lit(out, "}");
    }
    if (with_text && !placed) {
        jb_lit(out, first ? "\"selected\": {" : ", \"selected\": {");
        put_mode_text(out, 0, 0);
        jb_lit(out, "}");
    }
    jb_lit(out, "}");
    jdoc_free(doc);
    return out->failed ? api_fail(e, "out of memory") : 0;
}

/* ------------------------------------------------------------ GUI launch */

static pthread_mutex_t compute_lock = PTHREAD_MUTEX_INITIALIZER, mailbox_lock = PTHREAD_MUTEX_INITIALIZER;
static int64_t revision;
static char *launch_json;

/* gui_launch.publish_gui_launch */
static int publish(const jvalue *launch, jbuf *out, api_error *e)
{
    static const char *const keys[] = {"schema", "stage", "source_name", "state_request", "amplitudes"};
    int shape = jis(launch, J_OBJECT) && launch->count == 5;
    for (int i = 0; i < 5 && shape; i++) shape = jget(launch, keys[i]) != 0;
    if (!shape) return api_fail(e, "GUI launch has an invalid shape");
    if (!jstr_is(jget(launch, "schema"), "APOSTRUCT.gui.launch.v1")) return api_fail(e, "GUI launch has an unknown schema");
    const jvalue *stage = jget(launch, "stage");
    if (!jstr_is(stage, "kpoints") && !jstr_is(stage, "irreps") && !jstr_is(stage, "opds") && !jstr_is(stage, "modes"))
        return api_fail(e, "GUI launch has an unknown stage");
    const jvalue *source = jget(launch, "source_name");
    int blank = 1;
    for (int64_t i = 0; jis(source, J_STRING) && i < source->count && blank; i++) blank = is_space(source->u.s[i]);
    if (!jis(source, J_STRING) || blank) return api_fail(e, "GUI launch source_name must be non-empty text");
    const jvalue *request = jget(launch, "state_request"), *amplitudes = jget(launch, "amplitudes");
    if (!jis(request, J_OBJECT)) return api_fail(e, "GUI launch state_request must be an object");
    if (!jis(amplitudes, J_ARRAY)) return api_fail(e, "GUI launch amplitudes must be an array");
    jbuf normalized = {0};
    jb_reset(&normalized);
    jb_lit(&normalized, "[");
    for (int64_t i = 0; i < amplitudes->count; i++) {
        const jvalue *row = jat(amplitudes, i), *id = jget(row, "definition_id"), *value = jget(row, "value");
        if (!jis(row, J_OBJECT) || row->count != 2 || !id || !value) {
            jb_free(&normalized);
            return api_fail(e, "GUI launch amplitude has an invalid shape");
        }
        if (!jis(id, J_STRING) || !id->count) {
            jb_free(&normalized);
            return api_fail(e, "GUI launch amplitude definition_id must be non-empty text");
        }
        for (int64_t k = 0; k < i; k++) {
            const jvalue *other = jget(jat(amplitudes, k), "definition_id");
            if (other->count == id->count && !memcmp(other->u.s, id->u.s, (size_t)id->count)) {
                jb_free(&normalized);
                return api_fail(e, "duplicate amplitude for definition %.*s", (int)id->count, id->u.s);
            }
        }
        jbuf text = {0};
        jb_reset(&text);
        py_str_of(&text, value);
        rat r;
        int bad = py_fraction_of_text(text.data ? text.data : "", (int64_t)text.used, &r) != 0;
        jb_free(&text);
        if (bad) {
            jb_free(&normalized);
            return api_fail(e, "amplitude for definition %.*s must be an exact number", (int)id->count, id->u.s);
        }
        char exact[96];
        py_rat_text(r, exact);
        if (i) jb_lit(&normalized, ", ");
        jb_lit(&normalized, "{\"definition_id\": ");
        jb_value(&normalized, id);
        jb_lit(&normalized, ", \"value\": ");
        jb_cstring(&normalized, exact);
        jb_lit(&normalized, "}");
    }
    jb_lit(&normalized, "]");
    if (amplitudes->count && !jstr_is(stage, "modes")) {
        jb_free(&normalized);
        return api_fail(e, "GUI launch amplitudes require the modes stage");
    }
    jbuf state = {0};
    jb_reset(&state);
    int status = gui_state(request, jstr_is(stage, "kpoints"), &state, e);
    if (status) {
        jb_free(&normalized);
        jb_free(&state);
        return 1;
    }
    if (amplitudes->count) {
        jdoc *doc = jdoc_parse(state.data, state.used);
        const jvalue *details = member(member(jdoc_root(doc), "selected"), "mode_details");
        if (!jstr_is(member(details, "status"), "ok")) {
            jdoc_free(doc);
            jb_free(&normalized);
            jb_free(&state);
            return api_fail(e, "GUI launch returned no mode definitions for amplitudes");
        }
        static const char *const kinds[3][2] = {{"displacive", "displacive_definitions"}, {"magnetic", "magnetic_definitions"}, {"strain", "strain_definitions"}};
        int navailable = 0;
        for (int k = 0; k < 3; k++) {
            const jvalue *list = member(details, kinds[k][1]);
            navailable += jis(list, J_ARRAY) ? (int)list->count : 0;
        }
        char **available = calloc((size_t)navailable + 1, sizeof *available);
        int filled = 0;
        for (int k = 0; k < 3; k++) {
            const jvalue *list = member(details, kinds[k][1]);
            for (int64_t i = 0; jis(list, J_ARRAY) && i < list->count; i++) {
                available[filled] = malloc(48);
                snprintf(available[filled++], 48, "%s-%lld", kinds[k][0], (long long)i + 1);
            }
        }
        qsort(available, (size_t)filled, sizeof *available, text_order);
        jbuf unknown = {0};
        jb_reset(&unknown);
        int nunknown = 0;
        for (int64_t i = 0; i < amplitudes->count; i++) {
            const jvalue *id = jget(jat(amplitudes, i), "definition_id");
            int found = 0;
            for (int k = 0; k < filled && !found; k++) found = (int64_t)strlen(available[k]) == id->count && !memcmp(available[k], id->u.s, (size_t)id->count);
            if (found) continue;
            jb_lit(&unknown, nunknown++ ? ", " : "[");
            py_repr_text(&unknown, id->u.s, id->count);
        }
        if (nunknown) {
            jb_lit(&unknown, "]; available IDs: [");
            for (int k = 0; k < filled; k++) {
                if (k) jb_lit(&unknown, ", ");
                py_repr_text(&unknown, available[k], (int64_t)strlen(available[k]));
            }
            jb_lit(&unknown, "]");
            api_fail(e, "mode definitions are unavailable in this case: %.*s", (int)unknown.used, unknown.data);
        }
        for (int k = 0; k < filled; k++) free(available[k]);
        free(available);
        jb_free(&unknown);
        if (nunknown) {
            jdoc_free(doc);
            jb_free(&normalized);
            jb_free(&state);
            return 1;
        }
        jdoc_free(doc);
    }
    jbuf prepared = {0};
    jb_reset(&prepared);
    const char *s = source->u.s, *end = source->u.s + source->count;
    while (s < end && is_space(*s)) s++;
    while (end > s && is_space(end[-1])) end--;
    jb_lit(&prepared, "{\"schema\": \"APOSTRUCT.gui.launch.v1\", \"stage\": ");
    jb_value(&prepared, stage);
    jb_lit(&prepared, ", \"source_name\": ");
    jb_string(&prepared, s, (size_t)(end - s));
    jb_lit(&prepared, ", \"state_request\": ");
    jb_value(&prepared, request);
    jb_lit(&prepared, ", \"amplitudes\": ");
    jb_raw(&prepared, normalized.data, normalized.used);
    jb_lit(&prepared, ", \"state\": ");
    jb_raw(&prepared, state.data, state.used);
    jb_lit(&prepared, "}");
    jb_free(&normalized);
    jb_free(&state);
    pthread_mutex_lock(&mailbox_lock);
    free(launch_json);
    launch_json = plat_strndup(prepared.data, prepared.used);
    int64_t current = ++revision;
    pthread_mutex_unlock(&mailbox_lock);
    jb_free(&prepared);
    char text[96];
    snprintf(text, sizeof text, "{\"revision\": %lld, \"stage\": ", (long long)current);
    jb_lit(out, text);
    jb_value(out, stage);
    jb_lit(out, "}");
    return 0;
}

static void read_mailbox(int64_t after, jbuf *out)
{
    pthread_mutex_lock(&mailbox_lock);
    char text[64];
    snprintf(text, sizeof text, "{\"revision\": %lld, \"launch\": ", (long long)revision);
    jb_lit(out, text);
    if (launch_json && revision > after) jb_lit(out, launch_json);
    else jb_lit(out, "null");
    jb_lit(out, "}");
    pthread_mutex_unlock(&mailbox_lock);
}

/* ------------------------------------------------------------ HTTP */

static char frontend[PATH_MAX];
static volatile sig_atomic_t stopping;
static net_socket listener = NET_INVALID;

typedef struct {
    net_socket fd;
} connection;

static void send_all(net_socket fd, const char *data, size_t n)
{
    while (n) {
        ssize_t sent = send(fd, data, n, 0);
        if (sent <= 0) return;
        data += sent;
        n -= (size_t)sent;
    }
}

static void respond(net_socket fd, int status, const char *type, const char *body, size_t n)
{
    const char *reason = status == 200 ? "OK" : status == 400 ? "Bad Request" : status == 404 ? "Not Found" : status == 413 ? "Payload Too Large" : "Internal Server Error";
    char head[512];
    int len = snprintf(head, sizeof head,
                       "HTTP/1.0 %d %s\r\nContent-Type: %s\r\nContent-Length: %zu\r\nCache-Control: no-store\r\nConnection: close\r\n\r\n", status,
                       reason, type, n);
    send_all(fd, head, (size_t)len);
    send_all(fd, body, n);
}

static void respond_json(net_socket fd, int status, const jbuf *body) { respond(fd, status, "application/json; charset=utf-8", body->data, body->used); }

static void respond_error(net_socket fd, const api_error *e)
{
    jbuf body = {0}, text = {0};
    jb_reset(&body);
    jb_reset(&text);
    jb_lit(&text, e->type);
    jb_lit(&text, ": ");
    jb_lit(&text, e->message);
    jb_lit(&body, "{\"error\": ");
    jb_string(&body, text.data, text.used);
    jb_lit(&body, "}");
    static const char *const caught[] = {"ValueError", "KeyError", "IndexError", "TypeError", "OpdNotFoundError", "LinAlgError"};
    int status = 500;
    for (int i = 0; i < 6 && status == 500; i++)
        if (!strcmp(e->type, caught[i])) status = 400;
    respond_json(fd, status, &body);
    jb_free(&body);
    jb_free(&text);
}

static const char *content_type(const char *path)
{
    const char *dot = strrchr(path, '.');
    if (!dot) return "application/octet-stream";
    static const char *const types[][2] = {{".csv", "text/csv; charset=utf-8"},         {".css", "text/css; charset=utf-8"},
                                           {".html", "text/html; charset=utf-8"},       {".js", "text/javascript; charset=utf-8"},
                                           {".json", "application/json; charset=utf-8"}, {".svg", "image/svg+xml"}};
    for (size_t i = 0; i < sizeof types / sizeof types[0]; i++)
        if (!strcasecmp(dot, types[i][0])) return types[i][1];
    return "application/octet-stream";
}

/* server._asset_path */
static int serve_asset(net_socket fd, const char *path)
{
    const char *relative;
    if (!strcmp(path, "/") || !strcmp(path, "/index.html") || !strcmp(path, "/tools/APOSTRUCT") || !strcmp(path, "/tools/APOSTRUCT/")) relative = "index.html";
    else if (!strncmp(path, "/tools/APOSTRUCT/", 17)) relative = path + 17;
    else relative = path + 1;
    char joined[PATH_MAX * 2], real[PATH_MAX];
    snprintf(joined, sizeof joined, "%s/%s", frontend, relative);
    struct stat st;
    if (!plat_realpath(joined, real) || strncmp(real, frontend, strlen(frontend)) || (real[strlen(frontend)] && real[strlen(frontend)] != '/') ||
        stat(real, &st) || !S_ISREG(st.st_mode))
        return 0;
    FILE *f = fopen(real, "rb");
    if (!f) return 0;
    char *data = malloc((size_t)st.st_size + 1);
    size_t n = fread(data, 1, (size_t)st.st_size, f);
    fclose(f);
    respond(fd, 200, content_type(real), data, n);
    free(data);
    return 1;
}

static int query_after(const char *query, int64_t *after, api_error *e)
{
    *after = 0;
    if (!query || !*query) return 0;
    if (strncmp(query, "after=", 6)) {
        char name[128];
        snprintf(name, sizeof name, "%.*s", (int)strcspn(query, "="), query);
        return api_fail(e, "unknown gui_launch filters: ['%s']", name);
    }
    const char *value = query + 6;
    char *end;
    long long n = strtoll(value, &end, 10);
    if (end == value || (*end && *end != '&') || n < 0) return api_fail(e, "gui_launch after must be one non-negative integer");
    if (*end == '&') return api_fail(e, "gui_launch after must be one non-negative integer");
    *after = n;
    return 0;
}

static void handle(net_socket fd)
{
    jbuf in = {0};
    jb_reset(&in);
    char chunk[65536];
    size_t header_end = 0, content_length = 0;
    for (;;) {
        ssize_t got = recv(fd, chunk, sizeof chunk, 0);
        if (got <= 0) break;
        jb_raw(&in, chunk, (size_t)got);
        if (!header_end) {
            char *mark = in.data ? strstr(in.data, "\r\n\r\n") : 0;
            if (mark) {
                header_end = (size_t)(mark - in.data) + 4;
                for (char *line = in.data; line < mark; line = strstr(line, "\r\n") + 2)
                    if (!strncasecmp(line, "Content-Length:", 15)) content_length = strtoull(line + 15, 0, 10);
                if (content_length > (size_t)16 << 20) {
                    const char *body = "{\"error\": \"request body is too large\"}";
                    respond(fd, 413, "application/json; charset=utf-8", body, strlen(body));
                    jb_free(&in);
                    return;
                }
            }
        }
        if (header_end && in.used >= header_end + content_length) break;
    }
    if (!header_end) {
        jb_free(&in);
        return;
    }
    char method[16] = "", target[4096] = "";
    sscanf(in.data, "%15s %4095s", method, target);
    char *query = strchr(target, '?');
    if (query) *query++ = 0;
    const char *body = in.data + header_end;
    jbuf out = {0};
    jb_reset(&out);
    api_error e = {"", "ValueError"};
    if (!strncmp(target, "/api/APOSTRUCT/", 15)) {
        const char *endpoint = target + 15;
        int failed = 0;
        if (!strcmp(method, "GET")) {
            if (!strcmp(endpoint, "gui_launch")) {
                int64_t after;
                failed = query_after(query, &after, &e);
                if (!failed) read_mailbox(after, &out);
            } else {
                failed = api_fail(&e, "unknown endpoint: %s", endpoint);
            }
        } else if (!strcmp(method, "POST")) {
            pthread_mutex_lock(&compute_lock);
            jdoc *doc = jdoc_parse(content_length ? body : "{}", content_length ? content_length : 2);
            if (!doc) {
                failed = api_fail(&e, "request body is not JSON");
            } else if (!jis(jdoc_root(doc), J_OBJECT)) {
                failed = api_fail(&e, "request body must be a JSON object");
            } else if (!strcmp(endpoint, "state")) {
                failed = gui_state(jdoc_root(doc), 0, &out, &e);
            } else if (!strcmp(endpoint, "gui_launch")) {
                failed = query && *query ? api_fail(&e, "gui_launch does not accept query parameters") : publish(jdoc_root(doc), &out, &e);
            } else {
                failed = api_fail(&e, "unknown endpoint: %s", endpoint);
            }
            jdoc_free(doc);
            pthread_mutex_unlock(&compute_lock);
        } else {
            respond(fd, 404, "text/plain; charset=utf-8", "Not Found", 9);
            jb_free(&in);
            jb_free(&out);
            return;
        }
        if (failed) respond_error(fd, &e);
        else respond_json(fd, 200, &out);
    } else if (!strcmp(method, "POST") && !strcmp(target, "/shutdown")) {
        respond(fd, 200, "text/plain; charset=utf-8", "ok", 2);
        net_close(fd);
        jb_free(&in);
        _exit(0);
    } else if (!strcmp(method, "GET") && serve_asset(fd, target)) {
    } else {
        respond(fd, 404, "text/plain; charset=utf-8", "Not Found", 9);
    }
    jb_free(&in);
    jb_free(&out);
}

static void *connection_main(void *arg)
{
    net_socket fd = ((connection *)arg)->fd;
    free(arg);
    handle(fd);
    net_close(fd);
    return 0;
}

int apo_serve(const char *host, int port, int open_browser, const char *frontend_dir, char *error, size_t error_size)
{
    if (port < 0 || port > 65535) {
        snprintf(error, error_size, "port must be an integer from 0 through 65535");
        return 1;
    }
    if (!plat_realpath(frontend_dir, frontend)) {
        snprintf(error, error_size, "the frontend is missing: %s", frontend_dir);
        return 1;
    }
    if (net_init()) {
        snprintf(error, error_size, "cannot start the network library");
        return 1;
    }
    struct addrinfo hints = {0}, *address;
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;
    char port_text[16];
    snprintf(port_text, sizeof port_text, "%d", port);
    if (getaddrinfo(host, port_text, &hints, &address)) {
        snprintf(error, error_size, "[Errno 8] nodename nor servname provided, or not known");
        return 1;
    }
    listener = socket(address->ai_family, SOCK_STREAM, 0);
    int yes = 1;
    setsockopt(listener, SOL_SOCKET, SO_REUSEADDR, (const char *)&yes, sizeof yes);
    if (listener == NET_INVALID || bind(listener, address->ai_addr, address->ai_addrlen) || listen(listener, 64)) {
        snprintf(error, error_size, "[Errno %d] %s", errno, strerror(errno));
        freeaddrinfo(address);
        return 1;
    }
    freeaddrinfo(address);
    struct sockaddr_in bound;
    socklen_t length = sizeof bound;
    getsockname(listener, (struct sockaddr *)&bound, &length);
    char actual[INET_ADDRSTRLEN] = "";
    inet_ntop(AF_INET, &bound.sin_addr, actual, sizeof actual);
    char url[256];
    snprintf(url, sizeof url, "http://%s:%d/", actual, ntohs(bound.sin_port));
    printf("APOSTRUCT: %s\nStop: Ctrl+C or POST /shutdown\n", url);
    fflush(stdout);
    if (open_browser) {
        char command[600];
#ifdef _WIN32
        snprintf(command, sizeof command, "start \"\" \"%s\"", url);
#else
        snprintf(command, sizeof command, "(open '%s' || xdg-open '%s') >/dev/null 2>&1 &", url, url);
#endif
        if (system(command)) {
        }
    }
    while (!stopping) {
        net_socket fd = accept(listener, 0, 0);
        if (fd == NET_INVALID) {
            if (stopping) break;
            if (errno == EINTR) continue;
            break;
        }
        connection *c = malloc(sizeof *c);
        c->fd = fd;
        pthread_t thread;
        if (pthread_create(&thread, 0, connection_main, c)) {
            net_close(fd);
            free(c);
            continue;
        }
        pthread_detach(thread);
    }
    net_close(listener);
    return 0;
}

/* ------------------------------------------------------------ show */

/* gui_launch.send_gui_launch: POST the launch, return the reply or an error. */
int apo_post(const char *server, const char *path, const jbuf *body, jbuf *reply, char *base, size_t base_size, char *error, size_t error_size)
{
    const char *s = server, *e = server + strlen(server);
    while (s < e && is_space(*s)) s++;
    while (e > s && is_space(e[-1])) e--;
    while (e > s && e[-1] == '/') e--;
    snprintf(base, base_size, "%.*s", (int)(e - s), s);
    if (!*base) {
        snprintf(error, error_size, "--server must be a non-empty URL");
        return 1;
    }
    if (strncmp(base, "http://", 7)) {
        snprintf(error, error_size, "unknown url type: '%s%s'", base, path);
        return 1;
    }
    const char *rest = base + 7;
    char host[512] = "", port[16] = "80";
    const char *colon = strchr(rest, ':'), *slash = strchr(rest, '/');
    size_t host_len = colon && (!slash || colon < slash) ? (size_t)(colon - rest) : slash ? (size_t)(slash - rest) : strlen(rest);
    snprintf(host, sizeof host, "%.*s", (int)host_len, rest);
    if (colon && (!slash || colon < slash)) snprintf(port, sizeof port, "%.*s", (int)strcspn(colon + 1, "/"), colon + 1);
    char target[2048];
    snprintf(target, sizeof target, "%s%s", slash ? slash : "", path);
    struct addrinfo hints = {0}, *address;
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;
    net_socket fd = NET_INVALID;
    if (!net_init() && !getaddrinfo(host, port, &hints, &address)) {
        for (struct addrinfo *a = address; a && fd == NET_INVALID; a = a->ai_next) {
            fd = socket(a->ai_family, SOCK_STREAM, 0);
            if (fd != NET_INVALID && connect(fd, a->ai_addr, a->ai_addrlen)) {
                net_close(fd);
                fd = NET_INVALID;
            }
        }
        freeaddrinfo(address);
    }
    if (fd == NET_INVALID) {
        snprintf(error, error_size, "could not reach APOSTRUCT GUI server at %s", base);
        return 1;
    }
    char head[1024];
    int len = snprintf(head, sizeof head, "POST %s HTTP/1.0\r\nHost: %s\r\nContent-Type: application/json\r\nContent-Length: %zu\r\n\r\n", target, host,
                       body->used);
    send_all(fd, head, (size_t)len);
    send_all(fd, body->data, body->used);
    jbuf in = {0};
    jb_reset(&in);
    char chunk[65536];
    for (;;) {
        ssize_t got = recv(fd, chunk, sizeof chunk, 0);
        if (got <= 0) break;
        jb_raw(&in, chunk, (size_t)got);
    }
    net_close(fd);
    int code = 0;
    sscanf(in.data ? in.data : "", "HTTP/%*s %d", &code);
    char *mark = in.data ? strstr(in.data, "\r\n\r\n") : 0;
    const char *payload = mark ? mark + 4 : "";
    if (code != 200) {
        jdoc *doc = jdoc_parse(payload, strlen(payload));
        const jvalue *message = doc ? jget(jdoc_root(doc), "error") : 0;
        if (jis(message, J_STRING)) snprintf(error, error_size, "%.*s", (int)message->count, message->u.s);
        else snprintf(error, error_size, "GUI server returned HTTP %d", code);
        jdoc_free(doc);
        jb_free(&in);
        return 1;
    }
    jb_lit(reply, payload);
    jb_free(&in);
    return 0;
}
