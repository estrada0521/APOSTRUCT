
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "par.h"
#include "sel.h"

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

static int64_t member_int(const jvalue *object, const char *key)
{
    int64_t v = 0;
    const jvalue *m = member(object, key);
    if (jtruthy(m)) py_int_of(m, &v);
    return v;
}

typedef struct {
    int n;
    int *items;
} ints_t;
static int site_selection(const jvalue *requested, int site_count, int default_all, ints_t *out)
{
    out->n = 0;
    out->items = calloc((size_t)site_count + 1, sizeof(int));
    if (!out->items) return md_fail(MDE_MEMORY, "out of memory");
    if (!requested || jis(requested, J_NULL) || jstr_is(requested, "all")) {
        if (!requested || jis(requested, J_NULL) ? default_all : 1)
            for (int i = 1; i <= site_count; i++) out->items[out->n++] = i;
        return 0;
    }
    if (jstr_is(requested, "none")) return 0;
    int64_t m = jis(requested, J_ARRAY) ? requested->count : 0;
    int64_t *values = calloc((size_t)m + 1, sizeof(int64_t));
    if (!values) return md_fail(MDE_MEMORY, "out of memory");
    int64_t unique = 0;
    for (int64_t i = 0; i < m; i++) {
        int64_t v;
        if (py_int_of(jat(requested, i), &v)) {
            free(values);
            return md_fail(MDE_VALUE, "invalid literal for int() with base 10");
        }
        int known = 0;
        for (int64_t k = 0; k < unique && !known; k++) known = values[k] == v;
        if (!known) values[unique++] = v;
    }
    for (int64_t a = 0; a < unique; a++)
        for (int64_t b = a + 1; b < unique; b++)
            if (values[b] < values[a]) {
                int64_t t = values[a];
                values[a] = values[b];
                values[b] = t;
            }
    jbuf invalid = {0};
    jb_reset(&invalid);
    int bad = 0;
    for (int64_t i = 0; i < unique; i++)
        if (values[i] < 1 || values[i] > site_count) {
            jb_lit(&invalid, bad++ ? ", " : "[");
            jb_int(&invalid, values[i]);
        }
    if (bad) {
        jb_lit(&invalid, "]");
        int status = md_fail(MDE_VALUE, "site indexes out of range 1..%d: %.*s", site_count, (int)invalid.used, invalid.data);
        jb_free(&invalid);
        free(values);
        return status;
    }
    jb_free(&invalid);
    free(out->items);
    out->items = calloc((size_t)unique + 1, sizeof(int));
    for (int64_t i = 0; out->items && i < unique; i++) out->items[out->n++] = (int)values[i];
    free(values);
    return out->items ? 0 : md_fail(MDE_MEMORY, "out of memory");
}

static void put_ints(jbuf *b, const ints_t *v)
{
    jb_lit(b, "[");
    for (int i = 0; i < v->n; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, v->items[i]);
    }
    jb_lit(b, "]");
}

static void put_row_ids(jbuf *b, const jvalue *atoms, const ints_t *v)
{
    int64_t rows[512];
    int n = 0;
    for (int i = 0; i < v->n && n < 512; i++) {
        int64_t row = member_int(jat(atoms, v->items[i] - 1), "wyckoff_row_id");
        if (!row) continue;
        int known = 0;
        for (int k = 0; k < n && !known; k++) known = rows[k] == row;
        if (!known) rows[n++] = row;
    }
    for (int a = 0; a < n; a++)
        for (int c = a + 1; c < n; c++)
            if (rows[c] < rows[a]) {
                int64_t t = rows[a];
                rows[a] = rows[c];
                rows[c] = t;
            }
    jb_lit(b, "[");
    for (int i = 0; i < n; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, rows[i]);
    }
    jb_lit(b, "]");
}

static int mode_set_is(const jvalue *mode_set, const char *a, const char *b) { return jstr_is(mode_set, a) || (b && jstr_is(mode_set, b)); }
int par_distortion(const jvalue *state, const jvalue *mode_set, const jvalue *requested, int allow_empty, jbuf *out)
{
    const jvalue *atoms = member(member(state, "input"), "atom_sites");
    int n = jis(atoms, J_ARRAY) ? (int)atoms->count : 0;
    jbuf types = {0}, labels = {0};
    jb_reset(&types);
    jb_reset(&labels);
    jb_lit(&types, "[");
    jb_lit(&labels, "[");
    for (int i = 0; i < n; i++) {
        const jvalue *atom = jat(atoms, i), *type = member(atom, "type"), *label = member(atom, "label");
        if (i) {
            jb_lit(&types, ", ");
            jb_lit(&labels, ", ");
        }
        jbuf t = {0};
        jb_reset(&t);
        if (jtruthy(type)) py_str_of(&t, type);
        else if (jtruthy(label)) py_str_of(&t, label);
        else jb_int(&t, i + 1);
        jb_string(&types, t.data, t.used);
        jb_reset(&t);
        if (jtruthy(label)) py_str_of(&t, label);
        else if (jtruthy(type)) py_str_of(&t, type);
        else jb_int(&t, i + 1);
        jb_string(&labels, t.data, t.used);
        jb_free(&t);
    }
    jb_lit(&types, "]");
    jb_lit(&labels, "]");
    int include_strain;
    if (jis(requested, J_OBJECT) && jget(requested, "strain")) include_strain = jtruthy(jget(requested, "strain"));
    else include_strain = mode_set_is(mode_set, "strain_displacive", "all");
    ints_t displacive = {0}, magnetic = {0};
    const jvalue *dsites = member(requested, "displacive_sites"), *msites = member(requested, "magnetic_sites");
    int status = site_selection(dsites, n, 1, &displacive);
    if (!status) {
        if (jstr_is(msites, "all")) status = site_selection(msites, n, 0, &magnetic);
        else if ((!msites || jis(msites, J_NULL)) && mode_set_is(mode_set, "magnetic", "all")) status = site_selection(0, n, 1, &magnetic);
        else status = site_selection(msites, n, 0, &magnetic);
    }
    for (int pass = 0; pass < 2 && !status; pass++) {
        const ints_t *selected = pass ? &magnetic : &displacive;
        const char *name = pass ? "magnetic" : "displacive";
        for (int i = 0; i < selected->n; i++) {
            const jvalue *site = jat(atoms, selected->items[i] - 1);
            if (member_int(site, "wyckoff_row_id")) continue;
            jbuf label = {0}, repr = {0}, coordinates = {0}, reason = {0};
            jb_reset(&label);
            jb_reset(&repr);
            jb_reset(&coordinates);
            jb_reset(&reason);
            const jvalue *l = member(site, "label"), *t = member(site, "type");
            if (jtruthy(l)) py_str_of(&label, l);
            else if (jtruthy(t)) py_str_of(&label, t);
            else {
                char text[32];
                snprintf(text, sizeof text, "site %d", selected->items[i]);
                jb_lit(&label, text);
            }
            py_repr_text(&repr, label.data, (int64_t)label.used);
            const jvalue *fract = member(site, "fract");
            if (fract) py_repr_of(&coordinates, fract);
            else jb_lit(&coordinates, "None");
            const jvalue *error = member(site, "wyckoff_mapping_error");
            if (jtruthy(error)) py_str_of(&reason, error);
            else jb_lit(&reason, "Wyckoff position was not resolved");
            status = md_fail(MDE_VALUE, "cannot select %s site %.*s at %.*s: %.*s", name, (int)repr.used, repr.data, (int)coordinates.used,
                             coordinates.data, (int)reason.used, reason.data);
            jb_free(&label);
            jb_free(&repr);
            jb_free(&coordinates);
            jb_free(&reason);
            break;
        }
    }
    if (!status && !allow_empty && !include_strain && !displacive.n && !magnetic.n) status = md_fail(MDE_VALUE, "You must include some kind of distortion");
    if (!status) {
        ints_t selectable = {0};
        site_selection(0, n, 1, &selectable);
        jbuf selection = {0}, internal = {0}, a = {0}, b = {0}, c = {0}, d = {0}, mode = {0};
        jb_reset(&selection);
        jb_reset(&internal);
        jb_reset(&a);
        jb_reset(&b);
        jb_reset(&c);
        jb_reset(&d);
        jb_reset(&mode);
        if (mode_set) jb_value(&mode, mode_set);
        else jb_lit(&mode, "null");
        jb_lit(&selection, include_strain ? "{\"strain\": true, \"displacive_sites\": " : "{\"strain\": false, \"displacive_sites\": ");
        put_ints(&selection, &displacive);
        jb_lit(&selection, ", \"magnetic_sites\": ");
        put_ints(&selection, &magnetic);
        jb_lit(&selection, ", \"selectable_sites\": ");
        put_ints(&selection, &selectable);
        jb_lit(&selection, ", \"site_types\": ");
        jb_raw(&selection, types.data, types.used);
        jb_lit(&selection, ", \"site_labels\": ");
        jb_raw(&selection, labels.data, labels.used);
        jb_lit(&selection, "}");
        put_row_ids(&a, atoms, &displacive);
        put_row_ids(&b, atoms, &magnetic);
        put_ints(&c, &displacive);
        put_ints(&d, &magnetic);
        kv internal_items[6] = {{"mode_set", mode.data, mode.used},
                                {"distortion_selection", selection.data, selection.used},
                                {"displacive_row_ids", a.data, a.used},
                                {"magnetic_row_ids", b.data, b.used},
                                {"displacive_site_indexes", c.data, c.used},
                                {"magnetic_site_indexes", d.data, d.used}};
        sel_merge(&internal, member(state, "_internal"), internal_items, 6);
        kv items[3] = {{"mode_set", mode.data, mode.used},
                       {"distortion_selection", selection.data, selection.used},
                       {"_internal", internal.data, internal.used}};
        sel_merge(out, state, items, 3);
        free(selectable.items);
        jb_free(&selection);
        jb_free(&internal);
        jb_free(&a);
        jb_free(&b);
        jb_free(&c);
        jb_free(&d);
        jb_free(&mode);
    }
    free(displacive.items);
    free(magnetic.items);
    jb_free(&types);
    jb_free(&labels);
    return status;
}
