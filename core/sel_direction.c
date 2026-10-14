/* isotropy.selection.select_order_parameter_direction: the OPD rows of the
 * selected irrep (or coupled irreps), the selected row and its factors. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"
#include "sel.h"
#include "tables.h"

int k_coordinate_map_words(int sg, int kslot, int64_t *out);

static struct {
    section_view image_label, image_dimension, image_order, image_type, image_landau, image_count, image_pointer, irrep_image, mag2nonmag;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (sel_tables() || sel_mag_tables() || section_find("images.image_label", &T.image_label) ||
        section_find("images.image_dimension", &T.image_dimension) || section_find("images.image_order", &T.image_order) ||
        section_find("images.image_type", &T.image_type) || section_find("images.image_landau", &T.image_landau) ||
        section_find("images.image_subgroup_count", &T.image_count) || section_find("images.image_subgroup_pointer", &T.image_pointer) ||
        section_find("irreps.irrep_image", &T.irrep_image) || section_find("magnetic.mag_point_op_mag2nonmag", &T.mag2nonmag))
        return md_fail(MDE_VALUE, "selection Source tables are inconsistent");
    T.ready = 1;
    return 0;
}

/* ------------------------------------------------------------ JSON helpers */

static int key_is(const jvalue *key, const char *text) { return (size_t)key->count == strlen(text) && !memcmp(key->u.s, text, (size_t)key->count); }

/* {**base, key: value, ...} */
void sel_merge(jbuf *b, const jvalue *base, const kv *items, int n)
{
    jb_lit(b, "{");
    int first = 1, *used = calloc((size_t)n + 1, sizeof(int));
    for (int64_t i = 0; jis(base, J_OBJECT) && i < base->count; i++) {
        const jvalue *key = base->u.members + 2 * i, *value = key + 1;
        if (!first) jb_lit(b, ", ");
        first = 0;
        jb_string(b, key->u.s, (size_t)key->count);
        jb_lit(b, ": ");
        int replaced = -1;
        for (int k = 0; k < n; k++)
            if (key_is(key, items[k].key)) replaced = k;
        if (replaced >= 0) {
            jb_raw(b, items[replaced].json, items[replaced].length);
            if (used) used[replaced] = 1;
        } else {
            jb_value(b, value);
        }
    }
    for (int k = 0; k < n; k++) {
        if (used && used[k]) continue;
        if (!first) jb_lit(b, ", ");
        first = 0;
        jb_cstring(b, items[k].key);
        jb_lit(b, ": ");
        jb_raw(b, items[k].json, items[k].length);
    }
    jb_lit(b, "}");
    free(used);
}

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

/* str(value or "") */
static void str_or_empty(jbuf *b, const jvalue *v)
{
    if (jtruthy(v)) py_str_of(b, v);
}

static void strip_bounds(const char **a, const char **e)
{
    while (*a < *e && sel_is_space(**a)) (*a)++;
    while (*e > *a && sel_is_space((*e)[-1])) (*e)--;
}

static int text_equals(const jbuf *t, const char *s, size_t n) { return t->used == n && !memcmp(t->data, s, n); }

/* str(x or "").strip() equals label */
static int stripped_is(const jvalue *v, const char *label, size_t n)
{
    jbuf t = {0};
    jb_reset(&t);
    str_or_empty(&t, v);
    const char *a = t.data, *e = t.data + t.used;
    strip_bounds(&a, &e);
    int same = (size_t)(e - a) == n && !memcmp(a, label, n);
    jb_free(&t);
    return same;
}

static int truthy_int(const jvalue *v, int64_t fallback, int64_t *out)
{
    if (!jtruthy(v)) {
        *out = fallback;
        return 0;
    }
    return py_int_of(v, out) ? md_fail(MDE_VALUE, "invalid literal for int()") : 0;
}

/* reciprocal.k_points.normalize_k_params as JSON text */
void sel_normalize_params(jbuf *b, const jvalue *v)
{
    jb_lit(b, "{");
    if (jtruthy(v) && jis(v, J_OBJECT)) {
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_lit(b, ", ");
            const jvalue *key = v->u.members + 2 * i;
            jb_string(b, key->u.s, (size_t)key->count);
            jb_lit(b, ": ");
            jbuf t = {0};
            jb_reset(&t);
            py_str_of(&t, key + 1);
            jb_string(b, t.data, t.used);
            jb_free(&t);
        }
    } else if (jtruthy(v) && jis(v, J_ARRAY)) {
        int64_t n = v->count, count = 0;
        jbuf *texts = calloc((size_t)n + 1, sizeof(jbuf));
        const char **ka = calloc((size_t)n + 1, sizeof(char *)), **ke = calloc((size_t)n + 1, sizeof(char *)),
                   **va = calloc((size_t)n + 1, sizeof(char *)), **ve = calloc((size_t)n + 1, sizeof(char *));
        for (int64_t i = 0; texts && ka && ke && va && ve && i < n; i++) {
            jb_reset(&texts[i]);
            py_str_of(&texts[i], jat(v, i));
            char *eq = texts[i].used ? memchr(texts[i].data, '=', texts[i].used) : 0;
            if (!eq) continue;
            const char *k0 = texts[i].data, *k1 = eq, *v0 = eq + 1, *v1 = texts[i].data + texts[i].used;
            strip_bounds(&k0, &k1);
            strip_bounds(&v0, &v1);
            int64_t slot = -1;
            for (int64_t j = 0; j < count && slot < 0; j++)
                if (ke[j] - ka[j] == k1 - k0 && !memcmp(ka[j], k0, (size_t)(k1 - k0))) slot = j;
            if (slot < 0) {
                slot = count++;
                ka[slot] = k0;
                ke[slot] = k1;
            }
            va[slot] = v0;
            ve[slot] = v1;
        }
        for (int64_t j = 0; j < count; j++) {
            if (j) jb_lit(b, ", ");
            jb_string(b, ka[j], (size_t)(ke[j] - ka[j]));
            jb_lit(b, ": ");
            jb_string(b, va[j], (size_t)(ve[j] - va[j]));
        }
        for (int64_t i = 0; texts && i < n; i++) jb_free(&texts[i]);
        free(texts);
        free(ka);
        free(ke);
        free(va);
        free(ve);
    }
    jb_lit(b, "}");
}

static int params_empty(const jbuf *json) { return json->used <= 2; }

/* ------------------------------------------------------------ selection */

/* reciprocal.irreps.select_irrep */
static int select_irrep(const jvalue *irreps, const jvalue *label, int64_t index, const jvalue **out)
{
    int64_t n = jis(irreps, J_ARRAY) ? irreps->count : 0;
    if (jtruthy(label)) {
        jbuf want = {0};
        jb_reset(&want);
        py_str_of(&want, label);
        for (int64_t i = 0; i < n; i++) {
            jbuf have = {0};
            jb_reset(&have);
            py_str_of(&have, member(jat(irreps, i), "symbol"));
            int same = have.used == want.used;
            for (size_t k = 0; same && k < have.used; k++) {
                char x = have.data[k], y = want.data[k];
                if (x >= 'A' && x <= 'Z') x = (char)(x - 'A' + 'a');
                if (y >= 'A' && y <= 'Z') y = (char)(y - 'A' + 'a');
                same = x == y;
            }
            jb_free(&have);
            if (same) {
                jb_free(&want);
                *out = jat(irreps, i);
                return 0;
            }
        }
        jbuf repr = {0};
        jb_reset(&repr);
        py_repr_text(&repr, want.data, (int64_t)want.used);
        int status = md_fail(MDE_KEY, "irrep %.*s not found", (int)repr.used, repr.data);
        jb_free(&repr);
        jb_free(&want);
        return status;
    }
    if (index < 1 || index > n) return md_fail(MDE_INDEX, "irrep index %lld out of range 1..%lld", (long long)index, (long long)n);
    *out = jat(irreps, index - 1);
    return 0;
}

static int64_t member_int(const jvalue *object, const char *key)
{
    int64_t v = 0;
    const jvalue *m = member(object, key);
    if (jtruthy(m)) py_int_of(m, &v);
    return v;
}

typedef struct {
    int has_setting;
    int64_t setting_id;
    int has_cell;
    double cell[6];
} sel_context;

static int rows_for(int sg, const jvalue *irrep, const sel_params *k, const sel_params *display, const jvalue *selected, const sel_context *c,
                    jbuf *out)
{
    int64_t old_id = member_int(irrep, "old_id"), gid = member_int(irrep, "gid");
    int magnetic = jtruthy(member(irrep, "magnetic"));
    if (old_id <= 0) return sel_dynamic_rows(magnetic, sg, gid, k, display, selected, c->has_setting, c->setting_id, c->has_cell, c->cell, out);
    if (magnetic)
        return sel_magnetic_static_rows((int)old_id, 1, sg, 1, gid, display, selected, c->has_setting, c->setting_id, c->has_cell, c->cell, out);
    return sel_static_rows((int)old_id, k, display, c->has_setting, c->setting_id, c->has_cell, c->cell, out);
}

static void repr_value(jbuf *b, const jvalue *v)
{
    if (jis(v, J_STRING)) py_repr_text(b, v->u.s, v->count);
    else py_repr_of(b, v);
}

/* selection._select_orderparam: *out NULL for None */
static int select_orderparam(const jvalue *rows, const jvalue *orderparam, int64_t orderparam_index, const jvalue **out)
{
    *out = 0;
    int64_t n = jis(rows, J_ARRAY) ? rows->count : 0;
    int requested = orderparam && !jis(orderparam, J_NULL);
    if (jis(orderparam, J_STRING)) {
        const char *a = orderparam->u.s, *e = a + orderparam->count;
        strip_bounds(&a, &e);
        if (a < e) {
            for (int64_t i = 0; i < n; i++) {
                const jvalue *row = jat(rows, i);
                if (stripped_is(member(member(row, "isotropy"), "opd_label"), a, (size_t)(e - a)) ||
                    stripped_is(member(member(row, "direction"), "label"), a, (size_t)(e - a))) {
                    *out = row;
                    return 0;
                }
            }
        }
    }
    int has_id = 0;
    int64_t row_id = 0;
    if (jis(orderparam, J_INT)) {
        has_id = 1;
        row_id = orderparam->u.i;
    } else if (jis(orderparam, J_STRING)) {
        const char *a = orderparam->u.s, *e = a + orderparam->count;
        strip_bounds(&a, &e);
        int digits = a < e;
        for (const char *c = a; c < e && digits; c++) digits = *c >= '0' && *c <= '9';
        if (digits && e - a < 18) {
            has_id = 1;
            row_id = strtoll(a, 0, 10);
        }
    }
    if (has_id)
        for (int64_t i = 0; i < n; i++)
            if (member_int(member(jat(rows, i), "isotropy"), "row_id") == row_id) {
                *out = jat(rows, i);
                return 0;
            }
    if (requested) {
        jbuf message = {0};
        jb_reset(&message);
        jb_lit(&message, "OPD ");
        repr_value(&message, orderparam);
        jb_lit(&message, " not found; available labels: [");
        jbuf *labels = calloc((size_t)n + 1, sizeof(jbuf));
        int count = 0;
        for (int64_t i = 0; labels && i < n; i++) {
            jbuf t = {0};
            jb_reset(&t);
            str_or_empty(&t, member(member(jat(rows, i), "direction"), "label"));
            const char *a = t.data, *e = t.data + t.used;
            strip_bounds(&a, &e);
            int known = a == e;
            for (int k = 0; k < count && !known; k++) known = text_equals(&labels[k], a, (size_t)(e - a));
            if (!known) {
                jb_reset(&labels[count]);
                jb_raw(&labels[count], a, (size_t)(e - a));
                if (count) jb_lit(&message, ", ");
                py_repr_text(&message, labels[count].data ? labels[count].data : "", (int64_t)labels[count].used);
                count++;
            }
            jb_free(&t);
        }
        for (int k = 0; k < count; k++) jb_free(&labels[k]);
        free(labels);
        jb_lit(&message, "]");
        int status = md_fail(MDE_OPD_NOT_FOUND, "%.*s", (int)message.used, message.data ? message.data : "");
        jb_free(&message);
        return status;
    }
    if (orderparam_index >= 1 && orderparam_index <= n) *out = jat(rows, orderparam_index - 1);
    return 0;
}

/* selection._domain_display */
static void domain_display(jbuf *b, const jvalue *iso, const char *label, size_t label_length, int64_t domain, const char *display_opd,
                           size_t display_length, const jbuf *basis, const jbuf *origin)
{
    const jvalue *subgroup = member(iso, "subgroup");
    jbuf first = {0}, second = {0};
    jb_reset(&first);
    jb_reset(&second);
    const jvalue *display_label = member(subgroup, "display_label");
    str_or_empty(&first, jtruthy(display_label) ? display_label : member(subgroup, "number"));
    str_or_empty(&second, member(subgroup, "symbol"));
    jb_raw(b, label, label_length);
    char number[32];
    snprintf(number, sizeof number, "(%lld)", (long long)domain);
    jb_lit(b, number);
    if (display_length) {
        jb_lit(b, " ");
        for (size_t i = 0; i < display_length; i++) jb_raw(b, display_opd[i] == ';' ? "," : display_opd + i, 1);
    }
    if (first.used || second.used) {
        jb_lit(b, " ");
        jb_raw(b, first.data, first.used);
        if (first.used && second.used) jb_lit(b, " ");
        jb_raw(b, second.data, second.used);
    }
    jbuf basis_text = {0}, origin_text = {0};
    jb_reset(&basis_text);
    jb_reset(&origin_text);
    if (basis) jb_raw(&basis_text, basis->data, basis->used);
    else str_or_empty(&basis_text, member(iso, "basis_text"));
    if (origin) jb_raw(&origin_text, origin->data, origin->used);
    else str_or_empty(&origin_text, member(iso, "origin"));
    if (basis_text.used) {
        jb_lit(b, ", basis={");
        jb_raw(b, basis_text.data, basis_text.used);
        jb_lit(b, "}");
    }
    if (origin_text.used) {
        jb_lit(b, ", origin=");
        jb_raw(b, origin_text.data, origin_text.used);
    }
    jb_free(&first);
    jb_free(&second);
    jb_free(&basis_text);
    jb_free(&origin_text);
}

static int k_parameter_names(int sg, int64_t kslot, char *names, int *count)
{
    int64_t map[29];
    int status = k_coordinate_map_words(sg, (int)kslot, map);
    if (status) return md_fail(status == 3 ? MDE_ARITH : MDE_VALUE, "K-family display expression is not affine rational");
    *count = (int)map[1];
    for (int i = 0; i < *count; i++) names[i] = (char)map[2 + i];
    return 0;
}

static const jvalue *row_for_label(const jvalue *slot, const char *label, size_t n)
{
    const jvalue *rows = member(slot, "opd_rows");
    for (int64_t i = 0; jis(rows, J_ARRAY) && i < rows->count; i++) {
        const jvalue *row = jat(rows, i);
        jbuf t = {0};
        jb_reset(&t);
        str_or_empty(&t, member(member(row, "isotropy"), "opd_label"));
        int same = text_equals(&t, label, n);
        jb_free(&t);
        if (same) return row;
    }
    return 0;
}

/* selection._selected_opd_factors */
int sel_opd_factors(jbuf *b, int sg, const jvalue *slots, const jvalue *selected)
{
    const jvalue *iso = member(selected, "isotropy");
    if (!jis(selected, J_OBJECT) || !jis(iso, J_OBJECT)) {
        jb_lit(b, "[]");
        return 0;
    }
    int64_t n = slots->count;
    jbuf *labels = calloc((size_t)n + 1, sizeof(jbuf)), *displays = calloc((size_t)n + 1, sizeof(jbuf)),
         *domain_displays = calloc((size_t)n + 1, sizeof(jbuf));
    int64_t *domains = calloc((size_t)n + 1, sizeof(int64_t)), *counts = calloc((size_t)n + 1, sizeof(int64_t));
    int *valid = calloc((size_t)n + 1, sizeof(int));
    if (!labels || !displays || !domain_displays || !domains || !counts || !valid) return md_fail(MDE_MEMORY, "out of memory");
    int status = 0, ok = 1;
    if (n == 1) {
        jb_reset(&labels[0]);
        jb_reset(&displays[0]);
        jb_reset(&domain_displays[0]);
        str_or_empty(&labels[0], member(iso, "opd_label"));
        int64_t domain = 1, count = 0;
        truthy_int(member(iso, "selected_domain"), 1, &domain);
        truthy_int(member(iso, "free"), 0, &count);
        domains[0] = domain;
        counts[0] = count;
        valid[0] = 1;
        const jvalue *display = member(iso, "display_opd");
        str_or_empty(&displays[0], jtruthy(display) ? display : member(iso, "source_opd"));
        domain_display(&domain_displays[0], iso, labels[0].data ? labels[0].data : "", labels[0].used, domain,
                       displays[0].data ? displays[0].data : "", displays[0].used, 0, 0);
    } else {
        const jvalue *components = member(member(selected, "direction"), "components");
        int64_t m = jis(components, J_ARRAY) ? components->count : 0;
        for (int64_t i = 0; i < m && !status; i++) {
            const jvalue *component = jat(components, i);
            const jvalue *slot = member(component, "slot");
            if (!jis(component, J_OBJECT) || !jis(slot, J_INT) || slot->u.i != i + 1)
                status = md_fail(MDE_VALUE, "selected OPD components do not match primary slots");
        }
        if (!status && m != n) ok = 0;
        for (int64_t i = 0; !status && ok && i < n; i++) {
            const jvalue *component = jat(components, i);
            jb_reset(&labels[i]);
            jb_reset(&displays[i]);
            jb_reset(&domain_displays[i]);
            str_or_empty(&labels[i], member(component, "opd"));
            const jvalue *domain = member(component, "domain");
            valid[i] = jis(domain, J_INT);
            domains[i] = valid[i] ? domain->u.i : 0;
            const jvalue *row = row_for_label(jat(slots, i), labels[i].data ? labels[i].data : "", labels[i].used);
            const jvalue *free = row ? member(member(row, "isotropy"), "free") : 0;
            counts[i] = jis(free, J_INT) ? free->u.i : 0;
            if (!jis(free, J_INT)) valid[i] = 0;
            str_or_empty(&displays[i], member(component, "display_opd"));
            str_or_empty(&domain_displays[i], member(component, "domain_display"));
        }
    }
    for (int64_t i = 0; !status && ok && i < n; i++)
        if (!labels[i].used || !displays[i].used || !domain_displays[i].used || !valid[i] || domains[i] < 1 || counts[i] < 1) ok = 0;
    if (!status && !ok) status = md_fail(MDE_VALUE, "selected OPD does not retain one exact displayed domain per primary slot");
    if (!status) jb_lit(b, "[");
    for (int64_t i = 0; !status && i < n; i++) {
        const jvalue *slot = jat(slots, i), *irrep = member(slot, "irrep"), *kpoint = member(slot, "kpoint");
        if (i) jb_lit(b, ", ");
        jb_lit(b, "{\"slot\": ");
        jb_int(b, member_int(slot, "slot"));
        jb_lit(b, ", \"gid\": ");
        jb_int(b, member_int(irrep, "gid"));
        jb_lit(b, ", \"label\": ");
        jbuf t = {0};
        jb_reset(&t);
        str_or_empty(&t, member(irrep, "symbol"));
        jb_string(b, t.data, t.used);
        jb_free(&t);
        jb_lit(b, jtruthy(member(irrep, "magnetic")) ? ", \"magnetic\": true" : ", \"magnetic\": false");
        jb_lit(b, ", \"k_parameters\": [");
        int64_t dimension = member_int(kpoint, "dimension");
        char names[3];
        int count = 0;
        if (dimension > 0) {
            status = k_parameter_names(sg, member_int(kpoint, "kslot"), names, &count);
            if (!status && count != dimension) status = md_fail(MDE_VALUE, "selected K slot does not retain its parameters");
            const jvalue *values = member(slot, "display_k_params");
            if (!jtruthy(values)) values = member(slot, "k_params");
            for (int k = 0; !status && k < count; k++) {
                char name[2] = {names[k], 0};
                jbuf v = {0};
                jb_reset(&v);
                const jvalue *value = member(values, name);
                if (jtruthy(value)) py_str_of(&v, value);
                else jb_lit(&v, "0");
                const char *a = v.data, *e = v.data + v.used;
                strip_bounds(&a, &e);
                if (k) jb_lit(b, ", ");
                if (a < e) jb_string(b, a, (size_t)(e - a));
                else jb_lit(b, "\"0\"");
                jb_free(&v);
            }
        }
        if (status) break;
        jb_lit(b, "], \"opd\": ");
        jb_string(b, labels[i].data, labels[i].used);
        jb_lit(b, ", \"domain\": ");
        jb_int(b, domains[i]);
        jb_lit(b, ", \"domain_display\": ");
        jb_string(b, domain_displays[i].data, domain_displays[i].used);
        jb_lit(b, ", \"parameter_count\": ");
        jb_int(b, counts[i]);
        if (dimension <= 0) jb_lit(b, ", \"source_k_parameters\": []");
        if (dimension > 0) {
            const jvalue *source = member(slot, "source_kparam");
            if (!source || jis(source, J_NULL)) source = 0;
            if (!source && n == 1) source = member(iso, "source_kparam");
            if (source && jis(source, J_NULL)) source = 0;
            if (!source) {
                const jvalue *row = row_for_label(slot, labels[i].data ? labels[i].data : "", labels[i].used);
                source = member(member(row, "isotropy"), "source_kparam");
            }
            if (jis(source, J_ARRAY) && source->count == 4) {
                int64_t values[4];
                int bad = 0;
                for (int k = 0; k < 4; k++) bad |= py_int_of(jat(source, k), values + k);
                if (!bad && values[3]) {
                    jb_lit(b, ", \"source_k_parameters\": [");
                    for (int k = 0; k < dimension; k++) {
                        if (k) jb_lit(b, ", ");
                        rat v = rat_make(values[k], values[3]);
                        char text[96];
                        text[0] = '"';
                        {
                            int n = py_rat_text(v, text + 1);
                            text[n + 1] = '"';
                            text[n + 2] = 0;
                        }
                        jb_lit(b, text);
                    }
                    jb_lit(b, "]");
                }
            }
        }
        jb_lit(b, "}");
    }
    if (!status) jb_lit(b, "]");
    for (int64_t i = 0; i < n; i++) {
        jb_free(&labels[i]);
        jb_free(&displays[i]);
        jb_free(&domain_displays[i]);
    }
    free(labels);
    free(displays);
    free(domain_displays);
    free(domains);
    free(counts);
    free(valid);
    return status;
}

static void image_record(jbuf *b, int64_t old_id)
{
    int64_t image = section_int(&T.irrep_image, old_id - 1), i = image - 1;
    jb_lit(b, "{\"id\": ");
    jb_int(b, image);
    jb_lit(b, ", \"label\": ");
    const char *s;
    int64_t n;
    if (section_string(&T.image_label, i, &s, &n)) {
        s = "";
        n = 0;
    }
    while (n && (*s == ' ' || *s == '\t')) {
        s++;
        n--;
    }
    while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    jb_string(b, s, (size_t)n);
    jb_lit(b, ", \"dimension\": ");
    jb_int(b, section_int(&T.image_dimension, i));
    jb_lit(b, ", \"order\": ");
    jb_int(b, section_int(&T.image_order, i));
    jb_lit(b, ", \"type\": ");
    if (i < section_count(&T.image_type)) jb_int(b, section_int(&T.image_type, i));
    else jb_lit(b, "null");
    jb_lit(b, ", \"landau\": ");
    jb_int(b, section_int(&T.image_landau, i));
    jb_lit(b, ", \"direction_count\": ");
    jb_int(b, section_int(&T.image_count, i));
    jb_lit(b, ", \"subgroup_pointer\": ");
    jb_int(b, section_int(&T.image_pointer, i));
    jb_lit(b, "}");
}

static void context_of(const jvalue *state, sel_context *c)
{
    const jvalue *internal = member(state, "_internal");
    c->has_setting = 1;
    c->setting_id = member_int(internal, "parent_inter_setting_id");
    const jvalue *lattice = member(member(state, "input"), "lattice");
    static const char *keys[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
    c->has_cell = 1;
    for (int i = 0; i < 6; i++)
        if (float_cif_number(member(lattice, keys[i]), c->cell + i)) c->has_cell = 0;
}

static int single(const jvalue *state, const jvalue *opts, int sg, jbuf *out)
{
    jbuf kp = {0}, dp = {0};
    jb_reset(&kp);
    jb_reset(&dp);
    sel_normalize_params(&kp, member(opts, "k_params"));
    sel_normalize_params(&dp, member(opts, "display_k_params"));
    if (params_empty(&dp)) {
        jb_reset(&dp);
        jb_raw(&dp, kp.data, kp.used);
    }
    jdoc *kdoc = jdoc_parse(kp.data, kp.used), *ddoc = jdoc_parse(dp.data, dp.used);
    sel_params k_params, display;
    sel_params_of(jdoc_root(kdoc), &k_params);
    sel_params_of(jdoc_root(ddoc), &display);
    const jvalue *selections = member(opts, "selections"), *selection = 0;
    if (jtruthy(selections)) selection = jat(selections, 0);
    const jvalue *irrep_label = member(selection, "irrep");
    if (!selection || !jget(selection, "irrep")) irrep_label = member(opts, "irrep");
    int64_t default_index = 1, index;
    truthy_int(member(opts, "irrep_index"), 1, &default_index);
    int status = truthy_int(member(selection, "irrep_index"), default_index, &index);
    const jvalue *irrep = 0;
    jbuf rows = {0}, block = {0};
    jb_reset(&rows);
    jb_reset(&block);
    if (!status) status = select_irrep(member(member(state, "selected"), "irreps"), irrep_label, index, &irrep);
    sel_context c;
    context_of(state, &c);
    const jvalue *orderparam = member(opts, "orderparam");
    int selected_only = jtruthy(member(opts, "selected_opd_only"));
    if (!status) status = rows_for(sg, irrep, &k_params, &display, selected_only ? orderparam : 0, &c, &rows);
    jdoc *rdoc = status ? 0 : jdoc_parse(rows.data, rows.used);
    if (!status && !rdoc) status = md_fail(MDE_RUNTIME, "OPD rows are not JSON");
    const jvalue *selected = 0;
    int64_t orderparam_index = jis(member(opts, "orderparam_index"), J_INT) ? member(opts, "orderparam_index")->u.i : 1;
    if (!status) status = select_orderparam(jdoc_root(rdoc), orderparam, orderparam_index, &selected);
    jbuf slots = {0}, factors = {0}, image = {0}, sel_json = {0};
    jb_reset(&slots);
    jb_reset(&factors);
    jb_reset(&image);
    jb_reset(&sel_json);
    if (!status) {
        jb_lit(&slots, "[{\"slot\": 1, \"kpoint\": ");
        const jvalue *kpoint = member(member(state, "selected"), "kpoint");
        if (kpoint) jb_value(&slots, kpoint);
        else jb_lit(&slots, "null");
        jb_lit(&slots, ", \"irrep\": ");
        jb_value(&slots, irrep);
        jb_lit(&slots, ", \"k_params\": ");
        jb_raw(&slots, kp.data, kp.used);
        jb_lit(&slots, ", \"display_k_params\": ");
        jb_raw(&slots, dp.data, dp.used);
        jb_lit(&slots, "}]");
        jdoc *sdoc = jdoc_parse(slots.data, slots.used);
        status = sel_opd_factors(&factors, sg, jdoc_root(sdoc), selected);
        jdoc_free(sdoc);
    }
    if (!status) {
        int64_t old_id = member_int(irrep, "old_id");
        if (old_id <= 0) jb_lit(&image, "null");
        else image_record(&image, old_id);
        jbuf selected_text = {0};
        jb_reset(&selected_text);
        if (selected) jb_value(&selected_text, selected);
        else jb_lit(&selected_text, "null");
        jbuf irrep_text = {0};
        jb_reset(&irrep_text);
        jb_value(&irrep_text, irrep);
        kv items[6] = {{"irrep", irrep_text.data, irrep_text.used},       {"selection_slots", slots.data, slots.used},
                       {"opd_factors", factors.data, factors.used},      {"image", image.data, image.used},
                       {"opd_rows", rows.data, rows.used},               {"orderparam", selected_text.data, selected_text.used}};
        sel_merge(&sel_json, member(state, "selected"), items, 6);
        kv top[2] = {{"schema", "\"APOSTRUCT.opd_direction.v1\"", strlen("\"APOSTRUCT.opd_direction.v1\"")}, {"selected", sel_json.data, sel_json.used}};
        sel_merge(out, state, top, 2);
        jb_free(&selected_text);
        jb_free(&irrep_text);
    }
    jdoc_free(rdoc);
    jdoc_free(kdoc);
    jdoc_free(ddoc);
    jb_free(&kp);
    jb_free(&dp);
    jb_free(&rows);
    jb_free(&block);
    jb_free(&slots);
    jb_free(&factors);
    jb_free(&image);
    jb_free(&sel_json);
    return status;
}

/* selection._coupled_selected_labels: count of labels (into starts/lengths) or 0 for None */
static int coupled_labels(const jvalue *orderparam, int slot_count, int selected_only, const char **starts, int *lengths)
{
    if (!selected_only || !jis(orderparam, J_STRING) || slot_count < 2) return 0;
    const char *a = orderparam->u.s, *e = a + orderparam->count;
    strip_bounds(&a, &e);
    int parts = 0;
    const char *c = a;
    while (c < e) {
        const char *start = c;
        if (*c == 'P' || *c == 'C' || *c == 'S') {
            c++;
        } else if (c + 1 < e && *c >= '4' && *c <= '9' && c[1] == 'D') {
            c += 2;
        } else if (c + 2 < e && *c >= '1' && *c <= '9' && c[1] >= '0' && c[1] <= '9' && c[2] == 'D') {
            c += 3;
        } else {
            return 0;
        }
        if (c >= e || *c < '1' || *c > '9') return 0;
        while (c < e && *c >= '0' && *c <= '9') c++;
        const char *label_end = c;
        if (c >= e || *c != '(') return 0;
        c++;
        if (c >= e || *c < '1' || *c > '9') return 0;
        while (c < e && *c >= '0' && *c <= '9') c++;
        if (c >= e || *c != ')') return 0;
        c++;
        if (parts >= 16) return 0;
        starts[parts] = start;
        lengths[parts] = (int)(label_end - start);
        parts++;
    }
    return parts == slot_count ? parts : 0;
}

/* The coupled rows with their components' internal keys dropped and the
 * selected row's domain displays added. */
static void put_coupled_row(jbuf *b, const jvalue *row, const jbuf *displays, int n)
{
    jb_lit(b, "{");
    for (int64_t i = 0; i < row->count; i++) {
        const jvalue *key = row->u.members + 2 * i, *value = key + 1;
        if (i) jb_lit(b, ", ");
        jb_string(b, key->u.s, (size_t)key->count);
        jb_lit(b, ": ");
        if (!key_is(key, "direction")) {
            jb_value(b, value);
            continue;
        }
        jb_lit(b, "{");
        for (int64_t j = 0; j < value->count; j++) {
            const jvalue *dkey = value->u.members + 2 * j, *dvalue = dkey + 1;
            if (j) jb_lit(b, ", ");
            jb_string(b, dkey->u.s, (size_t)dkey->count);
            jb_lit(b, ": ");
            if (!key_is(dkey, "components")) {
                jb_value(b, dvalue);
                continue;
            }
            jb_lit(b, "[");
            for (int64_t c = 0; c < dvalue->count; c++) {
                const jvalue *component = jat(dvalue, c);
                if (c) jb_lit(b, ", ");
                jb_lit(b, "{");
                int first = 1;
                for (int64_t k = 0; k < component->count; k++) {
                    const jvalue *ckey = component->u.members + 2 * k;
                    if (key_is(ckey, "_domain_operation_record") || key_is(ckey, "_magnetic_domain_operation")) continue;
                    if (!first) jb_lit(b, ", ");
                    first = 0;
                    jb_string(b, ckey->u.s, (size_t)ckey->count);
                    jb_lit(b, ": ");
                    jb_value(b, ckey + 1);
                }
                if (displays && c < n) {
                    if (!first) jb_lit(b, ", ");
                    jb_lit(b, "\"domain_display\": ");
                    jb_string(b, displays[c].data, displays[c].used);
                }
                jb_lit(b, "}");
            }
            jb_lit(b, "]");
        }
        jb_lit(b, "}");
    }
    jb_lit(b, "}");
}

/* selection._complete_coupled_domain_displays for one component */
static int coupled_domain_display(jbuf *b, int sg, const sel_context *c, const jvalue *slot, const jvalue *component)
{
    jbuf label = {0};
    jb_reset(&label);
    str_or_empty(&label, member(component, "opd"));
    const jvalue *row = row_for_label(slot, label.data ? label.data : "", label.used);
    if (!row) {
        int status = md_fail(MDE_VALUE, "selected Source OPD row not found: %.*s", (int)label.used, label.data ? label.data : "");
        jb_free(&label);
        return status;
    }
    const jvalue *iso = member(row, "isotropy");
    int64_t domain = 0;
    if (py_int_of(member(component, "domain"), &domain)) {
        jb_free(&label);
        return md_fail(MDE_TYPE, "int() argument must be a string, a bytes-like object or a real number");
    }
    jbuf basis = {0}, origin = {0};
    jb_reset(&basis);
    jb_reset(&origin);
    str_or_empty(&basis, member(iso, "basis_text"));
    str_or_empty(&origin, member(iso, "origin"));
    int status = 0;
    if (domain > 1) {
        const jvalue *raw = member(component, "_domain_operation_record");
        int ok = jis(raw, J_ARRAY) && raw->count == 5;
        int64_t record[5];
        for (int k = 0; ok && k < 5; k++) {
            ok = jis(jat(raw, k), J_INT);
            if (ok) record[k] = jat(raw, k)->u.i;
        }
        if (!ok) status = md_fail(MDE_VALUE, "selected domain operation record is incomplete");
        if (!status && jtruthy(member(component, "_magnetic_domain_operation"))) record[4] = section_int(&T.mag2nonmag, record[4] - 1);
        int64_t sb[9], so[4];
        const jvalue *sbv = member(iso, "source_basis_values"), *sov = member(iso, "source_origin_values");
        for (int k = 0; !status && k < 9; k++)
            if (py_int_of(jat(sbv, k), sb + k)) status = md_fail(MDE_TYPE, "source basis values are not integers");
        for (int k = 0; !status && k < 4; k++)
            if (py_int_of(jat(sov, k), so + k)) status = md_fail(MDE_TYPE, "source origin values are not integers");
        int64_t point[9];
        if (!status && geo_pml_point_matrix(sg, (int)record[4], point))
            status = md_fail(MDE_KEY, "missing pml point operation SG%d/op%lld", sg, (long long)record[4]);
        if (!status) {
            rat transform[4][4], operation[4][4];
            for (int r = 0; r < 4; r++)
                for (int k = 0; k < 4; k++) transform[r][k] = operation[r][k] = rat_int(0);
            for (int r = 0; r < 3; r++)
                for (int k = 0; k < 3; k++) {
                    transform[r][k] = rat_int(sb[r * 3 + k]);
                    operation[r][k] = rat_int(point[r * 3 + k]);
                }
            for (int k = 0; k < 3; k++) {
                transform[3][k] = rat_make(so[k], so[3] ? so[3] : 1);
                operation[3][k] = rat_make(record[k], record[3]);
            }
            transform[3][3] = operation[3][3] = rat_int(1);
            rat basis_t[4][4], origin_t[4][4];
            for (int r = 0; r < 4; r++)
                for (int k = 0; k < 4; k++) {
                    rat x = rat_int(0), y = rat_int(0);
                    for (int m = 0; m < 4; m++) {
                        x = rat_add(x, rat_mul(operation[r][m], transform[m][k]));
                        y = rat_add(y, rat_mul(transform[r][m], operation[m][k]));
                    }
                    basis_t[r][k] = x;
                    origin_t[r][k] = y;
                }
            int64_t values[9], mod[4];
            for (int r = 0; r < 3; r++)
                for (int k = 0; k < 3; k++) values[r * 3 + k] = (int64_t)(basis_t[r][k].n / basis_t[r][k].d);
            i128 den = 1;
            for (int k = 0; k < 3; k++) den = den / rat_gcd(den, origin_t[3][k].d) * origin_t[3][k].d;
            i128 nums[3], g = den;
            for (int k = 0; k < 3; k++) {
                nums[k] = origin_t[3][k].n * (den / origin_t[3][k].d);
                g = rat_gcd(g, nums[k] < 0 ? -nums[k] : nums[k]);
            }
            if (g) {
                den /= g;
                for (int k = 0; k < 3; k++) nums[k] /= g;
            }
            for (int k = 0; k < 3; k++) mod[k] = (int64_t)(((nums[k] % den) + den) % den);
            mod[3] = (int64_t)den;
            const jvalue *subgroup = member(iso, "subgroup");
            int64_t ordinary = member_int(subgroup, "ordinary_number");
            if (!ordinary) ordinary = member_int(subgroup, "number");
            sel_setting displayed;
            int same = sg == ordinary;
            status = sel_change_setting(sg, (int)ordinary, values, mod, 1, c->setting_id, same, c->setting_id, &displayed);
            if (!status) {
                rat rows[3][3];
                for (int r = 0; r < 3; r++)
                    for (int k = 0; k < 3; k++) rows[r][k] = rat_make(displayed.basis[r][k], displayed.denominator ? displayed.denominator : 1);
                jbuf q = {0};
                jb_reset(&q);
                sel_put_basis_text(&q, (const rat (*)[3])rows);
                jdoc *doc = jdoc_parse(q.data, q.used);
                jb_reset(&basis);
                jb_raw(&basis, jdoc_root(doc)->u.s, (size_t)jdoc_root(doc)->count);
                jdoc_free(doc);
                jb_reset(&q);
                rat integer_rows[3][3];
                for (int r = 0; r < 3; r++)
                    for (int k = 0; k < 3; k++) integer_rows[r][k] = rat_int(displayed.basis[r][k]);
                sel_put_origin_text(&q, sg, (int)ordinary, (const rat (*)[3])integer_rows, displayed.origin);
                doc = jdoc_parse(q.data, q.used);
                jb_reset(&origin);
                jb_raw(&origin, jdoc_root(doc)->u.s, (size_t)jdoc_root(doc)->count);
                jdoc_free(doc);
                jb_free(&q);
            }
        }
    }
    if (!status) {
        jbuf display = {0};
        jb_reset(&display);
        str_or_empty(&display, member(component, "display_opd"));
        domain_display(b, iso, label.data ? label.data : "", label.used, domain, display.data ? display.data : "", display.used, &basis, &origin);
        jb_free(&display);
    }
    jb_free(&label);
    jb_free(&basis);
    jb_free(&origin);
    return status;
}

static int multi(const jvalue *state, const jvalue *opts, int sg, jbuf *out)
{
    const jvalue *selections = member(opts, "selections");
    const jvalue *irrep_slots = member(member(state, "selected"), "irrep_slots");
    int64_t n = selections->count;
    if (!jis(irrep_slots, J_ARRAY) || irrep_slots->count != n) return md_fail(MDE_VALUE, "K/IR slot count changed before coupled OPD generation");
    const jvalue *orderparam = member(opts, "orderparam");
    const char *starts[16];
    int lengths[16];
    int nlabels = coupled_labels(orderparam, (int)n, jtruthy(member(opts, "selected_opd_only")), starts, lengths);
    sel_context c;
    context_of(state, &c);
    jbuf slots = {0};
    jb_reset(&slots);
    jb_lit(&slots, "[");
    int status = 0;
    for (int64_t i = 0; i < n && !status; i++) {
        const jvalue *slot_state = jat(irrep_slots, i), *selection = jat(selections, i);
        int64_t index;
        status = truthy_int(member(selection, "irrep_index"), 1, &index);
        const jvalue *irrep = 0;
        if (!status) status = select_irrep(member(slot_state, "irreps"), member(selection, "irrep"), index, &irrep);
        if (status) break;
        jbuf kp = {0}, dp = {0};
        jb_reset(&kp);
        jb_reset(&dp);
        sel_normalize_params(&kp, member(selection, "k_params"));
        sel_normalize_params(&dp, member(selection, "display_k_params"));
        if (params_empty(&dp)) {
            jb_reset(&dp);
            jb_raw(&dp, kp.data, kp.used);
        }
        jdoc *kdoc = jdoc_parse(kp.data, kp.used), *ddoc = jdoc_parse(dp.data, dp.used);
        sel_params k_params, display;
        sel_params_of(jdoc_root(kdoc), &k_params);
        sel_params_of(jdoc_root(ddoc), &display);
        jbuf selector = {0}, rows = {0};
        jb_reset(&selector);
        jb_reset(&rows);
        jdoc *selector_doc = 0;
        if (nlabels) {
            jb_string(&selector, starts[i], (size_t)lengths[i]);
            selector_doc = jdoc_parse(selector.data, selector.used);
        }
        status = rows_for(sg, irrep, &k_params, &display, selector_doc ? jdoc_root(selector_doc) : 0, &c, &rows);
        int64_t kparam[4] = {0, 0, 0, 1};
        if (!status) {
            int found = sel_dynamic_kparam(member_int(irrep, "gid"), &k_params, kparam);
            if (found < 0) status = -found;
            else if (!found) {
                kparam[0] = kparam[1] = kparam[2] = 0;
                kparam[3] = 1;
            }
        }
        if (!status) {
            if (i) jb_lit(&slots, ", ");
            jb_lit(&slots, "{\"slot\": ");
            jb_int(&slots, i + 1);
            jb_lit(&slots, ", \"kpoint\": ");
            const jvalue *kpoint = member(slot_state, "kpoint");
            if (kpoint) jb_value(&slots, kpoint);
            else jb_lit(&slots, "null");
            jb_lit(&slots, ", \"irrep\": ");
            jb_value(&slots, irrep);
            jb_lit(&slots, ", \"k_params\": ");
            jb_raw(&slots, kp.data, kp.used);
            jb_lit(&slots, ", \"display_k_params\": ");
            jb_raw(&slots, dp.data, dp.used);
            jb_lit(&slots, ", \"source_kparam\": [");
            for (int k = 0; k < 4; k++) {
                if (k) jb_lit(&slots, ", ");
                jb_int(&slots, kparam[k]);
            }
            jb_lit(&slots, "], \"opd_rows\": ");
            if (nlabels) {
                jdoc *rdoc = jdoc_parse(rows.data, rows.used);
                const jvalue *list = jdoc_root(rdoc);
                jb_lit(&slots, "[");
                int first = 1;
                for (int64_t r = 0; list && r < list->count; r++) {
                    const jvalue *row = jat(list, r);
                    if (!stripped_is(member(member(row, "isotropy"), "opd_label"), starts[i], (size_t)lengths[i]) &&
                        !stripped_is(member(member(row, "direction"), "label"), starts[i], (size_t)lengths[i]))
                        continue;
                    if (!first) jb_lit(&slots, ", ");
                    first = 0;
                    jb_value(&slots, row);
                }
                jb_lit(&slots, "]");
                jdoc_free(rdoc);
            } else {
                jb_raw(&slots, rows.data, rows.used);
            }
            jb_lit(&slots, "}");
        }
        jdoc_free(selector_doc);
        jdoc_free(kdoc);
        jdoc_free(ddoc);
        jb_free(&kp);
        jb_free(&dp);
        jb_free(&selector);
        jb_free(&rows);
    }
    jb_lit(&slots, "]");
    jdoc *sdoc = status ? 0 : jdoc_parse(slots.data, slots.used);
    jbuf request = {0}, coupled = {0};
    jb_reset(&request);
    jb_reset(&coupled);
    if (!status) {
        jb_lit(&request, "{\"parent_sg\": ");
        jb_int(&request, sg);
        jb_lit(&request, ", \"slots\": ");
        jb_raw(&request, slots.data, slots.used);
        jb_lit(&request, ", \"parent_setting_id\": ");
        jb_int(&request, c.setting_id);
        jb_lit(&request, ", \"parent_cell\": ");
        if (c.has_cell) {
            jb_lit(&request, "[");
            for (int k = 0; k < 6; k++) {
                if (k) jb_lit(&request, ", ");
                jb_float(&request, c.cell[k]);
            }
            jb_lit(&request, "]");
        } else {
            jb_lit(&request, "null");
        }
        jb_lit(&request, ", \"row_limit\": ");
        const jvalue *limit = member(opts, "opd_row_limit");
        if (limit) jb_value(&request, limit);
        else jb_lit(&request, "null");
        jb_lit(&request, "}");
        jdoc *qdoc = jdoc_parse(request.data, request.used);
        status = sel_coupled_rows(jdoc_root(qdoc), &coupled);
        jdoc_free(qdoc);
    }
    jdoc *cdoc = status ? 0 : jdoc_parse(coupled.data, coupled.used);
    const jvalue *rows = cdoc ? jdoc_root(cdoc) : 0, *selected = 0;
    int64_t orderparam_index = jis(member(opts, "orderparam_index"), J_INT) ? member(opts, "orderparam_index")->u.i : 1;
    if (!status) status = select_orderparam(rows, orderparam, orderparam_index, &selected);
    jbuf *displays = 0;
    int ndisplays = 0;
    if (!status && selected) {
        const jvalue *components = member(member(selected, "direction"), "components");
        int64_t m = jis(components, J_ARRAY) ? components->count : 0;
        if (m != n) status = md_fail(MDE_VALUE, "selected OPD components do not match primary slots");
        displays = calloc((size_t)n + 1, sizeof(jbuf));
        if (!status && !displays) status = md_fail(MDE_MEMORY, "out of memory");
        for (int64_t i = 0; !status && i < n; i++) {
            jb_reset(&displays[i]);
            status = coupled_domain_display(&displays[i], sg, &c, jat(jdoc_root(sdoc), i), jat(components, i));
            ndisplays = (int)i + 1;
        }
    }
    jbuf rows_json = {0}, selected_json = {0}, factors = {0}, sel_json = {0}, irreps = {0}, first_irrep = {0};
    jb_reset(&rows_json);
    jb_reset(&selected_json);
    jb_reset(&factors);
    jb_reset(&sel_json);
    jb_reset(&irreps);
    jb_reset(&first_irrep);
    if (!status) {
        jb_lit(&rows_json, "[");
        for (int64_t i = 0; i < rows->count; i++) {
            if (i) jb_lit(&rows_json, ", ");
            put_coupled_row(&rows_json, jat(rows, i), jat(rows, i) == selected ? displays : 0, ndisplays);
        }
        jb_lit(&rows_json, "]");
        if (selected) put_coupled_row(&selected_json, selected, displays, ndisplays);
        else jb_lit(&selected_json, "null");
        jdoc *seldoc = jdoc_parse(selected_json.data, selected_json.used);
        status = sel_opd_factors(&factors, sg, jdoc_root(sdoc), jdoc_root(seldoc));
        jdoc_free(seldoc);
    }
    if (!status) {
        const jvalue *list = jdoc_root(sdoc);
        jb_value(&first_irrep, member(jat(list, 0), "irrep"));
        jb_lit(&irreps, "[");
        for (int64_t i = 0; i < list->count; i++) {
            if (i) jb_lit(&irreps, ", ");
            jb_value(&irreps, member(jat(list, i), "irrep"));
        }
        jb_lit(&irreps, "]");
        kv items[6] = {{"irrep", first_irrep.data, first_irrep.used}, {"irreps_selected", irreps.data, irreps.used},
                       {"selection_slots", slots.data, slots.used},    {"opd_factors", factors.data, factors.used},
                       {"opd_rows", rows_json.data, rows_json.used},   {"orderparam", selected_json.data, selected_json.used}};
        sel_merge(&sel_json, member(state, "selected"), items, 6);
        kv top[2] = {{"schema", "\"APOSTRUCT.opd_direction.multi.v1\"", strlen("\"APOSTRUCT.opd_direction.multi.v1\"")}, {"selected", sel_json.data, sel_json.used}};
        sel_merge(out, state, top, 2);
    }
    for (int i = 0; i < ndisplays; i++) jb_free(&displays[i]);
    free(displays);
    jdoc_free(cdoc);
    jdoc_free(sdoc);
    jb_free(&slots);
    jb_free(&request);
    jb_free(&coupled);
    jb_free(&rows_json);
    jb_free(&selected_json);
    jb_free(&factors);
    jb_free(&sel_json);
    jb_free(&irreps);
    jb_free(&first_irrep);
    return status;
}

/* select_order_parameter_direction(state, **opts) */
int sel_order_parameter_direction(const jvalue *state, const jvalue *opts, jbuf *out)
{
    MD_TRY(load());
    int64_t sg = member_int(member(member(state, "input"), "parent"), "number");
    const jvalue *selections = member(opts, "selections");
    if (jtruthy(selections) && jis(selections, J_ARRAY) && selections->count > 1) return multi(state, opts, (int)sg, out);
    return single(state, opts, (int)sg, out);
}
