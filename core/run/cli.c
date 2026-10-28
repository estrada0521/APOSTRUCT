/* The CLI entry from a normalized case to the printed result: case_input's
 * pipeline request, cli.execute_case, cli_output.compact_result and
 * json.dumps(sort_keys=True, ensure_ascii=False, indent=N). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "md.h"
#include "par.h"
#include "sel.h"

const char *cif_summary_json(const char *data, int64_t length, int from_file);
const char *k_catalog_json(int sg);
int k_source_kparam_of_visible(int sg, int kslot, const int64_t *visible_words, int64_t *out);
int st_state(const jvalue *request, jbuf *b);
int sel_k_points(const jvalue *state, const jvalue *opts, jbuf *out);
int64_t mode_details_vectors(const double **data);
int md_render_text(const jvalue *state, const double *vectors, int64_t nvectors, jbuf *b);

static jbuf out;
static int run_status;

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

static int ascii_alpha(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }

static int ieq(const char *a, size_t an, const char *b, size_t bn)
{
    if (an != bn) return 0;
    for (size_t i = 0; i < an; i++) {
        char x = a[i], y = b[i];
        if (x >= 'A' && x <= 'Z') x = (char)(x - 'A' + 'a');
        if (y >= 'A' && y <= 'Z') y = (char)(y - 'A' + 'a');
        if (x != y) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------ output */

typedef struct {
    const double *vectors;
    int64_t nvectors;
} expand_t;

static int is_f64(const jvalue *v)
{
    return jis(v, J_OBJECT) && v->count == 1 && jstr_is(v->u.members, "$f64") && jis(v->u.members + 1, J_ARRAY) &&
           (v->u.members[1].count == 2 || v->u.members[1].count == 3);
}

/* A compact copy of v; private keys dropped when public. */
static void put_copy(jbuf *b, const jvalue *v, int public_only, const expand_t *x)
{
    if (!v) {
        jb_lit(b, "null");
        return;
    }
    if (is_f64(v)) {
        jb_value(b, v);
        return;
    }
    if (jis(v, J_ARRAY)) {
        jb_lit(b, "[");
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_lit(b, ", ");
            put_copy(b, v->u.items + i, public_only, x);
        }
        jb_lit(b, "]");
        return;
    }
    if (jis(v, J_OBJECT)) {
        jb_lit(b, "{");
        int first = 1;
        for (int64_t i = 0; i < v->count; i++) {
            const jvalue *key = v->u.members + 2 * i;
            if (public_only && key->count && key->u.s[0] == '_') continue;
            if (!first) jb_lit(b, ", ");
            first = 0;
            jb_value(b, key);
            jb_lit(b, ": ");
            put_copy(b, key + 1, public_only, x);
        }
        jb_lit(b, "}");
        return;
    }
    jb_value(b, v);
}

/* json.dumps with ensure_ascii=False */
static void put_text(jbuf *b, const char *s, int64_t n)
{
    jb_raw(b, "\"", 1);
    int64_t start = 0;
    for (int64_t i = 0; i < n; i++) {
        unsigned char c = (unsigned char)s[i];
        const char *escape = 0;
        char unicode[8];
        if (c == '"') escape = "\\\"";
        else if (c == '\\') escape = "\\\\";
        else if (c == '\n') escape = "\\n";
        else if (c == '\r') escape = "\\r";
        else if (c == '\t') escape = "\\t";
        else if (c == '\b') escape = "\\b";
        else if (c == '\f') escape = "\\f";
        else if (c < 0x20) {
            snprintf(unicode, sizeof unicode, "\\u%04x", c);
            escape = unicode;
        }
        if (!escape) continue;
        jb_raw(b, s + start, (size_t)(i - start));
        jb_lit(b, escape);
        start = i + 1;
    }
    jb_raw(b, s + start, (size_t)(n - start));
    jb_raw(b, "\"", 1);
}

static int key_order(const void *a, const void *b)
{
    const jvalue *x = *(const jvalue *const *)a, *y = *(const jvalue *const *)b;
    int64_t n = x->count < y->count ? x->count : y->count;
    int c = memcmp(x->u.s, y->u.s, (size_t)n);
    if (c) return c;
    return x->count < y->count ? -1 : x->count > y->count;
}

static void newline(jbuf *b, int indent, int depth)
{
    jb_raw(b, "\n", 1);
    static const char spaces[] = "                                                                ";
    int64_t remaining = (int64_t)indent * depth;
    while (remaining > 0) {
        size_t n = remaining < (int64_t)(sizeof spaces - 1) ? (size_t)remaining : sizeof spaces - 1;
        jb_raw(b, spaces, n);
        remaining -= (int64_t)n;
    }
}

static int put_pretty(jbuf *b, const jvalue *v, int indent, int depth, const expand_t *x)
{
    if (is_f64(v)) {
        const jvalue *pair = v->u.members + 1;
        int64_t offset = jat(pair, 0)->u.i, atoms = jat(pair, 1)->u.i;
        jb_lit(b, "[");
        if (pair->count == 3) {
            for (int64_t a = 0; a < atoms; a++) {
                if (a) jb_lit(b, ",");
                newline(b, indent, depth + 1);
                double value = offset + a < x->nvectors ? x->vectors[offset + a] : 0.0;
                if (value != value || value - value != 0.0)
                    return md_fail(MDE_VALUE, "Out of range float values are not JSON compliant");
                jb_float(b, value);
            }
            if (atoms) newline(b, indent, depth);
            jb_lit(b, "]");
            return 0;
        }
        for (int64_t a = 0; a < atoms; a++) {
            if (a) jb_lit(b, ",");
            newline(b, indent, depth + 1);
            jb_lit(b, "[");
            for (int k = 0; k < 3; k++) {
                if (k) jb_lit(b, ",");
                newline(b, indent, depth + 2);
                int64_t at = offset + 3 * a + k;
                double value = at < x->nvectors ? x->vectors[at] : 0.0;
                if (value != value || value - value != 0.0)
                    return md_fail(MDE_VALUE, "Out of range float values are not JSON compliant");
                jb_float(b, value);
            }
            newline(b, indent, depth + 1);
            jb_lit(b, "]");
        }
        if (atoms) newline(b, indent, depth);
        jb_lit(b, "]");
        return 0;
    }
    switch (v->kind) {
    case J_STRING: put_text(b, v->u.s, v->count); return 0;
    case J_FLOAT:
        if (v->u.f != v->u.f || v->u.f - v->u.f != 0.0) return md_fail(MDE_VALUE, "Out of range float values are not JSON compliant");
        jb_float(b, v->u.f);
        return 0;
    case J_ARRAY:
        if (!v->count) {
            jb_lit(b, "[]");
            return 0;
        }
        jb_lit(b, "[");
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_lit(b, ",");
            newline(b, indent, depth + 1);
            MD_TRY(put_pretty(b, v->u.items + i, indent, depth + 1, x));
        }
        newline(b, indent, depth);
        jb_lit(b, "]");
        return 0;
    case J_OBJECT: {
        if (!v->count) {
            jb_lit(b, "{}");
            return 0;
        }
        const jvalue **keys = malloc(sizeof *keys * (size_t)v->count);
        if (!keys) return md_fail(MDE_MEMORY, "out of memory");
        for (int64_t i = 0; i < v->count; i++) keys[i] = v->u.members + 2 * i;
        qsort(keys, (size_t)v->count, sizeof *keys, key_order);
        jb_lit(b, "{");
        int status = 0;
        for (int64_t i = 0; i < v->count && !status; i++) {
            if (i) jb_lit(b, ",");
            newline(b, indent, depth + 1);
            put_text(b, keys[i]->u.s, keys[i]->count);
            jb_lit(b, ": ");
            status = put_pretty(b, keys[i] + 1, indent, depth + 1, x);
        }
        free(keys);
        if (status) return status;
        newline(b, indent, depth);
        jb_lit(b, "}");
        return 0;
    }
    default: jb_value(b, v); return 0;
    }
}

/* Parse metadata only; expand vector placeholders in the final formatter. */
static int print_json(jbuf *compact, int indent, jbuf *b, const expand_t *x)
{
    if (compact->failed) return md_fail(MDE_MEMORY, "out of memory");
    jdoc *doc = jdoc_parse(compact->data, compact->used);
    jb_free(compact);
    if (!doc) return md_fail(MDE_MEMORY, "out of memory");
    jb_reset(b);
    int status = put_pretty(b, jdoc_root(doc), indent < 0 ? 0 : indent, 0, x);
    jb_raw(b, "\n", 1);
    jdoc_free(doc);
    return status ? status : b->failed ? md_fail(MDE_MEMORY, "out of memory") : 0;
}

/* ------------------------------------------------------------ Python texts */

static void put_repr(jbuf *b, const char *s, size_t n) { py_repr_text(b, s, (int64_t)n); }

typedef struct {
    const char *s;
    size_t n;
} text;

static void put_repr_list(jbuf *b, const text *items, int count)
{
    jb_lit(b, "[");
    for (int i = 0; i < count; i++) {
        if (i) jb_lit(b, ", ");
        put_repr(b, items[i].s, items[i].n);
    }
    jb_lit(b, "]");
}

static int fail_text(jbuf *message)
{
    int status = md_fail(MDE_VALUE, "%.*s", (int)message->used, message->data ? message->data : "");
    jb_free(message);
    return status;
}

static text text_of(const jvalue *v) { return (text){v->u.s, (size_t)v->count}; }

static int text_cmp(const void *a, const void *b)
{
    const text *x = a, *y = b;
    size_t n = x->n < y->n ? x->n : y->n;
    int c = memcmp(x->s, y->s, n);
    return c ? c : x->n < y->n ? -1 : x->n > y->n;
}

static int text_in(const text *list, int count, text t)
{
    for (int i = 0; i < count; i++)
        if (!text_cmp(list + i, &t)) return 1;
    return 0;
}

/* ------------------------------------------------------------ request */

static void put_or_null(jbuf *b, const jvalue *v)
{
    if (v) jb_value(b, v);
    else jb_lit(b, "null");
}

static void put_str(jbuf *b, const jvalue *v)
{
    jbuf t = {0};
    jb_reset(&t);
    py_str_of(&t, v);
    jb_string(b, t.data ? t.data : "", t.used);
    jb_free(&t);
}

/* the parameter names (a, b, g) of a k point */
static int parameter_names(const jvalue *kpoint, text *names)
{
    const jvalue *expression = member(kpoint, "display_kvector");
    const char *s = jis(expression, J_STRING) ? expression->u.s : "";
    int64_t n = jis(expression, J_STRING) ? expression->count : 0;
    static const char *const order[] = {"a", "b", "g"};
    int count = 0;
    for (int k = 0; k < 3; k++)
        for (int64_t i = 0; i < n; i++)
            if (s[i] == order[k][0] && !(i > 0 && ascii_alpha(s[i - 1])) && !(i + 1 < n && ascii_alpha(s[i + 1]))) {
                names[count++] = (text){order[k], 1};
                break;
            }
    return count;
}

/* An object of k parameters for one k slot (case_input._slot_params). */
static int slot_params(const jvalue *slot, int ordinal, const jvalue *catalog, jbuf *b)
{
    const jvalue *label = member(slot, "label"), *params = member(slot, "params"), *values = member(slot, "values");
    const jvalue *match = 0;
    int matches = 0;
    for (int64_t i = 0; jis(catalog, J_ARRAY) && i < catalog->count; i++) {
        const jvalue *item = jat(catalog, i), *item_label = member(item, "label");
        const char *s = jis(item_label, J_STRING) ? item_label->u.s : "";
        size_t n = jis(item_label, J_STRING) ? (size_t)item_label->count : 0;
        if (ieq(s, n, label->u.s, (size_t)label->count)) {
            matches++;
            match = item;
        }
    }
    jbuf message = {0};
    jb_reset(&message);
    if (matches != 1) {
        char head[32];
        snprintf(head, sizeof head, "k[%d] label ", ordinal);
        jb_lit(&message, head);
        put_repr(&message, label->u.s, (size_t)label->count);
        jb_lit(&message, " is not uniquely available");
        return fail_text(&message);
    }
    text names[3];
    int nnames = parameter_names(match, names);
    int64_t nparams = jis(params, J_ARRAY) ? params->count : 0;
    char head[32];
    snprintf(head, sizeof head, "k[%d] ", ordinal);
    if (nparams) {
        int same = 1;
        text given[16];
        int ngiven = 0;
        for (int64_t i = 0; i < nparams; i++) {
            text t = text_of(jat(jat(params, i), 0));
            if (!text_in(given, ngiven, t) && ngiven < 16) given[ngiven++] = t;
        }
        if (ngiven != nnames) same = 0;
        for (int i = 0; same && i < nnames; i++) same = text_in(given, ngiven, names[i]);
        if (!same) {
            jb_lit(&message, head);
            jb_raw(&message, label->u.s, (size_t)label->count);
            jb_lit(&message, " parameters must be exactly ");
            put_repr_list(&message, names, nnames);
            return fail_text(&message);
        }
        jb_free(&message);
        jb_lit(b, "{");
        int first = 1;
        for (int64_t i = 0; i < nparams; i++) {
            const jvalue *pair = jat(params, i);
            int later = 0;
            for (int64_t k = i + 1; k < nparams && !later; k++) later = !text_cmp(&(text){jat(pair, 0)->u.s, (size_t)jat(pair, 0)->count},
                                                                             &(text){jat(jat(params, k), 0)->u.s, (size_t)jat(jat(params, k), 0)->count});
            if (later) continue;
            if (!first) jb_lit(b, ", ");
            first = 0;
            jb_value(b, jat(pair, 0));
            jb_lit(b, ": ");
            jb_value(b, jat(pair, 1));
        }
        jb_lit(b, "}");
        return 0;
    }
    int64_t nvalues = jis(values, J_ARRAY) ? values->count : 0;
    if (nvalues != nnames) {
        jb_lit(&message, head);
        jb_raw(&message, label->u.s, (size_t)label->count);
        char count[32];
        snprintf(count, sizeof count, " requires %d parameters ", nnames);
        jb_lit(&message, count);
        put_repr_list(&message, names, nnames);
        return fail_text(&message);
    }
    jb_free(&message);
    jb_lit(b, "{");
    for (int i = 0; i < nnames; i++) {
        if (i) jb_lit(b, ", ");
        jb_string(b, names[i].s, names[i].n);
        jb_lit(b, ": ");
        jb_value(b, jat(values, i));
    }
    jb_lit(b, "}");
    return 0;
}

static void put_index_list(jbuf *b, const int *items, int count)
{
    jb_lit(b, "[");
    for (int i = 0; i < count; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, items[i]);
    }
    jb_lit(b, "]");
}

static text site_name(const jvalue *atom, const char *first, const char *second, int index, char *scratch)
{
    const jvalue *a = member(atom, first), *b = member(atom, second);
    if (jtruthy(a) && jis(a, J_STRING)) return text_of(a);
    if (jtruthy(b) && jis(b, J_STRING)) return text_of(b);
    snprintf(scratch, 16, "%d", index);
    return (text){scratch, strlen(scratch)};
}

/* the selection options, as
 * object members without braces. */
static int pipeline_options(const jvalue *request, const jvalue *cif_info, int require_irreps, int allow_empty, jbuf *b)
{
    const jvalue *atoms = member(cif_info, "atom_sites"), *sites = member(request, "sites"), *slots = member(request, "k");
    int natoms = jis(atoms, J_ARRAY) ? (int)atoms->count : 0, nsites = jis(sites, J_ARRAY) ? (int)sites->count : 0;
    int nslots = jis(slots, J_ARRAY) ? (int)slots->count : 0, strain = jtruthy(member(request, "strain"));
    text *types = calloc((size_t)natoms + 1, sizeof(text)), *labels = calloc((size_t)natoms + 1, sizeof(text));
    text *selectors = calloc((size_t)(2 * natoms) + 1, sizeof(text)), *unknown = calloc((size_t)nsites + 1, sizeof(text));
    char (*scratch)[2][16] = calloc((size_t)natoms + 1, sizeof *scratch);
    int *displacive = calloc((size_t)natoms + 1, sizeof(int)), *magnetic = calloc((size_t)natoms + 1, sizeof(int));
    int *occupancy = calloc((size_t)natoms + 1, sizeof(int));
    int status = 0;
    if (!types || !labels || !selectors || !unknown || !scratch || !displacive || !magnetic || !occupancy) {
        status = md_fail(MDE_MEMORY, "out of memory");
        goto done;
    }
    int nselectors = 0, ndisp = 0, nmag = 0, nocc = 0, nunknown = 0;
    for (int i = 0; i < natoms; i++) {
        types[i] = site_name(jat(atoms, i), "type", "label", i + 1, scratch[i][0]);
        labels[i] = site_name(jat(atoms, i), "label", "type", i + 1, scratch[i][1]);
    }
    for (int i = 0; i < natoms; i++)
        if (!text_in(selectors, nselectors, types[i])) selectors[nselectors++] = types[i];
    for (int i = 0; i < natoms; i++)
        if (!text_in(selectors, nselectors, labels[i])) selectors[nselectors++] = labels[i];
    for (int s = 0; s < nsites; s++) {
        const jvalue *site = jat(sites, s);
        text name = text_of(jat(site, 0));
        if (!jtruthy(jat(site, 1)) && !jtruthy(jat(site, 2)) && !jtruthy(jat(site, 3))) continue;
        if (!text_in(selectors, nselectors, name) && !text_in(unknown, nunknown, name)) unknown[nunknown++] = name;
    }
    if (nunknown) {
        qsort(unknown, (size_t)nunknown, sizeof(text), text_cmp);
        jbuf message = {0};
        jb_reset(&message);
        jb_lit(&message, "unknown CIF site selectors: ");
        put_repr_list(&message, unknown, nunknown);
        jb_lit(&message, "; available selectors: ");
        put_repr_list(&message, selectors, nselectors);
        status = fail_text(&message);
        goto done;
    }
    for (int i = 0; i < natoms; i++) {
        int d = 0, m = 0, o = 0;
        for (int s = 0; s < nsites; s++) {
            const jvalue *site = jat(sites, s);
            text name = text_of(jat(site, 0));
            if (!jtruthy(jat(site, 1)) && !jtruthy(jat(site, 2)) && !jtruthy(jat(site, 3))) continue;
            if (!text_cmp(&name, types + i) || !text_cmp(&name, labels + i)) {
                int last = 1;
                for (int t = s + 1; t < nsites && last; t++) {
                    const jvalue *other = jat(sites, t);
                    text again = text_of(jat(other, 0));
                    if ((jtruthy(jat(other, 1)) || jtruthy(jat(other, 2)) || jtruthy(jat(other, 3))) && !text_cmp(&again, &name)) last = 0;
                }
                if (!last) continue;
                d |= jtruthy(jat(site, 1));
                m |= jtruthy(jat(site, 2));
                o |= jtruthy(jat(site, 3));
            }
        }
        if (d) displacive[ndisp++] = i + 1;
        if (m) magnetic[nmag++] = i + 1;
        if (o) occupancy[nocc++] = i + 1;
    }
    if (!allow_empty && !strain && !ndisp && !nmag && !nocc) {
        status = md_fail(MDE_VALUE, "case must include strain, displacive, magnetic, or occupancy distortion");
        goto done;
    }
    int pure_strain_irreps = !require_irreps && !allow_empty && strain && !ndisp && !nmag && !nocc;
    if (pure_strain_irreps)
        for (int k = 0; k < nslots; k++) {
            const jvalue *label = member(jat(slots, k), "label");
            if (!ieq(label->u.s, (size_t)label->count, "gm", 2)) {
                status = md_fail(MDE_VALUE, "strain-only irreps require k label GM or no k");
                goto done;
            }
        }
    if (!nslots && !allow_empty && !pure_strain_irreps) {
        status = md_fail(MDE_VALUE, "case requires one to four ordered k slots");
        goto done;
    }
    int64_t sg = 0;
    py_int_of(member(member(cif_info, "parent"), "number"), &sg);
    const char *catalog_text = k_catalog_json((int)sg);
    if (!catalog_text) {
        status = md_fail(MDE_VALUE, "k-point catalog Source tables are inconsistent");
        goto done;
    }
    jdoc *catalog = jdoc_parse(catalog_text, strlen(catalog_text));
    if (!catalog) {
        status = md_fail(MDE_MEMORY, "out of memory");
        goto done;
    }
    jbuf selections = {0}, first = {0};
    jb_reset(&selections);
    jb_reset(&first);
    jb_lit(&selections, "[");
    for (int k = 0; k < nslots && !status; k++) {
        const jvalue *slot = jat(slots, k), *irrep = member(slot, "ir");
        if (require_irreps && !jtruthy(irrep)) {
            status = md_fail(MDE_VALUE, "k[%d] requires an irrep", k + 1);
            break;
        }
        jbuf params = {0}, one = {0};
        jb_reset(&params);
        jb_reset(&one);
        status = slot_params(slot, k + 1, member(jdoc_root(catalog), "kpoints"), &params);
        if (!status) {
            jb_lit(&one, "{\"k_label\": ");
            jb_value(&one, member(slot, "label"));
            jb_lit(&one, ", \"k_index\": 1, \"k_params\": ");
            jb_raw(&one, params.data, params.used);
            jb_lit(&one, ", \"display_k_params\": ");
            jb_raw(&one, params.data, params.used);
            jb_lit(&one, ", \"irrep\": ");
            put_or_null(&one, irrep);
            jb_lit(&one, ", \"irrep_index\": 1}");
            if (k) jb_lit(&selections, ", ");
            jb_raw(&selections, one.data, one.used);
            if (!k) jb_raw(&first, one.data, one.used);
        }
        jb_free(&params);
        jb_free(&one);
    }
    jb_lit(&selections, "]");
    jdoc_free(catalog);
    if (!status) {
        jdoc *fdoc = nslots ? jdoc_parse(first.data, first.used) : 0;
        const jvalue *f = fdoc ? jdoc_root(fdoc) : 0;
        jb_lit(b, "\"mode_set\": ");
        jb_lit(b, nmag ? "\"all\"" : strain ? "\"strain_displacive\"" : "\"displacive\"");
        jb_lit(b, ", \"distortion_selection\": {\"strain\": ");
        jb_lit(b, strain ? "true" : "false");
        jb_lit(b, ", \"displacive_sites\": ");
        put_index_list(b, displacive, ndisp);
        jb_lit(b, ", \"magnetic_sites\": ");
        put_index_list(b, magnetic, nmag);
        if (nocc) {
            jb_lit(b, ", \"occupancy_sites\": ");
            put_index_list(b, occupancy, nocc);
        }
        jb_lit(b, "}, \"k_label\": ");
        if (f) jb_value(b, jget(f, "k_label"));
        else jb_lit(b, "\"GM\"");
        jb_lit(b, ", \"k_index\": 1, \"k_params\": ");
        if (f) jb_value(b, jget(f, "k_params"));
        else jb_lit(b, "{}");
        jb_lit(b, ", \"display_k_params\": ");
        if (f) jb_value(b, jget(f, "display_k_params"));
        else jb_lit(b, "{}");
        jb_lit(b, ", \"irrep\": ");
        if (f) jb_value(b, jget(f, "irrep"));
        else jb_lit(b, "null");
        jb_lit(b, ", \"irrep_index\": 1, \"orderparam\": ");
        const jvalue *opd = member(request, "opd");
        if (opd) jb_value(b, opd);
        else jb_lit(b, "null");
        jb_lit(b, ", \"orderparam_index\": 1, \"selections\": ");
        if (nslots) jb_raw(b, selections.data, selections.used);
        else jb_lit(b, "null");
        jdoc_free(fdoc);
    }
    jb_free(&selections);
    jb_free(&first);
done:
    free(types);
    free(labels);
    free(selectors);
    free(unknown);
    free(scratch);
    free(displacive);
    free(magnetic);
    free(occupancy);
    return status;
}

/* ------------------------------------------------------------ compact */

static void put_fields(jbuf *b, const jvalue *row, const char *const *fields, int nfields, int *first)
{
    for (int i = 0; i < nfields; i++) {
        const jvalue *v = member(row, fields[i]);
        if (!v) continue;
        if (!*first) jb_lit(b, ", ");
        *first = 0;
        jb_cstring(b, fields[i]);
        jb_lit(b, ": ");
        jb_value(b, v);
    }
}

static void put_key(jbuf *b, int *first, const char *key)
{
    if (!*first) jb_lit(b, ", ");
    *first = 0;
    jb_cstring(b, key);
    jb_lit(b, ": ");
}

static void kpoint_result(jbuf *b, const jvalue *row)
{
    static const char *const fields[] = {"slot",      "label",     "display",       "kid",       "kslot",  "dimension", "parameters",
                                         "display_parameters", "miller_love_parameters", "sg_specific", "star_size",
                                         "n_irreps",  "star"};
    jb_lit(b, "{");
    int first = 1;
    put_fields(b, row, fields, (int)(sizeof fields / sizeof fields[0]), &first);
    const jvalue *display = member(row, "display_kvector");
    if (display && !jis(display, J_NULL)) {
        put_key(b, &first, "kvector");
        jb_value(b, display);
    }
    put_key(b, &first, "parameter_names");
    text names[3];
    int n = parameter_names(row, names);
    jb_lit(b, "[");
    for (int i = 0; i < n; i++) {
        if (i) jb_lit(b, ", ");
        jb_string(b, names[i].s, names[i].n);
    }
    jb_lit(b, "]}");
}

static void irrep_result(jbuf *b, const jvalue *row)
{
    static const char *const fields[] = {"symbol",      "gid",      "k_label",  "kov_label",           "type",
                                         "little_dim", "full_dim", "displacive_mode_total", "occupancy_mode_total",
                                         "magnetic_mode_total"};
    jb_lit(b, "{");
    int first = 1;
    put_fields(b, row, fields, (int)(sizeof fields / sizeof fields[0]), &first);
    put_key(b, &first, "strain_visible");
    jb_lit(b, jtruthy(member(row, "strain_visible")) ? "true}" : "false}");
}

static void subgroup_result(jbuf *b, const jvalue *subgroup)
{
    const jvalue *ordinary = member(subgroup, "ordinary_number");
    int magnetic = ordinary && !jis(ordinary, J_NULL);
    if (!magnetic) ordinary = member(subgroup, "number");
    jb_lit(b, "{\"display_label\": ");
    const jvalue *label = member(subgroup, "display_label");
    if (jtruthy(label)) put_str(b, label);
    else if (jtruthy(ordinary)) put_str(b, ordinary);
    else jb_lit(b, "\"\"");
    jb_lit(b, ", \"ordinary_number\": ");
    put_or_null(b, ordinary);
    jb_lit(b, ", \"symbol\": ");
    const jvalue *symbol = member(subgroup, "symbol");
    put_or_null(b, symbol);
    jb_lit(b, magnetic ? ", \"magnetic\": true" : ", \"magnetic\": false");
    if (magnetic) {
        jb_lit(b, ", \"magnetic_id\": ");
        const jvalue *number = member(subgroup, "number");
        put_or_null(b, number);
    }
    jb_lit(b, "}");
}

/* one OPD row: the invariants of its isotropy subgroup apart from the choices (direction, frame, arms) that present it */
static void opd_result(jbuf *b, const jvalue *row)
{
    const jvalue *direction = member(row, "direction"), *iso = member(row, "isotropy");
    jb_lit(b, "{\"label\": ");
    put_or_null(b, member(direction, "label"));
    jb_lit(b, ", \"subgroup\": ");
    subgroup_result(b, member(iso, "subgroup"));
    jb_lit(b, ", \"invariants\": {\"parameter_count\": ");
    put_or_null(b, member(iso, "free"));
    jb_lit(b, ", \"index\": ");
    put_or_null(b, member(iso, "i"));
    jb_lit(b, ", \"cell_index\": ");
    put_or_null(b, member(iso, "s"));
    const jvalue *active = member(iso, "k_active_vectors");
    jb_lit(b, ", \"k_active_count\": ");
    if (jis(active, J_ARRAY)) jb_int(b, active->count);
    else jb_lit(b, "null");
    jb_lit(b, ", \"ferroic_properties\": ");
    const jvalue *ferroic = member(iso, "ferroic_properties");
    if (jtruthy(ferroic)) jb_value(b, ferroic);
    else jb_lit(b, "[]");
    jb_lit(b, "}, \"setting\": {\"opd\": ");
    put_or_null(b, member(direction, "opd"));
    jb_lit(b, ", \"basis\": ");
    put_or_null(b, member(iso, "basis_text"));
    jb_lit(b, ", \"origin\": ");
    put_or_null(b, member(iso, "origin"));
    jb_lit(b, ", \"k_active\": ");
    put_or_null(b, jtruthy(active) ? active : member(iso, "k_active"));
    const jvalue *components = member(direction, "components");
    if (components && !jis(components, J_NULL)) {
        jb_lit(b, ", \"components\": [");
        for (int64_t c = 0; jis(components, J_ARRAY) && c < components->count; c++) {
            const jvalue *component = jat(components, c);
            jb_lit(b, c ? ", {" : "{");
            int first = 1;
            for (int64_t m = 0; jis(component, J_OBJECT) && m < component->count; m++) {
                const jvalue *key = component->u.members + 2 * m;
                if (jstr_is(key, "domain_display")) continue;
                if (!first) jb_lit(b, ", ");
                first = 0;
                jb_value(b, key);
                jb_lit(b, ": ");
                jb_value(b, key + 1);
            }
            jb_lit(b, "}");
        }
        jb_lit(b, "]");
    }
    jb_lit(b, "}}");
}

static void selection_rows(jbuf *b, const jvalue *slots, const char *irrep_key, int irrep_list)
{
    jb_lit(b, "[");
    for (int64_t i = 0; jis(slots, J_ARRAY) && i < slots->count; i++) {
        const jvalue *slot = jat(slots, i);
        if (i) jb_lit(b, ", ");
        jb_lit(b, "{\"slot\": ");
        put_or_null(b, member(slot, "slot"));
        jb_lit(b, ", \"kpoint\": ");
        kpoint_result(b, member(slot, "kpoint"));
        jb_lit(b, ", \"");
        jb_lit(b, irrep_key);
        jb_lit(b, "\": ");
        if (irrep_list) {
            const jvalue *rows = member(slot, "irreps");
            jb_lit(b, "[");
            for (int64_t r = 0; jis(rows, J_ARRAY) && r < rows->count; r++) {
                if (r) jb_lit(b, ", ");
                irrep_result(b, jat(rows, r));
            }
            jb_lit(b, "]");
        } else {
            irrep_result(b, member(slot, "irrep"));
        }
        jb_lit(b, "}");
    }
    jb_lit(b, "]");
}

/* The public mode details with definition ids and the selected subgroup. */
static void put_field(jbuf *b, const char *key, const jvalue *value, const expand_t *x)
{
    jb_lit(b, "\"");
    jb_lit(b, key);
    jb_lit(b, "\": ");
    if (value) put_copy(b, value, 1, x);
    else jb_lit(b, "null");
}

/* one mode definition: what identifies the mode apart from the choices (irrep basis, star arm, signs, the basis inside
 * a degenerate space) that present it */
static void definition_result(jbuf *b, const jvalue *definition, const char *kind, int64_t number, const expand_t *x)
{
    const jvalue *mode = member(definition, "mode"), *factor = member(definition, "factor");
    char id[64];
    snprintf(id, sizeof id, "{\"definition_id\": \"%s-%lld\", ", kind, (long long)number);
    jb_lit(b, id);
    put_field(b, "role", member(definition, "role"), x);
    jb_lit(b, ", ");
    put_field(b, "slot", member(factor, "slot"), x);
    jb_lit(b, ", ");
    put_field(b, "atom_block", member(definition, "atom_block"), x);
    jb_lit(b, ", \"invariants\": {");
    static const char *const invariant[] = {"kind", "k", "irrep", "site_label", "site_type", "wyckoff"};
    for (size_t i = 0; i < sizeof invariant / sizeof invariant[0]; i++) {
        if (i) jb_lit(b, ", ");
        put_field(b, invariant[i], member(mode, invariant[i]), x);
    }
    jb_lit(b, ", \"site_irrep\": ");
    const jvalue *site_irrep = member(mode, "site_irrep");
    if (jis(site_irrep, J_STRING)) {
        size_t n = 0;
        while (n < (size_t)site_irrep->count && site_irrep->u.s[n] != '(') n++;
        jb_string(b, site_irrep->u.s, n);
    } else {
        jb_lit(b, "null");
    }
    jb_lit(b, "}, \"setting\": {");
    put_field(b, "label", member(definition, "label"), x);
    jb_lit(b, ", ");
    put_field(b, "direction", member(mode, "direction"), x);
    jb_lit(b, ", ");
    put_field(b, "kvector", member(mode, "kvector"), x);
    jb_lit(b, ", ");
    put_field(b, "site_irrep", site_irrep, x);
    jb_lit(b, ", ");
    put_field(b, "parameters", member(factor, "parameters"), x);
    jb_lit(b, ", ");
    put_field(b, "normfactor", member(definition, "normfactor"), x);
    jb_lit(b, ", ");
    put_field(b, "vectors", member(definition, "vectors"), x);
    jb_lit(b, "}}");
}

static void mode_details_result(jbuf *b, const jvalue *details, const jvalue *orderparam, const expand_t *x)
{
    jb_lit(b, "{");
    int first = 1;
    for (int64_t i = 0; jis(details, J_OBJECT) && i < details->count; i++) {
        const jvalue *key = details->u.members + 2 * i, *value = key + 1;
        if (key->count && key->u.s[0] == '_') continue;
        if (!first) jb_lit(b, ", ");
        first = 0;
        jb_value(b, key);
        jb_lit(b, ": ");
        int definitions = jstr_is(key, "displacive_definitions") || jstr_is(key, "magnetic_definitions") ||
                          jstr_is(key, "occupancy_definitions");
        if (definitions && jis(value, J_ARRAY)) {
            const char *kind = jstr_is(key, "displacive_definitions") ? "displacive"
                               : jstr_is(key, "magnetic_definitions") ? "magnetic"
                                                                      : "occupancy";
            jb_lit(b, "[");
            for (int64_t d = 0; d < value->count; d++) {
                const jvalue *definition = jat(value, d);
                if (d) jb_lit(b, ", ");
                if (!jis(definition, J_OBJECT)) {
                    put_copy(b, definition, 1, x);
                    continue;
                }
                definition_result(b, definition, kind, d + 1, x);
            }
            jb_lit(b, "]");
        } else if (jstr_is(key, "subgroup_details") && jis(value, J_OBJECT)) {
            jb_lit(b, "{\"subgroup\": ");
            subgroup_result(b, member(member(orderparam, "isotropy"), "subgroup"));
            jb_lit(b, ", \"setting\": {");
            put_field(b, "basis", member(member(orderparam, "isotropy"), "basis_text"), x);
            jb_lit(b, ", ");
            put_field(b, "origin", member(value, "origin"), x);
            jb_lit(b, ", ");
            put_field(b, "parent_setting", member(value, "parent_inter_setting_id"), x);
            jb_lit(b, "}}");
        } else {
            put_copy(b, value, 1, x);
        }
    }
    jb_lit(b, "}");
}

static void info_sites(jbuf *b, const jvalue *atoms)
{
    static const char *const fields[] = {"type",       "label",      "wyckoff_multiplicity", "wyckoff",         "wyckoff_formula", "wyckoff_parameters",
                                         "wyckoff_params", "coordinates_realized", "fract", "occupancy", "wyckoff_mapping_error"};
    jb_lit(b, "[");
    for (int64_t i = 0; jis(atoms, J_ARRAY) && i < atoms->count; i++) {
        if (i) jb_lit(b, ", ");
        jb_lit(b, "{");
        int first = 1;
        put_fields(b, jat(atoms, i), fields, (int)(sizeof fields / sizeof fields[0]), &first);
        jb_lit(b, "}");
    }
    jb_lit(b, "]");
}
static void compact(jbuf *b, const jvalue *value, const char *stage, const expand_t *x)
{
    const jvalue *input = member(value, "input"), *selected = member(value, "selected");
    const jvalue *space_group = member(value, "space_group");
    if (!jtruthy(space_group)) space_group = member(input, "parent");
    const jvalue *symbol = member(member(input, "parent"), "symbol");
    jbuf symbol_text = {0};
    jb_reset(&symbol_text);
    if (jtruthy(symbol)) py_str_of(&symbol_text, symbol);
    size_t a = 0, e = symbol_text.used;
    while (a < e && sel_is_space(symbol_text.data[a])) a++;
    while (e > a && sel_is_space(symbol_text.data[e - 1])) e--;
    jb_lit(b, "{\"schema\": \"APOSTRUCT.cli.");
    jb_lit(b, stage);
    jb_lit(b, "\", \"space_group\": {");
    int first = 1;
    for (int64_t i = 0; jtruthy(space_group) && jis(space_group, J_OBJECT) && i < space_group->count; i++) {
        const jvalue *key = space_group->u.members + 2 * i;
        if (e > a && jstr_is(key, "symbol")) continue;
        if (!first) jb_lit(b, ", ");
        first = 0;
        jb_value(b, key);
        jb_lit(b, ": ");
        jb_value(b, key + 1);
    }
    if (e > a) {
        if (!first) jb_lit(b, ", ");
        jb_lit(b, "\"symbol\": ");
        jb_string(b, symbol_text.data + a, e - a);
    }
    jb_free(&symbol_text);
    jb_lit(b, "}");
    if (!strcmp(stage, "info")) {
        jb_lit(b, ", \"formula\": ");
        put_or_null(b, member(input, "formula"));
        jb_lit(b, ", \"cell\": ");
        put_or_null(b, member(input, "lattice"));
        jb_lit(b, ", \"sites\": ");
        info_sites(b, member(input, "atom_sites"));
    } else if (!strcmp(stage, "kpoints")) {
        const jvalue *rows = member(value, "kpoints");
        jb_lit(b, ", \"kpoints\": [");
        for (int64_t i = 0; jis(rows, J_ARRAY) && i < rows->count; i++) {
            if (i) jb_lit(b, ", ");
            kpoint_result(b, jat(rows, i));
        }
        jb_lit(b, "]");
    } else {
        jb_lit(b, ", \"distortion_selection\": ");
        const jvalue *distortion = member(value, "distortion_selection");
        if (jtruthy(distortion)) jb_value(b, distortion);
        else jb_lit(b, "{}");
        jb_lit(b, ", \"selections\": ");
        if (!strcmp(stage, "irreps")) {
            selection_rows(b, member(selected, "irrep_slots"), "irreps", 1);
        } else {
            selection_rows(b, member(selected, "selection_slots"), "irrep", 0);
        }
        if (!strcmp(stage, "opds")) {
            const jvalue *rows = member(selected, "opd_rows");
            jb_lit(b, ", \"opds\": [");
            for (int64_t i = 0; jis(rows, J_ARRAY) && i < rows->count; i++) {
                if (i) jb_lit(b, ", ");
                opd_result(b, jat(rows, i));
            }
            jb_lit(b, "]");
        } else if (!strcmp(stage, "modes")) {
            const jvalue *orderparam = member(selected, "orderparam");
            jb_lit(b, ", \"opd\": ");
            opd_result(b, orderparam);
            jb_lit(b, ", \"mode_details\": ");
            mode_details_result(b, member(selected, "mode_details"), orderparam, x);
        }
    }
    jb_lit(b, "}");
}

/* ------------------------------------------------------------ entry */

/* The CIF summary of the case's structure file. */
static int read_structure(const jvalue *request, jdoc **summary)
{
    *summary = 0;
    const jvalue *structure = member(request, "structure");
    if (!jis(structure, J_STRING)) return md_fail(MDE_VALUE, "case has no structure");
    char path[4096];
    if (structure->count >= (int64_t)sizeof path) return md_fail(MDE_VALUE, "structure path is too long");
    memcpy(path, structure->u.s, (size_t)structure->count);
    path[structure->count] = 0;
    struct stat info;
    if (stat(path, &info) || !S_ISREG(info.st_mode)) return md_fail(MDE_VALUE, "structure CIF does not exist: %s", path);
    const char *summary_text = cif_summary_json(structure->u.s, structure->count, 1);
    if (!summary_text) return md_fail(MDE_MEMORY, "CIF reader ran out of memory");
    *summary = jdoc_parse(summary_text, strlen(summary_text));
    if (!*summary) return md_fail(MDE_MEMORY, "out of memory");
    const jvalue *cif_error = member(jdoc_root(*summary), "error");
    if (cif_error) {
        const jvalue *message = member(cif_error, "message");
        int status = md_fail(MDE_VALUE, "%.*s", message ? (int)message->count : 0, message ? message->u.s : "");
        jdoc_free(*summary);
        *summary = 0;
        return status;
    }
    return 0;
}

/* case_input.pipeline_request: {"cif_info": ..., the selection options} */
static int request_of(const jvalue *request, jbuf *b)
{
    jdoc *summary;
    MD_TRY(read_structure(request, &summary));
    jb_lit(b, "{\"cif_info\": ");
    jb_value(b, jdoc_root(summary));
    jb_lit(b, ", ");
    int status = pipeline_options(request, jdoc_root(summary), jtruthy(member(request, "require_irreps")),
                                  jtruthy(member(request, "allow_empty_distortion")), b);
    jb_lit(b, "}");
    jdoc_free(summary);
    return status ? status : b->failed ? md_fail(MDE_MEMORY, "out of memory") : 0;
}

int spec_gamma(const jvalue *cif_info, const int *sites, const int *kinds, int nsites, int strain, const rat *kvec, int sg, int setting,
               const void *names, jbuf *b);
void *spec_names_new(int count, const int *gids, const char *const *symbols, const int64_t *kparam);
const char *k_catalog_json(int sg);
int k_coordinate_map_words(int sg, int kslot, int64_t *out);

static int run(const jvalue *request, jbuf *b)
{
    const jvalue *stage_value = member(request, "stage");
    char stage[16] = {0};
    if (!jis(stage_value, J_STRING) || stage_value->count >= (int64_t)sizeof stage) return md_fail(MDE_VALUE, "unknown stage");
    memcpy(stage, stage_value->u.s, (size_t)stage_value->count);
    int full_state = jtruthy(member(request, "full_state")), text_format = jstr_is(member(request, "format"), "text");
    int64_t indent = 2;
    if (member(request, "indent")) py_int_of(member(request, "indent"), &indent);
    jdoc *summary;
    MD_TRY(read_structure(request, &summary));
    const jvalue *cif_info = jdoc_root(summary);
    jbuf state = {0}, compact_text = {0};
    jb_reset(&state);
    jb_reset(&compact_text);
    int status = 0;
    jdoc *sdoc = 0;
    if (!strcmp(stage, "spec")) {
        jbuf options = {0};
        jb_reset(&options);
        jb_lit(&options, "{");
        status = pipeline_options(request, cif_info, 0, 0, &options);
        jb_lit(&options, "}");
        jdoc *odoc = status ? 0 : jdoc_parse(options.data, options.used);
        /* the selected sites in CIF order, each with its kinds (1 displacive, 2 magnetic) */
        const jvalue *sel = odoc ? member(jdoc_root(odoc), "distortion_selection") : 0;
        const jvalue *lists[2] = {member(sel, "displacive_sites"), member(sel, "magnetic_sites")};
        int64_t natoms = jis(member(cif_info, "atom_sites"), J_ARRAY) ? member(cif_info, "atom_sites")->count : 0;
        int *flags = calloc((size_t)natoms + 1, sizeof(int)), nsites = 0;
        for (int kind = 0; kind < 2; kind++)
            for (int64_t i = 0; jis(lists[kind], J_ARRAY) && i < lists[kind]->count; i++) {
                int64_t v = 0;
                py_int_of(jat(lists[kind], i), &v);
                if (v >= 1 && v <= natoms) flags[v - 1] |= 1 << kind;
            }
        int *sites = calloc((size_t)natoms + 1, sizeof(int)), *kinds = calloc((size_t)natoms + 1, sizeof(int));
        for (int64_t i = 0; i < natoms; i++)
            if (flags[i]) {
                sites[nsites] = (int)i + 1;
                kinds[nsites++] = flags[i];
            }
        int strain = jtruthy(member(sel, "strain"));
        free(flags);
        rat kvec[3] = {rat_int(0), rat_int(0), rat_int(0)};
        int64_t spec_setting = 0, spec_sg = 0;
        void *spec_names = 0;
        py_int_of(member(member(cif_info, "parent"), "number"), &spec_sg);
        if (!status && odoc) {
            /* the selected k point's coordinates: the catalog map origin + parameters . columns */
            jbuf parent = {0}, kstate = {0};
            jb_reset(&parent);
            jb_reset(&kstate);
            status = par_parent_state(cif_info, &parent);
            jdoc *pdoc = status ? 0 : jdoc_parse(parent.data, parent.used);
            if (!status && pdoc) status = sel_k_points(jdoc_root(pdoc), jdoc_root(odoc), &kstate);
            jdoc *kdoc = status ? 0 : jdoc_parse(kstate.data, kstate.used);
            const jvalue *kp = kdoc ? member(member(jdoc_root(kdoc), "selected"), "kpoint") : 0;
            int64_t kslot = 0, sg = 0, words[29];
            if (pdoc) py_int_of(member(member(jdoc_root(pdoc), "_internal"), "parent_inter_setting_id"), &spec_setting);
            py_int_of(member(member(cif_info, "parent"), "number"), &sg);
            if (!status && (!kp || py_int_of(member(kp, "kslot"), &kslot) || k_coordinate_map_words((int)sg, (int)kslot, words)))
                status = md_fail(MDE_VALUE, "spec: selected k point has no coordinate map");
            const jvalue *params = kp ? member(kp, "parameters") : 0;
            rat_overflow = 0;
            for (int axis = 0; axis < 3 && !status; axis++) kvec[axis] = rat_make(words[5 + 2 * axis], words[6 + 2 * axis]);
            for (int i = 0; i < (int)words[1] && !status; i++) {
                char key[2] = {(char)words[2 + i], 0};
                const jvalue *value = member(params, key);
                rat v;
                if (!jis(value, J_STRING) || py_fraction_of_text(value->u.s, value->count, &v)) {
                    status = md_fail(MDE_VALUE, "spec: k parameter %s is missing", key);
                    break;
                }
                for (int axis = 0; axis < 3; axis++)
                    kvec[axis] = rat_add(kvec[axis], rat_mul(v, rat_make(words[11 + 6 * i + 2 * axis], words[12 + 6 * i + 2 * axis])));
            }
            if (rat_overflow) status = md_fail(MDE_ARITH, "spec: k out of range");
            /* naming candidates: the catalog irreps of this k point, and its parameters as (p1, p2, p3, den) */
            if (!status) {
                int64_t kparam[4] = {0, 0, 0, 1}, visible[6];
                for (int axis = 0; axis < 3; axis++) {
                    visible[2 * axis] = (int64_t)kvec[axis].n;
                    visible[2 * axis + 1] = (int64_t)kvec[axis].d;
                }
                if (k_source_kparam_of_visible((int)sg, (int)kslot, visible, kparam) < 0) status = md_fail(MDE_VALUE, "spec: k parameters");
                const char *catalog = k_catalog_json((int)sg);
                jdoc *cdoc = catalog ? jdoc_parse(catalog, strlen(catalog)) : 0;
                const jvalue *kps = cdoc ? member(jdoc_root(cdoc), "kpoints") : 0;
                int gids[64], count = 0;
                const char *symbols[64];
                char held[64][24];
                for (int64_t q = 0; jis(kps, J_ARRAY) && q < kps->count; q++) {
                    int64_t slot = 0;
                    py_int_of(member(jat(kps, q), "kslot"), &slot);
                    if (slot != kslot) continue;
                    const jvalue *list = member(jat(kps, q), "irreps");
                    for (int64_t r = 0; jis(list, J_ARRAY) && r < list->count && count < 64; r++) {
                        int64_t gid = 0;
                        const jvalue *sym = member(jat(list, r), "symbol");
                        py_int_of(member(jat(list, r), "gid"), &gid);
                        if (!jis(sym, J_STRING)) continue;
                        snprintf(held[count], sizeof held[count], "%.*s", (int)sym->count, sym->u.s);
                        symbols[count] = held[count];
                        gids[count++] = (int)gid;
                    }
                }
                spec_names = spec_names_new(count, gids, symbols, kparam);
                jdoc_free(cdoc);
            }
            jdoc_free(kdoc);
            jdoc_free(pdoc);
            jb_free(&parent);
            jb_free(&kstate);
        }
        if (!status) status = odoc ? spec_gamma(cif_info, sites, kinds, nsites, strain, kvec, (int)spec_sg, (int)spec_setting, spec_names, b) : md_fail(MDE_MEMORY, "out of memory");
        free(sites);
        free(kinds);
        free(spec_names);
        jdoc_free(odoc);
        jb_free(&options);
        jdoc_free(summary);
        return status;
    }
    if (!strcmp(stage, "info")) {
        status = par_parent_state(cif_info, &state);
    } else {
        jbuf options = {0};
        jb_reset(&options);
        int require_irreps = !strcmp(stage, "opds") || !strcmp(stage, "modes");
        jb_lit(&options, "{");
        status = pipeline_options(request, cif_info, require_irreps, !strcmp(stage, "kpoints"), &options);
        jb_lit(&options, ", \"allow_empty_distortion\": false, \"opd_row_limit\": null");
        jb_lit(&options, require_irreps ? ", \"include_opd\": true" : ", \"include_opd\": false");
        jb_lit(&options, !strcmp(stage, "modes") ? ", \"selected_opd_only\": true}" : ", \"selected_opd_only\": false}");
        jdoc *odoc = status ? 0 : jdoc_parse(options.data, options.used);
        if (!status && !odoc) status = md_fail(MDE_MEMORY, "out of memory");
        if (!status && !strcmp(stage, "kpoints")) {
            jbuf parent = {0};
            jb_reset(&parent);
            status = par_parent_state(cif_info, &parent);
            jdoc *pdoc = status ? 0 : jdoc_parse(parent.data, parent.used);
            if (!status && !pdoc) status = md_fail(MDE_MEMORY, "out of memory");
            if (!status) status = sel_k_points(jdoc_root(pdoc), jdoc_root(odoc), &state);
            jdoc_free(pdoc);
            jb_free(&parent);
        } else if (!status) {
            jbuf call = {0};
            jb_reset(&call);
            jb_lit(&call, "{\"call\": \"state\", \"cif_info\": ");
            jb_value(&call, cif_info);
            jb_lit(&call, ", \"options\": ");
            jb_raw(&call, options.data, options.used);
            jb_lit(&call, !strcmp(stage, "modes") ? ", \"include_mode_details\": true}" : ", \"include_mode_details\": false}");
            jdoc *cdoc = call.failed ? 0 : jdoc_parse(call.data, call.used);
            jb_free(&call);
            if (!cdoc) status = md_fail(MDE_MEMORY, "out of memory");
            else status = st_state(jdoc_root(cdoc), &state);
            jdoc_free(cdoc);
        }
        jdoc_free(odoc);
        jb_free(&options);
    }
    jdoc_free(summary);
    if (!status) {
        sdoc = jdoc_parse(state.data, state.used);
        jb_free(&state);
        if (!sdoc) status = md_fail(MDE_MEMORY, "out of memory");
    }
    if (!status) {
        const jvalue *root = jdoc_root(sdoc), *error = member(root, "error");
        if (error && root->count == 1) {
            const jvalue *message = member(error, "message");
            int64_t kind = 0, code = 0;
            py_int_of(member(error, "kind"), &kind);
            py_int_of(member(error, "status"), &code);
            md_fail((int)kind, "%.*s", message ? (int)message->count : 0, message ? message->u.s : "");
            md_err.status = (int)code;
            status = (int)kind;
        }
    }
    if (!status) {
        expand_t x;
        x.nvectors = mode_details_vectors(&x.vectors);
        const jvalue *root = jdoc_root(sdoc);
        if (text_format) {
            jb_reset(b);
            status = md_render_text(root, x.vectors, x.nvectors, b);
            jdoc_free(sdoc);
            jb_free(&compact_text);
            return status;
        }
        if (full_state) put_copy(&compact_text, root, 0, &x);
        else compact(&compact_text, root, stage, &x);
        const jvalue *details = member(member(root, "selected"), "mode_details");
        run_status = !strcmp(stage, "modes") && jstr_is(member(details, "status"), "error") ? 2 : 0;
        jdoc_free(sdoc);
        sdoc = 0;
        status = print_json(&compact_text, (int)indent, b, &x);
    }
    jdoc_free(sdoc);
    jb_free(&state);
    jb_free(&compact_text);
    return status;
}

/* The printed result of one normalized CLI case; cli_run_status() is 0, 2
 * for a failed mode-detail answer, or -1 when the text is an error record. */
const char *cli_run_json(const char *text, int64_t length)
{
    if (out.cap > ((size_t)16 << 20)) jb_free(&out);
    jb_reset(&out);
    run_status = 0;
    jdoc *doc = jdoc_parse(text, (size_t)length);
    int status;
    if (!doc) status = md_fail(MDE_VALUE, "run request is not JSON");
    else if (jstr_is(jget(jdoc_root(doc), "call"), "request")) status = request_of(jdoc_root(doc), &out);
    else status = run(jdoc_root(doc), &out);
    jdoc_free(doc);
    if (status) {
        run_status = -1;
        jb_reset(&out);
        jb_lit(&out, "{\"error\": {\"kind\": ");
        jb_int(&out, md_err.kind);
        jb_lit(&out, ", \"status\": ");
        jb_int(&out, md_err.status);
        jb_lit(&out, ", \"message\": ");
        jb_cstring(&out, md_err.message);
        jb_lit(&out, "}}");
    }
    return out.failed ? 0 : out.data;
}

int cli_run_status(void) { return run_status; }

int64_t cli_run_length(void) { return (int64_t)out.used; }
