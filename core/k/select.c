/* The k point and irrep selections of the selection stage. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "sel.h"
#include "tables.h"

const char *k_catalog_json(int sg);
int k_coordinate_map_words(int sg, int kslot, int64_t *out);
int k_source_kparam_of_visible(int sg, int kslot, const int64_t *visible_words, int64_t *out);
int k_source_little_order(int sg, int kslot, int count, const int64_t *parameter_words, int *order_out);
int vector_mode_total(int sg, int gid, int row_count, const int *row_ids, int axial, int64_t *total);
int strain_irrep(int sg, int gid);

static struct {
    section_view k_dim, mag_pointer;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (sel_tables() || section_find("k.slots.dim", &T.k_dim) ||
        section_find("isotropy.magnetic.by_old_start", &T.mag_pointer))
        return md_fail(MDE_VALUE, "k-point Source tables are inconsistent");
    T.ready = 1;
    return 0;
}

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

static int64_t member_int(const jvalue *object, const char *key)
{
    int64_t v = 0;
    const jvalue *m = member(object, key);
    if (jtruthy(m)) py_int_of(m, &v);
    return v;
}

static int casefold_equal(const jbuf *a, const jbuf *b)
{
    if (a->used != b->used) return 0;
    for (size_t i = 0; i < a->used; i++) {
        char x = a->data[i], y = b->data[i];
        if (x >= 'A' && x <= 'Z') x = (char)(x - 'A' + 'a');
        if (y >= 'A' && y <= 'Z') y = (char)(y - 'A' + 'a');
        if (x != y) return 0;
    }
    return 1;
}
static int select_k(const jvalue *kpoints, const jvalue *label, int64_t index, const jvalue **out)
{
    int64_t n = jis(kpoints, J_ARRAY) ? kpoints->count : 0;
    if (jtruthy(label)) {
        jbuf want = {0};
        jb_reset(&want);
        py_str_of(&want, label);
        for (int64_t i = 0; i < n; i++) {
            jbuf have = {0};
            jb_reset(&have);
            py_str_of(&have, member(jat(kpoints, i), "label"));
            int same = casefold_equal(&have, &want);
            jb_free(&have);
            if (same) {
                jb_free(&want);
                *out = jat(kpoints, i);
                return 0;
            }
        }
        jbuf repr = {0};
        jb_reset(&repr);
        if (jis(label, J_STRING)) py_repr_text(&repr, label->u.s, label->count);
        else py_repr_of(&repr, label);
        int status = md_fail(MDE_KEY, "k label %.*s not found", (int)repr.used, repr.data);
        jb_free(&repr);
        jb_free(&want);
        return status;
    }
    if (index < 1 || index > n) return md_fail(MDE_INDEX, "k index %lld out of range 1..%lld", (long long)index, (long long)n);
    *out = jat(kpoints, index - 1);
    return 0;
}

static void put_frac(jbuf *b, rat v)
{
    char text[96];
    text[0] = '"';
    {
        int n = py_rat_text(v, text + 1);
        text[n + 1] = '"';
        text[n + 2] = 0;
    }
    jb_lit(b, text);
}

static int is_space(char c) { return sel_is_space(c); }

/* the miller_love_parameters of one selected k point */
static int miller_love(jbuf *b, int sg, const jvalue *selected, const sel_params *params)
{
    int64_t dimension = member_int(selected, "dimension"), kslot = member_int(selected, "kslot");
    int64_t map[29];
    int status = k_coordinate_map_words(sg, (int)kslot, map);
    if (status) return md_fail(status == 3 ? MDE_ARITH : MDE_VALUE, "K-family display expression is not affine rational");
    int count = (int)map[1];
    jb_lit(b, "{");
    int all = dimension > 0;
    for (int i = 0; all && i < count; i++) {
        char name[2] = {(char)map[2 + i], 0};
        const char *raw = sel_param(params, name);
        int nonblank = 0;
        for (const char *c = raw; c && *c && !nonblank; c++) nonblank = !is_space(*c);
        all = nonblank;
    }
    if (all && params->count) {
        int lattice = geo_sg_lattice(sg);
        int64_t dim = section_int(&T.k_dim, (int64_t)(lattice - 1) * 27 + kslot - 1);
        int64_t visible[6], source[4];
        int has = dim > 0 ? sel_visible_k_vector(sg, (int)kslot, params, visible) : 0;
        if (has < 0) return md_fail(-has, "k-vector Source tables are inconsistent");
        int found = has ? k_source_kparam_of_visible(sg, (int)kslot, visible, source) : 0;
        if (found < 0) return md_fail(found == -2 ? MDE_ARITH : MDE_VALUE, "k-vector Source tables are inconsistent");
        if (found) {
            rat values[3];
            int64_t words[6];
            for (int i = 0; i < dimension && i < 3; i++) {
                values[i] = rat_make(source[i], source[3]);
                words[2 * i] = (int64_t)values[i].n;
                words[2 * i + 1] = (int64_t)values[i].d;
            }
            if (dim > 0) {
                /* the little-group order of the family at a generic parameter */
                static const int64_t generic[6] = {1, 1009, 1, 1013, 1, 1019};
                int order = 0, expected = 0;
                int code = k_source_little_order(sg, (int)kslot, (int)dimension, words, &order);
                if (!code) code = k_source_little_order(sg, (int)kslot, (int)dimension, generic, &expected);
                if (code == -3) return md_fail(MDE_VALUE, "K-family parameter count does not match its Source slot");
                if (code) return md_fail(code == -2 ? MDE_ARITH : MDE_VALUE, "k-vector Source tables are inconsistent");
                if (order > expected) {
                    jbuf label = {0};
                    jb_reset(&label);
                    py_str_of(&label, member(selected, "label"));
                    int fail = md_fail(MDE_VALUE, "%.*s parameters are at a special K point; select the corresponding fixed K label",
                                       (int)label.used, label.data ? label.data : "");
                    jb_free(&label);
                    return fail;
                }
            }
            static const char *names[3] = {"alpha", "beta", "gamma"};
            for (int i = 0; i < dimension && i < 3; i++) {
                if (i) jb_lit(b, ", ");
                jb_cstring(b, names[i]);
                jb_lit(b, ": ");
                put_frac(b, values[i]);
            }
        }
    }
    jb_lit(b, "}");
    return 0;
}

static void put_member(jbuf *b, const jvalue *object, const char *key, int first)
{
    if (!first) jb_lit(b, ", ");
    jb_cstring(b, key);
    jb_lit(b, ": ");
    const jvalue *v = member(object, key);
    if (v) jb_value(b, v);
    else jb_lit(b, "null");
}

int sel_k_points(const jvalue *state, const jvalue *opts, jbuf *out)
{
    MD_TRY(load());
    int sg = (int)member_int(member(member(state, "input"), "parent"), "number");
    jbuf kp = {0}, dp = {0};
    jb_reset(&kp);
    jb_reset(&dp);
    sel_normalize_params(&kp, member(opts, "k_params"));
    sel_normalize_params(&dp, member(opts, "display_k_params"));
    if (dp.used <= 2) {
        jb_reset(&dp);
        jb_raw(&dp, kp.data, kp.used);
    }
    const char *catalog_text = k_catalog_json(sg);
    if (!catalog_text) return md_fail(MDE_VALUE, "k-point catalog Source tables are inconsistent");
    jdoc *cdoc = jdoc_parse(catalog_text, strlen(catalog_text));
    const jvalue *catalog = jdoc_root(cdoc), *kpoints = member(catalog, "kpoints");
    const jvalue *k_label_value = member(opts, "k_label");
    int has_label_key = jis(opts, J_OBJECT) && jget(opts, "k_label");
    int64_t k_index = 1;
    const jvalue *index_value = member(opts, "k_index");
    if (jtruthy(index_value)) py_int_of(index_value, &k_index);
    else if (jis(index_value, J_INT)) k_index = index_value->u.i;
    jbuf specs_text = {0};
    jb_reset(&specs_text);
    const jvalue *selections = member(opts, "selections");
    jdoc *sdoc = 0;
    const jvalue *specs;
    if (jtruthy(selections)) {
        specs = selections;
    } else {
        jb_lit(&specs_text, "[{\"k_label\": ");
        if (has_label_key) jb_value(&specs_text, k_label_value);
        else jb_lit(&specs_text, "\"GM\"");
        jb_lit(&specs_text, ", \"k_index\": ");
        jb_int(&specs_text, k_index);
        jb_lit(&specs_text, ", \"k_params\": ");
        jb_raw(&specs_text, kp.data, kp.used);
        jb_lit(&specs_text, ", \"display_k_params\": ");
        jb_raw(&specs_text, dp.data, dp.used);
        jb_lit(&specs_text, "}]");
        sdoc = jdoc_parse(specs_text.data, specs_text.used);
        specs = jdoc_root(sdoc);
    }
    jbuf selected = {0}, first = {0};
    jb_reset(&selected);
    jb_reset(&first);
    jb_lit(&selected, "[");
    int status = 0;
    for (int64_t i = 0; i < specs->count && !status; i++) {
        const jvalue *spec = jat(specs, i), *label;
        if (jis(spec, J_OBJECT) && jget(spec, "k_label")) label = jget(spec, "k_label");
        else label = has_label_key ? k_label_value : 0;
        jdoc *gm = 0;
        if (!(jis(spec, J_OBJECT) && jget(spec, "k_label")) && !has_label_key) {
            gm = jdoc_parse("\"GM\"", 4);
            label = jdoc_root(gm);
        }
        int64_t index = k_index;
        if (jtruthy(member(spec, "k_index"))) py_int_of(member(spec, "k_index"), &index);
        const jvalue *row = 0;
        status = select_k(kpoints, label, index, &row);
        jdoc_free(gm);
        if (status) break;
        jbuf sp = {0}, sdp = {0};
        jb_reset(&sp);
        jb_reset(&sdp);
        sel_normalize_params(&sp, member(spec, "k_params"));
        sel_normalize_params(&sdp, member(spec, "display_k_params"));
        if (sdp.used <= 2) {
            jb_reset(&sdp);
            jb_raw(&sdp, sp.data, sp.used);
        }
        jdoc *pdoc = jdoc_parse(sp.data, sp.used);
        sel_params params;
        sel_params_of(jdoc_root(pdoc), &params);
        jbuf ml = {0};
        jb_reset(&ml);
        status = miller_love(&ml, sg, row, &params);
        if (!status) {
            jbuf entry = {0};
            jb_reset(&entry);
            jb_lit(&entry, "{\"slot\": ");
            jb_int(&entry, i + 1);
            jb_lit(&entry, ", ");
            size_t body = entry.used;
            static const char *keys[6] = {"kslot", "kid", "label", "display", "kvector", "display_kvector"};
            for (int k = 0; k < 6; k++) put_member(&entry, row, keys[k], k == 0);
            jb_lit(&entry, ", \"parameters\": ");
            jb_raw(&entry, sp.data, sp.used);
            jb_lit(&entry, ", \"display_parameters\": ");
            jb_raw(&entry, sdp.data, sdp.used);
            jb_lit(&entry, ", \"miller_love_parameters\": ");
            jb_raw(&entry, ml.data, ml.used);
            jb_lit(&entry, ", \"dimension\": ");
            jb_int(&entry, member_int(row, "dimension"));
            jb_lit(&entry, jtruthy(member(row, "sg_specific")) ? ", \"sg_specific\": true" : ", \"sg_specific\": false");
            put_member(&entry, row, "star", 0);
            jb_lit(&entry, "}");
            if (i) jb_lit(&selected, ", ");
            jb_raw(&selected, entry.data, entry.used);
            if (!i) {
                jb_lit(&first, "{");
                jb_raw(&first, entry.data + body, entry.used - body);
            }
            jb_free(&entry);
        }
        jdoc_free(pdoc);
        jb_free(&sp);
        jb_free(&sdp);
        jb_free(&ml);
    }
    jb_lit(&selected, "]");
    if (!status) {
        jbuf list = {0}, sel = {0}, internal = {0}, space_group = {0};
        jb_reset(&list);
        jb_reset(&sel);
        jb_reset(&internal);
        jb_reset(&space_group);
        jb_value(&space_group, member(catalog, "space_group"));
        jb_lit(&list, "[");
        for (int64_t i = 0; i < kpoints->count; i++) {
            const jvalue *k = jat(kpoints, i);
            if (i) jb_lit(&list, ", ");
            jb_lit(&list, "{");
            static const char *keys[7] = {"kslot", "kid", "label", "display", "kvector", "display_kvector", "dimension"};
            for (int j = 0; j < 7; j++) put_member(&list, k, keys[j], j == 0);
            jb_lit(&list, jtruthy(member(k, "sg_specific")) ? ", \"sg_specific\": true" : ", \"sg_specific\": false");
            put_member(&list, k, "star_size", 0);
            put_member(&list, k, "n_irreps", 0);
            jb_lit(&list, "}");
        }
        jb_lit(&list, "]");
        jb_lit(&sel, "{\"kpoint\": ");
        jb_raw(&sel, first.data, first.used);
        jb_lit(&sel, ", \"kpoints\": ");
        jb_raw(&sel, selected.data, selected.used);
        jb_lit(&sel, "}");
        kv add[1] = {{"k_catalog", catalog_text, strlen(catalog_text)}};
        sel_merge(&internal, member(state, "_internal"), add, 1);
        kv items[5] = {{"schema", "\"APOSTRUCT.k_points.v1\"", strlen("\"APOSTRUCT.k_points.v1\"")},
                       {"space_group", space_group.data, space_group.used},
                       {"kpoints", list.data, list.used},
                       {"selected", sel.data, sel.used},
                       {"_internal", internal.data, internal.used}};
        sel_merge(out, state, items, 5);
        jb_free(&list);
        jb_free(&sel);
        jb_free(&internal);
        jb_free(&space_group);
    }
    jdoc_free(cdoc);
    jdoc_free(sdoc);
    jb_free(&specs_text);
    jb_free(&selected);
    jb_free(&first);
    jb_free(&kp);
    jb_free(&dp);
    return status;
}

/* ------------------------------------------------------------ irreps */

static int strain_visible(int sg, const jvalue *kpoint, const jvalue *candidate)
{
    const jvalue *label = member(kpoint, "label");
    if (!jstr_is(label, "GM")) return 0;
    jbuf symbol = {0};
    jb_reset(&symbol);
    if (jtruthy(member(candidate, "symbol"))) py_str_of(&symbol, member(candidate, "symbol"));
    int gm1 = (symbol.used == 3 && !memcmp(symbol.data, "GM1", 3)) || (symbol.used == 4 && !memcmp(symbol.data, "GM1+", 4));
    jb_free(&symbol);
    if (gm1) return 1;
    int64_t old_id = member_int(candidate, "old_id");
    return old_id > 0 && strain_irrep(sg, (int)member_int(candidate, "gid"));
}

static int has_magnetic_variant(const jvalue *candidate)
{
    int64_t old_id = member_int(candidate, "old_id");
    if (old_id <= 0) return 1;
    if (old_id >= section_count(&T.mag_pointer)) return 0;
    return section_int(&T.mag_pointer, old_id) > section_int(&T.mag_pointer, old_id - 1);
}

static int mode_total(int sg, const jvalue *candidate, const int *rows, int nrows, int axial, int64_t *total)
{
    int status = vector_mode_total(sg, (int)member_int(candidate, "gid"), nrows, rows, axial, total);
    if (!status) return 0;
    md_err.kind = MDE_PROJECTION_ENTRY;
    md_err.status = status;
    md_err.message[0] = 0;
    return MDE_PROJECTION_ENTRY;
}

static int int_list(const jvalue *v, int **out, int *n)
{
    *n = jis(v, J_ARRAY) ? (int)v->count : 0;
    *out = calloc((size_t)*n + 1, sizeof(int));
    if (!*out) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < *n; i++) {
        int64_t x;
        if (py_int_of(jat(v, i), &x)) return md_fail(MDE_TYPE, "row ids must be integers");
        (*out)[i] = (int)x;
    }
    return 0;
}

int sel_irreps(const jvalue *state, jbuf *out)
{
    MD_TRY(load());
    int sg = (int)member_int(member(member(state, "input"), "parent"), "number");
    const jvalue *internal = member(state, "_internal"), *catalog = member(internal, "k_catalog");
    jbuf mode_set = {0};
    jb_reset(&mode_set);
    const jvalue *ms = member(internal, "mode_set");
    if (!jtruthy(ms)) ms = member(state, "mode_set");
    if (jtruthy(ms)) py_str_of(&mode_set, ms);
    else jb_lit(&mode_set, "displacive");
    const jvalue *selection = member(internal, "distortion_selection");
    int include_strain;
    if (jis(selection, J_OBJECT) && jtruthy(selection) && jget(selection, "strain")) {
        include_strain = jtruthy(jget(selection, "strain"));
    } else {
        include_strain = (mode_set.used == 17 && !memcmp(mode_set.data, "strain_displacive", 17)) ||
                         (mode_set.used == 3 && !memcmp(mode_set.data, "all", 3));
    }
    int *displacive, *magnetic_rows, nd, nm;
    MD_TRY(int_list(member(internal, "displacive_row_ids"), &displacive, &nd));
    MD_TRY(int_list(member(internal, "magnetic_row_ids"), &magnetic_rows, &nm));
    const jvalue *kpoints = member(catalog, "kpoints");
    int64_t nk = jis(kpoints, J_ARRAY) ? kpoints->count : 0;
    jbuf *irreps = calloc((size_t)nk + 1, sizeof(jbuf)), *kp_json = calloc((size_t)nk + 1, sizeof(jbuf));
    int64_t *counts = calloc((size_t)nk + 1, sizeof(int64_t));
    if (!irreps || !kp_json || !counts) return md_fail(MDE_MEMORY, "out of memory");
    int status = 0;
    for (int64_t k = 0; k < nk && !status; k++) {
        const jvalue *kpoint = jat(kpoints, k), *candidates = member(kpoint, "irreps");
        int64_t nc = jis(candidates, J_ARRAY) ? candidates->count : 0;
        jb_reset(&irreps[k]);
        jb_lit(&irreps[k], "[");
        int64_t count = 0;
        for (int64_t c = 0; c < nc && !status; c++) {
            const jvalue *candidate = jat(candidates, c);
            int64_t total = 0;
            status = mode_total(sg, candidate, displacive, nd, 0, &total);
            if (status) break;
            int strain = include_strain && strain_visible(sg, kpoint, candidate);
            if (total <= 0 && !strain) continue;
            char number[32];
            snprintf(number, sizeof number, "%lld", (long long)total);
            kv add[2];
            int n = 0;
            if (total > 0) add[n++] = (kv){"displacive_mode_total", number, strlen(number)};
            if (strain) add[n++] = (kv){"strain_visible", "true", 4};
            if (count++) jb_lit(&irreps[k], ", ");
            sel_merge(&irreps[k], candidate, add, n);
        }
        if (status) break;
        int64_t magnetic_count = 0;
        jbuf magnetic = {0};
        jb_reset(&magnetic);
        for (int64_t c = 0; nm && c < nc && !status; c++) {
            const jvalue *candidate = jat(candidates, c);
            if (!has_magnetic_variant(candidate)) continue;
            int64_t total = 0;
            status = mode_total(sg, candidate, magnetic_rows, nm, 1, &total);
            if (status || total <= 0) continue;
            jbuf symbol = {0}, kov = {0}, ordinary = {0}, text = {0}, m = {0};
            jb_reset(&symbol);
            jb_reset(&kov);
            jb_reset(&ordinary);
            jb_reset(&text);
            jb_reset(&m);
            py_str_of(&text, member(candidate, "symbol"));
            jb_lit(&m, "m");
            jb_raw(&m, text.data, text.used);
            jb_string(&symbol, m.data, m.used);
            jb_reset(&text);
            if (jtruthy(member(candidate, "kov_label"))) py_str_of(&text, member(candidate, "kov_label"));
            jb_reset(&m);
            if (text.used && text.data[0] != 'm') jb_lit(&m, "m");
            jb_raw(&m, text.data, text.used);
            jb_string(&kov, m.data, m.used);
            const jvalue *ordinary_symbol = member(candidate, "symbol");
            if (ordinary_symbol) jb_value(&ordinary, ordinary_symbol);
            else jb_lit(&ordinary, "null");
            char number[32];
            snprintf(number, sizeof number, "%lld", (long long)total);
            kv add[6] = {{"symbol", symbol.data, symbol.used},
                         {"kov_label", kov.data, kov.used},
                         {"magnetic", "true", 4},
                         {"ordinary_symbol", ordinary.data, ordinary.used},
                         {"order_parameter_type", "\"magnetic\"", 10},
                         {"magnetic_mode_total", number, strlen(number)}};
            if (magnetic_count++) jb_lit(&magnetic, ", ");
            sel_merge(&magnetic, candidate, add, 6);
            jb_free(&symbol);
            jb_free(&kov);
            jb_free(&ordinary);
            jb_free(&text);
            jb_free(&m);
        }
        if (magnetic_count) {
            if (count) jb_lit(&irreps[k], ", ");
            jb_raw(&irreps[k], magnetic.data, magnetic.used);
        }
        jb_free(&magnetic);
        jb_lit(&irreps[k], "]");
        counts[k] = count + magnetic_count;
        char number[32];
        snprintf(number, sizeof number, "%lld", (long long)counts[k]);
        kv items[2] = {{"irreps", irreps[k].data, irreps[k].used},
                       {"n_irreps", number, strlen(number)}};
        jb_reset(&kp_json[k]);
        sel_merge(&kp_json[k], kpoint, items, 2);
    }
    if (!status) {
        jbuf list = {0}, new_catalog = {0};
        jb_reset(&list);
        jb_reset(&new_catalog);
        jb_lit(&list, "[");
        for (int64_t k = 0; k < nk; k++) {
            if (k) jb_lit(&list, ", ");
            jb_raw(&list, kp_json[k].data, kp_json[k].used);
        }
        jb_lit(&list, "]");
        kv cat_items[1] = {{"kpoints", list.data, list.used}};
        sel_merge(&new_catalog, catalog, cat_items, 1);
        jdoc *ndoc = jdoc_parse(new_catalog.data, new_catalog.used);
        const jvalue *rows = member(jdoc_root(ndoc), "kpoints");
        jbuf irrep_catalog = {0}, slots = {0}, sel = {0}, internal_json = {0}, mode_json = {0};
        jb_reset(&irrep_catalog);
        jb_reset(&slots);
        jb_reset(&sel);
        jb_reset(&internal_json);
        jb_reset(&mode_json);
        jb_string(&mode_json, mode_set.data, mode_set.used);
        jb_lit(&irrep_catalog, "[");
        for (int64_t k = 0; k < nk; k++) {
            const jvalue *row = jat(rows, k);
            if (k) jb_lit(&irrep_catalog, ", ");
            jb_lit(&irrep_catalog, "{");
            put_member(&irrep_catalog, row, "kslot", 1);
            put_member(&irrep_catalog, row, "kid", 0);
            jb_lit(&irrep_catalog, ", \"k_label\": ");
            jb_value(&irrep_catalog, member(row, "label"));
            static const char *keys[3] = {"kvector", "display_kvector", "dimension"};
            for (int j = 0; j < 3; j++) put_member(&irrep_catalog, row, keys[j], 0);
            jb_lit(&irrep_catalog, jtruthy(member(row, "sg_specific")) ? ", \"sg_specific\": true" : ", \"sg_specific\": false");
            put_member(&irrep_catalog, row, "star_size", 0);
            put_member(&irrep_catalog, row, "irreps", 0);
            jb_lit(&irrep_catalog, "}");
        }
        jb_lit(&irrep_catalog, "]");
        const jvalue *selected_kpoints = member(member(state, "selected"), "kpoints");
        int64_t ns = jis(selected_kpoints, J_ARRAY) ? selected_kpoints->count : 0;
        const jvalue *first_row = 0;
        jb_lit(&slots, "[");
        for (int64_t i = 0; i < ns && !status; i++) {
            const jvalue *kpoint = jat(selected_kpoints, i);
            int64_t kslot = member_int(kpoint, "kslot");
            const jvalue *row = 0;
            for (int64_t k = 0; k < nk; k++)
                if (member_int(jat(rows, k), "kslot") == kslot) row = jat(rows, k);
            if (!row) {
                status = md_fail(MDE_KEY, "%lld", (long long)kslot);
                break;
            }
            if (!i) first_row = row;
            if (i) jb_lit(&slots, ", ");
            jb_lit(&slots, "{\"slot\": ");
            int64_t slot = member_int(kpoint, "slot");
            jb_int(&slots, slot ? slot : i + 1);
            jb_lit(&slots, ", \"kpoint\": ");
            jb_value(&slots, kpoint);
            jb_lit(&slots, ", \"irreps\": ");
            jb_value(&slots, member(row, "irreps"));
            jb_lit(&slots, "}");
        }
        jb_lit(&slots, "]");
        if (!status && !first_row) status = md_fail(MDE_INDEX, "list index out of range");
        if (!status) {
            jbuf first_irreps = {0};
            jb_reset(&first_irreps);
            jb_value(&first_irreps, member(first_row, "irreps"));
            kv sel_items[2] = {{"irreps", first_irreps.data, first_irreps.used}, {"irrep_slots", slots.data, slots.used}};
            sel_merge(&sel, member(state, "selected"), sel_items, 2);
            kv internal_items[2] = {{"k_catalog", new_catalog.data, new_catalog.used}, {"mode_set", mode_json.data, mode_json.used}};
            sel_merge(&internal_json, internal, internal_items, 2);
            kv items[5] = {{"schema", "\"APOSTRUCT.irreps.v1\"", strlen("\"APOSTRUCT.irreps.v1\"")},
                           {"mode_set", mode_json.data, mode_json.used},
                           {"irrep_catalog", irrep_catalog.data, irrep_catalog.used},
                           {"selected", sel.data, sel.used},
                           {"_internal", internal_json.data, internal_json.used}};
            sel_merge(out, state, items, 5);
            jb_free(&first_irreps);
        }
        jdoc_free(ndoc);
        jb_free(&list);
        jb_free(&new_catalog);
        jb_free(&irrep_catalog);
        jb_free(&slots);
        jb_free(&sel);
        jb_free(&internal_json);
        jb_free(&mode_json);
    }
    for (int64_t k = 0; k < nk; k++) {
        jb_free(&irreps[k]);
        jb_free(&kp_json[k]);
    }
    free(irreps);
    free(kp_json);
    free(counts);
    free(displacive);
    free(magnetic_rows);
    jb_free(&mode_set);
    return status;
}
