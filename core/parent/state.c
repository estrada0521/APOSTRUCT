/* The parent state: the International setting and every atom site's Wyckoff position. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "par.h"
#include "sel.h"
#include "tables.h"

int match_inter_wyckoff_images(int sg, int setting_id, int image_count, const double *fracts, const double *source_values, int wanted_mult,
                               double tol, int *rows_out, double *params_out, int *masks_out, double *representatives_out);
int solve_wyckoff_formula(const double *base, const double *cols, const double *target, double tol, int modulo, double *params, int *mask);

static struct {
    section_view label_bc, wyckoff_label, pg_label, inter_number, inter_version, inter_origin, inter_cell, inter_axis, inter_abc,
        inter_label_full, inter_label_short, lattice_label, point_group, group_order;
    int ready;
} T;

int par_tables(void)
{
    if (T.ready) return 0;
    struct {
        const char *name;
        section_view *view;
    } list[] = {{"group.space_groups.label_bc", &T.label_bc},
                {"settings.wyckoff.label", &T.wyckoff_label},
                {"group.point_groups.label", &T.pg_label},
                {"settings.inter_settings.number", &T.inter_number},
                {"settings.inter_settings.version", &T.inter_version},
                {"settings.inter_settings.origin", &T.inter_origin},
                {"settings.inter_settings.cell", &T.inter_cell},
                {"settings.inter_settings.axis", &T.inter_axis},
                {"settings.inter_settings.abc", &T.inter_abc},
                {"settings.inter_settings.full", &T.inter_label_full},
                {"settings.inter_settings.short", &T.inter_label_short},
                {"group.lattices.label", &T.lattice_label},
                {"group.space_groups.point_group", &T.point_group},
                {"group.point_groups.order", &T.group_order}};
    if (geo_load()) return 1;
    for (size_t i = 0; i < sizeof list / sizeof list[0]; i++)
        if (section_find(list[i].name, list[i].view)) return 1;
    T.ready = 1;
    return 0;
}

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

void par_put_stripped(jbuf *b, const section_view *view, int64_t index)
{
    const char *s;
    int64_t n;
    if (section_string(view, index, &s, &n)) {
        jb_lit(b, "\"\"");
        return;
    }
    while (n && sel_is_space(*s)) {
        s++;
        n--;
    }
    while (n && sel_is_space(s[n - 1])) n--;
    jb_string(b, s, (size_t)n);
}

static void stripped(const section_view *view, int64_t index, const char **s, int64_t *n)
{
    if (section_string(view, index, s, n)) {
        *s = "";
        *n = 0;
    }
    while (*n && sel_is_space(**s)) {
        (*s)++;
        (*n)--;
    }
    while (*n && sel_is_space((*s)[*n - 1])) (*n)--;
}

/* empty for "" */
void par_preference_text(jbuf *b, int64_t id)
{
    int64_t i = id - 1;
    int64_t origin = section_int(&T.inter_origin, i), cell = section_int(&T.inter_cell, i);
    const char *abc, *axis;
    int64_t abc_n, axis_n;
    stripped(&T.inter_abc, i, &abc, &abc_n);
    stripped(&T.inter_axis, i, &axis, &axis_n);
    int parts = 0;
    if (axis_n && !(axis_n == 1 && (axis[0] == ' ' || axis[0] == '0'))) {
        jb_lit(b, "axis ");
        jb_raw(b, axis, (size_t)axis_n);
        parts++;
    }
    if (abc_n) {
        if (parts++) jb_lit(b, ", ");
        jb_lit(b, "axes ");
        jb_raw(b, abc, (size_t)abc_n);
    }
    char text[64];
    if (cell) {
        if (parts++) jb_lit(b, ", ");
        snprintf(text, sizeof text, "cell choice %lld", (long long)cell);
        jb_lit(b, text);
    }
    if (origin) {
        if (parts++) jb_lit(b, ", ");
        snprintf(text, sizeof text, "origin choice %lld", (long long)origin);
        jb_lit(b, text);
    }
}

const char *PAR_DEFAULT_PREFERENCES = "monoclinic axes a(b)c, monoclinic cell choice 1, orthorhombic axes abc, origin choice 2, hexagonal axes, "
                                      "SSG standard setting";

int64_t par_wyckoff_multiplicity(int sg, int row_id)
{
    const char *symbol;
    int64_t n;
    stripped(&T.label_bc, sg - 1, &symbol, &n);
    int64_t centering = 1;
    if (n) {
        switch (symbol[0]) {
        case 'A': case 'B': case 'C': case 'I': centering = 2; break;
        case 'F': centering = 4; break;
        case 'R': centering = 3; break;
        default: centering = 1;
        }
    }
    int64_t pg = section_int(&T.point_group, sg - 1);
    int64_t parent_order = section_int(&T.group_order, pg - 1), site_order = section_int(&T.group_order, geo_wyckoff_site_pg(row_id) - 1);
    return centering * parent_order / site_order;
}

void par_wyckoff_label(jbuf *b, int row_id) { par_put_stripped(b, &T.wyckoff_label, row_id - 1); }


void par_site_pg_label(jbuf *b, int site_pg) { par_put_stripped(b, &T.pg_label, site_pg - 1); }

static void put_frac(jbuf *b, rat v)
{
    char text[96];
    py_rat_text(v, text);
    jb_lit(b, text);
}

/* formula text (JSON) and the used parameters */
void par_formula(const rat (*v)[3], jbuf *formula, int *used)
{
    static const char names[3] = {'x', 'y', 'z'};
    for (int i = 0; i < 3; i++) used[i] = !rat_zero(v[i + 1][0]) || !rat_zero(v[i + 1][1]) || !rat_zero(v[i + 1][2]);
    jbuf t = {0};
    jb_reset(&t);
    jb_lit(&t, "(");
    for (int axis = 0; axis < 3; axis++) {
        if (axis) jb_lit(&t, ",");
        int terms = 0;
        for (int i = 0; i < 3; i++) {
            rat c = v[i + 1][axis];
            if (rat_zero(c)) continue;
            jbuf term = {0};
            jb_reset(&term);
            if (rat_eq(c, rat_int(1))) jb_raw(&term, names + i, 1);
            else if (rat_eq(c, rat_int(-1))) {
                jb_lit(&term, "-");
                jb_raw(&term, names + i, 1);
            } else {
                put_frac(&term, c);
                jb_raw(&term, names + i, 1);
            }
            if (terms && term.data[0] != '-') jb_lit(&t, "+");
            jb_raw(&t, term.data, term.used);
            jb_free(&term);
            terms++;
        }
        rat base = v[0][axis];
        if (!rat_zero(base)) {
            jbuf term = {0};
            jb_reset(&term);
            put_frac(&term, base);
            if (terms && term.data[0] != '-') jb_lit(&t, "+");
            jb_raw(&t, term.data, term.used);
            jb_free(&term);
            terms++;
        }
        if (!terms) jb_lit(&t, "0");
    }
    jb_lit(&t, ")");
    jb_string(formula, t.data, t.used);
    jb_free(&t);
}
int par_inter_vectors(int sg, int offset0, int setting_id, rat (*out)[3])
{
    rat pml[4][3], matrix[3][3], origin[3];
    geo_inter_wyckoff_vectors(setting_id, offset0, pml);
    if (geo_setting_to_cinter(sg, 1, setting_id, matrix, origin)) return md_fail(MDE_VALUE, "setting Source tables are inconsistent");
    for (int i = 0; i < 4; i++) rat_row_mul(pml[i], (const rat (*)[3])matrix, out[i]);
    for (int axis = 0; axis < 3; axis++) out[0][axis] = rat_add(out[0][axis], origin[axis]);
    return 0;
}

/* ------------------------------------------------------------ helpers */

double par_mod1(double x)
{
    double m = fmod(x, 1.0);
    if (m != 0.0) {
        if (m < 0) m += 1.0;
    } else {
        m = 0.0;
    }
    return m;
}

static double round10(double x)
{
    char text[512];
    snprintf(text, sizeof text, "%.10f", x);
    return strtod(text, 0);
}

/* float_cif_number of the three fract values: 1 when all parse */
static int site_fract(const jvalue *site, double *out)
{
    const jvalue *fract = member(site, "fract");
    if (!jtruthy(fract) || !jis(fract, J_ARRAY) || fract->count != 3) return 0;
    for (int i = 0; i < 3; i++)
        if (float_cif_number(jat(fract, i), out + i)) return 0;
    return 1;
}

static int wanted_multiplicity(const jvalue *site)
{
    const jvalue *m = member(site, "multiplicity");
    if (!jtruthy(m)) return 0;
    jbuf t = {0};
    jb_reset(&t);
    py_str_of(&t, m);
    int value = 0, found = 0;
    for (size_t i = 0; i < t.used; i++) {
        if (t.data[i] >= '0' && t.data[i] <= '9') {
            found = 1;
            value = value * 10 + (t.data[i] - '0');
        } else if (found) {
            break;
        }
    }
    jb_free(&t);
    return value;
}

typedef struct {
    double fract[3];
    char *operation;
    int index; /* -1 when absent */
} image_t;

typedef struct {
    int n;
    image_t *items;
} images_t;

static void images_free(images_t *im)
{
    for (int i = 0; i < im->n; i++) free(im->items[i].operation);
    free(im->items);
    memset(im, 0, sizeof *im);
}
static int site_images(const jvalue *site, const jvalue *ops, images_t *out)
{
    memset(out, 0, sizeof *out);
    double v[3];
    if (!site_fract(site, v)) return 0;
    int64_t nops = jis(ops, J_ARRAY) ? ops->count : 0;
    out->items = calloc((size_t)nops + 1, sizeof(image_t));
    if (!out->items) return md_fail(MDE_MEMORY, "out of memory");
    double fract[3] = {par_mod1(v[0]), par_mod1(v[1]), par_mod1(v[2])};
    image_t *first = out->items + out->n++;
    memcpy(first->fract, fract, sizeof fract);
    first->operation = (char *)md_strdup("x,y,z", 5);
    first->index = -1;
    double (*seen)[3] = calloc((size_t)nops + 1, sizeof *seen);
    if (!seen) return md_fail(MDE_MEMORY, "out of memory");
    int nseen = 0;
    for (int a = 0; a < 3; a++) seen[0][a] = round10(fract[a]);
    nseen = 1;
    for (int64_t k = 0; k < nops; k++) {
        jbuf text = {0};
        jb_reset(&text);
        py_str_of(&text, jat(ops, k));
        int rot[3][3], tran[3];
        char notation;
        if (gemmi_parse_triplet(text.data, (int64_t)text.used, rot, tran, &notation)) {
            jb_free(&text);
            continue;
        }
        double moved[3], wrapped[3], key[3];
        gemmi_apply_to_xyz(rot, tran, fract, moved);
        for (int a = 0; a < 3; a++) {
            wrapped[a] = par_mod1(moved[a]);
            key[a] = round10(wrapped[a]);
        }
        int known = 0;
        for (int s = 0; s < nseen && !known; s++) known = seen[s][0] == key[0] && seen[s][1] == key[1] && seen[s][2] == key[2];
        if (!known) {
            memcpy(seen[nseen++], key, sizeof key);
            image_t *im = out->items + out->n++;
            memcpy(im->fract, wrapped, sizeof wrapped);
            im->operation = (char *)md_strdup(text.data, (int64_t)text.used);
            im->index = (int)k;
        }
        jb_free(&text);
    }
    free(seen);
    return 0;
}

/* One Wyckoff match with its representative image. */
typedef struct {
    int found;
    int row_id, offset0;
    int64_t multiplicity;
    jbuf formula;
    int used[3];
    double params[3];
    int has_inter;
    int64_t inter_setting;
    double representative[3];
    int has_image;
    const jvalue *source_site;
    double matched[3];
    char *operation;
    int operation_index;
    int has_translation;
    int64_t translation[3];
} match_t;

static void match_free(match_t *m)
{
    jb_free(&m->formula);
    free(m->operation);
    memset(m, 0, sizeof *m);
}

static int free_count(const match_t *m) { return (m->used[0] != 0) + (m->used[1] != 0) + (m->used[2] != 0); }

static int specificity_better(const match_t *a, int a_index, const match_t *b, int b_index)
{
    int ca = 3 - free_count(a), cb = 3 - free_count(b);
    if (ca < 0) ca = 0;
    if (cb < 0) cb = 0;
    if (ca != cb) return ca > cb;
    if (a->multiplicity != b->multiplicity) return a->multiplicity < b->multiplicity;
    return a_index < b_index;
}

static void set_image(match_t *m, const jvalue *site, const image_t *im)
{
    m->has_image = 1;
    m->source_site = site;
    memcpy(m->matched, im->fract, sizeof m->matched);
    free(m->operation);
    m->operation = (char *)md_strdup(im->operation, (int64_t)strlen(im->operation));
    m->operation_index = im->index;
}

/* _match_inter_wyckoff_site_with_symmetry */
static int inter_match(int sg, const jvalue *site, const images_t *im, int setting_id, match_t *best)
{
    memset(best, 0, sizeof *best);
    double source[3];
    if (!im->n || !site_fract(site, source)) return 0;
    int n = im->n;
    double *fracts = malloc(sizeof(double) * 3 * (size_t)n), *params = malloc(sizeof(double) * 3 * (size_t)n),
           *reps = malloc(sizeof(double) * 3 * (size_t)n);
    int *rows = malloc(sizeof(int) * (size_t)n), *masks = malloc(sizeof(int) * 3 * (size_t)n);
    if (!fracts || !params || !reps || !rows || !masks) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < n; i++) memcpy(fracts + 3 * i, im->items[i].fract, sizeof(double) * 3);
    int status = match_inter_wyckoff_images(sg, setting_id, n, fracts, source, wanted_multiplicity(site), 1e-5, rows, params, masks, reps);
    if (status < 0) status = md_fail(MDE_VALUE, "Wyckoff site match failed (%d)", status);
    int start, count;
    geo_wyckoff_row_range(sg, &start, &count);
    int best_index = -1;
    for (int i = 0; !status && i < n; i++) {
        if (rows[i] < 0) continue;
        match_t m = {0};
        m.found = 1;
        m.offset0 = rows[i];
        m.row_id = start + rows[i];
        m.multiplicity = par_wyckoff_multiplicity(sg, m.row_id);
        rat v[4][3];
        status = par_inter_vectors(sg, m.offset0, setting_id, v);
        if (status) break;
        int used[3];
        jb_reset(&m.formula);
        par_formula((const rat (*)[3])v, &m.formula, used);
        for (int a = 0; a < 3; a++) {
            m.used[a] = used[a] ? (masks[3 * i + a] ? 1 : 2) : 0;
            m.params[a] = params[3 * i + a];
            m.representative[a] = reps[3 * i + a];
        }
        if (best_index >= 0 && !specificity_better(&m, i, best, best_index)) {
            match_free(&m);
            continue;
        }
        m.has_inter = 1;
        m.inter_setting = setting_id;
        set_image(&m, site, im->items + i);
        if (best_index >= 0) match_free(best);
        *best = m;
        best_index = i;
    }
    free(fracts);
    free(params);
    free(reps);
    free(rows);
    free(masks);
    return status;
}

/* the Wyckoff match in the Source frame */
static int legacy_match(int sg, const jvalue *site, const double *fract, match_t *out)
{
    memset(out, 0, sizeof *out);
    int start, count;
    geo_wyckoff_row_range(sg, &start, &count);
    int wanted = wanted_multiplicity(site);
    for (int k = 0; k < count; k++) {
        int row_id = start + k;
        int64_t mult = par_wyckoff_multiplicity(sg, row_id);
        if (wanted && mult != wanted) continue;
        rat v[4][3];
        geo_wyckoff_vectors(row_id, v);
        double base[3], cols[9];
        for (int a = 0; a < 3; a++) {
            py_float_of_fraction(v[0][a], base + a);
            for (int i = 0; i < 3; i++) py_float_of_fraction(v[i + 1][a], cols + 3 * i + a);
        }
        double params[3];
        int mask[3];
        int found = solve_wyckoff_formula(base, cols, fract, 1e-5, 1, params, mask);
        if (found < 0) return md_fail(MDE_VALUE, "Wyckoff formula solve failed (%d)", found);
        if (!found) continue;
        out->found = 1;
        out->row_id = row_id;
        out->offset0 = k;
        out->multiplicity = mult;
        int used[3];
        jb_reset(&out->formula);
        par_formula((const rat (*)[3])v, &out->formula, used);
        for (int a = 0; a < 3; a++) {
            out->used[a] = used[a] && mask[a] ? 1 : used[a] ? 2 : 0;
            out->params[a] = params[a];
        }
        return 0;
    }
    return 0;
}
static int legacy_symmetry_match(int sg, const jvalue *site, const jvalue *ops, match_t *best)
{
    memset(best, 0, sizeof *best);
    images_t im;
    MD_TRY(site_images(site, ops, &im));
    int best_index = -1, status = 0;
    for (int i = 0; i < im.n && !status; i++) {
        double fract[3];
        for (int a = 0; a < 3; a++) {
            char text[64];
            py_float_repr(im.items[i].fract[a], text);
            py_float_of_text(text, (int64_t)strlen(text), fract + a);
        }
        match_t m;
        status = legacy_match(sg, site, fract, &m);
        if (status || !m.found) {
            match_free(&m);
            continue;
        }
        if (best_index >= 0 && !specificity_better(&m, i, best, best_index)) {
            match_free(&m);
            continue;
        }
        set_image(&m, site, im.items + i);
        if (best_index >= 0) match_free(best);
        *best = m;
        best_index = i;
    }
    images_free(&im);
    if (status || best->found) return status;
    double v[3];
    if (geo_sg_lattice(sg) == 2 && site_fract(site, v)) {
        double x = par_mod1(v[0]), y = par_mod1(v[1]), z = par_mod1(v[2]);
        double swapped[3] = {x, z, y};
        match_t m;
        MD_TRY(legacy_match(sg, site, swapped, &m));
        if (m.found) {
            image_t image = {{x, z, y}, "x,z,y", -1};
            set_image(&m, site, &image);
            *best = m;
        }
    }
    return 0;
}

/* onto m */
static void canonical_mapping(const jvalue *site, const jvalue *ops, match_t *m)
{
    if (!m->has_inter) return;
    double source[3];
    if (!site_fract(site, source)) return;
    int64_t nops = jis(ops, J_ARRAY) ? ops->count : 0;
    for (int64_t k = -1; k < nops; k++) {
        jbuf text = {0};
        jb_reset(&text);
        if (k < 0) jb_lit(&text, "x,y,z");
        else py_str_of(&text, jat(ops, k));
        int rot[3][3], tran[3];
        char notation;
        if (gemmi_parse_triplet(text.data, (int64_t)text.used, rot, tran, &notation)) {
            jb_free(&text);
            continue;
        }
        double applied[3];
        gemmi_apply_to_xyz(rot, tran, source, applied);
        int ok = 1;
        int64_t integer[3];
        for (int a = 0; a < 3 && ok; a++) {
            double translation = applied[a] - m->representative[a];
            integer[a] = (int64_t)nearbyint(translation);
            ok = fabs(translation - (double)integer[a]) <= 2e-5;
        }
        if (ok) {
            m->has_image = 1;
            m->source_site = site;
            memcpy(m->matched, m->representative, sizeof m->matched);
            free(m->operation);
            m->operation = (char *)md_strdup(text.data, (int64_t)text.used);
            m->operation_index = k < 0 ? -1 : (int)k;
            m->has_translation = 1;
            memcpy(m->translation, integer, sizeof integer);
            jb_free(&text);
            return;
        }
        jb_free(&text);
    }
}

typedef struct {
    int64_t score[5];
    int id;
    match_t *matches;
} scored_t;
static int setting_score(int sg, const jvalue *sites, int drop_multiplicity, images_t *images, int setting_id, scored_t *out)
{
    int64_t n = sites->count;
    out->id = setting_id;
    memset(out->score, 0, sizeof out->score);
    out->matches = calloc((size_t)n + 1, sizeof(match_t));
    if (!out->matches) return md_fail(MDE_MEMORY, "out of memory");
    jdoc *stripped_doc = 0;
    for (int64_t i = 0; i < n; i++) {
        const jvalue *site = jat(sites, i);
        if (!jtruthy(member(site, "fract"))) continue;
        const jvalue *probe = site;
        if (drop_multiplicity && member(site, "multiplicity")) {
            jbuf t = {0};
            jb_reset(&t);
            jb_lit(&t, "{");
            int first = 1;
            for (int64_t k = 0; k < site->count; k++) {
                const jvalue *key = site->u.members + 2 * k;
                if (key->count == 12 && !memcmp(key->u.s, "multiplicity", 12)) continue;
                if (!first) jb_lit(&t, ", ");
                first = 0;
                jb_string(&t, key->u.s, (size_t)key->count);
                jb_lit(&t, ": ");
                jb_value(&t, key + 1);
            }
            jb_lit(&t, "}");
            jdoc_free(stripped_doc);
            stripped_doc = jdoc_parse(t.data, t.used);
            jb_free(&t);
            probe = jdoc_root(stripped_doc);
        }
        match_t *m = out->matches + i;
        int status = inter_match(sg, probe, images + i, setting_id, m);
        if (m->found) m->source_site = site;
        if (status) {
            jdoc_free(stripped_doc);
            return status;
        }
        if (m->found && m->operation_index < 0) out->score[3]++;
        if (m->found) {
            out->score[0]++;
            int c = 3 - free_count(m);
            out->score[1] += c > 0 ? c : 0;
        }
    }
    jdoc_free(stripped_doc);
    int default_id = geo_default_setting(sg);
    out->score[4] = setting_id == default_id;
    int same = section_int(&T.inter_cell, setting_id - 1) == section_int(&T.inter_cell, default_id - 1);
    const char *a, *b;
    int64_t an, bn;
    stripped(&T.inter_axis, setting_id - 1, &a, &an);
    stripped(&T.inter_axis, default_id - 1, &b, &bn);
    same &= an == bn && !memcmp(a, b, (size_t)an);
    stripped(&T.inter_abc, setting_id - 1, &a, &an);
    stripped(&T.inter_abc, default_id - 1, &b, &bn);
    same &= an == bn && !memcmp(a, b, (size_t)an);
    out->score[2] = same;
    return 0;
}

static void scored_free(scored_t *s, int64_t n)
{
    for (int64_t i = 0; s->matches && i < n; i++) match_free(s->matches + i);
    free(s->matches);
    s->matches = 0;
}

static int better(const scored_t *a, const scored_t *b)
{
    for (int i = 0; i < 5; i++)
        if (a->score[i] != b->score[i]) return a->score[i] > b->score[i];
    return a->id < b->id;
}

static int best_setting(int sg, const jvalue *sites, int drop, images_t *images, const int *candidates, int ncand, scored_t *best)
{
    int64_t n = sites->count;
    memset(best, 0, sizeof *best);
    for (int c = 0; c < ncand; c++) {
        scored_t s;
        int status = setting_score(sg, sites, drop, images, candidates[c], &s);
        if (status) {
            scored_free(&s, n);
            scored_free(best, n);
            return status;
        }
        if (!c || better(&s, best)) {
            if (c) scored_free(best, n);
            *best = s;
        } else {
            scored_free(&s, n);
        }
    }
    return 0;
}

typedef struct {
    int64_t v[12];
} op_key;

static int key_cmp(const void *a, const void *b)
{
    const int64_t *x = ((const op_key *)a)->v, *y = ((const op_key *)b)->v;
    for (int i = 0; i < 12; i++)
        if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
    return 0;
}

static int unique_keys(op_key *keys, int n)
{
    qsort(keys, (size_t)n, sizeof *keys, key_cmp);
    int m = 0;
    for (int i = 0; i < n; i++)
        if (!m || key_cmp(keys + i, keys + m - 1)) keys[m++] = keys[i];
    return m;
}

static int key_of(const rat (*rot)[3], const rat *tran, op_key *out)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            rat v = rat_mul(rot[i][j], rat_int(24));
            if (v.d != 1) return 1;
            out->v[3 * i + j] = (int64_t)v.n;
        }
    for (int i = 0; i < 3; i++) {
        rat v = rat_mul(rat_mod1(tran[i]), rat_int(24));
        if (v.d != 1) return 1;
        out->v[9 + i] = (int64_t)v.n;
    }
    return 0;
}

/* The full operation set of one inter setting: Source records moved from the
 * cinter frame by the setting's pml matrix and origin, with every centering
 * translation. Rotations act on columns, as CIF xyz triplets do. */
static int setting_keys(int sg, int setting_id, op_key *keys, int capacity)
{
    rat m[3][3], origin[3], inverse[3][3];
    if (geo_setting_to_cinter(sg, 1, setting_id, m, origin) || rat_inverse3((const rat (*)[3])m, inverse)) return -1;
    rat centers[64][3];
    int ncenters = 1;
    for (int k = 0; k < 3; k++) centers[0][k] = rat_int(0);
    for (int done = 0; done < ncenters; done++)
        for (int row = 0; row < 3; row++) {
            rat c[3];
            for (int k = 0; k < 3; k++) c[k] = rat_mod1(rat_add(centers[done][k], m[row][k]));
            int known = 0;
            for (int i = 0; i < ncenters && !known; i++) known = rat_eq(c[0], centers[i][0]) && rat_eq(c[1], centers[i][1]) && rat_eq(c[2], centers[i][2]);
            if (known) continue;
            if (ncenters == 64) return -1;
            memcpy(centers[ncenters++], c, sizeof c);
        }
    int64_t records[GEO_MAX_OPS][5];
    int nrec = geo_sg_records(sg, records), n = 0;
    for (int r = 0; r < nrec; r++) {
        rat point[3][3], rm[3][3], r1[3][3];
        for (int i = 0; i < 3; i++) {
            rat unit[3] = {rat_int(i == 0), rat_int(i == 1), rat_int(i == 2)};
            geo_vrot(sg, (int)records[r][4], unit, point[i]);
        }
        for (int i = 0; i < 3; i++) rat_row_mul(point[i], (const rat (*)[3])m, rm[i]);
        for (int i = 0; i < 3; i++) rat_row_mul(inverse[i], (const rat (*)[3])rm, r1[i]);
        rat t[3], t1[3], shifted[3], rot[3][3];
        for (int i = 0; i < 3; i++) t[i] = rat_make(records[r][i], records[r][3]);
        rat_row_mul(t, (const rat (*)[3])m, t1);
        rat_row_mul(origin, (const rat (*)[3])r1, shifted);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) rot[i][j] = r1[j][i];
        for (int c = 0; c < ncenters; c++) {
            rat tran[3];
            for (int i = 0; i < 3; i++) tran[i] = rat_add(rat_sub(rat_add(t1[i], origin[i]), shifted[i]), centers[c][i]);
            if (n == capacity || key_of((const rat (*)[3])rot, tran, keys + n)) return -1;
            n++;
        }
    }
    if (geo_take_error()) return -1;
    return unique_keys(keys, n);
}

/* The candidate settings whose operation set equals the CIF's, written to out; returns their number, or -1 when the
 * CIF gives no operations or they cannot be read. */
static int matching_settings(int sg, const jvalue *ops, const int *candidates, int ncand, int *out)
{
    int64_t n = jis(ops, J_ARRAY) ? ops->count : 0;
    if (!n || n > 4096) return -1;
    op_key *want = calloc((size_t)n, sizeof *want), *have = calloc(4096, sizeof *have);
    int ok = want && have, nwant = 0, matches = 0;
    for (int64_t i = 0; ok && i < n; i++) {
        jbuf text = {0};
        jb_reset(&text);
        py_str_of(&text, jat(ops, i));
        int rot[3][3], tran[3];
        char notation;
        if (text.failed || gemmi_parse_triplet(text.data, (int64_t)text.used, rot, tran, &notation)) ok = 0;
        jb_free(&text);
        if (!ok) break;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) want[nwant].v[3 * a + b] = rot[a][b];
        for (int a = 0; a < 3; a++) want[nwant].v[9 + a] = ((tran[a] % 24) + 24) % 24;
        nwant++;
    }
    if (ok) {
        nwant = unique_keys(want, nwant);
        for (int c = 0; c < ncand; c++) {
            int nhave = setting_keys(sg, candidates[c], have, 4096);
            if (nhave == nwant && !memcmp(have, want, sizeof *have * (size_t)nwant)) out[matches++] = candidates[c];
        }
    }
    free(want);
    free(have);
    return ok ? matches : -1;
}

static int select_setting(int sg, const jvalue *sites, const jvalue *ops, scored_t *chosen)
{
    int64_t n = sites->count;
    int candidates[64], ncand = 0;
    int64_t total = section_count(&T.inter_number);
    for (int64_t i = 0; i < total && ncand < 64; i++)
        if (section_int(&T.inter_number, i) == sg) candidates[ncand++] = (int)i + 1;
    if (!ncand) candidates[ncand++] = geo_default_setting(sg);
    int matched[64], nmatched = matching_settings(sg, ops, candidates, ncand, matched);
    if (nmatched == 0)
        return md_fail(MDE_VALUE, "the CIF's symmetry operations match no setting of space group %d in the Source", sg);
    if (nmatched > 0) {
        memcpy(candidates, matched, sizeof(int) * (size_t)nmatched);
        ncand = nmatched;
        const char *axis;
        int64_t length;
        if (!section_string(&T.inter_axis, candidates[0] - 1, &axis, &length) && length == 1 && axis[0] == 'R')
            return md_fail(MDE_VALUE, "rhombohedral axes are not supported; give the structure in hexagonal axes");
    }
    images_t *images = calloc((size_t)n + 1, sizeof(images_t));
    if (!images) return md_fail(MDE_MEMORY, "out of memory");
    int status = 0;
    if (nmatched > 0) {
        int matches = ncand, operation_setting = candidates[0];
        if (matches == 1) {
            for (int64_t i = 0; i < n && !status; i++) status = site_images(jat(sites, i), ops, images + i);
            scored_t s;
            if (!status) status = setting_score(sg, sites, 0, images, operation_setting, &s);
            int all = !status;
            for (int64_t i = 0; all && i < n; i++)
                if (jtruthy(member(jat(sites, i), "fract")) && !s.matches[i].found) all = 0;
            if (all) {
                *chosen = s;
                for (int64_t i = 0; i < n; i++) images_free(images + i);
                free(images);
                return 0;
            }
            if (!status) scored_free(&s, n);
            for (int64_t i = 0; i < n; i++) images_free(images + i);
        }
    }
    for (int64_t i = 0; i < n && !status; i++) status = site_images(jat(sites, i), ops, images + i);
    if (!status) status = best_setting(sg, sites, 0, images, candidates, ncand, chosen);
    int64_t site_count = 0;
    for (int64_t i = 0; i < n; i++) site_count += jtruthy(member(jat(sites, i), "fract"));
    if (!status && chosen->score[0] < site_count) {
        scored_t coordinate;
        status = best_setting(sg, sites, 1, images, candidates, ncand, &coordinate);
        if (!status && coordinate.score[0] > chosen->score[0]) {
            scored_free(chosen, n);
            *chosen = coordinate;
        } else if (!status) {
            scored_free(&coordinate, n);
        }
    }
    for (int64_t i = 0; i < n; i++) images_free(images + i);
    free(images);
    return status;
}

/* ------------------------------------------------------------ formula */

static int element_of(const jvalue *site, char *out)
{
    const jvalue *v = member(site, "type");
    if (!jtruthy(v)) v = member(site, "label");
    jbuf t = {0};
    jb_reset(&t);
    if (jtruthy(v)) py_str_of(&t, v);
    size_t i = 0;
    while (i < t.used && sel_is_space(t.data[i])) i++;
    int ok = i < t.used && t.data[i] >= 'A' && t.data[i] <= 'Z';
    if (ok) {
        out[0] = t.data[i];
        out[1] = i + 1 < t.used && t.data[i + 1] >= 'a' && t.data[i + 1] <= 'z' ? t.data[i + 1] : 0;
        out[2] = 0;
    }
    jb_free(&t);
    return ok;
}

/* 0 with *out, 1 when absent */
static int formula_number(const jvalue *v, int has_default, rat *out)
{
    jbuf t = {0};
    jb_reset(&t);
    if (jtruthy(v)) py_str_of(&t, v);
    const char *a = t.data ? t.data : "", *e = a + t.used;
    while (a < e && sel_is_space(*a)) a++;
    while (e > a && sel_is_space(e[-1])) e--;
    int status;
    if (a == e || (e - a == 1 && (*a == '.' || *a == '?'))) {
        status = !has_default;
        if (has_default) *out = rat_int(1);
    } else {
        if (e[-1] == ')') {
            const char *leftmost = 0;
            for (const char *p = e - 2; p >= a && *p != ')'; p--)
                if (*p == '(') leftmost = p;
            if (leftmost) e = leftmost;
        }
        status = py_fraction_of_text(a, e - a, out) ? 1 : 0;
    }
    jb_free(&t);
    return status;
}

/* the formula over the output sites (wyckoff_multiplicity set); 0 when there is none */
static int formula_from_sites(const jvalue *sites, jbuf *out)
{
    char elements[64][3];
    rat totals[64];
    int n = 0;
    for (int64_t i = 0; i < sites->count; i++) {
        const jvalue *site = jat(sites, i);
        char element[3];
        rat multiplicity, occupancy;
        if (!element_of(site, element)) return 0;
        const jvalue *m = member(site, "wyckoff_multiplicity");
        if (!jtruthy(m)) m = member(site, "multiplicity");
        int bad = formula_number(m, 0, &multiplicity);
        if (bad || formula_number(member(site, "occupancy"), 1, &occupancy)) return 0;
        int k = 0;
        while (k < n && strcmp(elements[k], element)) k++;
        if (k == n) {
            if (n == 64) return 0;
            memcpy(elements[n], element, 3);
            totals[n++] = rat_int(0);
        }
        totals[k] = rat_add(totals[k], rat_mul(multiplicity, occupancy));
    }
    if (!n) return 0;
    for (int k = 0; k < n; k++)
        if (totals[k].n <= 0) return 0;
    i128 den = 1;
    for (int k = 0; k < n; k++) den = den / rat_gcd(den, totals[k].d) * totals[k].d;
    i128 counts[64], g = 0;
    for (int k = 0; k < n; k++) {
        counts[k] = totals[k].n * (den / totals[k].d);
        g = rat_gcd(g, counts[k]);
    }
    jbuf t = {0};
    jb_reset(&t);
    for (int k = 0; k < n; k++) {
        if (k) jb_lit(&t, " ");
        jb_lit(&t, elements[k]);
        i128 c = counts[k] / g;
        if (c != 1) {
            char number[64];
            snprintf(number, sizeof number, "%lld", (long long)c);
            jb_lit(&t, number);
        }
    }
    jb_string(out, t.data, t.used);
    jb_free(&t);
    return 1;
}

/* ------------------------------------------------------------ state */

static void put_match_site(jbuf *b, int sg, const jvalue *site, const match_t *m)
{
    static const char *popped[7] = {"wyckoff", "wyckoff_row_id", "site_pg", "wyckoff_formula", "wyckoff_params", "wyckoff_multiplicity",
                                    "wyckoff_mapping_error"};
    jbuf base = {0};
    jb_reset(&base);
    jb_lit(&base, "{");
    int first = 1;
    for (int64_t k = 0; k < site->count; k++) {
        const jvalue *key = site->u.members + 2 * k;
        int drop = 0;
        for (int p = 0; p < 7 && !drop; p++) drop = (size_t)key->count == strlen(popped[p]) && !memcmp(key->u.s, popped[p], (size_t)key->count);
        if (drop) continue;
        if (!first) jb_lit(&base, ", ");
        first = 0;
        jb_string(&base, key->u.s, (size_t)key->count);
        jb_lit(&base, ": ");
        jb_value(&base, key + 1);
    }
    jb_lit(&base, "}");
    jdoc *doc = jdoc_parse(base.data, base.used);
    jbuf texts[12] = {{0}};
    for (int i = 0; i < 12; i++) jb_reset(&texts[i]);
    kv items[12];
    int n = 0;
    if (m && m->found) {
        par_wyckoff_label(&texts[0], m->row_id);
        items[n++] = (kv){"wyckoff", texts[0].data, texts[0].used};
        jb_int(&texts[1], m->row_id);
        items[n++] = (kv){"wyckoff_row_id", texts[1].data, texts[1].used};
        par_site_pg_label(&texts[2], geo_wyckoff_site_pg(m->row_id));
        items[n++] = (kv){"site_pg", texts[2].data, texts[2].used};
        items[n++] = (kv){"wyckoff_formula", m->formula.data, m->formula.used};
        jb_lit(&texts[3], "{");
        static const char *names[3] = {"\"x\": ", "\"y\": ", "\"z\": "};
        int any = 0;
        for (int a = 0; a < 3; a++)
            if (m->used[a] == 1) {
                if (any++) jb_lit(&texts[3], ", ");
                jb_lit(&texts[3], names[a]);
                jb_float(&texts[3], m->params[a]);
            }
        jb_lit(&texts[3], "}");
        items[n++] = (kv){"wyckoff_params", texts[3].data, texts[3].used};
        jb_int(&texts[4], m->multiplicity);
        items[n++] = (kv){"wyckoff_multiplicity", texts[4].data, texts[4].used};
        if (m->has_image) {
            const jvalue *fract = member(m->source_site, "fract");
            if (jtruthy(fract)) jb_value(&texts[5], fract);
            else jb_lit(&texts[5], "[]");
            items[n++] = (kv){"representative_source_fract", texts[5].data, texts[5].used};
            jb_lit(&texts[6], "[");
            for (int a = 0; a < 3; a++) {
                if (a) jb_lit(&texts[6], ", ");
                jb_float(&texts[6], m->matched[a]);
            }
            jb_lit(&texts[6], "]");
            items[n++] = (kv){"representative_matched_fract", texts[6].data, texts[6].used};
        }
    } else {
        double v[3];
        const jvalue *fract = member(site, "fract");
        if (!jis(fract, J_ARRAY) || fract->count != 3) {
            jb_cstring(&texts[10], "fractional coordinates are missing");
        } else {
            int invalid = -1;
            for (int a = 0; a < 3 && invalid < 0; a++)
                if (float_cif_number(jat(fract, a), v + a)) invalid = a;
            jbuf message = {0};
            jb_reset(&message);
            if (invalid >= 0) {
                jb_lit(&message, "unsupported fractional coordinate value ");
                py_repr_of(&message, jat(fract, invalid));
            } else {
                char text[96];
                snprintf(text, sizeof text, "coordinates do not match a Wyckoff position in any setting of SG%d", sg);
                jb_lit(&message, text);
            }
            jb_string(&texts[10], message.data, message.used);
            jb_free(&message);
        }
        items[n++] = (kv){"wyckoff_mapping_error", texts[10].data, texts[10].used};
    }
    sel_merge(b, jdoc_root(doc), items, n);
    jdoc_free(doc);
    for (int i = 0; i < 12; i++) jb_free(&texts[i]);
    jb_free(&base);
}

static int64_t member_int(const jvalue *object, const char *key)
{
    int64_t v = 0;
    const jvalue *m = member(object, key);
    if (jtruthy(m)) py_int_of(m, &v);
    return v;
}

int par_parent_state(const jvalue *cif_info, jbuf *out)
{
    if (par_tables()) return md_fail(MDE_VALUE, "parent Source tables are inconsistent");
    int sg = (int)member_int(member(cif_info, "parent"), "number");
    const jvalue *sites = member(cif_info, "atom_sites"), *ops = member(cif_info, "symmetry_operations");
    jdoc *empty = 0;
    if (!jis(sites, J_ARRAY)) {
        empty = jdoc_parse("[]", 2);
        sites = jdoc_root(empty);
    }
    int64_t n = sites->count;
    scored_t chosen;
    int status = select_setting(sg, sites, ops, &chosen);
    if (status) {
        jdoc_free(empty);
        return status;
    }
    match_t *final = calloc((size_t)n + 1, sizeof(match_t));
    if (!final) return md_fail(MDE_MEMORY, "out of memory");
    for (int64_t i = 0; i < n && !status; i++) {
        const jvalue *site = jat(sites, i);
        if (chosen.matches[i].found) {
            final[i] = chosen.matches[i];
            memset(&chosen.matches[i], 0, sizeof(match_t));
            final[i].source_site = site;
        } else {
            status = legacy_symmetry_match(sg, site, ops, &final[i]);
        }
        if (!status && final[i].found) canonical_mapping(site, ops, &final[i]);
    }
    if (status) {
        for (int64_t i = 0; i < n; i++) match_free(final + i);
        free(final);
        scored_free(&chosen, n);
        jdoc_free(empty);
        return status;
    }
    jbuf site_list = {0};
    jb_reset(&site_list);
    jb_lit(&site_list, "[");
    for (int64_t i = 0; i < n; i++) {
        if (i) jb_lit(&site_list, ", ");
        jbuf one = {0};
        jb_reset(&one);
        put_match_site(&one, sg, jat(sites, i), final + i);
        jb_raw(&site_list, one.data, one.used);
        jb_free(&one);
    }
    jb_lit(&site_list, "]");
    jdoc *sdoc = jdoc_parse(site_list.data, site_list.used);
    jbuf formula = {0};
    jb_reset(&formula);
    kv input_items[2];
    int ninput = 0;
    if (jis(cif_info, J_OBJECT) && jget(cif_info, "atom_sites")) input_items[ninput++] = (kv){"atom_sites", site_list.data, site_list.used};
    jbuf formula_text = {0};
    jb_reset(&formula_text);
    if (jtruthy(member(cif_info, "formula"))) py_str_of(&formula_text, member(cif_info, "formula"));
    int blank = 1;
    for (size_t k = 0; k < formula_text.used && blank; k++) blank = sel_is_space(formula_text.data[k]);
    if (blank) {
        if (!formula_from_sites(jdoc_root(sdoc), &formula)) jb_lit(&formula, "null");
        input_items[ninput++] = (kv){"formula", formula.data, formula.used};
    }
    jbuf input = {0}, internal = {0}, space_group = {0}, rows = {0};
    jb_reset(&input);
    jb_reset(&internal);
    jb_reset(&space_group);
    jb_reset(&rows);
    sel_merge(&input, cif_info, input_items, ninput);
    int parent_id = chosen.id;
    int all_matched = 1;
    for (int64_t i = 0; i < n; i++) all_matched &= final[i].found && final[i].row_id;
    jb_lit(&rows, "[");
    if (all_matched) {
        int last = -1, first = 1;
        for (;;) {
            int next = -1;
            for (int64_t i = 0; i < n; i++)
                if (final[i].row_id > last && (next < 0 || final[i].row_id < next)) next = final[i].row_id;
            if (next < 0) break;
            if (!first) jb_lit(&rows, ", ");
            first = 0;
            jb_int(&rows, next);
            last = next;
        }
    }
    jb_lit(&rows, "]");
    int lattice = geo_sg_lattice(sg);
    jb_lit(out, "{\"schema\": \"APOSTRUCT.parent.v1\", \"input\": ");
    jb_raw(out, input.data, input.used);
    jb_lit(out, ", \"space_group\": {\"number\": ");
    jb_int(out, sg);
    jb_lit(out, ", \"symbol\": ");
    const jvalue *symbol = member(member(cif_info, "parent"), "symbol");
    if (symbol) jb_value(out, symbol);
    else jb_lit(out, "\"\"");
    jb_lit(out, ", \"lattice\": ");
    jb_int(out, lattice);
    jb_lit(out, ", \"lattice_type\": ");
    par_put_stripped(out, &T.lattice_label, lattice - 1);
    jb_lit(out, "}");
    jb_lit(out, ", \"_internal\": {\"parent_inter_setting_id\": ");
    jb_int(out, parent_id);
    jb_lit(out, ", \"displacive_row_ids\": ");
    jb_raw(out, rows.data, rows.used);
    jb_lit(out, "}}");
    for (int64_t i = 0; i < n; i++) match_free(final + i);
    free(final);
    scored_free(&chosen, n);
    jdoc_free(sdoc);
    jdoc_free(empty);
    jb_free(&site_list);
    jb_free(&formula);
    jb_free(&formula_text);
    jb_free(&input);
    jb_free(&internal);
    jb_free(&space_group);
    jb_free(&rows);
    return 0;
}
