/* What the OPD rows share: k parameters, OPD text, the displayed setting and basis of a subgroup, and the active k
 * vectors. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "sel.h"
#include "irreps.h"
#include "tables.h"

int opd_basis_presentation(int subgroup, const int64_t *basis, const double *cell, int parametric, int coupled, int magnetic_group,
                           int64_t *out);
const char *k_numeric_star_json(int parent_sg, int gid, int has_visible, const int64_t *visible_words, int has_source_visible,
                                const int64_t *source_visible_words);
const char *k_display_star_json(int sg, int kslot, int dim);
int k_coordinate_map_words(int sg, int kslot, int64_t *out);

static struct {
    section_view point_group, group_order, irr_sg, irr_k, k_dim;
    int ready;
} T;

int sel_tables(void)
{
    if (T.ready) return 0;
    struct {
        const char *name;
        section_view *view;
    } list[] = {{"group.space_groups.point_group", &T.point_group},
                {"group.point_groups.order", &T.group_order},
                {"irreps.little.sg", &T.irr_sg},
                {"irreps.little.k", &T.irr_k},
                {"k.slots.dim", &T.k_dim}};
    if (geo_load()) return 1;
    for (size_t i = 0; i < sizeof list / sizeof list[0]; i++)
        if (section_find(list[i].name, list[i].view)) return 1;
    T.ready = 1;
    return 0;
}

/* ------------------------------------------------------------ k params */

void sel_params_of(const jvalue *object, sel_params *out)
{
    out->count = 0;
    if (!jis(object, J_OBJECT)) return;
    for (int64_t i = 0; i < object->count && out->count < 16; i++) {
        const jvalue *key = object->u.members + 2 * i, *value = key + 1;
        if (!jis(value, J_STRING)) continue;
        out->names[out->count] = key->u.s;
        out->values[out->count] = value->u.s;
        out->count++;
    }
}

const char *sel_param(const sel_params *p, const char *name)
{
    const char *found = 0;
    for (int i = 0; i < p->count; i++)
        if (!strcmp(p->names[i], name)) found = p->values[i];
    return found;
}

/* ------------------------------------------------------------ formatting */

void sel_opd_from_rows(jbuf *b, char ***rows, int nrows, const int *widths)
{
    if (!nrows) return;
    int width = 0;
    for (int r = 0; r < nrows; r++)
        if (widths[r] > width) width = widths[r];
    jb_lit(b, "(");
    for (int col = 0; col < width; col++) {
        if (col) jb_lit(b, ",");
        int first = 1;
        for (int r = 0; r < nrows; r++) {
            const char *coef = col < widths[r] ? rows[r][col] : "0";
            if (!strcmp(coef, "0")) continue;
            int negative = coef[0] == '-';
            const char *mag = negative ? coef + 1 : coef;
            char param = (char)('a' + r);
            if (negative) jb_lit(b, "-");
            else if (!first) jb_lit(b, "+");
            if (strcmp(mag, "1")) jb_lit(b, mag);
            jb_raw(b, &param, 1);
            first = 0;
        }
        if (first) jb_lit(b, "0");
    }
    jb_lit(b, ")");
}

int sel_is_space(char c) { return c == ' ' || (c >= '\t' && c <= '\r'); }
void sel_display_opd(jbuf *b, const char *opd)
{
    const char *s = opd, *e = opd + strlen(opd);
    while (s < e && sel_is_space(*s)) s++;
    while (e > s && sel_is_space(e[-1])) e--;
    jbuf text = {0};
    jb_reset(&text);
    if (e - s >= 2 && *s == '(' && e[-1] == ')') {
        const char *body = s + 1, *end = e - 1;
        int parts = 1;
        for (const char *c = body; c < end; c++) parts += *c == ',';
        jb_lit(&text, "(");
        int index = 0;
        const char *a = body;
        for (const char *c = body;; c++) {
            if (c < end && *c != ',') continue;
            const char *x = a, *y = c;
            while (x < y && sel_is_space(*x)) x++;
            while (y > x && sel_is_space(y[-1])) y--;
            if (index) jb_lit(&text, parts >= 4 && parts % 2 == 0 && index % 2 == 0 ? ";" : ",");
            if (parts >= 4 && parts % 2 == 0) jb_raw(&text, x, (size_t)(y - x));
            else jb_raw(&text, a, (size_t)(c - a));
            index++;
            if (c >= end) break;
            a = c + 1;
        }
        jb_lit(&text, ")");
    } else {
        jb_raw(&text, s, (size_t)(e - s));
    }
    /* re.sub(r"-?\d+\.\d+(?=[A-Za-z])", round_coeff, text) */
    const char *t = text.data;
    size_t n = text.used, i = 0;
    while (i < n) {
        size_t start = i;
        if (t[i] == '-' && i + 1 < n && t[i + 1] >= '0' && t[i + 1] <= '9') i++;
        if (i < n && t[i] >= '0' && t[i] <= '9') {
            size_t d = i;
            while (d < n && t[d] >= '0' && t[d] <= '9') d++;
            if (d < n && t[d] == '.' && d + 1 < n && t[d + 1] >= '0' && t[d + 1] <= '9') {
                size_t f = d + 1;
                while (f < n && t[f] >= '0' && t[f] <= '9') f++;
                size_t g = f;
                /* backtrack the fraction digits until a letter follows */
                while (g > d + 1 && !(g < n && ((t[g] >= 'a' && t[g] <= 'z') || (t[g] >= 'A' && t[g] <= 'Z')))) g--;
                if (g < n && ((t[g] >= 'a' && t[g] <= 'z') || (t[g] >= 'A' && t[g] <= 'Z')) && g > d + 1) {
                    char number[128];
                    size_t length = g - start < sizeof number - 1 ? g - start : sizeof number - 1;
                    memcpy(number, t + start, length);
                    number[length] = 0;
                    double value = strtod(number, 0);
                    char rounded[64];
                    double nearest = nearbyint(value);
                    if (fabs(value - nearest) < 1e-8) snprintf(rounded, sizeof rounded, "%lld", (long long)nearest);
                    else snprintf(rounded, sizeof rounded, "%.3f", value);
                    jb_lit(b, strcmp(rounded, "-0") ? rounded : "0");
                    i = g;
                    continue;
                }
            }
            jb_raw(b, t + start, d - start);
            i = d;
            continue;
        }
        jb_raw(b, t + start, i - start + 1);
        i = start + 1;
    }
    jb_free(&text);
}

/* ------------------------------------------------------------ settings */

char *sel_default_symbol(int sg) { return md_default_symbol(sg); }

int64_t sel_point_group_order(int sg)
{
    int64_t pg = section_int(&T.point_group, sg - 1);
    return section_int(&T.group_order, pg - 1);
}

static void mat4_mul(const int64_t *a, const int64_t *b, int64_t *out)
{
    int64_t r[16];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            int64_t s = 0;
            for (int k = 0; k < 4; k++) s += a[4 * i + k] * b[4 * k + j];
            r[4 * i + j] = s;
        }
    memcpy(out, r, sizeof r);
}

/* last . middle . first */
static void trans4(const int64_t *first, const int64_t *middle, const int64_t *last, int64_t *out)
{
    int64_t inner[16];
    mat4_mul(middle, first, inner);
    mat4_mul(last, inner, out);
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

static int64_t floor_div(int64_t a, int64_t b)
{
    int64_t q = a / b;
    if ((a % b) && ((a < 0) != (b < 0))) q--;
    return q;
}

int sel_change_setting(int parent_sg, int subgroup, const int64_t *basis, const int64_t *origin, int has_parent_setting,
                       int64_t parent_setting, int has_subgroup_setting, int64_t subgroup_setting, sel_setting *out)
{
    int64_t den = origin[3] ? origin[3] : 1;
    int64_t m[16] = {basis[0] * den, basis[1] * den, basis[2] * den, 0, basis[3] * den, basis[4] * den, basis[5] * den, 0,
                     basis[6] * den, basis[7] * den, basis[8] * den, 0, origin[0], origin[1], origin[2], den};
    int64_t a[16], c[16];
    geo_cinter_base(parent_sg, 1, a);
    geo_cinter_base(subgroup, 0, c);
    trans4(a, m, c, m);
    geo_inter_block(has_parent_setting ? (int)parent_setting : geo_default_setting(parent_sg), 0, a);
    geo_inter_block(has_subgroup_setting ? (int)subgroup_setting : geo_default_setting(subgroup), 1, c);
    trans4(a, m, c, m);
    int64_t d = m[15];
    if (!d) return md_fail(MDE_VALUE, "zero denominator in subgroup basis");
    int64_t sign = d < 0 ? -1 : 1, values[9];
    for (int i = 0; i < 9; i++) values[i] = sign * m[4 * (i / 3) + i % 3];
    int64_t rden = sign * d, g = rden < 0 ? -rden : rden;
    for (int i = 0; i < 9; i++) g = gcd64(g, values[i]);
    if (!g) g = 1;
    for (int i = 0; i < 9; i++) out->basis[i / 3][i % 3] = g == 1 ? values[i] : floor_div(values[i], g);
    out->denominator = floor_div(rden, g);
    int64_t onum[3] = {m[12], m[13], m[14]}, oden = m[15];
    if (oden < 0) {
        for (int i = 0; i < 3; i++) onum[i] = -onum[i];
        oden = -oden;
    }
    int64_t og = oden;
    for (int i = 0; i < 3; i++) og = gcd64(og, onum[i]);
    if (!og) og = 1;
    for (int i = 0; i < 3; i++) out->origin[i] = floor_div(onum[i], og);
    out->origin[3] = floor_div(oden, og);
    const int64_t *r0 = out->basis[0], *r1 = out->basis[1], *r2 = out->basis[2];
    int64_t det = r0[0] * (r1[1] * r2[2] - r1[2] * r2[1]) - r0[1] * (r1[0] * r2[2] - r1[2] * r2[0]) + r0[2] * (r1[0] * r2[1] - r1[1] * r2[0]);
    if (det < 0) det = -det;
    char *symbol = sel_default_symbol(subgroup);
    int64_t centering = 1;
    if (symbol && symbol[0]) {
        switch (symbol[0]) {
        case 'A': case 'B': case 'C': case 'I': centering = 2; break;
        case 'F': centering = 4; break;
        case 'R': centering = 3; break;
        default: centering = 1;
        }
    }
    free(symbol);
    int64_t size = det % centering == 0 ? (det / centering > 1 ? det / centering : 1) : det;
    int64_t rd = out->denominator;
    if (rd && rd != 1 && size % (rd < 0 ? -rd : rd) == 0) size /= rd < 0 ? -rd : rd;
    out->size = size;
    return 0;
}

/* ------------------------------------------------------------ directions */

typedef struct {
    int64_t id;
    int valid;
    const char *label;
    int64_t label_length;
    int64_t dimension, subduction;
    int nrows;
    char ***rows;
    int *widths;
} direction;

/* ------------------------------------------------------------ active k vectors */

void strings_free(strings *s)
{
    for (int i = 0; i < s->count; i++) free(s->items[i]);
    free(s->items);
    memset(s, 0, sizeof *s);
}

static int strings_of_json(const char *json, strings *out)
{
    memset(out, 0, sizeof *out);
    if (!json) return md_fail(MDE_VALUE, "k-vector Source tables are inconsistent");
    jdoc *doc = jdoc_parse(json, strlen(json));
    if (!doc) return md_fail(MDE_VALUE, "k-vector Source tables are inconsistent");
    const jvalue *root = jdoc_root(doc);
    out->items = calloc((size_t)root->count + 1, sizeof(char *));
    for (int64_t i = 0; out->items && i < root->count; i++) out->items[out->count++] = (char *)md_strdup(jat(root, i)->u.s, jat(root, i)->count);
    jdoc_free(doc);
    return out->items ? 0 : md_fail(MDE_MEMORY, "out of memory");
}

/* 0 none, 1 with words */
int sel_visible_k_vector(int sg, int kslot, const sel_params *p, int64_t *words)
{
    int64_t map[29];
    int status = k_coordinate_map_words(sg, kslot, map);
    if (status) return -(status == 3 ? MDE_ARITH : MDE_VALUE);
    int count = (int)map[1];
    rat values[3];
    for (int i = 0; i < count; i++) {
        char name[2] = {(char)map[2 + i], 0};
        const char *raw = sel_param(p, name);
        if (!raw) raw = sel_param(p, "_single");
        if (!raw) return 0;
        const char *a = raw, *e = raw + strlen(raw);
        while (a < e && sel_is_space(*a)) a++;
        while (e > a && sel_is_space(e[-1])) e--;
        if (a == e) return 0;
        if (py_fraction_of_text(a, e - a, &values[i])) return -MDE_VALUE;
    }
    rat_overflow = 0;
    for (int axis = 0; axis < 3; axis++) {
        rat v = rat_make(map[5 + 2 * axis], map[6 + 2 * axis]);
        for (int i = 0; i < count; i++) v = rat_add(v, rat_mul(values[i], rat_make(map[11 + 6 * i + 2 * axis], map[12 + 6 * i + 2 * axis])));
        if (!rat_fits64(v)) return -MDE_ARITH;
        words[2 * axis] = (int64_t)v.n;
        words[2 * axis + 1] = (int64_t)v.d;
    }
    return rat_overflow ? -MDE_ARITH : 1;
}

/* 0 with *out, 1 when absent */
static int evaluate_k_component(const char *text, int64_t length, const sel_params *p, rat *out)
{
    char *s = malloc((size_t)length + 1);
    if (!s) return 1;
    int64_t n = 0;
    for (int64_t i = 0; i < length; i++)
        if (text[i] != ' ') s[n++] = text[i];
    s[n] = 0;
    rat total = rat_int(0);
    int terms = 0, status = 0;
    for (int64_t i = 0; i < n && !status;) {
        int64_t start = i;
        if (s[i] == '+' || s[i] == '-') i++;
        while (i < n && s[i] != '+' && s[i] != '-') i++;
        if (i == start + 1 && (s[start] == '+' || s[start] == '-')) {
            continue;
        }
        terms++;
        int sign = s[start] == '-' ? -1 : 1;
        int64_t body = s[start] == '+' || s[start] == '-' ? start + 1 : start;
        int64_t blen = i - body;
        char last = blen ? s[body + blen - 1] : 0;
        int is_param = 0;
        if (blen && (last == 'a' || last == 'b' || last == 'g')) {
            is_param = 1;
            int64_t coef_len = blen - 1, slash = -1;
            for (int64_t k = 0; k < coef_len; k++) {
                char c = s[body + k];
                if (c == '/' && slash < 0 && k > 0) slash = k;
                else if (c < '0' || c > '9') is_param = 0;
            }
            if (slash >= 0 && slash == coef_len - 1) is_param = 0;
            if (is_param) {
                char name[2] = {last, 0};
                const char *raw = sel_param(p, name);
                if (!raw) raw = sel_param(p, "_single");
                if (!raw) {
                    status = 1;
                    break;
                }
                const char *a = raw, *e = raw + strlen(raw);
                while (a < e && sel_is_space(*a)) a++;
                while (e > a && sel_is_space(e[-1])) e--;
                rat value, coef = rat_int(1);
                if (a == e || py_fraction_of_text(a, e - a, &value)) {
                    status = a == e ? 1 : 2;
                    break;
                }
                if (coef_len && py_fraction_of_text(s + body, coef_len, &coef)) {
                    status = 2;
                    break;
                }
                total = rat_add(total, rat_mul(rat_int(sign), rat_mul(coef, value)));
                continue;
            }
        }
        rat value;
        if (py_fraction_of_text(s + body, blen, &value)) {
            status = 1;
            break;
        }
        total = rat_add(total, rat_mul(rat_int(sign), value));
    }
    free(s);
    if (status == 2) return 2;
    if (status || !terms) return 1;
    *out = total;
    return 0;
}

/* the k vector with the parameters substituted, whitespace removed */
static char *substitute_k_vector(const char *vector, const sel_params *given)
{
    sel_params p = *given;
    if (!p.count) return (char *)md_strdup(vector, (int64_t)strlen(vector));
    if (!sel_param(&p, "_single")) {
        const char *values[3];
        int n = 0;
        static const char *keys[3] = {"a", "b", "g"};
        for (int i = 0; i < 3; i++) {
            const char *raw = sel_param(&p, keys[i]);
            if (!raw) continue;
            const char *a = raw, *e = raw + strlen(raw);
            while (a < e && sel_is_space(*a)) a++;
            while (e > a && sel_is_space(e[-1])) e--;
            if (a == e) continue;
            values[n++] = raw;
        }
        int same = n > 0;
        for (int i = 1; i < n; i++) {
            const char *x = values[0], *y = values[i];
            while (sel_is_space(*x)) x++;
            while (sel_is_space(*y)) y++;
            size_t lx = strlen(x), ly = strlen(y);
            while (lx && sel_is_space(x[lx - 1])) lx--;
            while (ly && sel_is_space(y[ly - 1])) ly--;
            same &= lx == ly && !memcmp(x, y, lx);
        }
        if (same && p.count < 16) {
            const char *x = values[0];
            p.names[p.count] = "_single";
            p.values[p.count] = x;
            p.count++;
        }
    }
    const char *s = vector, *e = vector + strlen(vector);
    while (s < e && sel_is_space(*s)) s++;
    while (e > s && sel_is_space(e[-1])) e--;
    jbuf out = {0};
    jb_reset(&out);
    if (!(e - s >= 2 && *s == '(' && e[-1] == ')')) {
        jb_raw(&out, s, (size_t)(e - s));
        return out.data;
    }
    jb_lit(&out, "(");
    const char *a = s + 1;
    int index = 0, failed = 0;
    for (const char *c = s + 1;; c++) {
        if (c < e - 1 && *c != ',') continue;
        rat value;
        if (evaluate_k_component(a, c - a, &p, &value)) {
            failed = 1;
            break;
        }
        if (index) jb_lit(&out, ",");
        char number[96];
        py_rat_text(value, number);
        jb_lit(&out, number);
        index++;
        if (c >= e - 1) break;
        a = c + 1;
    }
    if (failed) {
        jb_reset(&out);
        jb_raw(&out, s, (size_t)(e - s));
        return out.data;
    }
    jb_lit(&out, ")");
    return out.data;
}

static char *compact(const char *text)
{
    size_t n = strlen(text);
    char *out = malloc(n + 1);
    if (!out) return 0;
    size_t k = 0;
    for (size_t i = 0; i < n; i++)
        if (!sel_is_space(text[i])) out[k++] = text[i];
    out[k] = 0;
    return out;
}

static int nonzero_opd_value(const char *text)
{
    double value;
    if (!py_float_of_text(text, (int64_t)strlen(text), &value)) return fabs(value) > 1e-9;
    const char *s = text, *e = text + strlen(text);
    while (s < e && sel_is_space(*s)) s++;
    while (e > s && sel_is_space(e[-1])) e--;
    static const char *zeros[5] = {"", "0", "0.0", "-0", "-0.0"};
    for (int i = 0; i < 5; i++)
        if ((size_t)(e - s) == strlen(zeros[i]) && !memcmp(s, zeros[i], (size_t)(e - s))) return 0;
    return 1;
}

static int string_nonzero(const void *ctx, int row, int col) { return nonzero_opd_value(((char ***)ctx)[row][col]); }

int sel_active_k_vectors(int parent_sg, int64_t gid, int64_t full_dim, char ***rows, int nrows, const int *widths, const sel_params *p,
                         strings *out)
{
    return sel_active_k_vectors_by(parent_sg, gid, full_dim, nrows, widths, string_nonzero, rows, p, out);
}

int sel_active_k_vectors_by(int parent_sg, int64_t gid, int64_t full_dim, int nrows, const int *widths,
                            int (*nonzero)(const void *ctx, int row, int col), const void *ctx, const sel_params *p, strings *out)
{
    memset(out, 0, sizeof *out);
    int kslot = (int)section_int(&T.irr_k, gid - 1);
    int lattice = geo_sg_lattice(parent_sg);
    int dim = (int)section_int(&T.k_dim, (int64_t)(lattice - 1) * 27 + kslot - 1);
    strings star = {0};
    if (p->count && dim > 0) {
        int gid_sg = (int)section_int(&T.irr_sg, gid - 1);
        int64_t visible[6], source_visible[6];
        int has = sel_visible_k_vector(parent_sg, kslot, p, visible);
        if (has < 0) return md_fail(-has, "k-vector Source tables are inconsistent");
        int has_source = sel_visible_k_vector(gid_sg, kslot, p, source_visible);
        if (has_source < 0) return md_fail(-has_source, "k-vector Source tables are inconsistent");
        MD_TRY(strings_of_json(k_numeric_star_json(parent_sg, (int)gid, has, visible, has_source, source_visible), &star));
    }
    if (!star.count) {
        strings_free(&star);
        strings raw;
        MD_TRY(strings_of_json(k_display_star_json(parent_sg, kslot, dim), &raw));
        star.items = calloc((size_t)raw.count + 1, sizeof(char *));
        if (!star.items) return md_fail(MDE_MEMORY, "out of memory");
        for (int i = 0; i < raw.count; i++) {
            char *substituted = substitute_k_vector(raw.items[i], p);
            star.items[star.count++] = substituted ? compact(substituted) : 0;
            free(substituted);
        }
        strings_free(&raw);
    }
    if (!star.count || full_dim <= 0 || full_dim % star.count) {
        strings_free(&star);
        return 0;
    }
    int64_t arm = full_dim / star.count;
    out->items = calloc((size_t)star.count + 1, sizeof(char *));
    if (!out->items) return md_fail(MDE_MEMORY, "out of memory");
    for (int a = 0; a < star.count; a++) {
        int active = 0;
        for (int r = 0; r < nrows && !active; r++)
            for (int64_t c = a * arm; c < (a + 1) * arm && c < widths[r] && !active; c++) active = nonzero(ctx, r, (int)c);
        if (active) {
            out->items[out->count++] = star.items[a];
            star.items[a] = 0;
        }
    }
    strings_free(&star);
    return 0;
}

/* ------------------------------------------------------------ isotropy rows */

typedef struct {
    int64_t row_id, direction_id, subgroup, free, size, s;
    int has_index;
    int64_t index;
} iso_key;

static void put_frac_text(jbuf *b, rat v)
{
    char text[96];
    py_rat_text(v, text);
    jb_lit(b, text);
}

int sel_present_basis(int subgroup, const sel_setting *displayed, int has_cell, const double *cell, int parametric, int coupled,
                      int magnetic_group, rat (*rows)[3])
{
    int64_t words[36], out[36];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            rat v = rat_make(displayed->basis[r][c], displayed->denominator ? displayed->denominator : 1);
            write_i128(words + 12 * r + 4 * c, v.n);
            write_i128(words + 12 * r + 4 * c + 2, v.d);
        }
    int status = opd_basis_presentation(subgroup, words, has_cell ? cell : 0, parametric, coupled, magnetic_group, out);
    if (status == 1) return md_fail(MDE_ARITH, "OPD basis presentation value out of range");
    if (status == 2) return md_fail(MDE_ZERO, "float division by zero");
    if (status == 3) return md_fail(MDE_VALUE, "OPD basis presentation Source tables are inconsistent");
    if (status) return md_fail(MDE_RUNTIME, "OPD basis presentation failed (%d)", status);
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) rows[r][c] = (rat){read_i128(out + 12 * r + 4 * c), read_i128(out + 12 * r + 4 * c + 2)};
    return 0;
}

void sel_put_basis(jbuf *b, const rat (*rows)[3])
{
    jb_lit(b, "[");
    for (int r = 0; r < 3; r++) {
        if (r) jb_lit(b, ", ");
        jb_lit(b, "[");
        for (int c = 0; c < 3; c++) {
            if (c) jb_lit(b, ", ");
            if (rows[r][c].d == 1) jb_int(b, (int64_t)rows[r][c].n);
            else {
                double v;
                py_float_of_fraction(rows[r][c], &v);
                jb_float(b, v);
            }
        }
        jb_lit(b, "]");
    }
    jb_lit(b, "]");
}

void sel_put_basis_text(jbuf *b, const rat (*rows)[3])
{
    jbuf t = {0};
    jb_reset(&t);
    for (int r = 0; r < 3; r++) {
        if (r) jb_lit(&t, ",");
        jb_lit(&t, "(");
        for (int c = 0; c < 3; c++) {
            if (c) jb_lit(&t, ",");
            put_frac_text(&t, rows[r][c]);
        }
        jb_lit(&t, ")");
    }
    jb_string(b, t.data, t.used);
    jb_free(&t);
}
void sel_put_origin_text(jbuf *b, int parent_sg, int subgroup, const rat (*rows)[3], const int64_t *origin)
{
    int identity = 1;
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) identity &= rat_eq(rows[r][c], rat_int(r == c));
    if (parent_sg == subgroup && ((!origin[0] && !origin[1] && !origin[2]) || identity)) {
        jb_lit(b, "\"(0,0,0)\"");
        return;
    }
    jbuf t = {0};
    jb_reset(&t);
    int64_t den = origin[3] ? origin[3] : 1;
    jb_lit(&t, "(");
    for (int i = 0; i < 3; i++) {
        if (i) jb_lit(&t, ",");
        put_frac_text(&t, rat_make(origin[i], den));
    }
    jb_lit(&t, ")");
    jb_string(b, t.data, t.used);
    jb_free(&t);
}

int64_t sel_subgroup_index(int parent_sg, int subgroup, int64_t size, int *has)
{
    int64_t parent_order = sel_point_group_order(parent_sg), subgroup_order = sel_point_group_order(subgroup);
    int64_t value = size * parent_order;
    *has = subgroup_order > 0 && value % subgroup_order == 0;
    return *has ? value / subgroup_order : 0;
}

void sel_put_strings(jbuf *b, const strings *s)
{
    jb_lit(b, "[");
    for (int i = 0; i < s->count; i++) {
        if (i) jb_lit(b, ", ");
        jb_cstring(b, s->items[i]);
    }
    jb_lit(b, "]");
}

void sel_put_joined(jbuf *b, const strings *s)
{
    jbuf t = {0};
    jb_reset(&t);
    for (int i = 0; i < s->count; i++) {
        if (i) jb_lit(&t, ",");
        jb_lit(&t, s->items[i]);
    }
    jb_string(b, t.data, t.used);
    jb_free(&t);
}

/* The subgroup, basis, basis_text, basis_denominator, source values, det and
 * origin members shared by every OPD row kind. */
void sel_put_subgroup(jbuf *b, int parent_sg, int subgroup, const int64_t *basis, const int64_t *origin, const sel_setting *displayed,
                      const rat (*rows)[3], int64_t det, const jbuf *subgroup_json, const jbuf *before_det)
{
    jb_lit(b, ", \"subgroup\": ");
    if (subgroup_json) {
        jb_raw(b, subgroup_json->data, subgroup_json->used);
    } else {
        jb_lit(b, "{\"number\": ");
        jb_int(b, subgroup);
        jb_lit(b, ", \"symbol\": ");
        char *symbol = sel_default_symbol(subgroup);
        jb_cstring(b, symbol ? symbol : "");
        free(symbol);
        jb_lit(b, "}");
    }
    jb_lit(b, ", \"basis\": ");
    sel_put_basis(b, rows);
    jb_lit(b, ", \"basis_text\": ");
    sel_put_basis_text(b, rows);
    jb_lit(b, ", \"basis_denominator\": ");
    jb_int(b, displayed->denominator);
    jb_lit(b, ", \"source_basis_values\": [");
    for (int i = 0; i < 9; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, basis[i]);
    }
    jb_lit(b, "], \"source_origin_values\": [");
    for (int i = 0; i < 4; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, origin[i]);
    }
    jb_lit(b, "]");
    if (before_det) jb_raw(b, before_det->data, before_det->used);
    jb_lit(b, ", \"det\": ");
    jb_int(b, det);
    jb_lit(b, ", \"origin\": ");
    sel_put_origin_text(b, parent_sg, subgroup, rows, displayed->origin);
}

/* subgroup_change_setting_cinter and present_opd_basis_rows of one row */
int sel_display(int parent_sg, int subgroup, const int64_t *basis, const int64_t *origin, int has_setting, int64_t setting_id, int has_cell,
                const double *cell, int parametric, int coupled, int magnetic_group, sel_setting *displayed, rat (*rows)[3])
{
    int same = has_setting && parent_sg == subgroup;
    MD_TRY(sel_change_setting(parent_sg, subgroup, basis, origin, has_setting, setting_id, same, setting_id, displayed));
    return sel_present_basis(subgroup, displayed, has_cell, cell, parametric, coupled, magnetic_group, rows);
}
