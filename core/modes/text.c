/* The mode details text (the Web complete-modes page) of a full modes state. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "md.h"

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

static int64_t code_points(const char *s, size_t n)
{
    int64_t count = 0;
    for (size_t i = 0; i < n; i++) count += ((unsigned char)s[i] & 0xC0) != 0x80;
    return count;
}

/* f"{text:<width}" */
static void put_left(jbuf *b, const char *s, size_t n, int width)
{
    jb_raw(b, s, n);
    for (int64_t pad = width - code_points(s, n); pad > 0; pad--) jb_raw(b, " ", 1);
}

static void put_str_left(jbuf *b, const jvalue *v, int width)
{
    jbuf t = {0};
    jb_reset(&t);
    if (jtruthy(v)) py_str_of(&t, v);
    put_left(b, t.data ? t.data : "", t.used, width);
    jb_free(&t);
}

static void put_fixed(jbuf *b, double v, int width, int digits)
{
    char text[64];
    int n = snprintf(text, sizeof text, "%*.*f", width, digits, v);
    jb_raw(b, text, (size_t)n);
}

/* a CIF value: a number, or text with an optional uncertainty in parentheses (5.380(1)) */
static int cif_value(const jvalue *v, double *out)
{
    return v && v->kind == J_STRING ? float_cif_number(v, out) : py_float_of(v, out);
}

static int float_or_zero(const jvalue *v, double *out)
{
    if (!v) {
        *out = 0.0;
        return 0;
    }
    return cif_value(v, out) ? md_fail(MDE_VALUE, "could not convert string to float") : 0;
}

static int lattice_lines(jbuf *b, const jvalue *lattice, const double *values, int digits)
{
    static const char *const keys[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
    double v[6];
    for (int i = 0; i < 6; i++) {
        if (values) v[i] = values[i];
        else MD_TRY(float_or_zero(member(lattice, keys[i]), v + i));
    }
    for (int line = 0; line < 2; line++) {
        for (int i = 3 * line; i < 3 * line + 3; i++) {
            if (i % 3) jb_lit(b, ", ");
            jb_lit(b, keys[i]);
            jb_lit(b, "=");
            put_fixed(b, v[i], 0, digits);
        }
        jb_lit(b, "\n");
    }
    return 0;
}

/* the supercell lattice: 0 and the six values, or 1 to keep the parent's */
static int supercell(const jvalue *parent, const jvalue *basis, double *out)
{
    if (!jtruthy(basis) || !jis(parent, J_OBJECT)) return 1;
    static const char *const keys[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
    double p[6];
    for (int i = 0; i < 6; i++) {
        const jvalue *v = member(parent, keys[i]);
        if (!v || jis(v, J_NULL) || cif_value(v, p + i)) return 1;
    }
    double a = p[0], b = p[1], c = p[2];
    double cal = geo_cos_deg(p[3]), cbe = geo_cos_deg(p[4]), cga = geo_cos_deg(p[5]), sga = geo_sin_deg(p[5]);
    double cx = c * cbe, cy = c * (cal - cbe * cga) / sga;
    double square = c * c - cx * cx - cy * cy;
    double cz = sqrt(square > 0.0 ? square : 0.0);
    double P[3][3] = {{a, 0.0, 0.0}, {b * cga, b * sga, 0.0}, {cx, cy, cz}};
    double rows[3][3];
    for (int r = 0; r < 3; r++) {
        const jvalue *row = jat(basis, r);
        for (int axis = 0; axis < 3; axis++) {
            double terms[3];
            for (int k = 0; k < 3; k++) {
                double factor;
                if (py_float_of(jat(row, k), &factor)) {
                    md_fail(MDE_VALUE, "could not convert basis to float");
                    return -1;
                }
                terms[k] = factor * P[k][axis];
            }
            rows[r][axis] = py_sum(terms, 3);
        }
    }
    double length[3];
    for (int r = 0; r < 3; r++) {
        double squares[3] = {rows[r][0] * rows[r][0], rows[r][1] * rows[r][1], rows[r][2] * rows[r][2]};
        length[r] = sqrt(py_sum(squares, 3));
    }
    int pairs[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    for (int i = 0; i < 3; i++) out[i] = length[i];
    for (int i = 0; i < 3; i++) {
        int u = pairs[i][0], w = pairs[i][1];
        if (!length[u] || !length[w]) {
            out[3 + i] = 0.0;
            continue;
        }
        double products[3] = {rows[u][0] * rows[w][0], rows[u][1] * rows[w][1], rows[u][2] * rows[w][2]};
        double d = py_sum(products, 3) / (length[u] * length[w]);
        d = d < 1.0 ? d : 1.0;
        d = -1.0 > d ? -1.0 : d;
        out[3 + i] = geo_acos_deg(d);
    }
    return 0;
}

/* ------------------------------------------------------------ Wyckoff formula */

typedef struct {
    const char *at, *end;
    const jvalue *params;
    int failed;
} expr;

static void skip(expr *e)
{
    while (e->at < e->end && (*e->at == ' ' || *e->at == '\t')) e->at++;
}

static rat parse_sum(expr *e);

static int name_char(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_' || (c >= '0' && c <= '9'); }

static rat parse_atom(expr *e)
{
    skip(e);
    if (e->at >= e->end) {
        e->failed = 1;
        return rat_int(0);
    }
    char c = *e->at;
    if (c == '(') {
        e->at++;
        rat v = parse_sum(e);
        skip(e);
        if (e->at >= e->end || *e->at != ')') e->failed = 1;
        else e->at++;
        return v;
    }
    if (c == '-' || c == '+') {
        e->at++;
        rat v = parse_atom(e);
        return c == '-' ? rat_neg(v) : v;
    }
    if ((c >= '0' && c <= '9') || c == '.') {
        const char *start = e->at;
        while (e->at < e->end && ((*e->at >= '0' && *e->at <= '9') || *e->at == '.')) e->at++;
        const char *stop = e->at;
        /* the Source form 1/2x reads (1/2)*x, as _safe_eval_fraction_expr rewrites it */
        if (stop < e->end && *stop == '/') {
            const char *c = stop + 1;
            while (c < e->end && *c >= '0' && *c <= '9') c++;
            if (c > stop + 1 && c < e->end && (*c == 'x' || *c == 'y' || *c == 'z') && !(c + 1 < e->end && name_char(c[1]))) stop = c;
        }
        int coefficient = stop < e->end && (*stop == 'x' || *stop == 'y' || *stop == 'z') && !(stop + 1 < e->end && name_char(stop[1]));
        rat v;
        if (py_fraction_of_text(start, stop - start, &v)) e->failed = 1;
        e->at = stop;
        if (coefficient) {
            rat name = parse_atom(e);
            return rat_mul(v, name);
        }
        return v;
    }
    if (name_char(c)) {
        const char *start = e->at;
        while (e->at < e->end && name_char(*e->at)) e->at++;
        char key[32];
        size_t n = (size_t)(e->at - start);
        if (n >= sizeof key) {
            e->failed = 1;
            return rat_int(0);
        }
        memcpy(key, start, n);
        key[n] = 0;
        const jvalue *value = member(e->params, key);
        if (!value) {
            e->failed = 1;
            return rat_int(0);
        }
        jbuf text = {0};
        jb_reset(&text);
        py_str_of(&text, value);
        rat v;
        if (py_fraction_of_text(text.data ? text.data : "", (int64_t)text.used, &v)) e->failed = 1;
        jb_free(&text);
        return v;
    }
    e->failed = 1;
    return rat_int(0);
}

static rat parse_product(expr *e)
{
    rat v = parse_atom(e);
    for (;;) {
        skip(e);
        if (e->at >= e->end || (*e->at != '*' && *e->at != '/')) return v;
        char op = *e->at++;
        rat w = parse_atom(e);
        if (op == '/') {
            if (rat_zero(w)) {
                e->failed = 1;
                return rat_int(0);
            }
            v = rat_div(v, w);
        } else {
            v = rat_mul(v, w);
        }
    }
}

static rat parse_sum(expr *e)
{
    rat v = parse_product(e);
    for (;;) {
        skip(e);
        if (e->at >= e->end || (*e->at != '+' && *e->at != '-')) return v;
        char op = *e->at++;
        rat w = parse_product(e);
        v = op == '+' ? rat_add(v, w) : rat_sub(v, w);
    }
}

/* Python float % 1.0 */
static double mod1(double x)
{
    double m = fmod(x, 1.0);
    if (m) {
        if (m < 0) m += 1.0;
    } else {
        m = copysign(0.0, 1.0);
    }
    return m;
}

static int parent_xyz(const jvalue *site, double *xyz)
{
    const jvalue *formula = member(site, "wyckoff_formula"), *params = member(site, "wyckoff_params");
    if (jtruthy(formula) && jis(formula, J_STRING) && (!params || jis(params, J_OBJECT))) {
        const char *s = formula->u.s, *end = formula->u.s + formula->count;
        while (s < end && (*s == ' ' || *s == '\t' || *s == '\n')) s++;
        while (end > s && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\n')) end--;
        if (end - s >= 2 && *s == '(' && end[-1] == ')') {
            s++;
            end--;
        }
        const char *parts[4] = {s, 0, 0, 0};
        int nparts = 1;
        for (const char *c = s; c < end; c++)
            if (*c == ',') {
                if (nparts < 4) parts[nparts] = c + 1;
                nparts++;
            }
        if (nparts == 3) {
            parts[3] = end + 1;
            double out[3];
            int ok = 1;
            rat_overflow = 0;
            for (int k = 0; k < 3 && ok; k++) {
                expr e = {parts[k], parts[k + 1] - 1, params, 0};
                rat v = parse_sum(&e);
                skip(&e);
                if (e.failed || e.at != e.end || rat_overflow || py_float_of_fraction(v, out + k)) ok = 0;
            }
            rat_overflow = 0;
            if (ok) {
                for (int k = 0; k < 3; k++) xyz[k] = mod1(out[k]);
                return 0;
            }
        }
    }
    const jvalue *fract = member(site, "fract");
    for (int k = 0; k < 3; k++) {
        double v = 0.0;
        if (jtruthy(fract) && cif_value(jat(fract, k), &v)) return md_fail(MDE_VALUE, "could not convert string to float");
        xyz[k] = mod1(v);
    }
    return 0;
}

static void site_label(jbuf *b, const jvalue *site)
{
    jbuf w = {0}, m = {0};
    jb_reset(&w);
    jb_reset(&m);
    const jvalue *wyckoff = member(site, "wyckoff"), *mult = member(site, "wyckoff_multiplicity");
    if (!jtruthy(mult)) mult = member(site, "multiplicity");
    if (jtruthy(wyckoff)) py_str_of(&w, wyckoff);
    if (jtruthy(mult)) py_str_of(&m, mult);
    jbuf label = {0};
    jb_reset(&label);
    if (w.used && w.data[0] >= '0' && w.data[0] <= '9') {
        jb_raw(&label, w.data, w.used);
    } else {
        jb_raw(&label, m.data ? m.data : "", m.used);
        jb_raw(&label, w.data ? w.data : "", w.used);
    }
    put_left(b, label.data ? label.data : "", label.used, 5);
    jb_free(&w);
    jb_free(&m);
    jb_free(&label);
}

/* ------------------------------------------------------------ subgroup line */

static void put_py_str(jbuf *b, const jvalue *v)
{
    if (v) py_str_of(b, v);
}

static void put_number_text(jbuf *b, const jvalue *v)
{
    if (!v) return;
    if (jis(v, J_STRING)) {
        jb_raw(b, v->u.s, (size_t)v->count);
        return;
    }
    py_str_of(b, v);
}

static int subgroup_line(jbuf *b, const jvalue *row)
{
    if (!jis(row, J_OBJECT)) return 0;
    const jvalue *iso = jis(member(row, "isotropy"), J_OBJECT) ? member(row, "isotropy") : row;
    const jvalue *sg = member(iso, "subgroup");
    if (!jis(sg, J_OBJECT)) return 0;
    const jvalue *number = member(sg, "display_label");
    if (!jtruthy(number)) number = member(sg, "number");
    if (number) put_py_str(b, number);
    jb_lit(b, " ");
    const jvalue *symbol = member(sg, "symbol");
    if (symbol) put_py_str(b, symbol);
    jb_lit(b, ", basis={");
    const jvalue *basis_text = member(iso, "basis_text");
    if (jtruthy(basis_text)) {
        put_py_str(b, basis_text);
    } else {
        const jvalue *basis = member(iso, "basis");
        if (!jis(basis, J_ARRAY)) {
            if (jtruthy(basis)) put_py_str(b, basis);
        } else {
            jbuf rows = {0};
            jb_reset(&rows);
            int nrows = 0;
            for (int64_t r = 0; r < basis->count; r++) {
                const jvalue *row_value = jat(basis, r);
                if (!jis(row_value, J_ARRAY)) continue;
                if (nrows++) jb_lit(&rows, ",");
                jb_lit(&rows, "(");
                for (int64_t k = 0; k < row_value->count; k++) {
                    if (k) jb_lit(&rows, ",");
                    put_py_str(&rows, jat(row_value, k));
                }
                jb_lit(&rows, ")");
            }
            if (nrows) {
                jb_lit(b, "{");
                jb_raw(b, rows.data, rows.used);
                jb_lit(b, "}");
            } else {
                put_py_str(b, basis);
            }
            jb_free(&rows);
        }
    }
    jb_lit(b, "}, origin=");
    const jvalue *origin = member(iso, "origin");
    if (jis(origin, J_STRING)) {
        jb_raw(b, origin->u.s, (size_t)origin->count);
    } else if (jis(origin, J_ARRAY)) {
        int done = 0;
        if (origin->count == 4) {
            int64_t parts[4];
            int ok = 1;
            for (int k = 0; k < 4 && ok; k++) ok = !py_int_of(jat(origin, k), parts + k) && jis(jat(origin, k), J_INT);
            if (ok && parts[3]) {
                jb_lit(b, "(");
                for (int k = 0; k < 3; k++) {
                    char text[64];
                    int64_t q = parts[k] / parts[3], r = parts[k] % parts[3];
                    if (r && ((r < 0) != (parts[3] < 0))) q--;
                    if (k) jb_lit(b, ",");
                    if (parts[k] % parts[3] == 0) snprintf(text, sizeof text, "%lld", (long long)q);
                    else snprintf(text, sizeof text, "%lld/%lld", (long long)parts[k], (long long)parts[3]);
                    jb_lit(b, text);
                }
                jb_lit(b, ")");
                done = 1;
            }
        }
        if (!done) {
            jb_lit(b, "(");
            for (int64_t k = 0; k < origin->count; k++) {
                if (k) jb_lit(b, ",");
                put_py_str(b, jat(origin, k));
            }
            jb_lit(b, ")");
        }
    } else if (jtruthy(origin)) {
        put_py_str(b, origin);
    }
    jb_lit(b, ", s=");
    put_number_text(b, member(iso, "s"));
    jb_lit(b, ", i=");
    put_number_text(b, member(iso, "i"));
    jb_lit(b, ", k-active= ");
    const jvalue *active = member(iso, "k_active");
    if (jtruthy(active)) {
        put_py_str(b, active);
    } else {
        const jvalue *vectors = member(iso, "k_active_vectors");
        for (int64_t k = 0; jis(vectors, J_ARRAY) && k < vectors->count; k++) {
            if (k) jb_lit(b, ",");
            put_py_str(b, jat(vectors, k));
        }
    }
    jb_lit(b, "\n");
    return 1;
}

/* ------------------------------------------------------------ definitions */

typedef struct {
    const char *s;
    size_t n;
} span;

typedef struct {
    span id;
    int site, order;
} atom_ref;

static int ref_cmp(const void *a, const void *b)
{
    const atom_ref *x = a, *y = b;
    size_t n = x->id.n < y->id.n ? x->id.n : y->id.n;
    int c = memcmp(x->id.s, y->id.s, n);
    return c ? c : x->id.n < y->id.n ? -1 : x->id.n > y->id.n;
}

typedef struct {
    atom_ref *refs;
    int nrefs, nsites;
    const jvalue **site_names;
    const jvalue **site_ids;
} layout;

static void layout_free(layout *l)
{
    free(l->refs);
    free(l->site_names);
    free(l->site_ids);
}

static int is_text(const jvalue *v) { return jis(v, J_STRING) && v->count > 0; }

/* The undistorted child sites as atom id -> (site, member order). */
static int layout_of(const jvalue *atoms, layout *l)
{
    memset(l, 0, sizeof *l);
    int64_t n = atoms->count, total = 0;
    l->site_names = calloc((size_t)n + 1, sizeof *l->site_names);
    l->site_ids = calloc((size_t)n + 1, sizeof *l->site_ids);
    if (!l->site_names || !l->site_ids) return md_fail(MDE_MEMORY, "out of memory");
    for (int64_t i = 0; i < n; i++) {
        const jvalue *atom = jat(atoms, i), *child = member(atom, "child_site"), *ids = member(atom, "atom_ids"),
                     *mult = member(atom, "multiplicity"), *label = member(atom, "label");
        int ok = is_text(child) && jis(label, J_STRING) && label->count == child->count && !memcmp(label->u.s, child->u.s, (size_t)child->count) &&
                 jis(ids, J_ARRAY) && ids->count > 0 && jis(mult, J_INT) && mult->u.i == ids->count;
        for (int64_t k = 0; ok && k < ids->count; k++) ok = is_text(jat(ids, k));
        if (!ok) return md_fail(MDE_VALUE, "undistorted child site has no complete atom identity");
        l->site_names[i] = child;
        l->site_ids[i] = ids;
        total += ids->count;
    }
    l->nsites = (int)n;
    l->refs = calloc((size_t)total + 1, sizeof *l->refs);
    if (!l->refs) return md_fail(MDE_MEMORY, "out of memory");
    for (int64_t i = 0; i < n; i++)
        for (int64_t k = 0; k < l->site_ids[i]->count; k++) {
            const jvalue *id = jat(l->site_ids[i], k);
            l->refs[l->nrefs++] = (atom_ref){{id->u.s, (size_t)id->count}, (int)i, (int)k};
        }
    qsort(l->refs, (size_t)l->nrefs, sizeof *l->refs, ref_cmp);
    return 0;
}

static const atom_ref *find_ref(const layout *l, span id)
{
    atom_ref key = {id, 0, 0};
    return bsearch(&key, l->refs, (size_t)l->nrefs, sizeof *l->refs, ref_cmp);
}

typedef struct {
    const double *vectors;
    int64_t nvectors;
} vectors_t;

static int vector_of(const jvalue *vectors, int64_t index, const vectors_t *x, double *out)
{
    if (jis(vectors, J_OBJECT)) {
        const jvalue *pair = jget(vectors, "$f64");
        int64_t at = jat(pair, 0)->u.i + 3 * index;
        for (int k = 0; k < 3; k++) out[k] = at + k < x->nvectors ? x->vectors[at + k] : 0.0;
        return 0;
    }
    const jvalue *row = jat(vectors, index);
    for (int k = 0; k < 3; k++)
        if (py_float_of(jat(row, k), out + k)) out[k] = 0.0;
    return 0;
}

static int64_t vector_count(const jvalue *vectors)
{
    if (jis(vectors, J_OBJECT)) return jat(jget(vectors, "$f64"), 1)->u.i;
    return jis(vectors, J_ARRAY) ? vectors->count : 0;
}

/* One definition's rows, checked against the child-site layout
 * (_mode_rows_grouped_by_child_orbit), printed with its normfactor line. */
static int unique_row_ids(const jvalue *details, const jvalue *definition);

static int definition_block(jbuf *b, const jvalue *details, const jvalue *definition, layout *l, const vectors_t *x, int precision,
                            int vector_precision, int **touched, int **count)
{
    const jvalue *norm = member(definition, "normfactor"), *label = member(definition, "label");
    if (label) put_py_str(b, label);
    jb_lit(b, " normfactor = ");
    if (!norm || jis(norm, J_NULL)) {
        jb_lit(b, "?");
    } else {
        double v;
        MD_TRY(float_or_zero(norm, &v));
        put_fixed(b, v, 0, precision);
    }
    jb_lit(b, "\n");
    const jvalue *atoms = member(details, "atoms"), *blocks = member(details, "atom_blocks"), *vectors = member(definition, "vectors");
    int64_t block_index = 0;
    if (py_int_of(member(definition, "atom_block"), &block_index)) return md_fail(MDE_KEY, "'atom_block'");
    const jvalue *block = jat(blocks, block_index);
    if (!jis(block, J_ARRAY)) return md_fail(MDE_INDEX, "list index out of range");
    int64_t nrows = block->count;
    int64_t nvalues = vector_count(vectors);
    if (nvalues != nrows)
        return md_fail(MDE_VALUE, "zip() argument 2 is %s than argument 1", nvalues < nrows ? "shorter" : "longer");
    if (!nrows) return md_fail(MDE_VALUE, "mode definition has no complete atom identity");
    MD_TRY(unique_row_ids(details, definition));
    if (!l->refs) MD_TRY(layout_of(member(details, "undistorted_atoms"), l));
    if (!*touched) {
        *touched = calloc((size_t)l->nsites + 1, sizeof(int));
        *count = calloc((size_t)l->nsites + 1, sizeof(int));
        if (!*touched || !*count) return md_fail(MDE_MEMORY, "out of memory");
    }
    int *seen = *touched, *members = *count;
    memset(seen, 0, sizeof(int) * (size_t)l->nsites);
    memset(members, 0, sizeof(int) * (size_t)l->nsites);
    const jvalue **row_atoms = malloc(sizeof *row_atoms * (size_t)nrows);
    if (!row_atoms) return md_fail(MDE_MEMORY, "out of memory");
    int status = 0;
    for (int64_t r = 0; r < nrows && !status; r++) {
        int64_t index;
        if (py_int_of(jat(block, r), &index) || !jat(atoms, index)) {
            status = md_fail(MDE_INDEX, "list index out of range");
            break;
        }
        row_atoms[r] = jat(atoms, index);
        const jvalue *id = member(row_atoms[r], "atom_id");
        if (!is_text(id)) {
            status = md_fail(MDE_VALUE, "mode definition has no complete atom identity");
            break;
        }
        const atom_ref *ref = find_ref(l, (span){id->u.s, (size_t)id->count});
        if (ref) {
            seen[ref->site] = 1;
            members[ref->site]++;
        }
    }
    if (!status) {
        int64_t expected = 0;
        for (int s = 0; s < l->nsites && !status; s++) {
            if (!seen[s]) continue;
            if (members[s] != l->site_ids[s]->count) status = md_fail(MDE_VALUE, "mode definition splits an undistorted child site");
            expected += l->site_ids[s]->count;
        }
        if (!status && expected != nrows)
            status = md_fail(MDE_VALUE, "mode definition does not cover the undistorted child atoms: %lld != %lld", (long long)nrows,
                             (long long)expected);
    }
    int site = -1, member_order = 0;
    for (int64_t r = 0; r < nrows && !status; r++) {
        if (site < 0 || member_order == l->site_ids[site]->count) {
            do site++;
            while (site < l->nsites && !seen[site]);
            member_order = 0;
        }
        const jvalue *atom = row_atoms[r], *id = member(atom, "atom_id"), *child = member(atom, "child_site"),
                     *want = jat(l->site_ids[site], member_order), *name = l->site_names[site];
        if (!jis(child, J_STRING) || child->count != name->count || memcmp(child->u.s, name->u.s, (size_t)name->count) ||
            id->count != want->count || memcmp(id->u.s, want->u.s, (size_t)want->count)) {
            status = md_fail(MDE_VALUE, "mode definition child-site identity/order changed");
            break;
        }
        const jvalue *label_value = member(atom, "atom");
        int has_label = jtruthy(label_value);
        int label_ok = member_order == 0 ? has_label && jis(label_value, J_STRING) && label_value->count == name->count &&
                                               !memcmp(label_value->u.s, name->u.s, (size_t)name->count)
                                         : !has_label;
        if (!label_ok) {
            status = md_fail(MDE_VALUE, "mode definition child-site label run is malformed");
            break;
        }
        member_order++;
    }
    for (int64_t r = 0; r < nrows && !status; r++) {
        const jvalue *atom = row_atoms[r], *xyz = member(atom, "xyz");
        put_str_left(b, member(atom, "atom"), 6);
        double p[3], v[3];
        for (int k = 0; k < 3; k++)
            if (!jtruthy(xyz) || py_float_of(jat(xyz, k), p + k)) p[k] = 0.0;
        vector_of(vectors, r, x, v);
        for (int k = 0; k < 3; k++) {
            jb_lit(b, " ");
            put_fixed(b, p[k], precision + 3, precision);
        }
        for (int k = 0; k < 3; k++) {
            jb_lit(b, " ");
            put_fixed(b, v[k], vector_precision + 3, vector_precision);
        }
        jb_lit(b, "\n");
    }
    free(row_atoms);
    return status;
}

static int unique_row_ids(const jvalue *details, const jvalue *definition)
{
    const jvalue *atoms = member(details, "atoms"), *block = jat(member(details, "atom_blocks"), member(definition, "atom_block") ? member(definition, "atom_block")->u.i : 0);
    int64_t n = jis(block, J_ARRAY) ? block->count : 0;
    span *ids = malloc(sizeof *ids * (size_t)(n + 1));
    if (!ids) return md_fail(MDE_MEMORY, "out of memory");
    for (int64_t r = 0; r < n; r++) {
        const jvalue *id = member(jat(atoms, jat(block, r) ? jat(block, r)->u.i : -1), "atom_id");
        if (!is_text(id)) {
            free(ids);
            return md_fail(MDE_VALUE, "mode definition has no complete atom identity");
        }
        ids[r] = (span){id->u.s, (size_t)id->count};
    }
    atom_ref *refs = malloc(sizeof *refs * (size_t)(n + 1));
    if (!refs) {
        free(ids);
        return md_fail(MDE_MEMORY, "out of memory");
    }
    for (int64_t r = 0; r < n; r++) refs[r] = (atom_ref){ids[r], 0, 0};
    qsort(refs, (size_t)n, sizeof *refs, ref_cmp);
    int status = 0;
    for (int64_t r = 1; r < n && !status; r++)
        if (!ref_cmp(refs + r, refs + r - 1)) status = md_fail(MDE_VALUE, "mode definition has no complete atom identity");
    free(ids);
    free(refs);
    return status;
}

static int definitions(jbuf *b, const jvalue *details, const jvalue *list, const char *title, const char *header, layout *l,
                       const vectors_t *x, int precision, int vector_precision, int **touched, int **count)
{
    if (title) {
        jb_lit(b, title);
        jb_lit(b, "\n\n");
        jb_lit(b, header);
        jb_lit(b, "\n");
    }
    for (int64_t i = 0; jis(list, J_ARRAY) && i < list->count; i++) {
        MD_TRY(definition_block(b, details, jat(list, i), l, x, precision, vector_precision, touched, count));
        jb_lit(b, "\n");
    }
    return 0;
}

/* render_mode_detail_text(state).rstrip() + "\n" */
int md_render_text(const jvalue *state, const double *vectors, int64_t nvectors, jbuf *b)
{
    const int precision = 5, vector_precision = 4;
    const jvalue *input = member(state, "input"), *selected = member(state, "selected");
    const jvalue *parent = member(state, "space_group");
    if (!jtruthy(parent)) parent = member(input, "parent");
    const jvalue *details = member(selected, "mode_details"), *status = member(details, "status");
    if (status && !jis(status, J_NULL) && !jstr_is(status, "ok")) {
        const jvalue *reason = member(details, "reason");
        if (jtruthy(reason)) {
            jbuf t = {0};
            jb_reset(&t);
            py_str_of(&t, reason);
            int code = md_fail(MDE_VALUE, "%.*s", (int)t.used, t.data ? t.data : "");
            jb_free(&t);
            return code;
        }
        jbuf t = {0};
        jb_reset(&t);
        py_repr_of(&t, status);
        int code = md_fail(MDE_VALUE, "mode details returned status %.*s", (int)t.used, t.data ? t.data : "");
        jb_free(&t);
        return code;
    }
    const jvalue *lattice = member(input, "lattice"), *orderparam = member(selected, "orderparam");
    const jvalue *iso = jis(member(orderparam, "isotropy"), J_OBJECT) ? member(orderparam, "isotropy") : orderparam;
    double cell[6];
    int super = supercell(lattice, member(iso, "basis"), cell);
    if (super < 0) return md_err.kind;
    const jvalue *atoms = member(details, "undistorted_atoms");
    int all_objects = jis(atoms, J_ARRAY);
    for (int64_t i = 0; all_objects && i < atoms->count; i++) all_objects = jis(jat(atoms, i), J_OBJECT);
    if (!all_objects) return md_fail(MDE_VALUE, "mode details have no canonical undistorted child sites");
    jb_lit(b, "Parent structure\n");
    jbuf head = {0};
    jb_reset(&head);
    const jvalue *number = member(parent, "number"), *symbol = member(parent, "symbol");
    if (number) py_str_of(&head, number);
    jb_lit(&head, " ");
    if (symbol) py_str_of(&head, symbol);
    while (head.used && head.data[head.used - 1] == ' ') head.used--;
    jb_raw(b, head.data ? head.data : "", head.used);
    jb_free(&head);
    jb_lit(b, "\n");
    MD_TRY(lattice_lines(b, lattice, 0, precision));
    jb_lit(b, "\natom site    x         y         z\n");
    const jvalue *sites = member(input, "atom_sites");
    for (int64_t i = 0; jis(sites, J_ARRAY) && i < sites->count; i++) {
        const jvalue *site = jat(sites, i);
        if (!jis(site, J_OBJECT)) continue;
        double xyz[3];
        MD_TRY(parent_xyz(site, xyz));
        put_str_left(b, member(site, "label"), 6);
        jb_lit(b, " ");
        site_label(b, site);
        for (int k = 0; k < 3; k++) {
            jb_lit(b, " ");
            put_fixed(b, xyz[k], precision + 4, precision);
        }
        jb_lit(b, "\n");
    }
    jb_lit(b, "\nSubgroup details\n");
    if (!subgroup_line(b, orderparam)) jb_lit(b, "(not available yet)\n");
    jb_lit(b, "\nUndistorted superstructure\n");
    MD_TRY(lattice_lines(b, lattice, super ? 0 : cell, precision));
    jb_lit(b, "atom site    x         y         z\n");
    for (int64_t i = 0; i < atoms->count; i++) {
        const jvalue *atom = jat(atoms, i), *xyz = member(atom, "xyz");
        put_str_left(b, member(atom, "child_site"), 6);
        jb_lit(b, " ");
        put_str_left(b, member(atom, "site"), 5);
        for (int k = 0; k < 3; k++) {
            double v = 0.0;
            if (jtruthy(xyz)) MD_TRY(float_or_zero(jat(xyz, k), &v));
            jb_lit(b, " ");
            put_fixed(b, v, precision + 4, precision);
        }
        jb_lit(b, "\n");
    }
    jb_lit(b, "\n");
    layout l = {0};
    int *touched = 0, *count = 0;
    vectors_t x = {vectors, nvectors};
    int code = definitions(b, details, member(details, "displacive_definitions"), "Displacive mode definitions",
                           "atom    x        y        z        dx      dy      dz", &l, &x, precision, vector_precision, &touched, &count);
    const jvalue *magnetic = member(details, "magnetic_definitions");
    if (!code && jtruthy(magnetic))
        code = definitions(b, details, magnetic, "Magnetic mode definitions", "atom    x        y        z        dmx     dmy     dmz", &l, &x,
                           precision, vector_precision, &touched, &count);
    free(touched);
    free(count);
    layout_free(&l);
    if (code) return code;
    const jvalue *strain = member(details, "strain_definitions");
    if (jtruthy(strain)) {
        jb_lit(b, "Parent-cell strain mode definitions\n\n    e1      e2      e3      e4      e5      e6\n");
        for (int64_t i = 0; jis(strain, J_ARRAY) && i < strain->count; i++) {
            const jvalue *mode = jat(strain, i), *norm = member(mode, "normfactor"), *components = member(mode, "components");
            const jvalue *label = member(mode, "label");
            if (label) put_py_str(b, label);
            jb_lit(b, " normfactor = ");
            if (!norm || jis(norm, J_NULL)) {
                jb_lit(b, "?");
            } else {
                double v;
                MD_TRY(float_or_zero(norm, &v));
                put_fixed(b, v, 0, precision);
            }
            jb_lit(b, "\n");
            int64_t ncomponents = jtruthy(components) && jis(components, J_ARRAY) ? (components->count < 6 ? components->count : 6) : 6;
            for (int64_t k = 0; k < ncomponents; k++) {
                double v = 0.0;
                if (jtruthy(components)) MD_TRY(float_or_zero(jat(components, k), &v));
                jb_lit(b, " ");
                put_fixed(b, v, vector_precision + 3, vector_precision);
            }
            if (!ncomponents) jb_lit(b, " ");
            jb_lit(b, "\n");
        }
    }
    while (b->used && (b->data[b->used - 1] == '\n' || b->data[b->used - 1] == ' ' || b->data[b->used - 1] == '\t' ||
                       b->data[b->used - 1] == '\r'))
        b->used--;
    if (b->data) b->data[b->used] = 0;
    jb_lit(b, "\n");
    return b->failed ? md_fail(MDE_MEMORY, "out of memory") : 0;
}
