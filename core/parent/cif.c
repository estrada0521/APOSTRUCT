/* parent.read_cif_summary: a CIF 1.1 reader for one data block (pairs, loops,
 * quoted values, text fields, comments) and the summary the parent stage reads
 * from it. save_ frames, global_ and stop_ are rejected. */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "json.h"
#include "pyconv.h"

typedef struct {
    const char *s;
    size_t n;
} span;

typedef struct {
    int loop;
    int ntags;
    span *tags;
    size_t nvalues;
    span *values;
} item;

typedef struct {
    span name;
    int nitems;
    item *items;
} block;

typedef struct {
    const char *p, *end, *start;
    int line;
    char error[256];
} lexer;

enum { TOK_END, TOK_DATA, TOK_LOOP, TOK_TAG, TOK_VALUE, TOK_BAD };

static int is_space(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f'; }

static int keyword(span t, const char *word)
{
    size_t n = strlen(word);
    if (t.n < n) return 0;
    for (size_t i = 0; i < n; i++)
        if (tolower((unsigned char)t.s[i]) != word[i]) return 0;
    return 1;
}

static int fail(lexer *lx, const char *what)
{
    snprintf(lx->error, sizeof lx->error, "CIF syntax error on line %d: %s", lx->line, what);
    return TOK_BAD;
}

static int next_token(lexer *lx, span *t)
{
    for (;;) {
        while (lx->p < lx->end && is_space(*lx->p)) {
            if (*lx->p == '\n') lx->line++;
            lx->p++;
        }
        if (lx->p < lx->end && *lx->p == '#') {
            while (lx->p < lx->end && *lx->p != '\n') lx->p++;
            continue;
        }
        break;
    }
    if (lx->p >= lx->end) return TOK_END;
    const char *a = lx->p;
    int line_start = a == lx->start || a[-1] == '\n' || a[-1] == '\r';
    if (*a == ';' && line_start) {
        const char *c = a + 1;
        for (;;) {
            while (c < lx->end && *c != '\n') c++;
            if (c >= lx->end) return fail(lx, "unterminated text field");
            lx->line++;
            c++;
            if (c < lx->end && *c == ';') break;
        }
        t->s = a + 1;
        t->n = (size_t)(c - 1 - (a + 1));
        lx->p = c + 1;
        return TOK_VALUE;
    }
    if (*a == '\'' || *a == '"') {
        const char *c = a + 1;
        while (c < lx->end && !(*c == *a && (c + 1 == lx->end || is_space(c[1])))) {
            if (*c == '\n' || *c == '\r') return fail(lx, "unterminated quoted value");
            c++;
        }
        if (c >= lx->end) return fail(lx, "unterminated quoted value");
        t->s = a;
        t->n = (size_t)(c + 1 - a);
        lx->p = c + 1;
        return TOK_VALUE;
    }
    const char *c = a;
    while (c < lx->end && !is_space(*c)) c++;
    t->s = a;
    t->n = (size_t)(c - a);
    lx->p = c;
    if (*a == '_') return TOK_TAG;
    if (keyword(*t, "data_")) {
        if (t->n == 5) return fail(lx, "data_ without a block name");
        return TOK_DATA;
    }
    if (keyword(*t, "loop_") && t->n == 5) return TOK_LOOP;
    if (keyword(*t, "save_")) return fail(lx, "save frames are not supported");
    if (keyword(*t, "global_") && t->n == 7) return fail(lx, "global_ is not supported");
    if (keyword(*t, "stop_") && t->n == 5) return fail(lx, "stop_ is not supported");
    return TOK_VALUE;
}

static int span_ieq(span a, span b)
{
    if (a.n != b.n) return 0;
    for (size_t i = 0; i < a.n; i++)
        if (tolower((unsigned char)a.s[i]) != tolower((unsigned char)b.s[i])) return 0;
    return 1;
}

static int span_ieq_c(span a, const char *b) { return span_ieq(a, (span){b, strlen(b)}); }

static void block_free(block *b)
{
    for (int i = 0; i < b->nitems; i++) {
        free(b->items[i].tags);
        free(b->items[i].values);
    }
    free(b->items);
    memset(b, 0, sizeof *b);
}

static int push(void **array, size_t size, size_t *count, size_t *capacity)
{
    if (*count < *capacity) return 0;
    size_t grown = *capacity ? *capacity * 2 : 8;
    void *p = realloc(*array, grown * size);
    if (!p) return 1;
    *array = p;
    *capacity = grown;
    return 0;
}

static int has_tag(const block *b, span tag)
{
    for (int i = 0; i < b->nitems; i++)
        for (int k = 0; k < b->items[i].ntags; k++)
            if (span_ieq(b->items[i].tags[k], tag)) return 1;
    return 0;
}

/* Parse the text into its only data block. 0 on success, else lx->error. */
static int parse(lexer *lx, block *out)
{
    memset(out, 0, sizeof *out);
    span t;
    int tok = next_token(lx, &t), blocks = 0;
    size_t nitems = 0, cap = 0;
    while (tok != TOK_END) {
        if (tok == TOK_BAD) return 1;
        if (tok == TOK_DATA) {
            if (++blocks > 1) {
                snprintf(lx->error, sizeof lx->error, "CIF must contain exactly one data block");
                return 1;
            }
            out->name = (span){t.s + 5, t.n - 5};
            tok = next_token(lx, &t);
            continue;
        }
        if (!blocks) {
            fail(lx, "content before the first data block");
            return 1;
        }
        if (tok == TOK_VALUE) {
            fail(lx, "value without a tag");
            return 1;
        }
        if (push((void **)&out->items, sizeof(item), &nitems, &cap)) goto memory;
        item *it = out->items + nitems;
        memset(it, 0, sizeof *it);
        out->nitems = (int)++nitems;
        if (tok == TOK_TAG) {
            if (has_tag(out, t)) {
                fail(lx, "duplicate tag");
                return 1;
            }
            it->tags = malloc(sizeof(span));
            it->values = malloc(sizeof(span));
            if (!it->tags || !it->values) goto memory;
            it->tags[0] = t;
            it->ntags = 1;
            tok = next_token(lx, &t);
            if (tok == TOK_BAD) return 1;
            if (tok != TOK_VALUE) {
                fail(lx, "tag without a value");
                return 1;
            }
            it->values[0] = t;
            it->nvalues = 1;
            tok = next_token(lx, &t);
            continue;
        }
        it->loop = 1;
        size_t ntags = 0, tcap = 0, nvalues = 0, vcap = 0;
        tok = next_token(lx, &t);
        while (tok == TOK_TAG) {
            if (has_tag(out, t)) {
                fail(lx, "duplicate tag");
                return 1;
            }
            if (push((void **)&it->tags, sizeof(span), &ntags, &tcap)) goto memory;
            it->tags[ntags++] = t;
            it->ntags = (int)ntags;
            tok = next_token(lx, &t);
        }
        if (!ntags) {
            fail(lx, "loop_ without tags");
            return 1;
        }
        while (tok == TOK_VALUE) {
            if (push((void **)&it->values, sizeof(span), &nvalues, &vcap)) goto memory;
            it->values[nvalues++] = t;
            it->nvalues = nvalues;
            tok = next_token(lx, &t);
        }
        if (nvalues % ntags) {
            fail(lx, "loop values do not fill whole rows");
            return 1;
        }
    }
    if (!blocks) {
        snprintf(lx->error, sizeof lx->error, "CIF contains no data block");
        return 1;
    }
    return 0;
memory:
    snprintf(lx->error, sizeof lx->error, "out of memory");
    return 1;
}

/* Python str.isspace over ASCII */
static int py_space(unsigned char c) { return c == ' ' || (c >= '\t' && c <= '\r') || (c >= 0x1c && c <= 0x1f); }

static span py_strip(span v)
{
    while (v.n && py_space((unsigned char)v.s[0])) {
        v.s++;
        v.n--;
    }
    while (v.n && py_space((unsigned char)v.s[v.n - 1])) v.n--;
    return v;
}

/* parent._strip_cif */
static span strip_cif(span v)
{
    v = py_strip(v);
    if (v.n >= 2 && v.s[0] == v.s[v.n - 1] && (v.s[0] == '\'' || v.s[0] == '"')) {
        v.s++;
        v.n -= 2;
    }
    return v;
}

/* gemmi Block::find_value: a pair, else a column of a one-row loop */
static const span *find_value(const block *b, const char *tag)
{
    for (int i = 0; i < b->nitems; i++)
        if (!b->items[i].loop && span_ieq_c(b->items[i].tags[0], tag)) return b->items[i].values;
    for (int i = 0; i < b->nitems; i++) {
        const item *it = b->items + i;
        if (!it->loop || it->nvalues != (size_t)it->ntags) continue;
        for (int k = 0; k < it->ntags; k++)
            if (span_ieq_c(it->tags[k], tag)) return it->values + k;
    }
    return 0;
}

/* The loop of the first item carrying tag, or 0 when that item is a pair. */
static const item *find_loop(const block *b, const char *tag, int *column)
{
    for (int i = 0; i < b->nitems; i++)
        for (int k = 0; k < b->items[i].ntags; k++)
            if (span_ieq_c(b->items[i].tags[k], tag)) {
                *column = k;
                return b->items[i].loop ? b->items + i : 0;
            }
    return 0;
}

/* parent._cif_value over a tag list */
static int cif_value(const block *b, const char *const *tags, int ntags, span *out)
{
    for (int i = 0; i < ntags; i++) {
        const span *raw = find_value(b, tags[i]);
        if (!raw) continue;
        span v = strip_cif(*raw);
        if (v.n) {
            *out = v;
            return 1;
        }
    }
    return 0;
}

static int contains_ci(span text, const char *word)
{
    size_t n = strlen(word);
    for (size_t i = 0; i + n <= text.n; i++) {
        size_t k = 0;
        while (k < n && tolower((unsigned char)text.s[i + k]) == word[k]) k++;
        if (k == n) return 1;
    }
    return 0;
}

static void put_span(jbuf *b, span s) { jb_string(b, s.s, s.n); }

static void put_optional(jbuf *b, int has, span s)
{
    if (has) put_span(b, s);
    else jb_lit(b, "null");
}

static int column_of(const item *it, const char *name)
{
    for (int k = 0; k < it->ntags; k++)
        if (span_ieq_c(it->tags[k], name)) return k;
    return -1;
}

static int summary(jbuf *b, const block *blk, int has_path, span path, char *error, size_t error_size)
{
    static const char *number_tags[] = {"_symmetry_Int_Tables_number", "_space_group_IT_number", "_space_group.it_number",
                                        "_space_group_IT_number"};
    span number_text;
    if (!cif_value(blk, number_tags, 4, &number_text)) {
        snprintf(error, error_size, "CIF does not contain an explicit IT space-group number");
        return 1;
    }
    size_t digit = 0;
    while (digit < number_text.n && !isdigit((unsigned char)number_text.s[digit])) digit++;
    if (digit == number_text.n) {
        jbuf r = {0};
        jb_reset(&r);
        py_repr_text(&r, number_text.s, (int64_t)number_text.n);
        snprintf(error, error_size, "cannot parse IT space-group number from %.*s", (int)r.used, r.data ? r.data : "");
        jb_free(&r);
        return 1;
    }
    size_t end = digit;
    while (end < number_text.n && isdigit((unsigned char)number_text.s[end])) end++;
    while (digit + 1 < end && number_text.s[digit] == '0') digit++;
    span digits = {number_text.s + digit, end - digit};

    static const char *symbol_tags[] = {"_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt", "_space_group.name_H-M_alt"};
    span symbol = {"", 0};
    int has_symbol = cif_value(blk, symbol_tags, 3, &symbol);

    static const char *pref1[] = {"_space_group.reference_setting"}, *pref2[] = {"_space_group_reference_setting"},
                      *pref3[] = {"_symmetry_space_group_setting"};
    span preference = {"", 0};
    int has_preference = cif_value(blk, pref1, 1, &preference) || cif_value(blk, pref2, 1, &preference) || cif_value(blk, pref3, 1, &preference);
    if (!has_preference) {
        static const char *name_tags[] = {"_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt", "_space_group.name-H-M_alt",
                                          "_space_group.name_H-M_alt"};
        for (int t = 0; t < 4 && !has_preference; t++) {
            span value;
            if (!cif_value(blk, name_tags + t, 1, &value)) continue;
            for (size_t open = 0; open < value.n && !has_preference; open++) {
                if (value.s[open] != '(') continue;
                size_t close = open + 1;
                while (close < value.n && value.s[close] != ')') close++;
                if (close >= value.n) break;
                span inner = {value.s + open + 1, close - open - 1};
                if (contains_ci(inner, "origin") || contains_ci(inner, "choice") || contains_ci(inner, "axis") || contains_ci(inner, "setting")) {
                    preference = py_strip(inner);
                    has_preference = 1;
                }
            }
        }
    }

    static const char *formula_tags[] = {"_chemical_formula_structural", "_chemical_formula_sum", "_chemical_formula_moiety",
                                         "_chemical_formula.structural", "_chemical_formula.sum", "_chemical_formula.moiety"};
    span formula;
    int has_formula = cif_value(blk, formula_tags, 6, &formula);
    if (has_formula && formula.n == 1 && (formula.s[0] == '.' || formula.s[0] == '?')) has_formula = 0;

    jb_lit(b, "{\"path\": ");
    put_optional(b, has_path, path);
    jb_lit(b, ", \"block\": ");
    put_span(b, blk->name);
    jb_lit(b, ", \"formula\": ");
    put_optional(b, has_formula, formula);
    jb_lit(b, ", \"parent\": {\"number\": ");
    jb_raw(b, digits.s, digits.n);
    jb_lit(b, ", \"symbol\": ");
    put_span(b, has_symbol ? symbol : (span){"", 0});
    jb_lit(b, "}, \"explicit_parent_preference\": ");
    put_optional(b, has_preference, preference);
    jb_lit(b, ", \"lattice\": {");
    static const char *keys[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
    static const char *cell[6][2] = {{"_cell_length_a", "_cell.length_a"},     {"_cell_length_b", "_cell.length_b"},
                                     {"_cell_length_c", "_cell.length_c"},     {"_cell_angle_alpha", "_cell.angle_alpha"},
                                     {"_cell_angle_beta", "_cell.angle_beta"}, {"_cell_angle_gamma", "_cell.angle_gamma"}};
    for (int i = 0; i < 6; i++) {
        if (i) jb_lit(b, ", ");
        jb_cstring(b, keys[i]);
        jb_lit(b, ": ");
        span v;
        int has = cif_value(blk, cell[i], 2, &v);
        put_optional(b, has, v);
    }
    jb_lit(b, "}, \"atom_sites\": [");
    int column;
    const item *atoms = find_loop(blk, "_atom_site_label", &column);
    if (atoms) {
        int idx_type = column_of(atoms, "_atom_site_type_symbol"), idx_label = column_of(atoms, "_atom_site_label"),
            idx_mult = column_of(atoms, "_atom_site_symmetry_multiplicity"), idx_x = column_of(atoms, "_atom_site_fract_x"),
            idx_y = column_of(atoms, "_atom_site_fract_y"), idx_z = column_of(atoms, "_atom_site_fract_z"),
            idx_occ = column_of(atoms, "_atom_site_occupancy"), idx_wyckoff = -1;
        static const char *wyckoff_tags[] = {"_atom_site_Wyckoff_label", "_atom_site_wyckoff_label", "_atom_site_Wyckoff_symbol",
                                             "_atom_site_wyckoff_symbol"};
        for (int t = 0; t < 4 && idx_wyckoff < 0; t++) idx_wyckoff = column_of(atoms, wyckoff_tags[t]);
        size_t rows = atoms->nvalues / (size_t)atoms->ntags;
        for (size_t row = 0; row < rows; row++) {
            const span *v = atoms->values + row * (size_t)atoms->ntags;
            if (row) jb_lit(b, ", ");
            jb_lit(b, "{");
            int any = 0;
            const struct {
                const char *key;
                int idx;
            } fields[] = {{"type", idx_type}, {"label", idx_label}, {"multiplicity", idx_mult}, {"wyckoff", idx_wyckoff}};
            for (int f = 0; f < 4; f++) {
                if (fields[f].idx < 0) continue;
                if (any) jb_lit(b, ", ");
                any = 1;
                jb_cstring(b, fields[f].key);
                jb_lit(b, ": ");
                put_span(b, strip_cif(v[fields[f].idx]));
            }
            if (idx_x >= 0 && idx_y >= 0 && idx_z >= 0) {
                if (any) jb_lit(b, ", ");
                any = 1;
                jb_lit(b, "\"fract\": [");
                put_span(b, strip_cif(v[idx_x]));
                jb_lit(b, ", ");
                put_span(b, strip_cif(v[idx_y]));
                jb_lit(b, ", ");
                put_span(b, strip_cif(v[idx_z]));
                jb_lit(b, "]");
            }
            if (idx_occ >= 0) {
                if (any) jb_lit(b, ", ");
                jb_lit(b, "\"occupancy\": ");
                put_span(b, strip_cif(v[idx_occ]));
            }
            jb_lit(b, "}");
        }
    }
    jb_lit(b, "], \"symmetry_operations\": [");
    static const char *op_tags[] = {"_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz", "_space_group_symop.operation_xyz"};
    for (int t = 0; t < 3; t++) {
        const item *ops = find_loop(blk, op_tags[t], &column);
        if (!ops) continue;
        int first = 1;
        size_t rows = ops->nvalues / (size_t)ops->ntags;
        for (size_t row = 0; row < rows; row++) {
            span op = py_strip(strip_cif(ops->values[row * (size_t)ops->ntags + (size_t)column]));
            if (!op.n) continue;
            if (!first) jb_lit(b, ", ");
            first = 0;
            put_span(b, op);
        }
        break;
    }
    jb_lit(b, "]}");
    return 0;
}

static jbuf out;

static void put_error(int kind, const char *message)
{
    jb_reset(&out);
    jb_lit(&out, "{\"error\": {\"kind\": ");
    jb_int(&out, kind);
    jb_lit(&out, ", \"status\": 0, \"message\": ");
    jb_cstring(&out, message);
    jb_lit(&out, "}}");
}

/* read_cif_summary (path, from_file) or read_cif_summary_from_text as JSON.
 * Errors: kind 1 ValueError. */
const char *cif_summary_json(const char *data, int64_t length, int from_file)
{
    jb_reset(&out);
    char *text = 0;
    size_t size = (size_t)length;
    span path = {data, (size_t)length};
    if (from_file) {
        char *name = malloc(size + 1);
        if (!name) return 0;
        memcpy(name, data, size);
        name[size] = 0;
        FILE *f = fopen(name, "rb");
        free(name);
        if (!f) {
            char message[256];
            snprintf(message, sizeof message, "cannot read CIF file: %.*s", (int)(size > 200 ? 200 : size), data);
            put_error(1, message);
            return out.failed ? 0 : out.data;
        }
        size_t cap = 1 << 16, used = 0;
        text = malloc(cap);
        while (text) {
            size_t got = fread(text + used, 1, cap - used, f);
            used += got;
            if (used < cap) break;
            char *grown = realloc(text, cap * 2);
            if (!grown) {
                free(text);
                text = 0;
                break;
            }
            text = grown;
            cap *= 2;
        }
        fclose(f);
        if (!text) return 0;
        data = text;
        size = used;
    }
    lexer lx = {data, data + size, data, 1, {0}};
    block blk;
    char error[512] = {0};
    if (parse(&lx, &blk)) put_error(1, lx.error);
    else if (summary(&out, &blk, from_file, path, error, sizeof error)) put_error(1, error);
    block_free(&blk);
    free(text);
    return out.failed ? 0 : out.data;
}
