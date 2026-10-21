/* The label texts of the modes: the evaluated selected k vector, the OPD direction texts, each spec's k parameters and displayed k vector,
 * and whether a spec emits every Source OPD row. */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "pyconv.h"
#include "irreps.h"
#include "tables.h"

int case_display_kvector(int sg, int kslot, int param_count, const int64_t *params, int64_t *out);
int k_coordinate_map_words(int sg, int kslot, int64_t *out);
const char *k_numeric_star_json(int parent_sg, int gid, int has_visible, const int64_t *visible_words, int has_source_visible,
                                const int64_t *source_visible_words);
int sp_matrix_rank(const double *m, int rows, int cols, double tol, int *rank);

enum {
    MT_OK = 0,
    MT_TABLES = 1,
    MT_RANGE = 2,         /* ArithmeticError("k-vector value out of range") */
    MT_MAP_RANGE = 3,     /* ArithmeticError of an affine K map */
    MT_INDEX = 4,
    MT_KEY = 5,           /* spec has no direction_matrix */
    MT_REQUEST = 6,       /* the request k parameters do not parse */
    MT_MEMORY = 7,
    MT_RANK = 8,
};

/* ------------------------------------------------------------ text buffer */

static struct {
    char *data;
    size_t used, cap;
    int failed;
} B;

static void put_n(const char *text, size_t length)
{
    if (B.used + length + 1 > B.cap) {
        size_t cap = B.cap ? B.cap : 256;
        while (B.used + length + 1 > cap) cap *= 2;
        char *data = realloc(B.data, cap);
        if (!data) {
            B.failed = 1;
            return;
        }
        B.data = data;
        B.cap = cap;
    }
    memcpy(B.data + B.used, text, length);
    B.used += length;
    B.data[B.used] = 0;
}

static void put(const char *text) { put_n(text, strlen(text)); }

static void reset_text(void)
{
    B.used = 0;
    B.failed = 0;
    put_n("", 0);
}

static void put_i128(i128 value)
{
    char digits[48];
    int n = 0;
    u128 magnitude = value < 0 ? (u128)0 - (u128)value : (u128)value;
    do {
        digits[n++] = (char)('0' + (int)(magnitude % 10));
        magnitude /= 10;
    } while (magnitude);
    if (value < 0) put("-");
    while (n) put_n(&digits[--n], 1);
}

/* str(Fraction) */
static void put_frac(rat v)
{
    put_i128(v.n);
    if (v.d != 1) {
        put("/");
        put_i128(v.d);
    }
}

static int is_space(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v'; }
static int is_word(char c) { return isalnum((unsigned char)c) || c == '_'; }

static void trim(const char **start, const char **end)
{
    while (*start < *end && is_space(**start)) (*start)++;
    while (*end > *start && is_space((*end)[-1])) (*end)--;
}

static void trim_parens(const char **start, const char **end)
{
    while (*start < *end && (**start == '(' || **start == ')')) (*start)++;
    while (*end > *start && ((*end)[-1] == '(' || (*end)[-1] == ')')) (*end)--;
}

static int parse_fraction(const char *text, int64_t length, rat *out) { return py_fraction_of_text(text, length, out); }

/* ------------------------------------------------- restricted expressions */

typedef struct {
    const char *name;
    int64_t length;
    int valid, overflow;
    rat value;
} param;

typedef struct {
    const char *at, *end;
    const param *params;
    int nparams;
    int failed, overflow;
} parser;

static void skip(parser *p)
{
    while (p->at < p->end && (*p->at == ' ' || *p->at == '\t')) p->at++;
}

static rat expression(parser *p);

static int is_decimal(char c) { return c >= '0' && c <= '9'; }
static int is_octal(char c) { return c >= '0' && c <= '7'; }
static int is_binary(char c) { return c == '0' || c == '1'; }
static int is_hex(char c) { return isxdigit((unsigned char)c); }

static rat fail(parser *p)
{
    p->failed = 1;
    return rat_int(0);
}

static rat atom(parser *p)
{
    skip(p);
    if (p->at >= p->end) return fail(p);
    char c = *p->at;
    if (c == '(') {
        p->at++;
        rat value = expression(p);
        skip(p);
        if (p->at >= p->end || *p->at != ')') return fail(p);
        p->at++;
        return value;
    }
    if (isdigit((unsigned char)c) || c == '.') {
        /* Python's numeric literal grammar (underscores between digits). */
        char clean[512];
        size_t n = 0;
        const char *e = p->end;
#define DIGITS(test)                                                                                     \
    do {                                                                                                 \
        int seen = 0;                                                                                    \
        while (p->at < e) {                                                                              \
            if (test(*p->at)) {                                                                          \
                if (n + 1 >= sizeof clean) return fail(p);                                               \
                clean[n++] = *p->at++;                                                                   \
                seen = 1;                                                                                \
            } else if (*p->at == '_' && seen && p->at + 1 < e && test(p->at[1])) {                       \
                p->at++;                                                                                 \
            } else {                                                                                     \
                break;                                                                                   \
            }                                                                                            \
        }                                                                                                \
    } while (0)
        if (c == '0' && p->at + 1 < e && strchr("xXoObB", p->at[1])) {
            char kind = (char)(p->at[1] | 0x20);
            int base = kind == 'x' ? 16 : kind == 'o' ? 8 : 2;
            p->at += 2;
            if (p->at < e && *p->at == '_') p->at++;
            if (base == 16) DIGITS(is_hex);
            else if (base == 8) DIGITS(is_octal);
            else DIGITS(is_binary);
            if (!n || (p->at < e && is_word(*p->at))) return fail(p);
            i128 value = 0;
            rat_overflow = 0;
            for (size_t i = 0; i < n; i++) {
                int d = isdigit((unsigned char)clean[i]) ? clean[i] - '0' : (clean[i] | 0x20) - 'a' + 10;
                value = rat_add_i(rat_mul_i(value, base), d);
            }
            if (rat_overflow) {
                p->overflow = 1;
                return rat_int(1);
            }
            return rat_int(value);
        }
        int is_float = 0;
        DIGITS(is_decimal);
        size_t integer_digits = n;
        if (p->at < e && *p->at == '.') {
            is_float = 1;
            clean[n++] = *p->at++;
            DIGITS(is_decimal);
            if (n == 1) return fail(p);
        }
        if (p->at < e && (*p->at == 'e' || *p->at == 'E') && n) {
            const char *mark = p->at;
            size_t keep = n;
            clean[n++] = *p->at++;
            if (p->at < e && (*p->at == '+' || *p->at == '-')) clean[n++] = *p->at++;
            size_t before = n;
            DIGITS(is_decimal);
            if (n == before) {
                p->at = mark;
                n = keep;
            } else {
                is_float = 1;
            }
        }
        clean[n] = 0;
        if (p->at < e && is_word(*p->at)) return fail(p);
        rat value;
        if (is_float) {
            double parsed = strtod(clean, 0);
            if (!isfinite(parsed)) return fail(p);
            int status = py_fraction_of_float(parsed, &value);
            if (status == 2) {
                p->overflow = 1;
                return rat_int(1);
            }
            if (status) return fail(p);
            return value;
        }
        if (integer_digits > 1 && clean[0] == '0')
            for (size_t z = 0; z < integer_digits; z++)
                if (clean[z] != '0') return fail(p);
        int status = parse_fraction(clean, (int64_t)n, &value);
        if (status == 2) {
            p->overflow = 1;
            return rat_int(1);
        }
        if (status) return fail(p);
        return value;
#undef DIGITS
    }
    if (isalpha((unsigned char)c) || c == '_') {
        const char *start = p->at;
        while (p->at < p->end && is_word(*p->at)) p->at++;
        int64_t length = p->at - start;
        for (int i = p->nparams - 1; i >= 0; i--)
            if (p->params[i].length == length && !memcmp(p->params[i].name, start, (size_t)length)) {
                if (p->params[i].overflow) {
                    p->overflow = 1;
                    return rat_int(1);
                }
                if (!p->params[i].valid) return fail(p);
                return p->params[i].value;
            }
        return fail(p);
    }
    return fail(p);
}

static rat factor(parser *p)
{
    skip(p);
    if (p->at < p->end && (*p->at == '-' || *p->at == '+')) {
        char sign = *p->at++;
        rat value = factor(p);
        return sign == '-' ? rat_neg(value) : value;
    }
    rat value = atom(p);
    skip(p);
    if (p->at + 1 < p->end && p->at[0] == '*' && p->at[1] == '*') return fail(p);
    return value;
}

static rat term(parser *p)
{
    rat value = factor(p);
    for (;;) {
        skip(p);
        if (p->failed || p->at >= p->end || (*p->at != '*' && *p->at != '/')) return value;
        char op = *p->at++;
        if (op == '/' && p->at < p->end && *p->at == '/') return fail(p);
        rat right = factor(p);
        if (op == '*') {
            value = rat_mul(value, right);
        } else {
            if (rat_zero(right)) return fail(p);
            value = rat_div(value, right);
        }
    }
}

static rat expression(parser *p)
{
    rat value = term(p);
    for (;;) {
        skip(p);
        if (p->failed || p->at >= p->end || (*p->at != '+' && *p->at != '-')) return value;
        char op = *p->at++;
        rat right = term(p);
        value = op == '+' ? rat_add(value, right) : rat_sub(value, right);
    }
}

/* re.sub(r"(?<![A-Za-z0-9_])(\d+(?:/\d+)?)([abg])\b", r"\1*\2") and then
 * evaluate_fraction_expression; 0 with *out, 1 for any failure. */
static int evaluate(const char *text, int64_t length, const param *params, int nparams, rat *out)
{
    char *rewritten = malloc((size_t)(2 * length + 1));
    if (!rewritten) return 1;
    int64_t n = 0;
    for (int64_t i = 0; i < length; i++) {
        rewritten[n++] = text[i];
        if (!isdigit((unsigned char)text[i]) || i + 1 >= length) continue;
        char next = text[i + 1];
        if ((next != 'a' && next != 'b' && next != 'g') || (i + 2 < length && is_word(text[i + 2]))) continue;
        int64_t start = i;
        while (start > 0 && isdigit((unsigned char)text[start - 1])) start--;
        if (start == 0 || !is_word(text[start - 1])) rewritten[n++] = '*';
    }
    rat_overflow = 0;
    parser p = {rewritten, rewritten + n, params, nparams, 0, 0};
    rat value = expression(&p);
    skip(&p);
    int overflow = p.overflow || rat_overflow;
    int failed = p.failed || p.at != p.end;
    free(rewritten);
    if (failed) return 1;
    if (overflow) return 2;
    *out = value;
    return 0;
}

/* _evaluated_kvector_text.  names/values: nparams concatenated texts with
 * their lengths, later names overriding earlier ones (the merged dict). */
static int text_status;

int mode_text_status(void) { return text_status; }

const char *evaluated_kvector_text(const char *text, int64_t length, int nparams, const char *names, const int64_t *name_lengths,
                                   const char *values, const int64_t *value_lengths)
{
    reset_text();
    text_status = MT_MEMORY;
    param *params = calloc((size_t)nparams + 1, sizeof(param));
    if (!params) return 0;
    for (int i = 0, name_at = 0, value_at = 0; i < nparams; i++) {
        params[i].name = names + name_at;
        params[i].length = name_lengths[i];
        int status = parse_fraction(values + value_at, value_lengths[i], &params[i].value);
        params[i].valid = !status;
        params[i].overflow = status == 2;
        name_at += (int)name_lengths[i];
        value_at += (int)value_lengths[i];
    }
    const char *start = text, *end = text + length;
    trim(&start, &end);
    trim_parens(&start, &end);
    const char *commas[2];
    int count = 0;
    for (const char *c = start; c < end; c++)
        if (*c == ',') {
            if (count < 2) commas[count] = c;
            count++;
        }
    if (count != 2) {
        put_n(start, (size_t)(end - start));
    } else {
        const char *bounds[4] = {start, commas[0] + 1, commas[1] + 1, end + 1};
        for (int k = 0; k < 3; k++) {
            const char *a = bounds[k], *b = bounds[k + 1] - 1;
            trim(&a, &b);
            if (k) put(",");
            rat value;
            int status = evaluate(a, b - a, params, nparams, &value);
            if (status == 2) {
                free(params);
                text_status = MT_RANGE;
                return 0;
            }
            if (status)
                put_n(a, (size_t)(b - a));
            else
                put_frac(value);
        }
    }
    free(params);
    return B.failed ? 0 : B.data;
}

/* _opd_direction_text after its little_dim choice. */
const char *opd_direction_text(const char *text, int64_t length, int little_dim)
{
    reset_text();
    if (little_dim <= 0) {
        put_n(text, (size_t)length);
        return B.failed ? 0 : B.data;
    }
    const char *start = text, *end = text + length;
    trim(&start, &end);
    trim_parens(&start, &end);
    int64_t parts = 1;
    for (const char *c = start; c < end; c++) parts += *c == ',' || *c == ';';
    if (parts > little_dim && parts % little_dim) {
        put_n(text, (size_t)length);
        return B.failed ? 0 : B.data;
    }
    put("(");
    int64_t index = 0;
    const char *a = start;
    for (const char *c = start;; c++) {
        if (c < end && *c != ',' && *c != ';') continue;
        const char *x = a, *y = c;
        trim(&x, &y);
        if (index) put(parts > little_dim && index % little_dim == 0 ? ";" : ",");
        put_n(x, (size_t)(y - x));
        index++;
        if (c >= end) break;
        a = c + 1;
    }
    put(")");
    return B.failed ? 0 : B.data;
}

/* _direction_matrix_text: row r has row_lengths[r] coefficients from values
 * (rows concatenated); group_width 0 for none. */
const char *direction_matrix_text(int rows, const int64_t *row_lengths, const double *values, double tol, int group_width)
{
    reset_text();
    if (rows <= 0) return B.data;
    int64_t free_count = row_lengths[0];
    int64_t at = 0;
    put("(");
    int grouped = group_width > 0 && rows > group_width && rows % group_width == 0;
    for (int r = 0; r < rows; r++) {
        if (r) put(grouped && r % group_width == 0 ? ";" : ",");
        int terms = 0;
        int64_t width = row_lengths[r] < free_count ? row_lengths[r] : free_count;
        for (int64_t c = 0; c < width; c++) {
            double value = values[at + c];
            if (fabs(value) <= tol) continue;
            int negative = value < 0;
            double magnitude = fabs(value);
            if (negative) put("-");
            else if (terms) put("+");
            if (!(fabs(magnitude - 1.0) <= tol)) {
                int size = snprintf(0, 0, "%.3f", magnitude);
                char *number = malloc((size_t)size + 1);
                if (!number) return 0;
                snprintf(number, (size_t)size + 1, "%.3f", magnitude);
                int n = size;
                while (n && number[n - 1] == '0') n--;
                while (n && number[n - 1] == '.') n--;
                put_n(number, (size_t)n);
                free(number);
            }
            char letter = (char)('a' + c % 26);
            for (int64_t repeat = c / 26; repeat >= 0; repeat--) put_n(&letter, 1);
            terms++;
        }
        if (!terms) put("0");
        at += row_lengths[r];
    }
    put(")");
    return B.failed ? 0 : B.data;
}

/* ------------------------------------------------------------ Source rows */

static struct {
    section_view irr_sg, irr_k, irr_old, irr_type, full_dim, k_dim, k_label, star, star_count, conv2ml;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (geo_load() || section_find("little.little_irr_space_group", &T.irr_sg) || section_find("little.little_irr_k", &T.irr_k) ||
        section_find("little.little_irr_old", &T.irr_old) || section_find("little.little_irr_type", &T.irr_type) ||
        section_find("little.little_irr_full_dim", &T.full_dim) ||
        section_find("little.little_k_dim", &T.k_dim) || section_find("little.little_k_label", &T.k_label) ||
        section_find("little.little_k_star", &T.star) || section_find("little.little_k_star_count", &T.star_count) ||
        section_find("little.little_k_star_conv2ml_pointer", &T.conv2ml))
        return 1;
    T.ready = 1;
    return 0;
}

/* view[index] with Python's negative indexing; *bad on IndexError. */
static int64_t at(const section_view *view, int64_t index, int *bad)
{
    int64_t count = section_count(view);
    if (index < 0) index += count;
    if (index < 0 || index >= count) {
        *bad = 1;
        return 0;
    }
    return section_int(view, index);
}

static int label_is(int lattice, int kslot, const char *label, int64_t length)
{
    const char *text;
    int64_t n;
    if (section_string(&T.k_label, (int64_t)(lattice - 1) * 27 + kslot - 1, &text, &n)) return 0;
    while (n && is_space(*text)) {
        text++;
        n--;
    }
    while (n && is_space(text[n - 1])) n--;
    return n == length && !memcmp(text, label, (size_t)length);
}

typedef struct {
    int has;
    int64_t v[4];
} record4;

typedef struct {
    int count; /* -1 when absent */
    rat v[8];
} params;

/* _source_parameter_record_to_case_k_params */
static int record_params(int64_t gid, const record4 *source, params *out)
{
    out->count = -1;
    if (!source->has || !source->v[3]) return 0;
    int bad = 0;
    int sg = (int)at(&T.irr_sg, gid - 1, &bad), kslot = (int)at(&T.irr_k, gid - 1, &bad);
    if (bad) return MT_INDEX;
    int dim = (int)at(&T.k_dim, (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1, &bad);
    if (bad) return MT_INDEX;
    out->count = dim > 0 ? (dim < 3 ? dim : 3) : 0;
    for (int i = 0; i < out->count; i++) out->v[i] = rat_make(source->v[i], source->v[3]);
    return 0;
}

typedef struct {
    int count;
    rat origin[3], columns[3][3];
} kmap;

static int coordinate_map(int sg, int kslot, kmap *map, int *value_error)
{
    int64_t words[29];
    int status = k_coordinate_map_words(sg, kslot, words);
    if (status == 3) return MT_MAP_RANGE;
    if (status) {
        *value_error = 1;
        return 0;
    }
    map->count = (int)words[1];
    for (int axis = 0; axis < 3; axis++) map->origin[axis] = rat_make(words[5 + 2 * axis], words[6 + 2 * axis]);
    for (int i = 0; i < map->count; i++)
        for (int axis = 0; axis < 3; axis++) map->columns[i][axis] = rat_make(words[11 + 6 * i + 2 * axis], words[12 + 6 * i + 2 * axis]);
    return 0;
}

/* _solve_linear_rational over three equations. */
static int solve(rat m[3][4], int width, rat *out)
{
    int pivots[3], found = 0, pivot_row = 0;
    for (int col = 0; col < width; col++) {
        int pivot = -1;
        for (int row = pivot_row; row < 3 && pivot < 0; row++)
            if (!rat_zero(m[row][col])) pivot = row;
        if (pivot < 0) continue;
        for (int i = 0; i <= width; i++) {
            rat t = m[pivot_row][i];
            m[pivot_row][i] = m[pivot][i];
            m[pivot][i] = t;
        }
        rat f = m[pivot_row][col];
        for (int i = 0; i <= width; i++) m[pivot_row][i] = rat_div(m[pivot_row][i], f);
        for (int row = 0; row < 3; row++) {
            if (row == pivot_row || rat_zero(m[row][col])) continue;
            rat g = m[row][col];
            for (int i = 0; i <= width; i++) m[row][i] = rat_sub(m[row][i], rat_mul(g, m[pivot_row][i]));
        }
        pivots[found++] = col;
        pivot_row++;
    }
    for (int row = 0; row < 3; row++) {
        int zero = 1;
        for (int col = 0; col < width; col++) zero &= rat_zero(m[row][col]);
        if (zero && !rat_zero(m[row][width])) return 0;
    }
    if (found < width) return 0;
    for (int col = 0; col < width; col++) out[col] = rat_int(0);
    for (int i = 0; i < found; i++) out[pivots[i]] = m[i][width];
    return 1;
}

/* _source_record_to_case_k_params */
static int source_record_params(int64_t gid, const record4 *source, params *out)
{
    out->count = -1;
    if (gid <= 0 || !source->has || !source->v[3]) return 0;
    rat target[3];
    for (int axis = 0; axis < 3; axis++) target[axis] = rat_make(source->v[axis], source->v[3]);
    int bad = 0;
    int sg = (int)at(&T.irr_sg, gid - 1, &bad), kslot = (int)at(&T.irr_k, gid - 1, &bad);
    if (bad) return MT_INDEX;
    int dim = (int)at(&T.k_dim, (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1, &bad);
    if (bad) return MT_INDEX;
    if (dim <= 0) {
        out->count = 0;
        return 0;
    }
    int64_t pointer = at(&T.conv2ml, (int64_t)(sg - 1) * 27 + kslot - 1, &bad);
    if (bad) return MT_INDEX;
    if (pointer <= 0 || 16 * pointer > section_count(&T.star)) return 0;
    int64_t raw[16];
    for (int i = 0; i < 16; i++) raw[i] = section_int(&T.star, 16 * (pointer - 1) + i);
    int64_t den = raw[15] ? raw[15] : 1;
    rat rows[3][3];
    for (int c = 0; c < 3; c++)
        for (int v = 0; v < 3; v++) rows[c][v] = rat_make(raw[4 * v + c], den);
    kmap map;
    int value_error = 0;
    int status = coordinate_map(sg, kslot, &map, &value_error);
    if (status) return status;
    if (value_error || map.count != dim) return 0;
    rat m[3][4];
    for (int c = 0; c < 3; c++) {
        rat zero = rat_make(raw[12 + c], den);
        for (int axis = 0; axis < 3; axis++) zero = rat_add(zero, rat_mul(rows[c][axis], map.origin[axis]));
        for (int p = 0; p < dim; p++) {
            rat column = rat_int(0);
            for (int axis = 0; axis < 3; axis++) column = rat_add(column, rat_mul(rows[c][axis], map.columns[p][axis]));
            m[c][p] = column;
        }
        m[c][dim] = rat_sub(target[c], zero);
    }
    if (solve(m, dim, out->v)) out->count = dim;
    return 0;
}

/* KCoordinateMap.evaluate */
static void evaluate_map(const kmap *map, const rat *values, rat *out)
{
    for (int axis = 0; axis < 3; axis++) {
        rat sum = rat_int(0);
        for (int i = 0; i < map->count; i++) sum = rat_add(sum, rat_mul(values[i], map->columns[i][axis]));
        out[axis] = rat_add(map->origin[axis], sum);
    }
}

/* the displayed k vector of a case: 0 with out, -1 when the label or the parameters are unknown */
static int case_vector(int sg, const char *label, int64_t label_length, const params *k, rat *out)
{
    int lattice = geo_sg_lattice(sg), kslot = 0, matches = 0;
    for (int slot = 1; slot <= 27; slot++)
        if (little_k_has_irreps(sg, slot) && label_is(lattice, slot, label, label_length)) {
            kslot = slot;
            matches++;
        }
    if (matches != 1) return -1;
    int64_t words[12], visible[12];
    int count = k->count < 3 ? k->count : 3;
    for (int i = 0; i < count; i++) {
        write_i128(words + 4 * i, k->v[i].n);
        write_i128(words + 4 * i + 2, k->v[i].d);
    }
    int status = case_display_kvector(sg, kslot, count, words, visible);
    if (status == 4) return MT_RANGE;
    if (status) return -1;
    read_rat3(visible, out);
    return 0;
}

static void put_vector(const rat *v)
{
    for (int axis = 0; axis < 3; axis++) {
        if (axis) put(",");
        put_frac(v[axis]);
    }
}

static struct {
    params k;
} R;

enum {
    SL_COUPLED = 1,
    SL_COUPLED_TRUE = 2,
    SL_PRIMARY = 4,
    SL_MAGNETIC_FALSE = 8,
    SL_RECIPROCAL = 16,
    SL_LABEL_MATCHES = 32,
};

static void read_params(int count, const int64_t *words, params *out)
{
    out->count = count < 8 ? count : 8;
    for (int i = 0; i < out->count; i++) out->v[i] = rat_make(read_i128(words + 4 * i), read_i128(words + 4 * i + 2));
}

static record4 record_of(int has, const int64_t *v)
{
    record4 r = {has, {0, 0, 0, 1}};
    if (has) memcpy(r.v, v, sizeof r.v);
    return r;
}

/* The display-vector try block of _spec_display_kvector: 0 with text, 1 to
 * fall through to the case vector, 2 to return the fallback, or MT_. */
static int display_try(int sg, const char *k_label, int64_t k_label_length, int flags, int64_t gid, const record4 *source,
                       const params *spec_case, const params *selected_display, int source_case)
{
    int bad = 0;
    int gsg = (int)at(&T.irr_sg, gid - 1, &bad), kslot = (int)at(&T.irr_k, gid - 1, &bad);
    if (bad) return 1;
    int kdim = (int)at(&T.k_dim, (int64_t)(geo_sg_lattice(gsg) - 1) * 27 + kslot - 1, &bad);
    if (bad) return 1;
    if (kdim <= 0) return (flags & SL_PRIMARY) && !(flags & SL_COUPLED) ? 2 : 1;
    rat visible[3];
    if (source_case) {
        int status = case_vector(sg, k_label, k_label_length, &R.k, visible);
        if (status > 0) return status;
        if (status) return 1;
        put_vector(visible);
        return 0;
    }
    kmap map;
    int value_error = 0;
    int status = coordinate_map(gsg, kslot, &map, &value_error);
    if (status) return status;
    if (value_error) return 1;
    params requested = *selected_display;
    if (requested.count < 0) requested = *spec_case;
    if (requested.count < 0) {
        status = source_record_params(gid, source, &requested);
        if (status == MT_INDEX) return 1;
        if (status) return status;
    }
    if (selected_display->count >= 0 && selected_display->count >= map.count && map.count > 0) {
        evaluate_map(&map, selected_display->v, visible);
        int64_t words[6];
        for (int axis = 0; axis < 3; axis++) {
            if (!rat_fits64(visible[axis])) return MT_RANGE;
            words[2 * axis] = (int64_t)visible[axis].n;
            words[2 * axis + 1] = (int64_t)visible[axis].d;
        }
        const char *json = k_numeric_star_json(gsg, (int)gid, 1, words, 1, words);
        if (!json) return 1;
        const char *open = strchr(json, '"');
        if (open) {
            const char *close = strchr(open + 1, '"');
            const char *a = open + 1, *b = close;
            trim(&a, &b);
            trim_parens(&a, &b);
            put_n(a, (size_t)(b - a));
            return 0;
        }
    }
    if (requested.count >= 0 && requested.count >= map.count) {
        evaluate_map(&map, requested.v, visible);
    } else {
        rat raw[3];
        status = case_vector(sg, k_label, k_label_length, &R.k, raw);
        if (status > 0) return status;
        if (status) return 1;
        rat m[3][4], solved[3];
        for (int axis = 0; axis < 3; axis++) {
            for (int p = 0; p < map.count; p++) m[axis][p] = map.columns[p][axis];
            m[axis][map.count] = rat_sub(raw[axis], map.origin[axis]);
        }
        if (!solve(m, map.count, solved)) return 1;
        for (int p = 0; p < map.count; p++)
            if (!(rat_cmp(solved[p], rat_make(-1, 2)) >= 0 && rat_cmp(solved[p], rat_int(1)) < 0)) solved[p] = rat_mod1(solved[p]);
        evaluate_map(&map, solved, visible);
    }
    put_vector(visible);
    return 0;
}

/* _case_k_params and _spec_display_kvector of one render spec.  Records are
 * four ints (has_* 0 when absent or a non-integral record); parameter lists are
 * int128 rational words with count -1 when absent; request_count -1 when the
 * request k parameters failed to parse.  The text is spec_label_text(). */
int spec_label(int sg, const char *k_label, int64_t k_label_length, int flags, int has_gid, int64_t gid, int has_source,
               const int64_t *source_words, int has_carrier, const int64_t *carrier_words, int case_count, const int64_t *case_words,
               int display_count, const int64_t *display_words, int has_selected_gid, int64_t selected_gid, int has_selected_source,
               const int64_t *selected_source_words, int request_count, const int64_t *request_words, const char *fallback,
               int64_t fallback_length)
{
    reset_text();
    if (load()) return MT_TABLES;
    rat_overflow = 0;
    record4 source = record_of(has_source, source_words), carrier = record_of(has_carrier, carrier_words);
    record4 selected_source = record_of(has_selected_source, selected_source_words);
    params spec_case, selected_display = {-1, {{0, 1}}}, request;
    read_params(case_count, case_words, &spec_case);
    read_params(request_count, request_words, &request);
    if ((flags & SL_COUPLED_TRUE) && display_count >= 0) read_params(display_count, display_words, &selected_display);

    int status = 0, done = 0;
    if (!(flags & SL_COUPLED) && spec_case.count >= 0) {
        R.k = spec_case;
        done = 1;
    }
    if (!done && has_gid) {
        const record4 *chosen = (flags & SL_MAGNETIC_FALSE) && carrier.has ? &carrier : &source;
        if (chosen->has) {
            status = record_params(gid, chosen, &R.k);
            if (status) return status;
            done = R.k.count >= 0;
        }
    }
    if (!done && spec_case.count >= 0) {
        R.k = spec_case;
        done = 1;
    }
    if (!done && !(flags & SL_LABEL_MATCHES)) {
        R.k.count = 0;
        done = 1;
    }
    if (!done && (flags & SL_PRIMARY) && has_selected_gid) {
        status = record_params(selected_gid, &selected_source, &R.k);
        if (status) return status;
        done = R.k.count >= 0;
    }
    if (!done) {
        if (request_count < 0) return MT_REQUEST;
        R.k = request;
    }

    int source_case = 0;
    if (selected_display.count < 0 && has_gid && gid > 0 && source.has && source.v[3]) {
        int bad = 0;
        int ssg = (int)at(&T.irr_sg, gid - 1, &bad), kslot = (int)at(&T.irr_k, gid - 1, &bad);
        params source_params;
        if (!bad && ssg == sg && label_is(geo_sg_lattice(ssg), kslot, k_label, k_label_length) &&
            !record_params(gid, &source, &source_params) && source_params.count == R.k.count) {
            source_case = 1;
            for (int i = 0; i < R.k.count; i++) source_case &= rat_eq(source_params.v[i], R.k.v[i]);
        }
    }
    if ((flags & SL_PRIMARY) && !(flags & SL_COUPLED) && !source_case) {
        put_n(fallback, (size_t)fallback_length);
        return rat_overflow ? MT_RANGE : B.failed ? MT_MEMORY : 0;
    }
    if ((flags & SL_RECIPROCAL) && has_gid) {
        status = display_try(sg, k_label, k_label_length, flags, gid, &source, &spec_case, &selected_display, source_case);
        if (status == 0 || status == 2) {
            if (status == 2) {
                reset_text();
                put_n(fallback, (size_t)fallback_length);
            }
            return rat_overflow ? MT_RANGE : B.failed ? MT_MEMORY : 0;
        }
        if (status != 1) return status;
        reset_text();
    }
    rat visible[3];
    if (case_vector(sg, k_label, k_label_length, &R.k, visible))
        put_n(fallback, (size_t)fallback_length);
    else
        put_vector(visible);
    return rat_overflow ? MT_RANGE : B.failed ? MT_MEMORY : 0;
}

const char *spec_label_text(void) { return B.data; }

/* The k parameters of the last spec_label as int128 rational words. */
int spec_label_params(int64_t *words)
{
    for (int i = 0; i < R.k.count; i++) {
        write_i128(words + 4 * i, R.k.v[i].n);
        write_i128(words + 4 * i + 2, R.k.v[i].d);
    }
    return R.k.count;
}

/* 0 when absent */
int little_star_arm_width(int64_t gid)
{
    if (load()) return -MT_TABLES;
    int bad = 0;
    int64_t full_dim = at(&T.full_dim, gid - 1, &bad);
    int sg = (int)at(&T.irr_sg, gid - 1, &bad), kslot = (int)at(&T.irr_k, gid - 1, &bad);
    if (bad) return -MT_INDEX;
    int64_t count = at(&T.star_count, (int64_t)(sg - 1) * 27 + kslot - 1, &bad);
    if (bad) return -MT_INDEX;
    if (count <= 0 || full_dim % count) return 0;
    return (int)(full_dim / count);
}

/* _should_emit_all_source_opd_rows.  rows: row_count x width source rows;
 * free_count: int(source_free_count or 0), has_free: it is not None; flags
 * carry the request's coupled and the spec's primary; direction:
 * direction_rows x direction_width (direction_rows -1 for a missing key, -2
 * for a value that is not a matrix).  Returns 0/1 or -MT_. */
int emit_all_source_rows(int row_count, int width, const double *rows, int has_free, int free_count, int flags, int64_t old_id, int64_t gid,
                         int little_type, int direction_rows, int direction_width, const double *direction)
{
    if (row_count <= 0) return 0;
    int rank, status;
    if ((flags & SL_PRIMARY) && free_count == row_count) {
        if ((status = sp_matrix_rank(rows, row_count, width, 1e-10, &rank))) return -MT_RANK;
        if (rank == row_count) return 1;
    }
    if (old_id <= 0) return !(flags & SL_PRIMARY) && row_count > 1;
    if (!(flags & SL_COUPLED) && !(flags & SL_PRIMARY) && gid > 0) {
        if (load()) return -MT_TABLES;
        int bad = 0;
        int64_t little_old = at(&T.irr_old, gid - 1, &bad), type = at(&T.irr_type, gid - 1, &bad);
        if (bad) return -MT_INDEX;
        if (little_old == old_id && type == 3) {
            if (direction_rows == -1) return -MT_KEY;
            int complete = row_count >= 4 && row_count % 2 == 0 && width >= row_count && direction_rows == width &&
                           direction_width == row_count;
            for (int r = 0; r < row_count && complete; r++)
                for (int c = 0; c < width && complete; c++) complete = direction[c * row_count + r] == rows[r * width + c];
            if (complete) {
                if ((status = sp_matrix_rank(rows, row_count, width, 1e-10, &rank))) return -MT_RANK;
                if (rank == row_count && (!has_free || free_count == row_count)) return 1;
            }
        }
    }
    return (flags & SL_COUPLED) && little_type == 3 && row_count > 1;
}

/* the case k parameters of a Source k record for a positive gid and a four-int record: the count (-1 when absent) with out as int128 words, or
 * -100 - MT_. */
int source_record_case_k_params(int64_t gid, const int64_t *record, int64_t *out)
{
    if (load()) return -100 - MT_TABLES;
    rat_overflow = 0;
    record4 source = record_of(1, record);
    params values;
    int status = source_record_params(gid, &source, &values);
    if (status) return -100 - status;
    if (rat_overflow) return -100 - MT_RANGE;
    for (int i = 0; i < values.count; i++) {
        write_i128(out + 4 * i, values.v[i].n);
        write_i128(out + 4 * i + 2, values.v[i].d);
    }
    return values.count;
}
