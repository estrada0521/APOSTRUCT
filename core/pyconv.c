#include "pyconv.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int py_space(unsigned char c) { return c == ' ' || (c >= '\t' && c <= '\r') || (c >= 0x1c && c <= 0x1f); }

static void py_strip(const char **start, const char **end)
{
    while (*start < *end && py_space((unsigned char)**start)) (*start)++;
    while (*end > *start && py_space((unsigned char)(*end)[-1])) (*end)--;
}

static int digit(char c) { return c >= '0' && c <= '9'; }

/* Underscores only between digits, then dropped (Python's number grammar). */
static int drop_underscores(const char *start, const char *end, char *out, size_t capacity)
{
    size_t n = 0;
    for (const char *c = start; c < end; c++) {
        if (*c == '_') {
            if (c == start || c + 1 >= end || !digit(c[-1]) || !digit(c[1])) return 1;
            continue;
        }
        if (n + 1 >= capacity) return 1;
        out[n++] = *c;
    }
    out[n] = 0;
    return 0;
}

static int casefold_is(const char *text, const char *word)
{
    for (; *word; text++, word++)
        if (tolower((unsigned char)*text) != *word) return 0;
    return !*text;
}

int py_float_of_text(const char *text, int64_t length, double *out)
{
    const char *start = text, *end = text + length;
    py_strip(&start, &end);
    char buffer[512];
    if (start == end || drop_underscores(start, end, buffer, sizeof buffer)) return 1;
    const char *c = buffer;
    int negative = 0;
    if (*c == '+' || *c == '-') negative = *c++ == '-';
    if (casefold_is(c, "inf") || casefold_is(c, "infinity")) {
        *out = negative ? -INFINITY : INFINITY;
        return 0;
    }
    if (casefold_is(c, "nan")) {
        *out = negative ? -NAN : NAN;
        return 0;
    }
    int digits = 0;
    const char *d = c;
    while (digit(*d)) {
        d++;
        digits++;
    }
    if (*d == '.') {
        d++;
        while (digit(*d)) {
            d++;
            digits++;
        }
    }
    if (!digits) return 1;
    if (*d == 'e' || *d == 'E') {
        d++;
        if (*d == '+' || *d == '-') d++;
        if (!digit(*d)) return 1;
        while (digit(*d)) d++;
    }
    if (*d) return 1;
    *out = strtod(buffer, 0);
    return 0;
}

int py_fraction_of_text(const char *text, int64_t length, rat *out)
{
    const char *c = text, *end = text + length;
    py_strip(&c, &end);
    int negative = 0;
    if (c < end && (*c == '-' || *c == '+')) negative = *c++ == '-';
    if (!(c < end && (digit(*c) || (*c == '.' && c + 1 < end && digit(c[1]))))) return 1;
    rat_overflow = 0;
    i128 numerator = 0, denominator = 1;
    int64_t count;
    char buffer[512];
#define RUN(at, value, n)                                                                 \
    do {                                                                                  \
        (n) = 0;                                                                          \
        while ((at) < end) {                                                              \
            if (digit(*(at))) {                                                           \
                (value) = rat_add_i(rat_mul_i((value), 10), *(at) - '0');                 \
                (n)++;                                                                    \
                (at)++;                                                                   \
            } else if (*(at) == '_' && (at) + 1 < end && digit((at)[1]) && (n)) {         \
                (at)++;                                                                   \
            } else {                                                                      \
                break;                                                                    \
            }                                                                             \
        }                                                                                 \
    } while (0)
    (void)buffer;
    RUN(c, numerator, count);
    const char *mark = c;
    while (mark < end && py_space((unsigned char)*mark)) mark++;
    if (mark < end && *mark == '/') {
        mark++;
        while (mark < end && py_space((unsigned char)*mark)) mark++;
        i128 value = 0;
        RUN(mark, value, count);
        if (!count) return 1;
        c = mark;
        if (!value) return 1;
        denominator = value;
    } else {
        if (c < end && *c == '.') {
            c++;
            i128 fraction = 0;
            RUN(c, fraction, count);
            for (int64_t i = 0; i < count; i++) {
                numerator = rat_mul_i(numerator, 10);
                denominator = rat_mul_i(denominator, 10);
            }
            numerator = rat_add_i(numerator, fraction);
        }
        if (c < end && (*c == 'e' || *c == 'E')) {
            const char *e = c + 1;
            int exponent_negative = 0;
            if (e < end && (*e == '-' || *e == '+')) exponent_negative = *e++ == '-';
            i128 exponent = 0;
            RUN(e, exponent, count);
            if (!count) return 1;
            if (exponent > 400) return 2;
            for (i128 i = 0; i < exponent; i++) {
                if (exponent_negative) denominator = rat_mul_i(denominator, 10);
                else if (denominator % 10 == 0) denominator /= 10;
                else numerator = rat_mul_i(numerator, 10);
            }
            c = e;
        }
    }
#undef RUN
    while (c < end && py_space((unsigned char)*c)) c++;
    if (c != end) return 1;
    if (rat_overflow) return 2;
    *out = rat_make(negative ? -numerator : numerator, denominator);
    return rat_overflow ? 2 : 0;
}

int py_fraction_of_float(double value, rat *out)
{
    if (!isfinite(value)) return 1;
    char text[40];
    int n = py_float_repr(value, text);
    return py_fraction_of_text(text, n, out);
}

int py_i128_text(i128 value, char *out)
{
    char digits[48];
    int n = 0;
    u128 magnitude = value < 0 ? (u128)0 - (u128)value : (u128)value;
    do {
        digits[n++] = (char)('0' + (int)(magnitude % 10));
        magnitude /= 10;
    } while (magnitude);
    int length = 0;
    if (value < 0) out[length++] = '-';
    while (n) out[length++] = digits[--n];
    out[length] = 0;
    return length;
}

int py_rat_text(rat value, char *out)
{
    int length = py_i128_text(value.n, out);
    if (value.d != 1) {
        out[length++] = '/';
        length += py_i128_text(value.d, out + length);
    }
    return length;
}

int py_float_of_fraction(rat value, double *out)
{
    if (value.d <= 0) return 1;
    if (!value.n) {
        *out = 0.0;
        return 0;
    }
    int negative = value.n < 0;
    u128 a = negative ? (u128)0 - (u128)value.n : (u128)value.n, b = (u128)value.d;
    u128 q = a / b, r = a % b;
    int exponent = 0, sticky = 0;
    const u128 low = (u128)1 << 62, high = (u128)1 << 63;
    while (q >= high) {
        sticky |= (int)(q & 1);
        q >>= 1;
        exponent++;
    }
    while (q < low) {
        r <<= 1;
        q <<= 1;
        if (r >= b) {
            r -= b;
            q |= 1;
        }
        exponent--;
    }
    sticky |= r != 0;
    uint64_t mantissa = (uint64_t)(q >> 10), rest = (uint64_t)(q & 0x3FF);
    if (rest > 0x200 || (rest == 0x200 && (sticky || (mantissa & 1)))) mantissa++;
    double magnitude = ldexp((double)mantissa, exponent + 10);
    *out = negative ? -magnitude : magnitude;
    return 0;
}

static i128 floor_div(i128 a, i128 b)
{
    i128 q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) q--;
    return q;
}

rat py_limit_denominator(rat value, i128 max_denominator)
{
    if (value.d <= max_denominator) return value;
    i128 p0 = 0, q0 = 1, p1 = 1, q1 = 0, n = value.n, d = value.d;
    for (;;) {
        i128 a = floor_div(n, d);
        i128 q2 = q0 + a * q1;
        if (q2 > max_denominator) break;
        i128 p2 = p0 + a * p1;
        p0 = p1;
        q0 = q1;
        p1 = p2;
        q1 = q2;
        i128 r = n - a * d;
        n = d;
        d = r;
    }
    i128 k = floor_div(max_denominator - q0, q1);
    if (2 * d * (q0 + k * q1) <= value.d) return (rat){p1, q1};
    return (rat){p0 + k * p1, q0 + k * q1};
}

int py_float_of(const jvalue *v, double *out)
{
    if (!v) return 1;
    switch (v->kind) {
    case J_INT: *out = (double)v->u.i; return 0;
    case J_FLOAT: *out = v->u.f; return 0;
    case J_TRUE: *out = 1.0; return 0;
    case J_FALSE: *out = 0.0; return 0;
    case J_STRING: return py_float_of_text(v->u.s, v->count, out);
    default: return 1;
    }
}

int py_int_of(const jvalue *v, int64_t *out)
{
    if (!v) return 1;
    switch (v->kind) {
    case J_INT: *out = v->u.i; return 0;
    case J_TRUE: *out = 1; return 0;
    case J_FALSE: *out = 0; return 0;
    case J_FLOAT:
        if (!isfinite(v->u.f) || fabs(v->u.f) >= 9.2e18) return 1;
        *out = (int64_t)v->u.f;
        return 0;
    case J_STRING: {
        const char *start = v->u.s, *end = v->u.s + v->count;
        py_strip(&start, &end);
        char buffer[64];
        if (start == end || drop_underscores(start, end, buffer, sizeof buffer)) return 1;
        const char *c = buffer;
        if (*c == '+' || *c == '-') c++;
        if (!*c) return 1;
        for (const char *d = c; *d; d++)
            if (!digit(*d)) return 1;
        errno = 0;
        long long value = strtoll(buffer, 0, 10);
        if (errno) return 1;
        *out = value;
        return 0;
    }
    default: return 1;
    }
}

void py_str_of(jbuf *b, const jvalue *v)
{
    if (!v) {
        jb_lit(b, "None");
        return;
    }
    char text[40];
    switch (v->kind) {
    case J_NULL: jb_lit(b, "None"); break;
    case J_TRUE: jb_lit(b, "True"); break;
    case J_FALSE: jb_lit(b, "False"); break;
    case J_INT: jb_int(b, v->u.i); break;
    case J_FLOAT: jb_raw(b, text, (size_t)py_float_repr(v->u.f, text)); break;
    case J_STRING: jb_raw(b, v->u.s, (size_t)v->count); break;
    default: b->failed = 1; break;
    }
}

int float_cif_number(const jvalue *value, double *out)
{
    if (!jtruthy(value)) return 1;
    jbuf text = {0};
    py_str_of(&text, value);
    if (text.failed || !text.data) {
        jb_free(&text);
        return 1;
    }
    const char *start = text.data, *end = text.data + text.used;
    py_strip(&start, &end);
    int status = 1;
    if (start < end && !(end - start == 1 && (*start == '.' || *start == '?'))) {
        char *cleaned = malloc((size_t)(end - start) + 1);
        size_t n = 0;
        if (cleaned) {
            for (const char *c = start; c < end;) {
                if (*c == '(') {
                    const char *close = c + 1;
                    while (close < end && *close != ')') close++;
                    if (close < end) {
                        c = close + 1;
                        continue;
                    }
                }
                cleaned[n++] = *c++;
            }
            rat fraction;
            if (!py_float_of_text(cleaned, (int64_t)n, out)) status = 0;
            else if (!py_fraction_of_text(cleaned, (int64_t)n, &fraction) && !py_float_of_fraction(fraction, out)) status = 0;
            free(cleaned);
        }
    }
    jb_free(&text);
    return status;
}

/* ------------------------------------------------------------------ gemmi */

static const char *skip_space(const char *p, const char *end)
{
    while (p < end && (*p == ' ' || *p == '\t' || *p == '_')) ++p;
    return p;
}

enum { DEN = 24 };

static int interpret_letter(char c, char *notation, int *index)
{
    static const signed char table[] = {('a' & ~3) | 0, ('a' & ~3) | 1, ('a' & ~3) | 2, 0, 0, 0, 0, ('h' & ~3) | 0, 0, 0,
                                        ('h' & ~3) | 1, ('h' & ~3) | 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        ('x' & ~3) | 0, ('x' & ~3) | 1, ('x' & ~3) | 2};
    size_t idx = (size_t)(unsigned char)((c | 0x20) - 'a');
    if (idx >= sizeof table || table[idx] == 0) return 1;
    int value = table[idx], detected = value & ~3;
    if ((*notation | 0x20) == ' ') *notation = (char)detected;
    else if (((*notation | 0x20) & ~3) != detected) return 1;
    *index = value & 3;
    return 0;
}

/* strtol on [p, end): leading isspace, sign, digits; *stop at the first unread. */
static long bounded_strtol(const char *p, const char *end, const char **stop)
{
    const char *c = p;
    while (c < end && isspace((unsigned char)*c)) c++;
    int negative = 0;
    if (c < end && (*c == '+' || *c == '-')) negative = *c++ == '-';
    if (c >= end || !digit(*c)) {
        *stop = p;
        return 0;
    }
    long value = 0;
    while (c < end && digit(*c)) {
        if (value < LONG_MAX / 10) value = value * 10 + (*c - '0');
        else value = LONG_MAX;
        c++;
    }
    *stop = c;
    return negative ? -value : value;
}

static int parse_part(const char *s, const char *end, char *notation, int r[4])
{
    r[0] = r[1] = r[2] = r[3] = 0;
    int num = DEN;
    const char *c = s;
    while ((c = skip_space(c, end)) < end) {
        if (*c == '+' || *c == '-') {
            num = *c == '+' ? DEN : -DEN;
            c = skip_space(c + 1, end);
        }
        if (num == 0) return 1;
        int r_idx, den = 1;
        double fract = 0;
        if (c < end && (digit(*c) || *c == '.')) {
            const char *endptr;
            int n = (int)bounded_strtol(c, end, &endptr);
            if (endptr < end && *endptr == '.') {
                fract = n;
                for (double denom = 0.1; ++endptr < end && digit(*endptr); denom *= 0.1) fract = fma((double)(*endptr - '0'), denom, fract);
                double rounded = round(fract * num);
                if (fabs(fma(-fract, (double)num, rounded)) > 0.05) return 1;
                num = (int)rounded;
            } else {
                num *= n;
            }
            if (endptr < end && *endptr == '/') den = (int)bounded_strtol(endptr + 1, end, &endptr);
            if (endptr < end && *endptr == '*') {
                c = skip_space(endptr + 1, end);
                if (c >= end || interpret_letter(*c, notation, &r_idx)) return 1;
                ++c;
            } else {
                c = endptr;
                r_idx = 3;
            }
        } else {
            if (c >= end || interpret_letter(*c, notation, &r_idx)) return 1;
            c = skip_space(c + 1, end);
            if (c < end && *c == '/') {
                const char *endptr;
                den = (int)bounded_strtol(c + 1, end, &endptr);
                c = endptr;
            }
        }
        if (den != 1) {
            if (den <= 0 || DEN % den != 0 || fract != 0) return 1;
            num /= den;
        }
        r[r_idx] += num;
        num = 0;
    }
    return num != 0;
}

int gemmi_parse_triplet(const char *text, int64_t length, int rot[3][3], int tran[3], char *notation)
{
    const char *end = text + length, *commas[2];
    int count = 0;
    for (const char *c = text; c < end; c++)
        if (*c == ',') {
            if (count < 2) commas[count] = c;
            count++;
        }
    if (count != 2) return 1;
    for (const char *c = text; c < end; c++)
        if (!*c) return 1;
    *notation = (char)((' ' | 0x20) & ~3);
    int parts[3][4];
    const char *bounds[4] = {text, commas[0] + 1, commas[1] + 1, end + 1};
    for (int k = 0; k < 3; k++)
        if (parse_part(bounds[k], bounds[k + 1] - 1, notation, parts[k])) return 1;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) rot[i][j] = parts[i][j];
        tran[i] = parts[i][3];
    }
    if (*notation == ('h' & ~3)) {
        if (tran[0] || tran[1] || tran[2]) return 1;
        for (int i = 0; i < 3; i++)
            for (int j = i + 1; j < 3; j++) {
                int t = rot[i][j];
                rot[i][j] = rot[j][i];
                rot[j][i] = t;
            }
    }
    return 0;
}

void gemmi_apply_to_xyz(const int rot[3][3], const int tran[3], const double *xyz, double *out)
{
    for (int i = 0; i < 3; i++)
        out[i] = (fma((double)rot[i][2], xyz[2], fma((double)rot[i][0], xyz[0], (double)rot[i][1] * xyz[1])) + (double)tran[i]) / DEN;
}

double py_sum(const double *values, int count)
{
    if (count <= 0) return 0.0;
    double hi = 0.0 + values[0], lo = 0.0;
    for (int i = 1; i < count; i++) {
        double x = values[i], t = hi + x;
        if (fabs(hi) >= fabs(x)) lo += (hi - t) + x;
        else lo += (x - t) + hi;
        hi = t;
    }
    if (lo && isfinite(lo)) return hi + lo;
    return hi;
}

static void repr_string(jbuf *b, const char *s, int64_t n)
{
    int single = 0, dbl = 0;
    for (int64_t i = 0; i < n; i++) {
        single |= s[i] == '\'';
        dbl |= s[i] == '"';
    }
    char quote = single && !dbl ? '"' : '\'';
    jb_raw(b, &quote, 1);
    for (int64_t i = 0; i < n; i++) {
        unsigned char c = (unsigned char)s[i];
        char escape[8];
        if (c == '\\') jb_lit(b, "\\\\");
        else if (c == (unsigned char)quote) {
            escape[0] = '\\';
            escape[1] = quote;
            jb_raw(b, escape, 2);
        } else if (c == '\n') jb_lit(b, "\\n");
        else if (c == '\r') jb_lit(b, "\\r");
        else if (c == '\t') jb_lit(b, "\\t");
        else if (c < 0x20 || c == 0x7f) {
            snprintf(escape, sizeof escape, "\\x%02x", c);
            jb_lit(b, escape);
        } else jb_raw(b, (const char *)&c, 1);
    }
    jb_raw(b, &quote, 1);
}

void py_repr_of(jbuf *b, const jvalue *v)
{
    if (!v) {
        jb_lit(b, "None");
        return;
    }
    char text[40];
    switch (v->kind) {
    case J_NULL: jb_lit(b, "None"); break;
    case J_TRUE: jb_lit(b, "True"); break;
    case J_FALSE: jb_lit(b, "False"); break;
    case J_INT: jb_int(b, v->u.i); break;
    case J_FLOAT: jb_raw(b, text, (size_t)py_float_repr(v->u.f, text)); break;
    case J_STRING: repr_string(b, v->u.s, v->count); break;
    case J_ARRAY:
        jb_lit(b, "[");
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_lit(b, ", ");
            py_repr_of(b, v->u.items + i);
        }
        jb_lit(b, "]");
        break;
    case J_OBJECT:
        jb_lit(b, "{");
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_lit(b, ", ");
            py_repr_of(b, v->u.members + 2 * i);
            jb_lit(b, ": ");
            py_repr_of(b, v->u.members + 2 * i + 1);
        }
        jb_lit(b, "}");
        break;
    }
}

void py_repr_text(jbuf *b, const char *text, int64_t length) { repr_string(b, text, length); }

double py_pow(double x, double y)
{
    static double (*volatile libm_pow)(double, double) = pow;
    return libm_pow(x, y);
}
