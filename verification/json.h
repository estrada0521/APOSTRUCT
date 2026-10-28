/* A cursor over JSON text read in place, shared by the native checks under verification/. */
#ifndef VERIFICATION_JSON_H
#define VERIFICATION_JSON_H
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *p, *end;
} cur;

static void ws(cur *c)
{
    while (c->p < c->end && (*c->p == ' ' || *c->p == '\n' || *c->p == '\r' || *c->p == '\t')) c->p++;
}

static int skip(cur *c);

static int skip_string(cur *c)
{
    if (c->p >= c->end || *c->p != '"') return -1;
    for (c->p++; c->p < c->end; c->p++) {
        if (*c->p == '\\') c->p++;
        else if (*c->p == '"') {
            c->p++;
            return 0;
        }
    }
    return -1;
}

static int skip(cur *c)
{
    ws(c);
    if (c->p >= c->end) return -1;
    char ch = *c->p;
    if (ch == '"') return skip_string(c);
    if (ch == '{' || ch == '[') {
        char close = ch == '{' ? '}' : ']';
        c->p++;
        ws(c);
        if (c->p < c->end && *c->p == close) {
            c->p++;
            return 0;
        }
        for (;;) {
            if (ch == '{') {
                if (skip_string(c)) return -1;
                ws(c);
                if (c->p >= c->end || *c->p != ':') return -1;
                c->p++;
            }
            if (skip(c)) return -1;
            ws(c);
            if (c->p >= c->end) return -1;
            if (*c->p == ',') {
                c->p++;
                ws(c);
                continue;
            }
            if (*c->p == close) {
                c->p++;
                return 0;
            }
            return -1;
        }
    }
    while (c->p < c->end && *c->p != ',' && *c->p != ']' && *c->p != '}' && *c->p != ' ' && *c->p != '\n' && *c->p != '\r' && *c->p != '\t') c->p++;
    return 0;
}

/* object iteration: key and klen set, c at the value; returns 1 per member, 0 at the end, -1 on error */
static int obj_begin(cur *c)
{
    ws(c);
    if (c->p >= c->end || *c->p != '{') return -1;
    c->p++;
    return 0;
}

static int obj_next(cur *c, const char **key, int *klen, int *first)
{
    ws(c);
    if (c->p < c->end && *c->p == '}') {
        c->p++;
        return 0;
    }
    if (!*first) {
        if (c->p >= c->end || *c->p != ',') return -1;
        c->p++;
        ws(c);
    }
    *first = 0;
    const char *s = c->p + 1;
    if (skip_string(c)) return -1;
    *key = s;
    *klen = (int)(c->p - 1 - s);
    ws(c);
    if (c->p >= c->end || *c->p != ':') return -1;
    c->p++;
    ws(c);
    return 1;
}

static int arr_begin(cur *c)
{
    ws(c);
    if (c->p >= c->end || *c->p != '[') return -1;
    c->p++;
    return 0;
}

static int arr_next(cur *c, int *first)
{
    ws(c);
    if (c->p < c->end && *c->p == ']') {
        c->p++;
        return 0;
    }
    if (!*first) {
        if (c->p >= c->end || *c->p != ',') return -1;
        c->p++;
        ws(c);
    }
    *first = 0;
    return 1;
}

static int keyis(const char *key, int klen, const char *want) { return (int)strlen(want) == klen && !memcmp(key, want, klen); }

static int is_null(cur *c)
{
    ws(c);
    return c->end - c->p >= 4 && !memcmp(c->p, "null", 4);
}

/* the correctly rounded value of a numeral [+-]digits[.digits][(e|E)[+-]digits], 0 when it is not one: one exact
 * division or product when the digits fit in 2^53 and the power of ten is exact (Clinger), strtod otherwise */
static int decimal(const char *s, size_t n, double *out)
{
    static const double pow10[23] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                     1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
    size_t i = 0;
    int neg = 0, digits = 0, scale = 0, exp = 0, eneg = 0, any = 0;
    uint64_t m = 0;
    if (i < n && (s[i] == '+' || s[i] == '-')) neg = s[i++] == '-';
    for (; i < n && s[i] >= '0' && s[i] <= '9'; i++, any = 1)
        if (m || s[i] != '0') m = 10 * m + (uint64_t)(s[i] - '0'), digits++;
    if (i < n && s[i] == '.')
        for (i++; i < n && s[i] >= '0' && s[i] <= '9'; i++, any = 1) {
            if (m || s[i] != '0') m = 10 * m + (uint64_t)(s[i] - '0'), digits++;
            scale++;
        }
    if (!any) return 0;
    if (i < n && (s[i] == 'e' || s[i] == 'E')) {
        size_t j = i + 1;
        if (j < n && (s[j] == '+' || s[j] == '-')) eneg = s[j++] == '-';
        if (!(j < n && s[j] >= '0' && s[j] <= '9')) return 0;
        for (; j < n && s[j] >= '0' && s[j] <= '9'; j++)
            if (exp < 100000) exp = 10 * exp + (s[j] - '0');
        i = j;
    }
    if (i != n) return 0;
    int k = (eneg ? -exp : exp) - scale;
    if (digits <= 19 && m <= (1ULL << 53) && k >= -22 && k <= 22) {
        double v = k < 0 ? (double)m / pow10[-k] : (double)m * pow10[k];
        *out = neg ? -v : v;
        return 1;
    }
    char small[64], *buf = n < sizeof small ? small : malloc(n + 1);
    memcpy(buf, s, n);
    buf[n] = 0;
    *out = strtod(buf, 0);
    if (buf != small) free(buf);
    return 1;
}

static int number(cur *c, double *out)
{
    ws(c);
    const char *s = c->p;
    while (c->p < c->end && ((*c->p >= '0' && *c->p <= '9') || *c->p == '.' || *c->p == '-' || *c->p == '+' || *c->p == 'e' || *c->p == 'E')) c->p++;
    return decimal(s, (size_t)(c->p - s), out) && isfinite(*out) ? 0 : -1;
}

static int string_value(cur *c, const char **s, int *len)
{
    ws(c);
    const char *start = c->p + 1;
    if (skip_string(c)) return -1;
    *s = start;
    *len = (int)(c->p - 1 - start);
    return 0;
}

static inline int string_dup(cur *c, char **out)
{
    const char *s;
    int len;
    if (string_value(c, &s, &len)) return -1;
    char *o = malloc(len + 1), *q = o;
    for (int i = 0; i < len; i++) {
        if (s[i] != '\\') {
            *q++ = s[i];
            continue;
        }
        if (++i >= len) break;
        char e = s[i];
        if (e == 'n') *q++ = '\n';
        else if (e == 't') *q++ = '\t';
        else if (e == 'r') *q++ = '\r';
        else if (e == 'b') *q++ = '\b';
        else if (e == 'f') *q++ = '\f';
        else if (e == 'u' && i + 4 < len) {
            unsigned u = (unsigned)strtoul((char[5]){s[i + 1], s[i + 2], s[i + 3], s[i + 4], 0}, 0, 16);
            i += 4;
            if (u >= 0xd800 && u < 0xdc00 && i + 6 < len && s[i + 1] == '\\' && s[i + 2] == 'u') {
                unsigned v = (unsigned)strtoul((char[5]){s[i + 3], s[i + 4], s[i + 5], s[i + 6], 0}, 0, 16);
                u = 0x10000 + ((u - 0xd800) << 10) + (v - 0xdc00);
                i += 6;
            }
            if (u < 0x80) *q++ = (char)u;
            else if (u < 0x800) *q++ = (char)(0xc0 | u >> 6), *q++ = (char)(0x80 | (u & 0x3f));
            else if (u < 0x10000) *q++ = (char)(0xe0 | u >> 12), *q++ = (char)(0x80 | (u >> 6 & 0x3f)), *q++ = (char)(0x80 | (u & 0x3f));
            else *q++ = (char)(0xf0 | u >> 18), *q++ = (char)(0x80 | (u >> 12 & 0x3f)), *q++ = (char)(0x80 | (u >> 6 & 0x3f)), *q++ = (char)(0x80 | (u & 0x3f));
        } else *q++ = e;
    }
    *q = 0;
    *out = o;
    return 0;
}

#endif
