#include "json.h"

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

/* ------------------------------------------------------------------ big blocks */

/* Blocks of BIG bytes and more are mapped directly: macOS malloc keeps freed
 * large blocks resident, so a long-lived process grew by every large result. */
enum { BIG = 8 << 20 };

static void *big_alloc(size_t size)
{
    if (size < BIG) return malloc(size);
    void *p = mmap(0, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    return p == MAP_FAILED ? 0 : p;
}

static void big_free(void *p, size_t size)
{
    if (!p) return;
    if (size < BIG) free(p);
    else munmap(p, size);
}

/* ------------------------------------------------------------------ arena */

typedef struct block {
    struct block *next;
    size_t used, cap;
    _Alignas(16) unsigned char data[];
} block;

struct jdoc {
    block *blocks;
    jvalue root;
};

static void *arena_alloc(jdoc *doc, size_t size)
{
    size = (size + 15) & ~(size_t)15;
    block *b = doc->blocks;
    if (!b || b->used + size > b->cap) {
        size_t cap = size > 65536 ? size : 65536;
        b = big_alloc(sizeof(block) + cap);
        if (!b) return 0;
        b->next = doc->blocks;
        b->used = 0;
        b->cap = cap;
        doc->blocks = b;
    }
    void *out = b->data + b->used;
    b->used += size;
    return out;
}

void jdoc_free(jdoc *doc)
{
    if (!doc) return;
    for (block *b = doc->blocks; b;) {
        block *next = b->next;
        big_free(b, sizeof(block) + b->cap);
        b = next;
    }
    free(doc);
}

const jvalue *jdoc_root(const jdoc *doc) { return &doc->root; }

/* ----------------------------------------------------------------- parser */

typedef struct {
    const char *at, *end;
    jdoc *doc;
    int depth;
} parser;

static void ws(parser *p)
{
    while (p->at < p->end && (*p->at == ' ' || *p->at == '\t' || *p->at == '\n' || *p->at == '\r')) p->at++;
}

static int hex4(const char *s, unsigned *out)
{
    unsigned v = 0;
    for (int i = 0; i < 4; i++) {
        char c = s[i];
        v <<= 4;
        if (c >= '0' && c <= '9') v |= (unsigned)(c - '0');
        else if (c >= 'a' && c <= 'f') v |= (unsigned)(c - 'a' + 10);
        else if (c >= 'A' && c <= 'F') v |= (unsigned)(c - 'A' + 10);
        else return 1;
    }
    *out = v;
    return 0;
}

static size_t utf8(unsigned code, char *out)
{
    if (code < 0x80) {
        out[0] = (char)code;
        return 1;
    }
    if (code < 0x800) {
        out[0] = (char)(0xC0 | (code >> 6));
        out[1] = (char)(0x80 | (code & 0x3F));
        return 2;
    }
    if (code < 0x10000) {
        out[0] = (char)(0xE0 | (code >> 12));
        out[1] = (char)(0x80 | ((code >> 6) & 0x3F));
        out[2] = (char)(0x80 | (code & 0x3F));
        return 3;
    }
    out[0] = (char)(0xF0 | (code >> 18));
    out[1] = (char)(0x80 | ((code >> 12) & 0x3F));
    out[2] = (char)(0x80 | ((code >> 6) & 0x3F));
    out[3] = (char)(0x80 | (code & 0x3F));
    return 4;
}

static int parse_string(parser *p, jvalue *out)
{
    p->at++;
    const char *start = p->at;
    int escaped = 0;
    while (p->at < p->end && *p->at != '"') {
        if (*p->at == '\\') {
            escaped = 1;
            p->at++;
        }
        p->at++;
    }
    if (p->at >= p->end) return 1;
    size_t raw = (size_t)(p->at - start);
    char *s = arena_alloc(p->doc, raw + 1);
    if (!s) return 1;
    size_t n = 0;
    if (!escaped) {
        memcpy(s, start, raw);
        n = raw;
    } else {
        for (const char *c = start; c < p->at; c++) {
            if (*c != '\\') {
                s[n++] = *c;
                continue;
            }
            c++;
            switch (*c) {
            case '"': s[n++] = '"'; break;
            case '\\': s[n++] = '\\'; break;
            case '/': s[n++] = '/'; break;
            case 'b': s[n++] = '\b'; break;
            case 'f': s[n++] = '\f'; break;
            case 'n': s[n++] = '\n'; break;
            case 'r': s[n++] = '\r'; break;
            case 't': s[n++] = '\t'; break;
            case 'u': {
                unsigned code;
                if (p->at - c < 5 || hex4(c + 1, &code)) return 1;
                c += 4;
                if (code >= 0xD800 && code < 0xDC00 && p->at - c >= 7 && c[1] == '\\' && c[2] == 'u') {
                    unsigned low;
                    if (!hex4(c + 3, &low) && low >= 0xDC00 && low < 0xE000) {
                        code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
                        c += 6;
                    }
                }
                n += utf8(code, s + n);
                break;
            }
            default: return 1;
            }
        }
    }
    s[n] = 0;
    p->at++;
    out->kind = J_STRING;
    out->count = (int64_t)n;
    out->u.s = s;
    return 0;
}

static int parse_value(parser *p, jvalue *out);

typedef struct chunk {
    jvalue values[64];
    struct chunk *next;
} chunk;

/* Collect items in a temporary chain, then copy them into the arena. */
static int parse_items(parser *p, jvalue *out, int object)
{
    char close = object ? '}' : ']';
    p->at++;
    ws(p);
    chunk head, *tail = &head;
    head.next = 0;
    int64_t count = 0, per = object ? 2 : 1;
    int status = 0;
    if (p->at < p->end && *p->at == close) {
        p->at++;
    } else {
        for (;;) {
            ws(p);
            for (int part = 0; part < per; part++) {
                if (count % 64 == 0 && count) {
                    chunk *next = malloc(sizeof(chunk));
                    if (!next) {
                        status = 1;
                        goto done;
                    }
                    next->next = 0;
                    tail->next = next;
                    tail = next;
                }
                jvalue *slot = tail->values + count % 64;
                if (part == 0 && object) {
                    if (p->at >= p->end || *p->at != '"' || parse_string(p, slot)) {
                        status = 1;
                        goto done;
                    }
                    ws(p);
                    if (p->at >= p->end || *p->at != ':') {
                        status = 1;
                        goto done;
                    }
                    p->at++;
                    ws(p);
                } else if (parse_value(p, slot)) {
                    status = 1;
                    goto done;
                }
                count++;
            }
            ws(p);
            if (p->at < p->end && *p->at == ',') {
                p->at++;
                continue;
            }
            if (p->at < p->end && *p->at == close) {
                p->at++;
                break;
            }
            status = 1;
            goto done;
        }
    }
    {
        jvalue *items = arena_alloc(p->doc, sizeof(jvalue) * (size_t)(count ? count : 1));
        if (!items) {
            status = 1;
            goto done;
        }
        chunk *c = &head;
        for (int64_t i = 0; i < count; i++) {
            if (i && i % 64 == 0) c = c->next;
            items[i] = c->values[i % 64];
        }
        out->kind = object ? J_OBJECT : J_ARRAY;
        out->count = count / per;
        if (object) out->u.members = items;
        else out->u.items = items;
    }
done:
    for (chunk *c = head.next; c;) {
        chunk *next = c->next;
        free(c);
        c = next;
    }
    return status;
}

static int parse_number(parser *p, jvalue *out)
{
    const char *start = p->at;
    int is_float = 0;
    if (p->at < p->end && *p->at == '-') p->at++;
    while (p->at < p->end) {
        char c = *p->at;
        if (c >= '0' && c <= '9') {
            p->at++;
        } else if (c == '.' || c == 'e' || c == 'E' || c == '+' || (c == '-' && is_float)) {
            is_float = 1;
            p->at++;
        } else {
            break;
        }
    }
    char buffer[400];
    size_t n = (size_t)(p->at - start);
    if (!n || n >= sizeof buffer) return 1;
    memcpy(buffer, start, n);
    buffer[n] = 0;
    char *stop;
    if (!is_float) {
        errno = 0;
        long long v = strtoll(buffer, &stop, 10);
        if (*stop || errno) return 1;
        out->kind = J_INT;
        out->u.i = v;
        return 0;
    }
    out->kind = J_FLOAT;
    out->u.f = strtod(buffer, &stop);
    return *stop != 0;
}

static int literal(parser *p, const char *word, jkind kind, jvalue *out)
{
    size_t n = strlen(word);
    if ((size_t)(p->end - p->at) < n || memcmp(p->at, word, n)) return 1;
    p->at += n;
    out->kind = kind;
    out->count = 0;
    return 0;
}

static int parse_value(parser *p, jvalue *out)
{
    ws(p);
    if (p->at >= p->end || ++p->depth > 512) return 1;
    int status;
    switch (*p->at) {
    case '{': status = parse_items(p, out, 1); break;
    case '[': status = parse_items(p, out, 0); break;
    case '"': status = parse_string(p, out); break;
    case 't': status = literal(p, "true", J_TRUE, out); break;
    case 'f': status = literal(p, "false", J_FALSE, out); break;
    case 'n': status = literal(p, "null", J_NULL, out); break;
    case 'N':
        status = literal(p, "NaN", J_FLOAT, out);
        if (!status) out->u.f = NAN;
        break;
    case 'I':
        status = literal(p, "Infinity", J_FLOAT, out);
        if (!status) out->u.f = INFINITY;
        break;
    case '-':
        if (p->end - p->at >= 9 && !memcmp(p->at, "-Infinity", 9)) {
            p->at += 9;
            out->kind = J_FLOAT;
            out->u.f = -INFINITY;
            status = 0;
            break;
        }
        status = parse_number(p, out);
        break;
    default: status = parse_number(p, out);
    }
    p->depth--;
    return status;
}

jdoc *jdoc_parse(const char *text, size_t length)
{
    jdoc *doc = calloc(1, sizeof(jdoc));
    if (!doc) return 0;
    parser p = {text, text + length, doc, 0};
    if (parse_value(&p, &doc->root)) {
        jdoc_free(doc);
        return 0;
    }
    ws(&p);
    if (p.at != p.end) {
        jdoc_free(doc);
        return 0;
    }
    return doc;
}

const jvalue *jget(const jvalue *object, const char *key)
{
    if (!object || object->kind != J_OBJECT) return 0;
    size_t n = strlen(key);
    const jvalue *found = 0;
    for (int64_t i = 0; i < object->count; i++) {
        const jvalue *k = object->u.members + 2 * i;
        if ((size_t)k->count == n && !memcmp(k->u.s, key, n)) found = k + 1;
    }
    return found;
}

const jvalue *jat(const jvalue *array, int64_t index)
{
    if (!array || array->kind != J_ARRAY || index < 0 || index >= array->count) return 0;
    return array->u.items + index;
}

int jis(const jvalue *value, jkind kind) { return value && value->kind == kind; }

/* Python truthiness of the decoded value (a missing value is falsy). */
int jtruthy(const jvalue *value)
{
    if (!value) return 0;
    switch (value->kind) {
    case J_NULL:
    case J_FALSE: return 0;
    case J_TRUE: return 1;
    case J_INT: return value->u.i != 0;
    case J_FLOAT: return value->u.f != 0.0;
    default: return value->count != 0;
    }
}

int jstr_is(const jvalue *value, const char *text)
{
    return value && value->kind == J_STRING && (size_t)value->count == strlen(text) && !memcmp(value->u.s, text, (size_t)value->count);
}

/* ----------------------------------------------------------------- writer */

void jb_reset(jbuf *b)
{
    b->used = 0;
    b->failed = 0;
    jb_raw(b, "", 0);
}

void jb_free(jbuf *b)
{
    big_free(b->data, b->cap);
    memset(b, 0, sizeof *b);
}

void jb_raw(jbuf *b, const char *text, size_t length)
{
    if (b->failed) return;
    if (b->used + length + 1 > b->cap) {
        size_t cap = b->cap ? b->cap : 4096;
        while (b->used + length + 1 > cap) cap *= 2;
        char *data;
        if (cap < BIG) {
            data = realloc(b->data, cap);
        } else {
            data = big_alloc(cap);
            if (data && b->data) {
                memcpy(data, b->data, b->used + 1);
                big_free(b->data, b->cap);
            }
        }
        if (!data) {
            b->failed = 1;
            return;
        }
        b->data = data;
        b->cap = cap;
    }
    memcpy(b->data + b->used, text, length);
    b->used += length;
    b->data[b->used] = 0;
}

void jb_lit(jbuf *b, const char *text) { jb_raw(b, text, strlen(text)); }

void jb_int(jbuf *b, int64_t value)
{
    char text[32];
    int n = snprintf(text, sizeof text, "%lld", (long long)value);
    jb_raw(b, text, (size_t)n);
}

static int float_repr(double value, char *out);

static struct {
    uint64_t bits;
    unsigned char length;
    char text[31];
} repr_cache[1 << 16];

int py_float_repr(double value, char *out)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof bits);
    uint64_t mixed = bits * 0x9E3779B97F4A7C15ull;
    size_t slot = (size_t)(mixed >> 48);
    if (repr_cache[slot].length && repr_cache[slot].bits == bits) {
        memcpy(out, repr_cache[slot].text, repr_cache[slot].length + 1u);
        return repr_cache[slot].length;
    }
    int n = float_repr(value, out);
    if (n > 0 && n < (int)sizeof repr_cache[slot].text) {
        repr_cache[slot].bits = bits;
        repr_cache[slot].length = (unsigned char)n;
        memcpy(repr_cache[slot].text, out, (size_t)n + 1);
    }
    return n;
}

static int roundtrips(int negative, const char *digits, int nd, int exponent, double value)
{
    char text[48];
    int n = 0;
    if (negative) text[n++] = '-';
    text[n++] = digits[0];
    text[n++] = '.';
    for (int i = 1; i < nd; i++) text[n++] = digits[i];
    snprintf(text + n, sizeof text - (size_t)n, "e%d", exponent);
    return strtod(text, 0) == value;
}

static int float_repr(double value, char *out)
{
    if (isnan(value)) return snprintf(out, 32, "nan");
    if (isinf(value)) return snprintf(out, 32, value < 0 ? "-inf" : "inf");
    if (value == 0.0) return snprintf(out, 32, signbit(value) ? "-0.0" : "0.0");
    /* For a normal double a shortest repr of at most 15 digits is the 15-digit
     * rounding with its trailing zeros dropped; subnormals need every width.
     * Near a power of two the rounding interval is lopsided, so a width whose
     * correct rounding misses can still hold its neighbor (repr's choice). */
    char digits[24];
    int nd = 0, exponent = 0, negative = signbit(value) != 0, found = 0;
    for (int precision = fabs(value) < DBL_MIN ? 1 : 15; precision <= 17 && !found; precision++) {
        char sci[40];
        snprintf(sci, sizeof sci, "%.*e", precision - 1, value);
        const char *c = sci + negative;
        nd = 0;
        for (; *c && *c != 'e'; c++)
            if (*c != '.') digits[nd++] = *c;
        exponent = atoi(c + 1);
        if (roundtrips(negative, digits, nd, exponent, value)) {
            found = 1;
            break;
        }
        for (int step = -1; step <= 1 && !found; step += 2) {
            char near[24];
            int near_exponent = exponent;
            memcpy(near, digits, (size_t)nd);
            if (step > 0) {
                int i = nd - 1;
                while (i >= 0 && near[i] == '9') near[i--] = '0';
                if (i >= 0) {
                    near[i]++;
                } else {
                    near[0] = '1';
                    near_exponent++;
                }
            } else {
                int i = nd - 1;
                while (i >= 0 && near[i] == '0') near[i--] = '9';
                if (i >= 0) near[i]--;
                if (near[0] == '0') {
                    for (int k = 0; k < nd; k++) near[k] = '9';
                    near_exponent--;
                }
            }
            if (roundtrips(negative, near, nd, near_exponent, value)) {
                memcpy(digits, near, (size_t)nd);
                exponent = near_exponent;
                found = 1;
            }
        }
    }
    while (nd > 1 && digits[nd - 1] == '0') nd--;
    int n = 0;
    if (negative) out[n++] = '-';
    int decpt = exponent + 1;
    if (decpt > -4 && decpt <= 16) {
        if (decpt <= 0) {
            out[n++] = '0';
            out[n++] = '.';
            for (int i = 0; i < -decpt; i++) out[n++] = '0';
            for (int i = 0; i < nd; i++) out[n++] = digits[i];
        } else if (decpt >= nd) {
            for (int i = 0; i < nd; i++) out[n++] = digits[i];
            for (int i = nd; i < decpt; i++) out[n++] = '0';
            out[n++] = '.';
            out[n++] = '0';
        } else {
            for (int i = 0; i < nd; i++) {
                if (i == decpt) out[n++] = '.';
                out[n++] = digits[i];
            }
        }
    } else {
        out[n++] = digits[0];
        if (nd > 1) {
            out[n++] = '.';
            for (int i = 1; i < nd; i++) out[n++] = digits[i];
        }
        n += snprintf(out + n, 16, "e%c%02d", exponent < 0 ? '-' : '+', exponent < 0 ? -exponent : exponent);
    }
    out[n] = 0;
    return n;
}

void jb_float(jbuf *b, double value)
{
    if (isnan(value)) {
        jb_lit(b, "NaN");
        return;
    }
    if (isinf(value)) {
        jb_lit(b, value < 0 ? "-Infinity" : "Infinity");
        return;
    }
    char text[40];
    int n = py_float_repr(value, text);
    jb_raw(b, text, (size_t)n);
}

/* json.dumps with ensure_ascii: non-ASCII as \uXXXX (surrogate pairs). */
void jb_string(jbuf *b, const char *text, size_t length)
{
    jb_raw(b, "\"", 1);
    const unsigned char *s = (const unsigned char *)text;
    size_t start = 0;
    for (size_t i = 0; i < length;) {
        unsigned char c = s[i];
        const char *escape = 0;
        char unicode[16];
        size_t width = 1;
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
        } else if (c >= 0x80) {
            unsigned code = 0;
            if ((c & 0xE0) == 0xC0 && i + 1 < length) {
                code = ((c & 0x1Fu) << 6) | (s[i + 1] & 0x3Fu);
                width = 2;
            } else if ((c & 0xF0) == 0xE0 && i + 2 < length) {
                code = ((c & 0x0Fu) << 12) | ((s[i + 1] & 0x3Fu) << 6) | (s[i + 2] & 0x3Fu);
                width = 3;
            } else if (i + 3 < length) {
                code = ((c & 0x07u) << 18) | ((s[i + 1] & 0x3Fu) << 12) | ((s[i + 2] & 0x3Fu) << 6) | (s[i + 3] & 0x3Fu);
                width = 4;
            }
            if (code >= 0x10000) {
                code -= 0x10000;
                snprintf(unicode, sizeof unicode, "\\u%04x\\u%04x", 0xD800 + (code >> 10), 0xDC00 + (code & 0x3FF));
            } else {
                snprintf(unicode, sizeof unicode, "\\u%04x", code);
            }
            escape = unicode;
        }
        if (escape) {
            jb_raw(b, text + start, i - start);
            jb_lit(b, escape);
            i += width;
            start = i;
        } else {
            i++;
        }
    }
    jb_raw(b, text + start, length - start);
    jb_raw(b, "\"", 1);
}

void jb_cstring(jbuf *b, const char *text) { jb_string(b, text, strlen(text)); }

void jb_value(jbuf *b, const jvalue *v)
{
    switch (v->kind) {
    case J_NULL: jb_lit(b, "null"); break;
    case J_FALSE: jb_lit(b, "false"); break;
    case J_TRUE: jb_lit(b, "true"); break;
    case J_INT: jb_int(b, v->u.i); break;
    case J_FLOAT: jb_float(b, v->u.f); break;
    case J_STRING: jb_string(b, v->u.s, (size_t)v->count); break;
    case J_ARRAY:
        jb_raw(b, "[", 1);
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_raw(b, ", ", 2);
            jb_value(b, v->u.items + i);
        }
        jb_raw(b, "]", 1);
        break;
    case J_OBJECT:
        jb_raw(b, "{", 1);
        for (int64_t i = 0; i < v->count; i++) {
            if (i) jb_raw(b, ", ", 2);
            jb_value(b, v->u.members + 2 * i);
            jb_raw(b, ": ", 2);
            jb_value(b, v->u.members + 2 * i + 1);
        }
        jb_raw(b, "}", 1);
        break;
    }
}
