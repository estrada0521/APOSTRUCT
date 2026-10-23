/* Bake the Source tables into one blob (source_blob.h).
 *
 *   source_blob_build <Source dir> <output file>
 *
 * layout.txt names each blob section and how it is rendered from a table column; the rendered items are then typed
 * as a block of text: quoted strings only make a string section; any number outside quotes makes it numeric (real
 * when a number has '.', 'e' or 'E'); bare T/F tokens make it boolean; otherwise its quoted strings (possibly none). */

#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "source_blob.h"

static void die(const char *message, const char *detail)
{
    fprintf(stderr, "source_blob_build: %s%s%s\n", message, detail ? ": " : "", detail ? detail : "");
    exit(1);
}

static void *xmalloc(size_t size)
{
    void *p = malloc(size ? size : 1);
    if (!p) die("out of memory", 0);
    return p;
}

static void *xrealloc(void *p, size_t size)
{
    p = realloc(p, size ? size : 1);
    if (!p) die("out of memory", 0);
    return p;
}

static char *read_file(const char *path, size_t *size)
{
    FILE *f = fopen(path, "rb");
    if (!f) die("cannot open", path);
    fseek(f, 0, SEEK_END);
    long length = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *data = xmalloc((size_t)length + 1);
    if (fread(data, 1, (size_t)length, f) != (size_t)length) die("cannot read", path);
    fclose(f);
    data[length] = 0;
    *size = (size_t)length;
    return data;
}

/* ---- SHA-256 ---- */

static const uint32_t K256[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
};

#define ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))

static void sha256_block(uint32_t *h, const unsigned char *p)
{
    uint32_t w[64];
    for (int i = 0; i < 16; i++)
        w[i] = (uint32_t)p[4 * i] << 24 | (uint32_t)p[4 * i + 1] << 16 | (uint32_t)p[4 * i + 2] << 8 | p[4 * i + 3];
    for (int i = 16; i < 64; i++) {
        uint32_t s0 = ROTR(w[i - 15], 7) ^ ROTR(w[i - 15], 18) ^ (w[i - 15] >> 3);
        uint32_t s1 = ROTR(w[i - 2], 17) ^ ROTR(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (int i = 0; i < 64; i++) {
        uint32_t s1 = ROTR(e, 6) ^ ROTR(e, 11) ^ ROTR(e, 25);
        uint32_t t1 = hh + s1 + ((e & f) ^ (~e & g)) + K256[i] + w[i];
        uint32_t s0 = ROTR(a, 2) ^ ROTR(a, 13) ^ ROTR(a, 22);
        uint32_t t2 = s0 + ((a & b) ^ (a & c) ^ (b & c));
        hh = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
}

static void sha256_hex(const unsigned char *data, size_t size, char *hex)
{
    uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    size_t full = size / 64 * 64;
    for (size_t i = 0; i < full; i += 64) sha256_block(h, data + i);
    unsigned char tail[128] = {0};
    size_t rest = size - full;
    memcpy(tail, data + full, rest);
    tail[rest] = 0x80;
    size_t tail_size = rest + 9 <= 64 ? 64 : 128;
    uint64_t bits = (uint64_t)size * 8;
    for (int i = 0; i < 8; i++) tail[tail_size - 1 - i] = (unsigned char)(bits >> (8 * i));
    for (size_t i = 0; i < tail_size; i += 64) sha256_block(h, tail + i);
    for (int i = 0; i < 8; i++) sprintf(hex + 8 * i, "%08x", h[i]);
}

/* ---- sections ---- */

typedef struct {
    char name[BLOB_NAME_SIZE];
    uint32_t kind;
    uint64_t count;
    unsigned char *payload;
    uint64_t bytes;
} section;

static section *sections;
static int section_count, section_capacity;

static section *add_section(const char *name)
{
    for (int i = 0; i < section_count; i++)
        if (!strcmp(sections[i].name, name)) {
            free(sections[i].payload);
            return &sections[i];
        }
    if (section_count == section_capacity) {
        section_capacity = section_capacity ? 2 * section_capacity : 256;
        sections = xrealloc(sections, sizeof *sections * (size_t)section_capacity);
    }
    section *s = &sections[section_count++];
    memset(s, 0, sizeof *s);
    if (strlen(name) >= BLOB_NAME_SIZE) die("section name too long", name);
    strcpy(s->name, name);
    return s;
}

/* Python re [-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)? scanned left to right. */
static size_t match_number(const char *s, size_t i, size_t n, int *real)
{
    size_t j = i;
    if (j < n && (s[j] == '+' || s[j] == '-')) j++;
    size_t digits = j;
    while (j < n && isdigit((unsigned char)s[j])) j++;
    if (j == digits) return 0;
    *real = 0;
    if (j < n && s[j] == '.') {
        j++;
        while (j < n && isdigit((unsigned char)s[j])) j++;
        *real = 1;
    }
    if (j < n && (s[j] == 'e' || s[j] == 'E')) {
        size_t k = j + 1;
        if (k < n && (s[k] == '+' || s[k] == '-')) k++;
        size_t exp_digits = k;
        while (k < n && isdigit((unsigned char)s[k])) k++;
        if (k > exp_digits) {
            j = k;
            *real = 1;
        }
    }
    return j - i;
}

static void parse_block(const char *table, const char *name, const char *block, size_t n)
{
    /* quoted strings and the remaining text */
    size_t string_count = 0, string_bytes = 0;
    char *rest = xmalloc(n + 1);
    for (size_t i = 0; i < n; i++) {
        if (block[i] == '"') {
            size_t j = i + 1;
            while (j < n && block[j] != '"') j++;
            if (j < n) {
                string_count++;
                string_bytes += j - i - 1;
                rest[i] = ' ';
                for (size_t k = i + 1; k <= j; k++) rest[k] = ' ';
                i = j;
                continue;
            }
        }
        rest[i] = block[i];
    }
    rest[n] = 0;
    /* numbers outside quotes */
    size_t number_count = 0;
    int any_real = 0;
    for (size_t i = 0; i < n;) {
        int real = 0;
        size_t m = match_number(rest, i, n, &real);
        if (m) {
            number_count++;
            any_real |= real;
            i += m;
        } else {
            i++;
        }
    }
    char full[BLOB_NAME_SIZE];
    if ((size_t)snprintf(full, sizeof full, "%s.%s", table, name) >= sizeof full) die("section name too long", name);
    section *s = add_section(full);
    if (number_count) {
        if (any_real) {
            double *values = xmalloc(sizeof(double) * number_count);
            size_t k = 0;
            for (size_t i = 0; i < n;) {
                int real = 0;
                size_t m = match_number(rest, i, n, &real);
                if (m) {
                    char token[64];
                    if (m >= sizeof token) die("number too long", full);
                    memcpy(token, rest + i, m);
                    token[m] = 0;
                    values[k++] = strtod(token, 0);
                    i += m;
                } else {
                    i++;
                }
            }
            s->kind = BLOB_F64;
            s->count = number_count;
            s->payload = (unsigned char *)values;
            s->bytes = sizeof(double) * number_count;
        } else {
            int64_t *values = xmalloc(sizeof(int64_t) * number_count);
            size_t k = 0;
            int64_t low = 0, high = 0;
            for (size_t i = 0; i < n;) {
                int real = 0;
                size_t m = match_number(rest, i, n, &real);
                if (m) {
                    char token[64];
                    if (m >= sizeof token) die("number too long", full);
                    memcpy(token, rest + i, m);
                    token[m] = 0;
                    errno = 0;
                    long long v = strtoll(token, 0, 10);
                    if (errno) die("integer out of range", full);
                    values[k] = v;
                    if (!k || v < low) low = v;
                    if (!k || v > high) high = v;
                    k++;
                    i += m;
                } else {
                    i++;
                }
            }
            int width = low >= INT8_MIN && high <= INT8_MAX ? 1
                      : low >= INT16_MIN && high <= INT16_MAX ? 2
                      : low >= INT32_MIN && high <= INT32_MAX ? 4 : 8;
            unsigned char *payload = xmalloc((size_t)width * number_count);
            for (size_t i = 0; i < number_count; i++) {
                if (width == 1) ((int8_t *)payload)[i] = (int8_t)values[i];
                else if (width == 2) ((int16_t *)payload)[i] = (int16_t)values[i];
                else if (width == 4) ((int32_t *)payload)[i] = (int32_t)values[i];
                else ((int64_t *)payload)[i] = values[i];
            }
            free(values);
            s->kind = width == 1 ? BLOB_I8 : width == 2 ? BLOB_I16 : width == 4 ? BLOB_I32 : BLOB_I64;
            s->count = number_count;
            s->payload = payload;
            s->bytes = (uint64_t)width * number_count;
        }
        free(rest);
        return;
    }
    if (!string_count) {
        /* bare T/F tokens */
        size_t tokens = 0, bools = 0;
        for (size_t i = 0; i < n;) {
            while (i < n && isspace((unsigned char)rest[i])) i++;
            if (i >= n) break;
            size_t j = i;
            while (j < n && !isspace((unsigned char)rest[j])) j++;
            tokens++;
            if (j - i == 1 && (rest[i] == 'T' || rest[i] == 'F')) bools++;
            i = j;
        }
        if (tokens && tokens == bools) {
            unsigned char *payload = xmalloc(tokens);
            size_t k = 0;
            for (size_t i = 0; i < n; i++)
                if (rest[i] == 'T' || rest[i] == 'F') payload[k++] = rest[i] == 'T';
            s->kind = BLOB_BOOL;
            s->count = tokens;
            s->payload = payload;
            s->bytes = tokens;
            free(rest);
            return;
        }
    }
    /* string section: offsets[count + 1] then bytes */
    uint64_t bytes = sizeof(int64_t) * (string_count + 1) + string_bytes;
    unsigned char *payload = xmalloc((size_t)bytes);
    int64_t *offsets = (int64_t *)payload;
    unsigned char *pool = payload + sizeof(int64_t) * (string_count + 1);
    size_t k = 0, fill = 0;
    for (size_t i = 0; i < n; i++) {
        if (block[i] != '"') continue;
        size_t j = i + 1;
        while (j < n && block[j] != '"') j++;
        if (j >= n) break;
        offsets[k++] = (int64_t)fill;
        memcpy(pool + fill, block + i + 1, j - i - 1);
        fill += j - i - 1;
        i = j;
    }
    offsets[k] = (int64_t)fill;
    s->kind = BLOB_STR;
    s->count = string_count;
    s->payload = payload;
    s->bytes = bytes;
    free(rest);
}

/* Exact constants: a + b sqrt2 + c sqrt3 + d sqrt6 with rational parts, parsed
 * from the form column of constants.txt. Forms outside that field (nested
 * radicals, 22.5 degrees) parse as absent. */
typedef struct {
    int64_t n, d;
} q64;

typedef struct {
    q64 x[4];
    int ok;
} q23;

static int64_t g64(int64_t a, int64_t b)
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

static q64 qn(__int128 n, __int128 d, int *ok)
{
    if (d < 0) {
        n = -n;
        d = -d;
    }
    __int128 a = n < 0 ? -n : n, b = d;
    while (b) {
        __int128 t = a % b;
        a = b;
        b = t;
    }
    if (a > 1) {
        n /= a;
        d /= a;
    }
    if (!d || n > INT64_MAX || n < -INT64_MAX || d > INT64_MAX) *ok = 0;
    return (q64){(int64_t)n, (int64_t)d};
}

static q64 qadd(q64 a, q64 b, int *ok) { return qn((__int128)a.n * b.d + (__int128)b.n * a.d, (__int128)a.d * b.d, ok); }
static q64 qmul(q64 a, q64 b, int *ok) { return qn((__int128)a.n * b.n, (__int128)a.d * b.d, ok); }

static q23 q23_rat(q64 r)
{
    q23 v = {{{0, 1}, {0, 1}, {0, 1}, {0, 1}}, 1};
    v.x[0] = r;
    return v;
}

static q23 q23_add(q23 a, q23 b, int sign)
{
    q23 v = a;
    v.ok = a.ok && b.ok;
    for (int i = 0; i < 4; i++) v.x[i] = qadd(a.x[i], (q64){sign * b.x[i].n, b.x[i].d}, &v.ok);
    return v;
}

/* basis e0=1, e1=sqrt2, e2=sqrt3, e3=sqrt6: e_i e_j = m * e_k */
static q23 q23_mul(q23 a, q23 b)
{
    static const int K[4][4] = {{0, 1, 2, 3}, {1, 0, 3, 2}, {2, 3, 0, 1}, {3, 2, 1, 0}};
    static const int F[4][4] = {{1, 1, 1, 1}, {1, 2, 1, 2}, {1, 1, 3, 3}, {1, 2, 3, 6}};
    q23 v = q23_rat((q64){0, 1});
    v.ok = a.ok && b.ok;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            q64 t = qmul(a.x[i], b.x[j], &v.ok);
            v.x[K[i][j]] = qadd(v.x[K[i][j]], qmul(t, (q64){F[i][j], 1}, &v.ok), &v.ok);
        }
    return v;
}

static q23 q23_conj(q23 a, int which)
{
    /* which 1: sqrt2 -> -sqrt2, which 2: sqrt3 -> -sqrt3 */
    q23 v = a;
    if (which == 1) {
        v.x[1].n = -v.x[1].n;
        v.x[3].n = -v.x[3].n;
    } else {
        v.x[2].n = -v.x[2].n;
        v.x[3].n = -v.x[3].n;
    }
    return v;
}

static q23 q23_div(q23 a, q23 b)
{
    q23 num = q23_mul(a, q23_conj(b, 1)), den = q23_mul(b, q23_conj(b, 1));
    num = q23_mul(num, q23_conj(den, 2));
    den = q23_mul(den, q23_conj(den, 2));
    q23 v = num;
    v.ok = num.ok && den.ok && den.x[0].n && !den.x[1].n && !den.x[2].n && !den.x[3].n;
    if (!v.ok) return v;
    q64 inv = {den.x[0].d, den.x[0].n};
    for (int i = 0; i < 4; i++) v.x[i] = qmul(v.x[i], inv, &v.ok);
    return v;
}

static q23 q23_sqrt(q23 a)
{
    q23 bad = q23_rat((q64){0, 1});
    bad.ok = 0;
    if (!a.ok || a.x[1].n || a.x[2].n || a.x[3].n || a.x[0].n < 0) return bad;
    /* sqrt(n/d) = sqrt(n d) / d with n d = s^2 * f, f squarefree in {1,2,3,6} */
    __int128 m = (__int128)a.x[0].n * a.x[0].d;
    if (!m) return a;
    int64_t f = 1;
    for (int64_t p = 2; p <= 3; p++)
        while (m % ((__int128)p * p) == 0) m /= (__int128)p * p;
    for (int64_t p = 2; p <= 3; p++)
        if (m % p == 0) {
            f *= p;
            m /= p;
        }
    int64_t root = (int64_t)sqrt((double)m);
    while ((__int128)root * root < m) root++;
    while ((__int128)root * root > m) root--;
    if ((__int128)root * root != m) return bad;
    __int128 nd = (__int128)a.x[0].n * a.x[0].d, s2 = nd / f;
    int64_t s = (int64_t)sqrt((double)s2);
    while ((__int128)s * s < s2) s++;
    while ((__int128)s * s > s2) s--;
    if ((__int128)s * s != s2) return bad;
    q23 v = q23_rat((q64){0, 1});
    int slot = f == 1 ? 0 : f == 2 ? 1 : f == 3 ? 2 : 3;
    v.x[slot] = qn(s, a.x[0].d, &v.ok);
    return v;
}

/* cos and sin of a multiple of 15 degrees */
static q23 q23_trig(q64 degrees, int sine)
{
    q23 bad = q23_rat((q64){0, 1});
    bad.ok = 0;
    if (degrees.d != 1 || degrees.n % 15) return bad;
    int64_t step = ((degrees.n / 15) % 24 + 24) % 24;
    if (sine) step = (step + 18) % 24;
    /* cos(15 j) for j = 0..6 over e0..e3 in quarters */
    static const int C[7][4] = {{4, 0, 0, 0}, {0, 1, 0, 1}, {0, 0, 2, 0}, {0, 2, 0, 0}, {2, 0, 0, 0}, {0, -1, 0, 1}, {0, 0, 0, 0}};
    int sign = 1, j = (int)step;
    if (j > 12) j = 24 - j;
    if (j > 6) {
        j = 12 - j;
        sign = -1;
    }
    q23 v = q23_rat((q64){0, 1});
    for (int i = 0; i < 4; i++) v.x[i] = qn(sign * C[j][i], 4, &v.ok);
    return v;
}

typedef struct {
    const char *p;
} q23_parser;

static q23 q23_expr(q23_parser *s);

static void q23_skip(q23_parser *s)
{
    while (*s->p == ' ' || *s->p == '\t') s->p++;
}

static q23 q23_factor(q23_parser *s)
{
    q23 bad = q23_rat((q64){0, 1});
    bad.ok = 0;
    q23_skip(s);
    if (*s->p == '-') {
        s->p++;
        q23 v = q23_factor(s);
        return q23_add(q23_rat((q64){0, 1}), v, -1);
    }
    if (*s->p == '(') {
        s->p++;
        q23 v = q23_expr(s);
        q23_skip(s);
        if (*s->p != ')') return bad;
        s->p++;
        return v;
    }
    if (isdigit((unsigned char)*s->p)) {
        int64_t n = 0, d = 1;
        while (isdigit((unsigned char)*s->p)) n = n * 10 + (*s->p++ - '0');
        if (*s->p == '.') {
            s->p++;
            while (isdigit((unsigned char)*s->p)) {
                n = n * 10 + (*s->p++ - '0');
                d *= 10;
            }
        }
        int ok = 1;
        q23 v = q23_rat(qn(n, d, &ok));
        v.ok = ok;
        return v;
    }
    static const char *const names[] = {"sqrt", "cos", "sin", "tan"};
    for (int f = 0; f < 4; f++) {
        size_t len = strlen(names[f]);
        if (strncmp(s->p, names[f], len) || s->p[len] != '(') continue;
        s->p += len + 1;
        q23 arg = q23_expr(s);
        q23_skip(s);
        if (*s->p != ')') return bad;
        s->p++;
        if (f == 0) return q23_sqrt(arg);
        if (!arg.ok || arg.x[1].n || arg.x[2].n || arg.x[3].n) return bad;
        if (f == 3) return q23_div(q23_trig(arg.x[0], 1), q23_trig(arg.x[0], 0));
        return q23_trig(arg.x[0], f == 2);
    }
    return bad;
}

static q23 q23_term(q23_parser *s)
{
    q23 v = q23_factor(s);
    for (;;) {
        q23_skip(s);
        char op = *s->p;
        if (op != '*' && op != '/') return v;
        s->p++;
        q23 r = q23_factor(s);
        v = op == '*' ? q23_mul(v, r) : q23_div(v, r);
    }
}

static q23 q23_expr(q23_parser *s)
{
    q23 v = q23_term(s);
    for (;;) {
        q23_skip(s);
        char op = *s->p;
        if (op != '+' && op != '-') return v;
        s->p++;
        q23 r = q23_term(s);
        v = q23_add(v, r, op == '+' ? 1 : -1);
    }
}

static q23 q23_parse(const char *text)
{
    q23_parser s = {text};
    q23 v = q23_expr(&s);
    q23_skip(&s);
    if (*s.p) v.ok = 0;
    return v;
}

static double q23_value(q23 v)
{
    const double e[4] = {1.0, sqrt(2.0), sqrt(3.0), sqrt(6.0)};
    double sum = 0.0;
    for (int i = 0; i < 4; i++) sum += (double)v.x[i].n / (double)v.x[i].d * e[i];
    return sum;
}

/* constants.txt: a header line, then one constant per line: number (1, 2, ... in order), decimal value, exact form
 * (empty when there is none).  A reference k in the tables is constant k, -k its negative, 0 no reference. */
static void load_constants(const char *dir)
{
    char path[4096];
    snprintf(path, sizeof path, "%s/constants.txt", dir);
    size_t size;
    char *text = read_file(path, &size);
    int64_t max_code = 0;
    int64_t codes[4096];
    double values[4096];
    static q23 exact[4096];
    int count = 0;
    char *line = strchr(text, '\n');
    line = line ? line + 1 : text + size;
    while (*line) {
        char *end = strchr(line, '\n');
        if (end) *end = 0;
        if (*line) {
            if (count >= 4096) die("too many constants", 0);
            char *after, *tab;
            codes[count] = strtoll(line, &after, 10);
            if (codes[count] != count + 1 || *after != '\t') die("constants must be numbered 1, 2, ... in order", line);
            values[count] = strtod(after + 1, &tab);
            if (*tab != '\t') die("constant needs a value and a form column", line);
            exact[count] = q23_parse(tab + 1);
            if (exact[count].ok && fabs(q23_value(exact[count]) - values[count]) > 1e-13 * (1.0 + fabs(values[count])))
                die("constant form does not match its value", tab + 1);
            if (codes[count] > max_code) max_code = codes[count];
            count++;
        }
        if (!end) break;
        line = end + 1;
    }
    free(text);
    int64_t n = max_code + 1;
    double *dense = xmalloc(sizeof(double) * (size_t)n);
    unsigned char *present = xmalloc((size_t)n);
    memset(dense, 0, sizeof(double) * (size_t)n);
    memset(present, 0, (size_t)n);
    for (int i = 0; i < count; i++) {
        dense[codes[i]] = values[i];
        present[codes[i]] = 1;
    }
    int64_t *exact_words = xmalloc(sizeof(int64_t) * 8 * (size_t)n);
    memset(exact_words, 0, sizeof(int64_t) * 8 * (size_t)n);
    for (int i = 0; i < count; i++) {
        if (!exact[i].ok) continue;
        for (int j = 0; j < 4; j++) {
            exact_words[8 * codes[i] + 2 * j] = exact[i].x[j].n;
            exact_words[8 * codes[i] + 2 * j + 1] = exact[i].x[j].d;
        }
    }
    section *s = add_section("const.exact");
    s->kind = BLOB_I64;
    s->count = (uint64_t)n * 8;
    s->payload = (unsigned char *)exact_words;
    s->bytes = sizeof(int64_t) * 8 * (uint64_t)n;
    s = add_section("const.values");
    s->kind = BLOB_F64;
    s->count = (uint64_t)n;
    s->payload = (unsigned char *)dense;
    s->bytes = sizeof(double) * (uint64_t)n;
    s = add_section("const.present");
    s->kind = BLOB_BOOL;
    s->count = (uint64_t)n;
    s->payload = present;
    s->bytes = (uint64_t)n;
}


/* ---- entity tables and layout.txt ---- */

/* A table file: a header of tab-separated column names, then one record per line with tab-separated fields; a
 * numeric field is a space-separated list, a text field is the text itself. */
typedef struct {
    char path[512];
    int columns, rows;
    char **names;
    char ***fields; /* rows x columns, each a pointer into text */
    char *text;
} entity_table;

static entity_table *etables;
static int etable_count;

static entity_table *table_of(const char *dir, const char *file)
{
    for (int i = 0; i < etable_count; i++)
        if (!strcmp(etables[i].path, file)) return &etables[i];
    etables = xrealloc(etables, sizeof *etables * (size_t)(etable_count + 1));
    entity_table *t = &etables[etable_count++];
    memset(t, 0, sizeof *t);
    snprintf(t->path, sizeof t->path, "%s", file);
    char path[4096];
    snprintf(path, sizeof path, "%s/%s", dir, file);
    size_t size;
    t->text = read_file(path, &size);
    char *line = t->text, *end = t->text + size;
    int header = 1;
    while (line < end) {
        char *next = memchr(line, '\n', (size_t)(end - line));
        if (!next) next = end;
        *next = 0;
        if (*line != '#') {
            int n = 1;
            for (char *c = line; *c; c++) n += *c == '\t';
            char **cells = xmalloc(sizeof *cells * (size_t)n);
            int k = 0;
            cells[k++] = line;
            for (char *c = line; *c; c++)
                if (*c == '\t') *c = 0, cells[k++] = c + 1;
            if (header) {
                t->names = cells;
                t->columns = n;
                header = 0;
            } else {
                if (n != t->columns) die("wrong field count in", file);
                t->fields = xrealloc(t->fields, sizeof *t->fields * (size_t)(t->rows + 1));
                t->fields[t->rows++] = cells;
            }
        }
        line = next + 1;
    }
    return t;
}

static int column_of(const entity_table *t, const char *name)
{
    for (int i = 0; i < t->columns; i++)
        if (!strcmp(t->names[i], name)) return i;
    die("no such column", name);
    return -1;
}

typedef struct {
    char *data;
    size_t used, capacity;
} text_buf;

static void put(text_buf *b, const char *s, size_t n)
{
    if (b->used + n + 2 > b->capacity) {
        b->capacity = 2 * (b->used + n + 2);
        b->data = xrealloc(b->data, b->capacity);
    }
    memcpy(b->data + b->used, s, n);
    b->used += n;
}

static void put_token(text_buf *b, const char *s, size_t n)
{
    put(b, s, n);
    put(b, (b->used % 1024) > 960 ? "\n" : " ", 1);
}

static void put_int(text_buf *b, long long v)
{
    char x[32];
    int n = snprintf(x, sizeof x, "%lld", v);
    put_token(b, x, (size_t)n);
}

/* the space-separated tokens of a numeric field */
static int field_tokens(const char *field, const char **starts, size_t *lengths, int capacity)
{
    int n = 0;
    for (const char *c = field; *c;) {
        while (*c == ' ') c++;
        if (!*c) break;
        const char *a = c;
        while (*c && *c != ' ') c++;
        if (n == capacity) die("field too long", field);
        starts[n] = a, lengths[n] = (size_t)(c - a), n++;
    }
    return n;
}

static long long option(const char *options, const char *key, long long fallback)
{
    size_t k = strlen(key);
    for (const char *c = options; c && *c; c++)
        if ((c == options || c[-1] == ' ') && !strncmp(c, key, k) && (c[k] == '=' || c[k] == ' ' || !c[k]))
            return c[k] == '=' ? atoll(c + k + 1) : 1;
    return fallback;
}

enum { MAX_TOKENS = 1 << 16 };

/* the magnetic point operations: each family of point_ops.txt (a block starting at an identity matrix), then the same
 * family with time reversal */
static long long magnetic_number(const char *dir, long long k, int primed)
{
    entity_table *ops = table_of(dir, "group/point_ops.txt");
    int c = column_of(ops, "matrix");
    long long start = 0, offset = 0;
    for (int r = 0; r <= ops->rows; r++) {
        if (r == ops->rows || (r > 0 && !strcmp(ops->fields[r][c], "1 0 0 0 1 0 0 0 1"))) {
            long long size = r - start;
            if (k > start && k <= r) return offset + (primed ? size : 0) + (k - start);
            offset += 2 * size, start = r;
        }
    }
    die("no such point operation", 0);
    return 0;
}

/* one blob section from its layout line: SECTION FILE COLUMN MODE [OPTIONS] */
static char byu_names[1 << 16];
static size_t byu_used;
static int byu_count;

static void render_section(const char *dir, char *line)
{
    char section[BLOB_NAME_SIZE], file[512], column[128], mode[32];
    int consumed = 0;
    if (sscanf(line, "%95s %511s %127s %31s %n", section, file, column, mode, &consumed) < 4) die("bad layout line", line);
    const char *options = line + consumed;
    const char *byu = strstr(options, "byu=");
    if (byu) {
        size_t n = strcspn(byu + 4, " \t");
        int w = snprintf(byu_names + byu_used, sizeof byu_names - byu_used, "%s=%.*s", section, (int)n, byu + 4);
        if (w < 0 || (size_t)w >= sizeof byu_names - byu_used) die("too many BYU names", section);
        byu_used += (size_t)w + 1;
        byu_count++;
    }
    if (!strcmp(mode, "derived")) return;
    char *dot = strchr(section, '.');
    if (!dot) die("section without table", section);
    *dot = 0;
    entity_table *t = table_of(dir, file);
    int c = column_of(t, column);
    long long first = option(options, "first", t->rows), unit = option(options, "div", 1), lead = option(options, "lead", 0);
    long long pad = option(options, "pad", 0), total = option(options, "total", 0), keys = option(options, "keys", 0);
    long long op_every = option(options, "ops", 0);
    const char *per_name = strstr(options, "per=");
    int per = -1;
    if (per_name) {
        char name[128];
        sscanf(per_name + 4, "%127s", name);
        per = column_of(t, name);
    }
    int zero_empty = option(options, "zero_empty", 0) > 0, sentinel = option(options, "sentinel", 0) > 0;
    int rows = (int)(first < t->rows ? first : t->rows);
    long long grid = option(options, "grid", 0);
    int *place = 0;
    if (grid) {
        int slot = column_of(t, "slot");
        place = xmalloc(sizeof *place * (size_t)(keys * grid));
        for (long long i = 0; i < keys * grid; i++) place[i] = -1;
        for (int r = 0; r < rows; r++) {
            long long key = atoll(t->fields[r][0]), at = atoll(t->fields[r][slot]);
            if (key < 1 || key > keys || at < 1 || at > grid || place[(key - 1) * grid + at - 1] >= 0)
                die("bad grid position", t->fields[r][0]);
            place[(key - 1) * grid + at - 1] = r;
        }
        rows = (int)(keys * grid);
    }
#define FIELD(r) (place ? (place[r] < 0 ? "" : t->fields[place[r]][c]) : t->fields[r][c])
    static const char *starts[MAX_TOKENS];
    static size_t lengths[MAX_TOKENS];
    text_buf b = {0};
    put(&b, "", 0);
    if (!strcmp(mode, "text")) {
        for (int r = 0; r < rows; r++) {
            put(&b, "\"", 1);
            put(&b, FIELD(r), strlen(FIELD(r)));
            put_token(&b, "\"", 1);
        }
    } else if (!strcmp(mode, "each")) {
        for (int r = 0; r < rows; r++) {
            int n = field_tokens(FIELD(r), starts, lengths, MAX_TOKENS);
            if (n) put_token(&b, starts[0], lengths[0]);
            else put_int(&b, 0);
        }
    } else if (!strcmp(mode, "flat")) {
        long long written = 0;
        for (long long i = 0; i < lead; i++, written++) put_int(&b, 0);
        for (int r = 0; r < rows; r++) {
            int n = field_tokens(FIELD(r), starts, lengths, MAX_TOKENS);
            for (int i = 0; i < n; i++, written++) {
                if (op_every && i % op_every == op_every - 1) {
                    char *after;
                    long long k = strtoll(starts[i], &after, 10);
                    int primed = after < starts[i] + lengths[i] && *after == '\'';
                    if (after + primed != starts[i] + lengths[i]) die("bad operation", FIELD(r));
                    put_int(&b, magnetic_number(dir, k, primed));
                } else {
                    put_token(&b, starts[i], lengths[i]);
                }
            }
            for (long long i = n; i < pad; i++, written++) put_int(&b, 0);
        }
        for (; written < total; written++) put_int(&b, 0);
    } else if (!strcmp(mode, "start") || !strcmp(mode, "count")) {
        long long at = 1 + lead;
        for (int r = 0; r < rows; r++) {
            int n = field_tokens(FIELD(r), starts, lengths, MAX_TOKENS);
            if (mode[0] == 'c' && per >= 0) {
                long long by = place && place[r] < 0 ? 0 : atoll(place ? t->fields[place[r]][per] : t->fields[r][per]);
                put_int(&b, by ? n / by : 0);
            } else if (mode[0] == 'c') put_int(&b, n / unit);
            else put_int(&b, n == 0 && zero_empty ? 0 : (at - 1) / unit + 1);
            at += n;
        }
        if (sentinel && mode[0] == 's') put_int(&b, (at - 1) / unit + 1);
    } else if (!strcmp(mode, "group_start") || !strcmp(mode, "group_count")) {
        long long previous = 1;
        for (long long k = 1; k <= keys; k++) {
            long long start = 0, count = 0;
            for (int r = 0; r < rows; r++)
                if (atoll(t->fields[r][c]) == k) {
                    if (!count) start = r + 1;
                    count++;
                }
            if (!count) start = previous;
            previous = start;
            put_int(&b, mode[6] == 'c' ? count : start);
        }
        if (sentinel) put_int(&b, rows + 1);
    } else if (!strcmp(mode, "directions")) {
        for (int r = 0; r < rows; r++) {
            int n = field_tokens(FIELD(r), starts, lengths, MAX_TOKENS), d = 0;
            for (int v = 4; v + 3 <= n; v += 4)
                if (strncmp(starts[v], "0", lengths[v]) || strncmp(starts[v + 1], "0", lengths[v + 1]) ||
                    strncmp(starts[v + 2], "0", lengths[v + 2]))
                    d++;
            put_int(&b, d);
        }
    } else if (!strncmp(mode, "magnetic_", 9)) {
        for (int primed = 0; primed < 2 && mode[9] == 'f'; primed++)
            for (int r = 1; r <= t->rows; r++) put_int(&b, magnetic_number(dir, r, primed));
        for (long long m = 1; mode[9] != 'f'; m++) {
            int found = 0;
            for (int r = 1; r <= t->rows && !found; r++)
                for (int primed = 0; primed < 2 && !found; primed++)
                    if (magnetic_number(dir, r, primed) == m) {
                        found = 1;
                        if (mode[9] == 'o') put_int(&b, r);
                        else put_token(&b, primed ? "T" : "F", 1);
                    }
            if (!found) break;
        }
    } else if (!strcmp(mode, "records")) {
        put_int(&b, rows);
    } else {
        die("unknown layout mode", mode);
    }
#undef FIELD
    parse_block(section, dot + 1, b.data, b.used);
    free(b.data);
    free(place);
}

/* every section of layout.txt; returns the table files it reads */
static int load_layout(const char *dir)
{
    char path[4096];
    snprintf(path, sizeof path, "%s/layout.txt", dir);
    size_t size;
    char *text = read_file(path, &size);
    for (char *line = text; line < text + size;) {
        char *next = memchr(line, '\n', (size_t)(text + size - line));
        if (!next) next = text + size;
        *next = 0;
        char *a = line;
        while (*a == ' ') a++;
        if (*a && *a != '#') render_section(dir, a);
        line = next + 1;
    }
    free(text);
    return etable_count;
}

static void add_hashes(const char *dir, const char *const *files, int file_count)
{
    size_t total = 0;
    char lines[64][256];
    for (int i = 0; i < file_count; i++) {
        char path[4096];
        snprintf(path, sizeof path, "%s/%s", dir, files[i]);
        size_t size;
        char *data = read_file(path, &size);
        char hex[65];
        sha256_hex((unsigned char *)data, size, hex);
        free(data);
        snprintf(lines[i], sizeof lines[i], "%s=%s", files[i], hex);
        total += strlen(lines[i]);
    }
    section *s = add_section("source.sha256");
    uint64_t bytes = sizeof(int64_t) * (uint64_t)(file_count + 1) + total;
    unsigned char *payload = xmalloc((size_t)bytes);
    int64_t *offsets = (int64_t *)payload;
    unsigned char *pool = payload + sizeof(int64_t) * (size_t)(file_count + 1);
    size_t fill = 0;
    for (int i = 0; i < file_count; i++) {
        offsets[i] = (int64_t)fill;
        memcpy(pool + fill, lines[i], strlen(lines[i]));
        fill += strlen(lines[i]);
    }
    offsets[file_count] = (int64_t)fill;
    s->kind = BLOB_STR;
    s->count = (uint64_t)file_count;
    s->payload = payload;
    s->bytes = bytes;
}

static uint64_t align64(uint64_t value) { return (value + 63) & ~(uint64_t)63; }

int main(int argc, char **argv)
{
    if (argc != 3) die("usage: source_blob_build <Source dir> <output file>", 0);
    int tables = load_layout(argv[1]);
    text_buf names = {0};
    put(&names, "", 0);
    for (size_t at = 0; at < byu_used; at += strlen(byu_names + at) + 1) {
        put(&names, "\"", 1);
        put(&names, byu_names + at, strlen(byu_names + at));
        put_token(&names, "\"", 1);
    }
    parse_block("source", "byu_names", names.data, names.used);
    free(names.data);
    load_constants(argv[1]);
    const char *files[64];
    int file_count = 0;
    files[file_count++] = "layout.txt";
    files[file_count++] = "constants.txt";
    for (int i = 0; i < tables && file_count < 64; i++) files[file_count++] = etables[i].path;
    add_hashes(argv[1], files, file_count);

    uint64_t index_offset = align64(sizeof(blob_header));
    uint64_t offset = align64(index_offset + sizeof(blob_section) * (uint64_t)section_count);
    blob_section *index = xmalloc(sizeof *index * (size_t)section_count);
    memset(index, 0, sizeof *index * (size_t)section_count);
    for (int i = 0; i < section_count; i++) {
        memcpy(index[i].name, sections[i].name, BLOB_NAME_SIZE);
        index[i].kind = sections[i].kind;
        index[i].count = sections[i].count;
        index[i].offset = offset;
        index[i].bytes = sections[i].bytes;
        offset = align64(offset + sections[i].bytes);
    }
    blob_header header;
    memset(&header, 0, sizeof header);
    memcpy(header.magic, BLOB_MAGIC, 8);
    header.version = BLOB_VERSION;
    header.section_count = (uint32_t)section_count;
    header.index_offset = index_offset;
    header.file_size = offset;

    char temp[4096];
    snprintf(temp, sizeof temp, "%s.tmp", argv[2]);
    FILE *out = fopen(temp, "wb");
    if (!out) die("cannot write", temp);
    unsigned char zero[64] = {0};
    uint64_t at = 0;
    fwrite(&header, sizeof header, 1, out);
    at += sizeof header;
    fwrite(zero, 1, (size_t)(index_offset - at), out);
    at = index_offset;
    fwrite(index, sizeof *index, (size_t)section_count, out);
    at += sizeof *index * (uint64_t)section_count;
    for (int i = 0; i < section_count; i++) {
        fwrite(zero, 1, (size_t)(index[i].offset - at), out);
        at = index[i].offset;
        fwrite(sections[i].payload, 1, (size_t)sections[i].bytes, out);
        at += sections[i].bytes;
    }
    fwrite(zero, 1, (size_t)(header.file_size - at), out);
    if (fclose(out)) die("cannot write", temp);
    if (rename(temp, argv[2])) die("cannot rename", argv[2]);
    return 0;
}
