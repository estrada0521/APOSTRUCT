/* settings.inter_settings.wyckoff and wyckoff_start: the formula of each Wyckoff position in every International Tables
 * setting, 16 numbers (the constant, then the x, y and z directions, each three numerators and a denominator).
 * A version-2 setting keeps the default formula.  Otherwise the formula is taken into the setting's coordinates, its
 * directions are row-reduced (each pivot row becoming the parameter of its pivot coordinate), the constant is cleared
 * at the pivot coordinates and reduced to [0,1) at the others, and the result is taken back.  A version-1 setting starts
 * from the formula of the setting of the same group and origin in the abc orientation (version 3 first), when there is
 * one, instead of the default.  settings.wyckoff_exceptions.formula lists (setting, position, formula) records that
 * replace the result: the representatives International Tables prints. */

#include <stdlib.h>
#include <string.h>

#include "tables.h"

typedef struct {
    int64_t n, d;
} q64;

static int overflow;

static int64_t gcd64(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        int64_t t = a % b;
        a = b, b = t;
    }
    return a;
}

static q64 q_make(int64_t n, int64_t d)
{
    if (d < 0) n = -n, d = -d;
    int64_t g = gcd64(n, d);
    if (g > 1) n /= g, d /= g;
    if (!d) overflow = 1, d = 1;
    return (q64){n, d};
}

static q64 q_int(int64_t v) { return (q64){v, 1}; }
static int q_zero(q64 a) { return a.n == 0; }

static q64 q_mul(q64 a, q64 b)
{
    int64_t n, d;
    if (__builtin_mul_overflow(a.n, b.n, &n) || __builtin_mul_overflow(a.d, b.d, &d)) return overflow = 1, q_int(0);
    return q_make(n, d);
}

static q64 q_add(q64 a, q64 b)
{
    int64_t x, y, n, d;
    if (__builtin_mul_overflow(a.n, b.d, &x) || __builtin_mul_overflow(b.n, a.d, &y) || __builtin_add_overflow(x, y, &n) ||
        __builtin_mul_overflow(a.d, b.d, &d))
        return overflow = 1, q_int(0);
    return q_make(n, d);
}

static q64 q_sub(q64 a, q64 b) { return q_add(a, (q64){-b.n, b.d}); }
static q64 q_div(q64 a, q64 b) { return q_mul(a, q_make(b.d, b.n)); }

static q64 q_mod1(q64 a)
{
    int64_t r = a.n % a.d;
    return (q64){r < 0 ? r + a.d : r, a.d};
}

int section_find_stored(const char *name, section_view *view);

static struct {
    int ready, failed;
    int64_t *data[2];
    blob_section entries[2];
} W;

static section_view number, version, origin, axis, abc, cell, inter, lattice, ml, ml_den, wstart, wcount, wfract, exceptions;
static int64_t *slot;
static signed char *state;

static int text_is(const section_view *v, int64_t i, const char *s)
{
    const char *t;
    int64_t n;
    if (section_string(v, i, &t, &n)) return 0;
    return (size_t)n == strlen(s) && !memcmp(t, s, (size_t)n);
}

static void decode(const int64_t *raw, q64 (*f)[3])
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 3; c++) f[r][c] = raw[4 * r + 3] ? q_make(raw[4 * r + c], raw[4 * r + 3]) : q_int(0);
}

static void encode(q64 (*f)[3], int64_t *raw)
{
    for (int r = 0; r < 4; r++) {
        int64_t d = 1;
        int any = 0;
        for (int c = 0; c < 3; c++) {
            d = d / gcd64(d, f[r][c].d) * f[r][c].d;
            any |= !q_zero(f[r][c]);
        }
        if (r && !any) d = 0;
        for (int c = 0; c < 3; c++) raw[4 * r + c] = d ? f[r][c].n * (d / f[r][c].d) : 0;
        raw[4 * r + 3] = d;
    }
}

static void row_mul(const q64 *v, q64 (*m)[3], q64 *out)
{
    q64 t[3];
    for (int c = 0; c < 3; c++) {
        t[c] = q_int(0);
        for (int k = 0; k < 3; k++) t[c] = q_add(t[c], q_mul(v[k], m[k][c]));
    }
    memcpy(out, t, sizeof t);
}

static int inverse(q64 (*m)[3], q64 (*out)[3])
{
    q64 det = q_int(0);
    for (int c = 0; c < 3; c++) {
        q64 minor = q_sub(q_mul(m[1][(c + 1) % 3], m[2][(c + 2) % 3]), q_mul(m[1][(c + 2) % 3], m[2][(c + 1) % 3]));
        det = q_add(det, q_mul(m[0][c], minor));
    }
    if (q_zero(det)) return 1;
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            q64 cof = q_sub(q_mul(m[(c + 1) % 3][(r + 1) % 3], m[(c + 2) % 3][(r + 2) % 3]),
                              q_mul(m[(c + 1) % 3][(r + 2) % 3], m[(c + 2) % 3][(r + 1) % 3]));
            out[r][c] = q_div(cof, det);
        }
    return 0;
}

/* PML coordinates to the setting's: x m + o */
static int setting(int64_t st, int64_t sg, q64 (*m)[3], q64 *o)
{
    int64_t base = (st - 1) * 32, den = section_int(&inter, base + 15), l = section_int(&lattice, sg - 1) - 1;
    int64_t pden = section_int(&ml_den, l * 2);
    if (!den || !pden) return 1;
    q64 t[3][3], cml[3][3];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            t[r][c] = q_make(section_int(&inter, base + 4 * r + c), den);
            cml[r][c] = q_make(section_int(&ml, l * 36 + 9 + 3 * r + c), pden);
        }
    for (int r = 0; r < 3; r++) row_mul(cml[r], t, m[r]);
    for (int a = 0; a < 3; a++) o[a] = q_make(section_int(&inter, base + 12 + a), den);
    return 0;
}

static int canonical(q64 (*f)[3], q64 (*m)[3], const q64 *o, q64 (*im)[3], q64 (*out)[3])
{
    q64 d[4][3];
    for (int r = 0; r < 4; r++) row_mul(f[r], m, d[r]);
    for (int a = 0; a < 3; a++) d[0][a] = q_add(d[0][a], o[a]);
    q64 (*rows)[3] = d + 1;
    int pivots[3], rank = 0;
    for (int c = 0; c < 3 && rank < 3; c++) {
        int p = -1;
        for (int r = rank; r < 3 && p < 0; r++)
            if (!q_zero(rows[r][c])) p = r;
        if (p < 0) continue;
        q64 t[3];
        memcpy(t, rows[p], sizeof t);
        memcpy(rows[p], rows[rank], sizeof t);
        memcpy(rows[rank], t, sizeof t);
        q64 lead = rows[rank][c];
        for (int k = 0; k < 3; k++) rows[rank][k] = q_div(rows[rank][k], lead);
        for (int r = 0; r < 3; r++) {
            if (r == rank || q_zero(rows[r][c])) continue;
            q64 x = rows[r][c];
            for (int k = 0; k < 3; k++) rows[r][k] = q_sub(rows[r][k], q_mul(x, rows[rank][k]));
        }
        pivots[rank++] = c;
    }
    q64 p[3];
    memcpy(p, d[0], sizeof p);
    for (int i = 0; i < rank; i++) {
        q64 x = p[pivots[i]];
        for (int k = 0; k < 3; k++) p[k] = q_sub(p[k], q_mul(x, rows[i][k]));
    }
    for (int a = 0; a < 3; a++) p[a] = q_sub(q_mod1(p[a]), o[a]);
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 3; c++) out[r][c] = q_int(0);
    row_mul(p, im, out[0]);
    for (int i = 0; i < rank; i++) row_mul(rows[i], im, out[pivots[i] + 1]);
    return 0;
}

/* per group and origin, the first two abc settings of version 3 or 1, version 3 first */
static int64_t abc_settings[231][3][2];

static int better(int64_t s, int64_t than)
{
    return !than || (section_int(&version, s - 1) == 3 && section_int(&version, than - 1) != 3);
}

static void find_abc_settings(void)
{
    for (int64_t s = 1; s <= section_count(&number); s++) {
        int64_t sg = section_int(&number, s - 1), v = section_int(&version, s - 1), o = section_int(&origin, s - 1);
        if (sg < 1 || sg > 230 || o < 0 || o > 2 || (v != 1 && v != 3) || section_int(&cell, s - 1) ||
            !text_is(&axis, s - 1, "") || !(text_is(&abc, s - 1, "abc") || text_is(&abc, s - 1, "")))
            continue;
        int64_t *best = abc_settings[sg][o];
        if (better(s, best[0])) best[1] = best[0], best[0] = s;
        else if (better(s, best[1])) best[1] = s;
    }
}

static int64_t source_of(int64_t st)
{
    int64_t sg = section_int(&number, st - 1), o = section_int(&origin, st - 1);
    if (section_int(&version, st - 1) != 1 || o < 0 || o > 2) return 0;
    const int64_t *best = abc_settings[sg][o];
    return best[0] == st ? best[1] : best[0];
}

static int derive(int64_t st, int64_t *out)
{
    if (state[st] == 2) return 0;
    if (state[st] == 1) return 1;
    state[st] = 1;
    int64_t sg = section_int(&number, st - 1);
    if (sg >= 1 && sg <= 230) {
        int64_t src = source_of(st), count = section_int(&wcount, sg - 1), first = section_int(&wstart, sg - 1) - 1;
        if (src && derive(src, out)) return 1;
        q64 m[3][3], o[3], im[3][3];
        if (section_int(&version, st - 1) != 2 && (setting(st, sg, m, o) || inverse(m, im))) return 1;
        for (int64_t w = 0; w < count; w++) {
            int64_t *raw = out + 16 * (slot[st] + w), def[16];
            for (int k = 0; k < 16; k++) def[k] = section_int(&wfract, 16 * (first + w) + k);
            if (section_int(&version, st - 1) == 2) {
                memcpy(raw, def, sizeof def);
                continue;
            }
            q64 f[4][3], g[4][3];
            decode(src ? out + 16 * (slot[src] + w) : def, f);
            if (canonical(f, m, o, im, g)) return 1;
            encode(g, raw);
        }
        for (int64_t e = 0; e + 18 <= section_count(&exceptions); e += 18) {
            int64_t pos = section_int(&exceptions, e + 1);
            if (section_int(&exceptions, e) != st) continue;
            if (pos < 1 || pos > count) return 1;
            for (int k = 0; k < 16; k++) out[16 * (slot[st] + pos - 1) + k] = section_int(&exceptions, e + 2 + k);
        }
    }
    state[st] = 2;
    return 0;
}

static int build(void)
{
    if (section_find_stored("settings.inter_settings.number", &number) ||
        section_find_stored("settings.inter_settings.version", &version) ||
        section_find_stored("settings.inter_settings.origin", &origin) ||
        section_find_stored("settings.inter_settings.axis", &axis) || section_find_stored("settings.inter_settings.abc", &abc) ||
        section_find_stored("settings.inter_settings.cell", &cell) || section_find_stored("settings.inter_settings.inter", &inter) ||
        section_find_stored("group.space_groups.lattice", &lattice) || section_find_stored("group.lattices.ml", &ml) ||
        section_find_stored("group.lattices.ml_denom", &ml_den) || section_find_stored("settings.wyckoff.by_sg_start", &wstart) ||
        section_find_stored("settings.wyckoff.by_sg_count", &wcount) || section_find_stored("settings.wyckoff.fract", &wfract) ||
        section_find_stored("settings.wyckoff_exceptions.formula", &exceptions))
        return 1;
    int64_t n = section_count(&number), total = 0;
    slot = calloc((size_t)n + 1, sizeof *slot);
    state = calloc((size_t)n + 1, 1);
    W.data[0] = calloc((size_t)n + 1, sizeof(int64_t));
    if (!slot || !state || !W.data[0]) return 1;
    for (int64_t st = 1; st <= n; st++) {
        int64_t sg = section_int(&number, st - 1), count = sg >= 1 && sg <= 230 ? section_int(&wcount, sg - 1) : 0;
        slot[st] = total;
        W.data[0][st - 1] = count ? total + 1 : 0;
        total += count;
    }
    if (!(W.data[1] = calloc((size_t)(16 * total + 1), sizeof(int64_t)))) return 1;
    find_abc_settings();
    for (int64_t st = 1; st <= n; st++)
        if (derive(st, W.data[1]) || overflow) return 1;
    free(slot);
    free(state);
    const char *names[2] = {"settings.inter_settings.wyckoff_start", "settings.inter_settings.wyckoff"};
    int64_t counts[2] = {n, 16 * total};
    for (int i = 0; i < 2; i++) {
        strncpy(W.entries[i].name, names[i], BLOB_NAME_SIZE - 1);
        W.entries[i].kind = BLOB_I64;
        W.entries[i].count = (uint64_t)counts[i];
        W.entries[i].bytes = (uint64_t)counts[i] * sizeof(int64_t);
    }
    return 0;
}

/* section_find for the two sections above; *handled 0 for any other name */
int settings_derived_section(const char *name, section_view *view, int *handled)
{
    int which = !strncmp(name, "settings.inter_settings.wyckoff_start", BLOB_NAME_SIZE) ? 0
              : !strncmp(name, "settings.inter_settings.wyckoff", BLOB_NAME_SIZE)       ? 1
                                                                                         : -1;
    *handled = which >= 0;
    if (!*handled) return 1;
    if (!W.ready && !W.failed) {
        W.failed = build();
        W.ready = !W.failed;
    }
    if (!W.ready) return 1;
    view->entry = &W.entries[which];
    view->data = (const unsigned char *)W.data[which];
    return 0;
}
