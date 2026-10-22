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

#include "rational.h"
#include "tables.h"

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

static void decode(const int64_t *raw, rat (*f)[3])
{
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 3; c++) f[r][c] = raw[4 * r + 3] ? rat_make(raw[4 * r + c], raw[4 * r + 3]) : rat_int(0);
}

static void encode(rat (*f)[3], int64_t *raw)
{
    for (int r = 0; r < 4; r++) {
        i128 d = 1;
        int any = 0;
        for (int c = 0; c < 3; c++) {
            d = d / rat_gcd(d, f[r][c].d) * f[r][c].d;
            any |= !rat_zero(f[r][c]);
        }
        if (r && !any) d = 0;
        for (int c = 0; c < 3; c++) raw[4 * r + c] = d ? (int64_t)(f[r][c].n * (d / f[r][c].d)) : 0;
        raw[4 * r + 3] = (int64_t)d;
    }
}

static void row_mul(const rat *v, rat (*m)[3], rat *out)
{
    rat t[3];
    for (int c = 0; c < 3; c++) {
        t[c] = rat_int(0);
        for (int k = 0; k < 3; k++) t[c] = rat_add(t[c], rat_mul(v[k], m[k][c]));
    }
    memcpy(out, t, sizeof t);
}

static int inverse(rat (*m)[3], rat (*out)[3])
{
    rat det = rat_int(0);
    for (int c = 0; c < 3; c++) {
        rat minor = rat_sub(rat_mul(m[1][(c + 1) % 3], m[2][(c + 2) % 3]), rat_mul(m[1][(c + 2) % 3], m[2][(c + 1) % 3]));
        det = rat_add(det, rat_mul(m[0][c], minor));
    }
    if (rat_zero(det)) return 1;
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            rat cof = rat_sub(rat_mul(m[(c + 1) % 3][(r + 1) % 3], m[(c + 2) % 3][(r + 2) % 3]),
                              rat_mul(m[(c + 1) % 3][(r + 2) % 3], m[(c + 2) % 3][(r + 1) % 3]));
            out[r][c] = rat_div(cof, det);
        }
    return 0;
}

/* PML coordinates to the setting's: x m + o */
static int setting(int64_t st, int64_t sg, rat (*m)[3], rat *o)
{
    int64_t base = (st - 1) * 32, den = section_int(&inter, base + 15), l = section_int(&lattice, sg - 1) - 1;
    int64_t pden = section_int(&ml_den, l * 2);
    if (!den || !pden) return 1;
    rat t[3][3], cml[3][3];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            t[r][c] = rat_make(section_int(&inter, base + 4 * r + c), den);
            cml[r][c] = rat_make(section_int(&ml, l * 36 + 9 + 3 * r + c), pden);
        }
    for (int r = 0; r < 3; r++) row_mul(cml[r], t, m[r]);
    for (int a = 0; a < 3; a++) o[a] = rat_make(section_int(&inter, base + 12 + a), den);
    return 0;
}

static int canonical(rat (*f)[3], int64_t st, int64_t sg, rat (*out)[3])
{
    rat m[3][3], o[3], im[3][3], d[4][3];
    if (setting(st, sg, m, o) || inverse(m, im)) return 1;
    for (int r = 0; r < 4; r++) row_mul(f[r], m, d[r]);
    for (int a = 0; a < 3; a++) d[0][a] = rat_add(d[0][a], o[a]);
    rat (*rows)[3] = d + 1;
    int pivots[3], rank = 0;
    for (int c = 0; c < 3 && rank < 3; c++) {
        int p = -1;
        for (int r = rank; r < 3 && p < 0; r++)
            if (!rat_zero(rows[r][c])) p = r;
        if (p < 0) continue;
        rat t[3];
        memcpy(t, rows[p], sizeof t);
        memcpy(rows[p], rows[rank], sizeof t);
        memcpy(rows[rank], t, sizeof t);
        rat lead = rows[rank][c];
        for (int k = 0; k < 3; k++) rows[rank][k] = rat_div(rows[rank][k], lead);
        for (int r = 0; r < 3; r++) {
            if (r == rank || rat_zero(rows[r][c])) continue;
            rat x = rows[r][c];
            for (int k = 0; k < 3; k++) rows[r][k] = rat_sub(rows[r][k], rat_mul(x, rows[rank][k]));
        }
        pivots[rank++] = c;
    }
    rat p[3];
    memcpy(p, d[0], sizeof p);
    for (int i = 0; i < rank; i++) {
        rat x = p[pivots[i]];
        for (int k = 0; k < 3; k++) p[k] = rat_sub(p[k], rat_mul(x, rows[i][k]));
    }
    for (int a = 0; a < 3; a++) p[a] = rat_sub(rat_mod1(p[a]), o[a]);
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 3; c++) out[r][c] = rat_int(0);
    row_mul(p, im, out[0]);
    for (int i = 0; i < rank; i++) row_mul(rows[i], im, out[pivots[i] + 1]);
    return 0;
}

static int64_t source_of(int64_t st)
{
    int64_t n = section_count(&number), sg = section_int(&number, st - 1), best = 0;
    if (section_int(&version, st - 1) != 1) return 0;
    for (int64_t s = 1; s <= n; s++) {
        int64_t v = section_int(&version, s - 1);
        if (s == st || section_int(&number, s - 1) != sg || (v != 1 && v != 3) ||
            section_int(&origin, s - 1) != section_int(&origin, st - 1) || section_int(&cell, s - 1) ||
            !text_is(&axis, s - 1, "") || !(text_is(&abc, s - 1, "abc") || text_is(&abc, s - 1, "")))
            continue;
        if (!best || (v == 3 && section_int(&version, best - 1) != 3)) best = s;
    }
    return best;
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
        for (int64_t w = 0; w < count; w++) {
            int64_t *raw = out + 16 * (slot[st] + w), def[16];
            for (int k = 0; k < 16; k++) def[k] = section_int(&wfract, 16 * (first + w) + k);
            if (section_int(&version, st - 1) == 2) {
                memcpy(raw, def, sizeof def);
                continue;
            }
            rat f[4][3], g[4][3];
            decode(src ? out + 16 * (slot[src] + w) : def, f);
            if (canonical(f, st, sg, g)) return 1;
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
    for (int64_t st = 1; st <= n; st++)
        if (derive(st, W.data[1])) return 1;
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
