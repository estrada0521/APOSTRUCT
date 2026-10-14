/* definition_presentation._apply_dynamic_source_family_presentation: the
 * public axis order of Source carrier tensors.  Modes of one print family
 * (spec, site, site irrep and, unless the family is one multi-row carrier,
 * the site-irrep component) are reordered by the carrier's tensor rule. */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { KIND_DSP = 1, KIND_MAG = 2 };

/* Per mode, 18 ints: kind, has family, family, has component, print
 * component, family width, dynamic, phase, little type, source row count,
 * opd group count, opd row count, every group has one row, component equals
 * print component, spec order, site order, site pg, pg irrep, and the
 * site-irrep component id (-2 none, -1 empty).  opd row count is 0 unless
 * every opd group is a mapping. */
enum {
    F_KIND, F_HAS_FAMILY, F_FAMILY, F_HAS_COMPONENT, F_COMPONENT, F_WIDTH, F_DYNAMIC, F_PHASE, F_LITTLE_TYPE,
    F_SOURCE_ROWS, F_GROUPS, F_GROUP_ROWS, F_SINGLE_ROW_GROUPS, F_COMPONENT_SAME, F_SPEC, F_SITE, F_SITE_PG,
    F_PG_IRREP, F_SITE_COMPONENT, FIELDS
};

typedef struct {
    const int64_t *m;
    int64_t key[5];
    int index;
} entry;

static int64_t py_floordiv(int64_t a, int64_t b)
{
    int64_t q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) q--;
    return q;
}

static int64_t py_mod(int64_t a, int64_t b)
{
    int64_t r = a % b;
    if (r != 0 && ((r < 0) != (b < 0))) r += b;
    return r;
}

static int64_t type3_rows(const int64_t *m)
{
    int64_t rows = m[F_GROUP_ROWS], width = m[F_WIDTH];
    if (!m[F_GROUPS]) return 0;
    if (m[F_LITTLE_TYPE] == 3 && rows >= 2 && width >= rows && py_mod(width, rows) == 0) return rows;
    return 0;
}

static int64_t type2_rows(const int64_t *m)
{
    int64_t rows = m[F_SOURCE_ROWS], width = m[F_WIDTH];
    if (m[F_LITTLE_TYPE] == 2 && rows > 2 && rows % 2 == 0 && py_mod(width, rows) == 0 && m[F_GROUPS] == rows &&
        m[F_SINGLE_ROW_GROUPS])
        return rows;
    return 0;
}

static int64_t type1_rows(const int64_t *m)
{
    int64_t groups = m[F_GROUPS], width = m[F_WIDTH];
    if (m[F_LITTLE_TYPE] == 1 && groups > 1 && py_mod(width, groups) == 0 && m[F_SINGLE_ROW_GROUPS]) return groups;
    return 0;
}

static int cmp_entries(const void *left, const void *right)
{
    const entry *a = left, *b = right;
    for (int i = 0; i < 5; i++)
        if (a->key[i] != b->key[i]) return a->key[i] < b->key[i] ? -1 : 1;
    return a->index - b->index;
}

/* The distinct value count of one field over a set of modes. */
static int distinct(const entry *e, int count, int field, int64_t *first)
{
    int n = 0;
    for (int i = 0; i < count; i++) {
        int seen = 0;
        for (int j = 0; j < i && !seen; j++) seen = e[j].m[field] == e[i].m[field];
        if (!seen) {
            if (!n) *first = e[i].m[field];
            n++;
        }
    }
    return n;
}

static int64_t mode_order_key[3 * 4096];

static int cmp_sorted(const void *left, const void *right)
{
    const int *a = left, *b = right;
    for (int i = 0; i < 3; i++) {
        int64_t x = mode_order_key[3 * a[1] + i], y = mode_order_key[3 * b[1] + i];
        if (x != y) return x < y ? -1 : 1;
    }
    return a[1] - b[1];
}

static void order_group(const entry *e, int count, int *perm)
{
    int64_t width = 0;
    int widths = distinct(e, count, F_WIDTH, &width);
    int families = 0;
    int64_t family_set[4096];
    for (int i = 0; i < count; i++) {
        int64_t f = e[i].m[F_FAMILY];
        int seen = 0;
        for (int j = 0; j < families && !seen; j++) seen = family_set[j] == f;
        if (!seen && families < 4096) family_set[families++] = f;
    }
    int families_are_range = 1;
    for (int j = 0; j < families && families_are_range; j++) families_are_range = family_set[j] >= 0 && family_set[j] < families;
    int distinct_pairs = 0;
    for (int i = 0; i < count; i++) {
        int seen = 0;
        for (int j = 0; j < i && !seen; j++)
            seen = e[j].m[F_FAMILY] == e[i].m[F_FAMILY] && e[j].m[F_COMPONENT] == e[i].m[F_COMPONENT];
        if (!seen) distinct_pairs++;
    }
    int components_in_range = 1;
    for (int i = 0; i < count && components_in_range; i++)
        components_in_range = e[i].m[F_COMPONENT] >= 0 && e[i].m[F_COMPONENT] < width;
    int rectangular = widths == 1 && width > 1 && families_are_range && components_in_range &&
                      distinct_pairs == families * width;
    int all_dsp = 1, all_mag = 1, any_dynamic = 0, all_dynamic = 1, any_phase = 0, all_phase = 1, all_type2 = 1;
    int all_multi_group = 1, any_type3_little = 0, any_dynamic_dsp = 0, kinds_mixed = 0, all_component_same = 1;
    for (int i = 0; i < count; i++) {
        const int64_t *m = e[i].m;
        all_dsp &= m[F_KIND] == KIND_DSP;
        all_mag &= m[F_KIND] == KIND_MAG;
        if (m[F_KIND] != e[0].m[F_KIND]) kinds_mixed = 1;
        any_dynamic |= m[F_DYNAMIC] != 0;
        all_dynamic &= m[F_DYNAMIC] != 0;
        any_phase |= m[F_PHASE] != 0;
        all_phase &= m[F_PHASE] != 0;
        all_type2 &= type2_rows(m) > 0;
        any_dynamic_dsp |= m[F_KIND] == KIND_DSP && m[F_DYNAMIC];
        all_multi_group &= m[F_GROUPS] > 1;
        any_type3_little |= m[F_LITTLE_TYPE] == 3;
        all_component_same &= m[F_COMPONENT_SAME] != 0;
    }
    int one_kind = !kinds_mixed;
    int static_rect_dsp = rectangular && all_dsp && !any_dynamic;
    int static_rect_mag = rectangular && all_mag && !any_dynamic;
    int64_t t1 = type1_rows(e[0].m), t2 = type2_rows(e[0].m), t3 = type3_rows(e[0].m);
    int t1_same = 1, t2_same = 1, t3_same = 1, source_rows_same = 1;
    for (int i = 1; i < count; i++) {
        t1_same &= type1_rows(e[i].m) == t1;
        t2_same &= type2_rows(e[i].m) == t2;
        t3_same &= type3_rows(e[i].m) == t3;
        source_rows_same &= e[i].m[F_SOURCE_ROWS] == e[0].m[F_SOURCE_ROWS];
    }
    int type3_carrier = rectangular && t3_same && t3 != 0;
    int64_t type1_source_rows = t1_same && source_rows_same && t1 == e[0].m[F_SOURCE_ROWS] ? t1 : 0;
    int dynamic_type1 = rectangular && one_kind && all_dynamic && width >= 4 && width % 2 == 0 && type1_source_rows > 1 &&
                        (all_dsp || (all_mag && all_phase && (width != 4 || families > 1)));
    int type2_order = rectangular && one_kind && families > 1 && count == families * width && all_component_same &&
                      t2_same && t2 != 0;
    int dynamic_type2_mag = all_mag && all_dynamic && !any_phase && all_type2;
    int static_type1 = (static_rect_dsp || static_rect_mag) && t1_same && t1 != 0;
    int single_family_multi_group = all_dsp && !any_dynamic && !any_type3_little && families == 1 && widths == 1 &&
                                    components_in_range && distinct_pairs == width && all_multi_group;
    int mode = 0;
    int64_t divisor = 1;
    if (type3_carrier) {
        mode = 1;
        divisor = t3;
    } else if (static_type1) {
        mode = 2;
        divisor = width / t1;
    } else if (dynamic_type1) {
        mode = 3;
        divisor = width / type1_source_rows;
    } else if (type2_order) {
        mode = 1;
        divisor = t2;
    } else if (dynamic_type2_mag || single_family_multi_group) {
        mode = 0;
    } else if (static_rect_mag || static_rect_dsp || any_phase || any_dynamic_dsp) {
        mode = 4;
    }
    int order[4096][2];
    for (int i = 0; i < count; i++) {
        int64_t c = e[i].m[F_COMPONENT], f = e[i].m[F_FAMILY], *k = mode_order_key + 3 * i;
        if (mode == 1) {
            k[0] = py_floordiv(c, divisor);
            k[1] = f;
            k[2] = py_mod(c, divisor);
        } else if (mode == 2 || mode == 3) {
            k[0] = py_mod(c, divisor);
            k[1] = f;
            k[2] = py_floordiv(c, divisor);
        } else if (mode == 4) {
            k[0] = c;
            k[1] = f;
            k[2] = 0;
        } else {
            k[0] = k[1] = k[2] = 0;
        }
        order[i][0] = e[i].index;
        order[i][1] = i;
    }
    qsort(order, (size_t)count, sizeof order[0], cmp_sorted);
    for (int i = 0; i < count; i++) perm[e[i].index] = e[order[i][1]].index;
}

/* perm[target] = the mode whose vectors, plan, normfactor and identity the
 * target takes (perm[i] = i for modes that keep their own). */
int family_presentation_order(int count, const int64_t *modes, int *perm)
{
    for (int i = 0; i < count; i++) perm[i] = i;
    if (count <= 0) return 0;
    entry *bases = malloc(sizeof(entry) * (size_t)count);
    entry *groups = malloc(sizeof(entry) * (size_t)count);
    if (!bases || !groups) {
        free(bases);
        free(groups);
        return 1;
    }
    int base_count = 0;
    for (int i = 0; i < count; i++) {
        const int64_t *m = modes + (size_t)FIELDS * i;
        if (!m[F_HAS_FAMILY]) continue;
        bases[base_count++] = (entry){m, {m[F_SPEC], m[F_SITE], m[F_SITE_PG], m[F_PG_IRREP], m[F_FAMILY]}, i};
    }
    qsort(bases, (size_t)base_count, sizeof *bases, cmp_entries);
    int group_count = 0;
    for (int i = 0; i < count; i++) {
        const int64_t *m = modes + (size_t)FIELDS * i;
        if (!m[F_HAS_FAMILY] || !m[F_HAS_COMPONENT]) continue;
        int dynamic_type2_mag = m[F_KIND] == KIND_MAG && m[F_DYNAMIC] && !m[F_PHASE] && type2_rows(m) > 0;
        int type3 = type3_rows(m) > 0;
        int static_pair = (m[F_KIND] == KIND_DSP || m[F_KIND] == KIND_MAG) && !m[F_DYNAMIC] && m[F_WIDTH] > 1;
        if (!m[F_PHASE] && !(m[F_KIND] == KIND_DSP && m[F_DYNAMIC]) && !dynamic_type2_mag && !type3 && !static_pair)
            continue;
        if (m[F_SITE_COMPONENT] == -2) continue;
        int64_t families = 0;
        for (int b = 0; b < base_count; b++) {
            const entry *x = bases + b;
            if (x->key[0] != m[F_SPEC] || x->key[1] != m[F_SITE] || x->key[2] != m[F_SITE_PG] || x->key[3] != m[F_PG_IRREP])
                continue;
            if (b == 0 || memcmp(bases[b - 1].key, x->key, sizeof x->key)) families++;
        }
        int multi_group = m[F_GROUPS] > 1;
        int merged = (multi_group && (m[F_LITTLE_TYPE] == 3 || families == 1)) || (type3 && families > 1);
        groups[group_count++] = (entry){m, {m[F_SPEC], m[F_SITE], m[F_SITE_PG], m[F_PG_IRREP], merged ? -1 : m[F_SITE_COMPONENT]}, i};
    }
    qsort(groups, (size_t)group_count, sizeof *groups, cmp_entries);
    int status = 0;
    for (int start = 0; start < group_count && !status;) {
        int end = start + 1;
        while (end < group_count && !memcmp(groups[end].key, groups[start].key, sizeof groups[start].key)) end++;
        if (end - start > 4096) status = 2;
        else order_group(groups + start, end - start, perm);
        start = end;
    }
    free(bases);
    free(groups);
    return status;
}
