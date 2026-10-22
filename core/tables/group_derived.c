/* Class codes and orders computed from the point operation matrices:
 * - group.point_ops.code: 3^(tr R + 3 + 6 [det R = 1]), one power of 3 for each kind of rotation;
 * - group.point_groups.order and code: the closure of the group's generators, and the sum of its operations' codes;
 * - group.magnetic_ops.code: a number for each kind of the ordinary rotation, times 1000 with time reversal. */

#include <stdlib.h>
#include <string.h>

#include "tables.h"

int section_find_stored(const char *name, section_view *view);

enum { NAMES = 4 };
static const char *const names[NAMES] = {"group.point_ops.code", "group.point_groups.order", "group.point_groups.code",
                                         "group.magnetic_ops.code"};
/* by det R and tr R: det -1 with tr -3 .. 1, then det 1 with tr -1 .. 3 */
static const int64_t MAGNETIC_KIND[10] = {1, 24, 12, 16, 4, 2, 14, 6, 20, 500};

static struct {
    int ready, failed;
    int64_t *data[NAMES];
    blob_section entries[NAMES];
} C;

static int64_t det_trace(const int64_t *m, int64_t *trace)
{
    *trace = m[0] + m[4] + m[8];
    return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
}

static int build(void)
{
    section_view matrix, generators, generator_count, ordinary, reversed;
    if (section_find_stored("group.point_ops.matrix", &matrix) || section_find_stored("group.point_groups.generators", &generators) ||
        section_find_stored("group.point_groups.generator_count", &generator_count) ||
        section_find_stored("group.magnetic_ops.ordinary", &ordinary) || section_find_stored("group.magnetic_ops.reversed", &reversed))
        return 1;
    int64_t ops = section_count(&matrix) / 9, groups = section_count(&generators) / 5, mops = section_count(&ordinary);
    if (ops != 72) return 1;
    int64_t m[73][9], kind[73];
    for (int a = 1; a <= ops; a++) {
        int64_t tr;
        for (int i = 0; i < 9; i++) m[a][i] = section_int(&matrix, 9 * (a - 1) + i);
        int64_t det = det_trace(m[a], &tr);
        if ((det != 1 && det != -1) || tr < (det > 0 ? -1 : -3) || tr > (det > 0 ? 3 : 1)) return 1;
        kind[a] = det > 0 ? 5 + tr + 1 : tr + 3;
    }
    int64_t counts[NAMES] = {ops, groups, groups, mops};
    for (int i = 0; i < NAMES; i++) {
        if (!(C.data[i] = calloc((size_t)counts[i] + 1, sizeof(int64_t)))) return 1;
        strncpy(C.entries[i].name, names[i], BLOB_NAME_SIZE - 1);
        C.entries[i].kind = BLOB_I64;
        C.entries[i].count = (uint64_t)counts[i];
        C.entries[i].bytes = (uint64_t)counts[i] * sizeof(int64_t);
    }
    for (int a = 1; a <= ops; a++) {
        int64_t p = 1;
        for (int e = 0; e < kind[a] + (kind[a] >= 5 ? 3 : 0); e++) p *= 3;
        C.data[0][a - 1] = p;
    }
    for (int64_t g = 0; g < groups; g++) {
        int in[73] = {0}, list[73], size = 0, family = 0, count = (int)section_int(&generator_count, g);
        if (count < 1 || count > 5) return 1;
        for (int j = 0; j < count; j++) {
            int64_t op = section_int(&generators, 5 * g + j);
            if (op < 0 || op > ops) return 1;
            if (op) family = op > 48;
        }
        int identity = family ? 49 : 1;
        in[identity] = 1, list[size++] = identity;
        for (int j = 0; j < count; j++) {
            int op = (int)section_int(&generators, 5 * g + j);
            if (op && !in[op]) in[op] = 1, list[size++] = op;
        }
        for (int a = 0; a < size; a++)
            for (int b = 0; b < size; b++) {
                int64_t p[9];
                for (int i = 0; i < 3; i++)
                    for (int k = 0; k < 3; k++)
                        p[3 * i + k] = m[list[a]][3 * i] * m[list[b]][k] + m[list[a]][3 * i + 1] * m[list[b]][3 + k] +
                                       m[list[a]][3 * i + 2] * m[list[b]][6 + k];
                int found = 0;
                for (int c = family ? 49 : 1; c <= (family ? 72 : 48) && !found; c++)
                    if (!memcmp(p, m[c], sizeof p)) found = c;
                if (!found) return 1;
                if (!in[found]) in[found] = 1, list[size++] = found;
            }
        C.data[1][g] = size;
        for (int a = 0; a < size; a++) C.data[2][g] += C.data[0][list[a] - 1];
    }
    for (int64_t o = 0; o < mops; o++) {
        int64_t op = section_int(&ordinary, o);
        if (op < 1 || op > ops) return 1;
        C.data[3][o] = MAGNETIC_KIND[kind[op]] * (section_int(&reversed, o) ? 1000 : 1);
    }
    return 0;
}

/* section_find for the sections above; *handled 0 for any other name */
int group_derived_section(const char *name, section_view *view, int *handled)
{
    int which = -1;
    for (int i = 0; i < NAMES && which < 0; i++)
        if (!strncmp(name, names[i], BLOB_NAME_SIZE)) which = i;
    *handled = which >= 0;
    if (!*handled) return 1;
    if (!C.ready && !C.failed) {
        C.failed = build();
        C.ready = !C.failed;
    }
    if (!C.ready) return 1;
    view->entry = &C.entries[which];
    view->data = (const unsigned char *)C.data[which];
    return 0;
}
