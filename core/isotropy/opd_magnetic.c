/* Magnetic subgroups in the OPD rows: labels, index, ferroic classes and operation records. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"
#include "sel.h"
#include "irreps.h"
#include "tables.h"

int magnetic_records(int group, int setting, int capacity, int64_t *records_out);

static struct {
    section_view op_r, mag2nonmag, bns, nlabel, point_op;
    int ready;
} T;

int sel_mag_tables(void)
{
    if (T.ready) return 0;
    struct {
        const char *name;
        section_view *view;
    } list[] = {{"group.magnetic_ops.reversed", &T.op_r},
                {"group.magnetic_ops.ordinary", &T.mag2nonmag},
                {"group.magnetic_groups.bns", &T.bns},
                {"group.magnetic_groups.og", &T.nlabel},
                {"group.point_ops.matrix", &T.point_op}};
    if (sel_tables() || mag_tables()) return 1;
    for (size_t i = 0; i < sizeof list / sizeof list[0]; i++)
        if (section_find(list[i].name, list[i].view)) return 1;
    T.ready = 1;
    return 0;
}

static void put_label(jbuf *b, const section_view *view, int64_t index)
{
    const char *s;
    int64_t n;
    if (section_string(view, index, &s, &n)) {
        jb_lit(b, "\"\"");
        return;
    }
    while (n && (*s == ' ' || *s == '\t')) {
        s++;
        n--;
    }
    while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    jb_string(b, s, (size_t)n);
}

int sel_mag_ordinary(int group) { return mag_ordinary_group(group); }

/* The magnetic subgroup dict {number, symbol, ordinary_number, display_label}. */
void sel_mag_subgroup_json(jbuf *b, int group)
{
    jb_lit(b, "{\"number\": ");
    jb_int(b, group);
    jb_lit(b, ", \"symbol\": ");
    put_label(b, &T.bns, group - 1);
    jb_lit(b, ", \"ordinary_number\": ");
    jb_int(b, sel_mag_ordinary(group));
    jb_lit(b, ", \"display_label\": ");
    put_label(b, &T.nlabel, group - 1);
    jb_lit(b, "}");
}
int64_t sel_mag_index(int parent_group, int subgroup_group, int64_t size, int *has)
{
    int64_t parent_order = mag_pg_order(mag_group_pg(parent_group)), sub_order = mag_pg_order(mag_group_pg(subgroup_group));
    int64_t value = size * parent_order;
    *has = sub_order > 0 && value % sub_order == 0;
    return *has ? value / sub_order : 0;
}

static int rank3(int64_t (*m)[3], int rows)
{
    int rank = 0;
    for (int col = 0; col < 3 && rank < rows; col++) {
        int pivot = -1;
        for (int r = rank; r < rows; r++)
            if (m[r][col]) {
                pivot = r;
                break;
            }
        if (pivot < 0) continue;
        for (int c = 0; c < 3; c++) {
            int64_t t = m[rank][c];
            m[rank][c] = m[pivot][c];
            m[pivot][c] = t;
        }
        for (int r = 0; r < rows; r++) {
            if (r == rank || !m[r][col]) continue;
            int64_t a = m[rank][col], f = m[r][col];
            int64_t g = 0;
            for (int c = 0; c < 3; c++) {
                m[r][c] = m[r][c] * a - f * m[rank][c];
                int64_t v = m[r][c] < 0 ? -m[r][c] : m[r][c];
                while (v) {
                    int64_t t = g % v;
                    g = v;
                    v = t;
                }
            }
            if (g > 1)
                for (int c = 0; c < 3; c++) m[r][c] /= g;
        }
        rank++;
    }
    return rank;
}

int sel_mag_ferroic(jbuf *b, int group)
{
    static signed char cache[256];
    int64_t total = mag_group_total();
    if (group < 1 || group > total) return md_fail(MDE_VALUE, "unknown magnetic group: %d", group);
    int64_t pg = mag_group_pg(group);
    if (pg < 1 || pg > 255) return md_fail(MDE_VALUE, "unknown magnetic point group: %lld", (long long)pg);
    if (!cache[pg]) {
        int representative = 0;
        for (int64_t g = 0; g < total && !representative; g++)
            if (mag_group_pg((int)g + 1) == pg) representative = (int)g + 1;
        if (!representative) return md_fail(MDE_VALUE, "unknown magnetic point group: %lld", (long long)pg);
        int64_t records[192][5];
        int count = magnetic_records(representative, 1, 192, &records[0][0]);
        if (count < 0) {
            md_err.kind = MDE_MAGNETIC;
            md_err.status = -count;
            md_err.message[0] = 0;
            return MDE_MAGNETIC;
        }
        int64_t expected = mag_pg_order((int)pg);
        if (count != expected)
            return md_fail(MDE_VALUE, "magnetic point-group closure is incomplete: group=%d expected=%lld actual=%d", representative,
                           (long long)expected, count);
        int64_t (*constraints)[3] = malloc(sizeof(int64_t) * 3 * 3 * (size_t)(count ? count : 1));
        if (!constraints) return md_fail(MDE_MEMORY, "out of memory");
        for (int i = 0; i < count; i++) {
            int op = (int)records[i][4];
            int64_t nonmag = section_int(&T.mag2nonmag, op - 1);
            int64_t r[9];
            for (int k = 0; k < 9; k++) r[k] = section_int(&T.point_op, (nonmag - 1) * 9 + k);
            int64_t det = r[0] * (r[4] * r[8] - r[5] * r[7]) - r[1] * (r[3] * r[8] - r[5] * r[6]) + r[2] * (r[3] * r[7] - r[4] * r[6]);
            int64_t axial = (section_int(&T.op_r, op - 1) ? -1 : 1) * det;
            for (int row = 0; row < 3; row++)
                for (int col = 0; col < 3; col++) constraints[3 * i + row][col] = axial * r[3 * row + col] - (row == col);
        }
        cache[pg] = rank3(constraints, 3 * count) < 3 ? 1 : 2;
        free(constraints);
    }
    jb_lit(b, cache[pg] == 1 ? "[\"ferromagnetic\"]" : "[]");
    return 0;
}

static void put_ints(jbuf *b, const int64_t *values, int64_t n)
{
    jb_lit(b, "[");
    for (int64_t i = 0; i < n; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, values[i]);
    }
    jb_lit(b, "]");
}

void sel_put_records(jbuf *b, const int64_t *records, int64_t count)
{
    jb_lit(b, "[");
    for (int64_t i = 0; i < count; i++) {
        if (i) jb_lit(b, ", ");
        put_ints(b, records + 5 * i, 5);
    }
    jb_lit(b, "]");
}
