/* settings.wyckoff.site_pg: the point group of each Wyckoff position, from its formula and its space group. The
 * operations kept are those fixing the formula identically (the constant modulo translations, each parameter direction
 * exactly); the point group is the one whose code (group.point_groups.code) is the sum of their rotations' codes.
 * apo asks for single rows (wyckoff_site_pg); the whole section is built only when it is read as a section. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "tables.h"

static struct {
    int ready, failed;
    int64_t rows, *data;
    blob_section entry;
} W;

static int64_t row_count(void)
{
    if (W.rows) return W.rows;
    for (int sg = 1; sg <= 230; sg++) {
        int first, count;
        if (geo_wyckoff_row_range(sg, &first, &count)) return 0;
        W.rows += count;
    }
    if (!(W.data = calloc((size_t)W.rows + 1, sizeof(int64_t)))) W.rows = 0;
    return W.rows;
}

/* the point groups of every Wyckoff row of space group sg; 0 on failure */
static int derive_group(int sg)
{
    int first, count;
    int64_t ops[GEO_MAX_OPS][5], ncodes;
    int nops = geo_sg_records(sg, ops);
    const int64_t *codes = geo_point_group_codes(&ncodes);
    if (geo_wyckoff_row_range(sg, &first, &count) || nops <= 0 || !codes) return 0;
    for (int row = first; row < first + count; row++) {
        rat v[4][3], t[3];
        int64_t code = 0;
        geo_wyckoff_vectors(row, v);
        for (int j = 0; j < nops; j++) {
            int keep = 1;
            geo_vrot(sg, (int)ops[j][4], v[0], t);
            for (int a = 0; a < 3 && keep; a++)
                keep = rat_zero(rat_mod1(rat_sub(rat_sub(v[0][a], t[a]), rat_make(ops[j][a], ops[j][3]))));
            for (int d = 1; d < 4 && keep; d++) {
                geo_vrot(sg, (int)ops[j][4], v[d], t);
                for (int a = 0; a < 3 && keep; a++) keep = rat_eq(t[a], v[d][a]);
            }
            if (keep) code += geo_point_op_code((int)ops[j][4]);
        }
        int found = 0;
        for (int64_t g = 0; g < ncodes; g++)
            if (codes[g] == code) W.data[row - 1] = g + 1, found++;
        if (found != 1) return 0;
    }
    return 1;
}

/* the point group of Wyckoff row row_id (1-based), 0 when it cannot be derived */
int wyckoff_site_pg(int row_id)
{
    if (geo_load() || row_id < 1 || row_id > row_count()) return 0;
    if (W.data[row_id - 1]) return (int)W.data[row_id - 1];
    int first, count;
    for (int sg = 1; sg <= 230; sg++)
        if (!geo_wyckoff_row_range(sg, &first, &count) && row_id < first + count) return derive_group(sg) ? (int)W.data[row_id - 1] : 0;
    return 0;
}

/* section_find for the section above; *handled 0 for any other name */
int wyckoff_derived_section(const char *name, section_view *view, int *handled)
{
    *handled = !strncmp(name, "settings.wyckoff.site_pg", BLOB_NAME_SIZE);
    if (!*handled) return 1;
    if (!W.ready && !W.failed) {
        int64_t rows = row_count();
        W.failed = !rows;
        for (int64_t r = 1; r <= rows && !W.failed; r++) W.failed = !wyckoff_site_pg((int)r);
        strncpy(W.entry.name, "settings.wyckoff.site_pg", BLOB_NAME_SIZE - 1);
        W.entry.kind = BLOB_I64;
        W.entry.count = (uint64_t)rows;
        W.entry.bytes = (uint64_t)rows * sizeof(int64_t);
        W.ready = !W.failed;
    }
    if (!W.ready) return 1;
    view->entry = &W.entry;
    view->data = (const unsigned char *)W.data;
    return 0;
}
