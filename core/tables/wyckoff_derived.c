/* settings.wyckoff.site_pg: the point group of each Wyckoff position, from its formula and its space group. The
 * operations kept are those fixing the formula identically (the constant modulo translations, each parameter direction
 * exactly); the point group is the one whose code (group.point_groups.code) is the sum of their rotations' codes. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "tables.h"

int section_find_stored(const char *name, section_view *view);

static struct {
    int ready, failed;
    int64_t *data;
    blob_section entry;
} W;

static int build(void)
{
    section_view by_sg_count;
    if (geo_load() || section_find_stored("settings.wyckoff.by_sg_count", &by_sg_count)) return 1;
    int64_t rows = 0, ncodes;
    for (int sg = 1; sg <= 230; sg++) rows += section_int(&by_sg_count, sg - 1);
    const int64_t *codes = geo_point_group_codes(&ncodes);
    if (!codes || !(W.data = calloc((size_t)rows + 1, sizeof(int64_t)))) return 1;
    strncpy(W.entry.name, "settings.wyckoff.site_pg", BLOB_NAME_SIZE - 1);
    W.entry.kind = BLOB_I64;
    W.entry.count = (uint64_t)rows;
    W.entry.bytes = (uint64_t)rows * sizeof(int64_t);
    for (int sg = 1; sg <= 230; sg++) {
        int first, count;
        int64_t ops[GEO_MAX_OPS][5];
        int nops = geo_sg_records(sg, ops);
        if (geo_wyckoff_row_range(sg, &first, &count) || nops <= 0) return 1;
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
            if (found != 1) return 1;
        }
    }
    return 0;
}

/* section_find for the section above; *handled 0 for any other name */
int wyckoff_derived_section(const char *name, section_view *view, int *handled)
{
    *handled = !strncmp(name, "settings.wyckoff.site_pg", BLOB_NAME_SIZE);
    if (!*handled) return 1;
    if (!W.ready && !W.failed) {
        W.failed = build();
        W.ready = !W.failed;
    }
    if (!W.ready) return 1;
    view->entry = &W.entry;
    view->data = (const unsigned char *)W.data;
    return 0;
}
