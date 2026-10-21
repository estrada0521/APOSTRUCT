/* The print layout of one projected mode block: each mode's site-irrep
 * label (summary, occurrence suffix, component) in the plain and the
 * OPD-row-component conventions, its Source print identity, and the
 * Gram-Schmidt keys and carrier family rows of present_block. */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "tables.h"

enum {
    BL_OK = 0,
    BL_RESOLVE = 1,
    BL_NONPOSITIVE = 2,
    BL_DIVIDE = 3,
    BL_COUNT = 4,
    BL_MEMORY = 5,
    BL_TABLES = 6,
};

enum { MODE_STRIDE = 7, ROW_STRIDE = 19 };

typedef struct {
    int64_t pg_irrep, family, source_rows, groups, flat;
    int64_t row_count;
    int explicit_rows, single_rows;
    int base_class, dimension, display_dimension, occurrences;
} summary;

static struct {
    section_view pg_irrep, irrep_image, image_dimension, label_pg;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (section_find("settings.site_point_groups.irrep", &T.pg_irrep) || section_find("irreps.old.image", &T.irrep_image) ||
        section_find("images.images.dimension", &T.image_dimension) ||
        section_find("irreps.old.mulliken", &T.label_pg))
        return BL_TABLES;
    T.ready = 1;
    return 0;
}

static int resolve(int site_pg, int64_t pg_irrep, int *old_id, int *dimension)
{
    int64_t slot = (int64_t)(site_pg - 1) * 12 + pg_irrep - 1;
    if (slot < 0 || slot >= section_count(&T.pg_irrep)) return BL_RESOLVE;
    int64_t old = section_int(&T.pg_irrep, slot);
    if (old <= 0 || old > section_count(&T.irrep_image) || old > section_count(&T.label_pg)) return BL_RESOLVE;
    int64_t image = section_int(&T.irrep_image, old - 1);
    if (image <= 0 || image > section_count(&T.image_dimension)) return BL_RESOLVE;
    *old_id = (int)old;
    *dimension = (int)section_int(&T.image_dimension, image - 1);
    return 0;
}

static int same_label(int left, int right)
{
    const char *a, *b;
    int64_t la, lb;
    if (section_string(&T.label_pg, left - 1, &a, &la) || section_string(&T.label_pg, right - 1, &b, &lb)) return 0;
    return la == lb && !memcmp(a, b, (size_t)la);
}

static int layout_labels(const summary *s, int count, int little_type, int opd_rows, int mode_count, int *totals,
                         int *seen, summary *local, int64_t *out, int64_t *diagnostic)
{
    for (int i = 0; i < count; i++) totals[i] = seen[i] = 0;
    memcpy(local, s, sizeof *s * (size_t)count);
    for (int i = 0; i < count; i++) {
        summary *x = local + i;
        int64_t flat = x->flat, dimension = x->dimension;
        if (dimension <= 0 || flat <= 0) return BL_NONPOSITIVE;
        int64_t display = dimension;
        int eligible = (little_type == 3 || x->groups > 1) && x->row_count > 0 && flat % x->row_count == 0;
        if (eligible && (opd_rows || (flat % dimension != 0 && x->explicit_rows)))
            display = x->row_count;
        else if (x->source_rows > 1 && flat == x->source_rows && x->groups == x->source_rows)
            display = 1;
        else if (little_type == 3 && x->source_rows > 1 && x->groups != x->source_rows && flat % x->source_rows == 0)
            display = x->source_rows;
        if (flat % display) {
            if (flat != 1) {
                diagnostic[0] = flat;
                diagnostic[1] = display;
                diagnostic[2] = i;
                return BL_DIVIDE;
            }
            display = 1;
        }
        x->display_dimension = (int)display;
        x->occurrences = (int)(flat / display);
        totals[x->base_class] += x->occurrences;
    }
    int written = 0;
    for (int i = 0; i < count; i++) {
        const summary *x = local + i;
        if (!(opd_rows && (little_type == 3 || x->groups > 1) && x->row_count > 0) && x->source_rows > 1 &&
            x->flat == x->source_rows * x->display_dimension && (little_type == 3 || x->groups == x->source_rows) &&
            (x->occurrences == 1 || x->groups == x->source_rows)) {
            for (int64_t c = 0; c < x->flat; c++, written++)
                if (written < mode_count) {
                    out[MODE_STRIDE * written] = i;
                    out[MODE_STRIDE * written + 1] = 0;
                    out[MODE_STRIDE * written + 2] = c;
                }
            continue;
        }
        for (int o = 0; o < x->occurrences; o++) {
            int suffix = ++seen[x->base_class];
            for (int c = 0; c < x->display_dimension; c++, written++)
                if (written < mode_count) {
                    out[MODE_STRIDE * written] = i;
                    out[MODE_STRIDE * written + 1] = totals[x->base_class] > 1 ? suffix : 0;
                    out[MODE_STRIDE * written + 2] = c;
                }
        }
    }
    if (written != mode_count) {
        diagnostic[0] = written;
        diagnostic[1] = mode_count;
        return BL_COUNT;
    }
    return 0;
}

/* summaries: the mode_block_ints summary words (count summaries); direct:
 * the display plan was adopted, so every OPD group is one direction row of
 * full_dim values; otherwise group rows are source rows [start, start +
 * length) of row_total rows with row_lengths.  block: kind code, dynamic,
 * phase, Stokes columns, spec order, site order.  modes (mode_count x 7):
 * summary, suffix, component, OPD-row suffix, OPD-row component, print
 * component, component within its summary; keys (x 3) and rows (x 19) for
 * present_block when the print identity covers every mode.  status: the
 * plain and the OPD-row label statuses; flags: identity covers the modes,
 * OPD-row labels are wanted.  Returns 0 or BL_TABLES / BL_MEMORY. */
int block_layout(int site_pg, int little_type, int count, const int64_t *summaries, int direct, int full_dim,
                 int row_total, const int *row_lengths, int mode_count, const int64_t *block, int64_t *modes,
                 int64_t *keys, int64_t *rows, int *status, int *flags, int64_t *diagnostic)
{
    if (load()) return BL_TABLES;
    summary *s = calloc((size_t)count + 1, sizeof *s);
    int *totals = calloc((size_t)count + 1, sizeof *totals), *seen = calloc((size_t)count + 1, sizeof *seen);
    int *old_ids = calloc((size_t)count + 1, sizeof *old_ids);
    int64_t *display = malloc(sizeof(int64_t) * MODE_STRIDE * ((size_t)mode_count + 1));
    summary *local = malloc(sizeof *local * ((size_t)count + 1));
    if (!s || !totals || !seen || !old_ids || !display || !local) {
        free(local);
        free(s);
        free(totals);
        free(seen);
        free(old_ids);
        free(display);
        return BL_MEMORY;
    }
    const int64_t *at = summaries;
    int want_display = 0, resolved = 1;
    for (int i = 0; i < count; i++) {
        summary *x = s + i;
        x->pg_irrep = at[0];
        x->family = at[1];
        x->source_rows = at[2];
        x->groups = at[3];
        const int64_t *groups = at + 4;
        x->explicit_rows = x->groups > 0;
        x->single_rows = 1;
        for (int64_t g = 0; g < x->groups; g++) {
            int64_t length = 1;
            int nonempty = full_dim > 0;
            if (!direct) {
                int64_t start = groups[3 * g], stop = start + groups[3 * g + 1];
                if (stop > row_total) stop = row_total;
                length = stop > start ? stop - start : 0;
                nonempty = 1;
                for (int64_t r = start; r < stop; r++) nonempty &= row_lengths[r] > 0;
            }
            x->row_count += length;
            x->explicit_rows &= length > 0 && nonempty;
            x->single_rows &= length == 1;
        }
        at += 4 + 4 * x->groups;
        x->flat = at[0];
        at += 5;
        if (little_type == 3 || x->groups > 1) want_display = 1;
        if (resolved && resolve(site_pg, x->pg_irrep, old_ids + i, &x->dimension)) resolved = 0;
    }
    int64_t identity_count = 0;
    for (int i = 0; i < count; i++) identity_count += s[i].flat;
    int covered = identity_count == mode_count;
    flags[0] = covered;
    flags[1] = want_display;
    status[0] = status[1] = 0;
    if (!resolved) status[0] = status[1] = BL_RESOLVE;
    else {
        for (int i = 0; i < count; i++) {
            s[i].base_class = i;
            for (int j = 0; j < i; j++)
                if (same_label(old_ids[i], old_ids[j])) {
                    s[i].base_class = s[j].base_class;
                    break;
                }
        }
        status[0] = layout_labels(s, count, little_type, 0, mode_count, totals, seen, local, modes, diagnostic);
        if (want_display)
            status[1] = layout_labels(s, count, little_type, 1, mode_count, totals, seen, local, display, diagnostic + 3);
        for (int m = 0; m < mode_count; m++) {
            modes[MODE_STRIDE * m + 3] = want_display && !status[1] ? display[MODE_STRIDE * m + 1] : 0;
            modes[MODE_STRIDE * m + 4] = want_display && !status[1] ? display[MODE_STRIDE * m + 2] : 0;
        }
    }
    if (covered) {
        int m = 0;
        for (int i = 0; i < count; i++) {
            const summary *x = s + i;
            int64_t blocks = x->row_count > 0 && x->flat % x->row_count == 0 ? x->flat / x->row_count : 0;
            int projected = little_type == 3 && blocks > 0;
            const int64_t *groups = summaries;
            if (projected) {
                const int64_t *walk = summaries;
                for (int j = 0; j < i; j++) walk += 4 + 4 * walk[3] + 5;
                groups = walk + 4;
            }
            int64_t row_offset = 0, written = 0;
            for (int64_t g = 0; projected && g < x->groups; g++) {
                int64_t length = 1;
                if (!direct) {
                    int64_t start = groups[3 * g], stop = start + groups[3 * g + 1];
                    if (stop > row_total) stop = row_total;
                    length = stop > start ? stop - start : 0;
                }
                for (int64_t b = 0; b < blocks; b++)
                    for (int64_t r = 0; r < length; r++, written++)
                        if (m + written < mode_count)
                            modes[MODE_STRIDE * (m + written) + 5] = b * x->row_count + row_offset + r;
                row_offset += length;
            }
            int keep = projected && written == x->flat;
            for (int64_t c = 0; c < x->flat; c++, m++) {
                if (!keep) modes[MODE_STRIDE * m + 5] = c;
                modes[MODE_STRIDE * m + 6] = c;
            }
        }
        int *key_ids = malloc(sizeof(int) * ((size_t)mode_count + 1)), *row_ids = malloc(sizeof(int) * ((size_t)mode_count + 1));
        if (!key_ids || !row_ids) {
            free(key_ids);
            free(row_ids);
            free(local);
            free(s);
            free(totals);
            free(seen);
            free(old_ids);
            free(display);
            return BL_MEMORY;
        }
        for (int i = 0; i <= mode_count; i++) key_ids[i] = row_ids[i] = -1;
        int key_count = 0, row_count = 0;
        int dynamic = (int)block[1], phase = (int)block[2], stokes = (int)block[3];
        for (m = 0; m < mode_count && !status[0]; m++) {
            const summary *x = s + modes[MODE_STRIDE * m];
            int64_t component = modes[MODE_STRIDE * m + 2];
            int blank = !stokes && (little_type == 1 || little_type == 3) && x->source_rows > 1 &&
                        (dynamic || (little_type == 3 && x->groups > 1));
            if (component < 0 || component > mode_count) component = mode_count;
            if (!blank && key_ids[component] < 0) key_ids[component] = key_count++;
            if (row_ids[component] < 0) row_ids[component] = row_count++;
            int64_t key_id = blank ? -1 : key_ids[component], row_id = row_ids[component];
            keys[3 * m] = site_pg;
            keys[3 * m + 1] = x->pg_irrep;
            keys[3 * m + 2] = key_id;
            int64_t *row = rows + ROW_STRIDE * m, print_component = modes[MODE_STRIDE * m + 5];
            row[0] = block[0];
            row[1] = 1;
            row[2] = x->family;
            row[3] = 1;
            row[4] = print_component;
            row[5] = x->flat;
            row[6] = dynamic;
            row[7] = phase;
            row[8] = little_type;
            row[9] = x->source_rows;
            row[10] = x->groups;
            row[11] = x->groups ? x->row_count : 0;
            row[12] = x->single_rows;
            row[13] = modes[MODE_STRIDE * m + 6] == print_component;
            row[14] = block[4];
            row[15] = block[5];
            row[16] = site_pg;
            row[17] = x->pg_irrep;
            row[18] = row_id;
        }
        free(key_ids);
        free(row_ids);
    }
    free(local);
    free(s);
    free(totals);
    free(seen);
    free(old_ids);
    free(display);
    return 0;
}
