/* The OPD labels of the isotropy rows, labels of image directions (images.image_directions.label):
 * - isotropy.ordinary.label: the row's image direction (isotropy.ordinary.image_direction);
 * - isotropy.magnetic.label: the k-th magnetic row of an irrep takes the k-th direction of the image with time reversal
 *   added (images.images.add_r of the irrep's image), in that image's list (images.image_subgroups.subgroup). */

#include <stdlib.h>
#include <string.h>

#include "tables.h"

enum { NAMES = 2 };
static const char *const names[NAMES] = {"isotropy.ordinary.label", "isotropy.magnetic.label"};

static struct {
    int ready[NAMES], failed[NAMES];
    unsigned char *data[NAMES];
    blob_section entry[NAMES];
} L;

/* the image direction of each magnetic row; returns the row count, -1 when the tables disagree */
static int64_t magnetic_directions(int64_t **out)
{
    section_view start, image, add_r, sub_start, sub_count, subgroup;
    if (section_find("isotropy.magnetic.by_old_start", &start) || section_find("irreps.old.image", &image) ||
        section_find("images.images.add_r", &add_r) || section_find("images.images.subgroup_start", &sub_start) ||
        section_find("images.images.subgroup_count", &sub_count) || section_find("images.image_subgroups.subgroup", &subgroup))
        return -1;
    int64_t olds = section_count(&start) - 1, rows = section_int(&start, olds) - 1;
    int64_t *dirs = malloc(sizeof(int64_t) * (size_t)(rows + 1));
    if (!dirs) return -1;
    for (int64_t x = 0; x < olds; x++) {
        int64_t first = section_int(&start, x) - 1, n = section_int(&start, x + 1) - 1 - first;
        if (!n) continue;
        int64_t i = x < section_count(&image) ? section_int(&image, x) : 0;
        int64_t j = i >= 1 && i <= section_count(&add_r) ? section_int(&add_r, i - 1) : 0;
        if (j < 1 || j > section_count(&sub_count) || section_int(&sub_count, j - 1) != n) {
            free(dirs);
            return -1;
        }
        for (int64_t k = 0; k < n; k++) dirs[first + k] = section_int(&subgroup, section_int(&sub_start, j - 1) - 1 + k);
    }
    *out = dirs;
    return rows;
}

static int build(int which)
{
    section_view direction, label;
    int64_t *dirs = 0, rows;
    if (section_find("images.image_directions.label", &label)) return 1;
    if (which == 0) {
        if (section_find("isotropy.ordinary.image_direction", &direction)) return 1;
        rows = section_count(&direction);
        if (!(dirs = malloc(sizeof(int64_t) * (size_t)(rows + 1)))) return 1;
        for (int64_t r = 0; r < rows; r++) dirs[r] = section_int(&direction, r);
    } else if ((rows = magnetic_directions(&dirs)) < 0)
        return 1;
    int64_t bytes = 0;
    for (int64_t r = 0; r < rows; r++) {
        const char *text;
        int64_t length;
        if (section_string(&label, dirs[r] - 1, &text, &length)) {
            free(dirs);
            return 1;
        }
        bytes += length;
    }
    L.data[which] = malloc(sizeof(int64_t) * (size_t)(rows + 1) + (size_t)bytes + 1);
    if (!L.data[which]) {
        free(dirs);
        return 1;
    }
    int64_t *offsets = (int64_t *)L.data[which];
    char *pool = (char *)(offsets + rows + 1);
    int64_t at = 0;
    for (int64_t r = 0; r < rows; r++) {
        const char *text;
        int64_t length;
        section_string(&label, dirs[r] - 1, &text, &length);
        offsets[r] = at;
        memcpy(pool + at, text, (size_t)length);
        at += length;
    }
    offsets[rows] = at;
    free(dirs);
    strncpy(L.entry[which].name, names[which], BLOB_NAME_SIZE - 1);
    L.entry[which].kind = BLOB_STR;
    L.entry[which].count = (uint64_t)rows;
    L.entry[which].bytes = sizeof(int64_t) * (uint64_t)(rows + 1) + (uint64_t)at;
    return 0;
}

/* section_find for the labels above; *handled 0 for any other name */
int isotropy_derived_section(const char *name, section_view *view, int *handled)
{
    int which = -1;
    for (int i = 0; i < NAMES && which < 0; i++)
        if (!strncmp(name, names[i], BLOB_NAME_SIZE)) which = i;
    *handled = which >= 0;
    if (!*handled) return 1;
    if (!L.ready[which] && !L.failed[which]) {
        L.failed[which] = build(which);
        L.ready[which] = !L.failed[which];
    }
    if (!L.ready[which]) return 1;
    view->entry = &L.entry[which];
    view->data = L.data[which];
    return 0;
}
