/* isotropy.ordinary.label: each ordinary isotropy row's OPD label is the label of its image direction
 * (isotropy.ordinary.image_direction into images.image_directions.label). */

#include <stdlib.h>
#include <string.h>

#include "tables.h"

static struct {
    int ready, failed;
    unsigned char *data;
    blob_section entry;
} L;

static int build(void)
{
    section_view direction, label;
    if (section_find("isotropy.ordinary.image_direction", &direction) || section_find("images.image_directions.label", &label))
        return 1;
    int64_t rows = section_count(&direction), bytes = 0;
    for (int64_t r = 0; r < rows; r++) {
        const char *text;
        int64_t length;
        if (section_string(&label, section_int(&direction, r) - 1, &text, &length)) return 1;
        bytes += length;
    }
    L.data = malloc(sizeof(int64_t) * (size_t)(rows + 1) + (size_t)bytes + 1);
    if (!L.data) return 1;
    int64_t *offsets = (int64_t *)L.data;
    char *pool = (char *)(offsets + rows + 1);
    int64_t at = 0;
    for (int64_t r = 0; r < rows; r++) {
        const char *text;
        int64_t length;
        section_string(&label, section_int(&direction, r) - 1, &text, &length);
        offsets[r] = at;
        memcpy(pool + at, text, (size_t)length);
        at += length;
    }
    offsets[rows] = at;
    strncpy(L.entry.name, "isotropy.ordinary.label", BLOB_NAME_SIZE - 1);
    L.entry.kind = BLOB_STR;
    L.entry.count = (uint64_t)rows;
    L.entry.bytes = sizeof(int64_t) * (uint64_t)(rows + 1) + (uint64_t)at;
    return 0;
}

/* section_find for isotropy.ordinary.label; *handled 0 for any other name */
int isotropy_derived_section(const char *name, section_view *view, int *handled)
{
    *handled = !strncmp(name, "isotropy.ordinary.label", BLOB_NAME_SIZE);
    if (!*handled) return 1;
    if (!L.ready && !L.failed) {
        L.failed = build();
        L.ready = !L.failed;
    }
    if (!L.ready) return 1;
    view->entry = &L.entry;
    view->data = L.data;
    return 0;
}
