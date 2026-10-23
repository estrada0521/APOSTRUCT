/* The baked Source blob (source_blob.h), mmapped once per process and read
 * by the kernels through typed section views. */

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "platform.h"
#include "tables.h"

static const unsigned char *blob;
static size_t blob_size;
static const blob_section *blob_index;
static uint32_t blob_count;

/* 0 on success; the mapping lives for the rest of the process. */
int tables_open(const char *path)
{
    if (blob) return 0;
    size_t size = 0;
    const unsigned char *map = plat_map_file(path, &size);
    if (!map) return 1;
    const blob_header *header = (const blob_header *)map;
    if (size < sizeof(blob_header) || memcmp(header->magic, BLOB_MAGIC, 8) || header->version != BLOB_VERSION ||
        header->file_size != (uint64_t)size ||
        header->index_offset + (uint64_t)header->section_count * sizeof(blob_section) > header->file_size) {
        plat_unmap_file(map, size);
        return 1;
    }
    blob = map;
    blob_size = size;
    blob_index = (const blob_section *)(blob + header->index_offset);
    blob_count = header->section_count;
    return 0;
}

int star_rule_section(const char *name, section_view *view, int *handled);
int little_derived_section(const char *name, section_view *view, int *handled);
int isotropy_derived_section(const char *name, section_view *view, int *handled);
int settings_derived_section(const char *name, section_view *view, int *handled);
int group_derived_section(const char *name, section_view *view, int *handled);

/* a section as stored in the blob */
int section_find_stored(const char *name, section_view *view)
{
    for (uint32_t i = 0; i < blob_count; i++)
        if (!strncmp(blob_index[i].name, name, BLOB_NAME_SIZE)) {
            if (blob_index[i].offset + blob_index[i].bytes > blob_size) return 1;
            view->entry = &blob_index[i];
            view->data = blob + blob_index[i].offset;
            return 0;
        }
    return 1;
}

/* the name of stored section i (0-based), 0 past the last */
const char *section_stored_name(uint32_t i) { return i < blob_count ? blob_index[i].name : 0; }

/* a section by name; the star of k sections are built from their rule (k/star.c), the irrep dimensions and old ids
 * and the ordinary OPD labels from other tables (tables/little_derived.c, tables/isotropy_derived.c) */
int section_find(const char *name, section_view *view)
{
    int handled, status = star_rule_section(name, view, &handled);
    if (!handled) status = little_derived_section(name, view, &handled);
    if (!handled) status = isotropy_derived_section(name, view, &handled);
    if (!handled) status = settings_derived_section(name, view, &handled);
    if (!handled) status = group_derived_section(name, view, &handled);
    return handled ? status : section_find_stored(name, view);
}

int64_t section_int(const section_view *view, int64_t index)
{
    switch (view->entry->kind) {
    case BLOB_I8: return ((const int8_t *)view->data)[index];
    case BLOB_I16: return ((const int16_t *)view->data)[index];
    case BLOB_I32: return ((const int32_t *)view->data)[index];
    case BLOB_I64: return ((const int64_t *)view->data)[index];
    case BLOB_BOOL: return view->data[index];
    default: return 0;
    }
}

int section_is_int(const section_view *view)
{
    uint32_t kind = view->entry->kind;
    return kind == BLOB_I8 || kind == BLOB_I16 || kind == BLOB_I32 || kind == BLOB_I64 || kind == BLOB_BOOL;
}

int64_t section_count(const section_view *view)
{
    return (int64_t)view->entry->count;
}

int section_string(const section_view *view, int64_t index, const char **text, int64_t *length)
{
    if (view->entry->kind != BLOB_STR || index < 0 || (uint64_t)index >= view->entry->count) return 1;
    const int64_t *offsets = (const int64_t *)view->data;
    const char *pool = (const char *)(offsets + view->entry->count + 1);
    int64_t start = offsets[index], end = offsets[index + 1];
    while (end > start && (pool[end - 1] == ' ' || pool[end - 1] == '\t')) end--;
    while (start < end && (pool[start] == ' ' || pool[start] == '\t')) start++;
    *text = pool + start;
    *length = end - start;
    return 0;
}

/* One Source little-irrep sparse matrix (data_little little_irr_full_*)
 * into out (dim x dim).  Returns nonzero when the tables do not hold it. */
static section_view const_values, const_present, const_exact;

static int const_ready(void)
{
    static int ready;
    if (!ready && blob && !section_find("const.values", &const_values) && !section_find("const.present", &const_present) &&
        !section_find("const.exact", &const_exact) && const_values.entry->kind == BLOB_F64)
        ready = 1;
    return ready;
}

double const_ref_value(int64_t ref)
{
    int64_t k = ref < 0 ? -ref : ref;
    if (!const_ready() || k <= 0 || k >= (int64_t)const_values.entry->count || !const_present.data[k]) return NAN;
    double v = ((const double *)const_values.data)[k];
    return ref < 0 ? -v : v;
}

int const_ref_exact(int64_t ref, int64_t words[8])
{
    int64_t k = ref < 0 ? -ref : ref;
    if (!const_ready() || k <= 0 || k >= (int64_t)const_values.entry->count || !const_present.data[k]) return 1;
    int has = 0;
    for (int j = 0; j < 8; j++) words[j] = section_int(&const_exact, 8 * k + j);
    for (int j = 0; j < 4; j++) {
        has |= words[2 * j + 1] != 0;
        if (ref < 0) words[2 * j] = -words[2 * j];
    }
    return !has;
}

int little_codes(int gid, int dim, int position, int *codes);

int little_sparse_matrix(int gid, int dim, int position, double *out)
{
    if (dim <= 0 || dim > 48) return 1;
    int codes[48 * 48];
    if (little_codes(gid, dim, position, codes)) return 1;
    for (int i = 0; i < dim * dim; i++) {
        out[i] = const_ref_value(codes[i]);
        if (out[i] != out[i]) return 1;
    }
    return 0;
}

/* Every generated position's sparse matrix of one little irrep; available[p]
 * is 0 where the tables do not hold position p + 1. */
int little_sparse_context(int gid, int dim, int element_count, double *sparse, int *available)
{
    for (int position = 1; position <= element_count; position++) {
        double *out = sparse + (int64_t)(position - 1) * dim * dim;
        if (little_sparse_matrix(gid, dim, position, out)) {
            for (int i = 0; i < dim * dim; i++) out[i] = 0.0;
            available[position - 1] = 0;
        } else {
            available[position - 1] = 1;
        }
    }
    return 0;
}
