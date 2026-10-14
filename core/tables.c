/* The baked Source blob (source_blob.h), mmapped once per process and read
 * by the kernels through typed section views. */

#include <fcntl.h>
#include <stdint.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "tables.h"

static const unsigned char *blob;
static size_t blob_size;
static const blob_section *blob_index;
static uint32_t blob_count;

/* 0 on success; the mapping lives for the rest of the process. */
int tables_open(const char *path)
{
    if (blob) return 0;
    int fd = open(path, O_RDONLY);
    if (fd < 0) return 1;
    struct stat st;
    if (fstat(fd, &st) || st.st_size < (off_t)sizeof(blob_header)) {
        close(fd);
        return 1;
    }
    void *map = mmap(0, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (map == MAP_FAILED) return 1;
    const blob_header *header = map;
    if (memcmp(header->magic, BLOB_MAGIC, 8) || header->version != BLOB_VERSION
        || header->file_size != (uint64_t)st.st_size
        || header->index_offset + (uint64_t)header->section_count * sizeof(blob_section) > header->file_size) {
        munmap(map, (size_t)st.st_size);
        return 1;
    }
    blob = map;
    blob_size = (size_t)st.st_size;
    blob_index = (const blob_section *)(blob + header->index_offset);
    blob_count = header->section_count;
    return 0;
}

int section_find(const char *name, section_view *view)
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
int little_sparse_matrix(int gid, int dim, int position, double *out)
{
    static section_view irr_pointer, pointer, matrices, values, present;
    static int ready;
    if (!ready) {
        if (!blob || section_find("little.little_irr_full_matrices_irr_pointer", &irr_pointer)
            || section_find("little.little_irr_full_matrices_pointer", &pointer)
            || section_find("little.little_irr_full_matrices", &matrices)
            || section_find("const.values", &values) || section_find("const.present", &present)
            || !section_is_int(&irr_pointer) || !section_is_int(&pointer) || !section_is_int(&matrices)
            || values.entry->kind != BLOB_F64)
            return 1;
        ready = 1;
    }
    int64_t irr_length = (int64_t)irr_pointer.entry->count, pointer_length = (int64_t)pointer.entry->count;
    int64_t matrices_length = (int64_t)matrices.entry->count, constant_count = (int64_t)values.entry->count;
    if (gid < 1 || gid > irr_length || position < 1 || dim <= 0) return 1;
    int64_t matrix_index = section_int(&irr_pointer, gid - 1) + position - 1;
    if (matrix_index < 1 || matrix_index >= pointer_length) return 1;
    int64_t start = section_int(&pointer, matrix_index - 1), end = section_int(&pointer, matrix_index);
    if (start < 1 || end < start || (end - 1) * 3 > matrices_length) return 1;
    const double *constants = (const double *)values.data;
    for (int i = 0; i < dim * dim; i++) out[i] = 0.0;
    for (int64_t offset = (start - 1) * 3; offset < (end - 1) * 3; offset += 3) {
        int64_t row = section_int(&matrices, offset), col = section_int(&matrices, offset + 1);
        int64_t code = section_int(&matrices, offset + 2);
        if (row < 1 || row > dim || col < 1 || col > dim) return 1;
        if (code < 0 || code >= constant_count || !present.data[code]) return 1;
        out[(row - 1) * dim + col - 1] = constants[code];
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
