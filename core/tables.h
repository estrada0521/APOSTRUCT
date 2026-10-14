/* Typed views of the mmapped Source blob (tables.c). */

#ifndef APOSTRUCT_TABLES_H
#define APOSTRUCT_TABLES_H

#include <stdint.h>

#include "source_blob.h"

typedef struct {
    const blob_section *entry;
    const unsigned char *data;
} section_view;

int section_find(const char *name, section_view *view);
int section_is_int(const section_view *view);
int64_t section_int(const section_view *view, int64_t index);
int64_t section_count(const section_view *view);
/* The index-th string of a string section without its trailing blanks. */
int section_string(const section_view *view, int64_t index, const char **text, int64_t *length);

#endif
