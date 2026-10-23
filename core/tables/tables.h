/* Typed views of the mmapped Source blob (tables/blob.c). */

#ifndef APOSTRUCT_TABLES_H
#define APOSTRUCT_TABLES_H

#include <stdint.h>

#include "source_blob.h"

typedef struct {
    const blob_section *entry;
    const unsigned char *data;
} section_view;

int section_find(const char *name, section_view *view);
const char *section_stored_name(uint32_t i);
int section_is_int(const section_view *view);
int64_t section_int(const section_view *view, int64_t index);
int64_t section_count(const section_view *view);
/* The index-th string of a string section without its trailing blanks. */
int section_string(const section_view *view, int64_t index, const char **text, int64_t *length);
/* A constant reference (constants.txt): k is constant k, -k its negative, 0 no reference. The value, NAN when there
 * is none; the exact form as four (numerator, denominator) pairs of r + s2 sqrt2 + s3 sqrt3 + s6 sqrt6, nonzero when
 * the constant has none. */
double const_ref_value(int64_t ref);
int const_ref_exact(int64_t ref, int64_t words[8]);

#endif
