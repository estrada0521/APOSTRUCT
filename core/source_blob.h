/* The baked Source blob: one little-endian file the kernels mmap.
 *
 *   header   magic "APSTBLB1", u32 version, u32 section count,
 *            u64 index offset, u64 file size
 *   index    section_count entries of struct blob_section
 *   data     every section's payload at a 64-byte aligned offset
 *
 * Integer sections use the narrowest of int8/16/32/64 holding their values.
 * A string section stores int64 offsets[count + 1] into the UTF-8 bytes that
 * follow them.  "source.sha256" lists "file=hex" for every authority file. */

#ifndef APOSTRUCT_SOURCE_BLOB_H
#define APOSTRUCT_SOURCE_BLOB_H

#include <stdint.h>

#define BLOB_MAGIC "APSTBLB1"
#define BLOB_VERSION 1u
#define BLOB_NAME_SIZE 96

enum blob_kind {
    BLOB_I8 = 1,
    BLOB_I16 = 2,
    BLOB_I32 = 3,
    BLOB_I64 = 4,
    BLOB_F64 = 5,
    BLOB_BOOL = 6,
    BLOB_STR = 7,
};

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t section_count;
    uint64_t index_offset;
    uint64_t file_size;
} blob_header;

typedef struct {
    char name[BLOB_NAME_SIZE];
    uint32_t kind;
    uint32_t reserved;
    uint64_t count;
    uint64_t offset;
    uint64_t bytes;
} blob_section;

#endif
