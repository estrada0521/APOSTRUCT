/* The few operating-system services apo uses, for POSIX (macOS, Linux) and Windows (platform.c): anonymous memory
 * for large blocks, a read-only file mapping, canonical paths and directories. */

#ifndef APOSTRUCT_PLATFORM_H
#define APOSTRUCT_PLATFORM_H

#include <stddef.h>

void *plat_map_anon(size_t size);
void plat_unmap_anon(void *p, size_t size);
/* the read-only contents of a file for the rest of the process, 0 on failure */
const unsigned char *plat_map_file(const char *path, size_t *size);
void plat_unmap_file(const unsigned char *view, size_t size);
/* the absolute path of an existing path into out (PATH_MAX bytes), with '/' separators; 0 on failure */
char *plat_realpath(const char *path, char *out);
int plat_mkdir(const char *path);
char *plat_strndup(const char *s, size_t n);
#ifdef _WIN32
/* the path of the running executable, with '/' separators; 1 on failure */
int plat_exe_path(char *out, size_t size);
#endif

#endif
