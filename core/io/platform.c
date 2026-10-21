#include "platform.h"

#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#include <windows.h>

void *plat_map_anon(size_t size)
{
    return VirtualAlloc(0, size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
}

void plat_unmap_anon(void *p, size_t size)
{
    (void)size;
    VirtualFree(p, 0, MEM_RELEASE);
}

const unsigned char *plat_map_file(const char *path, size_t *size)
{
    HANDLE file = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
    if (file == INVALID_HANDLE_VALUE) return 0;
    LARGE_INTEGER length;
    if (!GetFileSizeEx(file, &length) || length.QuadPart <= 0) {
        CloseHandle(file);
        return 0;
    }
    HANDLE mapping = CreateFileMappingA(file, 0, PAGE_READONLY, 0, 0, 0);
    CloseHandle(file);
    if (!mapping) return 0;
    const unsigned char *view = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
    CloseHandle(mapping);
    if (view) *size = (size_t)length.QuadPart;
    return view;
}

void plat_unmap_file(const unsigned char *view, size_t size)
{
    (void)size;
    UnmapViewOfFile(view);
}

char *plat_realpath(const char *path, char *out)
{
    if (!_fullpath(out, path, 4096) || GetFileAttributesA(out) == INVALID_FILE_ATTRIBUTES) return 0;
    for (char *c = out; *c; c++)
        if (*c == '\\') *c = '/';
    return out;
}

int plat_mkdir(const char *path) { return _mkdir(path); }

int plat_exe_path(char *out, size_t size)
{
    DWORD n = GetModuleFileNameA(0, out, (DWORD)size);
    if (!n || n >= size) return 1;
    for (char *c = out; *c; c++)
        if (*c == '\\') *c = '/';
    return 0;
}

#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

void *plat_map_anon(size_t size)
{
    void *p = mmap(0, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    return p == MAP_FAILED ? 0 : p;
}

void plat_unmap_anon(void *p, size_t size) { munmap(p, size); }

const unsigned char *plat_map_file(const char *path, size_t *size)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) return 0;
    struct stat st;
    if (fstat(fd, &st) || st.st_size <= 0) {
        close(fd);
        return 0;
    }
    void *map = mmap(0, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (map == MAP_FAILED) return 0;
    *size = (size_t)st.st_size;
    return map;
}

void plat_unmap_file(const unsigned char *view, size_t size) { munmap((void *)view, size); }

char *plat_realpath(const char *path, char *out) { return realpath(path, out); }

int plat_mkdir(const char *path) { return mkdir(path, 0777); }
#endif

char *plat_strndup(const char *s, size_t n)
{
    size_t length = strnlen(s, n);
    char *out = malloc(length + 1);
    if (!out) return 0;
    memcpy(out, s, length);
    out[length] = 0;
    return out;
}
