/* Write the BYU-format tables (data_<table>.txt) that apo reads, from the compiled Source: every stored section and
 * every section apo builds at startup from the others (the star of k, the little-irrep dimensions, types, numbers,
 * old ids, matrix pointers and seed positions).  The sections apo computes in its kernels (group elements, products,
 * isotropy subgroups, ...) are not written here.  With --full-matrices, data_little_full.txt also holds the little
 * irrep matrices at every position (derived from the seeds) in the BYU layout: sparse (row, column, code) triples by
 * column, then row, with BYU's position and irrep pointers.
 *
 *   byu_export <source.blob> <output dir> [--full-matrices] */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "irreps.h"
#include "tables.h"

int tables_open(const char *path);
int little_codes(int gid, int dim, int position, int *codes);
int little_full_dim(int gid);
int irrep_context_k(int gid, int param_count, const int64_t *kparam, irrep_context *out);

static const char *const derived[] = {
    "little.little_k_star", "little.little_k_star_special", "little.little_k_star_ml_pointer", "little.little_k_star_count",
    "little.little_k_star_minusk", "little.little_k_star_conv2ml_pointer", "little.little_irr_full_dim", "little.little_irr_old",
    "little.little_irr_full_matrices_pointer", "little.little_irr_full_matrices_position", "little.little_irr_full_matrices_irr_pointer",
    "little.little_irr_irrnum", "little.little_irr_type",
};

static void write_section(FILE *out, const char *name)
{
    section_view v;
    if (section_find(name, &v)) {
        fprintf(stderr, "byu_export: no section %s\n", name);
        exit(1);
    }
    fprintf(out, "%s\n", strchr(name, '.') + 1);
    int64_t n = section_count(&v);
    for (int64_t i = 0; i < n; i++) {
        if (v.entry->kind == BLOB_STR) {
            const char *text;
            int64_t length;
            section_string(&v, i, &text, &length);
            fprintf(out, "\"%.*s\"", (int)length, text);
        } else if (v.entry->kind == BLOB_BOOL) {
            fputs(section_int(&v, i) ? "T" : "F", out);
        } else if (section_is_int(&v)) {
            fprintf(out, "%lld", (long long)section_int(&v, i));
        } else {
            fprintf(out, "%.17g", ((const double *)v.data)[i]);
        }
        fputc(i % 16 == 15 || i == n - 1 ? '\n' : ' ', out);
    }
}

static void put_list(FILE *out, const char *name, const int64_t *v, int64_t n)
{
    fprintf(out, "%s\n", name);
    for (int64_t i = 0; i < n; i++) fprintf(out, "%lld%c", (long long)v[i], i % 16 == 15 || i == n - 1 ? '\n' : ' ');
}

/* the matrices at every position of every little irrep, in the BYU layout */
static int full_matrices(const char *dir)
{
    section_view sg;
    if (section_find("little.little_irr_space_group", &sg)) return 1;
    int64_t count = section_count(&sg), last = 0;
    for (int64_t g = 0; g < count; g++)
        if (section_int(&sg, g) > 0) last = g + 1;
    int64_t *irr = calloc((size_t)last + 2, sizeof(int64_t)), *pointer = malloc(sizeof(int64_t)), *triples = malloc(sizeof(int64_t));
    int64_t positions = 0, entries = 0, pcap = 1, tcap = 1;
    int *codes = malloc(sizeof(int) * 48 * 48);
    for (int64_t g = 1; g <= last; g++) {
        irr[g - 1] = positions + 1;
        int dim = little_full_dim((int)g);
        irrep_context c;
        const int64_t none[1] = {0};
        if (section_int(&sg, g - 1) <= 0 || dim <= 0 || irrep_context_k((int)g, 0, none, &c)) continue;
        for (int p = 1; p <= c.element_count; p++) {
            if (little_codes((int)g, dim, p, codes)) {
                fprintf(stderr, "byu_export: no matrix for gid %lld position %d\n", (long long)g, p);
                return 1;
            }
            if (positions + 1 >= pcap) pointer = realloc(pointer, sizeof(int64_t) * (size_t)(pcap *= 2));
            pointer[positions++] = entries + 1;
            for (int col = 0; col < dim; col++)
                for (int row = 0; row < dim; row++) {
                    int code = codes[row * dim + col];
                    if (code == 1) continue;
                    if (3 * (entries + 1) >= tcap) triples = realloc(triples, sizeof(int64_t) * (size_t)(tcap = 2 * tcap + 3));
                    triples[3 * entries] = row + 1, triples[3 * entries + 1] = col + 1, triples[3 * entries + 2] = code;
                    entries++;
                }
        }
    }
    irr[last] = positions + 1;
    pointer[positions] = entries + 1;
    char path[4096];
    snprintf(path, sizeof path, "%s/data_little_full.txt", dir);
    FILE *out = fopen(path, "w");
    if (!out) return 1;
    put_list(out, "little_irr_full_matrices_irr_pointer", irr, last + 1);
    put_list(out, "little_irr_full_matrices_pointer", pointer, positions + 1);
    put_list(out, "little_irr_full_matrices", triples, 3 * entries);
    fclose(out);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 3 || tables_open(argv[1])) {
        fprintf(stderr, "usage: byu_export <source.blob> <output dir> [--full-matrices]\n");
        return 2;
    }
    if (argc > 3 && !strcmp(argv[3], "--full-matrices") && full_matrices(argv[2])) return 1;
    static const char *const tables[] = {"space", "irreps", "little", "isotropy", "wyckoff", "images", "magnetic"};
    for (size_t t = 0; t < sizeof tables / sizeof tables[0]; t++) {
        char path[4096], prefix[32];
        snprintf(path, sizeof path, "%s/data_%s.txt", argv[2], tables[t]);
        snprintf(prefix, sizeof prefix, "%s.", tables[t]);
        FILE *out = fopen(path, "w");
        if (!out) return 1;
        const char *name;
        for (uint32_t i = 0; (name = section_stored_name(i)); i++)
            if (!strncmp(name, prefix, strlen(prefix))) write_section(out, name);
        for (size_t d = 0; d < sizeof derived / sizeof derived[0]; d++)
            if (!strncmp(derived[d], prefix, strlen(prefix))) write_section(out, derived[d]);
        fclose(out);
    }
    return 0;
}
