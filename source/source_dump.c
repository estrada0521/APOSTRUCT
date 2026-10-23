/* What apo reads from the compiled Source, independent of how the tables encode it: every stored section and every
 * section apo builds at startup, read through section_find, with constant references written as their values, and
 * the little irrep matrices at every position.  Two Source layouts with the same content give the same dump.
 *
 *   source_dump <source.blob> <output file> */

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
    "k.star.arms", "k.star.special_arms", "k.star.arm_start", "k.star.arm_count", "k.star.minus_k", "k.star.conv2ml_start",
    "irreps.little.dim", "irreps.little.old", "irreps.little.matrix_start", "irreps.little.positions",
    "irreps.little.position_start", "irreps.little.number", "irreps.little.type", "isotropy.ordinary.label",
    "isotropy.magnetic.label", "group.magnetic_ops.code", "group.point_groups.code", "group.point_groups.order",
    "group.point_ops.code", "settings.inter_settings.wyckoff_start", "settings.inter_settings.wyckoff",
};

/* sections whose entries are constant references: every `stride` entries, the ones flagged in `mask` */
static const struct {
    const char *name;
    int stride, mask;
} references[] = {
    {"irreps.little.matrices", 3, 4},
    {"irreps.little.real_form", 3, 6},
    {"images.image_directions.vector", 1, 1},
    {"isotropy.ordinary.direction_codes", 1, 1},
    {"isotropy.magnetic.direction_codes", 1, 1},
    {"settings.site_settings.vector_basis", 1, 1},
};

/* a constant reference as its value: exact as [r, s2, s3, s6] (r + s2 sqrt2 + s3 sqrt3 + s6 sqrt6), 0 for zero, else
 * ~decimal; - for no reference */
static void put_constant(FILE *out, int64_t code)
{
    if (code == 0) {
        fputs("-", out);
        return;
    }
    int64_t w[8];
    if (const_ref_exact(code, w)) {
        double v = const_ref_value(code);
        if (v != v) fprintf(out, "?%lld", (long long)code);
        else fprintf(out, "~%.15g", v);
        return;
    }
    if (!w[0] && !w[2] && !w[4] && !w[6]) {
        fputc('0', out);
        return;
    }
    fputc('[', out);
    for (int j = 0; j < 4; j++) {
        int64_t n = w[2 * j], d = w[2 * j + 1];
        if (j) fputc(',', out);
        if (!n || !d) fputc('0', out);
        else if (d == 1) fprintf(out, "%lld", (long long)n);
        else fprintf(out, "%lld/%lld", (long long)n, (long long)d);
    }
    fputc(']', out);
}

static void put_section(FILE *out, const char *name)
{
    section_view v;
    if (section_find(name, &v)) {
        fprintf(stderr, "source_dump: no section %s\n", name);
        exit(1);
    }
    int stride = 0, mask = 0;
    for (size_t r = 0; r < sizeof references / sizeof references[0]; r++)
        if (!strcmp(references[r].name, name)) stride = references[r].stride, mask = references[r].mask;
    int64_t n = section_count(&v);
    fprintf(out, "# %s %lld\n", name, (long long)n);
    for (int64_t i = 0; i < n; i++) {
        if (v.entry->kind == BLOB_STR) {
            const char *text;
            int64_t length;
            section_string(&v, i, &text, &length);
            fprintf(out, "\"%.*s\"", (int)length, text);
        } else if (stride && (mask >> (i % stride) & 1)) {
            put_constant(out, section_int(&v, i));
        } else if (section_is_int(&v) || v.entry->kind == BLOB_BOOL) {
            fprintf(out, "%lld", (long long)section_int(&v, i));
        } else {
            fprintf(out, "%.17g", ((const double *)v.data)[i]);
        }
        fputc(i % 16 == 15 || i == n - 1 ? '\n' : ' ', out);
    }
}

static int compare_names(const void *a, const void *b) { return strcmp(*(const char *const *)a, *(const char *const *)b); }

/* the matrices at every position of every little irrep: per position, the nonzero (row, column, value) entries */
static int put_full_matrices(FILE *out)
{
    section_view sg;
    if (section_find("irreps.little.sg", &sg)) return 1;
    int *codes = malloc(sizeof(int) * 48 * 48);
    fputs("# little irrep matrices at every position\n", out);
    for (int64_t g = 1; g <= section_count(&sg); g++) {
        int dim = little_full_dim((int)g);
        irrep_context c;
        const int64_t none[1] = {0};
        if (section_int(&sg, g - 1) <= 0 || dim <= 0 || irrep_context_k((int)g, 0, none, &c)) continue;
        for (int p = 1; p <= c.element_count; p++) {
            if (little_codes((int)g, dim, p, codes)) {
                fprintf(stderr, "source_dump: no matrix for gid %lld position %d\n", (long long)g, p);
                return 1;
            }
            fprintf(out, "%lld %d:", (long long)g, p);
            for (int col = 0; col < dim; col++)
                for (int row = 0; row < dim; row++) {
                    int code = codes[row * dim + col];
                    if (code == 1) continue;
                    fprintf(out, " %d,%d,", row + 1, col + 1);
                    put_constant(out, code);
                }
            fputc('\n', out);
        }
    }
    free(codes);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc != 3 || tables_open(argv[1])) {
        fprintf(stderr, "usage: source_dump <source.blob> <output file>\n");
        return 2;
    }
    FILE *out = fopen(argv[2], "w");
    if (!out) return 1;
    const char *names[4096];
    int count = 0;
    for (uint32_t i = 0; section_stored_name(i); i++) {
        const char *name = section_stored_name(i);
        if (!strncmp(name, "const.", 6) || !strncmp(name, "source.", 7)) continue;
        if (count < 4096) names[count++] = name;
    }
    for (size_t i = 0; i < sizeof derived / sizeof derived[0] && count < 4096; i++) names[count++] = derived[i];
    qsort(names, (size_t)count, sizeof names[0], compare_names);
    for (int i = 0; i < count; i++)
        if (!i || strcmp(names[i], names[i - 1])) put_section(out, names[i]);
    int status = put_full_matrices(out);
    fclose(out);
    return status;
}
