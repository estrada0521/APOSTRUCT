/* structure_runtime's mode topology and row plan of one site block: every
 * Source row/centering representation assigned to its canonical child atom,
 * the raw vector slot each atom reads, its Source-equivalent slots, the
 * presentation inverse and the magnetic occurrence negation. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "tables.h"

int row_assignments(int child_sg, const int64_t *basis, int point_operation_count, int atom_count, const int64_t *atoms,
                    int row_count, const int64_t *rows, int *out);

enum {
    BP_TOPOLOGY = 1,     /* mode topology requires Source rows and centerings */
    BP_COVER = 2,        /* mode topology does not cover every canonical child atom */
    BP_SINGULAR = 3,     /* mode presentation basis is singular */
    BP_LAYOUT = 4,       /* mode rows do not cover the canonical presentation layout */
    BP_BIJECTION = 5,    /* mode rows do not bijectively cover canonical child atoms */
    BP_MEMORY = 6,
    BP_TABLES = 7,
    BP_ASSIGNMENT = 100, /* + row_assignments status */
};

static struct {
    int count, centering, vector_count, has_inverse;
    int64_t *raw_indices, *offsets, *equivalents;
    unsigned char *negated;
    double inverse[9];
} P;

static int float_inverse3(const double *m, double *out)
{
    double det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    if (fabs(det) < 1e-12) return 1;
    out[0] = (m[4] * m[8] - m[5] * m[7]) / det;
    out[1] = (m[2] * m[7] - m[1] * m[8]) / det;
    out[2] = (m[1] * m[5] - m[2] * m[4]) / det;
    out[3] = (m[5] * m[6] - m[3] * m[8]) / det;
    out[4] = (m[0] * m[8] - m[2] * m[6]) / det;
    out[5] = (m[2] * m[3] - m[0] * m[5]) / det;
    out[6] = (m[3] * m[7] - m[4] * m[6]) / det;
    out[7] = (m[1] * m[6] - m[0] * m[7]) / det;
    out[8] = (m[0] * m[4] - m[1] * m[3]) / det;
    return 0;
}

static int negative_determinant(int sg, int op, int *out)
{
    int64_t m[9];
    geo_pml_point_matrix(sg, op, m);
    int64_t det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    *out = det < 0;
    return 0;
}

static void reset(void)
{
    free(P.raw_indices);
    free(P.offsets);
    free(P.equivalents);
    free(P.negated);
    memset(&P, 0, sizeof P);
}

/* raw_vector_count < 0 checks the topology only.  presented: whether the block's raw vectors are presented through
 * presentation_basis (then every row/centering is read), else they already
 * follow the layout's presentation rows.  Returns 0 or a BP_ status. */
int block_plan(int sg, int child_sg, const int64_t *basis_pml, int point_operation_count, int atom_count, const int64_t *topology_atoms,
               int presentation_row_count, int layout_point_operation_count, int record_count, const int64_t *records,
               int raw_vector_count, int presented, const double *presentation_basis, int magnetic)
{
    reset();
    if (geo_load()) return BP_TABLES;
    rat centering[4][3];
    int cc = geo_centering(child_sg, centering);
    if (cc <= 0 || record_count <= 0) return BP_TOPOLOGY;
    int64_t row_count = (int64_t)record_count * cc;
    int64_t *rows = malloc(sizeof(int64_t) * 6 * (size_t)row_count);
    int *owners = malloc(sizeof(int) * (size_t)row_count);
    if (!rows || !owners) {
        free(rows);
        free(owners);
        return BP_MEMORY;
    }
    for (int r = 0; r < record_count; r++)
        for (int c = 0; c < cc; c++) {
            int64_t *row = rows + 6 * ((int64_t)r * cc + c);
            memcpy(row, records + 5 * r, 5 * sizeof(int64_t));
            row[5] = c;
        }
    int status = row_assignments(child_sg, basis_pml, point_operation_count, atom_count, topology_atoms, (int)row_count, rows, owners);
    free(rows);
    if (status) {
        free(owners);
        return BP_ASSIGNMENT + status;
    }
    int *counts = calloc((size_t)atom_count + 1, sizeof(int));
    if (!counts) {
        free(owners);
        return BP_MEMORY;
    }
    int distinct = 0;
    for (int64_t i = 0; i < row_count; i++) {
        if (owners[i] < 0 || owners[i] >= atom_count) {
            free(owners);
            free(counts);
            return BP_COVER;
        }
        distinct += !counts[owners[i]]++;
    }
    if (distinct != atom_count) {
        free(owners);
        free(counts);
        return BP_COVER;
    }
    if (raw_vector_count < 0) {
        free(owners);
        free(counts);
        return 0;
    }
    P.count = atom_count;
    P.centering = cc;
    P.raw_indices = malloc(sizeof(int64_t) * (size_t)(atom_count + 1));
    P.offsets = calloc((size_t)atom_count + 1, sizeof(int64_t));
    P.equivalents = malloc(sizeof(int64_t) * (size_t)(row_count + 1));
    P.negated = calloc((size_t)atom_count + 1, 1);
    if (!P.raw_indices || !P.offsets || !P.equivalents || !P.negated) {
        free(owners);
        free(counts);
        return BP_MEMORY;
    }
    if (presented) {
        if (float_inverse3(presentation_basis, P.inverse)) {
            free(owners);
            free(counts);
            return BP_SINGULAR;
        }
        P.has_inverse = 1;
        int64_t *starts = malloc(sizeof(int64_t) * (size_t)(atom_count + 1)), *fill = malloc(sizeof(int64_t) * (size_t)(atom_count + 1));
        int64_t *by_atom = malloc(sizeof(int64_t) * (size_t)(row_count + 1));
        if (!starts || !fill || !by_atom) {
            free(owners);
            free(counts);
            free(starts);
            free(fill);
            free(by_atom);
            return BP_MEMORY;
        }
        starts[0] = 0;
        for (int a = 0; a < atom_count; a++) starts[a + 1] = starts[a] + counts[a];
        memcpy(fill, starts, sizeof(int64_t) * (size_t)atom_count);
        for (int64_t i = 0; i < row_count; i++) by_atom[fill[owners[i]]++] = i / cc;
        int64_t equivalents = 0;
        for (int a = 0; a < atom_count; a++) {
            P.raw_indices[a] = by_atom[starts[a]];
            for (int64_t k = starts[a] + 1; k < starts[a + 1]; k++) P.equivalents[equivalents++] = by_atom[k];
            P.offsets[a + 1] = equivalents;
        }
        free(starts);
        free(fill);
        free(by_atom);
        P.vector_count = record_count;
    } else {
        if (layout_point_operation_count <= 0 || presentation_row_count != atom_count || raw_vector_count != presentation_row_count ||
            raw_vector_count > record_count) {
            free(owners);
            free(counts);
            return BP_LAYOUT;
        }
        unsigned char *seen = calloc((size_t)atom_count + 1, 1);
        if (!seen) {
            free(owners);
            free(counts);
            return BP_MEMORY;
        }
        int unique = 0;
        for (int r = 0; r < raw_vector_count; r++) unique += !seen[owners[(int64_t)r * cc]]++;
        free(seen);
        if (unique != raw_vector_count) {
            free(owners);
            free(counts);
            return BP_BIJECTION;
        }
        for (int r = 0; r < raw_vector_count; r++) P.raw_indices[owners[(int64_t)r * cc]] = r;
        P.vector_count = raw_vector_count;
    }
    free(owners);
    free(counts);
    if (magnetic)
        for (int a = 0; a < atom_count; a++) {
            int negative;
            if (negative_determinant(sg, (int)records[5 * P.raw_indices[a] + 4], &negative)) return BP_TABLES;
            P.negated[a] = (unsigned char)negative;
        }
    return 0;
}

/* The last plan: vector count, has inverse, inverse[9] and its arrays. */
int block_plan_result(int64_t *header, double *inverse, const int64_t **raw_indices, const int64_t **offsets,
                      const int64_t **equivalents, const unsigned char **negated)
{
    header[0] = P.vector_count;
    header[1] = P.has_inverse;
    header[2] = P.count ? P.offsets[P.count] : 0;
    memcpy(inverse, P.inverse, sizeof P.inverse);
    *raw_indices = P.raw_indices;
    *offsets = P.offsets;
    *equivalents = P.equivalents;
    *negated = P.negated;
    return P.count;
}
