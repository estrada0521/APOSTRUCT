/* One projected mode block of pipeline_trace: the project_ inputs of every
 * site irrep of a little irrep, the Source DISPLAY DISTORTION plan on the
 * occurrence cell when it is complete, else the recordwise bridge
 * projection on the baseline cell, with the block summaries. */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "kernels.h"

enum {
    MB_MEMORY = 1,
    MB_NO_MODES = 2,
    MB_NO_BASIS_ID = 3,
    MB_ENTRY = 100,
    MB_RECORDWISE = 200,
};

typedef struct {
    int pg_irrep, project_count, vector_dim, prep_status, basis_ok, family_count, active_count;
    int64_t basis_id;
    int64_t *out, *index;
    double *value;
    const int64_t *families; /* family, bridge count, loop count, members */
    int64_t source_at;
} prep;

static int64_t *ints;
static size_t ints_used, ints_cap;
static double *doubles;
static size_t doubles_used, doubles_cap;
static int *summed_flags;

static int64_t *grow_ints(size_t count)
{
    if (ints_used + count > ints_cap) {
        size_t cap = ints_cap ? ints_cap : 1024;
        while (ints_used + count > cap) cap *= 2;
        int64_t *grown = realloc(ints, cap * sizeof *grown);
        if (!grown) return 0;
        ints = grown;
        ints_cap = cap;
    }
    int64_t *at = ints + ints_used;
    ints_used += count;
    return at;
}

static double *grow_doubles(size_t count)
{
    if (doubles_used + count > doubles_cap) {
        size_t cap = doubles_cap ? doubles_cap : 4096;
        while (doubles_used + count > cap) cap *= 2;
        double *grown = realloc(doubles, cap * sizeof *grown);
        if (!grown) return 0;
        doubles = grown;
        doubles_cap = cap;
    }
    double *at = doubles + doubles_used;
    doubles_used += count;
    return at;
}

static int load_prep(int sg, int row_id, int site_pg, int gid, int pg_irrep, int param_count, const int64_t *params,
                     const double *inverse, int inverse_error, prep *p)
{
    memset(p, 0, sizeof *p);
    p->pg_irrep = pg_irrep;
    int status = 4;
    for (int64_t capacity = 1 << 13; status == 4 && capacity <= 1 << 22; capacity <<= 2) {
        free(p->out);
        free(p->index);
        free(p->value);
        p->out = malloc(sizeof(int64_t) * capacity);
        p->index = malloc(sizeof(int64_t) * capacity);
        p->value = malloc(sizeof(double) * capacity);
        if (!p->out || !p->index || !p->value) return MB_MEMORY;
        status = projection_entry(sg, row_id, site_pg, gid, pg_irrep, param_count, params, inverse, capacity, p->out,
                                  p->index, p->value);
    }
    if (status) return MB_ENTRY + status;
    const int64_t *q = p->out;
    p->project_count = (int)*q++;
    q++;
    p->vector_dim = (int)*q++;
    p->prep_status = (int)*q++;
    p->basis_id = *q++;
    q += *q + 1;
    int64_t basis_status = *q++;
    p->basis_ok = !basis_status && !inverse_error;
    if (!p->basis_ok) return 0;
    q += p->project_count;
    p->active_count = (int)*q++;
    p->family_count = (int)*q++;
    p->families = q;
    p->source_at = p->active_count;
    return 0;
}

static void free_prep(prep *p)
{
    free(p->out);
    free(p->index);
    free(p->value);
}

/* _freeze_projection_preps: every prep has a basis, a vector basis and one
 * bridge source per project family in family order. */
static int display_preps_complete(int site_pg, const prep *preps, int count)
{
    if (!count) return 0;
    for (int i = 0; i < count; i++) {
        const prep *p = preps + i;
        if (p->pg_irrep <= 0 || p->project_count <= 0 || site_pg <= 0 || p->prep_status || p->basis_id <= 0 ||
            p->vector_dim < 1 || p->vector_dim > 3 || !p->basis_ok || p->family_count != p->project_count)
            return 0;
        for (int f = 0; f < p->family_count; f++)
            if (p->families[4 * f] != f || p->families[4 * f + 1] < 0) return 0;
    }
    return 1;
}

static int display_block(int gid, int full_dim, int site_pg, int vector_setting, int param_count,
                         const int64_t *params, const prep *preps, int count, int direction_rows,
                         const double *direction, int lanes, const int *records, int transform_count,
                         const double *transforms, int negate, int *adopted)
{
    *adopted = 0;
    if (!display_preps_complete(site_pg, preps, count)) return 0;
    size_t total = 0;
    for (int i = 0; i < count; i++) total += (size_t)preps[i].project_count * direction_rows * lanes * 3;
    double *modes = malloc(sizeof(double) * (total + 1));
    if (!modes) return MB_MEMORY;
    size_t at = 0;
    for (int i = 0; i < count; i++) {
        const prep *p = preps + i;
        if (!display_families(gid, full_dim, site_pg, p->pg_irrep, (int)p->basis_id, p->vector_dim, vector_setting,
                              p->project_count, param_count, params, p->active_count, p->index, p->value,
                              direction_rows, direction, lanes, records, modes + at)) {
            free(modes);
            return 0;
        }
        at += (size_t)p->project_count * direction_rows * lanes * 3;
    }
    finalize_vectors((int64_t)(total / 3), modes, transform_count, transforms, negate);
    int64_t summaries = 0;
    for (int i = 0; i < count; i++) summaries += preps[i].project_count;
    int64_t *head = grow_ints(2);
    if (!head) {
        free(modes);
        return MB_MEMORY;
    }
    head[0] = 1;
    head[1] = summaries;
    for (int i = 0; i < count; i++) {
        const prep *p = preps + i;
        for (int f = 0; f < p->project_count; f++) {
            int64_t *s = grow_ints(9 + 4 * (size_t)direction_rows);
            if (!s) {
                free(modes);
                return MB_MEMORY;
            }
            *s++ = p->pg_irrep;
            *s++ = f;
            *s++ = direction_rows;
            *s++ = direction_rows;
            for (int r = 0; r < direction_rows; r++) {
                *s++ = r;
                *s++ = 1;
                *s++ = 1;
            }
            for (int r = 0; r < direction_rows; r++) *s++ = 1;
            *s++ = direction_rows;
            *s++ = p->project_count;
            *s++ = p->project_count;
            *s++ = p->families[4 * f + 1];
            *s++ = p->vector_dim;
        }
    }
    int64_t mode_count = (int64_t)(total / ((size_t)lanes * 3));
    int64_t *tail = grow_ints(1 + (size_t)mode_count);
    double *values = grow_doubles(total);
    if (!tail || !values) {
        free(modes);
        return MB_MEMORY;
    }
    *tail++ = mode_count;
    for (int64_t m = 0; m < mode_count; m++) *tail++ = 0;
    memcpy(values, modes, sizeof(double) * total);
    free(modes);
    *adopted = 1;
    return 0;
}

static int bridge_block(int gid, int irrep_type, int full_dim, int old_id, int site_pg, int vector_setting,
                        int param_count, const int64_t *params, const prep *preps, int count, int row_total,
                        int row_stride, const int *row_lengths, const double *rows, int has_rows, int source_row_id,
                        int freeparam, int flags, int lanes, const int *records, int transform_count,
                        const double *transforms, int negate)
{
    int64_t *head = grow_ints(2);
    if (!head) return MB_MEMORY;
    size_t head_at = ints_used - 2;
    ints[head_at] = 0;
    int64_t summaries = 0, mode_total = 0;
    size_t mode_size = (size_t)lanes * 3;
    int64_t mode_capacity = 0;
    for (int i = 0; i < count; i++) {
        int64_t c = (3 * (int64_t)(preps[i].project_count > 1 ? preps[i].project_count : 1) * 2 + 3) *
                        (row_total > 1 ? row_total : 1) * (preps[i].family_count > 1 ? preps[i].family_count : 1) +
                    16;
        if (c > mode_capacity) mode_capacity = c;
    }
    double *modes = malloc(sizeof(double) * ((size_t)mode_capacity * mode_size + 1));
    int *summed = malloc(sizeof(int) * ((size_t)mode_capacity + 1));
    int64_t summary[4096];
    int *sources = 0;
    int status = !modes || !summed ? MB_MEMORY : 0;
    for (int i = 0; i < count && !status; i++) {
        const prep *p = preps + i;
        if (irrep_type < 1 || irrep_type > 3 || p->vector_dim < 1 || p->vector_dim > 3) continue;
        if (!p->family_count || !has_rows) {
            status = MB_NO_MODES;
            break;
        }
        if (p->prep_status) {
            status = MB_NO_BASIS_ID;
            break;
        }
        int *grown = realloc(sources, sizeof(int) * (3 * (size_t)p->family_count + 3));
        if (!grown) {
            status = MB_MEMORY;
            break;
        }
        sources = grown;
        for (int f = 0; f < p->family_count; f++) {
            sources[3 * f] = (int)p->families[4 * f];
            sources[3 * f + 1] = (int)p->families[4 * f + 1];
            sources[3 * f + 2] = (int)p->families[4 * f + 2];
        }
        int64_t capacity = (3 * (int64_t)(p->project_count > 1 ? p->project_count : 1) * 2 + 3) *
                               (row_total > 1 ? row_total : 1) * (p->family_count > 1 ? p->family_count : 1) +
                           16;
        int written = recordwise_modes(gid, irrep_type, full_dim, old_id, param_count > 0, param_count, params,
                                       p->pg_irrep, p->project_count, p->vector_dim, (int)p->basis_id, site_pg,
                                       vector_setting, p->family_count, sources, p->active_count, p->index, p->value,
                                       row_total, row_stride, row_lengths, rows, source_row_id, lanes, records,
                                       transform_count, transforms, negate, flags & 1, freeparam, (flags >> 1) & 1,
                                       (flags >> 2) & 1, capacity, modes, summed, 4096, summary);
        if (written < 0) {
            status = MB_RECORDWISE - written;
            break;
        }
        if (!written) {
            status = MB_NO_MODES;
            break;
        }
        const int64_t *cursor = summary;
        for (int f = 0; f < p->family_count && !status; f++) {
            int64_t family = cursor[0], source_rows = cursor[1], groups = cursor[2];
            int64_t *s = grow_ints(9 + 4 * (size_t)groups);
            if (!s) {
                status = MB_MEMORY;
                break;
            }
            *s++ = p->pg_irrep;
            *s++ = family;
            *s++ = source_rows;
            *s++ = groups;
            memcpy(s, cursor + 3, sizeof(int64_t) * 3 * (size_t)groups);
            s += 3 * groups;
            int64_t flat = 0;
            for (int64_t g = 0; g < groups; g++) {
                *s++ = cursor[3 + 3 * groups + g];
                flat += cursor[3 + 3 * groups + g];
            }
            *s++ = flat;
            *s++ = -1;
            *s++ = p->project_count;
            int64_t bridge = 0;
            for (int k = 0; k < p->family_count; k++)
                if (p->families[4 * k] == family) bridge = p->families[4 * k + 1];
            *s++ = bridge;
            *s++ = p->vector_dim;
            cursor += 3 + 4 * groups;
            summaries++;
        }
        double *values = grow_doubles((size_t)written * mode_size);
        int *flags_grown = realloc(summed_flags, sizeof(int) * ((size_t)(mode_total + written) + 1));
        if (!values || !flags_grown) {
            status = MB_MEMORY;
            break;
        }
        summed_flags = flags_grown;
        memcpy(values, modes, sizeof(double) * (size_t)written * mode_size);
        memcpy(summed_flags + mode_total, summed, sizeof(int) * (size_t)written);
        mode_total += written;
    }
    if (!status) {
        ints[head_at + 1] = summaries;
        int64_t *tail = grow_ints(1 + (size_t)mode_total);
        if (!tail) status = MB_MEMORY;
        else {
            *tail++ = mode_total;
            for (int64_t m = 0; m < mode_total; m++) *tail++ = summed_flags[m];
        }
    }
    free(modes);
    free(summed);
    free(sources);
    return status;
}

/* little: gid, irrep type, full dim, old id; pg_irreps: the little's project
 * entries in order; inverse: the type-2/3 transform inverse (complex
 * full_dim x full_dim) or 0 with inverse_error set when it failed.  rows: the
 * order-parameter source rows, has_rows 0 without an order-parameter source;
 * flags: emit freeparam groups, all rows, all families.  direction (display
 * rows x full_dim) with candidate records tries the display plan; baseline
 * records (baseline_lanes < 0: none) carry the bridge projection.
 * Results (mode_block_ints): adopted, summary count, per summary pg irrep,
 * family, source row count, group count, groups (start, length, row count),
 * block mode counts, flat mode count, emitted count (-1 none), project
 * count, bridge project count, vector dim; then the mode count and each
 * mode's pair-sum flag; mode_block_doubles: modes x lanes x 3.  Returns 0,
 * 1 when the display plan declined and no baseline was given, or a status. */
int mode_block(int sg, int row_id, int site_pg, int gid, int irrep_type, int full_dim, int old_id, int entry_count,
               const int *pg_irreps, int param_count, const int64_t *params, const double *inverse, int inverse_error,
               int vector_setting, int row_total, int row_stride, const int *row_lengths, const double *rows,
               int has_rows, int source_row_id, int freeparam, int flags, int direction_rows, const double *direction,
               int candidate_lanes, const int *candidate_records, int baseline_lanes, const int *baseline_records,
               int transform_count, const double *transforms, int negate)
{
    ints_used = doubles_used = 0;
    prep *preps = calloc((size_t)entry_count + 1, sizeof *preps);
    if (!preps) return -MB_MEMORY;
    int status = 0, loaded = 0;
    for (; loaded < entry_count && !status; loaded++)
        status = load_prep(sg, row_id, site_pg, gid, pg_irreps[loaded], param_count, params, inverse, inverse_error,
                           preps + loaded);
    int adopted = 0;
    if (!status && direction_rows > 0 && candidate_lanes > 0)
        status = display_block(gid, full_dim, site_pg, vector_setting, param_count, params, preps, entry_count,
                               direction_rows, direction, candidate_lanes, candidate_records, transform_count,
                               transforms, negate, &adopted);
    int declined = 0;
    if (!status && !adopted) {
        ints_used = doubles_used = 0;
        if (baseline_lanes < 0)
            declined = 1;
        else
            status = bridge_block(gid, irrep_type, full_dim, old_id, site_pg, vector_setting, param_count, params,
                                  preps, entry_count, row_total, row_stride, row_lengths, rows, has_rows,
                                  source_row_id, freeparam, flags, baseline_lanes, baseline_records, transform_count,
                                  transforms, negate);
    }
    for (int i = 0; i < loaded; i++) free_prep(preps + i);
    free(preps);
    return status ? -status : declined;
}

int64_t mode_block_ints(const int64_t **data)
{
    *data = ints;
    return (int64_t)ints_used;
}

int64_t mode_block_doubles(const double **data)
{
    *data = doubles;
    return (int64_t)doubles_used;
}
