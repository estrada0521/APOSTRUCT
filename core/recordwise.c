/* project_vector_ for one project preparation, by the Source bridge
 * (pipeline_trace recordwise projection: the OPD row groups of every bridge
 * family, their bridge modes on the displayed atoms, and the Source rules
 * that interleave or pair-sum the group blocks) or by the Source DISPLAY
 * DISTORTION plan on an occurrence cell (mode_b0 then project_vector_). */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "irreps.h"
#include "pynum.h"
#include "tables.h"

void bridge_coefficients(int lanes, int dim, const double *weights, int active_length, const double *active,
                         int vector_loop_count, int row_count, double *out);
void display_mode_b0(int lanes, int dim, const double *irreps, int direction_rows, const double *direction,
                     int project_count, const double *basis, double *out);
int project_vectors(int lanes, const double *basis, int basis_length, int atom_count, const int *reps, int target_rep,
                    int vector_dim, const double *constants, const double *points, int project_count, double *output,
                    int output_length);
int bridge_source_modes(int lanes, const double *coefficients, int row_count, int row_stride, const int *row_lengths,
                        const double *rows, int vector_loop_count, int bridge_project_count, const int *reps,
                        int target_rep, int vector_dim, const double *constants, const double *points,
                        int project_count, int split_row_slots, int transform_count, const double *transforms,
                        int negate, int capacity, double *modes, int *mode_count);

enum {
    RW_OK = 0,
    RW_MEMORY = 1,
    RW_NO_ACTIVE = 2,
    RW_TABLES = 3,
    RW_INDEX = 4,
    RW_CAPACITY = 5,
    RW_STACK = 6,
    RW_MODES = 7,
};

#define COEFFICIENTS 181

typedef struct {
    int start, count, row_count;
} group;

static int projection_tables(int site_pg, int vector_basis_id, int vector_setting, int vector_dim, int target,
                             int lanes, const int *records, int *reps, double *constants, double *points)
{
    static section_view reps_table, basis_codes, point_table, values;
    static int ready;
    if (!ready) {
        if (section_find("wyckoff.iwyckoff_pg_vector_reps", &reps_table) ||
            section_find("wyckoff.iwyckoff_pg_vector_basis", &basis_codes) ||
            section_find("space.ipoint_op", &point_table) || section_find("const.values", &values))
            return RW_TABLES;
        ready = 1;
    }
    const double *constant = (const double *)values.data;
    int64_t constant_count = section_count(&values);
    int step = vector_dim > 1 ? vector_dim : 1;
    reps[0] = reps[1] = reps[2] = 0;
    for (int start = 0; start < 3; start += step) {
        int64_t index = (int64_t)(site_pg - 1) * 6 + (vector_setting - 1) * 3 + start;
        if (index < 0 || index >= section_count(&reps_table)) return RW_INDEX;
        reps[start] = (int)section_int(&reps_table, index);
    }
    memset(constants, 0, sizeof(double) * 3 * (size_t)vector_dim * 3);
    for (int start = 0; start < 3; start += step) {
        if (reps[start] != target) continue;
        for (int column = 0; column < vector_dim; column++)
            for (int axis = 0; axis < 3; axis++) {
                int64_t index = (int64_t)(vector_basis_id - 1) * 18 + (vector_setting - 1) * 9 + start * 3 + column * 3 + axis;
                if (index < 0 || index >= section_count(&basis_codes)) return RW_INDEX;
                int64_t code = section_int(&basis_codes, index);
                if (code < 0 || code >= constant_count) return RW_INDEX;
                constants[(start * vector_dim + column) * 3 + axis] = constant[code];
            }
    }
    for (int lane = 0; lane < lanes; lane++) {
        int op = records[5 * lane + 4];
        for (int i = 0; i < 9; i++) {
            int64_t index = 9 * (int64_t)(op - 1) + i;
            if (op < 1 || index >= section_count(&point_table)) return RW_INDEX;
            points[9 * lane + i] = (double)section_int(&point_table, index);
        }
    }
    return 0;
}

static int effective_rows(const double *rows, const int *lengths, int row_count, int stride)
{
    if (!row_count) return 0;
    for (int r = 1; r < row_count; r++) {
        if (lengths[r] != lengths[0]) return row_count;
        for (int c = 0; c < lengths[0]; c++)
            if (python_round(rows[r * stride + c], 12) != python_round(rows[c], 12)) return row_count;
    }
    return 1;
}

static int add_group(group *groups, int *count, int start, int length, int row_count)
{
    if (*count >= 96) return 1;
    groups[(*count)++] = (group){start, length, row_count};
    return 0;
}

static void normalize(double *mode, int lanes)
{
    double scale = 0.0;
    for (int i = 0; i < 3 * lanes; i++)
        if (fabs(mode[i]) > scale) scale = fabs(mode[i]);
    if (scale <= 1e-12 || fabs(scale - 1.0) <= 1e-12) return;
    for (int i = 0; i < 3 * lanes; i++) mode[i] = mode[i] / scale;
}

/* little: gid, irrep type, full dim, old id; sources: family, bridge
 * project count, vector loop count per family in family order; active:
 * the prep's basis (index, value) pairs; rows: the order-parameter source
 * rows (row_total x full_dim); records: the displayed atom records.
 * Writes modes (mode x lane x 3) with summed[mode] marking pair sums, and per
 * source: family, source row count, group count, groups (start, length,
 * row count), block mode counts.  Returns the mode count or a negative
 * status. */
int recordwise_modes(int gid, int irrep_type, int full_dim, int old_id, int k_params, int param_count,
                     const int64_t *params, int pg_irrep, int project_count, int vector_dim, int vector_basis_id,
                     int site_pg, int vector_setting, int source_count, const int *sources, int active_count,
                     const int64_t *active_index, const double *active_value, int row_total, int row_stride,
                     const int *row_lengths, const double *rows,
                     int source_row_id, int lanes, const int *records, int transform_count, const double *transforms,
                     int negate, int emit_freeparam, int freeparam, int emit_all_rows, int emit_all_families,
                     int64_t capacity, double *modes, int *summed, int64_t summary_capacity, int64_t *summary)
{
    if (!source_count || !row_total) return 0;
    int source_row_count = effective_rows(rows, row_lengths, row_total, row_stride);
    if (source_row_count <= 0) return 0;
    int identity_pair = irrep_type == 1 && k_params && full_dim == 4 && (vector_dim == 1 || vector_dim == 2) &&
                        source_row_count == 2 && row_lengths[0] >= 2 && row_lengths[1] >= 2;
    if (identity_pair) {
        identity_pair = fabs(rows[0] - 1.0) <= 1e-12 && fabs(rows[row_stride + 1] - 1.0) <= 1e-12;
        for (int c = 0; c < row_lengths[0] && identity_pair; c++)
            if (c != 0 && !(fabs(rows[c]) <= 1e-12)) identity_pair = 0;
        for (int c = 0; c < row_lengths[1] && identity_pair; c++)
            if (c != 1 && !(fabs(rows[row_stride + c]) <= 1e-12)) identity_pair = 0;
    }
    int status = 0;
    int64_t mode_count = 0, summary_used = 0;
    double *weights = malloc(sizeof(double) * (size_t)lanes * (size_t)full_dim * (size_t)full_dim + 1);
    double *coefficients = malloc(sizeof(double) * (size_t)lanes * COEFFICIENTS + 1);
    double *points = malloc(sizeof(double) * 9 * (size_t)lanes + 1);
    double constants[3 * 3 * 3];
    int reps[3];
    if (!weights || !coefficients || !points) status = RW_MEMORY;
    if (!status && case_irrep_stack(gid, param_count, params, lanes, records,
                                    (int64_t)lanes * full_dim * full_dim, weights))
        status = RW_STACK;
    int freeparam_groups = 0;
    if (emit_freeparam && source_row_id > 0 && freeparam >= 0) {
        int divisor = irrep_type == 1 ? 1 : irrep_type == 2 ? 4 : 2;
        freeparam_groups = freeparam / divisor;
        if (freeparam_groups < 0) freeparam_groups = 0;
    }
    for (int order = 0; order < source_count && !status; order++) {
        const int *source = sources + 3 * (identity_pair ? source_count - 1 - order : order);
        int family = source[0], bridge_count = source[1], loops = source[2];
        int two_row_slots = irrep_type == 3 && source_row_count >= 2 &&
                            (source_row_count == 2 || lanes <= 4 ||
                             (source_row_count >= 4 && (project_count == 1 || project_count == 2) &&
                              (vector_dim == 1 || vector_dim == 2)));
        group groups[96];
        int group_count = 0, n = source_row_count, over = 0;
        if (emit_freeparam && source_row_id > 0 && irrep_type == 2 && n > 1) {
            for (int i = 0; i < n; i++) over |= add_group(groups, &group_count, i, 1, 1);
        } else if (freeparam_groups > 0 && irrep_type == 1) {
            for (int i = 0; i < (n < freeparam_groups ? n : freeparam_groups); i++) over |= add_group(groups, &group_count, i, 1, 1);
        } else if (freeparam_groups > 0 && irrep_type == 3) {
            int limit = n < freeparam_groups * 2 ? n : freeparam_groups * 2;
            for (int i = 0; i < limit; i += 2)
                if (i + 1 < n) over |= add_group(groups, &group_count, i, 2, 2);
        } else if (freeparam_groups > 0) {
            for (int i = 0; i < (n < freeparam_groups ? n : freeparam_groups); i++) over |= add_group(groups, &group_count, i, 1, 1);
        } else if (emit_all_rows && irrep_type == 1) {
            for (int i = 0; i < n; i++) over |= add_group(groups, &group_count, i, 1, 1);
        } else if (emit_all_rows && irrep_type == 3 && n >= 2 && n % 2 == 0) {
            for (int i = 0; i < n; i += 2) over |= add_group(groups, &group_count, i, 2, 2);
        } else if (emit_all_rows) {
            for (int i = 0; i < n; i++) over |= add_group(groups, &group_count, i, 1, 1);
        } else if (source_row_id <= 0 && old_id > 0 && irrep_type == 1 && n > 1) {
            for (int i = 0; i < n; i++) over |= add_group(groups, &group_count, i, 1, 1);
        } else if (irrep_type == 1) {
            over |= add_group(groups, &group_count, 0, 1, 1);
        } else if (two_row_slots && n >= 4 && k_params && irrep_type == 3 && full_dim == 8) {
            int pairs = n / 2 > 1 ? n / 2 : 1;
            int start = 2 * (family % pairs);
            int length = start + 2 <= row_total ? 2 : row_total - start;
            over |= add_group(groups, &group_count, start, length < 0 ? 0 : length, 2);
        } else if (two_row_slots && n >= 4 && k_params) {
            for (int i = 0; i < n; i += 2) over |= add_group(groups, &group_count, i, i + 2 <= row_total ? 2 : row_total - i, 2);
        } else if (two_row_slots) {
            over |= add_group(groups, &group_count, 0, row_total < 2 ? row_total : 2, 2);
        } else if (irrep_type == 3 && k_params && n >= 4) {
            for (int i = 0; i < n; i += 2) over |= add_group(groups, &group_count, i, i + 2 <= row_total ? 2 : row_total - i, 2);
        } else if (irrep_type == 3 && n >= 4) {
            over |= add_group(groups, &group_count, 0, 2, 2);
        } else {
            for (int i = 0; i < n; i++) over |= add_group(groups, &group_count, i, 1, 1);
        }
        if (over) {
            status = RW_CAPACITY;
            break;
        }
        int active_length = 144 > (bridge_count + 1) * 3 ? 144 : (bridge_count + 1) * 3;
        double *active = calloc((size_t)active_length, sizeof(double));
        if (!active) {
            status = RW_MEMORY;
            break;
        }
        int found = 0;
        for (int i = 0; i < active_count; i++) {
            if (active_index[i] / 144 != family) continue;
            int64_t local = active_index[i] % 144;
            if (local < active_length) {
                active[local] = active_value[i];
                found = 1;
            }
        }
        if (!found) {
            free(active);
            status = RW_NO_ACTIVE;
            break;
        }
        bridge_coefficients(lanes, full_dim, weights, active_length, active, loops, bridge_count + 1, coefficients);
        free(active);
        status = projection_tables(site_pg, vector_basis_id, vector_setting, vector_dim, pg_irrep, lanes, records, reps,
                                   constants, points);
        if (status) break;
        int block_start[96], block_count[96];
        int64_t block_base = mode_count;
        int family_project_count = emit_all_families ? project_count : 1;
        for (int g = 0; g < group_count && !status; g++) {
            int lengths[2] = {0, 0}, stride = 1;
            double row_values[2 * 64] = {0};
            for (int r = 0; r < groups[g].row_count; r++)
                if (r < groups[g].count) {
                    lengths[r] = row_lengths[groups[g].start + r];
                    if (lengths[r] > stride) stride = lengths[r];
                }
            if (stride > 64) {
                status = RW_CAPACITY;
                break;
            }
            for (int r = 0; r < groups[g].row_count; r++)
                for (int c = 0; c < lengths[r]; c++) row_values[r * stride + c] = rows[(groups[g].start + r) * row_stride + c];
            int block_capacity = 3 * (family_project_count > 1 ? family_project_count : 1) * groups[g].row_count + 3;
            if (mode_count + block_capacity > capacity) {
                status = RW_CAPACITY;
                break;
            }
            int produced = 0;
            int bs = bridge_source_modes(lanes, coefficients, groups[g].row_count, stride, lengths, row_values, loops,
                                         bridge_count, reps, pg_irrep, vector_dim, constants, points,
                                         family_project_count, irrep_type == 3 && groups[g].row_count > 1,
                                         transform_count, transforms, negate, block_capacity,
                                         modes + (size_t)mode_count * lanes * 3, &produced);
            if (bs) {
                status = RW_MODES;
                break;
            }
            block_start[g] = (int)mode_count;
            block_count[g] = produced;
            mode_count += produced;
        }
        if (status) break;
        if (summary_used + 3 + 3 * group_count + group_count > summary_capacity) {
            status = RW_CAPACITY;
            break;
        }
        summary[summary_used++] = family;
        summary[summary_used++] = source_row_count;
        summary[summary_used++] = group_count;
        for (int g = 0; g < group_count; g++) {
            summary[summary_used++] = groups[g].start;
            summary[summary_used++] = groups[g].count;
            summary[summary_used++] = groups[g].row_count;
        }
        for (int g = 0; g < group_count; g++) summary[summary_used++] = block_count[g];
        /* Arrange this family's blocks into its emitted modes. */
        int64_t produced_total = mode_count - block_base;
        double *staging = malloc(sizeof(double) * (size_t)(produced_total + 1) * lanes * 3);
        if (!staging) {
            status = RW_MEMORY;
            break;
        }
        memcpy(staging, modes + (size_t)block_base * lanes * 3, sizeof(double) * (size_t)produced_total * lanes * 3);
        int64_t written = block_base;
        size_t mode_size = (size_t)lanes * 3;
#define BLOCK(g, m) (staging + (size_t)(block_start[g] - block_base + (m)) * mode_size)
        if (irrep_type == 2 && n > 2 && n % 2 == 0 && group_count == n && group_count) {
            int emitted = block_count[0];
            for (int g = 1; g < group_count; g++)
                if (block_count[g] < emitted) emitted = block_count[g];
            for (int m = 0; m < emitted; m++)
                for (int g = 0; g < group_count; g++) {
                    memcpy(modes + (size_t)written * mode_size, BLOCK(g, m), sizeof(double) * mode_size);
                    summed[written++] = 0;
                }
        } else if (irrep_type == 3 && k_params && n > 2 && n % 2 == 0 && group_count == n && group_count) {
            for (int g = 0; g < group_count; g++)
                for (int m = 0; m < block_count[g]; m++) {
                    memcpy(modes + (size_t)written * mode_size, BLOCK(g, m), sizeof(double) * mode_size);
                    summed[written++] = 0;
                }
        } else if (irrep_type != 1 && n > 2 && n % 2 == 0 && group_count == n) {
            for (int g = 0; g + 1 < group_count; g += 2) {
                int pairs = block_count[g] < block_count[g + 1] ? block_count[g] : block_count[g + 1];
                for (int m = 0; m < pairs; m++) {
                    double *target = modes + (size_t)written * mode_size;
                    const double *left = BLOCK(g, m), *right = BLOCK(g + 1, m);
                    for (size_t i = 0; i < mode_size; i++) target[i] = left[i] + right[i];
                    normalize(target, lanes);
                    summed[written++] = 1;
                }
            }
        } else if (irrep_type == 3 && n == 2 && group_count == 2 && block_count[0] == block_count[1] &&
                   block_count[0] > 1) {
            for (int m = 0; m < block_count[0]; m++)
                for (int g = 0; g < 2; g++) {
                    memcpy(modes + (size_t)written * mode_size, BLOCK(g, m), sizeof(double) * mode_size);
                    summed[written++] = 0;
                }
        } else {
            for (int g = 0; g < group_count; g++)
                for (int m = 0; m < block_count[g]; m++) {
                    memcpy(modes + (size_t)written * mode_size, BLOCK(g, m), sizeof(double) * mode_size);
                    summed[written++] = 0;
                }
        }
#undef BLOCK
        free(staging);
        mode_count = written;
    }
    free(weights);
    free(coefficients);
    free(points);
    return status ? -status : (int)mode_count;
}

/* display_distortion.source_display_projection_families: modes_out
 * (project_count x direction_rows x lanes x 3) when every family projects
 * to a nonzero mode; returns 1 then, 0 when the display plan declines. */
int display_families(int gid, int full_dim, int site_pg, int pg_irrep, int vector_basis_id, int vector_dim,
                     int vector_setting, int project_count, int param_count, const int64_t *params, int active_count,
                     const int64_t *active_index, const double *active_value, int direction_rows,
                     const double *direction, int lanes, const int *records, double *modes_out)
{
    if (gid <= 0 || full_dim <= 0 || full_dim > 48 || site_pg <= 0 || pg_irrep <= 0 || vector_basis_id <= 0 ||
        vector_dim < 1 || vector_dim > 3 || vector_setting < 1 || vector_setting > 2 || project_count <= 0 || lanes <= 0)
        return 0;
    if (direction_rows <= 0 || direction_rows > 48) return 0;
    for (int i = 0; i < direction_rows * full_dim; i++)
        if (!isfinite(direction[i])) return 0;
    int basis_length = project_count * 144;
    double *basis = calloc((size_t)basis_length, sizeof(double));
    char *seen = calloc((size_t)basis_length, 1);
    char *families = calloc((size_t)project_count, 1);
    double *weights = malloc(sizeof(double) * (size_t)lanes * full_dim * full_dim);
    double *mode_b0 = malloc(sizeof(double) * (size_t)lanes * basis_length);
    double *output = malloc(sizeof(double) * (size_t)lanes * basis_length);
    double *points = malloc(sizeof(double) * 9 * (size_t)lanes);
    int ok = basis && seen && families && weights && mode_b0 && output && points && active_count > 0;
    for (int i = 0; i < active_count && ok; i++) {
        int64_t index = active_index[i];
        double value = active_value[i];
        if (index < 0 || index >= basis_length || index / 144 >= project_count || index % 144 >= full_dim * 3 ||
            !isfinite(value) || seen[index])
            ok = 0;
        else {
            seen[index] = 1;
            basis[index] = value;
            if (value != 0.0) families[index / 144] = 1;
        }
    }
    for (int f = 0; f < project_count && ok; f++) ok = families[f];
    for (int r = 0; r < lanes && ok; r++) ok = records[5 * r + 3] > 0 && records[5 * r + 4] > 0;
    if (ok && case_irrep_stack(gid, param_count, params, lanes, records, (int64_t)lanes * full_dim * full_dim, weights))
        ok = 0;
    for (int i = 0; i < lanes * full_dim * full_dim && ok; i++) ok = isfinite(weights[i]);
    if (ok) {
        display_mode_b0(lanes, full_dim, weights, direction_rows, direction, project_count, basis, mode_b0);
        for (int i = 0; i < lanes * basis_length && ok; i++) ok = isfinite(mode_b0[i]);
    }
    int reps[3];
    double constants[27];
    if (ok && projection_tables(site_pg, vector_basis_id, vector_setting, vector_dim, pg_irrep, lanes, records, reps,
                                constants, points))
        ok = 0;
    if (ok) {
        int emitted = project_vectors(lanes, mode_b0, basis_length, direction_rows, reps, pg_irrep, vector_dim,
                                      constants, points, project_count, output, basis_length);
        ok = emitted == project_count;
        for (int i = 0; i < lanes * basis_length && ok; i++) ok = isfinite(output[i]);
    }
    for (int f = 0; f < project_count && ok; f++) {
        double scale = 0.0;
        for (int row = 0; row < direction_rows; row++)
            for (int lane = 0; lane < lanes; lane++)
                for (int axis = 0; axis < 3; axis++) {
                    double value = output[(size_t)lane * basis_length + f * 144 + row * 3 + axis];
                    modes_out[(((size_t)f * direction_rows + row) * lanes + lane) * 3 + axis] = value;
                    if (fabs(value) > scale) scale = fabs(value);
                }
        ok = scale > 1e-12;
    }
    free(basis);
    free(seen);
    free(families);
    free(weights);
    free(mode_b0);
    free(output);
    free(points);
    return ok;
}

/* pipeline_trace finalize_vector over count vectors in place: each
 * transform as 0.0 + x m0c + y m1c + z m2c, then the optional negation. */
void finalize_vectors(int64_t count, double *vectors, int transform_count, const double *transforms, int negate)
{
    for (int64_t i = 0; i < count; i++) {
        double *v = vectors + 3 * i;
        for (int t = 0; t < transform_count; t++) {
            const double *m = transforms + 9 * t;
            double x = v[0], y = v[1], z = v[2];
            for (int c = 0; c < 3; c++) v[c] = 0.0 + x * m[c] + y * m[3 + c] + z * m[6 + c];
        }
        if (negate)
            for (int c = 0; c < 3; c++) v[c] = -v[c];
    }
}
