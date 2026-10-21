/* The subduction rows of the k families, ordinary and magnetic: for every k family in order and every little irrep of that k label, the
 * invariant directions over the operations plus the supercell translations,
 * merged by (gid, reciprocal point-group orbit) with the preferred source k
 * parameter taking over a row's family. */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "tables.h"

int subduction_directions(int gid, int mode_case, int point_mode, int param_count, const int64_t *params,
                          int record_count, const int *records, double *directions_out);
int reciprocal_orbit_key_of(int sg, const rat *vector, int64_t *out);

enum {
    SR_TABLES = 1,
    SR_MEMORY = 2,
    SR_DIRECTIONS = 3,
    SR_ORBIT = 4,
    SR_ZERO = 5,
};

static struct {
    section_view sg, kslot, full_dim, k_label, k_dim;
    int ready;
} T;

static int64_t *row_ints, *occurrence_ints;
static double *direction_values;
static size_t row_cap, occurrence_cap, direction_cap, row_used, occurrence_used, direction_used;
static int64_t failed_gid, failed_status;

static int load(void)
{
    if (T.ready) return 0;
    if (geo_load() || section_find("little.little_irr_space_group", &T.sg) || section_find("little.little_irr_k", &T.kslot) ||
        section_find("little.little_irr_full_dim", &T.full_dim) || section_find("little.little_k_label", &T.k_label) ||
        section_find("little.little_k_dim", &T.k_dim))
        return SR_TABLES;
    T.ready = 1;
    return 0;
}

static int grow(void **buffer, size_t *cap, size_t need, size_t size)
{
    if (need <= *cap) return 0;
    size_t next = *cap ? *cap : 256;
    while (next < need) next *= 2;
    void *grown = realloc(*buffer, next * size);
    if (!grown) return SR_MEMORY;
    *buffer = grown;
    *cap = next;
    return 0;
}

static int label_of(int lattice, int gid, const char **text, int64_t *length)
{
    int64_t kslot = section_int(&T.kslot, gid - 1);
    return section_string(&T.k_label, (int64_t)(lattice - 1) * 27 + kslot - 1, text, length);
}

static int64_t gcd64(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        int64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static int same_mod1(int64_t an, int64_t ad, int64_t bn, int64_t bd)
{
    if (ad < 0) {
        an = -an;
        ad = -ad;
    }
    if (bd < 0) {
        bn = -bn;
        bd = -bd;
    }
    __int128 left = ((__int128)an % ad + ad) % ad * bd, right = ((__int128)bn % bd + bd) % bd * ad;
    return left == right;
}

/* families: count x (gid, kparam[4], vector (num, den) x 3); preferred:
 * entries with label text and kparam[4]; ops: count x 5 with the three
 * supercell translations already appended.  use_kparam: the irrep frame
 * takes the source k parameter record (magnetic, or an isotropy irrep
 * source); otherwise the case's evaluated k parameters.  Results: rows
 * (gid, family, occurrence start, occurrence count, direction count,
 * full dim), occurrence family indices, directions (count x full dim per
 * row).  Returns 0 or a status (failed gid and status through getter). */
int subduction_rows(int sg, int family_count, const int64_t *families, int preferred_count,
                    const char *const *preferred_labels, const int64_t *preferred_lengths, const int64_t *preferred_kparams,
                    int op_count, const int *ops, int magnetic, int use_kparam)
{
    row_used = occurrence_used = direction_used = 0;
    failed_gid = failed_status = 0;
    if (load()) return SR_TABLES;
    int lattice = geo_sg_lattice(sg);
    int64_t total = section_count(&T.sg);
    int *gids = malloc(sizeof(int) * ((size_t)total + 1));
    double *buffer = malloc(sizeof(double) * 48 * 48);
    int gid_count = 0, status = !gids || !buffer ? SR_MEMORY : 0, key_count = 0;
    for (int64_t i = 0; i < total && !status; i++)
        if (section_int(&T.sg, i) == sg) gids[gid_count++] = (int)i + 1;
    int key_cap = family_count * gid_count + 1;
    int64_t (*keys)[7] = malloc(sizeof *keys * (size_t)key_cap);
    int *preferred_rows = malloc(sizeof(int) * (size_t)key_cap);
    if (!keys || !preferred_rows) status = SR_MEMORY;
    /* per row: the occurrence family indices in merge order */
    int **lists = calloc((size_t)key_cap, sizeof *lists), *list_lengths = calloc((size_t)key_cap, sizeof(int));
    int *row_family = calloc((size_t)key_cap, sizeof(int)), *row_gid = calloc((size_t)key_cap, sizeof(int));
    double **row_directions = calloc((size_t)key_cap, sizeof *row_directions);
    int *row_count = calloc((size_t)key_cap, sizeof(int)), *row_dim = calloc((size_t)key_cap, sizeof(int));
    if (!lists || !list_lengths || !row_family || !row_gid || !row_directions || !row_count || !row_dim) status = SR_MEMORY;
    for (int f = 0; f < family_count && !status; f++) {
        const int64_t *family = families + 11 * f;
        const char *label;
        int64_t label_length;
        if (label_of(lattice, (int)family[0], &label, &label_length)) {
            status = SR_TABLES;
            break;
        }
        int64_t den = family[4];
        int dim = -1, preferred = 0;
        for (int p = 0; p < preferred_count; p++) {
            if (preferred_lengths[p] != label_length || memcmp(preferred_labels[p], label, (size_t)label_length)) continue;
            const int64_t *target = preferred_kparams + 4 * p;
            if (!target[3] || !den) {
                status = SR_ZERO;
                break;
            }
            int same = 1;
            for (int axis = 0; axis < 3 && same; axis++) same = same_mod1(family[1 + axis], den, target[axis], target[3]);
            if (same) preferred = 1;
        }
        if (status) break;
        rat vector[3];
        for (int axis = 0; axis < 3; axis++) vector[axis] = rat_make(family[5 + 2 * axis], family[6 + 2 * axis]);
        for (int g = 0; g < gid_count && !status; g++) {
            int gid = gids[g];
            const char *text;
            int64_t length;
            if (label_of(lattice, gid, &text, &length)) {
                status = SR_TABLES;
                break;
            }
            if (length != label_length || memcmp(text, label, (size_t)label_length)) continue;
            if (dim < 0) dim = (int)section_int(&T.k_dim, (int64_t)(lattice - 1) * 27 + section_int(&T.kslot, gid - 1) - 1);
            int full_dim = (int)section_int(&T.full_dim, gid - 1);
            int64_t params[6];
            int param_count, mode_case;
            if (use_kparam) {
                for (int i = 0; i < 4; i++) params[i] = family[1 + i];
                param_count = 4;
                mode_case = 0;
            } else {
                mode_case = 1;
                param_count = dim > 0 ? 3 : 0;
                for (int axis = 0; axis < 3; axis++) {
                    int64_t n = axis < dim ? family[1 + axis] : 0, d = axis < dim ? den : 1;
                    if (!d) {
                        status = SR_ZERO;
                        break;
                    }
                    int64_t gg = gcd64(n, d);
                    if (d < 0) gg = -gg;
                    params[2 * axis] = gg ? n / gg : 0;
                    params[2 * axis + 1] = gg ? d / gg : 1;
                }
                if (status) break;
            }
            int count = subduction_directions(gid, mode_case, magnetic, param_count, params, op_count, ops, buffer);
            if (count < 0) {
                failed_gid = gid;
                failed_status = -count;
                status = SR_DIRECTIONS;
                break;
            }
            if (count == 0) continue;
            int64_t key[7];
            key[0] = gid;
            int orbit = reciprocal_orbit_key_of(sg, vector, key + 1);
            if (orbit) {
                failed_status = orbit;
                status = SR_ORBIT;
                break;
            }
            int index = -1;
            for (int k = 0; k < key_count && index < 0; k++)
                if (!memcmp(keys[k], key, sizeof key)) index = k;
            double *directions = malloc(sizeof(double) * ((size_t)count * full_dim + 1));
            if (!directions) {
                status = SR_MEMORY;
                break;
            }
            for (int r = 0; r < count; r++)
                for (int c = 0; c < full_dim; c++) directions[r * full_dim + c] = buffer[r * full_dim + c];
            if (index < 0) {
                if (key_count == key_cap) {
                    free(directions);
                    status = SR_MEMORY;
                    break;
                }
                index = key_count++;
                memcpy(keys[index], key, sizeof key);
                preferred_rows[index] = preferred;
                row_family[index] = f;
                row_gid[index] = gid;
                row_directions[index] = directions;
                row_count[index] = count;
                row_dim[index] = full_dim;
            } else if (preferred && !preferred_rows[index]) {
                preferred_rows[index] = 1;
                row_family[index] = f;
                free(row_directions[index]);
                row_directions[index] = directions;
                row_count[index] = count;
                row_dim[index] = full_dim;
            } else {
                free(directions);
            }
            int *grown = realloc(lists[index], sizeof(int) * ((size_t)list_lengths[index] + 1));
            if (!grown) {
                status = SR_MEMORY;
                break;
            }
            lists[index] = grown;
            lists[index][list_lengths[index]++] = f;
        }
    }
    for (int k = 0; k < key_count && !status; k++) {
        if (grow((void **)&row_ints, &row_cap, row_used + 6, sizeof(int64_t)) ||
            grow((void **)&occurrence_ints, &occurrence_cap, occurrence_used + (size_t)list_lengths[k], sizeof(int64_t)) ||
            grow((void **)&direction_values, &direction_cap, direction_used + (size_t)row_count[k] * row_dim[k], sizeof(double))) {
            status = SR_MEMORY;
            break;
        }
        int64_t *row = row_ints + row_used;
        row[0] = row_gid[k];
        row[1] = row_family[k];
        row[2] = (int64_t)occurrence_used;
        row[3] = list_lengths[k];
        row[4] = row_count[k];
        row[5] = row_dim[k];
        row_used += 6;
        for (int i = 0; i < list_lengths[k]; i++) occurrence_ints[occurrence_used++] = lists[k][i];
        memcpy(direction_values + direction_used, row_directions[k], sizeof(double) * (size_t)row_count[k] * row_dim[k]);
        direction_used += (size_t)row_count[k] * row_dim[k];
    }
    for (int k = 0; k < key_count; k++) {
        free(lists ? lists[k] : 0);
        free(row_directions ? row_directions[k] : 0);
    }
    free(lists);
    free(list_lengths);
    free(row_family);
    free(row_gid);
    free(row_directions);
    free(row_count);
    free(row_dim);
    free(gids);
    free(keys);
    free(preferred_rows);
    free(buffer);
    return status;
}

int64_t subduction_rows_ints(const int64_t **rows, const int64_t **occurrences, int64_t *occurrence_count)
{
    *rows = row_ints;
    *occurrences = occurrence_ints;
    *occurrence_count = (int64_t)occurrence_used;
    return (int64_t)(row_used / 6);
}

const double *subduction_rows_directions(void) { return direction_values; }

void subduction_rows_failure(int64_t *gid, int64_t *status)
{
    *gid = failed_gid;
    *status = failed_status;
}
