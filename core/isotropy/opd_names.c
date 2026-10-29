/* The isotropy table as a catalogue of names: the row of an irrep with tabulated rows whose direction space lies in the
 * G-orbit of a computed one. The free count and the subgroup number decide most rows; the rest are told apart by the
 * orbit of the computed space under the image group (the irrep matrices of the parent operations and translations;
 * for the magnetic table the gray group, a time-reversed element acting as minus its ordinary part). */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "irreps.h"
#include "magnetic.h"
#include "tables.h"

int isotropy_row_subgroup(int64_t row_id);
int isotropy_row_free(int64_t row_id);
int mag_isotropy_row_subgroup(int64_t row_id);
int mag_isotropy_row_free(int64_t row_id);
int little_old_id(int gid);
int little_full_dim(int gid);

enum { GROUP_CAP = 20000 };

static int rref(double *B, int f, int D)
{
    int rank = 0;
    for (int col = 0; col < D && rank < f; col++) {
        int best = -1;
        double bv = 1e-7;
        for (int r = rank; r < f; r++)
            if (fabs(B[r * D + col]) > bv) bv = fabs(B[r * D + col]), best = r;
        if (best < 0) continue;
        for (int j = 0; j < D; j++) {
            double t = B[rank * D + j];
            B[rank * D + j] = B[best * D + j];
            B[best * D + j] = t;
        }
        double pv = B[rank * D + col];
        for (int j = 0; j < D; j++) B[rank * D + j] /= pv;
        for (int r = 0; r < f; r++)
            if (r != rank && B[r * D + col] != 0) {
                double x = B[r * D + col];
                for (int j = 0; j < D; j++) B[r * D + j] -= x * B[rank * D + j];
            }
        rank++;
    }
    return rank;
}

static int same_rows(const double *a, const double *b, int n)
{
    for (int i = 0; i < n; i++)
        if (fabs(a[i] - b[i]) > 1e-6) return 0;
    return 1;
}


static uint64_t matrix_key(const double *m, int D)
{
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < D * D; i++) h = (h ^ (uint64_t)llround(m[i] * 1e5)) * 1099511628211ULL;
    return h;
}

/* the image group's elements (closure of the generator matrices), at most GROUP_CAP; returns the order or -1 */
static int image_group(int gid, int magnetic, int D, double **out)
{
    int sg = little_irrep_sg(gid);
    int64_t gens[200][5], sign[200];
    int ng;
    if (magnetic) {
        ng = magnetic_records_iso(mag_orderparam_group(sg), 192, gens[0]);
        for (int i = 0; i < ng; i++) {
            int s = 1;
            gens[i][4] = mag_point_action(1, (int)gens[i][4], &s);
            sign[i] = s;
        }
    } else {
        ng = geo_sg_records(sg, gens);
        for (int i = 0; i < ng; i++) sign[i] = 1;
    }
    if (ng <= 0) return -1;
    int64_t ids[GEO_MAX_OPS][5];
    geo_sg_records(sg, ids);
    for (int a = 0; a < 3; a++) {
        sign[ng] = 1;
        int64_t *t = gens[ng++];
        t[0] = a == 0, t[1] = a == 1, t[2] = a == 2, t[3] = 1, t[4] = ids[0][4];
    }
    int *r32 = malloc(sizeof(int) * 5 * (size_t)ng);
    double *gm = malloc(sizeof(double) * (size_t)ng * D * D);
    for (int i = 0; i < ng; i++)
        for (int j = 0; j < 5; j++) r32[5 * i + j] = (int)gens[i][j];
    const int64_t none[1] = {0};
    int status = case_irrep_stack(gid, 0, none, ng, r32, (int64_t)ng * D * D, gm);
    free(r32);
    if (status) {
        free(gm);
        return -1;
    }
    for (int i = 0; i < ng; i++)
        for (int j = 0; j < D * D; j++) gm[(size_t)i * D * D + j] *= (double)sign[i];
    double *mats = malloc(sizeof(double) * (size_t)GROUP_CAP * D * D), *p = malloc(sizeof(double) * D * D);
    enum { SLOTS = 4 * GROUP_CAP + 1 };
    int *slot = malloc(sizeof(int) * SLOTS), n = 1;
    for (int i = 0; i < SLOTS; i++) slot[i] = -1;
    memset(mats, 0, sizeof(double) * D * D);
    for (int i = 0; i < D; i++) mats[i * D + i] = 1;
    slot[matrix_key(mats, D) % SLOTS] = 0;
    for (int i = 0; i < n; i++)
        for (int g = 0; g < ng; g++) {
            for (int r = 0; r < D; r++)
                for (int c = 0; c < D; c++) {
                    double s = 0;
                    for (int k = 0; k < D; k++) s += mats[(size_t)i * D * D + r * D + k] * gm[(size_t)g * D * D + k * D + c];
                    p[r * D + c] = s;
                }
            uint64_t s = matrix_key(p, D) % SLOTS;
            int known = 0;
            while (slot[s] >= 0 && !(known = same_rows(mats + (size_t)slot[s] * D * D, p, D * D))) s = (s + 1) % SLOTS;
            if (known) continue;
            if (n == GROUP_CAP) {
                free(mats), free(p), free(gm), free(slot);
                return -1;
            }
            slot[s] = n;
            memcpy(mats + (size_t)n++ * D * D, p, sizeof(double) * D * D);
        }
    free(p), free(gm), free(slot);
    *out = mats;
    return n;
}

/* the table row (1-based) of gid whose direction space is in the orbit of the f rows of `space`, with that free count
 * and subgroup; 0 when none, -1 when the orbit cannot be computed */
int opd_table_row(int gid, int magnetic, int64_t free_count, int64_t subgroup, const double *space)
{
    section_view ip, opp, op;
    const char *family = magnetic ? "isotropy.magnetic." : "isotropy.ordinary.";
    char name[3][64];
    static const char *const part[3] = {"by_old_start", "direction_codes_start", "direction_codes"};
    for (int i = 0; i < 3; i++) {
        strcpy(name[i], family);
        strcat(name[i], part[i]);
    }
    int old = little_old_id(gid), D = little_full_dim(gid);
    if (old <= 0 || D <= 0 || section_find(name[0], &ip) || section_find(name[1], &opp) || section_find(name[2], &op)) return 0;
    int64_t candidates[256];
    int nc = 0;
    for (int64_t row = section_int(&ip, old - 1); row < section_int(&ip, old) && nc < 256; row++) {
        int64_t f = magnetic ? mag_isotropy_row_free(row) : isotropy_row_free(row);
        int64_t s = magnetic ? mag_isotropy_row_subgroup(row) : isotropy_row_subgroup(row);
        if (f == free_count && s == subgroup) candidates[nc++] = row;
    }
    if (nc <= 1) return nc ? (int)candidates[0] : 0;
    int f = (int)free_count;
    double *mine = malloc(sizeof(double) * f * D), *theirs = malloc(sizeof(double) * f * D), *image = malloc(sizeof(double) * f * D);
    double *group = 0;
    int order = image_group(gid, magnetic, D, &group), found = order < 0 ? -1 : 0;
    for (int c = 0; c < nc && found == 0; c++) {
        int64_t p = section_int(&opp, candidates[c] - 1) - 1;
        for (int r = 0; r < f; r++)
            for (int j = 0; j < D; j++) theirs[r * D + j] = const_ref_value(section_int(&op, p + r * D + j));
        if (rref(theirs, f, D) != f) continue;
        for (int e = 0; e < order && !found; e++) {
            const double *E = group + (size_t)e * D * D;
            for (int r = 0; r < f; r++)
                for (int j = 0; j < D; j++) {
                    double x = 0;
                    for (int k = 0; k < D; k++) x += space[r * D + k] * E[j * D + k];
                    image[r * D + j] = x;
                }
            if (rref(image, f, D) == f && same_rows(image, theirs, f * D)) found = (int)candidates[c];
        }
    }
    free(mine), free(theirs), free(image), free(group);
    return found;
}
