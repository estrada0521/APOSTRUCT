/* Ferroelectric and ferroelastic classes of an isotropy subgroup: present when the subgroup fixes more of the polar
 * vector (the symmetric strain) than the parent's point group, proper when the primary irrep is at Gamma and occurs in
 * that representation. Codes: 0 none, 1 proper, 2 improper. */

#include <stdlib.h>

#include "geometry.h"
#include "tables.h"

int parent_operation_records(int parent_sg, int child_sg, const int64_t *basis, const int64_t *origin, int64_t (*out)[5]);
int little_sparse_matrix(int gid, int dim, int position, double *out);

static section_view point_op, little_k, full_dim, k_dim;
static int ready;

static void point_matrix(int op, int64_t *m)
{
    for (int k = 0; k < 9; k++) m[k] = section_int(&point_op, (int64_t)(op - 1) * 9 + k);
}

/* fixed dimensions of the polar vector and the symmetric strain under the distinct point parts of ops */
static int fixed(const int64_t (*ops)[5], int n, int64_t *polar, int64_t *strain)
{
    int seen[GEO_MAX_OPS * 4], unique = 0;
    int64_t p = 0, s = 0;
    for (int i = 0; i < n; i++) {
        int dup = 0;
        for (int j = 0; j < unique && !dup; j++) dup = seen[j] == (int)ops[i][4];
        if (dup) continue;
        seen[unique++] = (int)ops[i][4];
        int64_t m[9], squared = 0;
        point_matrix((int)ops[i][4], m);
        int64_t trace = m[0] + m[4] + m[8];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) squared += m[3 * r + c] * m[3 * c + r];
        p += trace;
        s += trace * trace + squared;
    }
    if (!unique || p % unique || s % (2 * unique)) return 1;
    *polar = p / unique;
    *strain = s / (2 * unique);
    return 0;
}

int ferroic_codes(int sg, int subgroup, const int64_t *basis, const int64_t *origin, int gid, int *fe, int *fa)
{
    *fe = *fa = 0;
    if (!ready) {
        if (geo_load() || section_find("space.ipoint_op", &point_op) || section_find("little.little_irr_k", &little_k) ||
            section_find("little.little_irr_full_dim", &full_dim) || section_find("little.little_k_dim", &k_dim))
            return 1;
        ready = 1;
    }
    int64_t G[GEO_MAX_OPS][5], H[GEO_MAX_OPS * 4][5], gp, gs, hp, hs;
    int ng = geo_sg_records(sg, G), nh = parent_operation_records(sg, subgroup, basis, origin, H);
    if (ng <= 0 || nh <= 0 || fixed((const int64_t (*)[5])G, ng, &gp, &gs) || fixed((const int64_t (*)[5])H, nh, &hp, &hs)) return 1;
    if (hp <= gp && hs <= gs) return 0;
    /* multiplicities of the primary irrep in the two representations, by characters, at Gamma */
    double mv = 0, ms = 0, norm = 0;
    int kslot = gid > 0 ? (int)section_int(&little_k, gid - 1) : 0;
    if (kslot == 1 && section_int(&k_dim, (int64_t)(geo_sg_lattice(sg) - 1) * 27) == 0) {
        int dim = (int)section_int(&full_dim, gid - 1);
        double *m = malloc(sizeof(double) * (size_t)dim * dim);
        if (!m) return 1;
        for (int q = 0; q < ng; q++) {
            if (little_sparse_matrix(gid, dim, q + 1, m)) {
                free(m);
                return 1;
            }
            double chi = 0, squared = 0;
            for (int a = 0; a < dim; a++) chi += m[a * dim + a];
            int64_t r[9];
            point_matrix((int)G[q][4], r);
            double trace = (double)(r[0] + r[4] + r[8]);
            for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++) squared += (double)(r[3 * x + y] * r[3 * y + x]);
            mv += chi * trace;
            ms += chi * (trace * trace + squared) / 2;
            norm += chi * chi;
        }
        free(m);
        if (norm > 0) mv /= norm, ms /= norm;
    }
    if (hp > gp) *fe = mv > 0.5 ? 1 : 2;
    if (hs > gs) *fa = ms > 0.5 ? 1 : 2;
    return 0;
}
