/* A subgroup's Source operations mapped into the selected raw parent embedding, integer lattice
 * translations kept. */

#include <stdint.h>
#include <string.h>

#include "geometry.h"

static void rotation_of(int sg, int op, rat (*out)[3])
{
    for (int row = 0; row < 3; row++) {
        rat unit[3] = {rat_int(row == 0), rat_int(row == 1), rat_int(row == 2)};
        geo_vrot(sg, op, unit, out[row]);
    }
}

static void matmul(const rat (*left)[3], const rat (*right)[3], rat (*out)[3])
{
    for (int row = 0; row < 3; row++) rat_row_mul(left[row], right, out[row]);
}

static int same3(const rat (*a)[3], const rat (*b)[3])
{
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
            if (!rat_eq(a[row][col], b[row][col])) return 0;
    return 1;
}

/* count >= 0 (0 for a singular basis), -10 - op when child operation op has
 * no parent point operation, -2 on table or range failure. */
int parent_operation_records(int parent_sg, int child_sg, const int64_t *basis, const int64_t *origin, int64_t (*out)[5])
{
    if (geo_load()) return -2;
    geo_take_error();
    rat_overflow = 0;
    rat b[3][3], inverse[3][3];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) b[row][col] = rat_int(basis[3 * row + col]);
    if (rat_inverse3((const rat (*)[3])b, inverse)) return 0;
    rat o[3];
    for (int axis = 0; axis < 3; axis++) o[axis] = rat_make(origin[axis], origin[3]);
    int64_t parent[GEO_MAX_OPS][5], child[GEO_MAX_OPS][5];
    int parent_count = geo_sg_records(parent_sg, parent), child_count = geo_sg_records(child_sg, child);
    if (parent_count <= 0 || child_count <= 0) return -2;
    rat parent_rotations[GEO_MAX_OPS][3][3];
    for (int i = 0; i < parent_count; i++) rotation_of(parent_sg, (int)parent[i][4], parent_rotations[i]);
    for (int i = 0; i < child_count; i++) {
        rat rotation[3][3], step[3][3], mapped[3][3];
        rotation_of(child_sg, (int)child[i][4], rotation);
        matmul((const rat (*)[3])rotation, (const rat (*)[3])b, step);
        matmul((const rat (*)[3])inverse, (const rat (*)[3])step, mapped);
        int point_op = -1;
        for (int j = 0; j < parent_count; j++)
            if (same3((const rat (*)[3])mapped, (const rat (*)[3])parent_rotations[j])) point_op = (int)parent[j][4];
        if (point_op < 0) return -10 - (int)child[i][4];
        rat t[3], ot[3], tb[3], result[3];
        for (int axis = 0; axis < 3; axis++) t[axis] = rat_make(child[i][axis], child[i][3]);
        rat_row_mul(o, (const rat (*)[3])mapped, ot);
        rat_row_mul(t, (const rat (*)[3])b, tb);
        i128 den = 1;
        for (int axis = 0; axis < 3; axis++) {
            result[axis] = rat_add(rat_sub(o[axis], ot[axis]), tb[axis]);
            den = den / rat_gcd(den, result[axis].d) * result[axis].d;
        }
        for (int axis = 0; axis < 3; axis++) out[i][axis] = (int64_t)(result[axis].n * (den / result[axis].d));
        out[i][3] = (int64_t)den;
        out[i][4] = point_op;
    }
    if (rat_overflow || geo_take_error()) return -2;
    return child_count;
}
