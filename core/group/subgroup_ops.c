/* Operation-record algebra of embedded subgroups (coupled_subgroups2_,
 * same_subgroup_, map_subgroup2_, generate_domains2_ and its magnetic
 * variant): kernel cosets of space operations, lattice membership, subgroup
 * intersection and the domain coset scan.  group > 0 selects the records of
 * that magnetic group in place of the ordinary space group's. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"

enum {
    SO_KEY = 1,
    SO_VALUE = 2,
    SO_RUNTIME = 3,
    SO_RANGE = 4,
    SO_CAPACITY = 5,
    SO_ZERO = 6,
    SO_TABLES = 7,
};

#define SO_MAX_RECORDS 192

static int so_bad;

static void so_fail(int status)
{
    if (!so_bad) so_bad = status;
}

static int records_of(int sg, int group, int64_t (*out)[5])
{
    if (!group) return geo_sg_records(sg, out);
    int count = magnetic_records_iso(group, SO_MAX_RECORDS, out[0]);
    if (count < 0) so_fail(SO_TABLES);
    return count < 0 ? 0 : count;
}

static void values_of(const int64_t *record, rat *out)
{
    if (!record[3]) {
        so_fail(SO_ZERO);
        out[0] = out[1] = out[2] = rat_int(0);
        return;
    }
    for (int axis = 0; axis < 3; axis++) out[axis] = rat_make(record[axis], record[3]);
}

/* (x, y, z, den) over the least common denominator. */
static void vector_record(const rat *v, int64_t *out)
{
    i128 den = 1;
    for (int axis = 0; axis < 3; axis++) den = rat_mul_i(den / rat_gcd(den, v[axis].d), v[axis].d);
    for (int axis = 0; axis < 3; axis++) {
        i128 n = rat_mul_i(v[axis].n, den / v[axis].d);
        if (n > INT64_MAX || n < INT64_MIN) rat_overflow = 1;
        out[axis] = (int64_t)n;
    }
    if (den > INT64_MAX) rat_overflow = 1;
    out[3] = (int64_t)den;
}

/* SourceData._vadd_fraction_operation */
static void vadd(const rat *fraction, const int64_t *operation, int64_t *out)
{
    rat t[3], sum[3];
    values_of(operation, t);
    for (int axis = 0; axis < 3; axis++) sum[axis] = rat_add(fraction[axis], t[axis]);
    vector_record(sum, out);
    out[4] = operation[4];
}

static int64_t kernel_of(const int64_t *basis, int64_t **out)
{
    int64_t count = geo_kernel_fractions(basis, out);
    if (count <= 0) so_fail(count == -GEO_RANGE ? SO_VALUE : SO_TABLES);
    return count < 0 ? 0 : count;
}

static void kernel_values(const int64_t *kernel, int64_t index, rat *out)
{
    for (int axis = 0; axis < 3; axis++) out[axis] = rat_make(kernel[4 * index + axis], kernel[4 * index + 3]);
}

static int finish(int64_t count)
{
    int status = so_bad;
    if (!status && geo_take_error()) status = SO_TABLES;
    if (!status && rat_overflow) status = SO_RANGE;
    return status ? -status : (int)count;
}

static int start(void)
{
    so_bad = 0;
    rat_overflow = 0;
    if (geo_load() || mag_tables()) return SO_TABLES;
    geo_take_error();
    return 0;
}

/* kernel_operation_records: every space operation times every kernel coset,
 * the coset rotated by the operation (magnetic records: unrotated). */
int subgroup_kernel_records(int sg, int group, const int64_t *basis, int64_t capacity, int64_t *out)
{
    if (start()) return -SO_TABLES;
    int64_t records[SO_MAX_RECORDS][5], *kernel = 0;
    int record_count = records_of(sg, group, records);
    int64_t kernel_count = kernel_of(basis, &kernel), count = 0;
    for (int r = 0; r < record_count && !so_bad; r++)
        for (int64_t f = 0; f < kernel_count; f++) {
            if (count == capacity) {
                so_fail(SO_CAPACITY);
                break;
            }
            rat fraction[3], rotated[3];
            kernel_values(kernel, f, fraction);
            if (group)
                memcpy(rotated, fraction, sizeof rotated);
            else
                geo_vrot(sg, (int)records[r][4], fraction, rotated);
            vadd(rotated, records[r], out + 5 * count++);
        }
    free(kernel);
    return finish(count);
}

typedef struct {
    i128 adj[9], det;
} lattice_test;

static void lattice_test_of(const int64_t *basis, lattice_test *t)
{
    i128 b[9];
    for (int i = 0; i < 9; i++) b[i] = basis[i];
    i128 adj[9] = {
        b[4] * b[8] - b[5] * b[7], b[2] * b[7] - b[1] * b[8], b[1] * b[5] - b[2] * b[4],
        b[5] * b[6] - b[3] * b[8], b[0] * b[8] - b[2] * b[6], b[2] * b[3] - b[0] * b[5],
        b[3] * b[7] - b[4] * b[6], b[1] * b[6] - b[0] * b[7], b[0] * b[4] - b[1] * b[3],
    };
    memcpy(t->adj, adj, sizeof adj);
    t->det = b[0] * adj[0] + b[1] * adj[3] + b[2] * adj[6];
}

static int in_tested_lattice(const int64_t *record, const lattice_test *t)
{
    const i128 *adj = t->adj;
    i128 det = t->det;
    if (!det) return 0;
    i128 scale = rat_mul_i(record[3] ? record[3] : 1, det);
    for (int col = 0; col < 3; col++) {
        i128 value = rat_add_i(rat_add_i(rat_mul_i(record[0], adj[col]), rat_mul_i(record[1], adj[3 + col])),
                               rat_mul_i(record[2], adj[6 + col]));
        if (value % scale) return 0;
    }
    return 1;
}

static int in_lattice(const int64_t *record, const int64_t *basis)
{
    i128 b[9];
    for (int i = 0; i < 9; i++) b[i] = basis[i];
    i128 adj[9] = {
        b[4] * b[8] - b[5] * b[7], b[2] * b[7] - b[1] * b[8], b[1] * b[5] - b[2] * b[4],
        b[5] * b[6] - b[3] * b[8], b[0] * b[8] - b[2] * b[6], b[2] * b[3] - b[0] * b[5],
        b[3] * b[7] - b[4] * b[6], b[1] * b[6] - b[0] * b[7], b[0] * b[4] - b[1] * b[3],
    };
    i128 det = b[0] * adj[0] + b[1] * adj[3] + b[2] * adj[6];
    if (!det) return 0;
    i128 scale = rat_mul_i(record[3] ? record[3] : 1, det);
    for (int col = 0; col < 3; col++) {
        i128 value = rat_add_i(rat_add_i(rat_mul_i(record[0], adj[col]), rat_mul_i(record[1], adj[3 + col])),
                               rat_mul_i(record[2], adj[6 + col]));
        if (value % scale) return 0;
    }
    return 1;
}

/* coupled._operation_in_subgroup */
static int in_subgroup(const int64_t *candidate, const int64_t *basis, int count, const int64_t *operations)
{
    const int64_t *representative = 0;
    for (int i = 0; i < count && !representative; i++)
        if (operations[5 * i + 4] == candidate[4]) representative = operations + 5 * i;
    if (!representative) return 0;
    rat c[3], r[3], delta[3];
    values_of(candidate, c);
    values_of(representative, r);
    for (int axis = 0; axis < 3; axis++) delta[axis] = rat_sub(c[axis], r[axis]);
    int64_t record[4];
    vector_record(delta, record);
    return in_lattice(record, basis);
}

/* One subgroup of subgroup_intersection: its lattice test and, per point
 * operation, the first operation record carrying it (values made on first
 * use, as in_subgroup would). */
typedef struct {
    lattice_test lattice;
    int count;
    const int64_t *operations;
    int first[256];
    int ready[256];
    rat values[256][3];
} subgroup_test;

static int in_tested_subgroup(const int64_t *candidate, const rat *candidate_values, subgroup_test *t)
{
    int64_t op = candidate[4];
    int index = -1;
    if (op >= 0 && op < 256) {
        index = t->first[op];
    } else {
        for (int i = 0; i < t->count && index < 0; i++)
            if (t->operations[5 * i + 4] == op) index = i;
    }
    if (index < 0) return 0;
    const rat *r;
    rat local[3];
    if (op >= 0 && op < 256) {
        if (!t->ready[op]) {
            values_of(t->operations + 5 * index, t->values[op]);
            t->ready[op] = 1;
        }
        r = t->values[op];
    } else {
        values_of(t->operations + 5 * index, local);
        r = local;
    }
    rat delta[3];
    for (int axis = 0; axis < 3; axis++) delta[axis] = rat_sub(candidate_values[axis], r[axis]);
    int64_t record[4];
    vector_record(delta, record);
    return in_tested_lattice(record, &t->lattice);
}

/* The first kernel coset of each space operation lying in every subgroup. */
int subgroup_intersection(int sg, int group, const int64_t *basis, int subgroup_count, const int64_t *bases,
                          const int *operation_counts, const int64_t *operations, int64_t capacity, int64_t *out)
{
    if (start()) return -SO_TABLES;
    int64_t records[SO_MAX_RECORDS][5], *kernel = 0;
    int record_count = records_of(sg, group, records);
    int64_t kernel_count = kernel_of(basis, &kernel), count = 0;
    subgroup_test *tests = malloc(sizeof *tests * ((size_t)subgroup_count + 1));
    if (!tests) {
        free(kernel);
        return -SO_TABLES;
    }
    const int64_t *ops = operations;
    for (int s = 0; s < subgroup_count; s++) {
        subgroup_test *t = tests + s;
        lattice_test_of(bases + 9 * s, &t->lattice);
        t->count = operation_counts[s];
        t->operations = ops;
        for (int i = 0; i < 256; i++) t->first[i] = -1, t->ready[i] = 0;
        for (int i = t->count - 1; i >= 0; i--) {
            int64_t op = ops[5 * i + 4];
            if (op >= 0 && op < 256) t->first[op] = i;
        }
        ops += 5 * (size_t)operation_counts[s];
    }
    for (int r = 0; r < record_count && !so_bad; r++) {
        int64_t point = records[r][4];
        int possible = 1;
        for (int s = 0; s < subgroup_count && possible; s++) {
            if (point >= 0 && point < 256) {
                possible = tests[s].first[point] >= 0;
            } else {
                possible = 0;
                for (int i = 0; i < tests[s].count && !possible; i++) possible = tests[s].operations[5 * i + 4] == point;
            }
        }
        if (!possible) continue;
        for (int64_t f = 0; f < kernel_count; f++) {
            rat fraction[3], shifted_values[3];
            int64_t shifted[5];
            kernel_values(kernel, f, fraction);
            vadd(fraction, records[r], shifted);
            int included = 1, valued = 0;
            for (int s = 0; s < subgroup_count && included; s++) {
                int64_t op = shifted[4];
                int has = 0;
                if (op >= 0 && op < 256)
                    has = tests[s].first[op] >= 0;
                else
                    for (int i = 0; i < tests[s].count && !has; i++) has = tests[s].operations[5 * i + 4] == op;
                if (!has) {
                    included = 0;
                    break;
                }
                if (!valued) {
                    values_of(shifted, shifted_values);
                    valued = 1;
                }
                included = in_tested_subgroup(shifted, shifted_values, tests + s);
            }
            if (!included) continue;
            if (count == capacity) {
                so_fail(SO_CAPACITY);
                break;
            }
            memcpy(out + 5 * count++, shifted, sizeof shifted);
            break;
        }
    }
    free(tests);
    free(kernel);
    return finish(count);
}

static int same_point_set(int left_count, const int64_t *left, int right_count, const int64_t *right)
{
    for (int pass = 0; pass < 2; pass++) {
        const int64_t *a = pass ? right : left, *b = pass ? left : right;
        int a_count = pass ? right_count : left_count, b_count = pass ? left_count : right_count;
        for (int i = 0; i < a_count; i++) {
            int found = 0;
            for (int j = 0; j < b_count && !found; j++) found = a[5 * i + 4] == b[5 * j + 4];
            if (!found) return 0;
        }
    }
    return 1;
}

/* same_subgroup_: equal lattices, equal point sets, mutual membership. */
int same_subgroup(const int64_t *left_basis, int left_count, const int64_t *left, const int64_t *right_basis,
                  int right_count, const int64_t *right)
{
    so_bad = 0;
    rat_overflow = 0;
    int same = 1;
    for (int pass = 0; pass < 2 && same; pass++) {
        const int64_t *source = pass ? right_basis : left_basis, *target = pass ? left_basis : right_basis;
        for (int row = 0; row < 3 && same; row++) {
            int64_t record[4] = {source[3 * row], source[3 * row + 1], source[3 * row + 2], 1};
            same = in_lattice(record, target);
        }
    }
    if (same) same = same_point_set(left_count, left, right_count, right);
    for (int i = 0; i < left_count && same; i++) same = in_subgroup(left + 5 * i, right_basis, right_count, right);
    for (int i = 0; i < right_count && same; i++) same = in_subgroup(right + 5 * i, left_basis, left_count, left);
    if (so_bad || rat_overflow) return -(so_bad ? so_bad : SO_RANGE);
    return same;
}

static void mat_mul(const rat (*a)[3], const rat (*b)[3], rat (*out)[3])
{
    rat t[3][3];
    for (int row = 0; row < 3; row++) rat_row_mul(a[row], b, t[row]);
    memcpy(out, t, sizeof t);
}

/* map_subgroup2_: the child's operations as parent operation records. */
int mapped_subgroup_records(int parent_sg, int child_sg, const int64_t *basis, const int64_t *origin,
                            int64_t capacity, int64_t *out)
{
    if (start()) return -SO_TABLES;
    rat b[3][3], inverse[3][3], o[3], shift[3];
    for (int i = 0; i < 9; i++) b[i / 3][i % 3] = rat_int(basis[i]);
    for (int axis = 0; axis < 3; axis++) o[axis] = rat_make(origin[axis], origin[3] ? origin[3] : 1);
    if (rat_inverse3((const rat (*)[3])b, inverse)) return -SO_VALUE;
    rat_row_mul(o, (const rat (*)[3])inverse, shift);
    for (int axis = 0; axis < 3; axis++) shift[axis] = rat_neg(shift[axis]);
    int64_t parents[SO_MAX_RECORDS][5], children[SO_MAX_RECORDS][5], parent_points[SO_MAX_RECORDS][9];
    int parent_count = geo_sg_records(parent_sg, parents), child_count = geo_sg_records(child_sg, children);
    if (child_count > capacity) return -SO_CAPACITY;
    for (int c = 0; c < child_count; c++) {
        int64_t point_values[9];
        if (geo_pml_point_matrix(child_sg, (int)children[c][4], point_values)) return -SO_KEY;
        rat point[3][3], product[3][3], t[3], moved[3], translation[3];
        for (int i = 0; i < 9; i++) point[i / 3][i % 3] = rat_int(point_values[i]);
        mat_mul((const rat (*)[3])inverse, (const rat (*)[3])point, product);
        mat_mul((const rat (*)[3])product, (const rat (*)[3])b, product);
        values_of(children[c], t);
        rat_row_mul(shift, (const rat (*)[3])point, moved);
        for (int axis = 0; axis < 3; axis++) moved[axis] = rat_add(moved[axis], t[axis]);
        rat_row_mul(moved, (const rat (*)[3])b, translation);
        for (int axis = 0; axis < 3; axis++) translation[axis] = rat_add(translation[axis], o[axis]);
        int64_t integral[9];
        for (int i = 0; i < 9; i++) {
            if (product[i / 3][i % 3].d != 1) return -SO_KEY;
            integral[i] = (int64_t)product[i / 3][i % 3].n;
        }
        int point_op = 0;
        for (int p = 0; p < parent_count; p++) {
            if (geo_pml_point_matrix(parent_sg, (int)parents[p][4], parent_points[p])) return -SO_KEY;
            if (!memcmp(parent_points[p], integral, sizeof integral)) point_op = (int)parents[p][4];
        }
        if (!point_op) return -SO_KEY;
        rat wrapped[3];
        for (int axis = 0; axis < 3; axis++) wrapped[axis] = rat_mod1(translation[axis]);
        int found = -1;
        for (int p = 0; p < parent_count && found < 0; p++) {
            if (parents[p][4] != point_op) continue;
            rat r[3];
            values_of(parents[p], r);
            int same = 1;
            for (int axis = 0; axis < 3 && same; axis++) same = rat_zero(rat_mod1(rat_sub(wrapped[axis], r[axis])));
            if (same) found = p;
        }
        if (found < 0) return -SO_KEY;
        memcpy(out + 5 * c, parents[found], sizeof parents[found]);
    }
    return finish(child_count);
}

static int index_of_point(const int64_t (*records)[5], int count, int64_t op)
{
    int index = 0;
    for (int i = 0; i < count; i++)
        if (records[i][4] == op) index = i + 1;
    return index;
}

/* generate_domains2_ (group 0) or generate_domains2_magnetic_: the domain
 * operations of a subgroup in the Source coset scan order. */
int subgroup_domain_records(int sg, int group, const int64_t *basis, int subgroup_count,
                            const int64_t *subgroup_operations, int64_t capacity, int64_t *out)
{
    if (start()) return -SO_TABLES;
    int64_t parents[SO_MAX_RECORDS][5];
    int parent_count = records_of(sg, group, parents);
    if (so_bad) return finish(0);
    int64_t points[SO_MAX_RECORDS];
    int point_count = 0;
    for (int i = 0; i < subgroup_count; i++) {
        int seen = 0;
        for (int j = 0; j < point_count && !seen; j++) seen = points[j] == subgroup_operations[5 * i + 4];
        if (!seen) points[point_count++] = subgroup_operations[5 * i + 4];
    }
    if (!point_count) return -SO_VALUE;
    int identity_index = 0;
    for (int i = 0; i < parent_count && !identity_index; i++) {
        int64_t m[9];
        if (group) {
            geo_point_op_matrix(mag_nonmagnetic_op((int)parents[i][4]), m);
            if (mag_time_reversed((int)parents[i][4])) continue;
        } else if (geo_pml_point_matrix(sg, (int)parents[i][4], m))
            return -SO_KEY;
        if (m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 && m[4] == 1 && m[5] == 0 && m[6] == 0 && m[7] == 0 &&
            m[8] == 1)
            identity_index = i + 1;
    }
    if (!identity_index) return -SO_RUNTIME;
    char marked[SO_MAX_RECORDS + 1] = {0};
    int cosets[SO_MAX_RECORDS], coset_count = 0;
    cosets[coset_count++] = identity_index;
    marked[identity_index] = 1;
    if (group) {
        for (int j = 0; j < point_count; j++) {
            int index = index_of_point((const int64_t (*)[5])parents, parent_count, points[j]);
            if (!index) return -SO_KEY;
            marked[index] = 1;
        }
    } else
        for (int i = 0; i < parent_count; i++)
            for (int j = 0; j < point_count; j++)
                if (parents[i][4] == points[j]) marked[i + 1] = 1;
    int target = parent_count / point_count;
    while (coset_count < target) {
        int chosen = 0;
        for (int i = 1; i <= parent_count && !chosen; i++)
            if (!marked[i]) chosen = i;
        if (!chosen) return -SO_RUNTIME;
        cosets[coset_count++] = chosen;
        for (int j = 0; j < point_count; j++) {
            int64_t product = group ? mag_point_product((int)parents[chosen - 1][4], (int)points[j])
                                    : geo_point_product((int)points[j], (int)parents[chosen - 1][4]);
            int index = index_of_point((const int64_t (*)[5])parents, parent_count, product);
            if (!index) return -SO_KEY;
            marked[index] = 1;
        }
    }
    int64_t *kernel = 0;
    int64_t kernel_count = kernel_of(basis, &kernel), count = 0;
    for (int64_t f = 0; f < kernel_count && !so_bad; f++)
        for (int c = 0; c < coset_count; c++) {
            if (count == capacity) {
                so_fail(SO_CAPACITY);
                break;
            }
            const int64_t *record = parents[cosets[c] - 1];
            int op = group ? mag_nonmagnetic_op((int)record[4]) : (int)record[4];
            rat fraction[3], rotated[3];
            kernel_values(kernel, f, fraction);
            geo_vrot(sg, op, fraction, rotated);
            vadd(rotated, record, out + 5 * count++);
        }
    free(kernel);
    return finish(count);
}
