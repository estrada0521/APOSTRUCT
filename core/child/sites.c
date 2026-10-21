/* The child Wyckoff sites of one parent Wyckoff orbit (get_new_wyckoff_):
 * the parent branches split by the subgroup basis/origin into child orbits
 * with their Formula15 rows (child/wyckoff_split.c), each occurrence evaluated at
 * the parent parameters and expanded by the child centering into the atoms
 * of the selected child cell. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"

int wyckoff_split_provenance(int branch_count, const int64_t *parent_branches, const int64_t *basis,
                             const int64_t *origin, int kernel_count, const int64_t *kernel,
                             int op_count, const int64_t *ops, const int64_t *rotations,
                             const int64_t *rotation_dens, int row_count, const int64_t *row_info,
                             const int64_t *row_records, int point_op_code_count, const int64_t *point_op_codes,
                             int point_group_code_count, const int64_t *point_group_codes, int capacity,
                             int64_t *out_rows, int64_t *out_formula, int64_t *out_orbits, int *group_count);

enum {
    CS_OK = 0,
    CS_GEOMETRY = 1,
    CS_RANGE = 2,
    CS_PARENT_ROW = 3,
    CS_NO_ROWS = 4,
    CS_PARTITION = 5,
    CS_MEMORY = 6,
    CS_CAPACITY = 7,
    CS_ZERO_DENOMINATOR = 8,
    CS_SPLIT = 100,
};

enum {
    SITE_WIDTH = 5,
    ATOM_WIDTH = 37,
};

static void canonical_record(const rat *v, int64_t *out)
{
    i128 den = 1;
    for (int axis = 0; axis < 3; axis++) den = rat_mul_i(den / rat_gcd(den, v[axis].d), v[axis].d);
    i128 n[3], g = den;
    for (int axis = 0; axis < 3; axis++) {
        n[axis] = rat_mul_i(v[axis].n, den / v[axis].d);
        g = rat_gcd(g, n[axis]);
    }
    if (g > 1) {
        for (int axis = 0; axis < 3; axis++) n[axis] /= g;
        den /= g;
    }
    for (int axis = 0; axis < 3; axis++) {
        if (n[axis] > INT32_MAX || n[axis] < INT32_MIN) rat_overflow = 1;
        out[axis] = (int64_t)n[axis];
    }
    if (den > INT32_MAX) rat_overflow = 1;
    out[3] = (int64_t)den;
}

static void record_values(const int64_t *record, rat *out)
{
    for (int axis = 0; axis < 3; axis++) out[axis] = rat_make(record[axis], record[3]);
}

/* parent_params: x, y, z as (num, den).  sites_out per child site: child row
 * offset (1-based), multiplicity, first atom, atom count, representative
 * atom; atoms_out per atom: flat index, kernel ordinal, kernel fraction[4],
 * parent branch ordinal, parent coset record[5], centering ordinal,
 * parent xyz[12], child xyz[12], each coordinate an int128 (num, den) as
 * (high, low) words; parent_params likewise. */
int wyckoff_child_sites(int parent_sg, int child_sg, int parent_row_id, const int64_t *basis, const int64_t *origin,
                        const int64_t *parent_params, int site_capacity, int atom_capacity, int64_t *sites_out,
                        int64_t *atoms_out, int *site_count_out, int *atom_count_out)
{
    if (geo_load()) return CS_GEOMETRY;
    geo_take_error();
    rat_overflow = 0;
    int status = CS_OK;
    int64_t *kernel = 0, *flat_rows = 0, *flat_formula = 0, *flat_orbits = 0;
    rat (*used_points)[3] = 0;

    int start, row_count;
    geo_wyckoff_row_range(parent_sg, &start, &row_count);
    if (parent_row_id < start || parent_row_id >= start + row_count) return CS_PARENT_ROW;
    if (!origin[3]) return CS_ZERO_DENOMINATOR;

    int parent_site_pg = geo_wyckoff_site_pg(parent_row_id);
    int64_t site_records[GEO_MAX_OPS][5], cosets[GEO_MAX_OPS][5];
    int site_count = geo_site_records(parent_sg, parent_row_id, parent_site_pg, site_records);
    int branch_count = geo_cosets(parent_sg, (const int64_t (*)[5])site_records, site_count, cosets);
    if (geo_take_error() || !branch_count) return CS_GEOMETRY;
    rat parent_vectors[4][3];
    geo_wyckoff_vectors(parent_row_id, parent_vectors);
    int64_t branches[GEO_MAX_OPS][16];
    for (int b = 0; b < branch_count; b++) {
        for (int piece = 0; piece < 4; piece++) {
            rat v[3];
            if (b == 0)
                memcpy(v, parent_vectors[piece], sizeof v);
            else
                geo_vrot(parent_sg, (int)cosets[b][4], parent_vectors[piece], v);
            canonical_record(v, branches[b] + 4 * piece);
        }
        if (b > 0) {
            rat translation[3], first[3], sum[3];
            record_values(cosets[b], translation);
            record_values(branches[b], first);
            for (int axis = 0; axis < 3; axis++) sum[axis] = rat_add(translation[axis], first[axis]);
            canonical_record(sum, branches[b]);
        }
    }

    int64_t kernel_count = geo_kernel_fractions(basis, &kernel);
    if (kernel_count <= 0) return CS_RANGE;

    int64_t child_ops[GEO_MAX_OPS][5], rotations[GEO_MAX_OPS][9], rotation_dens[GEO_MAX_OPS];
    int op_count = geo_sg_records(child_sg, child_ops);
    for (int op = 0; op < op_count; op++) {
        rat columns[3][3];
        for (int col = 0; col < 3; col++) {
            rat unit[3] = {rat_int(col == 0), rat_int(col == 1), rat_int(col == 2)};
            geo_vrot(child_sg, (int)child_ops[op][4], unit, columns[col]);
        }
        i128 den = 1;
        for (int col = 0; col < 3; col++)
            for (int axis = 0; axis < 3; axis++) den = rat_mul_i(den / rat_gcd(den, columns[col][axis].d), columns[col][axis].d);
        rotation_dens[op] = (int64_t)den;
        for (int axis = 0; axis < 3; axis++)
            for (int col = 0; col < 3; col++)
                rotations[op][axis * 3 + col] = (int64_t)(columns[col][axis].n * (den / columns[col][axis].d));
    }
    int child_start, child_rows;
    geo_wyckoff_row_range(child_sg, &child_start, &child_rows);
    int64_t *row_info = malloc(sizeof(int64_t) * 2 * (size_t)(child_rows + 1));
    int64_t *row_records = malloc(sizeof(int64_t) * 16 * (size_t)(child_rows + 1));
    if (!row_info || !row_records) {
        status = CS_MEMORY;
        goto done;
    }
    for (int row = 0; row < child_rows; row++) {
        row_info[2 * row] = geo_wyckoff_site_pg(child_start + row);
        row_info[2 * row + 1] = row;
        rat vectors[4][3];
        geo_wyckoff_vectors(child_start + row, vectors);
        for (int piece = 0; piece < 4; piece++) canonical_record(vectors[piece], row_records + 16 * row + 4 * piece);
    }
    if (geo_take_error() || rat_overflow) {
        status = CS_GEOMETRY;
        goto done;
    }

    int flat_count = (int)kernel_count * branch_count;
    flat_rows = malloc(sizeof(int64_t) * 5 * (size_t)flat_count);
    flat_formula = malloc(sizeof(int64_t) * 15 * (size_t)flat_count);
    flat_orbits = malloc(sizeof(int64_t) * (size_t)flat_count);
    if (!flat_rows || !flat_formula || !flat_orbits) {
        status = CS_MEMORY;
        goto done;
    }
    int64_t origin_record[4] = {origin[0], origin[1], origin[2], origin[3]};
    if (origin_record[3] < 0)
        for (int i = 0; i < 4; i++) origin_record[i] = -origin_record[i];
    i128 g = rat_gcd(rat_gcd(origin_record[0], origin_record[1]), rat_gcd(origin_record[2], origin_record[3]));
    if (g > 1)
        for (int i = 0; i < 4; i++) origin_record[i] /= (int64_t)g;
    int64_t op_code_count, group_code_count;
    const int64_t *op_codes = geo_point_op_codes(&op_code_count);
    const int64_t *group_codes = geo_point_group_codes(&group_code_count);
    int group_count = 0;
    int split = wyckoff_split_provenance(branch_count, (const int64_t *)branches, basis, origin_record, (int)kernel_count,
                                         kernel, op_count, (const int64_t *)child_ops, (const int64_t *)rotations,
                                         rotation_dens, child_rows, row_info, row_records, (int)op_code_count, op_codes,
                                         (int)group_code_count, group_codes, flat_count, flat_rows, flat_formula,
                                         flat_orbits, &group_count);
    if (split) {
        status = CS_SPLIT + split;
        goto done;
    }
    if (!group_count) {
        status = CS_NO_ROWS;
        goto done;
    }

    rat params[3];
    for (int axis = 0; axis < 3; axis++)
        params[axis] = rat_make(read_i128(parent_params + 4 * axis), read_i128(parent_params + 4 * axis + 2));
    rat result_basis[3][3], result_inverse[3][3], result_origin[3];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) result_basis[row][col] = rat_int(basis[3 * row + col]);
    if (origin[3] <= 0 || rat_inverse3((const rat (*)[3])result_basis, result_inverse)) {
        status = CS_PARTITION;
        goto done;
    }
    for (int axis = 0; axis < 3; axis++) result_origin[axis] = rat_make(origin[axis], origin[3]);
    rat child_to_cinter[3][3], child_origin[3], child_to_pml[3][3];
    if (geo_setting_to_cinter(child_sg, 1, 0, child_to_cinter, child_origin) ||
        rat_inverse3((const rat (*)[3])child_to_cinter, child_to_pml)) {
        status = CS_PARTITION;
        goto done;
    }
    rat centering[4][3];
    int centering_count = geo_centering(child_sg, centering);

    used_points = malloc(sizeof(rat[3]) * (size_t)(atom_capacity + 1));
    if (!used_points) {
        status = CS_MEMORY;
        goto done;
    }
    int atoms = 0, sites = 0;
    for (int group = 0; group < group_count && !status; group++) {
        const int64_t *row = flat_rows + 5 * group;
        int64_t offset = row[0];
        if (offset < 1 || offset > child_rows) {
            status = CS_NO_ROWS;
            break;
        }
        int64_t multiplicity = geo_multiplicity(child_sg, geo_wyckoff_site_pg(child_start + (int)offset - 1));
        int64_t representative_flat = row[1], orbit_start = row[2], orbit_size = row[3];
        if (multiplicity <= 0 || !orbit_size || multiplicity % orbit_size) {
            status = CS_PARTITION;
            break;
        }
        int64_t expansion = multiplicity / orbit_size;
        int translation_count = expansion == 1 ? 1 : centering_count;
        if (expansion != 1 && centering_count != expansion) {
            status = CS_PARTITION;
            break;
        }
        rat shifts[4][3];
        for (int t = 0; t < translation_count; t++) {
            rat pml[3];
            rat_row_mul(centering[t], (const rat (*)[3])child_to_pml, pml);
            rat_row_mul(pml, (const rat (*)[3])result_basis, shifts[t]);
        }
        if (sites >= site_capacity) {
            status = CS_CAPACITY;
            break;
        }
        int64_t *site = sites_out + SITE_WIDTH * sites++;
        site[0] = offset;
        site[1] = multiplicity;
        site[2] = atoms;
        site[4] = -1;
        int site_first = atoms;
        for (int64_t o = 0; o < orbit_size && !status; o++) {
            int64_t flat = flat_orbits[orbit_start + o];
            int64_t kernel_ordinal = flat / branch_count, branch_ordinal = flat % branch_count;
            rat branch[4][3], kernel_value[3], parent_point[3];
            for (int piece = 0; piece < 4; piece++) record_values(branches[branch_ordinal] + 4 * piece, branch[piece]);
            record_values(kernel + 4 * kernel_ordinal, kernel_value);
            for (int axis = 0; axis < 3; axis++) {
                rat value = rat_add(branch[0][axis], kernel_value[axis]);
                for (int i = 0; i < 3; i++) value = rat_add(value, rat_mul(params[i], branch[i + 1][axis]));
                parent_point[axis] = value;
            }
            rat relative[3], child_pml[3], child_cinter[3];
            for (int axis = 0; axis < 3; axis++) relative[axis] = rat_sub(parent_point[axis], result_origin[axis]);
            rat_row_mul(relative, (const rat (*)[3])result_inverse, child_pml);
            rat_row_mul(child_pml, (const rat (*)[3])child_to_cinter, child_cinter);
            for (int axis = 0; axis < 3; axis++) child_cinter[axis] = rat_add(child_cinter[axis], child_origin[axis]);
            for (int t = 0; t < translation_count; t++) {
                rat xyz[3], parent_xyz[3];
                for (int axis = 0; axis < 3; axis++) {
                    xyz[axis] = rat_mod1(rat_add(child_cinter[axis], centering[t][axis]));
                    parent_xyz[axis] = rat_add(parent_point[axis], shifts[t][axis]);
                }
                for (int a = 0; a < atoms && !status; a++)
                    if (rat_eq(used_points[a][0], xyz[0]) && rat_eq(used_points[a][1], xyz[1]) && rat_eq(used_points[a][2], xyz[2]))
                        status = CS_PARTITION;
                for (int a = 0; a < atoms && !status; a++) {
                    const int64_t *other = atoms_out + ATOM_WIDTH * a;
                    if (other[0] == flat && other[12] == t) status = CS_PARTITION;
                }
                if (status) break;
                if (atoms >= atom_capacity) {
                    status = CS_CAPACITY;
                    break;
                }
                memcpy(used_points[atoms], xyz, sizeof xyz);
                if (flat == representative_flat && t == 0) site[4] = atoms - site_first;
                int64_t *atom = atoms_out + ATOM_WIDTH * atoms++;
                atom[0] = flat;
                atom[1] = kernel_ordinal;
                memcpy(atom + 2, kernel + 4 * kernel_ordinal, 4 * sizeof(int64_t));
                atom[6] = branch_ordinal;
                memcpy(atom + 7, cosets[branch_ordinal], 5 * sizeof(int64_t));
                atom[12] = t;
                write_rat3(atom + 13, parent_xyz);
                write_rat3(atom + 25, xyz);
            }
        }
        site[3] = atoms - site_first;
        if (!status && (site[3] != multiplicity || site[4] < 0)) status = CS_PARTITION;
    }
    if (!status && rat_overflow) status = CS_RANGE;
    *site_count_out = sites;
    *atom_count_out = atoms;
done:
    free(kernel);
    free(row_info);
    free(row_records);
    free(flat_rows);
    free(flat_formula);
    free(flat_orbits);
    free(used_points);
    return status;
}
