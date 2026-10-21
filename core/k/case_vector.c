/* The k vector of a mode case and the records find_isotropy_subgroup_
 * builds from it: the newlat_ input records, the order-parameter entry
 * records, and the newlat_ integer superlattice basis. */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "tables.h"

enum {
    KV_TABLES = 1,
    KV_PARAMS = 2,
    KV_POINTER = 3,
    KV_RANGE = 4,
    KV_NEWLAT = 5,
};

static int64_t at(const section_view *view, int64_t index, int *bad)
{
    if (index < 0 || index >= section_count(view)) {
        *bad = 1;
        return 0;
    }
    return section_int(view, index);
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

/* params: the case's k parameters as rationals (int128 words, 4 per value).
 * newlat_out: 4 x 4 newlat_ input records; entry_out: 8 x 4 entry records
 * (entry_count 0 for fixed k). */
int case_k_records(int sg, int kslot, int param_count, const int64_t *params, int64_t *newlat_out,
                   int64_t *entry_out, int *entry_count)
{
    static section_view little_k, k_dim, star, conv2ml;
    static int ready;
    if (!ready) {
        if (geo_load() || section_find("little.little_k", &little_k) || section_find("little.little_k_dim", &k_dim) ||
            section_find("little.little_k_star", &star) ||
            section_find("little.little_k_star_conv2ml_pointer", &conv2ml))
            return KV_TABLES;
        ready = 1;
    }
    rat_overflow = 0;
    int bad = 0;
    int64_t slot = (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1;
    int dim = (int)at(&k_dim, slot, &bad);
    if (bad) return KV_TABLES;
    if (param_count < dim) return KV_PARAMS;
    rat vectors[4][3], p[3];
    for (int v = 0; v < 4; v++) {
        int64_t den = at(&little_k, slot * 16 + 4 * v + 3, &bad);
        for (int axis = 0; axis < 3; axis++)
            vectors[v][axis] = den ? rat_make(at(&little_k, slot * 16 + 4 * v + axis, &bad), den) : rat_int(0);
    }
    for (int i = 0; i < dim && i < 3; i++) p[i] = rat_make(read_i128(params + 4 * i), read_i128(params + 4 * i + 2));
    rat evaluated[3];
    for (int axis = 0; axis < 3; axis++) {
        evaluated[axis] = vectors[0][axis];
        for (int i = 0; i < dim; i++) evaluated[axis] = rat_add(evaluated[axis], rat_mul(p[i], vectors[i + 1][axis]));
    }
    vector_record(evaluated, newlat_out);
    for (int i = 0; i < 3; i++) {
        int64_t *row = newlat_out + 4 * (i + 1);
        if (i < dim) {
            rat contribution[3];
            for (int axis = 0; axis < 3; axis++) contribution[axis] = rat_mul(p[i], vectors[i + 1][axis]);
            vector_record(contribution, row);
        } else {
            row[0] = row[1] = row[2] = 0;
            row[3] = 1;
        }
    }
    *entry_count = 0;
    if (dim > 0) {
        int64_t pointer = at(&conv2ml, (int64_t)(sg - 1) * 27 + kslot - 1, &bad);
        if (bad || pointer <= 0) return KV_POINTER;
        if (16 * pointer > section_count(&star)) return KV_POINTER;
        rat transformed[4] = {rat_int(0), rat_int(0), rat_int(0), rat_int(0)};
        for (int i = 0; i < dim; i++) {
            const int64_t *record = newlat_out + 4 * (i + 1);
            rat vector[4] = {rat_make(record[0], record[3]), rat_make(record[1], record[3]), rat_make(record[2], record[3]),
                             rat_int(0)};
            for (int col = 0; col < 4; col++) {
                rat total = rat_int(0);
                for (int row = 0; row < 4; row++)
                    total = rat_add(total, rat_mul(vector[row], rat_int(at(&star, 16 * (pointer - 1) + 4 * row + col, &bad))));
                transformed[col] = rat_add(transformed[col], total);
            }
        }
        rat scale = rat_zero(transformed[3]) ? rat_int(1) : transformed[3], first[3];
        for (int axis = 0; axis < 3; axis++) first[axis] = rat_div(transformed[axis], scale);
        memset(entry_out, 0, sizeof(int64_t) * 32);
        vector_record(first, entry_out);
        memcpy(entry_out + 4, newlat_out, sizeof(int64_t) * 4);
        memcpy(entry_out + 12, newlat_out, sizeof(int64_t) * 4);
        for (int i = 0; i < dim; i++) memcpy(entry_out + 16 + 4 * i, newlat_out + 4 * (i + 1), sizeof(int64_t) * 4);
        *entry_count = 8;
    }
    if (bad) return KV_TABLES;
    return rat_overflow ? KV_RANGE : 0;
}

static int64_t max_abs(const int64_t *row)
{
    int64_t m = 0;
    for (int axis = 0; axis < 3; axis++) {
        int64_t v = row[axis] < 0 ? -row[axis] : row[axis];
        if (v > m) m = v;
    }
    return m;
}

/* newlat_: rows with integral dot products against the first mode_count
 * records, then three passes of max-norm reduction. */
int newlat_basis(int mode_count, int record_count, const int64_t *records, int64_t *out)
{
    if (mode_count <= 0 || record_count < mode_count) return KV_NEWLAT;
    int64_t max_den = records[3];
    if (max_den <= 0) return KV_NEWLAT;
    for (int r = 0; r < mode_count; r++)
        if (!records[4 * r + 3]) return KV_NEWLAT;
    int64_t rows[3][3];
    for (int pivot_index = 0; pivot_index < 3; pivot_index++) {
        int found = 0;
        for (int64_t pivot = 1; pivot <= max_den && !found; pivot++) {
            int64_t combos = 1;
            for (int i = 0; i < pivot_index; i++) combos *= max_den + 1;
            for (int64_t serial = 0; serial < combos && !found; serial++) {
                int64_t candidate[3] = {0, 0, 0}, value = serial;
                for (int axis = 0; axis < pivot_index; axis++) {
                    candidate[axis] = value % (max_den + 1);
                    value /= max_den + 1;
                }
                candidate[pivot_index] = pivot;
                int ok = 1;
                for (int r = 0; r < mode_count && ok; r++) {
                    i128 dot = 0;
                    for (int axis = 0; axis <= pivot_index; axis++) dot += (i128)candidate[axis] * records[4 * r + axis];
                    ok = dot % records[4 * r + 3] == 0;
                }
                if (ok) {
                    memcpy(rows[pivot_index], candidate, sizeof candidate);
                    found = 1;
                }
            }
        }
        if (!found) return KV_NEWLAT;
    }
    for (int pass = 0; pass < 3; pass++)
        for (int target = 0; target < 3; target++)
            for (int other = 0; other < 3; other++) {
                if (other == target) continue;
                for (int sign = -1; sign <= 1; sign += 2) {
                    int64_t candidate[3];
                    for (int axis = 0; axis < 3; axis++) candidate[axis] = rows[target][axis] + sign * rows[other][axis];
                    if (max_abs(candidate) < max_abs(rows[target])) memcpy(rows[target], candidate, sizeof candidate);
                }
            }
    for (int i = 0; i < 9; i++) out[i] = rows[i / 3][i % 3];
    return 0;
}

/* display_k_vector_from_case: the first Miller-Love star arm of the case's
 * SG/k slot (parametric arms with the case parameters substituted), in the
 * reciprocal basis of the default International setting.  out: three
 * rationals as int128 words. */
int case_display_kvector(int sg, int kslot, int param_count, const int64_t *params, int64_t *out)
{
    static section_view little_sg, little_kslot, star, special, star_count, star_pointer, k_dim;
    static int ready;
    if (!ready) {
        if (geo_load() || section_find("little.little_irr_space_group", &little_sg) ||
            section_find("little.little_irr_k", &little_kslot) || section_find("little.little_k_star", &star) ||
            section_find("little.little_k_star_special", &special) ||
            section_find("little.little_k_star_count", &star_count) ||
            section_find("little.little_k_star_ml_pointer", &star_pointer) || section_find("little.little_k_dim", &k_dim))
            return KV_TABLES;
        ready = 1;
    }
    rat_overflow = 0;
    int bad = 0, gid = 0;
    int64_t total = section_count(&little_sg);
    for (int64_t i = 0; i < total && !gid; i++)
        if (section_int(&little_sg, i) == sg && section_int(&little_kslot, i) == kslot) gid = (int)i + 1;
    if (!gid) return KV_POINTER;
    int64_t slot = (int64_t)(sg - 1) * 27 + kslot - 1;
    int64_t count = at(&star_count, slot, &bad), pointer = at(&star_pointer, slot, &bad);
    int dim = (int)at(&k_dim, (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1, &bad);
    if (bad || count <= 0 || pointer <= 0) return KV_POINTER;
    rat arm[3];
    if (dim > 0) {
        rat pieces[4][3], p[3] = {rat_int(0), rat_int(0), rat_int(0)};
        for (int i = 0; i < param_count && i < 3; i++) p[i] = rat_make(read_i128(params + 4 * i), read_i128(params + 4 * i + 2));
        int64_t start = 16 * (pointer - 1);
        if (start + 16 > section_count(&star)) return KV_POINTER;
        for (int piece = 0; piece < 4; piece++) {
            int64_t den = at(&star, start + 4 * piece + 3, &bad);
            for (int axis = 0; axis < 3; axis++)
                pieces[piece][axis] = den ? rat_make(at(&star, start + 4 * piece + axis, &bad), den) : rat_int(0);
        }
        for (int axis = 0; axis < 3; axis++) {
            arm[axis] = pieces[0][axis];
            for (int i = 0; i < 3; i++) arm[axis] = rat_add(arm[axis], rat_mul(p[i], pieces[i + 1][axis]));
        }
    } else {
        int64_t start = 4 * (pointer - 1), den = at(&special, start + 3, &bad);
        if (!den) return KV_POINTER;
        for (int axis = 0; axis < 3; axis++) arm[axis] = rat_make(at(&special, start + axis, &bad), den);
    }
    rat matrix[3][3], origin[3], inverse[3][3], visible[3];
    geo_setting_to_cinter(sg, 1, 0, matrix, origin);
    if (rat_inverse3((const rat (*)[3])matrix, inverse)) return KV_TABLES;
    for (int col = 0; col < 3; col++) {
        visible[col] = rat_int(0);
        for (int row = 0; row < 3; row++) visible[col] = rat_add(visible[col], rat_mul(arm[row], inverse[col][row]));
    }
    write_rat3(out, visible);
    if (bad || geo_take_error()) return KV_TABLES;
    return rat_overflow ? KV_RANGE : 0;
}

/* the orbit key of a k vector: the least
 * ((x mod 1, den), (y mod 1, den), (z mod 1, den)) over the images of a
 * reciprocal PML vector under the inverses of the parent point operations.
 * vector: three rationals as int128 words; out: (numerator, denominator) x 3. */
static int orbit_key(int sg, const rat *source, int64_t *out);

/* reciprocal_orbit_key for a vector already held as rationals. */
int reciprocal_orbit_key_of(int sg, const rat *vector, int64_t *out)
{
    if (geo_load()) return KV_TABLES;
    rat_overflow = 0;
    return orbit_key(sg, vector, out);
}

static int orbit_key(int sg, const rat *source, int64_t *out)
{
    rat best[3];
    int64_t records[GEO_MAX_OPS][5];
    int count = geo_sg_records(sg, records), seen[145] = {0}, found = 0;
    for (int r = 0; r < count; r++) {
        int op = (int)records[r][4];
        if (op < 1 || op > 144 || seen[op]) continue;
        seen[op] = 1;
        rat direct[3][3], inverse[3][3], image[3];
        for (int col = 0; col < 3; col++) {
            rat unit[3] = {rat_int(col == 0), rat_int(col == 1), rat_int(col == 2)}, rotated[3];
            geo_vrot(sg, op, unit, rotated);
            for (int row = 0; row < 3; row++) direct[row][col] = rotated[row];
        }
        if (rat_inverse3((const rat (*)[3])direct, inverse)) return KV_RANGE;
        rat_row_mul(source, (const rat (*)[3])inverse, image);
        for (int axis = 0; axis < 3; axis++) image[axis] = rat_mod1(image[axis]);
        int less = !found;
        for (int axis = 0; axis < 3 && !less; axis++) {
            if (image[axis].n != best[axis].n) {
                less = image[axis].n < best[axis].n;
                break;
            }
            if (image[axis].d != best[axis].d) {
                less = image[axis].d < best[axis].d;
                break;
            }
        }
        if (less) memcpy(best, image, sizeof best);
        found = 1;
    }
    if (!found)
        for (int axis = 0; axis < 3; axis++) best[axis] = rat_mod1(source[axis]);
    for (int axis = 0; axis < 3; axis++) {
        if (best[axis].n > INT64_MAX || best[axis].d > INT64_MAX) return KV_RANGE;
        out[2 * axis] = (int64_t)best[axis].n;
        out[2 * axis + 1] = (int64_t)best[axis].d;
    }
    if (geo_take_error()) return KV_TABLES;
    return rat_overflow ? KV_RANGE : 0;
}
