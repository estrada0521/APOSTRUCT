/* Space-group, setting and Wyckoff geometry decoded from the Source blob
 * (data_space, data_wyckoff). */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "tables.h"

static struct {
    section_view sg_point_group, point_group_order, generators, lattice, lattice_ml,
        lattice_ml_denom, point_op, inter_choice, settings_inter, label_bc, inter_label_short,
        point_op_code, point_group_code, wyckoff_pointer, wyckoff_count, wyckoff_fract,
        inter_wyckoff_pointer, inter_wyckoff_fract,
        pg_elements_pointer, pg_settings_count, pg_elements, settings,
        group_gen_count, group_gen, inter_number;
    int64_t *op_codes, *group_codes, gen_start[231];
    int ready, order[231];
    short products[72 * 72];
    int64_t (*records)[GEO_MAX_OPS][5];
} G;

int geo_load(void)
{
    if (G.ready) return 0;
    struct {
        const char *name;
        section_view *view;
        int text;
    } items[] = {
        {"group.space_groups.point_group", &G.sg_point_group, 0},
        {"group.point_groups.order", &G.point_group_order, 0},
        {"group.space_groups.generators", &G.generators, 0},
        {"group.space_groups.lattice", &G.lattice, 0},
        {"group.lattices.ml", &G.lattice_ml, 0},
        {"group.lattices.ml_denom", &G.lattice_ml_denom, 0},
        {"group.point_ops.matrix", &G.point_op, 0},
        {"group.space_groups.inter_choice", &G.inter_choice, 0},
        {"settings.inter_settings.inter", &G.settings_inter, 0},
        {"group.space_groups.label_bc", &G.label_bc, 1},
        {"settings.inter_settings.short", &G.inter_label_short, 1},
        {"group.point_ops.code", &G.point_op_code, 0},
        {"group.point_groups.code", &G.point_group_code, 0},
        {"settings.inter_settings.wyckoff_start", &G.inter_wyckoff_pointer, 0},
        {"settings.inter_settings.wyckoff", &G.inter_wyckoff_fract, 0},
        {"group.space_groups.cinter", &G.settings, 0},
        {"settings.inter_settings.number", &G.inter_number, 0},
        {"group.point_groups.generator_count", &G.group_gen_count, 0},
        {"group.point_groups.generators", &G.group_gen, 0},
        {"settings.wyckoff.by_sg_start", &G.wyckoff_pointer, 0},
        {"settings.wyckoff.by_sg_count", &G.wyckoff_count, 0},
        {"settings.wyckoff.fract", &G.wyckoff_fract, 0},
        {"settings.site_settings.by_pg_start", &G.pg_elements_pointer, 0},
        {"settings.site_settings.by_pg_count", &G.pg_settings_count, 0},
        {"settings.site_settings.elements", &G.pg_elements, 0},
    };
    for (size_t i = 0; i < sizeof items / sizeof items[0]; i++) {
        if (section_find(items[i].name, items[i].view)) return GEO_TABLES;
        if (!items[i].text && !section_is_int(items[i].view)) return GEO_TABLES;
    }
    int64_t ops = section_count(&G.point_op_code), groups = section_count(&G.point_group_code);
    G.op_codes = malloc(sizeof(int64_t) * (size_t)(ops + 1));
    G.group_codes = malloc(sizeof(int64_t) * (size_t)(groups + 1));
    if (!G.op_codes || !G.group_codes) return GEO_TABLES;
    for (int64_t i = 0; i < ops; i++) G.op_codes[i] = section_int(&G.point_op_code, i);
    for (int64_t i = 0; i < groups; i++) G.group_codes[i] = section_int(&G.point_group_code, i);
    /* products of the point operations by their matrices, within the cubic (1-48) and the hexagonal (49-72) family */
    if (section_count(&G.point_op) != 72 * 9) return GEO_TABLES;
    for (int a = 0; a < 72; a++)
        for (int b = 0; b < 72; b++) {
            int family = a >= 48, product = 0;
            if ((b >= 48) == family) {
                int64_t m[9];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
                        m[3 * i + j] = 0;
                        for (int k = 0; k < 3; k++) m[3 * i + j] += section_int(&G.point_op, 9 * a + 3 * i + k) * section_int(&G.point_op, 9 * b + 3 * k + j);
                    }
                for (int c = family ? 48 : 0; c < (family ? 72 : 48) && !product; c++) {
                    int same = 1;
                    for (int i = 0; i < 9 && same; i++) same = section_int(&G.point_op, 9 * c + i) == m[i];
                    if (same) product = c + 1;
                }
            }
            G.products[72 * a + b] = (short)product;
        }
    G.records = calloc(231, sizeof *G.records);
    if (!G.records) return GEO_TABLES;
    for (int sg = 1; sg < 230; sg++)
        G.gen_start[sg + 1] = G.gen_start[sg] + section_int(&G.group_gen_count, section_int(&G.sg_point_group, sg - 1) - 1);
    G.ready = 1;
    return 0;
}

static int geo_bad;

int geo_take_error(void)
{
    int status = geo_bad ? geo_bad : (rat_overflow ? GEO_RANGE : 0);
    geo_bad = 0;
    rat_overflow = 0;
    return status;
}

static int64_t at(const section_view *view, int64_t index)
{
    if (index < 0 || index >= section_count(view)) {
        if (!geo_bad) geo_bad = GEO_TABLES;
        return 0;
    }
    return section_int(view, index);
}

const int64_t *geo_point_op_codes(int64_t *count)
{
    *count = section_count(&G.point_op_code);
    return G.op_codes;
}

const int64_t *geo_point_group_codes(int64_t *count)
{
    *count = section_count(&G.point_group_code);
    return G.group_codes;
}

int64_t geo_point_group_order(int point_group) { return at(&G.point_group_order, point_group - 1); }

int64_t geo_point_product(int left, int right)
{
    if (geo_load() || left < 1 || left > 72 || right < 1 || right > 72) {
        if (!geo_bad) geo_bad = GEO_TABLES;
        return 0;
    }
    return G.products[72 * (left - 1) + right - 1];
}

int geo_inter_setting_ids(int sg, int capacity, int *out)
{
    int count = 0;
    int64_t total = section_count(&G.inter_number);
    for (int64_t i = 0; i < total; i++)
        if (section_int(&G.inter_number, i) == sg) {
            if (count == capacity) {
                if (!geo_bad) geo_bad = GEO_RANGE;
                return count;
            }
            out[count++] = (int)i + 1;
        }
    return count;
}

/* The 4x4 row-major affine blocks of the SG-indexed cinter base table
 * (half 0 or 1) and of one inter setting. */
void geo_cinter_base(int sg, int half, int64_t *out)
{
    for (int i = 0; i < 16; i++) out[i] = at(&G.settings, (int64_t)(sg - 1) * 32 + half * 16 + i);
}

void geo_inter_block(int setting_id, int half, int64_t *out)
{
    for (int i = 0; i < 16; i++) out[i] = at(&G.settings_inter, (int64_t)(setting_id - 1) * 32 + half * 16 + i);
}

static signed char pml_matrix[15][72][9], pml_present[15][72], pml_ready[15];

/* The pml point matrices of a lattice: ml . R . ml^-1 for the operations that occur in its space groups. */
static int pml_matrices(int lattice)
{
    if (lattice < 1 || lattice > 14) return 1;
    if (pml_ready[lattice]) return 0;
    for (int sg = 1; sg <= 230; sg++) {
        int64_t records[GEO_MAX_OPS][5];
        if (at(&G.lattice, sg - 1) != lattice) continue;
        int n = geo_sg_records(sg, records);
        for (int i = 0; i < n; i++)
            if (records[i][4] >= 1 && records[i][4] <= 72) pml_present[lattice][records[i][4] - 1] = 1;
    }
    int64_t ml[36], den = geo_ml_setting(lattice, ml);
    for (int op = 0; op < 72; op++) {
        if (!pml_present[lattice][op]) continue;
        int64_t m[9], t[9];
        geo_point_op_matrix(op + 1, m);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                t[i + 3 * j] = 0;
                for (int k = 0; k < 3; k++) t[i + 3 * j] += ml[i + 3 * k] * m[k + 3 * j];
            }
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                int64_t v = 0;
                for (int k = 0; k < 3; k++) v += t[i + 3 * k] * ml[9 + k + 3 * j];
                if (!den || v % den || v / den < -1 || v / den > 1) return 1;
                pml_matrix[lattice][op][i + 3 * j] = (signed char)(v / den);
            }
    }
    pml_ready[lattice] = 1;
    return 0;
}

/* the pml matrix of op in the lattice of sg; 1, with every entry -1, when the lattice's space groups do not use op */
int geo_pml_point_matrix(int sg, int op, int64_t *out)
{
    int lattice = !geo_load() && sg >= 1 && sg <= 230 ? (int)at(&G.lattice, sg - 1) : 0;
    int found = op >= 1 && op <= 72 && !pml_matrices(lattice) && pml_present[lattice][op - 1];
    for (int i = 0; i < 9; i++) out[i] = found ? pml_matrix[lattice][op - 1][i] : -1;
    return !found;
}

int geo_site_settings(int site_pg, int capacity, int (*out)[48], int *lengths)
{
    int64_t pointer = at(&G.pg_elements_pointer, site_pg - 1);
    int count = (int)at(&G.pg_settings_count, site_pg - 1);
    if (count > capacity) {
        if (!geo_bad) geo_bad = GEO_RANGE;
        count = capacity;
    }
    for (int setting = 0; setting < count; setting++) {
        lengths[setting] = 0;
        for (int i = 0; i < 48; i++) {
            int64_t value = at(&G.pg_elements, (pointer - 1 + setting) * 48 + i);
            if (value) out[setting][lengths[setting]++] = (int)value;
        }
    }
    return count;
}

int64_t geo_site_settings_pointer(int site_pg) { return at(&G.pg_elements_pointer, site_pg - 1); }

int geo_sg_point_group(int sg) { return (int)at(&G.sg_point_group, sg - 1); }

int geo_sg_lattice(int sg) { return (int)at(&G.lattice, sg - 1); }

int geo_point_op_count(void) { return geo_load() ? 0 : (int)(section_count(&G.point_op) / 9); }

int64_t geo_point_op_code(int op) { return at(&G.point_op_code, op - 1); }

void geo_point_op_matrix(int op, int64_t *out)
{
    for (int i = 0; i < 9; i++) out[i] = at(&G.point_op, 9 * (int64_t)(op - 1) + i);
}

int64_t geo_ml_setting(int lattice, int64_t *record)
{
    for (int i = 0; i < 36; i++) record[i] = at(&G.lattice_ml, (int64_t)(lattice - 1) * 36 + i);
    return at(&G.lattice_ml_denom, (int64_t)(lattice - 1) * 2);
}

static int reduce_record(int64_t *record)
{
    int64_t den = record[3];
    if (!den) return 0;
    i128 v[3], g = den < 0 ? -den : den;
    for (int axis = 0; axis < 3; axis++) {
        i128 r = record[axis] % den;
        v[axis] = r && (r < 0) != (den < 0) ? r + den : r;
    }
    for (int axis = 0; axis < 3; axis++) g = rat_gcd(g, v[axis]);
    i128 d = den;
    if (d < 0) {
        d = -d;
        for (int axis = 0; axis < 3; axis++) v[axis] = -v[axis];
    }
    d /= g;
    for (int axis = 0; axis < 3; axis++) {
        i128 r = (v[axis] / g) % d;
        record[axis] = (int64_t)(r < 0 ? r + d : r);
    }
    record[3] = (int64_t)d;
    return 0;
}

int geo_pml_generators(int sg, int64_t (*out)[5])
{
    int point_group = geo_sg_point_group(sg);
    int count = (int)at(&G.group_gen_count, point_group - 1);
    if (count < 0 || count > 5) {
        if (!geo_bad) geo_bad = GEO_TABLES;
        return 0;
    }
    int identity = 1;
    for (int i = 0; i < 32; i++)
        if (at(&G.settings, (int64_t)(sg - 1) * 32 + i) != ((i % 16) % 5 == 0)) identity = 0;
    if (identity) {
        for (int index = 0; index < count; index++) {
            for (int part = 0; part < 5; part++) out[index][part] = at(&G.generators, (G.gen_start[sg] + index) * 5 + part);
            reduce_record(out[index]);
        }
        return count;
    }
    int64_t generated[GEO_MAX_OPS][5];
    int total = geo_sg_records(sg, generated), used[73] = {0};
    for (int index = 0; index < count; index++) {
        int op = (int)at(&G.group_gen, (int64_t)(point_group - 1) * 5 + index), found = -1;
        for (int g = 0; g < total; g++)
            if (generated[g][4] == op) found = g;
        if (found < 0) {
            int64_t code = geo_point_op_code(op);
            for (int g = 0; g < total; g++) {
                int candidate = (int)generated[g][4];
                if (candidate < 1 || candidate > 72 || used[candidate] || geo_point_op_code(candidate) != code) continue;
                found = g;
                if (point_group == 18 && op == 58) break;
            }
        }
        if (found < 0) return -1;
        memcpy(out[index], generated[found], sizeof out[index]);
        if (out[index][4] >= 1 && out[index][4] <= 72) used[out[index][4]] = 1;
        reduce_record(out[index]);
    }
    return count;
}

static void vmlt(const int64_t *m, const rat *v, rat *out);

/* The elements of a space group: the closure of its cml generators, by ascending point operation, the translations
 * written in the ml basis and reduced into [0, 1). */
static int build_records(int sg)
{
    int point_group = (int)at(&G.sg_point_group, sg - 1), count = (int)geo_point_group_order(point_group);
    int ng = (int)at(&G.group_gen_count, point_group - 1);
    if (count < 1 || count > GEO_MAX_OPS || ng < 0 || ng > 5) return 1;
    int op[GEO_MAX_OPS], gop[5], m = 1;
    rat t[GEO_MAX_OPS][3], gt[5][3];
    op[0] = 1;
    for (int g = 0; g < ng; g++) {
        int64_t r[5];
        for (int k = 0; k < 5; k++) r[k] = at(&G.generators, (G.gen_start[sg] + g) * 5 + k);
        gop[g] = (int)r[4];
        for (int a = 0; a < 3; a++) gt[g][a] = rat_mod1(rat_make(r[a], r[3]));
        if (gop[g] > 48) op[0] = 49;
    }
    for (int a = 0; a < 3; a++) t[0][a] = rat_int(0);
    for (int i = 0; i < m; i++)
        for (int g = 0; g < ng; g++) {
            int product = (int)geo_point_product(gop[g], op[i]), seen = 0;
            for (int j = 0; j < m && !seen; j++) seen = op[j] == product;
            if (seen) continue;
            if (m == count) return 1;
            int64_t point[9];
            geo_point_op_matrix(op[i], point);
            vmlt(point, gt[g], t[m]);
            for (int a = 0; a < 3; a++) t[m][a] = rat_mod1(rat_add(t[m][a], t[i][a]));
            op[m++] = product;
        }
    if (m != count) return 1;
    int64_t ml[36];
    geo_ml_setting((int)at(&G.lattice, sg - 1), ml);
    for (int i = 1; i < m; i++)
        for (int j = i; j > 0 && op[j - 1] > op[j]; j--) {
            int o = op[j];
            op[j] = op[j - 1], op[j - 1] = o;
            for (int a = 0; a < 3; a++) {
                rat x = t[j][a];
                t[j][a] = t[j - 1][a], t[j - 1][a] = x;
            }
        }
    for (int i = 0; i < m; i++) {
        rat v[3];
        vmlt(ml, t[i], v);
        i128 den = 1;
        for (int a = 0; a < 3; a++) {
            v[a] = rat_mod1(v[a]);
            den = den / rat_gcd(den, v[a].d) * v[a].d;
        }
        int64_t *record = G.records[sg][i];
        for (int a = 0; a < 3; a++) record[a] = (int64_t)(v[a].n * (den / v[a].d));
        record[3] = (int64_t)den;
        record[4] = op[i];
    }
    G.order[sg] = m;
    return 0;
}

int geo_sg_records(int sg, int64_t (*out)[5])
{
    if (geo_load() || sg < 1 || sg > 230 || (!G.order[sg] && build_records(sg))) {
        if (!geo_bad) geo_bad = GEO_TABLES;
        return 0;
    }
    memcpy(out, G.records[sg], sizeof(int64_t) * 5 * (size_t)G.order[sg]);
    return G.order[sg];
}

/* vmlt_ with a column-major integer matrix. */
static void vmlt(const int64_t *m, const rat *v, rat *out)
{
    rat t[3];
    for (int axis = 0; axis < 3; axis++)
        t[axis] = rat_add(rat_add(rat_mul(rat_int(m[axis]), v[0]), rat_mul(rat_int(m[axis + 3]), v[1])),
                          rat_mul(rat_int(m[axis + 6]), v[2]));
    memcpy(out, t, sizeof t);
}

void geo_vrot(int sg, int op, const rat *v, rat *out)
{
    int64_t lattice = at(&G.lattice, sg - 1);
    int64_t back[9], forward[9], point[9];
    for (int i = 0; i < 9; i++) {
        back[i] = at(&G.lattice_ml, (lattice - 1) * 36 + i);
        forward[i] = at(&G.lattice_ml, (lattice - 1) * 36 + 9 + i);
        point[i] = at(&G.point_op, 9 * (int64_t)(op - 1) + i);
    }
    int64_t den = at(&G.lattice_ml_denom, (lattice - 1) * 2);
    if (!den) {
        if (!geo_bad) geo_bad = GEO_TABLES;
        return;
    }
    rat t[3];
    vmlt(forward, v, t);
    vmlt(point, t, t);
    vmlt(back, t, t);
    for (int axis = 0; axis < 3; axis++) out[axis] = rat_div(t[axis], rat_int(den));
}

int geo_wyckoff_row_range(int sg, int *start, int *count)
{
    *start = (int)at(&G.wyckoff_pointer, sg - 1);
    *count = (int)at(&G.wyckoff_count, sg - 1);
    return geo_bad;
}

int wyckoff_site_pg(int row_id);

int geo_wyckoff_site_pg(int row_id) { return wyckoff_site_pg(row_id); }

static void vectors_from(const section_view *view, int64_t base, rat (*out)[3])
{
    for (int piece = 0; piece < 4; piece++) {
        int64_t den = at(view, base + 4 * piece + 3);
        for (int axis = 0; axis < 3; axis++)
            out[piece][axis] = den ? rat_make(at(view, base + 4 * piece + axis), den) : rat_int(0);
    }
}

void geo_wyckoff_vectors(int row_id, rat (*out)[3]) { vectors_from(&G.wyckoff_fract, (int64_t)(row_id - 1) * 16, out); }

void geo_inter_wyckoff_vectors(int setting_id, int offset0, rat (*out)[3])
{
    int64_t pointer = at(&G.inter_wyckoff_pointer, setting_id - 1);
    vectors_from(&G.inter_wyckoff_fract, (pointer - 1 + offset0) * 16, out);
}

int geo_default_setting(int sg) { return (int)at(&G.inter_choice, sg - 1); }

static void operation_record(const rat *v, int op, int64_t *out)
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
    for (int axis = 0; axis < 3; axis++) out[axis] = (int64_t)n[axis];
    out[3] = (int64_t)den;
    out[4] = op;
}

int geo_site_records(int sg, int row_id, int site_pg, int64_t (*out)[5])
{
    rat vectors[4][3];
    geo_wyckoff_vectors(row_id, vectors);
    int64_t records[GEO_MAX_OPS][5];
    int count = geo_sg_records(sg, records), found = 0;
    for (int index = 0; index < count; index++) {
        int op = (int)records[index][4];
        rat rotated[3], shift[3];
        geo_vrot(sg, op, vectors[0], rotated);
        int keep = 1;
        for (int axis = 0; axis < 3 && keep; axis++) {
            shift[axis] = rat_sub(vectors[0][axis], rotated[axis]);
            keep = rat_zero(rat_mod1(rat_sub(shift[axis], rat_make(records[index][axis], records[index][3]))));
        }
        for (int p = 1; p < 4 && keep; p++) {
            rat image[3];
            geo_vrot(sg, op, vectors[p], image);
            keep = rat_eq(image[0], vectors[p][0]) && rat_eq(image[1], vectors[p][1]) && rat_eq(image[2], vectors[p][2]);
        }
        if (keep) operation_record(shift, op, out[found++]);
    }
    int64_t pointer = at(&G.pg_elements_pointer, site_pg - 1);
    int64_t settings = at(&G.pg_settings_count, site_pg - 1);
    for (int64_t setting = 0; setting < settings; setting++) {
        int64_t ops[48];
        int length = 0;
        for (int i = 0; i < 48; i++) {
            int64_t value = at(&G.pg_elements, (pointer - 1 + setting) * 48 + i);
            if (value) ops[length++] = value;
        }
        if (length != found) continue;
        int same = 1;
        for (int i = 0; i < length && same; i++) {
            int hit = 0;
            for (int j = 0; j < found && !hit; j++) hit = out[j][4] == ops[i];
            same = hit;
        }
        for (int j = 0; j < found && same; j++) {
            int hit = 0;
            for (int i = 0; i < length && !hit; i++) hit = out[j][4] == ops[i];
            same = hit;
        }
        if (!same) continue;
        int64_t ordered[GEO_MAX_OPS][5];
        for (int i = 0; i < length; i++)
            for (int j = 0; j < found; j++)
                if (out[j][4] == ops[i]) memcpy(ordered[i], out[j], sizeof ordered[i]);
        memcpy(out, ordered, sizeof(int64_t) * 5 * (size_t)length);
        return length;
    }
    return found;
}

int geo_cosets(int sg, const int64_t (*site)[5], int site_count, int64_t (*out)[5])
{
    int64_t generated[GEO_MAX_OPS][5];
    int count = geo_sg_records(sg, generated), reps = 0;
    int visited[GEO_MAX_OPS] = {0};
    for (int index = 0; index < count; index++) {
        if (visited[index]) continue;
        memcpy(out[reps++], generated[index], sizeof generated[index]);
        for (int s = 0; s < site_count; s++) {
            int64_t product = geo_point_product((int)site[s][4], (int)generated[index][4]);
            int mark = -1;
            for (int j = 0; j < count; j++)
                if (generated[j][4] == product) mark = j;
            if (mark < 0 || visited[mark]) {
                geo_bad = GEO_COSETS;
                return 0;
            }
            visited[mark] = 1;
        }
    }
    if (reps * site_count != count) {
        geo_bad = GEO_COSETS;
        return 0;
    }
    return reps;
}

static const char *first_char(const section_view *view, int64_t index)
{
    const char *text;
    int64_t length;
    if (section_string(view, index, &text, &length) || !length) return 0;
    return text;
}

int64_t geo_multiplicity(int sg, int site_pg)
{
    const char *symbol = first_char(&G.label_bc, sg - 1);
    int64_t centering = 1;
    if (symbol) switch (symbol[0]) {
        case 'A': case 'B': case 'C': case 'I': centering = 2; break;
        case 'F': centering = 4; break;
        case 'R': centering = 3; break;
        default: break;
        }
    int64_t parent = geo_point_group_order((int)at(&G.sg_point_group, sg - 1));
    int64_t site = geo_point_group_order(site_pg);
    if (!site) {
        if (!geo_bad) geo_bad = GEO_TABLES;
        return 0;
    }
    return centering * parent / site;
}

int geo_setting_to_cinter(int sg, int pml, int setting_id, rat (*matrix)[3], rat *origin)
{
    int64_t choice = setting_id ? setting_id : at(&G.inter_choice, sg - 1);
    int64_t base = (choice - 1) * 32, den = at(&G.settings_inter, base + 15);
    if (!den) return geo_bad = GEO_TABLES;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) matrix[row][col] = rat_make(at(&G.settings_inter, base + 4 * row + col), den);
    for (int axis = 0; axis < 3; axis++) origin[axis] = rat_make(at(&G.settings_inter, base + 12 + axis), den);
    if (pml) {
        int64_t lattice = at(&G.lattice, sg - 1);
        int64_t pml_den = at(&G.lattice_ml_denom, (lattice - 1) * 2);
        if (!pml_den) return geo_bad = GEO_TABLES;
        rat to_cml[3][3], product[3][3];
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
                to_cml[row][col] = rat_make(at(&G.lattice_ml, (lattice - 1) * 36 + 9 + 3 * row + col), pml_den);
        for (int row = 0; row < 3; row++) rat_row_mul(to_cml[row], (const rat (*)[3])matrix, product[row]);
        memcpy(matrix, product, sizeof product);
    }
    return geo_bad;
}

int64_t geo_kernel_fractions(const int64_t *b, int64_t **out)
{
    int64_t adj[9] = {
        b[4] * b[8] - b[5] * b[7], b[2] * b[7] - b[1] * b[8], b[1] * b[5] - b[2] * b[4],
        b[5] * b[6] - b[3] * b[8], b[0] * b[8] - b[2] * b[6], b[2] * b[3] - b[0] * b[5],
        b[3] * b[7] - b[4] * b[6], b[1] * b[6] - b[0] * b[7], b[0] * b[4] - b[1] * b[3],
    };
    int64_t det = b[0] * adj[0] + b[1] * adj[3] + b[2] * adj[6];
    if (det < 0) det = -det;
    *out = 0;
    if (!det) return -GEO_RANGE;
    int64_t factor = det;
    for (int i = 0; i < 9; i++) factor = (int64_t)rat_gcd(factor, adj[i]);
    int64_t den = det / factor;
    int64_t *images = malloc(sizeof(int64_t) * 4 * (size_t)det);
    if (!images) return -GEO_RANGE;
    int64_t count = 0;
    for (int64_t x = 0; x < den; x++)
        for (int64_t y = 0; y < den; y++)
            for (int64_t z = 0; z < den; z++) {
                int64_t n[3];
                int integral = 1;
                for (int col = 0; col < 3 && integral; col++) {
                    n[col] = x * b[col] + y * b[3 + col] + z * b[6 + col];
                    integral = n[col] % den == 0;
                }
                if (!integral) continue;
                if (count >= det) {
                    free(images);
                    return -GEO_RANGE;
                }
                for (int col = 0; col < 3; col++) images[4 * count + col] = n[col] / den;
                images[4 * count + 3] = 1;
                count++;
            }
    if (count != det) {
        free(images);
        return -GEO_RANGE;
    }
    *out = images;
    return count;
}

int geo_centering(int sg, rat (*out)[3])
{
    int64_t choice = at(&G.inter_choice, sg - 1);
    const char *symbol = first_char(&G.inter_label_short, choice - 1);
    rat z = rat_int(0), h = rat_make(1, 2);
    for (int axis = 0; axis < 3; axis++) out[0][axis] = z;
    switch (symbol ? symbol[0] : 'P') {
    case 'C': out[1][0] = h; out[1][1] = h; out[1][2] = z; return 2;
    case 'A': out[1][0] = z; out[1][1] = h; out[1][2] = h; return 2;
    case 'B': out[1][0] = h; out[1][1] = z; out[1][2] = h; return 2;
    case 'I': out[1][0] = h; out[1][1] = h; out[1][2] = h; return 2;
    case 'F':
        out[1][0] = z; out[1][1] = h; out[1][2] = h;
        out[2][0] = h; out[2][1] = z; out[2][2] = h;
        out[3][0] = h; out[3][1] = h; out[3][2] = z;
        return 4;
    case 'R':
        out[1][0] = rat_make(2, 3); out[1][1] = rat_make(1, 3); out[1][2] = rat_make(1, 3);
        out[2][0] = rat_make(1, 3); out[2][1] = rat_make(2, 3); out[2][2] = rat_make(2, 3);
        return 3;
    default: return 1;
    }
}
