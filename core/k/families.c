/* Source k families of a supercell's reciprocal mesh.
 *
 * The mesh of a basis is scanned in Source first-occurrence order; every
 * vector is identified with the little-k families of the parent space group
 * (data_little) whose star holds it, with each family's parameters brought to
 * the kvec_standard_ representative.  All arithmetic is exact over reduced
 * int64 rationals with int128 intermediates; anything wider is an error. */

#include <stdlib.h>
#include <string.h>

#include "intmath.h"
#include "geometry.h"
#include "tables.h"

typedef struct {
    int64_t n, d;
} q64;

enum {
    KF_OK = 0,
    KF_OVERFLOW = 1,
    KF_TABLES = 2,
    KF_MEMORY = 3,
    KF_CAPACITY = 4,
    KF_NO_IDENTITY = 5,
    KF_ZERO_DENOMINATOR = 6,
    KF_NO_SPECIAL_VECTOR = 7,
};

enum {
    KF_DYNAMIC = 0,
    KF_DISPLAY = 1,
};

#define MAX_FAMILIES 32
#define MAX_OPS 96
#define MAX_ARMS 48

int kstar_solutions(int mode, int kdim, int arm_count, const int64_t *origin,
                    const int *valid, const int64_t *left, const int64_t *consistency,
                    const int64_t *scale, const int64_t *target, int capacity, int64_t *out);

static int bad;

static void fail(int status)
{
    if (!bad) bad = status;
}

static i128 gcd128(i128 a, i128 b) { return gcd_i128(a, b); }

static int64_t lcm64(int64_t a, int64_t b)
{
    i128 g = gcd128(a, b);
    if (!g) return 0;
    i128 value = (i128)a / g * b;
    if (value < 0) value = -value;
    if (value > ((i128)1 << 62)) {
        fail(KF_OVERFLOW);
        return 1;
    }
    return (int64_t)value;
}

static q64 qmake(i128 n, i128 d)
{
    if (!d) {
        fail(KF_ZERO_DENOMINATOR);
        return (q64){0, 1};
    }
    if (d < 0) {
        n = -n;
        d = -d;
    }
    i128 g = gcd128(n, d);
    if (g > 1) {
        n /= g;
        d /= g;
    }
    const i128 limit = (i128)1 << 62;
    if (n > limit || n < -limit || d > limit) {
        fail(KF_OVERFLOW);
        return (q64){0, 1};
    }
    return (q64){(int64_t)n, (int64_t)d};
}

static q64 qint(int64_t n) { return (q64){n, 1}; }
static q64 qadd(q64 a, q64 b) { return qmake((i128)a.n * b.d + (i128)b.n * a.d, (i128)a.d * b.d); }
static q64 qsub(q64 a, q64 b) { return qmake((i128)a.n * b.d - (i128)b.n * a.d, (i128)a.d * b.d); }
static q64 qmul(q64 a, q64 b) { return qmake((i128)a.n * b.n, (i128)a.d * b.d); }
static q64 qdiv(q64 a, q64 b) { return qmake((i128)a.n * b.d, (i128)a.d * b.n); }
static q64 qneg(q64 a) { return (q64){-a.n, a.d}; }
static q64 qabs(q64 a) { return (q64){a.n < 0 ? -a.n : a.n, a.d}; }
static int qeq(q64 a, q64 b) { return a.n == b.n && a.d == b.d; }

static int qcmp(q64 a, q64 b)
{
    i128 left = (i128)a.n * b.d, right = (i128)b.n * a.d;
    return (left > right) - (left < right);
}

static int64_t pymod(int64_t a, int64_t m)
{
    int64_t r = a % m;
    if (r && ((r < 0) != (m < 0))) r += m;
    return r;
}

static q64 qmod1(q64 a) { return (q64){pymod(a.n, a.d), a.d}; }

/* Python's _key: (numerator mod denominator, denominator) per axis. */
static int key_cmp(const q64 *a, const q64 *b)
{
    for (int axis = 0; axis < 3; axis++) {
        int64_t an = pymod(a[axis].n, a[axis].d), bn = pymod(b[axis].n, b[axis].d);
        if (an != bn) return an < bn ? -1 : 1;
        if (a[axis].d != b[axis].d) return a[axis].d < b[axis].d ? -1 : 1;
    }
    return 0;
}

/* (|v|, v < 0) per axis, lexicographic. */
static int cinter_cmp(const q64 *a, const q64 *b)
{
    for (int axis = 0; axis < 3; axis++) {
        int c = qcmp(qabs(a[axis]), qabs(b[axis]));
        if (c) return c;
        int na = a[axis].n < 0, nb = b[axis].n < 0;
        if (na != nb) return na - nb;
    }
    return 0;
}

static struct {
    section_view irr_sg, irr_k, k_label, k_dim, star_count, star_pointer, star, star_special, little_k;
    section_view lattice, lattice_ml, lattice_ml_denom, point_op;
    section_view sg_point_group, point_group_order, inter_choice, settings_inter;
    int ready;
} T;

static int load_tables(void)
{
    if (T.ready) return 0;
    struct {
        const char *name;
        section_view *view;
    } items[] = {
        {"irreps.little.sg", &T.irr_sg},
        {"irreps.little.k", &T.irr_k},
        {"k.slots.label", &T.k_label},
        {"k.slots.dim", &T.k_dim},
        {"k.star.arm_count", &T.star_count},
        {"k.star.arm_start", &T.star_pointer},
        {"k.star.arms", &T.star},
        {"k.star.special_arms", &T.star_special},
        {"k.slots.vector", &T.little_k},
        {"group.space_groups.lattice", &T.lattice},
        {"group.lattices.ml", &T.lattice_ml},
        {"group.lattices.ml_denom", &T.lattice_ml_denom},
        {"group.point_ops.matrix", &T.point_op},
        {"group.space_groups.point_group", &T.sg_point_group},
        {"group.point_groups.order", &T.point_group_order},
        {"group.space_groups.inter_choice", &T.inter_choice},
        {"settings.inter_settings.inter", &T.settings_inter},
    };
    for (size_t i = 0; i < sizeof items / sizeof items[0]; i++) {
        if (section_find(items[i].name, items[i].view)) return KF_TABLES;
        if (items[i].view != &T.k_label && !section_is_int(items[i].view)) return KF_TABLES;
    }
    T.ready = 1;
    return 0;
}

static int64_t at(const section_view *view, int64_t index)
{
    if (index < 0 || index >= section_count(view)) {
        fail(KF_TABLES);
        return 0;
    }
    return section_int(view, index);
}

/* Process-lifetime memo of pure table functions: fixed-width int64 keys to
 * owned int64 values. */
typedef struct {
    int width;
    size_t cap, count;
    int64_t *keys;
    int64_t **values;
} memo;

static uint64_t memo_hash(const int64_t *key, int width)
{
    uint64_t h = 1469598103934665603ull;
    for (int i = 0; i < width; i++) {
        h ^= (uint64_t)key[i];
        h *= 1099511628211ull;
        h ^= h >> 29;
    }
    return h;
}

static int64_t *memo_get(const memo *m, const int64_t *key)
{
    if (!m->cap) return 0;
    for (size_t i = memo_hash(key, m->width) & (m->cap - 1);; i = (i + 1) & (m->cap - 1)) {
        if (!m->values[i]) return 0;
        if (!memcmp(m->keys + i * m->width, key, sizeof(int64_t) * m->width)) return m->values[i];
    }
}

static int memo_put(memo *m, const int64_t *key, int64_t *value)
{
    if (2 * (m->count + 1) > m->cap) {
        size_t cap = m->cap ? 2 * m->cap : 1024;
        int64_t *keys = malloc(sizeof(int64_t) * m->width * cap);
        int64_t **values = calloc(cap, sizeof(int64_t *));
        if (!keys || !values) {
            free(keys);
            free(values);
            return KF_MEMORY;
        }
        for (size_t j = 0; j < m->cap; j++) {
            if (!m->values[j]) continue;
            size_t i = memo_hash(m->keys + j * m->width, m->width) & (cap - 1);
            while (values[i]) i = (i + 1) & (cap - 1);
            memcpy(keys + i * m->width, m->keys + j * m->width, sizeof(int64_t) * m->width);
            values[i] = m->values[j];
        }
        free(m->keys);
        free(m->values);
        m->keys = keys;
        m->values = values;
        m->cap = cap;
    }
    size_t i = memo_hash(key, m->width) & (m->cap - 1);
    while (m->values[i]) i = (i + 1) & (m->cap - 1);
    memcpy(m->keys + i * m->width, key, sizeof(int64_t) * m->width);
    m->values[i] = value;
    m->count++;
    return 0;
}

static void write_vector(int64_t *out, const q64 *v)
{
    for (int axis = 0; axis < 3; axis++) {
        out[2 * axis] = v[axis].n;
        out[2 * axis + 1] = v[axis].d;
    }
}

static memo standard_memo = {6, 0, 0, 0, 0};
static memo orbit_memo = {8, 0, 0, 0, 0};
static memo candidate_memo = {7, 0, 0, 0, 0};

typedef struct {
    int sg, lattice;
    int64_t reciprocal[36], reciprocal_den;
    int lattice_count;
    int64_t lattice_rot[MAX_OPS][9];
    int standard_count;
    int standard_op[MAX_OPS];
    int64_t standard_point[MAX_OPS][9];
    int64_t standard_rot[MAX_OPS][9];
    int identity_order;
    q64 to_cinter[3][3];
    int family_count;
    int family_gid[MAX_FAMILIES], family_kdim[MAX_FAMILIES];
    const char *family_label[MAX_FAMILIES];
    int64_t family_label_length[MAX_FAMILIES];
} sg_context;

static int sg_records(int sg, int (*records)[5])
{
    int64_t generated[GEO_MAX_OPS][5];
    int count = geo_load() ? 0 : geo_sg_records(sg, generated);
    if (count <= 0 || count > MAX_OPS) {
        fail(KF_TABLES);
        return 0;
    }
    for (int index = 0; index < count; index++)
        for (int part = 0; part < 5; part++) records[index][part] = (int)generated[index][part];
    return count;
}

static void point_matrix(int op, int64_t *out)
{
    for (int i = 0; i < 9; i++) out[i] = at(&T.point_op, 9 * (op - 1) + i);
}

/* record[18:27] . point . record[27:36] with column-major 3x3 blocks, as a
 * row-major matrix. */
static void compose_rotation(const int64_t *record, const int64_t *point, int64_t *out)
{
    int64_t left[3][3], middle[3][3], right[3][3], tmp[3][3];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            left[row][col] = record[18 + row + 3 * col];
            middle[row][col] = point[row + 3 * col];
            right[row][col] = record[27 + row + 3 * col];
        }
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            tmp[row][col] = 0;
            for (int k = 0; k < 3; k++) tmp[row][col] += left[row][k] * middle[k][col];
        }
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            int64_t value = 0;
            for (int k = 0; k < 3; k++) value += tmp[row][k] * right[k][col];
            out[3 * row + col] = value;
        }
}

static void rotate(const sg_context *ctx, const int64_t *rot, const q64 *v, int reduce, q64 *out)
{
    int64_t common = 1;
    for (int axis = 0; axis < 3; axis++) common = lcm64(common, v[axis].d);
    i128 x[3];
    for (int axis = 0; axis < 3; axis++) x[axis] = (i128)v[axis].n * (common / v[axis].d);
    i128 den = (i128)common * ctx->reciprocal_den;
    for (int row = 0; row < 3; row++) {
        i128 num = rot[3 * row] * x[0] + rot[3 * row + 1] * x[1] + rot[3 * row + 2] * x[2];
        if (reduce) {
            i128 r = num % den;
            if (r && ((r < 0) != (den < 0))) r += den;
            num = r;
        }
        out[row] = qmake(num, den);
    }
}

static int label_equal(const char *a, int64_t alen, const char *b, int64_t blen)
{
    return alen == blen && !memcmp(a, b, (size_t)alen);
}

static int load_sg_context(int sg, sg_context *ctx)
{
    memset(ctx, 0, sizeof *ctx);
    ctx->sg = sg;
    ctx->lattice = (int)at(&T.lattice, sg - 1);
    int lattice = ctx->lattice;
    for (int i = 0; i < 36; i++) ctx->reciprocal[i] = at(&T.lattice_ml, (lattice - 1) * 36 + i);
    ctx->reciprocal_den = at(&T.lattice_ml_denom, (lattice - 1) * 2 + 1);
    if (!ctx->reciprocal_den) return KF_ZERO_DENOMINATOR;

    int records[MAX_OPS][5];
    int representative = 0;
    for (int candidate = 1; candidate <= section_count(&T.lattice); candidate++)
        if (at(&T.lattice, candidate - 1) == lattice) representative = candidate;
    if (!representative) return KF_TABLES;
    ctx->lattice_count = sg_records(representative, records);
    for (int i = 0; i < ctx->lattice_count; i++) {
        int64_t point[9];
        point_matrix(records[i][4], point);
        compose_rotation(ctx->reciprocal, point, ctx->lattice_rot[i]);
    }

    int count = sg_records(sg, records);
    for (int i = 0; i < count; i++) {
        int op = records[i][4], seen = 0;
        for (int j = 0; j < ctx->standard_count; j++) seen |= ctx->standard_op[j] == op;
        if (!seen) ctx->standard_op[ctx->standard_count++] = op;
    }
    if (!ctx->standard_count) return KF_TABLES;
    int inversion = ctx->standard_op[0] == 1 ? 25 : 61, has_inversion = 0;
    for (int j = 0; j < ctx->standard_count; j++) has_inversion |= ctx->standard_op[j] == inversion;
    if (!has_inversion) {
        int original = ctx->standard_count;
        for (int j = 0; j < original; j++) {
            int mate = (int)geo_point_product(inversion, ctx->standard_op[j]), seen = 0;
            for (int i = 0; i < ctx->standard_count; i++) seen |= ctx->standard_op[i] == mate;
            if (!seen) {
                if (ctx->standard_count >= MAX_OPS) return KF_TABLES;
                ctx->standard_op[ctx->standard_count++] = mate;
            }
        }
    }
    ctx->identity_order = -1;
    for (int j = 0; j < ctx->standard_count; j++) {
        point_matrix(ctx->standard_op[j], ctx->standard_point[j]);
        compose_rotation(ctx->reciprocal, ctx->standard_point[j], ctx->standard_rot[j]);
        const int64_t *p = ctx->standard_point[j];
        if (ctx->identity_order < 0 && p[0] == 1 && p[4] == 1 && p[8] == 1 && !p[1] && !p[2] && !p[3]
            && !p[5] && !p[6] && !p[7])
            ctx->identity_order = j;
    }

    int64_t pml_den = at(&T.lattice_ml_denom, (lattice - 1) * 2);
    int64_t choice = at(&T.inter_choice, sg - 1);
    int64_t inter_den = at(&T.settings_inter, (choice - 1) * 32 + 15);
    if (!pml_den || !inter_den) return KF_ZERO_DENOMINATOR;
    q64 pml[3][3], inter[3][3], forward[3][3];
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            pml[row][col] = qmake(at(&T.lattice_ml, (lattice - 1) * 36 + 9 + 3 * row + col), pml_den);
            inter[row][col] = qmake(at(&T.settings_inter, (choice - 1) * 32 + 4 * row + col), inter_den);
        }
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            q64 acc = qint(0);
            for (int k = 0; k < 3; k++) acc = qadd(acc, qmul(pml[row][k], inter[k][col]));
            forward[row][col] = acc;
        }
    q64 det = qint(0);
    for (int col = 0; col < 3; col++) {
        q64 minor = qsub(qmul(forward[1][(col + 1) % 3], forward[2][(col + 2) % 3]),
                         qmul(forward[1][(col + 2) % 3], forward[2][(col + 1) % 3]));
        det = qadd(det, qmul(forward[0][col], minor));
    }
    if (!det.n) return KF_ZERO_DENOMINATOR;
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) {
            int r1 = (col + 1) % 3, r2 = (col + 2) % 3, c1 = (row + 1) % 3, c2 = (row + 2) % 3;
            q64 cofactor = qsub(qmul(forward[r1][c1], forward[r2][c2]), qmul(forward[r1][c2], forward[r2][c1]));
            ctx->to_cinter[row][col] = qdiv(cofactor, det);
        }

    int64_t gid_count = section_count(&T.irr_sg);
    for (int64_t gid = 1; gid <= gid_count; gid++) {
        if (at(&T.irr_sg, gid - 1) != sg) continue;
        int64_t kslot = at(&T.irr_k, gid - 1);
        int64_t slot = (lattice - 1) * 27 + kslot - 1;
        const char *label;
        int64_t length;
        if (section_string(&T.k_label, slot, &label, &length)) return KF_TABLES;
        int seen = 0;
        for (int f = 0; f < ctx->family_count && !seen; f++)
            seen = label_equal(label, length, ctx->family_label[f], ctx->family_label_length[f]);
        if (seen) continue;
        if (ctx->family_count >= MAX_FAMILIES) return KF_TABLES;
        int f = ctx->family_count++;
        ctx->family_gid[f] = (int)gid;
        ctx->family_kdim[f] = (int)at(&T.k_dim, slot);
        ctx->family_label[f] = label;
        ctx->family_label_length[f] = length;
    }
    return bad;
}

static void to_cinter(const sg_context *ctx, const q64 *v, q64 *out)
{
    for (int row = 0; row < 3; row++) {
        q64 acc = qint(0);
        for (int col = 0; col < 3; col++) acc = qadd(acc, qmul(ctx->to_cinter[row][col], v[col]));
        out[row] = acc;
    }
}

static int k_dimension(const sg_context *ctx, int gid)
{
    int64_t kslot = at(&T.irr_k, gid - 1);
    return (int)at(&T.k_dim, (ctx->lattice - 1) * 27 + kslot - 1);
}

static void star_location(int gid, int *count, int64_t *pointer)
{
    int64_t sg = at(&T.irr_sg, gid - 1), kslot = at(&T.irr_k, gid - 1);
    int64_t slot = (sg - 1) * 27 + kslot - 1;
    *count = (int)at(&T.star_count, slot);
    *pointer = at(&T.star_pointer, slot);
}

static void star_piece(int64_t pointer, int arm, int piece, q64 *out)
{
    int64_t base = 16 * (pointer - 1 + arm) + 4 * piece;
    int64_t den = at(&T.star, base + 3);
    for (int axis = 0; axis < 3; axis++) out[axis] = den ? qmake(at(&T.star, base + axis), den) : qint(0);
}

/* The fixed special-k vectors of a gid (little_k_star_special). */
static int special_vectors(int gid, q64 (*out)[3])
{
    int count;
    int64_t pointer;
    star_location(gid, &count, &pointer);
    if (count < 0 || count > MAX_ARMS || (count && pointer <= 0)) {
        fail(KF_TABLES);
        return 0;
    }
    for (int arm = 0; arm < count; arm++) {
        int64_t base = 4 * (pointer - 1 + arm);
        int64_t den = at(&T.star_special, base + 3);
        if (!den) {
            fail(KF_ZERO_DENOMINATOR);
            return 0;
        }
        for (int axis = 0; axis < 3; axis++) out[arm][axis] = qmake(at(&T.star_special, base + axis), den);
    }
    return count;
}

static void kparam_params(const int64_t *kparam, q64 *params)
{
    for (int i = 0; i < 3; i++) params[i] = kparam[3] ? qmake(kparam[i], kparam[3]) : qint(kparam[i]);
}

static void affine_at(const q64 (*pieces)[3], int used, const q64 *params, q64 *out)
{
    for (int axis = 0; axis < 3; axis++) {
        q64 acc = pieces[0][axis];
        for (int i = 0; i < used; i++) acc = qadd(acc, qmul(params[i], pieces[i + 1][axis]));
        out[axis] = acc;
    }
}

/* Every star arm of a gid at kparam. */
static int star_at(int gid, const int64_t *kparam, q64 (*out)[3])
{
    int count;
    int64_t pointer;
    star_location(gid, &count, &pointer);
    if (count < 0 || count > MAX_ARMS) {
        fail(KF_TABLES);
        return 0;
    }
    q64 params[3];
    kparam_params(kparam, params);
    for (int arm = 0; arm < count; arm++) {
        q64 pieces[4][3];
        for (int piece = 0; piece < 4; piece++) star_piece(pointer, arm, piece, pieces[piece]);
        affine_at((const q64 (*)[3])pieces, 3, params, out[arm]);
    }
    return count;
}

static int first_arm_at(int gid, const int64_t *kparam, q64 *out)
{
    int count;
    int64_t pointer;
    star_location(gid, &count, &pointer);
    if (count <= 0) return 0;
    q64 params[3], pieces[4][3];
    kparam_params(kparam, params);
    for (int piece = 0; piece < 4; piece++) star_piece(pointer, 0, piece, pieces[piece]);
    affine_at((const q64 (*)[3])pieces, 3, params, out);
    return 1;
}

/* Per-arm integer solvers of origin + params . directions = target (mod the
 * reciprocal shifts), cached per gid. */
typedef struct {
    int ready, kdim, count;
    int64_t origin[MAX_ARMS * 6];
    int valid[MAX_ARMS];
    int64_t left[MAX_ARMS * 9], consistency[MAX_ARMS * 9], scale[MAX_ARMS];
} star_solver;

static star_solver **solvers;
static int64_t solver_capacity;

static const star_solver *solver_for(int gid, int kdim)
{
    if (!solvers) {
        solver_capacity = section_count(&T.irr_sg) + 1;
        solvers = calloc((size_t)solver_capacity, sizeof *solvers);
        if (!solvers) {
            fail(KF_MEMORY);
            return 0;
        }
    }
    if (gid < 1 || gid >= solver_capacity) {
        fail(KF_TABLES);
        return 0;
    }
    if (solvers[gid] && solvers[gid]->kdim == kdim) return solvers[gid];
    star_solver *s = calloc(1, sizeof *s);
    if (!s) {
        fail(KF_MEMORY);
        return 0;
    }
    int64_t pointer;
    star_location(gid, &s->count, &pointer);
    if (s->count < 0 || s->count > MAX_ARMS || kdim < 0 || kdim > 3) {
        free(s);
        fail(KF_TABLES);
        return 0;
    }
    s->kdim = kdim;
    for (int arm = 0; arm < s->count; arm++) {
        q64 pieces[4][3];
        for (int piece = 0; piece < 4; piece++) star_piece(pointer, arm, piece, pieces[piece]);
        for (int axis = 0; axis < 3; axis++) {
            s->origin[6 * arm + 2 * axis] = pieces[0][axis].n;
            s->origin[6 * arm + 2 * axis + 1] = pieces[0][axis].d;
        }
        s->scale[arm] = 1;
        q64 m[3][6];
        int width = kdim + 3;
        for (int axis = 0; axis < 3; axis++) {
            for (int index = 0; index < kdim; index++) m[axis][index] = pieces[index + 1][axis];
            for (int col = 0; col < 3; col++) m[axis][kdim + col] = qint(axis == col);
        }
        int pivot_row = 0;
        for (int col = 0; col < kdim; col++) {
            int pivot = -1;
            for (int row = pivot_row; row < 3 && pivot < 0; row++)
                if (m[row][col].n) pivot = row;
            if (pivot < 0) break;
            for (int j = 0; j < width; j++) {
                q64 t = m[pivot_row][j];
                m[pivot_row][j] = m[pivot][j];
                m[pivot][j] = t;
            }
            q64 divisor = m[pivot_row][col];
            for (int j = 0; j < width; j++) m[pivot_row][j] = qdiv(m[pivot_row][j], divisor);
            for (int row = 0; row < 3; row++) {
                if (row == pivot_row || !m[row][col].n) continue;
                q64 factor = m[row][col];
                for (int j = 0; j < width; j++) m[row][j] = qsub(m[row][j], qmul(factor, m[pivot_row][j]));
            }
            pivot_row++;
        }
        if (pivot_row != kdim) continue;
        int64_t scale = 1;
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) scale = lcm64(scale, m[row][kdim + col].d);
        s->valid[arm] = 1;
        s->scale[arm] = scale;
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) {
                int64_t value = (int64_t)((i128)m[row][kdim + col].n * (scale / m[row][kdim + col].d));
                if (row < kdim)
                    s->left[9 * arm + 3 * row + col] = value;
                else
                    s->consistency[9 * arm + 3 * (row - kdim) + col] = value;
            }
    }
    free(solvers[gid]);
    solvers[gid] = s;
    s->ready = 1;
    return s;
}

/* Star-arm parameter solutions of target: the first one, or all distinct
 * ones mod 1. */
static int star_solve(int gid, int kdim, const q64 *target, int first, q64 (*out)[3])
{
    const star_solver *s = solver_for(gid, kdim);
    if (!s) return 0;
    int64_t target_values[6];
    for (int axis = 0; axis < 3; axis++) {
        target_values[2 * axis] = target[axis].n;
        target_values[2 * axis + 1] = target[axis].d;
    }
    static int64_t buffer[27 * MAX_ARMS * 3 * 2];
    int capacity = first ? 1 : (s->count > 0 ? 27 * s->count : 1);
    int count = kstar_solutions(first ? 0 : 1, kdim, s->count, s->origin, s->valid, s->left, s->consistency,
                                s->scale, target_values, capacity, buffer);
    if (count < 0) {
        fail(KF_OVERFLOW);
        return 0;
    }
    for (int solution = 0; solution < count; solution++)
        for (int row = 0; row < kdim; row++)
            out[solution][row] = (q64){buffer[(solution * kdim + row) * 2], buffer[(solution * kdim + row) * 2 + 1]};
    return count;
}

static int64_t kparam_from(const q64 *values, int used, int64_t *out)
{
    int64_t den = 1;
    for (int i = 0; i < used; i++) den = lcm64(den, values[i].d);
    for (int i = 0; i < 3; i++) out[i] = i < used ? (int64_t)((i128)values[i].n * (den / values[i].d)) : 0;
    out[3] = den;
    return den;
}

/* The Source parameters of the gid's star arm holding vector, if any. */
static int source_kparam(const sg_context *ctx, int gid, const q64 *vector, int64_t *out)
{
    int kdim = k_dimension(ctx, gid);
    if (kdim <= 0) {
        q64 vectors[MAX_ARMS][3];
        int count = special_vectors(gid, vectors);
        for (int arm = 0; arm < count; arm++)
            if (!key_cmp(vectors[arm], vector)) {
                out[0] = out[1] = out[2] = 0;
                out[3] = 1;
                return 1;
            }
        return 0;
    }
    q64 solutions[1][3];
    if (!star_solve(gid, kdim, vector, 1, solutions)) return 0;
    kparam_from(solutions[0], kdim, out);
    return 1;
}

typedef struct {
    int kdim;
    int64_t periods[3];
    int has_affine;
    q64 affine[4][3];
    int shift_count;
    int shifts[27][3];
    q64 deltas[27][3];
} standard_context;

static void load_standard_context(const sg_context *ctx, int gid, standard_context *sc)
{
    memset(sc, 0, sizeof *sc);
    sc->kdim = k_dimension(ctx, gid);
    if (sc->kdim < 0 || sc->kdim > 3) {
        fail(KF_TABLES);
        return;
    }
    int64_t kslot = at(&T.irr_k, gid - 1);
    int64_t slot = (ctx->lattice - 1) * 27 + kslot - 1;
    for (int index = 0; index < sc->kdim; index++) {
        int64_t start = slot * 16 + 4 * (index + 1);
        int64_t den = at(&T.little_k, start + 3);
        if (den < 0) den = -den;
        int64_t period = 1;
        if (den) {
            for (int axis = 0; axis < 3; axis++) {
                int64_t num = at(&T.little_k, start + axis);
                if (num < 0) num = -num;
                period = lcm64(period, den / (int64_t)gcd128(den, num));
            }
        }
        sc->periods[index] = period > 1 ? period : 1;
    }
    int count;
    int64_t pointer;
    star_location(gid, &count, &pointer);
    if (count > 0) {
        sc->has_affine = 1;
        for (int piece = 0; piece < 4; piece++) {
            q64 raw[3];
            star_piece(pointer, 0, piece, raw);
            to_cinter(ctx, raw, sc->affine[piece]);
        }
    }
    sc->shift_count = 1;
    for (int i = 0; i < sc->kdim; i++) sc->shift_count *= 3;
    for (int s = 0; s < sc->shift_count; s++) {
        int rest = s;
        for (int i = sc->kdim - 1; i >= 0; i--) {
            sc->shifts[s][i] = rest % 3 - 1;
            rest /= 3;
        }
        for (int axis = 0; axis < 3; axis++) {
            q64 acc = qint(0);
            if (sc->has_affine)
                for (int i = 0; i < sc->kdim; i++)
                    acc = qadd(acc, qmul(qint(sc->shifts[s][i] * sc->periods[i]), sc->affine[i + 1][axis]));
            sc->deltas[s][axis] = acc;
        }
    }
}

typedef struct {
    int found, updates;
    q64 key[3];
    q64 solved[3];
    int shift;
} standard_best;

static void standard_param(const standard_context *sc, const standard_best *best, int64_t *out)
{
    q64 signed_values[3];
    for (int i = 0; i < sc->kdim; i++)
        signed_values[i] = qadd(best->solved[i], qint(sc->shifts[best->shift][i] * sc->periods[i]));
    kparam_from(signed_values, sc->kdim, out);
}

/* Candidates of one rotated vector, in Source order; improve when strictly
 * below the running best (or the bound, when given). */
static void standard_scan(int gid, const standard_context *sc, const q64 *rotated, standard_best *best)
{
    static q64 solutions[27 * MAX_ARMS][3];
    int count = sc->kdim <= 0 ? 1 : star_solve(gid, sc->kdim, rotated, 0, solutions);
    for (int solution = 0; solution < count; solution++) {
        q64 solved_cinter[3];
        affine_at((const q64 (*)[3])sc->affine, sc->kdim, solutions[solution], solved_cinter);
        for (int s = 0; s < sc->shift_count; s++) {
            q64 candidate[3];
            for (int axis = 0; axis < 3; axis++) candidate[axis] = qadd(solved_cinter[axis], sc->deltas[s][axis]);
            if (best->found && cinter_cmp(candidate, best->key) >= 0) continue;
            best->found = 1;
            memcpy(best->key, candidate, sizeof candidate);
            for (int i = 0; i < sc->kdim; i++) best->solved[i] = solutions[solution][i];
            best->shift = s;
            best->updates++;
        }
    }
}

/* kvec_standard_: the parameters of the minimal cinter first arm over the
 * point orbit, when it is below the original's. */
static void standard_kparam(const sg_context *ctx, int gid, const int64_t *kparam, int64_t *out)
{
    int64_t key[6] = {ctx->sg, gid, kparam[0], kparam[1], kparam[2], kparam[3]};
    const int64_t *hit = memo_get(&standard_memo, key);
    if (hit) {
        memcpy(out, hit, 4 * sizeof(int64_t));
        return;
    }
    memcpy(out, kparam, 4 * sizeof(int64_t));
    standard_context sc;
    load_standard_context(ctx, gid, &sc);
    q64 original[3];
    if (!bad && sc.has_affine && first_arm_at(gid, kparam, original)) {
        q64 rotated[MAX_OPS][3];
        int rotated_count = 0, canonical = 0;
        for (int j = 0; j < ctx->standard_count; j++) {
            q64 r[3];
            rotate(ctx, ctx->standard_rot[j], original, 1, r);
            int seen = 0;
            for (int i = 0; i < rotated_count && !seen; i++)
                seen = qeq(rotated[i][0], r[0]) && qeq(rotated[i][1], r[1]) && qeq(rotated[i][2], r[2]);
            if (seen) continue;
            memcpy(rotated[rotated_count], r, sizeof r);
            if (key_cmp(r, rotated[canonical]) < 0) canonical = rotated_count;
            rotated_count++;
        }
        int64_t orbit_key[8] = {ctx->sg, gid};
        write_vector(orbit_key + 2, rotated[canonical]);
        int64_t *orbit = memo_get(&orbit_memo, orbit_key);
        if (!orbit) {
            standard_best best = {0};
            for (int i = 0; i < rotated_count; i++) standard_scan(gid, &sc, rotated[i], &best);
            orbit = malloc(11 * sizeof(int64_t));
            if (!orbit) {
                fail(KF_MEMORY);
                return;
            }
            orbit[0] = best.found;
            if (best.found) {
                write_vector(orbit + 1, best.key);
                standard_param(&sc, &best, orbit + 7);
            }
            if (bad || memo_put(&orbit_memo, orbit_key, orbit)) {
                free(orbit);
                fail(KF_MEMORY);
                return;
            }
        }
        if (orbit[0]) {
            q64 best_key[3], current[3];
            for (int axis = 0; axis < 3; axis++) best_key[axis] = (q64){orbit[1 + 2 * axis], orbit[2 + 2 * axis]};
            to_cinter(ctx, original, current);
            if (cinter_cmp(best_key, current) < 0) memcpy(out, orbit + 7, 4 * sizeof(int64_t));
        }
    }
    if (bad) return;
    int64_t *value = malloc(4 * sizeof(int64_t));
    if (!value) {
        fail(KF_MEMORY);
        return;
    }
    memcpy(value, out, 4 * sizeof(int64_t));
    if (memo_put(&standard_memo, key, value)) {
        free(value);
        fail(KF_MEMORY);
    }
}

/* Families whose standardized star holds target: [count, (family, kparam[4])...]. */
static const int64_t *family_candidates(const sg_context *ctx, const q64 *target)
{
    int64_t key[7] = {ctx->sg, target[0].n, target[0].d, target[1].n, target[1].d, target[2].n, target[2].d};
    const int64_t *hit = memo_get(&candidate_memo, key);
    if (hit) return hit;
    int64_t *value = malloc(sizeof(int64_t) * (1 + 5 * MAX_FAMILIES));
    if (!value) {
        fail(KF_MEMORY);
        return 0;
    }
    value[0] = 0;
    for (int f = 0; f < ctx->family_count && !bad; f++) {
        int gid = ctx->family_gid[f];
        int64_t solved[4], kparam[4];
        if (!source_kparam(ctx, gid, target, solved)) continue;
        q64 arms[MAX_ARMS][3];
        int count;
        if (ctx->family_kdim[f] <= 0) {
            kparam[0] = kparam[1] = kparam[2] = 0;
            kparam[3] = 1;
            count = special_vectors(gid, arms);
        } else {
            standard_kparam(ctx, gid, solved, kparam);
            count = star_at(gid, kparam, arms);
        }
        int holds = 0;
        for (int arm = 0; arm < count && !holds; arm++) holds = !key_cmp(arms[arm], target);
        if (!holds) continue;
        int64_t *row = value + 1 + 5 * value[0]++;
        row[0] = f;
        memcpy(row + 1, kparam, 4 * sizeof(int64_t));
    }
    if (bad || memo_put(&candidate_memo, key, value)) {
        free(value);
        fail(KF_MEMORY);
        return 0;
    }
    return value;
}

typedef struct {
    int count;
    const char *const *labels;
    const int64_t *label_lengths;
    const int64_t *kparams;
} preferred_set;

static int is_preferred(const sg_context *ctx, const preferred_set *preferred, int f, const int64_t *kparam)
{
    if (!preferred) return 0;
    for (int i = 0; i < preferred->count; i++) {
        if (!label_equal(preferred->labels[i], preferred->label_lengths[i], ctx->family_label[f],
                         ctx->family_label_length[f]))
            continue;
        const int64_t *target = preferred->kparams + 4 * i;
        int same = 1;
        for (int axis = 0; axis < 3 && same; axis++)
            same = qeq(qmod1(qmake(kparam[axis], kparam[3])), qmod1(qmake(target[axis], target[3])));
        if (same) return 1;
    }
    return 0;
}

typedef struct {
    int family;
    int64_t kparam[4];
    q64 vector[3];
} family_match;

/* The k families of one mesh vector (reduced mod 1), in Source order. */
static int family_matches(const sg_context *ctx, const q64 *vector, const preferred_set *preferred,
                          family_match *out)
{
    family_match found[MAX_FAMILIES];
    int have[MAX_FAMILIES] = {0};
    q64 targets[MAX_OPS + 1][3];
    int target_count = 0;
    for (int j = -1; j < ctx->lattice_count; j++) {
        q64 r[3];
        if (j < 0)
            memcpy(r, vector, sizeof r);
        else
            rotate(ctx, ctx->lattice_rot[j], vector, 1, r);
        int seen = 0;
        for (int i = 0; i < target_count && !seen; i++) seen = !key_cmp(targets[i], r);
        if (seen) continue;
        memcpy(targets[target_count++], r, sizeof r);
        const int64_t *candidates = family_candidates(ctx, r);
        if (!candidates) return 0;
        for (int c = 0; c < candidates[0]; c++) {
            const int64_t *row = candidates + 1 + 5 * c;
            int f = (int)row[0];
            if (have[f]) continue;
            have[f] = 1;
            found[f].family = f;
            memcpy(found[f].kparam, row + 1, 4 * sizeof(int64_t));
            memcpy(found[f].vector, r, sizeof r);
        }
    }
    int best_unpreferred = 2, best_kdim = 0;
    for (int f = 0; f < ctx->family_count; f++) {
        if (!have[f]) continue;
        int unpreferred = !is_preferred(ctx, preferred, f, found[f].kparam);
        int kdim = ctx->family_kdim[f];
        if (unpreferred < best_unpreferred || (unpreferred == best_unpreferred && kdim < best_kdim)) {
            best_unpreferred = unpreferred;
            best_kdim = kdim;
        }
    }
    int count = 0;
    for (int f = 0; f < ctx->family_count; f++) {
        if (!have[f]) continue;
        if ((!is_preferred(ctx, preferred, f, found[f].kparam)) != best_unpreferred || ctx->family_kdim[f] != best_kdim)
            continue;
        family_match match = found[f];
        if (ctx->family_kdim[f] > 0) {
            int gid = ctx->family_gid[f];
            int64_t occurrence[4], standardized[4];
            if (source_kparam(ctx, gid, vector, occurrence)) {
                standard_kparam(ctx, gid, occurrence, standardized);
                q64 arms[MAX_ARMS][3];
                int arm_count = star_at(gid, standardized, arms), holds = 0;
                for (int arm = 0; arm < arm_count && !holds; arm++) {
                    q64 negated[3] = {qneg(arms[arm][0]), qneg(arms[arm][1]), qneg(arms[arm][2])};
                    holds = !key_cmp(arms[arm], vector) || !key_cmp(negated, vector);
                }
                if (holds) memcpy(match.kparam, standardized, sizeof standardized);
            }
        }
        out[count++] = match;
    }
    return bad ? 0 : count;
}

static void lattice_orbit_key(const sg_context *ctx, const q64 *vector, q64 *out)
{
    memcpy(out, vector, 3 * sizeof(q64));
    for (int j = 0; j < ctx->lattice_count; j++) {
        q64 r[3];
        rotate(ctx, ctx->lattice_rot[j], vector, 1, r);
        if (!j || key_cmp(r, out) < 0) memcpy(out, r, sizeof r);
    }
    for (int axis = 0; axis < 3; axis++) out[axis] = qmod1(out[axis]);
}

/* Offsets 0, -g, -2g, ... (mod det) up to the first one inside subgroup. */
static int cyclic_offsets(const int64_t *generator, const int64_t *subgroup, int subgroup_count, int64_t det,
                          int64_t *out, int capacity)
{
    int64_t negative[3], current[3];
    for (int axis = 0; axis < 3; axis++) negative[axis] = pymod(-generator[axis], det);
    int count = 1;
    out[0] = out[1] = out[2] = 0;
    memcpy(current, negative, sizeof current);
    for (;;) {
        int inside = 0;
        for (int i = 0; i < subgroup_count && !inside; i++) inside = !memcmp(subgroup + 3 * i, current, sizeof current);
        if (inside) return count;
        if (count >= capacity) return -1;
        memcpy(out + 3 * count++, current, sizeof current);
        for (int axis = 0; axis < 3; axis++) current[axis] = pymod(current[axis] + negative[axis], det);
    }
}

/* The reciprocal quotient of basis in Source scan order, as numerators over
 * det; returns the point count or a negative error. */
static int64_t reciprocal_mesh(const int64_t *b, int64_t **points, int64_t *det_out)
{
    int64_t adj[9] = {
        b[4] * b[8] - b[5] * b[7], b[2] * b[7] - b[1] * b[8], b[1] * b[5] - b[2] * b[4],
        b[5] * b[6] - b[3] * b[8], b[0] * b[8] - b[2] * b[6], b[2] * b[3] - b[0] * b[5],
        b[3] * b[7] - b[4] * b[6], b[1] * b[6] - b[0] * b[7], b[0] * b[4] - b[1] * b[3],
    };
    int64_t signed_det = b[0] * adj[0] + b[1] * adj[3] + b[2] * adj[6];
    int64_t det = signed_det < 0 ? -signed_det : signed_det;
    *det_out = det;
    *points = 0;
    if (!det) return 0;
    int64_t *result = malloc(sizeof(int64_t) * 3 * (size_t)det);
    int64_t *third = malloc(sizeof(int64_t) * 3 * (size_t)det);
    int64_t *second = malloc(sizeof(int64_t) * 3 * (size_t)det);
    int64_t *first = malloc(sizeof(int64_t) * 3 * (size_t)det);
    int64_t *joined = malloc(sizeof(int64_t) * 3 * (size_t)det);
    if (!result || !third || !second || !first || !joined) {
        free(result);
        free(third);
        free(second);
        free(first);
        free(joined);
        return -KF_MEMORY;
    }
    int64_t sign = signed_det > 0 ? 1 : -1;
    int64_t generators[3][3];
    for (int axis = 0; axis < 3; axis++)
        for (int row = 0; row < 3; row++) generators[axis][row] = pymod(sign * adj[3 * row + axis], det);
    int64_t zero[3] = {0, 0, 0};
    int64_t count = 0;
    int third_count = cyclic_offsets(generators[2], zero, 1, det, third, (int)det);
    int second_count = third_count < 0 ? -1 : cyclic_offsets(generators[1], third, third_count, det, second, (int)det);
    int joined_count = 0;
    if (second_count > 0 && (int64_t)second_count * third_count <= det)
        for (int s = 0; s < second_count; s++)
            for (int t = 0; t < third_count; t++) {
                int64_t v[3];
                for (int axis = 0; axis < 3; axis++) v[axis] = pymod(second[3 * s + axis] + third[3 * t + axis], det);
                int seen = 0;
                for (int i = 0; i < joined_count && !seen; i++) seen = !memcmp(joined + 3 * i, v, sizeof v);
                if (!seen) memcpy(joined + 3 * joined_count++, v, sizeof v);
            }
    int first_count = second_count <= 0 ? -1 : cyclic_offsets(generators[0], joined, joined_count, det, first, (int)det);
    if (first_count < 0 || (int64_t)first_count * second_count * third_count > det) {
        count = -KF_CAPACITY;
    } else {
        for (int f = 0; f < first_count; f++)
            for (int s = 0; s < second_count; s++)
                for (int t = 0; t < third_count; t++) {
                    for (int axis = 0; axis < 3; axis++)
                        result[3 * count + axis] =
                            pymod(first[3 * f + axis] + second[3 * s + axis] + third[3 * t + axis], det);
                    count++;
                }
    }
    free(third);
    free(second);
    free(first);
    free(joined);
    if (count < 0) {
        free(result);
        return count;
    }
    *points = result;
    return count;
}

enum {
    FAMILY_WIDTH = 13,
    OCCURRENCE_WIDTH = 17,
};

typedef struct {
    int mesh_order;
    int family;
    int group;
    int64_t kparam[4];
    q64 vector[3];
    q64 orbit[3];
} occurrence;

typedef struct {
    int family;
    int first_mesh;
    family_match match;
    q64 orbit[3];
} family_group;

static int group_order(const void *a, const void *b)
{
    const family_group *x = a, *y = b;
    if (x->family != y->family) return x->family - y->family;
    return x->first_mesh - y->first_mesh;
}

/* families_out: gid, kparam[4], vector[6], occurrence start, occurrence count;
 * occurrences_out: mesh order, kparam[4], vector[6], orbit key[6].  Returns
 * the family count, or a negative error. */
int k_family_records(int mode, int sg, const int64_t *basis, int preferred_count, const char *const *preferred_labels,
                     const int64_t *preferred_label_lengths, const int64_t *preferred_kparams, int family_capacity,
                     int occurrence_capacity, int64_t *families_out, int64_t *occurrences_out,
                     int64_t *occurrence_count_out)
{
    bad = 0;
    if (load_tables()) return -KF_TABLES;
    sg_context ctx;
    int loaded = load_sg_context(sg, &ctx);
    if (loaded) return -loaded;
    if (ctx.identity_order < 0) return -KF_NO_IDENTITY;
    preferred_set preferred = {preferred_count, preferred_labels, preferred_label_lengths, preferred_kparams};
    int64_t *mesh, det;
    int64_t mesh_count = reciprocal_mesh(basis, &mesh, &det);
    if (mesh_count < 0) return (int)mesh_count;
    occurrence *occurrences = malloc(sizeof(occurrence) * (size_t)(mesh_count * MAX_FAMILIES + 1));
    family_group *groups = malloc(sizeof(family_group) * (size_t)(mesh_count * MAX_FAMILIES + 1));
    if (!occurrences || !groups) {
        free(mesh);
        free(occurrences);
        free(groups);
        return -KF_MEMORY;
    }
    int occurrence_count = 0, group_count = 0;
    for (int64_t order = 0; order < mesh_count && !bad; order++) {
        q64 vector[3];
        for (int axis = 0; axis < 3; axis++) vector[axis] = qmake(mesh[3 * order + axis], det);
        family_match matches[MAX_FAMILIES];
        int match_count = family_matches(&ctx, vector, mode == KF_DYNAMIC ? &preferred : 0, matches);
        if (!match_count) continue;
        q64 orbit[3];
        lattice_orbit_key(&ctx, vector, orbit);
        int used = mode == KF_DYNAMIC ? match_count : 1;
        for (int m = 0; m < used; m++) {
            int group = -1;
            for (int g = 0; g < group_count && group < 0; g++) {
                if (mode == KF_DYNAMIC && groups[g].family != matches[m].family) continue;
                if (!key_cmp(groups[g].orbit, orbit)) group = g;
            }
            if (group < 0) {
                group = group_count++;
                groups[group].family = matches[m].family;
                groups[group].first_mesh = (int)order;
                groups[group].match = matches[m];
                memcpy(groups[group].orbit, orbit, sizeof orbit);
            } else if (mode == KF_DISPLAY) {
                continue;
            }
            occurrence *o = &occurrences[occurrence_count++];
            o->mesh_order = (int)order;
            o->family = matches[m].family;
            o->group = group;
            memcpy(o->kparam, matches[m].kparam, sizeof o->kparam);
            memcpy(o->vector, vector, sizeof vector);
            memcpy(o->orbit, orbit, sizeof orbit);
        }
    }
    free(mesh);
    int status = bad;
    int *order = malloc(sizeof(int) * (size_t)(group_count + 1));
    if (!status && !order) status = KF_MEMORY;
    int family_count = 0;
    int64_t written = 0;
    if (!status) {
        {
            for (int g = 0; g < group_count; g++) order[g] = g;
            for (int i = 1; i < group_count; i++) {
                int g = order[i], j = i - 1;
                while (j >= 0 && group_order(&groups[order[j]], &groups[g]) > 0) {
                    order[j + 1] = order[j];
                    j--;
                }
                order[j + 1] = g;
            }
            unsigned char *taken = calloc((size_t)occurrence_count + 1, 1);
            if (!taken) status = KF_MEMORY;
            for (int i = 0; i < group_count && !status; i++) {
                const family_group *g = &groups[order[i]];
                int gid = ctx.family_gid[g->family];
                for (int pass = 0; pass < 2 && !status; pass++)
                    for (int o = 0; o < occurrence_count && !status; o++) {
                        if (taken[o] || occurrences[o].group != order[i]) continue;
                        int selected = !memcmp(occurrences[o].kparam, g->match.kparam, sizeof g->match.kparam);
                        if (pass == 0 && !selected) continue;
                        if (family_count >= family_capacity) {
                            status = KF_CAPACITY;
                            break;
                        }
                        int64_t *family = families_out + FAMILY_WIDTH * family_count++;
                        family[0] = gid;
                        memcpy(family + 1, occurrences[o].kparam, 4 * sizeof(int64_t));
                        q64 vector[3];
                        if (mode == KF_DISPLAY) {
                            if (ctx.family_kdim[g->family] <= 0) {
                                q64 specials[MAX_ARMS][3];
                                if (!special_vectors(gid, specials)) {
                                    status = bad ? bad : KF_NO_SPECIAL_VECTOR;
                                    break;
                                }
                                memcpy(vector, specials[0], sizeof vector);
                            } else if (!first_arm_at(gid, g->match.kparam, vector)) {
                                memcpy(vector, occurrences[o].vector, sizeof vector);
                            }
                        } else {
                            memcpy(vector, selected ? g->match.vector : occurrences[o].vector, sizeof vector);
                        }
                        write_vector(family + 5, vector);
                        family[11] = written;
                        for (int p = o; p < occurrence_count; p++) {
                            if (taken[p] || occurrences[p].group != order[i]
                                || memcmp(occurrences[p].kparam, occurrences[o].kparam, sizeof occurrences[o].kparam))
                                continue;
                            if (written >= occurrence_capacity) {
                                status = KF_CAPACITY;
                                break;
                            }
                            taken[p] = 1;
                            int64_t *row = occurrences_out + OCCURRENCE_WIDTH * written++;
                            row[0] = occurrences[p].mesh_order;
                            memcpy(row + 1, occurrences[p].kparam, 4 * sizeof(int64_t));
                            write_vector(row + 5, occurrences[p].vector);
                            write_vector(row + 11, occurrences[p].orbit);
                        }
                        family[12] = written - family[11];
                        if (bad) status = bad;
                    }
            }
            free(taken);
        }
    }
    free(order);
    free(occurrences);
    free(groups);
    if (status) return -status;
    *occurrence_count_out = written;
    return family_count;
}

/* The Source kparam of the gid's star arm holding vector: 1 and out[4], 0
 * when none holds it, or a negative error. */
int k_family_source_kparam(int sg, int gid, const int64_t *vector, int64_t *out)
{
    bad = 0;
    if (load_tables()) return -KF_TABLES;
    sg_context ctx;
    int loaded = load_sg_context(sg, &ctx);
    if (loaded) return -loaded;
    q64 v[3];
    for (int axis = 0; axis < 3; axis++) v[axis] = qmake(vector[2 * axis], vector[2 * axis + 1]);
    int found = source_kparam(&ctx, gid, v, out);
    return bad ? -bad : found;
}

/* kvec_standard_'s full selection: out = standard kparam[4], selected source
 * order, selected point op. */
int k_family_standard_selection(int sg, int gid, const int64_t *kparam, int64_t *out)
{
    bad = 0;
    if (load_tables()) return -KF_TABLES;
    sg_context ctx;
    int loaded = load_sg_context(sg, &ctx);
    if (loaded) return -loaded;
    if (ctx.identity_order < 0) return -KF_NO_IDENTITY;
    memcpy(out, kparam, 4 * sizeof(int64_t));
    out[4] = ctx.identity_order;
    out[5] = ctx.standard_op[ctx.identity_order];
    standard_context sc;
    load_standard_context(&ctx, gid, &sc);
    q64 original[3];
    if (bad || !sc.has_affine || !first_arm_at(gid, kparam, original)) return bad ? -bad : 0;
    standard_best best = {0};
    to_cinter(&ctx, original, best.key);
    best.found = 1;
    q64 seen[MAX_OPS][3];
    int seen_count = 0;
    for (int j = 0; j < ctx.standard_count && !bad; j++) {
        q64 r[3];
        rotate(&ctx, ctx.standard_rot[j], original, 1, r);
        int repeated = 0;
        for (int i = 0; i < seen_count && !repeated; i++)
            repeated = qeq(seen[i][0], r[0]) && qeq(seen[i][1], r[1]) && qeq(seen[i][2], r[2]);
        if (repeated) continue;
        memcpy(seen[seen_count++], r, sizeof r);
        int updates = best.updates;
        standard_scan(gid, &sc, r, &best);
        if (best.updates != updates) {
            standard_param(&sc, &best, out);
            out[4] = j;
            out[5] = ctx.standard_op[j];
        }
    }
    return bad ? -bad : 0;
}
