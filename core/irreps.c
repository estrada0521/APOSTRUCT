/* The get_irreps_ context of one little irrep at one k parameter: its star
 * arms, the generated operations of its space group and the decoded Source
 * matrices at their positions (SourceData._native_irrep_context). */

#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "irreps.h"
#include "tables.h"

int little_sparse_context(int gid, int dim, int element_count, double *sparse, int *available);

typedef struct {
    int ready, dim, subtract, element_count, arms;
    int elements[GEO_MAX_OPS * 5];
    int available[GEO_MAX_OPS];
    double *sparse;
    int64_t key[4];
    int key_count;
    int64_t k[48 * 6];
} gid_context;

static struct {
    section_view sg, dim, kslot, star_count, star_pointer, k_dim, star_special, star;
    int ready;
    gid_context *contexts;
    int64_t count;
} L;

static int load(void)
{
    if (L.ready) return 0;
    if (geo_load()) return IR_TABLES;
    struct {
        const char *name;
        section_view *view;
    } items[] = {
        {"little.little_irr_space_group", &L.sg},     {"little.little_irr_full_dim", &L.dim},
        {"little.little_irr_k", &L.kslot},            {"little.little_k_star_count", &L.star_count},
        {"little.little_k_star_ml_pointer", &L.star_pointer}, {"little.little_k_dim", &L.k_dim},
        {"little.little_k_star_special", &L.star_special},   {"little.little_k_star", &L.star},
    };
    for (size_t i = 0; i < sizeof items / sizeof items[0]; i++)
        if (section_find(items[i].name, items[i].view) || !section_is_int(items[i].view)) return IR_TABLES;
    L.count = section_count(&L.sg);
    L.contexts = calloc((size_t)L.count + 1, sizeof *L.contexts);
    if (!L.contexts) return IR_TABLES;
    L.ready = 1;
    return 0;
}

int irrep_space_group(int gid)
{
    if (load() || gid < 1 || gid > L.count) return 0;
    return (int)section_int(&L.sg, gid - 1);
}

static int64_t at(const section_view *view, int64_t index, int *bad)
{
    if (index < 0 || index >= section_count(view)) {
        *bad = 1;
        return 0;
    }
    return section_int(view, index);
}

/* _little_k_star_vectors_uncached as (numerator, denominator) per axis. */
static int star_vectors(int gid, int param_count, const int64_t *kparam, gid_context *g)
{
    int bad = 0;
    int sg = (int)at(&L.sg, gid - 1, &bad), kslot = (int)at(&L.kslot, gid - 1, &bad);
    int64_t slot = (int64_t)(sg - 1) * 27 + kslot - 1;
    int count = (int)at(&L.star_count, slot, &bad);
    int64_t pointer = at(&L.star_pointer, slot, &bad);
    rat params[3] = {rat_int(0), rat_int(0), rat_int(0)};
    if (param_count >= 4 && kparam[3])
        for (int i = 0; i < 3; i++) params[i] = rat_make(kparam[i], kparam[3]);
    int64_t lattice_slot = (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1;
    int special = at(&L.k_dim, lattice_slot, &bad) <= 0;
    if (bad || count < 0 || count > 48) return IR_TABLES;
    rat_overflow = 0;
    for (int arm = 0; arm < count; arm++) {
        rat out[3];
        if (special) {
            int64_t start = 4 * (pointer - 1 + arm), den = at(&L.star_special, start + 3, &bad);
            for (int axis = 0; axis < 3; axis++) out[axis] = rat_make(at(&L.star_special, start + axis, &bad), den ? den : 1);
        } else {
            int64_t start = 16 * (pointer - 1 + arm);
            rat pieces[4][3];
            for (int piece = 0; piece < 4; piece++) {
                int64_t den = at(&L.star, start + 4 * piece + 3, &bad);
                for (int axis = 0; axis < 3; axis++)
                    pieces[piece][axis] = den ? rat_make(at(&L.star, start + 4 * piece + axis, &bad), den) : rat_int(0);
            }
            for (int axis = 0; axis < 3; axis++) {
                out[axis] = pieces[0][axis];
                for (int p = 0; p < 3; p++) out[axis] = rat_add(out[axis], rat_mul(params[p], pieces[p + 1][axis]));
            }
        }
        for (int axis = 0; axis < 3; axis++) {
            if (!rat_fits64(out[axis])) rat_overflow = 1;
            g->k[6 * arm + 2 * axis] = (int64_t)out[axis].n;
            g->k[6 * arm + 2 * axis + 1] = (int64_t)out[axis].d;
        }
    }
    if (bad) return IR_TABLES;
    if (rat_overflow) return IR_RANGE;
    g->arms = count;
    return 0;
}

int irrep_context_for(int gid, int param_count, const int64_t *kparam, irrep_context *out)
{
    if (load()) return IR_TABLES;
    if (gid < 1 || gid > L.count) return IR_NO_TABLES;
    gid_context *g = &L.contexts[gid];
    if (!g->ready) {
        int bad = 0;
        int sg = (int)at(&L.sg, gid - 1, &bad), kslot = (int)at(&L.kslot, gid - 1, &bad);
        g->dim = (int)at(&L.dim, gid - 1, &bad);
        int64_t records[GEO_MAX_OPS][5];
        g->element_count = geo_sg_records(sg, records);
        if (bad || geo_take_error()) return IR_TABLES;
        if (g->element_count != geo_point_group_order(geo_sg_point_group(sg)) || g->dim <= 0 || g->dim > 48)
            return IR_NO_TABLES;
        for (int i = 0; i < g->element_count; i++)
            for (int j = 0; j < 5; j++) g->elements[5 * i + j] = (int)records[i][j];
        g->sparse = malloc(sizeof(double) * (size_t)g->element_count * (size_t)g->dim * (size_t)g->dim);
        if (!g->sparse) return IR_TABLES;
        little_sparse_context(gid, g->dim, g->element_count, g->sparse, g->available);
        g->subtract = at(&L.k_dim, (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1, &bad) <= 0;
        g->key_count = -1;
        g->ready = 1;
    }
    int count = param_count < 4 ? param_count : 4;
    if (count != g->key_count || memcmp(g->key, kparam, sizeof(int64_t) * (size_t)count)) {
        int status = star_vectors(gid, param_count, kparam, g);
        if (status) {
            g->key_count = -1;
            return status;
        }
        memcpy(g->key, kparam, sizeof(int64_t) * (size_t)count);
        g->key_count = count;
    }
    out->dim = g->dim;
    out->arms = g->arms;
    out->k = g->k;
    out->subtract = g->subtract;
    out->element_count = g->element_count;
    out->elements = g->elements;
    out->sparse = g->sparse;
    out->available = g->available;
    out->identity_op = g->elements[4];
    return 0;
}

/* the context of gid as it stands (its last k parameters), without changing them; nonzero when none is ready */
int irrep_context_current(int gid, irrep_context *out)
{
    if (load() || gid < 1 || gid > L.count) return 1;
    gid_context *g = &L.contexts[gid];
    if (!g->ready || g->key_count < 0) return 1;
    out->dim = g->dim;
    out->arms = g->arms;
    out->k = g->k;
    out->subtract = g->subtract;
    out->element_count = g->element_count;
    out->elements = g->elements;
    out->sparse = g->sparse;
    out->available = g->available;
    out->identity_op = g->elements[4];
    return 0;
}

int bridge_weight_matrices(int dim, int arms, const int64_t *k, int subtract_generated, int element_count,
                           const int *elements, int op_limit, const int *canonical_position, const double *sparse,
                           const int *available, int record_count, const int *records, double *out);

/* get_irrep4_ matrices for a mode case: parametric k (params: three
 * rationals as (numerator, denominator), param_count 3) substituted into
 * little_k_star, else the fixed little_k_star_special arms; a record whose
 * translation is not generated falls back to the first generated operation
 * of its point op.  Returns a bridge-weight status, or 100 + IR_* for the
 * context (IR_RANGE also for bad k-star pointers). */
static int case_irrep_stack_uncached(int gid, int param_count, const int64_t *params, int record_count,
                                     const int *records, int64_t capacity, double *out);

/* The last few stacks, keyed by their exact inputs. */
#define STACK_MEMO 8

static struct {
    int gid, param_count, record_count, status;
    int64_t params[6];
    int *records;
    double *stack;
    size_t size;
} memo[STACK_MEMO];
static int memo_next;

int case_irrep_stack(int gid, int param_count, const int64_t *params, int record_count, const int *records,
                     int64_t capacity, double *out)
{
    int words = param_count == 3 ? 6 : 0;
    for (int i = 0; i < STACK_MEMO; i++)
        if (memo[i].stack && memo[i].gid == gid && memo[i].param_count == param_count &&
            memo[i].record_count == record_count && !memcmp(memo[i].params, params, sizeof(int64_t) * (size_t)words) &&
            !memcmp(memo[i].records, records, sizeof(int) * 5 * (size_t)record_count)) {
            if ((int64_t)memo[i].size > capacity) return 100 + IR_RANGE;
            memcpy(out, memo[i].stack, sizeof(double) * memo[i].size);
            return memo[i].status;
        }
    int status = case_irrep_stack_uncached(gid, param_count, params, record_count, records, capacity, out);
    if (status) return status;
    irrep_context c;
    const int64_t none[1] = {0};
    if (irrep_context_for(gid, 0, none, &c)) return status;
    size_t size = (size_t)record_count * c.dim * c.dim;
    int slot = memo_next;
    memo_next = (memo_next + 1) % STACK_MEMO;
    free(memo[slot].records);
    free(memo[slot].stack);
    memo[slot].records = malloc(sizeof(int) * 5 * (size_t)(record_count + 1));
    memo[slot].stack = malloc(sizeof(double) * (size + 1));
    if (!memo[slot].records || !memo[slot].stack) {
        free(memo[slot].records);
        free(memo[slot].stack);
        memo[slot].records = 0;
        memo[slot].stack = 0;
        return status;
    }
    memcpy(memo[slot].records, records, sizeof(int) * 5 * (size_t)record_count);
    memcpy(memo[slot].stack, out, sizeof(double) * size);
    memcpy(memo[slot].params, params, sizeof(int64_t) * (size_t)words);
    memo[slot].gid = gid;
    memo[slot].param_count = param_count;
    memo[slot].record_count = record_count;
    memo[slot].size = size;
    memo[slot].status = status;
    return status;
}

static int case_irrep_stack_uncached(int gid, int param_count, const int64_t *params, int record_count,
                                     const int *records, int64_t capacity, double *out)
{
    irrep_context c;
    const int64_t none[1] = {0};
    int status = irrep_context_for(gid, 0, none, &c);
    if (status) return 100 + status;
    if ((int64_t)record_count * c.dim * c.dim > capacity) return 100 + IR_RANGE;
    int bad = 0;
    int sg = (int)at(&L.sg, gid - 1, &bad), kslot = (int)at(&L.kslot, gid - 1, &bad);
    int64_t slot = (int64_t)(sg - 1) * 27 + kslot - 1;
    int count = (int)at(&L.star_count, slot, &bad);
    int64_t pointer = at(&L.star_pointer, slot, &bad);
    if (bad || count < 0 || count > 48 || (count && pointer <= 0)) return 100 + IR_RANGE;
    int64_t k[48 * 6];
    int subtract = 1;
    rat_overflow = 0;
    for (int arm = 0; arm < count; arm++) {
        rat out_k[3];
        if (param_count == 3) {
            int64_t start = 16 * (pointer - 1 + arm);
            if (start + 16 > section_count(&L.star)) return 100 + IR_RANGE;
            rat pieces[4][3];
            for (int piece = 0; piece < 4; piece++) {
                int64_t den = at(&L.star, start + 4 * piece + 3, &bad);
                for (int axis = 0; axis < 3; axis++)
                    pieces[piece][axis] = den ? rat_make(at(&L.star, start + 4 * piece + axis, &bad), den) : rat_int(0);
            }
            for (int axis = 0; axis < 3; axis++) {
                out_k[axis] = pieces[0][axis];
                for (int p = 0; p < 3; p++)
                    out_k[axis] = rat_add(out_k[axis], rat_mul(rat_make(params[2 * p], params[2 * p + 1]), pieces[p + 1][axis]));
            }
        } else {
            int64_t start = 4 * (pointer - 1 + arm);
            if (start + 4 > section_count(&L.star_special)) return 100 + IR_RANGE;
            int64_t den = at(&L.star_special, start + 3, &bad);
            if (!den) return 100 + IR_RANGE;
            for (int axis = 0; axis < 3; axis++) out_k[axis] = rat_make(at(&L.star_special, start + axis, &bad), den);
        }
        for (int axis = 0; axis < 3; axis++) {
            if (!rat_fits64(out_k[axis])) rat_overflow = 1;
            k[6 * arm + 2 * axis] = (int64_t)out_k[axis].n;
            k[6 * arm + 2 * axis + 1] = (int64_t)out_k[axis].d;
        }
    }
    if (param_count == 3) subtract = at(&L.k_dim, (int64_t)(geo_sg_lattice(sg) - 1) * 27 + kslot - 1, &bad) <= 0;
    if (bad || rat_overflow) return 100 + IR_RANGE;
    int op_limit = 0;
    for (int i = 0; i < c.element_count; i++)
        if (c.elements[5 * i + 4] > op_limit) op_limit = c.elements[5 * i + 4];
    int canonical[GEO_MAX_OPS * 4 + 1];
    if (op_limit > GEO_MAX_OPS * 4) return 100 + IR_RANGE;
    memset(canonical, 0, sizeof(int) * (size_t)(op_limit + 1));
    for (int i = c.element_count - 1; i >= 0; i--) canonical[c.elements[5 * i + 4]] = i + 1;
    return bridge_weight_matrices(c.dim, count, k, subtract, c.element_count, c.elements, op_limit, canonical,
                                  c.sparse, c.available, record_count, records, out);
}
