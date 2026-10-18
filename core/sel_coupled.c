/* isotropy.coupled.coupled_opd_rows: relative domains of every selected
 * order-parameter subspace, their common stabilizers and identified
 * subgroups. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "magnetic.h"
#include "sel.h"
#include "tables.h"

int orderparam_to_subgroup_rows(int gid, int param_count, const int64_t *kparam, const double *values, int64_t value_count, int row_count,
                                int candidate_count, const int64_t *candidates, const int *candidate_matrix_ops, const int *candidate_signs,
                                int max_iter, int64_t *basis_out, int64_t *kept_out, int *kept_count);
int coupled_newlat_step(const int64_t *slot, const double *row, int has_previous, const int64_t *previous, int64_t *basis_out,
                        int *code_out);
int subgroup_intersection(int sg, int group, const int64_t *basis, int subgroup_count, const int64_t *bases, const int *operation_counts,
                          const int64_t *operations, int64_t capacity, int64_t *out);
int subgroup_kernel_records(int sg, int group, const int64_t *basis, int64_t capacity, int64_t *out);
int same_subgroup(const int64_t *left_basis, int left_count, const int64_t *left, const int64_t *right_basis, int right_count,
                  const int64_t *right);
int mapped_subgroup_records(int parent_sg, int child_sg, const int64_t *basis, const int64_t *origin, int64_t capacity, int64_t *out);
int subgroup_domain_records(int sg, int group, const int64_t *basis, int subgroup_count, const int64_t *subgroup_operations,
                            int64_t capacity, int64_t *out);
int coupled_domain_orbit(int gid, int param_count, const int64_t *kparam, int mode, int free_rows, const double *values, int candidate_count,
                         const int64_t *candidates, int op_count, const int64_t *op_records, int numbered_count,
                         const int64_t *numbered_records, int use_numbered, int generator_count, const int64_t *generator_records);
int64_t coupled_domain_orbit_ints(const int64_t **data);
int64_t coupled_domain_orbit_doubles(const double **data);
int relative_domain_tuples(int slot_count, const int *sizes, int op_count, const int32_t *actions, int64_t capacity, int32_t *out);
void coupled_tuples_reset(void);
int coupled_tuples(int sg, int group, int slot_count, const int64_t *slots, const int64_t *kernel_basis, const double *values,
                   const int64_t *bases, const int *op_counts, const int64_t *ops, const int64_t *orders, int magnetic, int mag_count,
                   const int *mag2nonmag, int tuple_count, const int32_t *tuples);
int64_t coupled_tuples_ints(const int64_t **data);
void coupled_tuples_failure(int64_t *out);
int id_subgroup_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t *out);
int id_subgroup_magnetic_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t capacity,
                                  int64_t *out);
int magnetic_records_iso(int group, int capacity, int64_t *records_out);

typedef int64_t rec[5];

/* Every allocation of one request, released together. */
static struct {
    void **items;
    size_t used, cap;
} pool;

static void *take(size_t bytes)
{
    void *p = calloc(1, bytes ? bytes : 1);
    if (!p) return 0;
    if (pool.used == pool.cap) {
        size_t cap = pool.cap ? 2 * pool.cap : 256;
        void **grown = realloc(pool.items, cap * sizeof *grown);
        if (!grown) {
            free(p);
            return 0;
        }
        pool.items = grown;
        pool.cap = cap;
    }
    pool.items[pool.used++] = p;
    return p;
}

static void pool_release(void)
{
    for (size_t i = 0; i < pool.used; i++) free(pool.items[i]);
    pool.used = 0;
}

#define TAKE(ptr, count)                                                      \
    do {                                                                      \
        (ptr) = take(sizeof *(ptr) * (size_t)(count));                        \
        if (!(ptr)) return md_fail(MDE_MEMORY, "out of memory");              \
    } while (0)

static int kernel_error(int kind, int status)
{
    md_err.kind = kind;
    md_err.status = status;
    md_err.message[0] = 0;
    return kind;
}

static struct {
    section_view irr_sg, point_op, mag2nonmag;
    int ready;
} T;

typedef struct {
    int64_t basis[9];
    int n;
    rec *ops;
} subgroup_t;

typedef struct {
    const jvalue *row, *iso;
    int dynamic;
    int free, nvalues;
    double *values;
    int block_ready;
    subgroup_t sub;
    rec *domain_ops;
    int ndomain_ops, source_order;
    int has_stored, nstored;
    rec *stored;
    int child;
    int64_t sbasis[9], sorigin[4];
    int orbit_ready, ndomains, ngen;
    double **domains;
    int64_t *numbers;
    int32_t *actions;
    subgroup_t *dsubs;
} srow_t;

typedef struct {
    const jvalue *slot;
    int gid, gid_sg, full_dim, little_dim, magnetic, mode, nk;
    int64_t kparam[6];
    sel_params active;
    int has_kparams;
    int nrows;
    srow_t *rows;
    int64_t kernel_basis[9];
    int ncand;
    rec *cand;
    int *cops, *csigns;
} cslot_t;

static struct {
    int sg, group, magnetic_request, mixed_dynamic, nslots;
    cslot_t *slots;
    rec *kernel_ops, *gens, *parent_ops;
    int nkernel_ops, ngens, nparent_ops;
} R;

static const int64_t IDENTITY[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

static int64_t det3(const int64_t *m)
{
    return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
}

static int64_t abs64(int64_t v) { return v < 0 ? -v : v; }

/* coupled._number */
static int number_of(const jvalue *v, double *out)
{
    if (jis(v, J_INT)) {
        *out = (double)v->u.i;
        return 0;
    }
    if (jis(v, J_FLOAT)) {
        *out = v->u.f;
        return 0;
    }
    jbuf text = {0};
    jb_reset(&text);
    py_str_of(&text, v);
    const char *a = text.data, *e = text.data + text.used;
    while (a < e && sel_is_space(*a)) a++;
    while (e > a && sel_is_space(e[-1])) e--;
    rat value;
    int status = 0;
    if (!py_fraction_of_text(a, e - a, &value) && !py_float_of_fraction(value, out)) status = 0;
    else if (py_float_of_text(a, e - a, out)) status = md_fail(MDE_VALUE, "could not convert string to float: '%.*s'", (int)(e - a), a);
    jb_free(&text);
    return status;
}

static int int_of(const jvalue *v, int64_t *out)
{
    if (!v || jis(v, J_NULL)) {
        *out = 0;
        return 0;
    }
    return py_int_of(v, out) ? md_fail(MDE_TYPE, "int() argument must be a string, a bytes-like object or a real number") : 0;
}

static const jvalue *member(const jvalue *object, const char *key) { return jis(object, J_OBJECT) ? jget(object, key) : 0; }

static int64_t member_int(const jvalue *object, const char *key)
{
    int64_t v = 0;
    const jvalue *m = member(object, key);
    if (jtruthy(m)) py_int_of(m, &v);
    return v;
}

static void put_str_member(jbuf *b, const jvalue *v)
{
    if (!jtruthy(v)) return;
    py_str_of(b, v);
}

/* ------------------------------------------------------------ kernels */

static int orderparam_subgroup(const cslot_t *s, const double *values, int nvalues, int free, subgroup_t *out)
{
    rec *kept;
    TAKE(kept, s->ncand + 1);
    int count = 0;
    int code = orderparam_to_subgroup_rows(s->gid, s->nk, s->kparam, values, nvalues, free, s->ncand, &s->cand[0][0], s->cops, s->csigns,
                                           256, out->basis, &kept[0][0], &count);
    if (code) return kernel_error(MDE_ORDERPARAM, code);
    out->n = count;
    out->ops = kept;
    return 0;
}

static int subgroup_ops_status(int status) { return status < 0 ? kernel_error(MDE_SUBGROUP_OPS, -status) : 0; }

static int combined_basis(cslot_t *const *slots, int n, const int *frees, double *const *values, const int *nvalues, int64_t *basis)
{
    subgroup_t *subs;
    TAKE(subs, n);
    for (int i = 0; i < n; i++) MD_TRY(orderparam_subgroup(slots[i], values[i], nvalues[i], frees[i], subs + i));
    int has = 0;
    int64_t current[9];
    for (int i = 0; i < n; i++) {
        const cslot_t *s = slots[i];
        int64_t words[9] = {s->gid, s->nk, s->kparam[0], s->kparam[1], s->kparam[2], s->kparam[3], s->kparam[4], s->kparam[5], s->full_dim};
        double row[48];
        for (int r = 0; r < frees[i]; r++) {
            for (int c = 0; c < s->full_dim && c < 48; c++) row[c] = values[i][r * 48 + c];
            int64_t next[9];
            int code = 0;
            int kind = coupled_newlat_step(words, row, has, current, next, &code);
            if (kind == 1) return kernel_error(MDE_ORDERPARAM, code);
            if (kind == 2) return md_fail(MDE_VALUE, "singular basis");
            if (kind == 5) return md_fail(MDE_MEMORY, "out of memory");
            if (kind) return md_fail(MDE_OVERFLOW, "Python int too large to convert to C long");
            memcpy(current, next, sizeof current);
            has = 1;
        }
    }
    memcpy(basis, has ? current : IDENTITY, sizeof current);
    int64_t *bases;
    int *counts, total = 0;
    TAKE(bases, 9 * n + 9);
    TAKE(counts, n + 1);
    for (int i = 0; i < n; i++) total += subs[i].n;
    rec *ops, *out;
    TAKE(ops, total + 1);
    TAKE(out, 192);
    for (int i = 0, at = 0; i < n; i++) {
        memcpy(bases + 9 * i, subs[i].basis, sizeof subs[i].basis);
        counts[i] = subs[i].n;
        memcpy(ops + at, subs[i].ops, sizeof(rec) * (size_t)subs[i].n);
        at += subs[i].n;
    }
    return subgroup_ops_status(subgroup_intersection(R.sg, R.group, basis, n, bases, counts, &ops[0][0], 192, &out[0][0]));
}

static int kernel_records(const int64_t *basis, int group, rec **out, int *count)
{
    int64_t det = abs64(det3(basis));
    int64_t capacity = 192 * (det > 1 ? det : 1);
    TAKE(*out, capacity);
    int status = subgroup_kernel_records(R.sg, group, basis, capacity, &(*out)[0][0]);
    MD_TRY(subgroup_ops_status(status));
    *count = status;
    return 0;
}

int exact_refine_rows(int rows, int cols, double *values);

/* coupled._orderparam */
static int row_orderparam(const cslot_t *s, srow_t *r)
{
    const jvalue *direction = member(r->row, "direction");
    const jvalue *vectors = 0;
    if (r->dynamic) vectors = member(r->iso, "source_display_rows");
    else vectors = member(r->iso, "source_numeric_rows");
    if (!jtruthy(vectors)) vectors = member(r->iso, "source_display_rows");
    if (!jtruthy(vectors)) vectors = member(direction, "vectors");
    int free = jis(vectors, J_ARRAY) ? (int)vectors->count : 0;
    r->free = free;
    r->nvalues = (free > 1 ? free : 1) * 48;
    TAKE(r->values, r->nvalues);
    for (int row = 0; row < free; row++) {
        const jvalue *values = jat(vectors, row);
        int64_t n = jis(values, J_ARRAY) ? values->count : 0;
        for (int64_t c = 0; c < n && c < s->full_dim; c++) MD_TRY(number_of(jat(values, c), r->values + row * 48 + c));
    }
    /* six-decimal dynamic rows: their exact values */
    int dim = s->full_dim;
    if (free > 0 && dim > 0 && dim <= 48) {
        double *compact;
        TAKE(compact, free * dim + 1);
        for (int row = 0; row < free; row++)
            for (int c = 0; c < dim; c++) compact[row * dim + c] = r->values[row * 48 + c];
        if (exact_refine_rows(free, dim, compact))
            for (int row = 0; row < free; row++)
                for (int c = 0; c < dim; c++) r->values[row * 48 + c] = compact[row * dim + c];
    }
    return 0;
}

static int records_of(const jvalue *values, int64_t *out, int n)
{
    if (!jis(values, J_ARRAY)) return 0;
    int count = 0;
    for (int64_t i = 0; i < values->count; i++) {
        int64_t v;
        MD_TRY(int_of(jat(values, i), &v));
        if (count < n) out[count] = v;
        count++;
    }
    return count;
}

/* coupled_opd_rows._slot_domain_block */
static int slot_block(const cslot_t *s, srow_t *r)
{
    if (r->block_ready) return 0;
    MD_TRY(orderparam_subgroup(s, r->values, r->nvalues, r->free, &r->sub));
    int64_t sbasis[16], sorigin[16];
    int nb = records_of(member(r->iso, "source_basis_values"), sbasis, 16);
    int no = records_of(member(r->iso, "source_origin_values"), sorigin, 16);
    if (nb < 0 || no < 0) return md_err.kind;
    const jvalue *subgroup = member(r->iso, "subgroup");
    int64_t child = member_int(subgroup, s->magnetic ? "ordinary_number" : "number");
    int use_source = (R.magnetic_request || R.mixed_dynamic) && nb == 9;
    const int64_t *numbering = use_source ? sbasis : r->sub.basis;
    int binary = 0;
    int64_t det_basis = abs64(det3(nb == 9 ? sbasis : IDENTITY));
    rec *mapped;
    TAKE(mapped, 192);
    if (!R.magnetic_request && child > 0 && nb == 9 && no == 4) {
        int status = mapped_subgroup_records(R.sg, (int)child, sbasis, sorigin, 192, &mapped[0][0]);
        if (status >= 0) {
            int64_t capacity = 192 * (det_basis > 1 ? det_basis : 1);
            rec *domains;
            TAKE(domains, capacity);
            int count = subgroup_domain_records(R.sg, 0, sbasis, status, &mapped[0][0], capacity, &domains[0][0]);
            if (count >= 0) {
                r->domain_ops = domains;
                r->ndomain_ops = count;
                numbering = sbasis;
                binary = 1;
            } else if (count == -4 || count == -6) {
                return kernel_error(MDE_SUBGROUP_OPS, -count);
            }
        } else if (status == -4 || status == -6) {
            return kernel_error(MDE_SUBGROUP_OPS, -status);
        }
    }
    if (!binary) {
        int64_t det = abs64(det3(numbering));
        int64_t capacity = 192 * (det > 1 ? det : 1);
        rec *domains;
        TAKE(domains, capacity);
        int count = subgroup_domain_records(R.sg, R.magnetic_request ? R.group : 0, numbering, r->sub.n, &r->sub.ops[0][0], capacity,
                                            &domains[0][0]);
        MD_TRY(subgroup_ops_status(count));
        r->domain_ops = domains;
        r->ndomain_ops = count;
        if (R.magnetic_request) binary = 1;
    }
    r->source_order = binary || use_source;
    r->has_stored = 0;
    if (child > 0 && nb == 9 && no == 4) {
        rec *stored;
        TAKE(stored, 192);
        int status = mapped_subgroup_records(R.sg, (int)child, sbasis, sorigin, 192, &stored[0][0]);
        if (status >= 0) {
            r->has_stored = 1;
            r->stored = stored;
            r->nstored = status;
            r->child = (int)child;
            memcpy(r->sbasis, sbasis, sizeof r->sbasis);
            memcpy(r->sorigin, sorigin, sizeof r->sorigin);
        } else if (status == -4 || status == -6) {
            return kernel_error(MDE_SUBGROUP_OPS, -status);
        }
    }
    r->block_ready = 1;
    return 0;
}

static int slot_orbit(const cslot_t *s, srow_t *r)
{
    if (r->orbit_ready) return 0;
    int dim = s->full_dim;
    double *compact;
    TAKE(compact, r->free * dim + 1);
    for (int row = 0; row < r->free; row++)
        for (int c = 0; c < dim; c++) compact[row * dim + c] = r->values[row * 48 + c];
    int status = coupled_domain_orbit(s->gid, s->nk, s->kparam, s->mode, r->free, compact, s->ncand, &s->cand[0][0], R.nkernel_ops,
                                      &R.kernel_ops[0][0], r->ndomain_ops, &r->domain_ops[0][0], r->source_order, R.ngens, &R.gens[0][0]);
    if (status) return kernel_error(MDE_DOMAIN_ORBIT, status);
    const int64_t *ints;
    const double *doubles;
    coupled_domain_orbit_ints(&ints);
    int64_t ndoubles = coupled_domain_orbit_doubles(&doubles);
    int dc = (int)ints[0], gc = (int)ints[1];
    int64_t cells = dc ? ndoubles / dc : 0;
    r->ndomains = dc;
    r->ngen = gc;
    TAKE(r->domains, dc + 1);
    TAKE(r->numbers, dc + 1);
    TAKE(r->actions, (int64_t)gc * dc + 1);
    TAKE(r->dsubs, dc + 1);
    for (int d = 0; d < dc; d++) {
        TAKE(r->domains[d], r->free * 48 + 1);
        for (int row = 0; row < r->free; row++)
            for (int c = 0; c < dim && row * dim + c < cells; c++) r->domains[d][row * 48 + c] = doubles[d * cells + row * dim + c];
        r->numbers[d] = ints[2 + d];
    }
    int64_t at = 2 + dc;
    for (int64_t i = 0; i < (int64_t)gc * dc; i++) r->actions[i] = (int32_t)ints[at + i];
    at += (int64_t)gc * dc;
    for (int d = 0; d < dc; d++) {
        memcpy(r->dsubs[d].basis, ints + at, sizeof r->dsubs[d].basis);
        int kept = (int)ints[at + 9];
        r->dsubs[d].n = kept;
        TAKE(r->dsubs[d].ops, kept + 1);
        memcpy(r->dsubs[d].ops, ints + at + 10, sizeof(rec) * (size_t)kept);
        at += 10 + 5 * kept;
    }
    r->orbit_ready = 1;
    return 0;
}

/* ------------------------------------------------------------ texts */

/* the OPD parameter name of index: a..z, then aa, bb, .., zz, aaa (as the Web lists print) */
static void put_variable(jbuf *b, int index)
{
    char c = (char)('a' + index % 26);
    for (int repeat = index / 26; repeat >= 0; repeat--) jb_raw(b, &c, 1);
}

static int rename_variables(jbuf *b, const char *text, size_t n, int offset)
{
    for (size_t i = 0; i < n; i++) {
        char c = text[i];
        if (c >= 'a' && c <= 'z') put_variable(b, offset + (c - 'a'));
        else jb_raw(b, &c, 1);
    }
    return 0;
}

/* _rename_opd_variables(_group_source_opd(opd, group), offset).strip("()") */
static int grouped_source_opd(jbuf *b, const char *opd, size_t length, int group, int offset)
{
    const char *a = opd, *e = opd + length;
    while (a < e && sel_is_space(*a)) a++;
    while (e > a && sel_is_space(e[-1])) e--;
    while (a < e && (*a == '(' || *a == ')')) a++;
    while (e > a && (e[-1] == '(' || e[-1] == ')')) e--;
    int columns = 1;
    for (const char *c = a; c < e; c++) columns += *c == ',' || *c == ';';
    if (group <= 0) group = columns;
    jbuf joined = {0};
    jb_reset(&joined);
    jb_lit(&joined, "(");
    int index = 0;
    const char *start = a;
    for (const char *c = a;; c++) {
        if (c < e && *c != ',' && *c != ';') continue;
        const char *x = start, *y = c;
        while (x < y && sel_is_space(*x)) x++;
        while (y > x && sel_is_space(y[-1])) y--;
        if (index) jb_lit(&joined, index % group == 0 ? ";" : ",");
        jb_raw(&joined, x, (size_t)(y - x));
        index++;
        if (c >= e) break;
        start = c + 1;
    }
    jb_lit(&joined, ")");
    jbuf renamed = {0};
    jb_reset(&renamed);
    int status = rename_variables(&renamed, joined.data, joined.used, offset);
    if (!status) {
        const char *x = renamed.data, *y = renamed.data + renamed.used;
        while (x < y && (*x == '(' || *x == ')')) x++;
        while (y > x && (y[-1] == '(' || y[-1] == ')')) y--;
        jb_raw(b, x, (size_t)(y - x));
    }
    jb_free(&joined);
    jb_free(&renamed);
    return status;
}

/* _symbolic_opd(...).strip("()") */
static int symbolic_opd(jbuf *b, int dim, int free, const double *values, int offset, int group)
{
    if (!group) group = dim;
    jbuf text = {0};
    jb_reset(&text);
    jb_lit(&text, "(");
    for (int col = 0; col < dim; col++) {
        if (col) jb_lit(&text, col % group == 0 ? ";" : ",");
        jbuf expr = {0};
        jb_reset(&expr);
        int any = 0;
        for (int row = 0; row < free; row++) {
            jbuf var = {0};
            jb_reset(&var);
            put_variable(&var, offset + row);
            double v = values[row * 48 + col];
            char term[96];
            if (fabs(v) < 5e-7) {
                jb_free(&var);
                continue;
            }
            if (fabs(v - 1.0) < 5e-7) snprintf(term, sizeof term, "%s", var.data);
            else if (fabs(v + 1.0) < 5e-7) snprintf(term, sizeof term, "-%s", var.data);
            else snprintf(term, sizeof term, "%.3f%s", v, var.data);
            jb_free(&var);
            if (any) jb_lit(&expr, "+");
            jb_lit(&expr, term);
            any = 1;
        }
        if (!any) {
            jb_lit(&text, "0");
        } else {
            for (size_t i = 0; i < expr.used; i++) {
                if (expr.data[i] == '+' && i + 1 < expr.used && expr.data[i + 1] == '-') continue;
                jb_raw(&text, expr.data + i, 1);
            }
        }
        jb_free(&expr);
    }
    jb_lit(&text, ")");
    const char *x = text.data, *y = text.data + text.used;
    while (x < y && (*x == '(' || *x == ')')) x++;
    while (y > x && (y[-1] == '(' || y[-1] == ')')) y--;
    jb_raw(b, x, (size_t)(y - x));
    jb_free(&text);
    return 0;
}

typedef struct {
    const double *values;
} float_rows;

static int float_nonzero(const void *ctx, int row, int col) { return fabs(((const float_rows *)ctx)->values[row * 48 + col]) > 1e-9; }

/* coupled._fixed_dimensions_for_point_actions over the distinct point matrices */
static int fixed_dimensions(const rec *ops, int n, int64_t *polar, int64_t *strain)
{
    int64_t (*seen)[9] = calloc((size_t)(n ? n : 1), sizeof *seen);
    if (!seen) return md_fail(MDE_MEMORY, "out of memory");
    int unique = 0;
    int64_t polar_c = 0, strain_c = 0;
    for (int i = 0; i < n; i++) {
        int64_t m[9];
        for (int k = 0; k < 9; k++) m[k] = section_int(&T.point_op, (ops[i][4] - 1) * 9 + k);
        int dup = 0;
        for (int j = 0; j < unique && !dup; j++) dup = !memcmp(seen[j], m, sizeof m);
        if (dup) continue;
        memcpy(seen[unique++], m, sizeof m);
        int64_t trace = m[0] + m[4] + m[8], squared = 0;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) squared += m[3 * r + c] * m[3 * c + r];
        polar_c += trace;
        strain_c += trace * trace + squared;
    }
    free(seen);
    if (!unique) return md_fail(MDE_VALUE, "empty point-action group");
    if (polar_c % unique || strain_c % (2 * unique)) return md_fail(MDE_VALUE, "point actions do not define integral fixed-space dimensions");
    *polar = polar_c / unique;
    *strain = strain_c / (2 * unique);
    return 0;
}

/* coupled._ordinary_coupled_ferroic_properties: 0 classified into b, 1 None */
static int ordinary_ferroic(jbuf *b, const rec *child_ops, int nchild, srow_t *const *rows, int *none)
{
    *none = 0;
    for (int i = 0; i < R.nslots; i++) {
        const jvalue *ferroic = member(rows[i]->iso, "ferroic");
        if (rows[i]->dynamic || !jis(ferroic, J_OBJECT)) {
            *none = 1;
            return 0;
        }
    }
    int64_t pp, ps, cp, cs;
    MD_TRY(fixed_dimensions(R.kernel_ops, R.nkernel_ops, &pp, &ps));
    MD_TRY(fixed_dimensions(child_ops, nchild, &cp, &cs));
    jb_lit(b, "[");
    if (cp > pp) {
        int proper = 0;
        for (int i = 0; i < R.nslots; i++) {
            const jvalue *f = member(rows[i]->iso, "ferroic");
            proper |= member_int(f, "raw") == 1 && member_int(f, "ferroelectric") == 1;
        }
        jb_lit(b, proper ? "\"proper ferroelectric\"" : "\"improper ferroelectric\"");
    }
    if (cs > ps) {
        int proper = 0;
        for (int i = 0; i < R.nslots; i++) {
            const jvalue *f = member(rows[i]->iso, "ferroic");
            proper |= member_int(f, "raw") == 1 && member_int(f, "ferroelastic") == 1;
        }
        if (cp > pp) jb_lit(b, ", ");
        jb_lit(b, proper ? "\"proper ferroelastic\"" : "\"improper ferroelastic\"");
    }
    jb_lit(b, "]");
    return 0;
}

/* coupled._unique_static_factor_presentation: 1 with basis/origin */
static int static_presentation(int ordinary, const int64_t *final_basis, const rec *final_ops, int nfinal, const int32_t *domains,
                               srow_t *const *rows, int64_t *basis_out, int64_t *origin_out, int *found)
{
    *found = 0;
    for (int i = 0; i < R.nslots; i++)
        if (domains[i]) return 0;
    int candidates = 0;
    for (int i = 0; i < R.nslots; i++) {
        const srow_t *r = rows[i];
        if (!r->has_stored || r->dynamic) return 0;
        if (r->child == ordinary) {
            int same = same_subgroup(r->sbasis, r->nstored, &r->stored[0][0], final_basis, nfinal, &final_ops[0][0]);
            MD_TRY(subgroup_ops_status(same));
            if (same) {
                if (!candidates || memcmp(basis_out, r->sbasis, sizeof r->sbasis) || memcmp(origin_out, r->sorigin, sizeof r->sorigin)) {
                    if (candidates) {
                        candidates = 2;
                    } else {
                        memcpy(basis_out, r->sbasis, sizeof r->sbasis);
                        memcpy(origin_out, r->sorigin, sizeof r->sorigin);
                        candidates = 1;
                    }
                }
                continue;
            }
        }
        if (r->child != R.sg) return 0;
        int same = same_subgroup(r->sbasis, r->nstored, &r->stored[0][0], IDENTITY, R.nparent_ops, &R.parent_ops[0][0]);
        MD_TRY(subgroup_ops_status(same));
        if (!same) return 0;
    }
    *found = candidates == 1;
    return 0;
}

/* ------------------------------------------------------------ dedup */

typedef struct {
    char *key;
    jbuf direction, rest;
    int64_t ordinal;
} out_row;

static struct {
    out_row *rows;
    int n, cap;
    char **keys;
    size_t key_cap;
} U;

static uint64_t hash_text(const char *s)
{
    uint64_t h = 1469598103934665603ull;
    for (; *s; s++) h = (h ^ (unsigned char)*s) * 1099511628211ull;
    return h;
}

static int seen_key(const char *key, int *seen)
{
    if ((size_t)U.n * 2 + 2 > U.key_cap) {
        size_t cap = U.key_cap ? U.key_cap * 2 : 256;
        char **grown = calloc(cap, sizeof *grown);
        if (!grown) return md_fail(MDE_MEMORY, "out of memory");
        for (size_t i = 0; i < U.key_cap; i++)
            if (U.keys[i]) {
                size_t at = hash_text(U.keys[i]) & (cap - 1);
                while (grown[at]) at = (at + 1) & (cap - 1);
                grown[at] = U.keys[i];
            }
        free(U.keys);
        U.keys = grown;
        U.key_cap = cap;
    }
    size_t at = hash_text(key) & (U.key_cap - 1);
    while (U.keys[at]) {
        if (!strcmp(U.keys[at], key)) {
            *seen = 1;
            return 0;
        }
        at = (at + 1) & (U.key_cap - 1);
    }
    *seen = 0;
    U.keys[at] = (char *)key;
    return 0;
}

static void unique_release(void)
{
    for (int i = 0; i < U.n; i++) {
        free(U.rows[i].key);
        jb_free(&U.rows[i].direction);
        jb_free(&U.rows[i].rest);
    }
    free(U.rows);
    free(U.keys);
    memset(&U, 0, sizeof U);
}

/* ------------------------------------------------------------ rows */

static int emit_tuple(const int32_t *domains, const int64_t *kept_basis, const rec *kept_ops, int nkept, const int64_t *binary,
                      const int64_t *coupled, srow_t *const *rows, int64_t ordinal, int has_setting, int64_t setting_id, int has_cell,
                      const double *cell, int parametric, int64_t parent_group, int64_t *flat, int *emitted)
{
    *emitted = 0;
    int64_t group, ordinary, basis[9], origin[4];
    if (R.magnetic_request) {
        int length = id_subgroup_magnetic_identify(R.sg, kept_basis, nkept, &kept_ops[0][0], 0, 8192, flat);
        if (length < 0) return kernel_error(MDE_ID_SUBGROUP, length);
        if (!length) return 0;
        group = flat[0];
        ordinary = flat[1];
        memcpy(basis, flat + 2, sizeof basis);
        memcpy(origin, flat + 11, sizeof origin);
    } else {
        int64_t out[14];
        int status = id_subgroup_identify(R.sg, kept_basis, nkept, &kept_ops[0][0], 0, out);
        if (status < 0) return kernel_error(MDE_ID_SUBGROUP, status);
        if (!status) return 0;
        group = ordinary = out[0];
        memcpy(basis, out + 1, sizeof basis);
        memcpy(origin, out + 10, sizeof origin);
    }
    *emitted = 1;
    int n = R.nslots;
    jbuf label = {0}, opd = {0}, active = {0};
    jb_reset(&label);
    jb_reset(&opd);
    jb_reset(&active);
    strings seen = {0};
    seen.items = calloc(64, sizeof(char *));
    int seen_cap = 64, active_parts = 0;
    jbuf *parts = calloc((size_t)n, sizeof(jbuf)), *component_opd = calloc((size_t)n, sizeof(jbuf));
    if (!seen.items || !parts || !component_opd) return md_fail(MDE_MEMORY, "out of memory");
    int status = 0, offset = 0, free_total = 0;
    jb_lit(&active, "[");
    for (int i = 0; i < n && !status; i++) {
        cslot_t *s = R.slots + i;
        srow_t *r = rows[i];
        int d = domains[i];
        jb_reset(&parts[i]);
        jb_reset(&component_opd[i]);
        const jvalue *opd_label = member(r->iso, "opd_label");
        if (jtruthy(opd_label)) put_str_member(&parts[i], opd_label);
        else put_str_member(&parts[i], member(member(r->row, "direction"), "label"));
        char number[32];
        snprintf(number, sizeof number, "(%lld)", (long long)r->numbers[d]);
        jb_lit(&parts[i], number);
        jb_raw(&label, parts[i].data, parts[i].used);
        if (d == 0) {
            const jvalue *text = member(r->iso, "display_opd");
            if (!jtruthy(text)) text = member(r->iso, "source_opd");
            jbuf t = {0};
            jb_reset(&t);
            put_str_member(&t, text);
            status = grouped_source_opd(&component_opd[i], t.data, t.used, s->little_dim, offset);
            jb_free(&t);
        } else {
            status = symbolic_opd(&component_opd[i], s->full_dim, r->free, r->domains[d], offset, s->little_dim);
        }
        if (status) break;
        if (i) jb_lit(&opd, "|");
        jb_raw(&opd, component_opd[i].data, component_opd[i].used);
        strings vectors = {0};
        int widths[48];
        for (int row = 0; row < r->free && row < 48; row++) widths[row] = s->full_dim;
        float_rows ctx = {r->domains[d]};
        status = sel_active_k_vectors_by(R.sg, s->gid, s->full_dim, r->free, widths, float_nonzero, &ctx, &s->active, &vectors);
        if (status) break;
        int fresh = 0;
        jbuf joined = {0};
        jb_reset(&joined);
        for (int v = 0; v < vectors.count; v++) {
            int known = 0;
            for (int k = 0; k < seen.count && !known; k++) known = !strcmp(seen.items[k], vectors.items[v]);
            if (known) continue;
            if (fresh) jb_lit(&joined, ",");
            jb_lit(&joined, vectors.items[v]);
            fresh++;
            if (seen.count == seen_cap) {
                seen_cap *= 2;
                char **grown = realloc(seen.items, sizeof(char *) * (size_t)seen_cap);
                if (!grown) {
                    status = md_fail(MDE_MEMORY, "out of memory");
                    break;
                }
                seen.items = grown;
            }
            seen.items[seen.count++] = vectors.items[v];
            vectors.items[v] = 0;
        }
        if (fresh) {
            if (active_parts) jb_lit(&active, ", ");
            jb_string(&active, joined.data, joined.used);
            active_parts++;
        }
        jb_free(&joined);
        strings_free(&vectors);
        offset += r->free;
        free_total += r->free;
    }
    jb_lit(&active, "]");
    int found_presentation = 0;
    int64_t pbasis[9], porigin[4];
    if (!status && !R.magnetic_request)
        status = static_presentation((int)ordinary, kept_basis, kept_ops, nkept, domains, rows, pbasis, porigin, &found_presentation);
    if (!found_presentation) {
        memcpy(pbasis, basis, sizeof pbasis);
        memcpy(porigin, origin, sizeof porigin);
    }
    int64_t det = abs64(det3(basis));
    sel_setting displayed;
    rat brows[3][3];
    if (!status)
        status = sel_display(R.sg, (int)ordinary, pbasis, porigin, has_setting, setting_id, has_cell, cell, parametric, 1,
                             R.magnetic_request ? (int)group : 0, &displayed, brows);
    jbuf ferroic = {0};
    jb_reset(&ferroic);
    int ferroic_none = 0;
    if (!status) {
        if (R.magnetic_request) status = sel_mag_ferroic(&ferroic, (int)group);
        else status = ordinary_ferroic(&ferroic, kept_ops, nkept, rows, &ferroic_none);
    }
    if (status) {
        for (int i = 0; i < n; i++) {
            jb_free(&parts[i]);
            jb_free(&component_opd[i]);
        }
        free(parts);
        free(component_opd);
        strings_free(&seen);
        jb_free(&label);
        jb_free(&opd);
        jb_free(&active);
        jb_free(&ferroic);
        return status;
    }
    out_row row = {0};
    jb_reset(&row.direction);
    jb_reset(&row.rest);
    row.ordinal = ordinal;
    jbuf *b = &row.direction;
    jb_lit(b, "{\"direction\": {\"label\": ");
    jb_string(b, label.data, label.used);
    jb_lit(b, ", \"opd\": ");
    jbuf opd_text = {0};
    jb_reset(&opd_text);
    jb_lit(&opd_text, "(");
    jb_raw(&opd_text, opd.data, opd.used);
    jb_lit(&opd_text, ")");
    jb_string(b, opd_text.data, opd_text.used);
    jb_lit(b, ", \"components\": [");
    for (int i = 0; i < n; i++) {
        srow_t *r = rows[i];
        int d = domains[i];
        if (i) jb_lit(b, ", ");
        jb_lit(b, "{\"slot\": ");
        jb_int(b, i + 1);
        jb_lit(b, ", \"label\": ");
        jb_string(b, parts[i].data, parts[i].used);
        jb_lit(b, ", \"opd\": ");
        jbuf t = {0};
        jb_reset(&t);
        put_str_member(&t, member(r->iso, "opd_label"));
        jb_string(b, t.data, t.used);
        jb_free(&t);
        jb_lit(b, ", \"domain\": ");
        jb_int(b, r->numbers[d]);
        jb_lit(b, ", \"display_opd\": ");
        jb_reset(&t);
        jb_lit(&t, "(");
        jb_raw(&t, component_opd[i].data, component_opd[i].used);
        jb_lit(&t, ")");
        jb_string(b, t.data, t.used);
        jb_free(&t);
        jb_lit(b, "}");
    }
    jb_lit(b, "]}, \"isotropy\": {\"row_id\": ");
    b = &row.rest;
    jb_lit(b, ", \"opd_label\": ");
    jb_string(b, label.data, label.used);
    jb_lit(b, ", \"direction_label\": ");
    jb_string(b, label.data, label.used);
    jb_lit(b, ", \"direction_domains\": [");
    for (int i = 0; i < n; i++) {
        if (i) jb_lit(b, ", ");
        jb_int(b, rows[i]->numbers[domains[i]]);
    }
    jb_lit(b, "], \"_slot_source_numeric_rows\": [");
    for (int i = 0; i < n; i++) {
        srow_t *r = rows[i];
        if (i) jb_lit(b, ", ");
        jb_lit(b, "[");
        for (int row_index = 0; row_index < r->free; row_index++) {
            if (row_index) jb_lit(b, ", ");
            jb_lit(b, "[");
            for (int c = 0; c < R.slots[i].full_dim; c++) {
                if (c) jb_lit(b, ", ");
                jb_float(b, r->domains[domains[i]][row_index * 48 + c]);
            }
            jb_lit(b, "]");
        }
        jb_lit(b, "]");
    }
    jb_lit(b, "], \"display_opd\": ");
    jb_string(b, opd_text.data, opd_text.used);
    jb_lit(b, ", \"source_opd\": ");
    jb_string(b, opd_text.data, opd_text.used);
    jb_lit(b, ", \"op_rows\": ");
    jb_int(b, free_total);
    jb_lit(b, ", \"free\": ");
    jb_int(b, free_total);
    jb_lit(b, ", \"s\": ");
    jb_int(b, det);
    jb_lit(b, ", \"i\": ");
    int has_index;
    int64_t index = R.magnetic_request ? sel_mag_index((int)parent_group, (int)group, det, &has_index)
                                       : sel_subgroup_index(R.sg, (int)group, det, &has_index);
    if (has_index) jb_int(b, index);
    else jb_lit(b, "null");
    jb_lit(b, ", \"k_active_vectors\": ");
    jb_raw(b, active.data, active.used);
    jb_lit(b, ", \"k_active\": ");
    jbuf joined = {0};
    jb_reset(&joined);
    {
        jdoc *doc = jdoc_parse(active.data, active.used);
        const jvalue *list = doc ? jdoc_root(doc) : 0;
        for (int64_t i = 0; list && i < list->count; i++) {
            if (i) jb_lit(&joined, ";");
            jb_raw(&joined, jat(list, i)->u.s, (size_t)jat(list, i)->count);
        }
        jdoc_free(doc);
    }
    jb_string(b, joined.data, joined.used);
    jb_free(&joined);
    jbuf subgroup_json = {0}, extra = {0};
    jb_reset(&subgroup_json);
    jb_reset(&extra);
    jb_lit(&subgroup_json, "{\"number\": ");
    jb_int(&subgroup_json, group);
    jb_lit(&subgroup_json, ", \"symbol\": ");
    if (R.magnetic_request) {
        jbuf full = {0};
        jb_reset(&full);
        sel_mag_subgroup_json(&full, (int)group);
        jdoc *doc = jdoc_parse(full.data, full.used);
        const jvalue *m = doc ? jdoc_root(doc) : 0;
        jb_value(&subgroup_json, jget(m, "symbol"));
        jb_lit(&subgroup_json, ", \"ordinary_number\": ");
        jb_int(&subgroup_json, ordinary);
        jb_lit(&subgroup_json, ", \"display_label\": ");
        jb_value(&subgroup_json, jget(m, "display_label"));
        jdoc_free(doc);
        jb_free(&full);
    } else {
        char *symbol = sel_default_symbol((int)group);
        jb_cstring(&subgroup_json, symbol ? symbol : "");
        free(symbol);
    }
    jb_lit(&subgroup_json, "}");
    if (R.magnetic_request) {
        jb_lit(&extra, ", \"magnetic_parent_group\": ");
        jb_int(&extra, parent_group);
        jb_lit(&extra, ", \"magnetic_subgroup_selection\": ");
        sel_put_selection(&extra, flat);
    }
    jb_lit(&extra, ", \"source_operation_records\": ");
    sel_put_records(&extra, &kept_ops[0][0], nkept);
    jb_lit(&extra, ", \"ferroic_properties\": ");
    if (R.magnetic_request || !ferroic_none) jb_raw(&extra, ferroic.data, ferroic.used);
    else jb_lit(&extra, "[]");
    jb_lit(&extra, ", \"subduction\": {\"expected\": [");
    for (int i = 0; i < n; i++) {
        if (i) jb_lit(&extra, ", ");
        jb_int(&extra, rows[i]->free);
    }
    jb_lit(&extra, "], \"binary_slot_kernel\": [");
    for (int i = 0; i < n; i++) {
        if (i) jb_lit(&extra, ", ");
        jb_int(&extra, binary[i]);
    }
    jb_lit(&extra, "], \"coupled_common_kernel\": [");
    for (int i = 0; i < n; i++) {
        if (i) jb_lit(&extra, ", ");
        jb_int(&extra, coupled[i]);
    }
    jb_lit(&extra, "]}");
    sel_put_subgroup(b, R.sg, (int)ordinary, basis, origin, &displayed, (const rat (*)[3])brows, det, &subgroup_json, &extra);
    jb_lit(b, "}}");
    jb_free(&subgroup_json);
    jb_free(&extra);
    jbuf key = {0};
    jb_reset(&key);
    for (int i = 0; i < n; i++) {
        put_str_member(&key, member(rows[i]->iso, "opd_label"));
        jb_raw(&key, "\x1f", 1);
    }
    int64_t den = origin[3] ? origin[3] : 1;
    char words[256];
    snprintf(words, sizeof words, "%lld|%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld|%lld,%lld,%lld|%lld", (long long)group,
             (long long)basis[0], (long long)basis[1], (long long)basis[2], (long long)basis[3], (long long)basis[4], (long long)basis[5],
             (long long)basis[6], (long long)basis[7], (long long)basis[8], (long long)(((origin[0] % den) + den) % den),
             (long long)(((origin[1] % den) + den) % den), (long long)(((origin[2] % den) + den) % den), (long long)den);
    jb_lit(&key, words);
    row.key = key.data;
    for (int i = 0; i < n; i++) {
        jb_free(&parts[i]);
        jb_free(&component_opd[i]);
    }
    free(parts);
    free(component_opd);
    strings_free(&seen);
    jb_free(&label);
    jb_free(&opd);
    jb_free(&opd_text);
    jb_free(&active);
    jb_free(&ferroic);
    if (status) {
        free(row.key);
        jb_free(&row.direction);
        jb_free(&row.rest);
        return status;
    }
    int duplicate;
    MD_TRY(seen_key(row.key, &duplicate));
    if (duplicate) {
        free(row.key);
        jb_free(&row.direction);
        jb_free(&row.rest);
        *emitted = 2;
        return 0;
    }
    if (U.n == U.cap) {
        int cap = U.cap ? 2 * U.cap : 64;
        out_row *grown = realloc(U.rows, sizeof *grown * (size_t)cap);
        if (!grown) return md_fail(MDE_MEMORY, "out of memory");
        U.rows = grown;
        U.cap = cap;
    }
    U.rows[U.n++] = row;
    return 0;
}

static int load_tables(void)
{
    if (T.ready) return 0;
    if (sel_tables() || sel_mag_tables() || geo_load() || section_find("little.little_irr_space_group", &T.irr_sg) ||
        section_find("space.ipoint_op", &T.point_op) || section_find("magnetic.mag_point_op_mag2nonmag", &T.mag2nonmag))
        return md_fail(MDE_VALUE, "coupled Source tables are inconsistent");
    T.ready = 1;
    return 0;
}

static int row_id_of(const srow_t *r) { return (int)member_int(r->iso, "row_id"); }

static int compare_rows(const void *a, const void *b)
{
    int x = row_id_of(a), y = row_id_of(b);
    return (x > y) - (x < y);
}

static int coupled(const jvalue *request, jbuf *out)
{
    MD_TRY(load_tables());
    const jvalue *slots = jget(request, "slots");
    int n = jis(slots, J_ARRAY) ? (int)slots->count : 0;
    if (n < 2) return md_fail(MDE_VALUE, "coupled_opd_rows requires at least two slots");
    const jvalue *limit_value = jget(request, "row_limit");
    int has_limit = jis(limit_value, J_INT);
    int64_t limit = has_limit ? limit_value->u.i : 0;
    if (has_limit && limit < 1) return md_fail(MDE_VALUE, "row_limit must be positive");
    int64_t sg;
    MD_TRY(int_of(jget(request, "parent_sg"), &sg));
    int64_t setting_id;
    MD_TRY(int_of(jget(request, "parent_setting_id"), &setting_id));
    const jvalue *cell_value = jget(request, "parent_cell");
    double cell[6];
    int has_cell = jis(cell_value, J_ARRAY) && cell_value->count == 6;
    for (int i = 0; has_cell && i < 6; i++)
        if (py_float_of(jat(cell_value, i), cell + i)) return md_fail(MDE_TYPE, "parent cell must be numeric");
    memset(&R, 0, sizeof R);
    R.sg = (int)sg;
    R.nslots = n;
    TAKE(R.slots, n);
    int parametric = 0;
    int64_t groups[16];
    int ngroups = 0;
    for (int i = 0; i < n; i++) {
        cslot_t *s = R.slots + i;
        const jvalue *slot = jat(slots, i), *irrep = member(slot, "irrep");
        s->slot = slot;
        s->gid = (int)member_int(irrep, "gid");
        s->full_dim = (int)member_int(irrep, "full_dim");
        s->little_dim = (int)member_int(irrep, "little_dim");
        if (!s->little_dim) s->little_dim = s->full_dim;
        s->magnetic = jtruthy(member(irrep, "magnetic"));
        s->gid_sg = (int)section_int(&T.irr_sg, s->gid - 1);
        const jvalue *kparam = member(slot, "source_kparam");
        s->nk = jis(kparam, J_ARRAY) ? (int)kparam->count : 0;
        if (s->nk > 6) return md_fail(MDE_VALUE, "source kparam is too long");
        for (int k = 0; k < s->nk; k++) MD_TRY(int_of(jat(kparam, k), s->kparam + k));
        const jvalue *kp = member(slot, "k_params"), *display = member(slot, "display_k_params");
        s->has_kparams = jtruthy(kp);
        parametric |= s->has_kparams;
        sel_params_of(jtruthy(display) ? display : kp, &s->active);
        if (s->magnetic) {
            int64_t g = mag_orderparam_group(s->gid_sg);
            int known = 0;
            for (int k = 0; k < ngroups; k++) known |= groups[k] == g;
            if (!known && ngroups < 16) groups[ngroups++] = g;
        }
    }
    R.mixed_dynamic = parametric;
    R.magnetic_request = ngroups > 0;
    if (ngroups > 1) {
        for (int a = 0; a < ngroups; a++)
            for (int b = a + 1; b < ngroups; b++)
                if (groups[b] < groups[a]) {
                    int64_t t = groups[a];
                    groups[a] = groups[b];
                    groups[b] = t;
                }
        char text[256] = "";
        size_t at = 0;
        for (int k = 0; k < ngroups; k++) at += (size_t)snprintf(text + at, sizeof text - at, k ? ", %lld" : "%lld", (long long)groups[k]);
        return md_fail(MDE_VALUE, "magnetic coupled slots use different parent groups: [%s]", text);
    }
    R.group = R.magnetic_request ? (int)groups[0] : 0;
    for (int i = 0; i < n; i++) R.slots[i].mode = R.magnetic_request ? (R.slots[i].magnetic ? 1 : 2) : 0;
    for (int i = 0; i < n; i++) {
        cslot_t *s = R.slots + i;
        TAKE(s->cand, 193);
        TAKE(s->cops, 193);
        TAKE(s->csigns, 193);
        if (R.magnetic_request) {
            int count = magnetic_records_iso(R.group, 192, &s->cand[0][0]);
            if (count < 0) return kernel_error(MDE_MAGNETIC, -count);
            s->ncand = count;
        } else {
            s->ncand = geo_sg_records(s->gid_sg, s->cand);
        }
        for (int k = 0; k < s->ncand; k++) s->cops[k] = mag_point_action(s->mode, (int)s->cand[k][4], &s->csigns[k]);
    }
    coupled_tuples_reset();

    double **kvalues;
    int *kfree, *knvalues;
    cslot_t **all;
    TAKE(kvalues, n);
    TAKE(kfree, n);
    TAKE(knvalues, n);
    TAKE(all, n);
    for (int i = 0; i < n; i++) {
        int dim = R.slots[i].full_dim;
        all[i] = R.slots + i;
        kfree[i] = dim;
        knvalues[i] = dim * 48;
        TAKE(kvalues[i], dim * 48 + 1);
        for (int r = 0; r < dim; r++) kvalues[i][r * 48 + r] = 1.0;
    }
    for (int i = 0; i < n; i++) MD_TRY(combined_basis(all + i, 1, kfree + i, kvalues + i, knvalues + i, R.slots[i].kernel_basis));
    int64_t kernel_basis[9];
    MD_TRY(combined_basis(all, n, kfree, kvalues, knvalues, kernel_basis));
    MD_TRY(kernel_records(kernel_basis, R.group, &R.kernel_ops, &R.nkernel_ops));
    int64_t fractions = abs64(det3(kernel_basis));
    TAKE(R.gens, R.nkernel_ops + 1);
    if (fractions <= 0 || R.nkernel_ops % fractions) {
        memcpy(R.gens, R.kernel_ops, sizeof(rec) * (size_t)R.nkernel_ops);
        R.ngens = R.nkernel_ops;
    } else {
        int count = 0;
        for (int pass = 0; pass < 2; pass++)
            for (int64_t i = 0; pass ? i < fractions : i < R.nkernel_ops; i += pass ? 1 : fractions) {
                int dup = 0;
                for (int k = 0; k < count && !dup; k++) dup = !memcmp(R.gens[k], R.kernel_ops[i], sizeof(rec));
                if (!dup) memcpy(R.gens[count++], R.kernel_ops[i], sizeof(rec));
            }
        R.ngens = count;
    }
    MD_TRY(kernel_records(IDENTITY, 0, &R.parent_ops, &R.nparent_ops));

    int64_t combos = 1;
    for (int i = 0; i < n; i++) {
        cslot_t *s = R.slots + i;
        const jvalue *rows = member(s->slot, "opd_rows");
        s->nrows = jis(rows, J_ARRAY) ? (int)rows->count : 0;
        TAKE(s->rows, s->nrows + 1);
        int all_static = s->nrows > 0;
        for (int k = 0; k < s->nrows; k++) {
            srow_t *r = s->rows + k;
            r->row = jat(rows, k);
            r->iso = member(r->row, "isotropy");
            r->dynamic = jtruthy(member(r->iso, "dynamic"));
            all_static &= !r->dynamic;
        }
        if (all_static) {
            for (int a = 1; a < s->nrows; a++) {
                srow_t t = s->rows[a];
                int b = a - 1;
                while (b >= 0 && compare_rows(&s->rows[b], &t) > 0) {
                    s->rows[b + 1] = s->rows[b];
                    b--;
                }
                s->rows[b + 1] = t;
            }
        }
        for (int k = 0; k < s->nrows; k++) MD_TRY(row_orderparam(s, s->rows + k));
        combos *= s->nrows;
    }
    int *mag2nonmag;
    int mag_count = R.magnetic_request ? (int)section_count(&T.mag2nonmag) : 0;
    TAKE(mag2nonmag, mag_count + 1);
    for (int i = 0; i < mag_count; i++) mag2nonmag[i] = (int)section_int(&T.mag2nonmag, i);
    int64_t *flat;
    TAKE(flat, 8192);
    int *index;
    srow_t **rows;
    TAKE(index, n);
    TAKE(rows, n);
    int64_t ordinal = 0;
    int done = 0;
    for (int64_t combo = 0; combo < combos && !done; combo++) {
        int64_t rest = combo;
        for (int i = 0; i < n; i++) {
            index[i] = (int)(rest % R.slots[i].nrows);
            rest /= R.slots[i].nrows;
            rows[i] = R.slots[i].rows + index[i];
        }
        for (int i = 0; i < n; i++) MD_TRY(slot_block(R.slots + i, rows[i]));
        for (int i = 0; i < n; i++) MD_TRY(slot_orbit(R.slots + i, rows[i]));
        int *sizes;
        TAKE(sizes, n);
        int64_t capacity = 1;
        for (int i = 0; i < n; i++) {
            if (!rows[i]->ngen) return md_fail(MDE_INDEX, "tuple index out of range");
            sizes[i] = rows[i]->ndomains;
            capacity *= sizes[i];
        }
        int width = 0;
        for (int i = 0; i < n; i++) width += sizes[i];
        int ngen = rows[0]->ngen;
        int32_t *joint, *tuples;
        TAKE(joint, (int64_t)ngen * width + 1);
        for (int g = 0; g < ngen; g++)
            for (int i = 0, at = 0; i < n; i++) {
                for (int d = 0; d < sizes[i]; d++) joint[(int64_t)g * width + at + d] = g < rows[i]->ngen ? rows[i]->actions[(int64_t)g * sizes[i] + d] : 0;
                at += sizes[i];
            }
        TAKE(tuples, capacity * n + n);
        int count = relative_domain_tuples(n, sizes, ngen, joint, capacity, tuples);
        if (count < 0) return md_fail(MDE_VALUE, "coupled domain tuples failed (%d)", -count);
        int64_t *words;
        TAKE(words, 22 * n);
        int64_t total_values = 0, total_domains = 0, total_ops = 0, total_orders = 0;
        for (int i = 0; i < n; i++) {
            total_values += (int64_t)rows[i]->ndomains * rows[i]->free * 48;
            total_domains += rows[i]->ndomains;
            for (int d = 0; d < rows[i]->ndomains; d++) total_ops += rows[i]->dsubs[d].n;
            if (rows[i]->has_stored) total_orders += rows[i]->nstored;
        }
        double *values;
        int64_t *bases, *orders;
        int *op_counts;
        rec *ops;
        TAKE(values, total_values + 1);
        TAKE(bases, 9 * total_domains + 9);
        TAKE(op_counts, total_domains + 1);
        TAKE(ops, total_ops + 1);
        TAKE(orders, total_orders + 1);
        int64_t vat = 0, dat = 0, oat = 0, rat_at = 0;
        for (int i = 0; i < n; i++) {
            cslot_t *s = R.slots + i;
            srow_t *r = rows[i];
            int64_t *w = words + 22 * i;
            w[0] = s->gid;
            w[1] = s->nk;
            for (int k = 0; k < 6; k++) w[2 + k] = k < s->nk ? s->kparam[k] : 0;
            w[8] = s->full_dim;
            w[9] = r->free;
            w[10] = s->mode;
            memcpy(w + 11, s->kernel_basis, sizeof s->kernel_basis);
            w[20] = r->ndomains;
            w[21] = r->has_stored ? r->nstored : -1;
            for (int d = 0; d < r->ndomains; d++) {
                memcpy(values + vat, r->domains[d], sizeof(double) * (size_t)(r->free * 48));
                vat += r->free * 48;
                memcpy(bases + 9 * dat, r->dsubs[d].basis, sizeof r->dsubs[d].basis);
                op_counts[dat++] = r->dsubs[d].n;
                memcpy(ops + oat, r->dsubs[d].ops, sizeof(rec) * (size_t)r->dsubs[d].n);
                oat += r->dsubs[d].n;
            }
            if (r->has_stored)
                for (int k = 0; k < r->nstored; k++) orders[rat_at++] = r->stored[k][4];
        }
        int status = coupled_tuples(R.sg, R.group, n, words, kernel_basis, values, bases, op_counts, &ops[0][0], orders, R.magnetic_request,
                                    mag_count, mag2nonmag, count, tuples);
        if (status == 5) return md_fail(MDE_MEMORY, "coupled tuple scan ran out of memory");
        if (status) return md_fail(MDE_ARITH, "coupled tuple scan value out of range");
        const int64_t *kept_words;
        int64_t used = coupled_tuples_ints(&kept_words);
        int64_t kept_count = used ? kept_words[0] : 0;
        int64_t *kept_copy;
        TAKE(kept_copy, used + 1);
        memcpy(kept_copy, kept_words, sizeof(int64_t) * (size_t)used);
        int64_t failure[3];
        coupled_tuples_failure(failure);
        int64_t at = 1;
        for (int64_t k = 0; k < kept_count && !done; k++) {
            int64_t tuple_index = kept_copy[at], op_count = kept_copy[at + 1];
            const int64_t *basis = kept_copy + at + 2;
            const rec *records = (const rec *)(kept_copy + at + 11);
            at += 11 + 5 * op_count;
            const int64_t *binary = kept_copy + at, *coupled_kernel = kept_copy + at + n;
            at += 2 * n;
            int emitted;
            int64_t next = ordinal + 1;
            MD_TRY(emit_tuple(tuples + tuple_index * n, basis, records, (int)op_count, binary, coupled_kernel, rows, next, 1, setting_id,
                              has_cell, cell, parametric, R.group, flat, &emitted));
            if (emitted) ordinal = next;
            if (has_limit && U.n >= limit) done = 1;
        }
        if (done) break;
        if (failure[0] >= 0) {
            int kind = (int)failure[1], code = (int)failure[2];
            if (kind == 1) return kernel_error(MDE_ORDERPARAM, code);
            if (kind == 2) return md_fail(MDE_VALUE, "singular basis");
            if (kind == 3) return kernel_error(MDE_SUBGROUP_OPS, code);
            if (kind == 4) return code >= 100 ? kernel_error(MDE_ORDERPARAM, code - 100) : kernel_error(MDE_BRIDGE, code);
            return md_fail(MDE_OVERFLOW, "Python int too large to convert to C long");
        }
    }
    jb_lit(out, "[");
    for (int i = 0; i < U.n; i++) {
        if (i) jb_lit(out, ", ");
        jb_raw(out, U.rows[i].direction.data, U.rows[i].direction.used);
        jb_int(out, i + 1);
        jb_lit(out, ", \"source_ordinal\": ");
        jb_int(out, U.rows[i].ordinal);
        jb_lit(out, i ? ", \"coupled\": true, \"canonical\": false" : ", \"coupled\": true, \"canonical\": true");
        jb_raw(out, U.rows[i].rest.data, U.rows[i].rest.used);
    }
    jb_lit(out, "]");
    return 0;
}

int sel_coupled_rows(const jvalue *request, jbuf *out)
{
    int status = coupled(request, out);
    unique_release();
    pool_release();
    return status;
}
