/* The row layer of the spec stage: dynamic_subduction rows over a frame,
 * occurrence alias rows, kvec_standard validity, and the numpy linear
 * algebra the spec decisions read. */

#include "linalg.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "spec_core.h"
#include "tables.h"

const char *irrep_full_label(int gid);
const char *irrep_old_label(int old_id);

int k_family_records(int mode, int sg, const int64_t *basis, int preferred_count, const char *const *preferred_labels,
                     const int64_t *preferred_label_lengths, const int64_t *preferred_kparams, int family_capacity,
                     int occurrence_capacity, int64_t *families_out, int64_t *occurrences_out,
                     int64_t *occurrence_count_out);
int k_family_source_kparam(int sg, int gid, const int64_t *vector, int64_t *out);
int k_family_standard_selection(int sg, int gid, const int64_t *kparam, int64_t *out);
int subduction_rows(int sg, int family_count, const int64_t *families, int preferred_count,
                    const char *const *preferred_labels, const int64_t *preferred_lengths, const int64_t *preferred_kparams,
                    int op_count, const int *ops, int magnetic, int use_kparam);
int64_t subduction_rows_ints(const int64_t **rows, const int64_t **occurrences, int64_t *occurrence_count);
const double *subduction_rows_directions(void);
void subduction_rows_failure(int64_t *gid, int64_t *status);
int subduction_directions(int gid, int mode_case, int point_mode, int param_count, const int64_t *params,
                          int record_count, const int *records, double *directions_out);

sp_error sp_err;

int sp_fail(int kind, const char *message)
{
    if (!sp_err.kind) {
        sp_err.kind = kind;
        snprintf(sp_err.message, sizeof sp_err.message, "%s", message);
    }
    return 1;
}

typedef struct chunk {
    struct chunk *next;
    size_t used, cap;
    _Alignas(16) unsigned char data[];
} chunk;

static chunk *arena;

void *sp_alloc(size_t size)
{
    size = (size + 15) & ~(size_t)15;
    if (!arena || arena->used + size > arena->cap) {
        size_t cap = size > (1u << 20) ? size : (1u << 20);
        chunk *next = malloc(sizeof(chunk) + cap);
        if (!next) {
            sp_fail(PY_MEMORY, "spec stage ran out of memory");
            return 0;
        }
        next->next = arena;
        next->used = 0;
        next->cap = cap;
        arena = next;
    }
    void *out = arena->data + arena->used;
    arena->used += size;
    return out;
}

sp_occ *sp_occ_pool;
int sp_occ_count;
static int occ_cap;

void sp_reset(void)
{
    while (arena) {
        chunk *next = arena->next;
        free(arena);
        arena = next;
    }
    sp_err.kind = PY_NONE;
    sp_err.message[0] = 0;
    sp_occ_pool = 0;
    sp_occ_count = occ_cap = 0;
}

int sp_str_eq(sp_str a, sp_str b) { return a.length == b.length && !memcmp(a.text, b.text, (size_t)a.length); }
static struct {
    section_view full_dim, old, type, k_label, irr_k, irr_sg, star, star_count, star_pointer, point_op, lattice_ml,
        lattice_ml_denom, nonmag2mag;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (geo_load() || section_find("little.little_irr_full_dim", &T.full_dim) || section_find("little.little_irr_old", &T.old) ||
        section_find("little.little_irr_type", &T.type) || section_find("little.little_k_label", &T.k_label) ||
        section_find("little.little_irr_k", &T.irr_k) || section_find("little.little_irr_space_group", &T.irr_sg) ||
        section_find("little.little_k_star", &T.star) || section_find("little.little_k_star_count", &T.star_count) ||
        section_find("little.little_k_star_ml_pointer", &T.star_pointer) || section_find("space.ipoint_op", &T.point_op) ||
        section_find("space.lattice_ml", &T.lattice_ml) || section_find("space.lattice_ml_denom", &T.lattice_ml_denom) ||
        section_find("magnetic.mag_point_op_nonmag2mag", &T.nonmag2mag))
        return sp_fail(PY_VALUE, "spec stage Source tables are missing");
    T.ready = 1;
    return 0;
}

static sp_str table_string(const section_view *view, int64_t index)
{
    sp_str out = {"", 0};
    if (section_string(view, index, &out.text, &out.length)) return (sp_str){"", 0};
    while (out.length && (out.text[0] == ' ' || out.text[0] == '\t')) {
        out.text++;
        out.length--;
    }
    return out;
}

int sp_little_old_id(int gid) { return load() ? 0 : (int)section_int(&T.old, gid - 1); }
int sp_little_type(int gid) { return load() ? 0 : (int)section_int(&T.type, gid - 1); }
int sp_little_full_dim(int gid) { return load() ? 0 : (int)section_int(&T.full_dim, gid - 1); }
sp_str sp_irrep_full_label(int gid) { return (sp_str){irrep_full_label(gid), (int64_t)strlen(irrep_full_label(gid))}; }

static sp_str gid_k_label(int gid)
{
    int sg = (int)section_int(&T.irr_sg, gid - 1);
    int lattice = geo_sg_lattice(sg);
    return table_string(&T.k_label, (int64_t)(lattice - 1) * 27 + section_int(&T.irr_k, gid - 1) - 1);
}

int sp_kp_exact_eq(const int64_t *a, const int64_t *b) { return !memcmp(a, b, 4 * sizeof(int64_t)); }

int sp_kp_identity_eq(const int64_t *a, const int64_t *b)
{
    for (int axis = 0; axis < 3; axis++)
        if (!rat_eq(rat_make(a[axis], a[3]), rat_make(b[axis], b[3]))) return 0;
    return 1;
}

int sp_kp_periodic_eq(const int64_t *a, const int64_t *b)
{
    for (int axis = 0; axis < 3; axis++)
        if (!rat_eq(rat_mod1(rat_make(a[axis], a[3])), rat_mod1(rat_make(b[axis], b[3])))) return 0;
    return 1;
}

int sp_pml_eq(const rat *a, const rat *b)
{
    for (int axis = 0; axis < 3; axis++)
        if (!rat_eq(a[axis], b[axis])) return 0;
    return 1;
}

int sp_pml_mod1_eq(const rat *a, const rat *b)
{
    for (int axis = 0; axis < 3; axis++)
        if (!rat_eq(rat_mod1(a[axis]), rat_mod1(b[axis]))) return 0;
    return 1;
}

int sp_rows_push(sp_rows *rows, const sp_row *row)
{
    if (rows->count == rows->cap) {
        int cap = rows->cap ? rows->cap * 2 : 16;
        sp_row *grown = sp_alloc(sizeof(sp_row) * (size_t)cap);
        if (!grown) return 1;
        if (rows->count) memcpy(grown, rows->items, sizeof(sp_row) * (size_t)rows->count);
        rows->items = grown;
        rows->cap = cap;
    }
    rows->items[rows->count++] = *row;
    return 0;
}

static int occ_push(const sp_occ *occ)
{
    if (sp_occ_count == occ_cap) {
        int cap = occ_cap ? occ_cap * 2 : 256;
        sp_occ *grown = sp_alloc(sizeof(sp_occ) * (size_t)cap);
        if (!grown) return -1;
        if (sp_occ_count) memcpy(grown, sp_occ_pool, sizeof(sp_occ) * (size_t)sp_occ_count);
        sp_occ_pool = grown;
        occ_cap = cap;
    }
    sp_occ_pool[sp_occ_count] = *occ;
    return sp_occ_count++;
}

static const struct {
    int kind;
    const char *message;
} K_FAMILY_ERRORS[] = {
    {PY_RUNTIME, "k family search failed"},
    {PY_ARITH, "k family arithmetic out of range"},
    {PY_VALUE, "k family Source tables are inconsistent"},
    {PY_MEMORY, "k family search ran out of memory"},
    {PY_VALUE, "k family output capacity exceeded"},
    {PY_VALUE, "lattice holohedry has no identity operation"},
    {PY_ZERO, "zero denominator in k family Source tables"},
    {PY_INDEX, "special k family has no star vector"},
};

static int k_family_fail(int status)
{
    int code = -status;
    if (code < 1 || code > 7) {
        char message[160];
        snprintf(message, sizeof message, "k family search failed (%d)", status);
        return sp_fail(PY_RUNTIME, message);
    }
    return sp_fail(K_FAMILY_ERRORS[code].kind, K_FAMILY_ERRORS[code].message);
}

static int identity_op(const sp_frame *frame)
{
    int64_t records[GEO_MAX_OPS][5];
    if (geo_sg_records(frame->sg, records) <= 0) return 0;
    int op = (int)records[0][4];
    return frame->magnetic ? (int)section_int(&T.nonmag2mag, op - 1) : op;
}

static int *constrained_ops(const sp_frame *frame, int *count)
{
    int op = identity_op(frame);
    int *records = sp_alloc(sizeof(int) * 5 * (size_t)(frame->nops + 3));
    if (!records) return 0;
    for (int i = 0; i < frame->nops; i++)
        for (int k = 0; k < 5; k++) records[5 * i + k] = (int)frame->ops[i][k];
    for (int row = 0; row < 3; row++) {
        int *r = records + 5 * (frame->nops + row);
        r[0] = (int)frame->basis[3 * row];
        r[1] = (int)frame->basis[3 * row + 1];
        r[2] = (int)frame->basis[3 * row + 2];
        r[3] = 1;
        r[4] = op;
    }
    *count = frame->nops + 3;
    return records;
}

static int preferred_arrays(const sp_preferred *preferred, int n, const char ***labels, int64_t **lengths, int64_t **kparams)
{
    *labels = sp_alloc(sizeof(char *) * (size_t)(n + 1));
    *lengths = sp_alloc(sizeof(int64_t) * (size_t)(n + 1));
    *kparams = sp_alloc(sizeof(int64_t) * 4 * (size_t)(n + 1));
    if (!*labels || !*lengths || !*kparams) return 1;
    for (int i = 0; i < n; i++) {
        (*labels)[i] = preferred[i].label.text;
        (*lengths)[i] = preferred[i].label.length;
        memcpy(*kparams + 4 * i, preferred[i].kp, 4 * sizeof(int64_t));
        if (!preferred[i].kp[3]) return sp_fail(PY_ZERO, "preferred kparam has a zero denominator");
    }
    if (!n) {
        (*lengths)[0] = 0;
        int64_t unit[4] = {0, 0, 0, 1};
        memcpy(*kparams, unit, sizeof unit);
    }
    return 0;
}

int sp_frame_rows(const sp_frame *frame, const sp_preferred *preferred, int npreferred, sp_rows *out)
{
    SP_TRY(load());
    const int64_t *b = frame->basis;
    int64_t det = b[0] * (b[4] * b[8] - b[5] * b[7]) - b[1] * (b[3] * b[8] - b[5] * b[6]) + b[2] * (b[3] * b[7] - b[4] * b[6]);
    if (det < 0) det = -det;
    int64_t capacity = (det + 1) * 32;
    const char **labels;
    int64_t *lengths, *kparams;
    SP_TRY(preferred_arrays(preferred, npreferred, &labels, &lengths, &kparams));
    int64_t *families = malloc(sizeof(int64_t) * 13 * (size_t)capacity);
    int64_t *occurrences = malloc(sizeof(int64_t) * 17 * (size_t)capacity);
    if (!families || !occurrences) {
        free(families);
        free(occurrences);
        return sp_fail(PY_MEMORY, "k family search ran out of memory");
    }
    int64_t occurrence_count = 0;
    int family_count = k_family_records(0, frame->sg, frame->basis, npreferred, labels, lengths, kparams, (int)capacity,
                                        (int)capacity, families, occurrences, &occurrence_count);
    if (family_count < 0) {
        free(families);
        free(occurrences);
        return k_family_fail(family_count);
    }
    int *family_occ_start = sp_alloc(sizeof(int) * (size_t)(family_count + 1));
    int64_t *sub_families = sp_alloc(sizeof(int64_t) * 11 * (size_t)(family_count + 1));
    if (!family_occ_start || !sub_families) {
        free(families);
        free(occurrences);
        return 1;
    }
    for (int f = 0; f < family_count; f++) {
        const int64_t *family = families + 13 * f;
        int gid = (int)family[0];
        sp_str label = gid_k_label(gid);
        int64_t *s = sub_families + 11 * f;
        s[0] = gid;
        memcpy(s + 1, family + 1, 4 * sizeof(int64_t));
        for (int axis = 0; axis < 3; axis++) {
            rat value = rat_make(family[5 + 2 * axis], family[6 + 2 * axis]);
            s[5 + 2 * axis] = (int64_t)value.n;
            s[6 + 2 * axis] = (int64_t)value.d;
        }
        family_occ_start[f] = sp_occ_count;
        for (int64_t o = family[11]; o < family[11] + family[12]; o++) {
            const int64_t *row = occurrences + 17 * o;
            sp_occ occ;
            occ.mesh = row[0];
            memcpy(occ.kp, row + 1, 4 * sizeof(int64_t));
            for (int axis = 0; axis < 3; axis++) occ.pml[axis] = rat_make(row[5 + 2 * axis], row[6 + 2 * axis]);
            memcpy(occ.bucket, row + 11, 6 * sizeof(int64_t));
            occ.k_label = label;
            if (occ_push(&occ) < 0) {
                free(families);
                free(occurrences);
                return 1;
            }
        }
    }
    int *family_occ_count = sp_alloc(sizeof(int) * (size_t)(family_count + 1));
    if (!family_occ_count) {
        free(families);
        free(occurrences);
        return 1;
    }
    for (int f = 0; f < family_count; f++) family_occ_count[f] = (int)families[13 * f + 12];
    free(families);
    free(occurrences);
    int op_count;
    int *ops = constrained_ops(frame, &op_count);
    if (!ops) return 1;
    int status = subduction_rows(frame->sg, family_count, sub_families, npreferred, labels, lengths, kparams, op_count, ops,
                                 frame->magnetic, 1);
    if (status) {
        int64_t gid, code;
        subduction_rows_failure(&gid, &code);
        char message[160];
        if (status == 3) {
            snprintf(message, sizeof message, "subduction directions failed for gid %lld (%lld)", (long long)gid, (long long)code);
            return sp_fail(PY_VALUE, message);
        }
        if (status == 4) {
            static const struct {
                int kind;
                const char *message;
            } K_RECORD[] = {{PY_RUNTIME, "reciprocal orbit key failed"},
                            {PY_VALUE, "k-vector Source tables are inconsistent"},
                            {PY_VALUE, "k label expects more parameters"},
                            {PY_VALUE, "missing little_k_star_conv2ml pointer"},
                            {PY_ARITH, "k-vector value out of range"},
                            {PY_VALUE, "newlat failed to find a basis"}};
            if (code >= 1 && code <= 5) return sp_fail(K_RECORD[code].kind, K_RECORD[code].message);
            snprintf(message, sizeof message, "reciprocal orbit key failed (%lld)", (long long)code);
            return sp_fail(PY_RUNTIME, message);
        }
        if (status == 5) return sp_fail(PY_ZERO, "k parameter has a zero denominator");
        snprintf(message, sizeof message, "subduction rows failed (%d)", status);
        return sp_fail(status == 2 ? PY_MEMORY : PY_VALUE, message);
    }
    const int64_t *rows, *occ_families;
    int64_t occ_family_count;
    int64_t count = subduction_rows_ints(&rows, &occ_families, &occ_family_count);
    const double *directions = subduction_rows_directions();
    int64_t at = 0;
    for (int64_t r = 0; r < count; r++) {
        const int64_t *row = rows + 6 * r;
        int gid = (int)row[0], family = (int)row[1], dcount = (int)row[4], dim = (int)row[5];
        sp_row out_row;
        memset(&out_row, 0, sizeof out_row);
        out_row.gid = gid;
        out_row.irrep_name = sp_irrep_full_label(gid);
        out_row.k_label = gid_k_label((int)sub_families[11 * family]);
        memcpy(out_row.kp, sub_families + 11 * family + 1, 4 * sizeof(int64_t));
        for (int axis = 0; axis < 3; axis++)
            out_row.pml[axis] = rat_make(sub_families[11 * family + 5 + 2 * axis], sub_families[11 * family + 6 + 2 * axis]);
        out_row.dim = dim;
        out_row.count = dcount;
        double *dir = sp_alloc(sizeof(double) * (size_t)(dim * dcount + 1));
        if (!dir) return 1;
        for (int c = 0; c < dim; c++)
            for (int f = 0; f < dcount; f++) dir[c * dcount + f] = directions[at + (int64_t)f * dim + c];
        at += (int64_t)dcount * dim;
        out_row.dir = dir;
        int total = 0;
        for (int64_t i = row[2]; i < row[2] + row[3]; i++) total += family_occ_count[occ_families[i]];
        int *occ = sp_alloc(sizeof(int) * (size_t)(total + 1));
        if (!occ) return 1;
        int n = 0;
        for (int64_t i = row[2]; i < row[2] + row[3]; i++) {
            int f = (int)occ_families[i];
            for (int k = 0; k < family_occ_count[f]; k++) occ[n++] = family_occ_start[f] + k;
        }
        out_row.occ = occ;
        out_row.nocc = n;
        SP_TRY(sp_rows_push(out, &out_row));
    }
    return 0;
}

static int alias_seen(const int64_t (*seen)[4], int n, const int64_t *kp)
{
    for (int i = 0; i < n; i++)
        if (sp_kp_identity_eq(seen[i], kp)) return 1;
    return 0;
}

typedef struct {
    const sp_frame *frame;
    const sp_row *row;
    int *ops, op_count, full_dim;
    int64_t (*seen)[4];
    int nseen;
    sp_rows *out;
} alias_ctx;

static int append_alias(alias_ctx *c, const int64_t *kp, int occurrence)
{
    if (alias_seen((const int64_t (*)[4])c->seen, c->nseen, kp)) return 0;
    memcpy(c->seen[c->nseen++], kp, 4 * sizeof(int64_t));
    double *buffer = sp_alloc(sizeof(double) * (size_t)(c->full_dim * c->full_dim + 1));
    if (!buffer) return 1;
    int count = subduction_directions(c->row->gid, 0, c->frame->magnetic, 4, kp, c->op_count, c->ops, buffer);
    if (count < 0) {
        char message[160];
        snprintf(message, sizeof message, "subduction directions failed for gid %d (%d)", c->row->gid, -count);
        return sp_fail(PY_VALUE, message);
    }
    if (count <= 0) return 0;
    sp_row alias;
    memset(&alias, 0, sizeof alias);
    alias.gid = c->row->gid;
    alias.irrep_name = c->row->irrep_name;
    alias.k_label = c->row->k_label;
    memcpy(alias.kp, kp, 4 * sizeof(int64_t));
    memcpy(alias.pml, sp_occ_pool[occurrence].pml, sizeof alias.pml);
    alias.dim = c->full_dim;
    alias.count = count;
    double *dir = sp_alloc(sizeof(double) * (size_t)(c->full_dim * count + 1));
    int *occ = sp_alloc(sizeof(int));
    if (!dir || !occ) return 1;
    for (int coordinate = 0; coordinate < c->full_dim; coordinate++)
        for (int f = 0; f < count; f++) dir[coordinate * count + f] = buffer[f * c->full_dim + coordinate];
    alias.dir = dir;
    occ[0] = occurrence;
    alias.occ = occ;
    alias.nocc = 1;
    return sp_rows_push(c->out, &alias);
}

static int occurrence_source_kparam(int sg, int gid, int occurrence, int64_t *out, int *found)
{
    int64_t words[6];
    for (int axis = 0; axis < 3; axis++) {
        words[2 * axis] = (int64_t)sp_occ_pool[occurrence].pml[axis].n;
        words[2 * axis + 1] = (int64_t)sp_occ_pool[occurrence].pml[axis].d;
    }
    int status = k_family_source_kparam(sg, gid, words, out);
    if (status < 0) return k_family_fail(status);
    *found = status;
    return 0;
}

int sp_aliases(const sp_frame *frame, const sp_row *row, const int64_t *requested, sp_rows *out)
{
    SP_TRY(load());
    alias_ctx c = {frame, row, 0, 0, sp_little_full_dim(row->gid), 0, 0, out};
    c.ops = constrained_ops(frame, &c.op_count);
    c.seen = sp_alloc(sizeof(int64_t) * 4 * (size_t)(row->nocc + 3));
    if (!c.ops || !c.seen) return 1;
    memcpy(c.seen[c.nseen++], row->kp, 4 * sizeof(int64_t));
    if (!frame->magnetic) {
        for (int i = 0; i < row->nocc; i++) {
            int64_t kp[4];
            int found;
            SP_TRY(occurrence_source_kparam(frame->sg, row->gid, row->occ[i], kp, &found));
            if (!found) continue;
            SP_TRY(append_alias(&c, kp, row->occ[i]));
        }
        if (requested) {
            for (int i = 0; i < row->nocc; i++)
                if (sp_kp_exact_eq(sp_occ_pool[row->occ[i]].kp, requested)) {
                    SP_TRY(append_alias(&c, requested, row->occ[i]));
                    break;
                }
        }
        return 0;
    }
    int exact = -1;
    if (requested)
        for (int i = 0; i < row->nocc && exact < 0; i++)
            if (sp_kp_exact_eq(sp_occ_pool[row->occ[i]].kp, requested)) exact = row->occ[i];
    int64_t (*occ_kp)[4] = sp_alloc(sizeof(int64_t) * 4 * (size_t)(row->nocc + 1));
    int *occ_index = sp_alloc(sizeof(int) * (size_t)(row->nocc + 1));
    if (!occ_kp || !occ_index) return 1;
    int n = 0;
    for (int i = 0; i < row->nocc; i++) {
        int64_t kp[4];
        int found;
        SP_TRY(occurrence_source_kparam(frame->sg, row->gid, row->occ[i], kp, &found));
        if (!found) continue;
        memcpy(occ_kp[n], kp, sizeof kp);
        occ_index[n++] = row->occ[i];
        if (exact >= 0 && sp_kp_identity_eq(kp, requested)) {
            SP_TRY(append_alias(&c, requested, exact));
        } else {
            SP_TRY(append_alias(&c, kp, row->occ[i]));
        }
    }
    if (requested && !alias_seen((const int64_t (*)[4])c.seen, c.nseen, requested)) {
        int match = exact;
        if (match < 0)
            for (int i = 0; i < n && match < 0; i++)
                if (sp_kp_periodic_eq(occ_kp[i], requested)) match = occ_index[i];
        if (match >= 0) SP_TRY(append_alias(&c, requested, match));
    }
    return 0;
}

static void star_vectors(int gid, const int64_t *kp, rat (*out)[3], int *count)
{
    int sg = (int)section_int(&T.irr_sg, gid - 1);
    int64_t slot = (int64_t)(sg - 1) * 27 + section_int(&T.irr_k, gid - 1) - 1;
    int n = (int)section_int(&T.star_count, slot);
    int64_t pointer = section_int(&T.star_pointer, slot);
    rat params[3];
    for (int i = 0; i < 3; i++) params[i] = kp[3] ? rat_make(kp[i], kp[3]) : rat_int(kp[i]);
    for (int arm = 0; arm < n; arm++) {
        rat pieces[4][3];
        for (int piece = 0; piece < 4; piece++) {
            int64_t base = 16 * (pointer - 1 + arm) + 4 * piece, den = section_int(&T.star, base + 3);
            for (int axis = 0; axis < 3; axis++) pieces[piece][axis] = den ? rat_make(section_int(&T.star, base + axis), den) : rat_int(0);
        }
        for (int axis = 0; axis < 3; axis++) {
            rat value = pieces[0][axis];
            for (int i = 0; i < 3; i++) value = rat_add(value, rat_mul(params[i], pieces[i + 1][axis]));
            out[arm][axis] = value;
        }
    }
    *count = n;
}

static void column_major(const int64_t *values, int64_t (*out)[3])
{
    for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) out[row][col] = values[row + 3 * col];
}

static void rotate_raw(int sg, const rat *vector, const int64_t *point, rat *out)
{
    int lattice = geo_sg_lattice(sg);
    int64_t record[36];
    for (int i = 0; i < 36; i++) record[i] = section_int(&T.lattice_ml, (int64_t)(lattice - 1) * 36 + i);
    int64_t recip_den = section_int(&T.lattice_ml_denom, (int64_t)(lattice - 1) * 2 + 1);
    int64_t a[3][3], p[3][3], b[3][3], ap[3][3], combined[3][3];
    column_major(record + 18, a);
    column_major(point, p);
    column_major(record + 27, b);
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            ap[r][c] = 0;
            for (int k = 0; k < 3; k++) ap[r][c] += a[r][k] * p[k][c];
        }
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            combined[r][c] = 0;
            for (int k = 0; k < 3; k++) combined[r][c] += ap[r][k] * b[k][c];
        }
    i128 common = 1;
    for (int axis = 0; axis < 3; axis++) common = common / rat_gcd(common, vector[axis].d) * vector[axis].d;
    i128 xyz[3];
    for (int axis = 0; axis < 3; axis++) xyz[axis] = vector[axis].n * (common / vector[axis].d);
    for (int axis = 0; axis < 3; axis++) {
        i128 num = combined[axis][0] * xyz[0] + combined[axis][1] * xyz[1] + combined[axis][2] * xyz[2];
        out[axis] = rat_make(num, common * recip_den);
    }
}

typedef struct {
    uint64_t *keys;
    size_t cap;
} dead_set;

static int dead_has(const dead_set *d, uint64_t key)
{
    for (size_t i = (size_t)(key * 0x9E3779B97F4A7C15ull) & (d->cap - 1);; i = (i + 1) & (d->cap - 1)) {
        if (!d->keys[i]) return 0;
        if (d->keys[i] == key) return 1;
    }
}

static void dead_add(dead_set *d, uint64_t key)
{
    for (size_t i = (size_t)(key * 0x9E3779B97F4A7C15ull) & (d->cap - 1);; i = (i + 1) & (d->cap - 1)) {
        if (!d->keys[i] || d->keys[i] == key) {
            d->keys[i] = key;
            return;
        }
    }
}

static int bijective(int arm, int n, const unsigned char (*match)[48], uint64_t used, dead_set *dead, size_t *dead_count)
{
    if (arm == n) return 1;
    uint64_t key = (used << 6 | (uint64_t)arm) + 1;
    if (dead_has(dead, key)) return 0;
    for (int s = 0; s < n; s++) {
        if (!match[arm][s] || (used & (1ull << s))) continue;
        if (bijective(arm + 1, n, match, used | (1ull << s), dead, dead_count)) return 1;
    }
    if (*dead_count * 2 < dead->cap) {
        dead_add(dead, key);
        (*dead_count)++;
    }
    return 0;
}

static char standard_message[160];

const char *sp_standard_message(void) { return standard_message; }

int sp_standard_kparam(int sg, int gid, const int64_t *kp, int64_t *out, int *ok, int *raised_kind)
{
    *ok = 0;
    *raised_kind = PY_NONE;
    SP_TRY(load());
    int64_t selection[6];
    int status = k_family_standard_selection(sg, gid, kp, selection);
    if (status < 0) {
        int code = -status;
        *raised_kind = code >= 1 && code <= 7 ? K_FAMILY_ERRORS[code].kind : PY_RUNTIME;
        if (code >= 1 && code <= 7)
            snprintf(standard_message, sizeof standard_message, "%s", K_FAMILY_ERRORS[code].message);
        else
            snprintf(standard_message, sizeof standard_message, "k family search failed (%d)", status);
        return 0;
    }
    int64_t point[9];
    for (int i = 0; i < 9; i++) point[i] = section_int(&T.point_op, 9 * (selection[5] - 1) + i);
    rat source[48][3], standard[48][3];
    int ns, nt;
    star_vectors(gid, kp, source, &ns);
    star_vectors(gid, selection, standard, &nt);
    if (ns != nt || ns > 48) {
        *raised_kind = PY_VALUE;
        return 0;
    }
    unsigned char match[48][48];
    memset(match, 0, sizeof match);
    for (int arm = 0; arm < ns; arm++) {
        rat rotated[3];
        rotate_raw(sg, source[arm], point, rotated);
        int any = 0;
        for (int conj = 0; conj < 2; conj++)
            for (int s = 0; s < nt; s++) {
                int integral = 1;
                for (int axis = 0; axis < 3; axis++) {
                    rat delta = conj ? rat_add(rotated[axis], standard[s][axis]) : rat_sub(rotated[axis], standard[s][axis]);
                    integral &= delta.d == 1;
                }
                if (integral) {
                    match[arm][s] = 1;
                    any = 1;
                }
            }
        if (!any) {
            *raised_kind = PY_VALUE;
            return 0;
        }
    }
    dead_set dead = {calloc(1u << 16, sizeof(uint64_t)), 1u << 16};
    if (!dead.keys) return sp_fail(PY_MEMORY, "spec stage ran out of memory");
    size_t dead_count = 0;
    int found = bijective(0, ns, (const unsigned char (*)[48])match, 0, &dead, &dead_count);
    free(dead.keys);
    if (!found) {
        *raised_kind = PY_VALUE;
        return 0;
    }
    memcpy(out, selection, 4 * sizeof(int64_t));
    *ok = 1;
    return 0;
}

int sp_matrix_rank(const double *m, int rows, int cols, double tol, int *rank)
{
    *rank = 0;
    if (rows <= 0 || cols <= 0) return 0;
    int k = rows < cols ? rows : cols;
    double *s = malloc(sizeof(double) * (size_t)(k + 1));
    if (!s || la_singular_values(rows, cols, m, s)) {
        free(s);
        return sp_fail(PY_MEMORY, "spec stage ran out of memory");
    }
    int out = 0;
    for (int i = 0; i < k; i++) out += s[i] > tol;
    free(s);
    *rank = out;
    return 0;
}

int sp_lstsq(const double *a_in, int rows, int acols, const double *b_in, int bcols, double *x, int *rank)
{
    double rcond = DBL_EPSILON * (rows > acols ? rows : acols);
    if (la_lstsq(rows, acols, bcols, a_in, b_in, rcond, x, rank)) return sp_fail(PY_MEMORY, "spec stage ran out of memory");
    return 0;
}

void sp_matmul(const double *a, int rows, int inner, const double *b, int cols, double *out)
{
    if (rows <= 0 || cols <= 0) return;
    la_dgemm(0, 0, rows, cols, inner, a, inner, b, cols, out, cols);
}
