/* modes.engine.trace.pipeline_trace.project_mode_block: the selected little
 * irrep's Wyckoff row, project entries, isotropy row and basis, Source
 * order-parameter rows, presentation setting frame and materialized atom
 * layouts around one mode_block call. */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "tables.h"
int isotropy_row_subgroup(int64_t row_id);

int64_t isotropy_row_count(void);
int isotropy_row_dim(int64_t row_id);
int isotropy_row_free(int64_t row_id);

const char *irrep_full_label(int gid);
const char *irrep_old_label(int old_id);

int mode_block(int sg, int row_id, int site_pg, int gid, int irrep_type, int full_dim, int old_id, int entry_count,
               const int *pg_irreps, int param_count, const int64_t *params, const double *inverse, int inverse_error,
               int vector_setting, int row_total, int row_stride, const int *row_lengths, const double *rows,
               int has_rows, int source_row_id, int freeparam, int flags, int direction_rows, const double *direction,
               int candidate_lanes, const int *candidate_records, int baseline_lanes, const int *baseline_records,
               int transform_count, const double *transforms, int negate);
int row_site_irreps(int gid, int row_id, int *out);
int64_t site_atom_table(int sg, int row_id, int site_pg, const int64_t *param_words, const int64_t *basis,
                        const int64_t *representative_words, int64_t capacity, int64_t *points, int64_t *records,
                        int64_t *order);
int case_k_records(int sg, int kslot, int param_count, const int64_t *params, int64_t *newlat_out, int64_t *entry_out,
                   int *entry_count);
int newlat_basis(int mode_count, int record_count, const int64_t *records, int64_t *out);
int irrep_transform_inverse(int gid, double complex *inverse);

/* Python exception classes (as spec_core.h). */
enum {
    PY_VALUE = 1,
    PY_KEY = 2,
    PY_INDEX = 3,
    PY_TYPE = 4,
    PY_ZERO = 5,
    PY_ARITH = 6,
    PY_MEMORY = 7,
    PY_RUNTIME = 8,
};

static int fail_kind;
static char fail_message[200];

static int fail(int kind, const char *message)
{
    if (!fail_kind) {
        fail_kind = kind;
        snprintf(fail_message, sizeof fail_message, "%s", message);
    }
    return 1;
}

static struct {
    section_view w_pointer, w_count, w_label, w_pg, vector_reps, irr_sg, irr_old, full_dim, irr_type, irr_k, k_label,
        irr_count, iso_irrep_pointer, iso_direction, iso_basis, iso_origin,
        iso_op_pointer, iso_op, const_values, sg_point_group, pg_order, pg_gen_count, pg_gen, lattice;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (geo_load() || section_find("wyckoff.iwyckoff_pointer", &T.w_pointer) || section_find("wyckoff.iwyckoff_count", &T.w_count) ||
        section_find("wyckoff.wyckoff_label", &T.w_label) || section_find("wyckoff.iwyckoff_pg", &T.w_pg) ||
        section_find("wyckoff.iwyckoff_pg_vector_reps", &T.vector_reps) ||
        section_find("little.little_irr_space_group", &T.irr_sg) || section_find("little.little_irr_old", &T.irr_old) ||
        section_find("little.little_irr_full_dim", &T.full_dim) ||
        section_find("little.little_irr_type", &T.irr_type) || section_find("little.little_irr_k", &T.irr_k) ||
        section_find("little.little_k_label", &T.k_label) || section_find("little.little_irr_count", &T.irr_count) ||
        section_find("isotropy.isotropy_irrep_pointer", &T.iso_irrep_pointer) ||
        
        section_find("isotropy.isotropy_direction", &T.iso_direction) || section_find("isotropy.isotropy_basis", &T.iso_basis) ||
        section_find("isotropy.isotropy_origin", &T.iso_origin) ||
        section_find("isotropy.isotropy_orderparam_pointer", &T.iso_op_pointer) ||
        section_find("isotropy.isotropy_orderparam", &T.iso_op) || section_find("const.values", &T.const_values) ||
        section_find("space.ispace_point_group", &T.sg_point_group) || section_find("space.ipoint_group_order", &T.pg_order) ||
        section_find("space.ipoint_group_gen_count", &T.pg_gen_count) || section_find("space.ipoint_group_gen", &T.pg_gen) ||
        section_find("space.ispace_lattice", &T.lattice))
        return fail(PY_VALUE, "mode frame Source tables are missing");
    T.ready = 1;
    return 0;
}

typedef struct {
    const char *text;
    int64_t length;
} str;

static str strip(const section_view *view, int64_t index)
{
    str out;
    if (section_string(view, index, &out.text, &out.length)) return (str){"", 0};
    while (out.length && (out.text[0] == ' ' || out.text[0] == '\t')) {
        out.text++;
        out.length--;
    }
    return out;
}

static int str_eq(str a, const char *b, int64_t blen) { return a.length == blen && !memcmp(a.text, b, (size_t)blen); }

/* k_label_from_irrep: the leading upper-case run, 0 when there is none. */
static int k_prefix(str label, str *out)
{
    int64_t n = 0;
    while (n < label.length && label.text[n] >= 'A' && label.text[n] <= 'Z') n++;
    *out = (str){label.text, n};
    return n > 0;
}

typedef struct {
    int row_id, subgroup, direction;
    int64_t basis[9], origin[4];
} iso_row;

static int iso_rows(int gid, iso_row **out, int *count)
{
    *out = 0;
    *count = 0;
    int old_id = (int)section_int(&T.irr_old, gid - 1);
    if (old_id <= 0) return 0;
    int64_t pointers = section_count(&T.iso_irrep_pointer);
    if (old_id < 1 || old_id >= pointers) {
        char message[120];
        snprintf(message, sizeof message, "old irrep id out of isotropy range: %d", old_id);
        return fail(PY_KEY, message);
    }
    int start = (int)section_int(&T.iso_irrep_pointer, old_id - 1), end = (int)section_int(&T.iso_irrep_pointer, old_id);
    iso_row *rows = calloc((size_t)(end - start + 1), sizeof *rows);
    if (!rows) return fail(PY_MEMORY, "mode frame ran out of memory");
    for (int row_id = start; row_id < end; row_id++) {
        iso_row *r = &rows[row_id - start];
        int64_t idx = row_id - 1;
        r->row_id = row_id;
        r->subgroup = isotropy_row_subgroup(row_id);
        r->direction = (int)section_int(&T.iso_direction, idx);
        for (int i = 0; i < 9; i++) r->basis[i] = section_int(&T.iso_basis, idx * 9 + i);
        for (int i = 0; i < 4; i++) r->origin[i] = section_int(&T.iso_origin, idx * 4 + i);
    }
    *out = rows;
    *count = end - start;
    return 0;
}

static int64_t abs_det(const int64_t *b)
{
    int64_t det = b[0] * (b[4] * b[8] - b[5] * b[7]) - b[1] * (b[3] * b[8] - b[5] * b[6]) + b[2] * (b[3] * b[7] - b[4] * b[6]);
    return det < 0 ? -det : det;
}

static int canonical_row(const iso_row *rows, int count)
{
    int best = -1;
    int64_t best_key[5] = {0};
    for (int i = 0; i < count; i++) {
        const iso_row *r = &rows[i];
        int pg = (int)section_int(&T.sg_point_group, r->subgroup - 1);
        int gens = (int)section_int(&T.pg_gen_count, pg - 1), inversion = 0;
        for (int g = 0; g < gens; g++) {
            int64_t op = section_int(&T.pg_gen, (int64_t)(pg - 1) * 5 + g);
            inversion |= op == 25 || op == 61;
        }
        int64_t key[5] = {isotropy_row_free(r->row_id), abs_det(r->basis), -section_int(&T.pg_order, pg - 1),
                          inversion ? 0 : 1, r->row_id};
        int less = best < 0;
        for (int k = 0; k < 5 && !less; k++) {
            if (key[k] != best_key[k]) {
                less = key[k] < best_key[k];
                break;
            }
        }
        if (less) {
            best = i;
            memcpy(best_key, key, sizeof key);
        }
    }
    return best;
}

/* ------------------------------------------------------------ atom layouts */

typedef struct {
    int count;
    int64_t *points; /* count x 3 rationals as int128 (n, d) words after the setting transform */
    int *records;    /* count x 5 */
    int64_t serial;  /* identifies this materialization to the Python boundary's decoded copy */
} layout;

typedef struct {
    int64_t key[80];
    int key_length;
    layout value;
} layout_entry;

enum { LAYOUT_CACHE = 4096, LAYOUT_BUCKETS = 8192 };
static layout_entry layout_cache[LAYOUT_CACHE];
static int layout_cache_count, layout_cache_next;
static int layout_bucket[LAYOUT_BUCKETS];
static int64_t layout_serial;

static uint64_t key_hash(const int64_t *key, int n)
{
    uint64_t h = 1469598103934665603ull;
    for (int i = 0; i < n; i++) h = (h ^ (uint64_t)key[i]) * 1099511628211ull;
    return h;
}

static int layout_find(const int64_t *key, int n)
{
    for (uint64_t i = key_hash(key, n) & (LAYOUT_BUCKETS - 1);; i = (i + 1) & (LAYOUT_BUCKETS - 1)) {
        int slot = layout_bucket[i] - 1;
        if (slot < 0) return -1;
        if (layout_cache[slot].key_length == n && !memcmp(layout_cache[slot].key, key, sizeof(int64_t) * (size_t)n)) return slot;
    }
}

static void layout_rehash(void)
{
    memset(layout_bucket, 0, sizeof layout_bucket);
    for (int slot = 0; slot < layout_cache_count; slot++) {
        uint64_t i = key_hash(layout_cache[slot].key, layout_cache[slot].key_length) & (LAYOUT_BUCKETS - 1);
        while (layout_bucket[i]) i = (i + 1) & (LAYOUT_BUCKETS - 1);
        layout_bucket[i] = slot + 1;
    }
}

typedef struct {
    int has;
    i128 ints[3][3], origin[3];
    i128 scale;
} setting_frame;

static int materialize(int sg, int row_id, int site_pg, const int64_t *site_words, int nsite, const int64_t *basis,
                       const setting_frame *frame, const int64_t *setting_key, layout *out, int *raised)
{
    *raised = 0;
    int64_t key[80];
    int n = 0;
    key[n++] = sg;
    key[n++] = row_id;
    key[n++] = nsite;
    for (int i = 0; i < 4 * nsite; i++) key[n++] = site_words[i];
    for (int i = 0; i < 9; i++) key[n++] = basis[i];
    for (int i = 0; i < 24; i++) key[n++] = setting_key[i];
    int hit = layout_find(key, n);
    if (hit >= 0) {
        *out = layout_cache[hit].value;
        return 0;
    }
    int64_t param_words[12] = {0};
    for (int i = 0; i < 3; i++) {
        if (i < nsite) {
            memcpy(param_words + 4 * i, site_words + 4 * i, 4 * sizeof(int64_t));
        } else {
            write_i128(param_words + 4 * i, 0);
            write_i128(param_words + 4 * i + 2, 1);
        }
    }
    int64_t det = abs_det(basis);
    int64_t capacity = 48 * (det > 1 ? det : 1);
    int64_t *points = malloc(sizeof(int64_t) * 12 * (size_t)capacity), *records = malloc(sizeof(int64_t) * 5 * (size_t)capacity);
    if (!points || !records) {
        free(points);
        free(records);
        return fail(PY_MEMORY, "mode frame ran out of memory");
    }
    int64_t count = site_atom_table(sg, row_id, site_pg, param_words, basis, 0, capacity, points, records, 0);
    if (count < 0) {
        static const struct {
            int kind;
            const char *message;
        } errors[] = {{PY_RUNTIME, "site atom table failed"},
                      {PY_VALUE, "site atom table Source tables are inconsistent"},
                      {PY_ARITH, "site atom table value out of range"},
                      {PY_VALUE, "singular basis matrix"},
                      {PY_VALUE, "parent representative has no Source orbit transport"},
                      {PY_VALUE, "Source site cosets overlap"},
                      {PY_VALUE, "Source atom operations are not unique in the supercell"},
                      {PY_VALUE, "Source representative order is not an atom permutation"},
                      {PY_KEY, "rebased atom operation has no Source atom"}};
        int code = (int)-count;
        free(points);
        free(records);
        *raised = 1;
        if (code >= 1 && code <= 8) return fail(errors[code].kind, errors[code].message);
        char message[80];
        snprintf(message, sizeof message, "site atom table failed (%d)", code);
        return fail(PY_RUNTIME, message);
    }
    layout value = {(int)count, malloc(sizeof(int64_t) * 12 * (size_t)(count + 1)), malloc(sizeof(int) * 5 * (size_t)(count + 1)),
                    ++layout_serial};
    if (!value.points || !value.records) {
        free(points);
        free(records);
        return fail(PY_MEMORY, "mode frame ran out of memory");
    }
    for (int64_t a = 0; a < count; a++) {
        rat p[3];
        read_rat3(points + 12 * a, p);
        rat q[3];
        if (frame->has) {
            i128 den = 1;
            for (int axis = 0; axis < 3; axis++) den = den / rat_gcd(den, p[axis].d) * p[axis].d;
            i128 nums[3];
            for (int axis = 0; axis < 3; axis++) nums[axis] = p[axis].n * (den / p[axis].d);
            for (int col = 0; col < 3; col++) {
                i128 sum = 0;
                for (int row = 0; row < 3; row++) sum += nums[row] * frame->ints[row][col];
                q[col] = rat_make(sum + frame->origin[col] * den, den * frame->scale);
            }
        } else {
            memcpy(q, p, sizeof q);
        }
        write_rat3(value.points + 12 * a, q);
        for (int k = 0; k < 5; k++) value.records[5 * a + k] = (int)records[5 * a + k];
    }
    free(points);
    free(records);
    if (rat_overflow) return fail(PY_ARITH, "mode frame atom coordinate out of range");
    layout_entry *slot = &layout_cache[layout_cache_next];
    if (layout_cache_count == LAYOUT_CACHE) {
        free(slot->value.points);
        free(slot->value.records);
    } else {
        layout_cache_count++;
    }
    memcpy(slot->key, key, sizeof(int64_t) * (size_t)n);
    slot->key_length = n;
    slot->value = value;
    layout_cache_next = (layout_cache_next + 1) % LAYOUT_CACHE;
    layout_rehash();
    *out = value;
    return 0;
}

/* irrep_transform_inverse per gid, memoized like the Python decoder's. */
typedef struct {
    int gid, status;
    double complex *inverse;
} inverse_entry;

enum { INVERSE_CACHE = 1024 };
static inverse_entry inverse_cache[INVERSE_CACHE];
static int inverse_cache_count, inverse_cache_next;

static int cached_transform_inverse(int gid, int dim, double complex *out)
{
    for (int i = 0; i < inverse_cache_count; i++)
        if (inverse_cache[i].gid == gid) {
            if (!inverse_cache[i].status) memcpy(out, inverse_cache[i].inverse, sizeof(double complex) * (size_t)(dim * dim));
            return inverse_cache[i].status;
        }
    int status = irrep_transform_inverse(gid, out);
    inverse_entry *e = &inverse_cache[inverse_cache_next];
    if (inverse_cache_count == INVERSE_CACHE) free(e->inverse);
    else inverse_cache_count++;
    inverse_cache_next = (inverse_cache_next + 1) % INVERSE_CACHE;
    e->gid = gid;
    e->status = status;
    e->inverse = 0;
    if (!status) {
        e->inverse = malloc(sizeof(double complex) * (size_t)(dim * dim + 1));
        if (e->inverse) memcpy(e->inverse, out, sizeof(double complex) * (size_t)(dim * dim));
        else e->gid = 0;
    }
    return status;
}

/* ------------------------------------------------------------ k-family cache */

typedef struct {
    int sg;
    char label[24];
    int64_t label_length;
    int count, failed, fail_kind;
    char message[200];
    int *littles;
    int has_canonical;
    int64_t canonical[9];
} k_family;

enum { K_FAMILY_CACHE = 512 };
static k_family k_families[K_FAMILY_CACHE];
static int k_family_count, k_family_next;

/* little_records_for_k(sg, k_label) and the first of its isotropy bases. */
static int k_family_of(int sg, const char *k_label, int64_t k_label_length, const k_family **out)
{
    for (int i = 0; i < k_family_count; i++) {
        const k_family *f = &k_families[i];
        if (f->sg == sg && f->label_length == k_label_length && !memcmp(f->label, k_label, (size_t)k_label_length)) {
            if (f->failed) return fail(f->fail_kind, f->message);
            *out = f;
            return 0;
        }
    }
    if (k_label_length >= (int64_t)sizeof k_families[0].label) return fail(PY_VALUE, "k label is too long");
    k_family *f = &k_families[k_family_next];
    if (k_family_count == K_FAMILY_CACHE) free(f->littles);
    else k_family_count++;
    k_family_next = (k_family_next + 1) % K_FAMILY_CACHE;
    memset(f, 0, sizeof *f);
    f->sg = sg;
    memcpy(f->label, k_label, (size_t)k_label_length);
    f->label_length = k_label_length;
    int64_t total = section_count(&T.irr_sg);
    f->littles = malloc(sizeof(int) * 4096);
    int *seen = malloc(sizeof(int) * 4096);
    if (!f->littles || !seen) {
        free(seen);
        f->sg = -1;
        return fail(PY_MEMORY, "mode frame ran out of memory");
    }
    int nseen = 0;
    for (int64_t index = 0; index < total && !f->failed; index++) {
        if (section_int(&T.irr_sg, index) != sg) continue;
        int gid = (int)index + 1, old = (int)section_int(&T.irr_old, index), key = old > 0 ? old : -gid, dup = 0;
        for (int i = 0; i < nseen && !dup; i++) dup = seen[i] == key;
        if (dup) continue;
        str label = {irrep_full_label(gid), (int64_t)strlen(irrep_full_label(gid))}, prefix;
        if (!k_prefix(label, &prefix)) {
            f->failed = 1;
            f->fail_kind = PY_VALUE;
            snprintf(f->message, sizeof f->message, "cannot infer k label from irrep label '%.*s'", (int)label.length, label.text);
            break;
        }
        if (!str_eq(prefix, k_label, k_label_length)) continue;
        if (nseen < 4096) seen[nseen++] = key;
        if (f->count < 4096) f->littles[f->count++] = gid;
    }
    free(seen);
    for (int i = 0; i < f->count && !f->failed && !f->has_canonical; i++) {
        iso_row *rows;
        int count;
        if (iso_rows(f->littles[i], &rows, &count)) {
            f->failed = 1;
            f->fail_kind = fail_kind;
            snprintf(f->message, sizeof f->message, "%s", fail_message);
            fail_kind = 0;
            break;
        }
        if (count) {
            memcpy(f->canonical, rows[0].basis, sizeof f->canonical);
            f->has_canonical = 1;
        }
        free(rows);
    }
    if (f->failed) return fail(f->fail_kind, f->message);
    *out = f;
    return 0;
}

/* find_isotropy_subgroup_'s initial newlat_ basis of a case, memoized. */
typedef struct {
    int sg, kslot, nparams, status;
    int64_t params[6], basis[9];
} basis_entry;

enum { BASIS_CACHE = 256 };
static basis_entry basis_cache[BASIS_CACHE];
static int basis_cache_count, basis_cache_next;

static int find_isotropy_basis(int sg, int kslot, int nparams, const int64_t *param_words, int64_t *words, int64_t *newlat,
                               int64_t *entry_records, int64_t *basis)
{
    int n = nparams < 3 ? nparams : 3;
    for (int i = 0; i < basis_cache_count; i++) {
        const basis_entry *e = &basis_cache[i];
        if (e->sg == sg && e->kslot == kslot && e->nparams == nparams && !memcmp(e->params, param_words, sizeof(int64_t) * 2 * (size_t)n)) {
            memcpy(basis, e->basis, sizeof e->basis);
            return e->status;
        }
    }
    int entry_count = 0;
    int status = case_k_records(sg, kslot, nparams, words, newlat, entry_records, &entry_count);
    if (!status) status = newlat_basis(1, 4, newlat, basis);
    basis_entry *e = &basis_cache[basis_cache_next];
    if (basis_cache_count < BASIS_CACHE) basis_cache_count++;
    basis_cache_next = (basis_cache_next + 1) % BASIS_CACHE;
    e->sg = sg;
    e->kslot = kslot;
    e->nparams = nparams;
    e->status = status;
    memset(e->params, 0, sizeof e->params);
    memcpy(e->params, param_words, sizeof(int64_t) * 2 * (size_t)n);
    memcpy(e->basis, basis, sizeof e->basis);
    return status;
}

/* ------------------------------------------------------------ results */

static struct {
    int site_pg, gid, has_selected, row_id, subgroup, direction, has_origin, empty, adopted, lanes, source_rows, source_row_id;
    int64_t basis[9], origin[4];
    layout atoms;
    double *rows;
    int *row_lengths, row_stride;
} R;

static void reset_result(void)
{
    free(R.rows);
    free(R.row_lengths);
    memset(&R, 0, sizeof R);
}

/* out: site_pg, gid, has_selected, row_id, subgroup, direction,
 * has_origin, origin[4], basis[9], empty, adopted, lanes, atom count,
 * source row count, source row stride, source row id. */
void mode_frame_result(int64_t *out)
{
    int64_t values[] = {R.site_pg, R.gid, R.has_selected, R.row_id, R.subgroup, R.direction, R.has_origin, R.origin[0],
                        R.origin[1], R.origin[2], R.origin[3], R.basis[0], R.basis[1], R.basis[2], R.basis[3], R.basis[4],
                        R.basis[5], R.basis[6], R.basis[7], R.basis[8], R.empty, R.adopted, R.lanes, R.atoms.count,
                        R.source_rows, R.row_stride, R.source_row_id, R.atoms.serial};
    memcpy(out, values, sizeof values);
}

const int64_t *mode_frame_points(void) { return R.atoms.points; }
const int *mode_frame_records(void) { return R.atoms.records; }
const double *mode_frame_rows(void) { return R.rows; }
const int *mode_frame_row_lengths(void) { return R.row_lengths; }

const char *mode_frame_error(int *kind)
{
    *kind = fail_kind;
    return fail_message;
}

/* ------------------------------------------------------------ entry */

/* k_params / site_params / setting values arrive as (numerator, denominator)
 * int64 pairs (site_params as int128 words, 4 per value); override rows as row_count x stride with lengths. Returns
 * 0, a Python exception class, or -(100 + code) for a mode_block failure. */
int mode_frame(int sg, const char *wyckoff, int64_t wyckoff_length, const char *k_label, int64_t k_label_length, int nparams,
               const int64_t *param_words, int nsite, const int64_t *site_words, int selected_gid, int has_row_id,
               int selected_row_id, int has_basis_override, const int64_t *basis_override, int emit_all_rows,
               int emit_freeparam, int vector_setting, int override_count, int override_stride, const int *override_lengths,
               const double *override_rows, int has_occurrence_basis, const int64_t *occurrence_basis, int has_setting,
               const int64_t *setting_words, const int64_t *setting_origin_words)
{
    reset_result();
    fail_kind = 0;
    fail_message[0] = 0;
    if (load()) return fail_kind;
    rat_overflow = 0;
    char message[200];
    int w_start = (int)section_int(&T.w_pointer, sg - 1), w_count = (int)section_int(&T.w_count, sg - 1);
    int matches = 0, offset0 = -1;
    for (int o = 0; o < w_count; o++)
        if (str_eq(strip(&T.w_label, w_start + o - 1), wyckoff, wyckoff_length)) {
            matches++;
            offset0 = o;
        }
    if (matches != 1) {
        snprintf(message, sizeof message, "expected one Wyckoff row for SG%d %.*s, got %d", sg, (int)wyckoff_length, wyckoff, matches);
        return fail(PY_KEY, message), fail_kind;
    }
    int row_id = w_start + offset0, site_pg = (int)section_int(&T.w_pg, row_id - 1);
    R.site_pg = site_pg;
    const k_family *family;
    if (k_family_of(sg, k_label, k_label_length, &family)) return fail_kind;
    const int *littles = family->littles;
    int nlittles = family->count;
    int found = 0;
    for (int i = 0; i < nlittles; i++) found += littles[i] == selected_gid;
    if (selected_gid <= 0) {
        return fail(PY_VALUE, "selected gid must be a positive exact integer"), fail_kind;
    }
    if (found != 1) {
        snprintf(message, sizeof message, "gid %d does not identify one SG%d %.*s irrep", selected_gid, sg, (int)k_label_length,
                 k_label);
        return fail(PY_VALUE, message), fail_kind;
    }
    if (vector_setting != 1 && vector_setting != 2) {
        snprintf(message, sizeof message, "unsupported vector setting: %d", vector_setting);
        return fail(PY_VALUE, message), fail_kind;
    }
    int gid = selected_gid;
    R.gid = gid;
    int irrep_type = (int)section_int(&T.irr_type, gid - 1), full_dim = (int)section_int(&T.full_dim, gid - 1),
        old_id = (int)section_int(&T.irr_old, gid - 1);
    /* project entries: the selected little's site-pg irreps in the vector slots */
    int vector_slots[3];
    for (int i = 0; i < 3; i++) vector_slots[i] = (int)section_int(&T.vector_reps, (int64_t)(site_pg - 1) * 6 + (vector_setting - 1) * 3 + i);
    int *entries = malloc(sizeof(int) * 64);
    int nentries = 0, target_direction = 0, has_target = 0;
    if (!entries) {
        return fail(PY_MEMORY, "mode frame ran out of memory"), fail_kind;
    }
    {
        int present[12], count = row_site_irreps(gid, row_id, present);
        if (count < 0) {
            free(entries);
            snprintf(message, sizeof message, "no subduction for gid=%d", gid);
            return fail(PY_KEY, message), fail_kind;
        }
        for (int j = 0; j < count && nentries < 64; j++) {
            int pg_irrep = present[j];
            int in = 0;
            for (int i = 0; i < 3; i++) in |= vector_slots[i] == pg_irrep;
            if (!in) continue;
            entries[nentries++] = pg_irrep;
            target_direction = offset0 + 1;
            has_target = 1;
        }
    }
    int64_t find_basis[9];
    int has_find_basis = 0;
    if (nparams) {
        int lattice = (int)section_int(&T.lattice, sg - 1), kslot = 0, kmatches = 0;
        for (int slot = 1; slot <= 27; slot++)
            if (section_int(&T.irr_count, (int64_t)(sg - 1) * 27 + slot - 1) > 0 &&
                str_eq(strip(&T.k_label, (int64_t)(lattice - 1) * 27 + slot - 1), k_label, k_label_length)) {
                kslot = slot;
                kmatches++;
            }
        if (kmatches != 1) {
            char list[120] = "";
            for (int slot = 1, at = 0; slot <= 27; slot++)
                if (section_int(&T.irr_count, (int64_t)(sg - 1) * 27 + slot - 1) > 0 &&
                    str_eq(strip(&T.k_label, (int64_t)(lattice - 1) * 27 + slot - 1), k_label, k_label_length))
                    at += snprintf(list + at, sizeof list - (size_t)at, at ? ", %d" : "%d", slot);
                free(entries);
            snprintf(message, sizeof message, "expected one k label '%.*s' for SG%d, got [%s]", (int)k_label_length, k_label, sg, list);
            return fail(PY_KEY, message), fail_kind;
        }
        int64_t words[12], newlat[16], entry_records[32];
        for (int i = 0; i < nparams && i < 3; i++) {
            write_i128(words + 4 * i, param_words[2 * i]);
            write_i128(words + 4 * i + 2, param_words[2 * i + 1]);
        }
        int k_status = find_isotropy_basis(sg, kslot, nparams, param_words, words, newlat, entry_records, find_basis);
        if (k_status) {
            static const struct {
                int kind;
                const char *message;
            } errors[] = {{PY_RUNTIME, "k records failed"},
                          {PY_VALUE, "k-vector Source tables are inconsistent"},
                          {PY_VALUE, "k label expects more parameters"},
                          {PY_VALUE, "missing little_k_star_conv2ml pointer"},
                          {PY_ARITH, "k-vector value out of range"},
                          {PY_VALUE, "newlat failed to find a basis"}};
                free(entries);
            if (k_status >= 1 && k_status <= 5) return fail(errors[k_status].kind, errors[k_status].message), fail_kind;
            return fail(PY_RUNTIME, "k records failed"), fail_kind;
        }
        has_find_basis = 1;
    }
    int64_t canonical_basis[9];
    int has_canonical = family->has_canonical;
    memcpy(canonical_basis, family->canonical, sizeof canonical_basis);
    if (!has_canonical && has_find_basis) {
        memcpy(canonical_basis, find_basis, sizeof canonical_basis);
        has_canonical = 1;
    }
    iso_row *rows;
    int nrows;
    if (iso_rows(gid, &rows, &nrows)) {
        free(entries);
        return fail_kind;
    }
    int selected = -1, explicit_row = -1;
    if (has_row_id)
        for (int i = 0; i < nrows && explicit_row < 0; i++)
            if (rows[i].row_id == selected_row_id) explicit_row = i;
    if (explicit_row >= 0) {
        selected = explicit_row;
    } else if (!nparams) {
        selected = canonical_row(rows, nrows);
    } else if (has_target) {
        for (int i = 0; i < nrows && selected < 0; i++)
            if (rows[i].direction == target_direction) selected = i;
    }
    if (selected < 0) {
        for (int i = 0; i < nrows && selected < 0; i++)
            if (rows[i].origin[0] == 0 && rows[i].origin[1] == 0 && rows[i].origin[2] == 0 && rows[i].origin[3] == 1) selected = i;
        if (selected < 0 && nrows) selected = 0;
    }
    int has_override = override_count >= 0;
    int64_t selected_basis[9];
    int has_selected_basis = 1;
    if (has_basis_override && (explicit_row >= 0 || has_override))
        memcpy(selected_basis, basis_override, sizeof selected_basis);
    else if (selected >= 0)
        memcpy(selected_basis, rows[selected].basis, sizeof selected_basis);
    else if (has_find_basis)
        memcpy(selected_basis, find_basis, sizeof selected_basis);
    else
        has_selected_basis = 0;
    R.has_selected = has_selected_basis;
    if (has_selected_basis) {
        memcpy(R.basis, selected_basis, sizeof R.basis);
        if (selected >= 0) {
            R.row_id = rows[selected].row_id;
            R.subgroup = rows[selected].subgroup;
            R.direction = rows[selected].direction;
            R.has_origin = 1;
            memcpy(R.origin, rows[selected].origin, sizeof R.origin);
        }
    }
    /* Source order-parameter rows */
    int source_row_id = selected >= 0 ? rows[selected].row_id : 0, has_source = 1;
    if (has_override) {
        R.source_rows = override_count;
        R.row_stride = override_stride;
        R.rows = malloc(sizeof(double) * (size_t)(override_count * override_stride + 1));
        R.row_lengths = malloc(sizeof(int) * (size_t)(override_count + 1));
        if (!R.rows || !R.row_lengths) {
            free(rows);
            free(entries);
            return fail(PY_MEMORY, "mode frame ran out of memory"), fail_kind;
        }
        memcpy(R.rows, override_rows, sizeof(double) * (size_t)(override_count * override_stride));
        memcpy(R.row_lengths, override_lengths, sizeof(int) * (size_t)override_count);
        source_row_id = 0;
    } else if (source_row_id > 0) {
        int64_t idx = source_row_id - 1;
        if (idx >= isotropy_row_count()) {
            free(rows);
            free(entries);
            snprintf(message, sizeof message, "isotropy row out of range: %d", source_row_id);
            return fail(PY_KEY, message), fail_kind;
        }
        int count = isotropy_row_dim(source_row_id);
        int64_t pointer = section_int(&T.iso_op_pointer, idx) - 1;
        const double *constants = (const double *)T.const_values.data;
        R.source_rows = count;
        R.row_stride = full_dim > 0 ? full_dim : 1;
        R.rows = malloc(sizeof(double) * (size_t)(count * R.row_stride + 1));
        R.row_lengths = malloc(sizeof(int) * (size_t)(count + 1));
        if (!R.rows || !R.row_lengths) {
            free(rows);
            free(entries);
            return fail(PY_MEMORY, "mode frame ran out of memory"), fail_kind;
        }
        for (int r = 0; r < count; r++) {
            R.row_lengths[r] = full_dim;
            for (int c = 0; c < full_dim; c++)
                R.rows[r * R.row_stride + c] = constants[section_int(&T.iso_op, pointer + (int64_t)r * full_dim + c)];
        }
    } else if (nparams) {
        free(rows);
        free(entries);
        return fail(PY_RUNTIME, "parametric mode block has no Source order-parameter rows"), fail_kind;
    } else {
        has_source = 0;
    }
    R.source_row_id = source_row_id;
    /* atom bases */
    int64_t baseline_basis[9];
    int has_baseline = has_selected_basis;
    if (has_baseline) {
        memcpy(baseline_basis, selected_basis, sizeof baseline_basis);
        if (has_canonical) {
            int same = 1;
            for (int r = 0; r < 3 && same; r++) {
                int in = 0;
                for (int c = 0; c < 3 && !in; c++) in = !memcmp(selected_basis + 3 * r, canonical_basis + 3 * c, 3 * sizeof(int64_t));
                same = in;
            }
            for (int r = 0; r < 3 && same; r++) {
                int in = 0;
                for (int c = 0; c < 3 && !in; c++) in = !memcmp(canonical_basis + 3 * r, selected_basis + 3 * c, 3 * sizeof(int64_t));
                same = in;
            }
            if (same) memcpy(baseline_basis, canonical_basis, sizeof baseline_basis);
        }
    }
    int has_occurrence = has_occurrence_basis && has_override;
    free(rows);
    /* presentation setting frame and the final float transforms */
    setting_frame frame = {0};
    double transforms[18];
    int transform_count = 0;
    {
        rat cml[3][3], cml_origin[3];
        if (geo_setting_to_cinter(sg, 0, 0, cml, cml_origin)) {
            free(entries);
            return fail(PY_VALUE, "setting Source tables are inconsistent"), fail_kind;
        }
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) transforms[3 * r + c] = (double)cml[r][c].n / (double)cml[r][c].d;
        transform_count = 1;
    }
    int negate = 0;
    if (has_setting) {
        rat m[3][3], o[3];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) m[r][c] = rat_make(setting_words[2 * (3 * r + c)], setting_words[2 * (3 * r + c) + 1]);
        for (int c = 0; c < 3; c++) o[c] = rat_make(setting_origin_words[2 * c], setting_origin_words[2 * c + 1]);
        rat det = rat_sub(rat_add(rat_mul(m[0][0], rat_sub(rat_mul(m[1][1], m[2][2]), rat_mul(m[1][2], m[2][1]))),
                                  rat_mul(m[0][2], rat_sub(rat_mul(m[1][0], m[2][1]), rat_mul(m[1][1], m[2][0])))),
                          rat_mul(m[0][1], rat_sub(rat_mul(m[1][0], m[2][2]), rat_mul(m[1][2], m[2][0]))));
        i128 scale = 1;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) scale = scale / rat_gcd(scale, m[r][c].d) * m[r][c].d;
        for (int c = 0; c < 3; c++) scale = scale / rat_gcd(scale, o[c].d) * o[c].d;
        frame.has = 1;
        frame.scale = scale;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) frame.ints[r][c] = m[r][c].n * (scale / m[r][c].d);
        for (int c = 0; c < 3; c++) frame.origin[c] = o[c].n * (scale / o[c].d);
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) transforms[9 + 3 * r + c] = (double)m[r][c].n / (double)m[r][c].d;
        transform_count = 2;
        negate = vector_setting == 2 && det.n < 0;
    }
    int64_t setting_key[24] = {0};
    if (has_setting) {
        memcpy(setting_key, setting_words, sizeof(int64_t) * 18);
        memcpy(setting_key + 18, setting_origin_words, sizeof(int64_t) * 6);
    }
    if (rat_overflow) {
        free(entries);
        return fail(PY_ARITH, "mode frame setting value out of range"), fail_kind;
    }
    if (!has_occurrence && !has_baseline) {
        R.empty = 1;
        free(entries);
        return 0;
    }
    int eligible = irrep_type == 1 && nparams && vector_setting == 1 && has_override && override_count > 0;
    layout candidate = {0}, baseline = {0};
    int have_candidate = 0, have_baseline = 0, baseline_failed = 0, baseline_kind = 0;
    char baseline_message[200] = "";
    if (has_occurrence && eligible) {
        int raised;
        if (materialize(sg, row_id, site_pg, site_words, nsite, occurrence_basis, &frame, setting_key, &candidate, &raised)) {
            if (!raised || fail_kind == PY_ARITH || fail_kind == PY_MEMORY || fail_kind == PY_RUNTIME) {
                free(entries);
                return fail_kind;
            }
            fail_kind = 0;
        } else {
            have_candidate = 1;
        }
    }
    if (has_baseline) {
        int raised;
        if (materialize(sg, row_id, site_pg, site_words, nsite, baseline_basis, &frame, setting_key, &baseline, &raised)) {
            if (fail_kind == PY_MEMORY) {
                free(entries);
                return fail_kind;
            }
            baseline_failed = 1;
            baseline_kind = fail_kind;
            snprintf(baseline_message, sizeof baseline_message, "%s", fail_message);
            fail_kind = 0;
        } else {
            have_baseline = 1;
        }
    }
    int freeparam = -1;
    if (emit_freeparam && source_row_id > 0) freeparam = isotropy_row_free(source_row_id);
    double complex *inverse = 0;
    int inverse_error = 0;
    if (irrep_type == 2 || irrep_type == 3) {
        inverse = malloc(sizeof(double complex) * (size_t)(full_dim * full_dim + 1));
        if (!inverse) {
            free(entries);
            return fail(PY_MEMORY, "mode frame ran out of memory"), fail_kind;
        }
        int t = cached_transform_inverse(gid, full_dim, inverse);
        if (t < 0) {
            inverse_error = 1;
            free(inverse);
            inverse = 0;
        }
    }
    int64_t params[6] = {0, 1, 0, 1, 0, 1};
    for (int i = 0; i < nparams && i < 3; i++) {
        rat v = rat_make(param_words[2 * i], param_words[2 * i + 1]);
        params[2 * i] = (int64_t)v.n;
        params[2 * i + 1] = (int64_t)v.d;
    }
    double zero_direction[48] = {0};
    int direction_rows = 0;
    double *direction = 0;
    if (have_candidate && override_count > 0) {
        int ok = 1;
        for (int r = 0; r < override_count; r++) ok &= override_lengths[r] == full_dim;
        if (ok) {
            direction = malloc(sizeof(double) * (size_t)(override_count * full_dim + 1));
            if (!direction) {
                free(entries);
                free(inverse);
                return fail(PY_MEMORY, "mode frame ran out of memory"), fail_kind;
            }
            for (int r = 0; r < override_count; r++)
                for (int c = 0; c < full_dim; c++) direction[r * full_dim + c] = override_rows[r * override_stride + c];
            direction_rows = override_count;
        }
    }
    double *row_values = R.rows;
    int row_stride = R.row_stride > 0 ? R.row_stride : 1;
    double one_row[1] = {0};
    int zero_length[1] = {0};
    int block = mode_block(sg, row_id, site_pg, gid, irrep_type, full_dim, old_id, nentries, entries, nparams ? 3 : 0, params,
                           (const double *)inverse, inverse_error, vector_setting, R.source_rows, row_stride,
                           R.source_rows ? R.row_lengths : zero_length, R.source_rows ? row_values : one_row,
                           has_source && R.source_rows > 0, source_row_id, freeparam,
                           (emit_freeparam ? 1 : 0) | (emit_all_rows ? 2 : 0), direction_rows,
                           direction_rows ? direction : zero_direction, direction_rows ? candidate.count : 0,
                           direction_rows ? candidate.records : 0, have_baseline ? baseline.count : -1,
                           have_baseline ? baseline.records : 0, transform_count, transforms, negate);
    free(entries);
    free(inverse);
    free(direction);
    if (block < 0) return -100 + block;
    if (block) {
        if (baseline_failed) return fail(baseline_kind, baseline_message), fail_kind;
        R.empty = 1;
        return 0;
    }
    extern int64_t mode_block_ints(const int64_t **out);
    const int64_t *ints;
    mode_block_ints(&ints);
    R.adopted = (int)ints[0];
    if (R.adopted) {
        R.atoms = candidate;
        R.lanes = direction_rows ? candidate.count : 0;
    } else {
        if (!have_baseline) return fail(PY_TYPE, "cannot unpack non-iterable NoneType object"), fail_kind;
        R.atoms = baseline;
        R.lanes = baseline.count;
    }
    return 0;
}
