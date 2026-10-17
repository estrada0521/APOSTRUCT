/* build_mode_details' site x spec loop and its output: every site block's
 * mode frame, layout, row plan and presented modes, the occurrence-alias
 * admission, the coupled independent union and the public definitions. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "md.h"
#include "tables.h"

int mode_frame(int sg, const char *wyckoff, int64_t wyckoff_length, const char *k_label, int64_t k_label_length, int nparams,
               const int64_t *param_words, int nsite, const int64_t *site_words, int selected_gid, int has_row_id,
               int selected_row_id, int has_basis_override, const int64_t *basis_override, int emit_all_rows, int emit_freeparam,
               int vector_setting, int override_count, int override_stride, const int *override_lengths, const double *override_rows,
               int has_occurrence_basis, const int64_t *occurrence_basis, int has_setting, const int64_t *setting_words,
               const int64_t *setting_origin_words);
void mode_frame_result(int64_t *out);
const int64_t *mode_frame_points(void);
const int *mode_frame_records(void);
const double *mode_frame_rows(void);
const int *mode_frame_row_lengths(void);
const char *mode_frame_error(int *kind);
int64_t mode_block_ints(const int64_t **data);
int64_t mode_block_doubles(const double **data);
int block_layout(int site_pg, int little_type, int count, const int64_t *summaries, int direct, int full_dim, int row_total,
                 const int *row_lengths, int mode_count, const int64_t *block, int64_t *modes, int64_t *keys, int64_t *rows,
                 int *status, int *flags, int64_t *diagnostic);
int block_plan(int sg, int child_sg, const int64_t *basis_pml, int point_operation_count, int atom_count, const int64_t *topology_atoms,
               int presentation_row_count, int layout_point_operation_count, int record_count, const int64_t *records,
               int raw_vector_count, int presented, const double *presentation_basis, int magnetic);
int block_plan_result(int64_t *header, double *inverse, const int64_t **raw_indices, const int64_t **offsets,
                      const int64_t **equivalents, const unsigned char **negated);
int mode_rows_from_plan(int mode_count, int raw_count, const double *raw, int has_pre, const double *pre, int has_inverse,
                        const double *inverse, int atom_count, const int64_t *raw_indices, const int64_t *equivalent_offsets,
                        const int64_t *equivalents, const unsigned char *negated, double *scratch, double *out);
int present_block(int mode_count, int64_t rows, double *dxyz, const int64_t *keys, const int64_t *family_rows, const double *basis,
                  int has_centering, double centering, double *normfactor, int *has_normfactor, int *perm);
int quantize4(int count, const double *values, int64_t *out);
int independent_columns(int row_count, int col_count, const int64_t *offsets, const int *rows, const int64_t *values, int *pivots_out,
                        int *pivot_count);
int spec_label(int sg, const char *k_label, int64_t k_label_length, int flags, int has_gid, int64_t gid, int has_source,
               const int64_t *source_words, int has_carrier, const int64_t *carrier_words, int case_count, const int64_t *case_words,
               int display_count, const int64_t *display_words, int has_selected_gid, int64_t selected_gid, int has_selected_source,
               const int64_t *selected_source_words, int request_count, const int64_t *request_words, const char *fallback,
               int64_t fallback_length);
const char *spec_label_text(void);
int spec_label_params(int64_t *words);
const char *direction_matrix_text(int rows, const int64_t *row_lengths, const double *values, double tol, int group_width);
int little_star_arm_width(int64_t gid);
int emit_all_source_rows(int row_count, int width, const double *rows, int has_free, int free_count, int flags, int64_t old_id,
                         int64_t gid, int little_type, int direction_rows, int direction_width, const double *direction);

/* ------------------------------------------------------------ tables */

static struct {
    section_view irr_type, irr_sg, irr_k, star_count, full_dim, pg_irrep, irrep_label_pg;
    int ready;
} T;

static int tables(void)
{
    if (T.ready) return 0;
    if (section_find("little.little_irr_type", &T.irr_type) || section_find("little.little_irr_space_group", &T.irr_sg) ||
        section_find("little.little_irr_k", &T.irr_k) || section_find("little.little_k_star_count", &T.star_count) ||
        section_find("little.little_irr_full_dim", &T.full_dim) || section_find("wyckoff.iwyckoff_pg_irrep", &T.pg_irrep) ||
        section_find("irreps.irrep_label_pg", &T.irrep_label_pg))
        return 1;
    T.ready = 1;
    return 0;
}

static int little_type_of(int64_t gid) { return (int)section_int(&T.irr_type, gid - 1); }

static int transform_block_count(int64_t gid)
{
    int64_t sg = section_int(&T.irr_sg, gid - 1), kslot = section_int(&T.irr_k, gid - 1);
    return (int)section_int(&T.star_count, (sg - 1) * 27 + kslot - 1);
}

static char *site_pg_irrep_label(int site_pg, int pg_irrep)
{
    int64_t old_id = section_int(&T.pg_irrep, (int64_t)(site_pg - 1) * 12 + pg_irrep - 1);
    if (old_id <= 0) {
        md_fail(MDE_KEY, "'PG%d irrep %d is not defined'", site_pg, pg_irrep);
        return 0;
    }
    const char *s;
    int64_t n;
    if (section_string(&T.irrep_label_pg, old_id - 1, &s, &n)) return (char *)md_strdup("", 0);
    while (n && (*s == ' ' || *s == '\t')) {
        s++;
        n--;
    }
    while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    return (char *)md_strdup(s, n);
}

/* ------------------------------------------------------------ records */

typedef struct {
    int index, row_count;
    md_matrix rows; /* ragged rows padded; lengths kept per row */
    int *lengths;
} md_group;

typedef struct {
    int64_t gid, pg_irrep, family, source_row_count, flat, little_type;
    int ngroups;
    md_group *groups;
} md_summary;

typedef struct {
    int spec_order, site_order, magnetic;
    char *label_prefix, *kvector, *direction;
    const char *k, *irrep;
    int64_t gid;
    const md_site *site;
    int subspace_replaced;
    const md_spec *spec;
    int stokes, dynamic, phase, nkp;
    int64_t kp[4];
    md_summary *summaries;
    int nsummaries;
} md_block;

typedef struct {
    int present;
    int64_t gid, site_pg, pg_irrep, family, component, print_component, family_width, source_row_count, little_type;
    int summary;
} md_identity;

typedef struct {
    int block;
    char *site_irrep, *display;
    md_identity identity;
    double *vectors;
    int plan;
    int has_normfactor;
    double normfactor;
} md_mode;

typedef struct {
    int site;
    int nrecords;
    int64_t *records;
    int raw_count, presented, magnetic;
    int vector_count, has_inverse;
    double inverse[9];
    int natoms;
    int64_t *raw_indices, *offsets, *equivalents;
    unsigned char *negated;
    int has_negated;
    int atom_block;
} md_plan;

typedef struct {
    int site, nrecords;
    int64_t *records;
} md_topology;

typedef struct {
    int spec, site, row_id;
    char *signature;
    int mode_count, all_nonzero;
    int nrecords;
    int64_t *records;
    int nfractionals;
    int64_t *fractionals; /* (num % den, den) x 3 per atom */
} md_emission;

typedef struct {
    int spec, block, first, count;
    int64_t *keys, *rows;
} md_pending;

typedef struct {
    md_block *blocks;
    int nblocks, block_cap;
    md_mode *modes;
    int nmodes, mode_cap;
    md_plan *plans;
    int nplans, plan_cap;
    md_topology *topologies;
    int ntopologies, topology_cap;
    md_emission *emissions;
    int nemissions, emission_cap;
    md_pending *pending;
    int npending, pending_cap;
    int *site_definitions, *site_atoms;
    int64_t *pending_counts; /* spec x site */
    int64_t *eligible;       /* spec * nsites + site, in observation order */
    int neligible, eligible_cap;
    /* spec texts */
    int *text_done;
    int (*text_k_count);
    rat (*text_k)[8];
    char **text_kvector, **text_direction;
    /* selected isotropy for subgroup_details */
    int site_isotropy;
    int64_t iso_row_id, iso_subgroup, iso_arms, iso_direction, iso_basis[9], iso_origin[4];
    int iso_has_origin;
} md_state;

#define GROW(array, count, cap)                                                              \
    do {                                                                                     \
        if ((count) == (cap)) {                                                              \
            int grown_cap_ = (cap) ? 2 * (cap) : 16;                                         \
            void *grown_ = realloc((array), sizeof(*(array)) * (size_t)grown_cap_);          \
            if (!grown_) return md_fail(MDE_MEMORY, "out of memory");                        \
            (array) = grown_;                                                                \
            (cap) = grown_cap_;                                                              \
        }                                                                                    \
    } while (0)

/* ------------------------------------------------------------ occurrence aliases */

/* _exact_source_kparam of an integral record */
static int exact_kp(const int64_t *kp) { return kp[3] != 0; }

typedef struct {
    int *candidate;  /* spec order -> 1 when it carries an anchor */
    int *anchor_of;  /* candidate -> anchor order, -1 none */
    int *observation;
    int *observed;
    int *distinct;
} md_aliases;

static int alias_orders(const md_request *r, md_aliases *a)
{
    int n = r->nspecs;
    a->candidate = calloc((size_t)n + 1, sizeof(int));
    a->anchor_of = malloc(sizeof(int) * ((size_t)n + 1));
    a->observation = calloc((size_t)n + 1, sizeof(int));
    a->observed = calloc((size_t)n + 1, sizeof(int));
    a->distinct = calloc((size_t)n + 1, sizeof(int));
    if (!a->candidate || !a->anchor_of || !a->observation || !a->observed || !a->distinct) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < n; i++) {
        a->candidate[i] = r->specs[i].has_anchor;
        a->anchor_of[i] = -1;
        a->observation[i] = r->specs[i].obs;
    }
    for (int c = 0; c < n; c++) {
        if (!a->candidate[c]) continue;
        const md_spec *cs = r->specs + c;
        if (!exact_kp(cs->anchor_kp)) continue;
        int matches = 0, match = -1;
        for (int i = 0; i < n; i++) {
            const md_spec *s = r->specs + i;
            if (a->candidate[i] || s->magnetic_kind != cs->magnetic_kind) continue;
            int64_t gid = s->gid_state == MD_SET ? s->gid : 0;
            if (gid <= 0 || s->kp_state != MD_SET || !exact_kp(s->kp)) continue;
            if (gid != cs->anchor_gid || memcmp(s->kp, cs->anchor_kp, sizeof s->kp)) continue;
            matches++;
            match = i;
        }
        if (matches == 1) a->anchor_of[c] = match;
    }
    for (int c = 0; c < n; c++) {
        if (a->anchor_of[c] < 0) continue;
        a->observed[c] = a->observed[a->anchor_of[c]] = 1;
        const md_spec *x = r->specs + c, *y = r->specs + a->anchor_of[c];
        int same = x->dir_state == y->dir_state || (x->dir_state != MD_SET && y->dir_state != MD_SET);
        if (same && x->dir_state == MD_SET)
            same = x->dir.rows == y->dir.rows && x->dir.cols == y->dir.cols &&
                   !memcmp(x->dir.values, y->dir.values, sizeof(double) * (size_t)(x->dir.rows * x->dir.cols));
        if (same && x->dir_state == MD_SET)
            for (int k = 0; k < x->dir.rows * x->dir.cols; k++) same &= x->dir.values[k] == y->dir.values[k];
        a->distinct[c] = !same;
    }
    return 0;
}

/* ------------------------------------------------------------ helpers */

static void put_rat_words(int64_t *out, const rat *values, int count)
{
    for (int i = 0; i < count; i++) {
        write_i128(out + 4 * i, values[i].n);
        write_i128(out + 4 * i + 2, values[i].d);
    }
}

/* request_context._spec_label: case k parameters and the displayed k vector */
static int spec_texts(md_request *r, md_state *st, int index)
{
    const md_spec *s = r->specs + index;
    rat request[8];
    int request_count = 0, request_error = 0;
    {
        static const char *keys[6] = {"a", "b", "g", "alpha", "beta", "gamma"};
        const jvalue *params = r->k_params;
        if (jtruthy(params) && jis(params, J_OBJECT))
            for (int i = 0; i < 6 && !request_error; i++) {
                const jvalue *value = jget(params, keys[i]);
                if (!value || jis(value, J_NULL)) continue;
                char *text = md_py_str(value);
                if (!text) return md_fail(MDE_MEMORY, "out of memory");
                const char *a = text, *e = text + strlen(text);
                while (a < e && (*a == ' ' || (*a >= '\t' && *a <= '\r'))) a++;
                while (e > a && (e[-1] == ' ' || (e[-1] >= '\t' && e[-1] <= '\r'))) e--;
                if (a < e) {
                    rat value_rat;
                    int status = py_fraction_of_text(a, e - a, &value_rat);
                    if (status) {
                        jbuf b = {0};
                        jb_reset(&b);
                        py_repr_text(&b, a, e - a);
                        if (status == 2) md_fail(MDE_ARITH, "k parameter out of range");
                        else md_fail(MDE_VALUE, "Invalid literal for Fraction: %s", b.data);
                        jb_free(&b);
                        request_error = 1;
                    } else if (request_count < 8) {
                        request[request_count++] = value_rat;
                    }
                }
                free(text);
            }
    }
    md_error saved = md_err;
    int flags = (s->coupled ? 1 : 0) | (s->coupled ? 2 : 0) | (s->primary ? 4 : 0) | (s->magnetic_state == 1 ? 8 : 0) |
                (s->pml_state == MD_SET ? 16 : 0);
    char *selected_label = md_py_str(jget(r->selected_k, "label"));
    if (!selected_label) return md_fail(MDE_MEMORY, "out of memory");
    if (!strcmp(s->k_label, selected_label)) flags |= 32;
    free(selected_label);
    int64_t case_words[32], display_words[32], request_words[32], kp[4] = {0, 0, 0, 1}, carrier[4] = {0, 0, 0, 1};
    int case_count = -1, display_count = -1;
    if (s->case_state == MD_SET) {
        case_count = s->case_n;
        put_rat_words(case_words, s->case_k, s->case_n);
    }
    if (s->coupled && s->dcase_state == MD_SET) {
        display_count = s->dcase_n;
        put_rat_words(display_words, s->dcase_k, s->dcase_n);
    }
    put_rat_words(request_words, request, request_count);
    int has_source = s->kp_state == MD_SET, has_carrier = s->carrier_state == MD_SET;
    if (has_source) memcpy(kp, s->kp, sizeof kp);
    if (has_carrier) memcpy(carrier, s->carrier, sizeof carrier);
    int has_selected_source = 0;
    int64_t selected_source[4] = {0, 0, 0, 1};
    const jvalue *raw = r->iso ? jget(r->iso, "source_kparam") : 0;
    if (raw && !jis(raw, J_NULL) && jis(raw, J_ARRAY) && raw->count == 4) {
        has_selected_source = 1;
        for (int i = 0; i < 4; i++)
            if (py_int_of(jat(raw, i), &selected_source[i])) has_selected_source = 0;
    }
    int status = spec_label(r->sg, s->k_label, (int64_t)strlen(s->k_label), flags, s->gid_state == MD_SET, s->gid, has_source, kp,
                            has_carrier, carrier, case_count, case_words, display_count, display_words, r->has_selected_gid,
                            r->selected_gid, has_selected_source, selected_source, request_error ? -1 : request_count, request_words,
                            r->display_kvector, (int64_t)strlen(r->display_kvector));
    if (status == 6 && request_error) {
        md_err = saved;
        return md_err.kind;
    }
    if (status) {
        static const int kinds[] = {0, MDE_VALUE, MDE_ARITH, MDE_ARITH, MDE_INDEX, MDE_KEY, MDE_RUNTIME, MDE_MEMORY, MDE_LINALG};
        static const char *messages[] = {"", "mode label Source tables are inconsistent", "k-vector value out of range",
                                         "K-family display expression is not affine rational", "index out of range", "'direction_matrix'",
                                         "mode label failed (6)", "mode label text ran out of memory", "SVD did not converge"};
        return md_fail(status < 9 ? kinds[status] : MDE_RUNTIME, "%s", status < 9 ? messages[status] : "mode label failed");
    }
    int64_t words[32];
    int count = spec_label_params(words);
    st->text_k_count[index] = count;
    for (int i = 0; i < count && i < 8; i++) st->text_k[index][i] = (rat){read_i128(words + 4 * i), read_i128(words + 4 * i + 2)};
    const char *text = spec_label_text();
    st->text_kvector[index] = (char *)md_strdup(text, (int64_t)strlen(text));
    /* _spec_opd_direction_text */
    if (s->dir_state == MD_SET && s->dir.rows > 0 && (!s->primary || s->coupled)) {
        int width = 0;
        if (s->gid_state == MD_SET) {
            width = little_star_arm_width(s->gid);
            if (width < 0) return md_fail(width == -4 ? MDE_INDEX : MDE_VALUE, width == -4 ? "index out of range" : "mode label Source tables are inconsistent");
        }
        int64_t *lengths = malloc(sizeof(int64_t) * (size_t)s->dir.rows);
        if (!lengths) return md_fail(MDE_MEMORY, "out of memory");
        for (int i = 0; i < s->dir.rows; i++) lengths[i] = s->dir.cols;
        const char *direction = direction_matrix_text(s->dir.rows, lengths, s->dir.values, 1e-7, width);
        free(lengths);
        if (!direction) return md_fail(MDE_MEMORY, "direction matrix text ran out of memory");
        st->text_direction[index] = (char *)md_strdup(direction, (int64_t)strlen(direction));
    } else {
        st->text_direction[index] = (char *)md_strdup(r->opd_direction, (int64_t)strlen(r->opd_direction));
    }
    if (!st->text_kvector[index] || !st->text_direction[index]) return md_fail(MDE_MEMORY, "out of memory");
    st->text_done[index] = 1;
    return 0;
}

/* request_context._should_emit_all_source_opd_rows */
static int emit_full_rows(const md_request *r, const md_spec *s, int *out)
{
    *out = 0;
    if (s->rows_state != MD_SET || s->rows.rows == 0) return 0;
    int direction_rows = -1, direction_width = 0;
    const double *direction = 0;
    double zero = 0.0;
    if (s->dir_state == MD_NONE || (s->dir_state == MD_SET && s->dir.rows == 0)) {
        direction_rows = -2;
    } else if (s->dir_state == MD_SET) {
        direction_rows = s->dir.rows;
        direction_width = s->dir.cols;
        direction = s->dir.values;
    }
    int flags = (r->coupled ? 1 : 0) | (s->primary ? 4 : 0);
    int result = emit_all_source_rows(s->rows.rows, s->rows.cols, s->rows.values, s->free_state == MD_SET, (int)(s->free_state == MD_SET ? s->free_count : 0),
                                      flags, s->old_id, s->gid_state == MD_SET ? s->gid : 0, (int)(s->little_type_state == MD_SET ? s->little_type : 0),
                                      direction_rows, direction_width, direction ? direction : &zero);
    if (result < 0) {
        int code = -result;
        if (code == 5) return md_fail(MDE_KEY, "'direction_matrix'");
        if (code == 4) return md_fail(MDE_INDEX, "index out of range");
        if (code == 8) return md_fail(MDE_LINALG, "SVD did not converge");
        return md_fail(MDE_VALUE, "mode label Source tables are inconsistent");
    }
    *out = result;
    return 0;
}

/* display_distortion.occurrence_sublattice(...) is not None */
static int occurrence_sublattice(const int64_t *source, const int64_t *occurrence)
{
    int64_t ds = source[0] * (source[4] * source[8] - source[5] * source[7]) - source[1] * (source[3] * source[8] - source[5] * source[6]) +
                 source[2] * (source[3] * source[7] - source[4] * source[6]);
    int64_t dc = occurrence[0] * (occurrence[4] * occurrence[8] - occurrence[5] * occurrence[7]) -
                 occurrence[1] * (occurrence[3] * occurrence[8] - occurrence[5] * occurrence[6]) +
                 occurrence[2] * (occurrence[3] * occurrence[7] - occurrence[4] * occurrence[6]);
    if (!ds || !dc) return 0;
    rat src[3][3], inverse[3][3], occ[3], out[3];
    for (int i = 0; i < 9; i++) src[i / 3][i % 3] = rat_int(source[i]);
    if (rat_inverse3((const rat (*)[3])src, inverse)) return 0;
    int64_t transition[9];
    for (int row = 0; row < 3; row++) {
        for (int k = 0; k < 3; k++) occ[k] = rat_int(occurrence[3 * row + k]);
        rat_row_mul(occ, (const rat (*)[3])inverse, out);
        for (int col = 0; col < 3; col++) {
            if (out[col].d != 1) return 0;
            transition[3 * row + col] = (int64_t)out[col].n;
        }
    }
    int64_t dt = transition[0] * (transition[4] * transition[8] - transition[5] * transition[7]) -
                 transition[1] * (transition[3] * transition[8] - transition[5] * transition[6]) +
                 transition[2] * (transition[3] * transition[7] - transition[4] * transition[6]);
    int64_t index = dt < 0 ? -dt : dt;
    int64_t ads = ds < 0 ? -ds : ds, adc = dc < 0 ? -dc : dc;
    return index > 1 && adc == ads * index;
}

static void put_rat_pairs(int64_t *out, const rat *values, int count)
{
    for (int i = 0; i < count; i++) {
        out[2 * i] = (int64_t)values[i].n;
        out[2 * i + 1] = (int64_t)values[i].d;
    }
}

/* canonical text of a frozen Source signature (equal iff Python's tuples are) */
static void signature_float(jbuf *b, double value)
{
    if (value == 0.0) value = 0.0;
    jb_float(b, value);
}

/* ------------------------------------------------------------ loop */

typedef struct {
    int empty, adopted, lanes, site_pg, atom_count, source_count, source_stride;
    int64_t gid;
    int has_selected;
    int64_t row_id, subgroup, arms, direction, basis[9], origin[4];
    int has_origin;
    int64_t *points;  /* atom_count x 12 words */
    int *records;     /* atom_count x 5 */
    double *source_rows;
    int *source_lengths;
    int nsummaries;
    md_summary *summaries;
    int64_t *summary_words;
    int64_t summary_word_count;
    int mode_count;
    double *modes; /* mode_count x lanes x 3 */
} md_frame;

static void summaries_free(md_summary *summaries, int count)
{
    for (int k = 0; k < count; k++) {
        for (int g = 0; g < summaries[k].ngroups; g++) {
            free(summaries[k].groups[g].rows.values);
            free(summaries[k].groups[g].lengths);
        }
        free(summaries[k].groups);
    }
    free(summaries);
}

static void frame_free(md_frame *f)
{
    summaries_free(f->summaries, f->nsummaries);
    free(f->points);
    free(f->records);
    free(f->source_rows);
    free(f->source_lengths);
    free(f->summary_words);
    free(f->modes);
    memset(f, 0, sizeof *f);
}

static int frame_error(int status)
{
    if (status < -100) {
        int code = -100 - status;
        md_err.kind = code > 200 ? 110 : code > 100 ? 111 : 112;
        md_err.status = code > 200 ? code - 200 : code > 100 ? code - 100 : code;
        md_err.message[0] = 0;
        return md_err.kind;
    }
    int kind = 0;
    const char *message = mode_frame_error(&kind);
    return md_fail(kind ? kind : MDE_RUNTIME, "%s", message);
}

static int run_frame(md_request *r, const md_site *site, const md_spec *spec, const rat *k_params, int k_count, int64_t trace_gid,
                     int raw_basis, int emit_all, int emit_freeparam, const md_matrix *override, int has_occurrence_basis, md_frame *f)
{
    memset(f, 0, sizeof *f);
    char *wyckoff = md_py_str(jget(site->site, "wyckoff"));
    if (!wyckoff) return md_fail(MDE_MEMORY, "out of memory");
    int64_t param_words[16], site_words[12];
    put_rat_pairs(param_words, k_params, k_count);
    rat site_params[3];
    for (int i = 0; i < site->site_param_count; i++)
        if (py_fraction_of_float(site->site_params[i], &site_params[i])) {
            free(wyckoff);
            return md_fail(MDE_ARITH, "site parameter out of range");
        }
    put_rat_words(site_words, site_params, site->site_param_count > 0 ? site->site_param_count : 0);
    int override_count = override ? override->rows : -1, stride = 1;
    if (override && override->cols > stride) stride = override->cols;
    int *lengths = malloc(sizeof(int) * (size_t)(override && override->rows ? override->rows : 1));
    double *rows = calloc((size_t)((override && override->rows ? override->rows : 1) * stride), sizeof(double));
    if (!lengths || !rows) {
        free(wyckoff);
        free(lengths);
        free(rows);
        return md_fail(MDE_MEMORY, "out of memory");
    }
    if (override)
        for (int i = 0; i < override->rows; i++) {
            lengths[i] = override->cols;
            memcpy(rows + (size_t)i * stride, override->values + (size_t)i * override->cols, sizeof(double) * (size_t)override->cols);
        }
    int64_t setting_words[18], origin_words[6];
    for (int i = 0; i < 9; i++) {
        setting_words[2 * i] = r->has_bridge ? (int64_t)r->bridge_matrix[i / 3][i % 3].n : 0;
        setting_words[2 * i + 1] = r->has_bridge ? (int64_t)r->bridge_matrix[i / 3][i % 3].d : 1;
    }
    for (int i = 0; i < 3; i++) {
        origin_words[2 * i] = r->has_bridge ? (int64_t)r->bridge_shift[i].n : 0;
        origin_words[2 * i + 1] = r->has_bridge ? (int64_t)r->bridge_shift[i].d : 1;
    }
    int64_t zero9[9] = {0};
    int status = mode_frame(r->sg, wyckoff, (int64_t)strlen(wyckoff), spec->k_label, (int64_t)strlen(spec->k_label), k_count, param_words,
                            site->site_param_count > 0 ? site->site_param_count : 0, site_words, (int)trace_gid,
                            spec->row_id_state == MD_SET, (int)(spec->row_id_state == MD_SET ? spec->row_id : 0), 1,
                            raw_basis ? r->source_basis_override : r->presentation_pml, emit_all, emit_freeparam, spec->vector_setting,
                            override_count, stride, lengths, rows, has_occurrence_basis, has_occurrence_basis ? r->presentation_pml : zero9,
                            r->has_bridge, setting_words, origin_words);
    free(wyckoff);
    free(lengths);
    free(rows);
    if (status) return frame_error(status);
    int64_t header[29];
    mode_frame_result(header);
    f->site_pg = (int)header[0];
    f->gid = header[1];
    f->has_selected = (int)header[2];
    f->row_id = header[3];
    f->subgroup = header[4];
    f->arms = header[5];
    f->direction = header[6];
    f->has_origin = (int)header[7];
    memcpy(f->origin, header + 8, sizeof f->origin);
    memcpy(f->basis, header + 12, sizeof f->basis);
    f->empty = (int)header[21];
    f->adopted = (int)header[22];
    f->lanes = (int)header[23];
    f->atom_count = (int)header[24];
    f->source_count = (int)header[25];
    f->source_stride = (int)header[26];
    if (f->empty) return 0;
    f->points = malloc(sizeof(int64_t) * 12 * (size_t)(f->atom_count + 1));
    f->records = malloc(sizeof(int) * 5 * (size_t)(f->atom_count + 1));
    f->source_rows = malloc(sizeof(double) * (size_t)(f->source_count * f->source_stride + 1));
    f->source_lengths = malloc(sizeof(int) * (size_t)(f->source_count + 1));
    if (!f->points || !f->records || !f->source_rows || !f->source_lengths) return md_fail(MDE_MEMORY, "out of memory");
    if (f->atom_count) {
        memcpy(f->points, mode_frame_points(), sizeof(int64_t) * 12 * (size_t)f->atom_count);
        memcpy(f->records, mode_frame_records(), sizeof(int) * 5 * (size_t)f->atom_count);
    }
    if (f->source_count) {
        memcpy(f->source_rows, mode_frame_rows(), sizeof(double) * (size_t)(f->source_count * f->source_stride));
        memcpy(f->source_lengths, mode_frame_row_lengths(), sizeof(int) * (size_t)f->source_count);
    }
    const int64_t *ints;
    const double *doubles;
    int64_t int_count = mode_block_ints(&ints), double_count = mode_block_doubles(&doubles);
    int64_t at = 2;
    f->nsummaries = (int)ints[1];
    f->summaries = calloc((size_t)f->nsummaries + 1, sizeof(md_summary));
    if (!f->summaries) return md_fail(MDE_MEMORY, "out of memory");
    for (int k = 0; k < f->nsummaries; k++) {
        md_summary *s = f->summaries + k;
        s->gid = f->gid;
        s->pg_irrep = ints[at];
        s->family = ints[at + 1];
        s->source_row_count = ints[at + 2];
        int group_count = (int)ints[at + 3];
        at += 4;
        s->ngroups = group_count;
        s->groups = calloc((size_t)group_count + 1, sizeof(md_group));
        if (!s->groups) return md_fail(MDE_MEMORY, "out of memory");
        const int64_t *groups = ints + at;
        at += 3 * group_count;
        at += group_count;
        s->flat = ints[at];
        at += 5;
        s->little_type = little_type_of(f->gid);
        if (f->adopted) {
            /* direct groups: one per override row */
            int rows_n = override ? override->rows : 0;
            s->ngroups = rows_n;
            free(s->groups);
            s->groups = calloc((size_t)rows_n + 1, sizeof(md_group));
            if (!s->groups) return md_fail(MDE_MEMORY, "out of memory");
            for (int g = 0; g < rows_n; g++) {
                md_group *grp = s->groups + g;
                grp->index = g;
                grp->row_count = 1;
                grp->rows.rows = 1;
                grp->rows.cols = override->cols;
                grp->rows.values = malloc(sizeof(double) * (size_t)(override->cols + 1));
                grp->lengths = malloc(sizeof(int));
                if (!grp->rows.values || !grp->lengths) return md_fail(MDE_MEMORY, "out of memory");
                memcpy(grp->rows.values, override->values + (size_t)g * override->cols, sizeof(double) * (size_t)override->cols);
                grp->lengths[0] = override->cols;
            }
        } else {
            for (int g = 0; g < group_count; g++) {
                md_group *grp = s->groups + g;
                int64_t start = groups[3 * g], length = groups[3 * g + 1];
                grp->index = g;
                grp->row_count = (int)groups[3 * g + 2];
                grp->rows.rows = (int)length;
                grp->rows.cols = f->source_stride;
                grp->rows.values = malloc(sizeof(double) * (size_t)(length * f->source_stride + 1));
                grp->lengths = malloc(sizeof(int) * (size_t)(length + 1));
                if (!grp->rows.values || !grp->lengths) return md_fail(MDE_MEMORY, "out of memory");
                for (int64_t q = 0; q < length; q++) {
                    int64_t src = start + q;
                    if (src < 0 || src >= f->source_count) {
                        grp->rows.rows = (int)q;
                        break;
                    }
                    grp->lengths[q] = f->source_lengths[src];
                    memcpy(grp->rows.values + q * f->source_stride, f->source_rows + src * f->source_stride, sizeof(double) * (size_t)f->source_stride);
                }
            }
        }
    }
    f->mode_count = (int)ints[at];
    f->summary_word_count = at - 2;
    f->summary_words = malloc(sizeof(int64_t) * (size_t)(f->summary_word_count + 1));
    f->modes = malloc(sizeof(double) * (size_t)(double_count + 1));
    if (!f->summary_words || !f->modes) return md_fail(MDE_MEMORY, "out of memory");
    if (f->summary_word_count) memcpy(f->summary_words, ints + 2, sizeof(int64_t) * (size_t)f->summary_word_count);
    else f->summary_words[0] = 0;
    if (double_count) memcpy(f->modes, doubles, sizeof(double) * (size_t)double_count);
    (void)int_count;
    if (f->lanes < 1) f->lanes = 1;
    return 0;
}

static int modes_nonzero(const md_frame *f)
{
    if (!f->mode_count) return 0;
    for (int m = 0; m < f->mode_count; m++) {
        const double *v = f->modes + (size_t)m * f->lanes * 3;
        int n = f->lanes * 3;
        if (!n) return 0;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            if (!isfinite(v[i])) return 0;
            sum += v[i] * v[i];
        }
        if (!(sqrt(sum) > 1e-12)) return 0;
    }
    return 1;
}

static void summary_group_json(jbuf *b, const md_summary *s)
{
    jb_lit(b, "[");
    for (int g = 0; g < s->ngroups; g++) {
        const md_group *grp = s->groups + g;
        if (g) jb_lit(b, ", ");
        jb_lit(b, "{\"index\": ");
        jb_int(b, grp->index);
        jb_lit(b, ", \"row_count\": ");
        jb_int(b, grp->row_count);
        jb_lit(b, ", \"rows\": [");
        for (int q = 0; q < grp->rows.rows; q++) {
            if (q) jb_lit(b, ", ");
            jb_lit(b, "[");
            for (int c = 0; c < grp->lengths[q]; c++) {
                if (c) jb_lit(b, ", ");
                jb_float(b, grp->rows.values[q * grp->rows.cols + c]);
            }
            jb_lit(b, "]");
        }
        jb_lit(b, "]}");
    }
    jb_lit(b, "]");
}

/* assembly._occurrence_source_signature as canonical text, NULL for None */
static char *source_signature(const md_frame *f, const int64_t *layout_modes, int identity_covered, int site_pg)
{
    if (!identity_covered || f->mode_count <= 0) return 0;
    jbuf b = {0};
    jb_reset(&b);
    for (int m = 0; m < f->mode_count; m++) {
        const int64_t *row = layout_modes + 7 * m;
        const md_summary *s = f->summaries + row[0];
        if (!s->ngroups) {
            jb_free(&b);
            return 0;
        }
        int64_t values[9] = {s->gid, site_pg, s->pg_irrep, s->family, row[6], row[5], s->flat, s->source_row_count, s->little_type};
        if (values[0] <= 0 || values[1] <= 0 || values[2] <= 0 || values[6] <= 0 || values[7] <= 0 || values[8] <= 0 || values[3] < 0 ||
            values[4] < 0 || values[5] < 0) {
            jb_free(&b);
            return 0;
        }
        jb_lit(&b, "(");
        for (int k = 0; k < 9; k++) {
            jb_int(&b, values[k]);
            jb_lit(&b, ",");
        }
        for (int g = 0; g < s->ngroups; g++) {
            const md_group *grp = s->groups + g;
            jb_lit(&b, "{index:");
            jb_int(&b, grp->index);
            jb_lit(&b, ",row_count:");
            jb_int(&b, grp->row_count);
            jb_lit(&b, ",rows:");
            for (int q = 0; q < grp->rows.rows; q++) {
                jb_lit(&b, "[");
                for (int c = 0; c < grp->lengths[q]; c++) {
                    double v = grp->rows.values[q * grp->rows.cols + c];
                    if (!isfinite(v)) {
                        jb_free(&b);
                        return 0;
                    }
                    signature_float(&b, v);
                    jb_lit(&b, ",");
                }
                jb_lit(&b, "]");
            }
            jb_lit(&b, "}");
        }
        jb_lit(&b, ")");
    }
    if (b.failed) {
        jb_free(&b);
        return 0;
    }
    return b.data;
}

static int label_error(md_request *r, int status, const int64_t *diagnostic, int site_pg, const md_frame *f, int mode_count)
{
    (void)r;
    if (status == 1)
        return md_fail(MDE_VALUE, "cannot resolve site-irrep labels from the mode block: mode_count=%d, summaries=%d, site_pg=%d", mode_count,
                       f->nsummaries, site_pg);
    if (status == 2) return md_fail(MDE_VALUE, "site-irrep label layout has a nonpositive dimension or mode count");
    if (status == 3) {
        char *base = site_pg_irrep_label(site_pg, (int)f->summaries[diagnostic[2]].pg_irrep);
        if (!base) return md_err.kind;
        jbuf b = {0};
        jb_reset(&b);
        py_repr_text(&b, base, (int64_t)strlen(base));
        md_fail(MDE_VALUE, "site-irrep label dimension does not divide the emitted mode count: count=%lld, dimension=%lld, base=%s",
                (long long)diagnostic[0], (long long)diagnostic[1], b.data);
        jb_free(&b);
        free(base);
        return MDE_VALUE;
    }
    return md_fail(MDE_VALUE, "site-irrep label count does not match the emitted mode count: labels=%lld, modes=%lld", (long long)diagnostic[0],
                   (long long)diagnostic[1]);
}

static int block_plan_error(int status)
{
    md_err.kind = 113;
    md_err.status = status;
    md_err.message[0] = 0;
    return md_err.kind;
}

static char *component_text(int64_t component)
{
    return component < 26 ? md_format("%c", (char)('a' + component)) : md_format("%lld", (long long)(component + 1));
}

static int topology_for(md_request *r, md_state *st, int site_index, const md_frame *f)
{
    const md_site *site = r->sites + site_index;
    for (int t = 0; t < st->ntopologies; t++) {
        const md_topology *top = st->topologies + t;
        if (top->site == site_index && top->nrecords == f->atom_count) {
            int same = 1;
            for (int i = 0; i < 5 * f->atom_count && same; i++) same = top->records[i] == f->records[i];
            if (same) return 0;
        }
    }
    int64_t *records = malloc(sizeof(int64_t) * 5 * (size_t)(f->atom_count + 1));
    if (!records) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < 5 * f->atom_count; i++) records[i] = f->records[i];
    int status = block_plan(0, r->child_sg, r->presentation_pml, r->point_operation_count, site->natoms, site->topology, site->presentation_rows,
                            r->point_operation_count, f->atom_count, records, -1, 0, 0, 0);
    if (status) {
        free(records);
        return block_plan_error(status);
    }
    GROW(st->topologies, st->ntopologies, st->topology_cap);
    st->topologies[st->ntopologies++] = (md_topology){site_index, f->atom_count, records};
    return 0;
}

static int plan_for(md_request *r, md_state *st, int site_index, const md_frame *f, int raw_count, int presented, int magnetic, int *out)
{
    const md_site *site = r->sites + site_index;
    for (int p = 0; p < st->nplans; p++) {
        const md_plan *plan = st->plans + p;
        if (plan->site != site_index || plan->nrecords != f->atom_count || plan->raw_count != raw_count || plan->presented != presented ||
            plan->magnetic != magnetic)
            continue;
        int same = 1;
        for (int i = 0; i < 5 * f->atom_count && same; i++) same = plan->records[i] == f->records[i];
        if (same) {
            *out = p;
            return 0;
        }
    }
    int64_t *records = malloc(sizeof(int64_t) * 5 * (size_t)(f->atom_count + 1));
    if (!records) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < 5 * f->atom_count; i++) records[i] = f->records[i];
    double presentation[9];
    for (int i = 0; i < 9; i++) presentation[i] = r->display_basis[i / 3][i % 3];
    int status = block_plan(r->sg, r->child_sg, r->presentation_pml, r->point_operation_count, site->natoms, site->topology,
                            site->presentation_rows, r->point_operation_count, f->atom_count, records, raw_count, presented, presentation,
                            magnetic);
    if (status) {
        free(records);
        return block_plan_error(status);
    }
    int64_t header[3];
    double inverse[9];
    const int64_t *raw_indices, *offsets, *equivalents;
    const unsigned char *negated;
    int count = block_plan_result(header, inverse, &raw_indices, &offsets, &equivalents, &negated);
    GROW(st->plans, st->nplans, st->plan_cap);
    md_plan *plan = st->plans + st->nplans;
    memset(plan, 0, sizeof *plan);
    plan->site = site_index;
    plan->nrecords = f->atom_count;
    plan->records = records;
    plan->raw_count = raw_count;
    plan->presented = presented;
    plan->magnetic = magnetic;
    plan->vector_count = (int)header[0];
    plan->has_inverse = (int)header[1];
    memcpy(plan->inverse, inverse, sizeof inverse);
    plan->natoms = count;
    plan->raw_indices = malloc(sizeof(int64_t) * (size_t)(count + 1));
    plan->offsets = malloc(sizeof(int64_t) * (size_t)(count + 1));
    plan->equivalents = malloc(sizeof(int64_t) * (size_t)(header[2] + 1));
    plan->negated = malloc((size_t)count + 1);
    plan->atom_block = -1;
    if (!plan->raw_indices || !plan->offsets || !plan->equivalents || !plan->negated) return md_fail(MDE_MEMORY, "out of memory");
    if (count) {
        memcpy(plan->raw_indices, raw_indices, sizeof(int64_t) * (size_t)count);
        memcpy(plan->offsets, offsets, sizeof(int64_t) * (size_t)(count + 1));
    } else {
        plan->offsets[0] = 0;
    }
    if (header[2]) memcpy(plan->equivalents, equivalents, sizeof(int64_t) * (size_t)header[2]);
    plan->has_negated = magnetic && count;
    if (plan->has_negated) memcpy(plan->negated, negated, (size_t)count);
    *out = st->nplans++;
    return 0;
}

static int md_present(md_request *r, md_state *st, int first, int count, const int64_t *keys, const int64_t *family_rows);

static int site_loop(md_request *r, md_state *st, md_aliases *al)
{
    const jvalue *allowed_sites[2] = {jget(r->request, "displacive_site_indexes"), jget(r->request, "magnetic_site_indexes")};
    const jvalue *allowed_rows[2] = {jget(r->request, "displacive_row_ids"), jget(r->request, "magnetic_row_ids")};
    int raw_override_proof = r->has_source_basis_override && occurrence_sublattice(r->source_basis_override, r->presentation_pml);
    for (int site_index = 0; site_index < r->nsites; site_index++) {
        md_site *site = r->sites + site_index;
        int site_definitions = 0, site_atoms = 0;
        for (int spec_index = 0; spec_index < r->nspecs; spec_index++) {
            md_spec *spec = r->specs + spec_index;
            int mag = spec->magnetic_kind;
            const jvalue *allowed = allowed_sites[mag];
            if (allowed && !jis(allowed, J_NULL)) {
                int found = 0;
                for (int64_t i = 0; i < allowed->count; i++) {
                    int64_t value;
                    if (!py_int_of(jat(allowed, i), &value) && value == site_index + 1) found = 1;
                }
                if (!found) continue;
            }
            int64_t site_row_id = 0;
            const jvalue *row_value = jget(site->site, "wyckoff_row_id");
            if (jtruthy(row_value) && py_int_of(row_value, &site_row_id)) return md_fail(MDE_VALUE, "invalid literal for int()");
            allowed = allowed_rows[mag];
            if (allowed && !jis(allowed, J_NULL)) {
                int found = 0;
                for (int64_t i = 0; i < allowed->count; i++) {
                    int64_t value;
                    if (!py_int_of(jat(allowed, i), &value) && value == site_row_id) found = 1;
                }
                if (!found) continue;
            }
            int candidate = al->candidate[spec_index], observation_only = al->observation[spec_index];
            if (al->observed[spec_index]) {
                GROW(st->eligible, st->neligible, st->eligible_cap);
                st->eligible[st->neligible++] = (int64_t)spec_index * r->nsites + site_index;
            }
            if (!st->text_done[spec_index]) MD_TRY(spec_texts(r, st, spec_index));
            const rat *k_params = st->text_k[spec_index];
            int k_count = st->text_k_count[spec_index];
            int raw_basis = r->has_source_basis_override && (mag || (spec->rows_state == MD_SET && spec->rows.rows > 0));
            int emit_full;
            MD_TRY(emit_full_rows(r, spec, &emit_full));
            md_matrix override = {0};
            int has_override = 0;
            if (spec->rows_state == MD_SET && spec->rows.rows > 0) {
                override = spec->rows;
                has_override = 1;
            }
            md_matrix iso_rows = {0};
            if (!has_override && mag && r->iso) {
                const jvalue *raw = jget(r->iso, "source_numeric_rows");
                if (jtruthy(raw)) {
                    int ndim;
                    MD_TRY(md_float_rows(raw, &iso_rows, &ndim));
                    if (iso_rows.rows) {
                        override = iso_rows;
                        has_override = 1;
                    }
                }
            }
            int64_t trace_gid = spec->gid_state == MD_SET && spec->gid ? spec->gid : 0;
            if (!trace_gid) {
                const jvalue *gid = jget(r->irrep, "gid");
                if (jtruthy(gid) && py_int_of(gid, &trace_gid)) return md_fail(MDE_VALUE, "invalid literal for int()");
            }
            if (trace_gid <= 0) {
                free(iso_rows.values);
                return md_fail(MDE_VALUE, "mode spec has no Source gid");
            }
            int occurrence_basis = raw_basis && !mag && k_count && has_override && raw_override_proof;
            md_frame f;
            int status = run_frame(r, site, spec, k_params, k_count, trace_gid, raw_basis, r->full_parameter_opd || emit_full,
                                   !r->full_parameter_opd && !emit_full, has_override ? &override : 0, occurrence_basis, &f);
            if (status) {
                free(iso_rows.values);
                frame_free(&f);
                return status;
            }
            free(iso_rows.values);
            if (spec->primary && f.has_selected) {
                st->site_isotropy = 1;
                st->iso_row_id = f.row_id;
                st->iso_subgroup = f.subgroup;
                st->iso_arms = f.arms;
                st->iso_direction = f.direction;
                memcpy(st->iso_basis, f.basis, sizeof f.basis);
                memcpy(st->iso_origin, f.origin, sizeof f.origin);
                st->iso_has_origin = f.has_origin;
            }
            int64_t spec_gid = spec->gid_state == MD_SET ? spec->gid : 0;
            int dynamic = spec_gid > 0 && spec->old_id <= 0;
            int phase = !observation_only && f.mode_count && dynamic && little_type_of(spec_gid) == 1 &&
                        (mag || transform_block_count(spec_gid) > 1);
            int mode_count = f.empty ? 0 : f.mode_count;
            int count = mode_count > 0 ? mode_count : 1;
            int64_t *layout_modes = calloc((size_t)count * 7, sizeof(int64_t)), *keys = calloc((size_t)count * 3, sizeof(int64_t));
            int64_t *family_rows = calloc((size_t)count * 19, sizeof(int64_t)), diagnostic[6] = {0};
            int label_status[2] = {0, 0}, flags[2] = {0, 0};
            int *row_lengths = malloc(sizeof(int) * (size_t)(f.source_count + 1));
            if (!layout_modes || !keys || !family_rows || !row_lengths) {
                frame_free(&f);
                return md_fail(MDE_MEMORY, "out of memory");
            }
            for (int i = 0; i < f.source_count; i++) row_lengths[i] = f.source_lengths[i];
            int64_t block_values[6] = {mag ? 2 : 1, dynamic, phase, spec->stokes, spec_index, site_index};
            int64_t zero = 0;
            int code = block_layout(f.site_pg, f.empty ? 0 : (int)little_type_of(f.gid), f.empty ? 0 : f.nsummaries,
                                    f.empty ? &zero : f.summary_words, f.empty ? 0 : f.adopted,
                                    f.empty ? 0 : (int)section_int(&T.full_dim, f.gid - 1), f.empty ? 0 : f.source_count, row_lengths,
                                    mode_count, block_values, layout_modes, keys, family_rows, label_status, flags, diagnostic);
            free(row_lengths);
            if (code) {
                frame_free(&f);
                return code == 5 ? md_fail(MDE_MEMORY, "block layout ran out of memory")
                                 : md_fail(MDE_VALUE, "block layout Source tables are inconsistent");
            }
            int identity_covered = flags[0], display_wanted = flags[1];
            if (al->observed[spec_index]) {
                char *signature = source_signature(&f, layout_modes, identity_covered, f.site_pg);
                for (int e = 0; e < st->nemissions; e++)
                    if (st->emissions[e].spec == spec_index && st->emissions[e].site == site_index) st->emissions[e].spec = -1;
                if (signature) {
                    GROW(st->emissions, st->nemissions, st->emission_cap);
                    md_emission *em = st->emissions + st->nemissions++;
                    em->spec = spec_index;
                    em->site = site_index;
                    em->row_id = (int)site_row_id;
                    em->signature = signature;
                    em->mode_count = mode_count;
                    em->all_nonzero = modes_nonzero(&f);
                    em->nrecords = f.atom_count;
                    em->records = malloc(sizeof(int64_t) * 5 * (size_t)(f.atom_count + 1));
                    em->nfractionals = f.atom_count;
                    em->fractionals = malloc(sizeof(int64_t) * 6 * (size_t)(f.atom_count + 1));
                    if (!em->records || !em->fractionals) return md_fail(MDE_MEMORY, "out of memory");
                    for (int i = 0; i < 5 * f.atom_count; i++) em->records[i] = f.records[i];
                    for (int a = 0; a < f.atom_count; a++)
                        for (int axis = 0; axis < 3; axis++) {
                            i128 n = read_i128(f.points + 12 * a + 4 * axis), d = read_i128(f.points + 12 * a + 4 * axis + 2);
                            rat v = rat_make(n, d), m = rat_mod1(v);
                            em->fractionals[6 * a + 2 * axis] = (int64_t)m.n;
                            em->fractionals[6 * a + 2 * axis + 1] = (int64_t)m.d;
                        }
                }
            }
            if (observation_only) {
                free(layout_modes);
                free(keys);
                free(family_rows);
                frame_free(&f);
                continue;
            }
            if (label_status[0] || (display_wanted && label_status[1])) {
                int s = label_status[0] ? label_status[0] : label_status[1];
                int status_code = label_error(r, s, label_status[0] ? diagnostic : diagnostic + 3, f.site_pg, &f, mode_count);
                frame_free(&f);
                return status_code;
            }
            int nbases = f.nsummaries;
            char **bases = calloc((size_t)nbases + 1, sizeof(char *));
            if (!bases) return md_fail(MDE_MEMORY, "out of memory");
            for (int k = 0; k < f.nsummaries; k++) {
                bases[k] = site_pg_irrep_label(f.site_pg, (int)f.summaries[k].pg_irrep);
                if (!bases[k]) return md_err.kind;
            }
            MD_TRY(topology_for(r, st, site_index, &f));
            int plan = -1;
            double *block_vectors = 0;
            int natoms = 0;
            if (mode_count) {
                int presented = raw_basis && !f.adopted;
                MD_TRY(plan_for(r, st, site_index, &f, f.lanes, presented, mag, &plan));
                const md_plan *p = st->plans + plan;
                if (f.lanes != p->vector_count)
                    return md_fail(MDE_VALUE, p->has_inverse ? "raw mode vectors and Source operation records differ"
                                                             : "mode rows do not cover the canonical presentation layout");
                natoms = p->natoms;
                block_vectors = malloc(sizeof(double) * (size_t)(mode_count * natoms * 3 + 1));
                double *scratch = malloc(sizeof(double) * (size_t)(3 * (f.lanes > 0 ? f.lanes : 1)));
                if (!block_vectors || !scratch) return md_fail(MDE_MEMORY, "out of memory");
                double pre[9];
                int has_pre = !presented && r->has_presentation_inverse;
                for (int i = 0; i < 9; i++) pre[i] = r->presentation_inverse[i / 3][i % 3];
                int rows_status = mode_rows_from_plan(mode_count, f.lanes, f.modes, has_pre, pre, p->has_inverse, p->inverse, natoms,
                                                      p->raw_indices, p->offsets, p->equivalents, p->has_negated ? p->negated : 0,
                                                      scratch, block_vectors);
                free(scratch);
                if (rows_status == 1) return md_fail(MDE_VALUE, "Source-equivalent rows carry different mode vectors");
                if (rows_status) return md_fail(MDE_VALUE, "mode atom row contains a nonfinite position or vector");
            }
            if (mode_count) {
                GROW(st->blocks, st->nblocks, st->block_cap);
                md_block *b = st->blocks + st->nblocks;
                memset(b, 0, sizeof *b);
                b->spec_order = spec_index;
                b->site_order = site_index;
                b->magnetic = mag;
                b->site = site;
                const jvalue *label = jget(site->site, "label"), *type = jget(site->site, "type"), *wyck = jget(site->site, "wyckoff");
                char *atom_label = jtruthy(label) ? md_py_str(label) : jtruthy(type) ? md_py_str(type) : (char *)md_strdup("Atom", 4);
                char *wyckoff = jtruthy(wyck) ? md_py_str(wyck) : (char *)md_strdup("?", 1);
                b->kvector = md_format("(%s)", st->text_kvector[spec_index]);
                b->label_prefix = md_format("%s[%s]%s%s[%s:%s:%s]", r->parent_symbol, st->text_kvector[spec_index], spec->irrep_label,
                                            st->text_direction[spec_index], atom_label, wyckoff, mag ? "mag" : "dsp");
                free(atom_label);
                free(wyckoff);
                b->direction = st->text_direction[spec_index];
                b->k = spec->k_label;
                b->irrep = spec->irrep_label;
                b->gid = spec_gid;
                b->subspace_replaced = spec->replaced_state == MD_SET && spec->replaced;
                b->spec = spec;
                b->stokes = spec->stokes;
                b->dynamic = dynamic;
                b->phase = phase;
                const int64_t *kp = spec->carrier_state == MD_SET ? spec->carrier : spec->kp_state == MD_SET ? spec->kp : 0;
                b->nkp = kp ? 4 : 0;
                if (kp) memcpy(b->kp, kp, sizeof b->kp);
                b->summaries = f.summaries;
                b->nsummaries = f.nsummaries;
                f.summaries = 0;
                f.nsummaries = 0;
                if (!b->kvector || !b->label_prefix) return md_fail(MDE_MEMORY, "out of memory");
                int block_index = st->nblocks++;
                int first = st->nmodes;
                for (int position = 0; position < mode_count; position++) {
                    const int64_t *row = layout_modes + 7 * position;
                    GROW(st->modes, st->nmodes, st->mode_cap);
                    md_mode *m = st->modes + st->nmodes++;
                    memset(m, 0, sizeof *m);
                    m->block = block_index;
                    char *component = component_text(row[2]);
                    const char *base = bases[row[0]];
                    m->site_irrep = row[1] ? md_format("%s_%lld(%s)", base, (long long)row[1], component) : md_format("%s(%s)", base, component);
                    free(component);
                    if (display_wanted) {
                        char *display_component = component_text(row[4]);
                        m->display = row[3] ? md_format("%s_%lld(%s)", base, (long long)row[3], display_component)
                                            : md_format("%s(%s)", base, display_component);
                        free(display_component);
                    }
                    if (identity_covered) {
                        const md_summary *s = b->summaries + row[0];
                        m->identity = (md_identity){1, s->gid, f.site_pg, s->pg_irrep, s->family, row[6], row[5], s->flat, s->source_row_count,
                                                    s->little_type, (int)row[0]};
                    }
                    m->vectors = malloc(sizeof(double) * (size_t)(natoms * 3 + 1));
                    if (!m->site_irrep || !m->vectors) return md_fail(MDE_MEMORY, "out of memory");
                    memcpy(m->vectors, block_vectors + (size_t)position * natoms * 3, sizeof(double) * (size_t)(natoms * 3));
                    m->plan = plan;
                }
                if (candidate) {
                    GROW(st->pending, st->npending, st->pending_cap);
                    md_pending *pd = st->pending + st->npending++;
                    pd->spec = spec_index;
                    pd->block = block_index;
                    pd->first = first;
                    pd->count = mode_count;
                    pd->keys = keys;
                    pd->rows = family_rows;
                    keys = family_rows = 0;
                } else {
                    MD_TRY(md_present(r, st, first, mode_count, keys, family_rows));
                }
            }
            free(block_vectors);
            for (int k = 0; k < nbases; k++) free(bases[k]);
            free(bases);
            if (candidate) st->pending_counts[(size_t)spec_index * r->nsites + site_index] = mode_count;
            else {
                site_definitions += mode_count;
                if (f.atom_count > site_atoms) site_atoms = f.atom_count;
            }
            free(layout_modes);
            free(keys);
            free(family_rows);
            frame_free(&f);
        }
        st->site_definitions[site_index] = site_definitions;
        st->site_atoms[site_index] = site_atoms;
    }
    return 0;
}

/* definition_table.present_block over modes[first .. first + count) */
static int md_present(md_request *r, md_state *st, int first, int count, const int64_t *keys, const int64_t *family_rows)
{
    if (!count) return 0;
    const md_plan *plan = st->plans + st->modes[first].plan;
    if (!plan->natoms) return md_fail(MDE_VALUE, "mode definition has no complete atom identity");
    for (int i = 0; i < count; i++)
        if (!st->modes[first + i].identity.present) return md_fail(MDE_KEY, "'site_pg'");
    int natoms = plan->natoms;
    double *vectors = malloc(sizeof(double) * (size_t)(count * natoms * 3 + 1));
    double *normfactor = malloc(sizeof(double) * (size_t)count);
    int *has = malloc(sizeof(int) * (size_t)count), *perm = malloc(sizeof(int) * (size_t)count);
    md_identity *identities = malloc(sizeof(md_identity) * (size_t)count);
    if (!vectors || !normfactor || !has || !perm || !identities) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < count; i++) {
        memcpy(vectors + (size_t)i * natoms * 3, st->modes[first + i].vectors, sizeof(double) * (size_t)(natoms * 3));
        identities[i] = st->modes[first + i].identity;
    }
    double basis[9];
    for (int i = 0; i < 9; i++) basis[i] = r->presentation_basis[i / 3][i % 3];
    int status = present_block(count, natoms, vectors, keys, family_rows, basis, r->has_centering, r->centering_factor, normfactor, has, perm);
    if (status) {
        free(vectors);
        free(normfactor);
        free(has);
        free(perm);
        free(identities);
        return status == 1 ? md_fail(MDE_VALUE, "mode family atom identity/order changed") : md_fail(MDE_MEMORY, "native block presentation failed");
    }
    for (int j = 0; j < count; j++) {
        md_mode *m = st->modes + first + j;
        int source = perm[j];
        m->has_normfactor = has[source];
        m->normfactor = normfactor[source];
        m->identity = identities[source];
        memcpy(m->vectors, vectors + (size_t)source * natoms * 3, sizeof(double) * (size_t)(natoms * 3));
    }
    free(vectors);
    free(normfactor);
    free(has);
    free(perm);
    free(identities);
    return 0;
}

/* ------------------------------------------------------------ post */

static const md_emission *emission_of(const md_state *st, int spec, int site)
{
    for (int e = 0; e < st->nemissions; e++)
        if (st->emissions[e].spec == spec && st->emissions[e].site == site) return st->emissions + e;
    return 0;
}

static int sites_of(const md_state *st, const md_request *r, int spec, int *out)
{
    int n = 0;
    for (int i = 0; i < st->neligible; i++)
        if (st->eligible[i] / r->nsites == spec) out[n++] = (int)(st->eligible[i] % r->nsites);
    return n;
}

static int emissions_equal(const md_emission *a, const md_emission *b)
{
    return a->row_id == b->row_id && !strcmp(a->signature, b->signature) && a->nrecords == b->nrecords &&
           !memcmp(a->records, b->records, sizeof(int64_t) * 5 * (size_t)a->nrecords) && a->nfractionals == b->nfractionals &&
           !memcmp(a->fractionals, b->fractionals, sizeof(int64_t) * 6 * (size_t)a->nfractionals) && a->mode_count == b->mode_count;
}

static int admitted_orders(const md_request *r, const md_state *st, const md_aliases *al, int *admitted)
{
    int *a_sites = malloc(sizeof(int) * (size_t)(st->neligible + 1)), *b_sites = malloc(sizeof(int) * (size_t)(st->neligible + 1));
    if (!a_sites || !b_sites) return md_fail(MDE_MEMORY, "out of memory");
    for (int c = 0; c < r->nspecs; c++) {
        int anchor = al->anchor_of[c];
        if (anchor < 0) continue;
        if (anchor != c) {
            int na = sites_of(st, r, c, a_sites), nb = sites_of(st, r, anchor, b_sites);
            if (na && na == nb && !memcmp(a_sites, b_sites, sizeof(int) * (size_t)na)) {
                int emitted = 0, ok = 1;
                for (int i = 0; i < na && ok; i++) {
                    const md_emission *x = emission_of(st, c, a_sites[i]), *y = emission_of(st, anchor, a_sites[i]);
                    if (!x || !y || !emissions_equal(x, y) || !x->all_nonzero) ok = 0;
                    else emitted += x->mode_count;
                }
                if (ok && emitted > 0) admitted[c] = 1;
            }
        }
        if (al->distinct[c]) {
            int n = sites_of(st, r, c, a_sites);
            if (n) {
                int ok = 1;
                for (int i = 0; i < n && ok; i++) {
                    const md_emission *x = emission_of(st, c, a_sites[i]);
                    ok = x && x->mode_count > 0 && x->all_nonzero && x->signature[0];
                }
                if (ok) admitted[c] = 1;
            }
        }
    }
    free(a_sites);
    free(b_sites);
    return 0;
}

static int occurrences_equal(const md_spec *a, const md_spec *b)
{
    if (a->nocc != b->nocc) return 0;
    for (int i = 0; i < a->nocc; i++) {
        const md_occurrence *x = a->occ + i, *y = b->occ + i;
        if (x->mesh != y->mesh || strcmp(x->k_label, y->k_label) || memcmp(x->kp, y->kp, sizeof x->kp) || memcmp(x->bucket, y->bucket, sizeof x->bucket))
            return 0;
        for (int axis = 0; axis < 3; axis++)
            if (!rat_eq(x->pml[axis], y->pml[axis])) return 0;
    }
    return 1;
}

/* definition_table.independent_child_mode_union over the ordered mode list */
static int independent_union(md_request *r, md_state *st, int *list, int *count)
{
    int n = *count;
    int *affected = calloc((size_t)n + 1, sizeof(int));
    if (!affected) return md_fail(MDE_MEMORY, "out of memory");
    int any = 0;
    for (int i = 0; i < n; i++) {
        const md_block *b = st->blocks + st->modes[list[i]].block;
        if (!b->subspace_replaced || !b->spec->nocc) continue;
        int specs = 0;
        int *seen = calloc((size_t)r->nspecs + 1, sizeof(int));
        if (!seen) return md_fail(MDE_MEMORY, "out of memory");
        for (int j = 0; j < n; j++) {
            const md_block *o = st->blocks + st->modes[list[j]].block;
            if (o->site_order == b->site_order && o->magnetic == b->magnetic && o->spec->nocc && occurrences_equal(o->spec, b->spec) &&
                !seen[o->spec_order]) {
                seen[o->spec_order] = 1;
                specs++;
            }
        }
        free(seen);
        if (specs > 1) affected[i] = 1;
    }
    for (int i = 0; i < n; i++) any |= affected[i];
    if (!any) {
        free(affected);
        return 0;
    }
    int *dependent = calloc((size_t)n + 1, sizeof(int));
    int *done = calloc((size_t)n + 1, sizeof(int));
    if (!dependent || !done) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < n; i++) {
        if (!affected[i] || done[i]) continue;
        const md_block *b = st->blocks + st->modes[list[i]].block;
        /* every mode of this (site, kind) group */
        int *members = malloc(sizeof(int) * (size_t)n);
        int nmembers = 0;
        if (!members) return md_fail(MDE_MEMORY, "out of memory");
        for (int j = 0; j < n; j++) {
            const md_block *o = st->blocks + st->modes[list[j]].block;
            int in_group = 0;
            for (int k = 0; k < n && !in_group; k++) {
                const md_block *q = st->blocks + st->modes[list[k]].block;
                in_group = affected[k] && q->site_order == o->site_order && q->magnetic == o->magnetic;
            }
            if (in_group && o->site_order == b->site_order && o->magnetic == b->magnetic) {
                members[nmembers++] = j;
                done[j] = 1;
            }
        }
        int natoms = st->plans[st->modes[list[members[0]]].plan].natoms;
        int64_t *row_of = malloc(sizeof(int64_t) * (size_t)(3 * natoms + 1));
        int64_t *offsets = malloc(sizeof(int64_t) * (size_t)(nmembers + 1));
        int *rows = malloc(sizeof(int) * (size_t)(nmembers * 3 * natoms + 1));
        int64_t *values = malloc(sizeof(int64_t) * (size_t)(nmembers * 3 * natoms + 1)), *quantized = malloc(sizeof(int64_t) * (size_t)(3 * natoms + 1));
        int *pivots = malloc(sizeof(int) * (size_t)(nmembers + 1));
        if (!row_of || !offsets || !rows || !values || !quantized || !pivots) return md_fail(MDE_MEMORY, "out of memory");
        for (int k = 0; k < 3 * natoms; k++) row_of[k] = -1;
        int64_t used = 0, nrows = 0;
        offsets[0] = 0;
        for (int m = 0; m < nmembers; m++) {
            const md_mode *mode = st->modes + list[members[m]];
            if (st->plans[mode->plan].natoms != natoms) return md_fail(MDE_VALUE, "mode definition vectors do not match its atoms");
            if (quantize4(3 * natoms, mode->vectors, quantized)) return md_fail(MDE_VALUE, "value cannot be quantized");
            for (int k = 0; k < 3 * natoms; k++) {
                if (!quantized[k]) continue;
                if (row_of[k] < 0) row_of[k] = nrows++;
                rows[used] = (int)row_of[k];
                values[used++] = quantized[k];
            }
            offsets[m + 1] = used;
        }
        int pivot_count = 0;
        int status = independent_columns((int)nrows, nmembers, offsets, rows, values, pivots, &pivot_count);
        if (status) return status == 1 ? md_fail(MDE_MEMORY, "independent mode columns are out of range")
                                       : md_fail(MDE_ARITH, "independent mode columns are out of range");
        int *is_pivot = calloc((size_t)nmembers + 1, sizeof(int));
        if (!is_pivot) return md_fail(MDE_MEMORY, "out of memory");
        for (int p = 0; p < pivot_count; p++) is_pivot[pivots[p]] = 1;
        for (int m = 0; m < nmembers; m++)
            if (!is_pivot[m]) dependent[members[m]] = 1;
        free(is_pivot);
        free(members);
        free(row_of);
        free(offsets);
        free(rows);
        free(values);
        free(quantized);
        free(pivots);
    }
    int kept = 0;
    for (int i = 0; i < n; i++)
        if (!dependent[i]) list[kept++] = list[i];
    *count = kept;
    free(affected);
    free(dependent);
    free(done);
    return 0;
}

static int compare_modes_state_key(const md_state *st, int a, int b)
{
    const md_block *x = st->blocks + st->modes[a].block, *y = st->blocks + st->modes[b].block;
    if (x->spec_order != y->spec_order) return x->spec_order < y->spec_order ? -1 : 1;
    if (x->site_order != y->site_order) return x->site_order < y->site_order ? -1 : 1;
    return 0;
}

static void stable_sort(const md_state *st, int *list, int n)
{
    for (int i = 1; i < n; i++) {
        int value = list[i], j = i;
        while (j > 0 && compare_modes_state_key(st, list[j - 1], value) > 0) {
            list[j] = list[j - 1];
            j--;
        }
        list[j] = value;
    }
}

static void put_frac_text(jbuf *b, rat v)
{
    char text[96];
    if (v.d == 1) snprintf(text, sizeof text, "%lld", (long long)v.n);
    else snprintf(text, sizeof text, "%lld/%lld", (long long)v.n, (long long)v.d);
    jb_cstring(b, text);
}

static void put_value_or_null(jbuf *b, const jvalue *v)
{
    if (v) jb_value(b, v);
    else jb_lit(b, "null");
}

static void put_display_basis(jbuf *b, const md_request *r)
{
    if (!r->has_display_basis) {
        jb_lit(b, "null");
        return;
    }
    jb_lit(b, "[");
    for (int row = 0; row < 3; row++) {
        if (row) jb_lit(b, ", ");
        jb_lit(b, "[");
        for (int col = 0; col < 3; col++) {
            if (col) jb_lit(b, ", ");
            jb_float(b, r->display_basis[row][col]);
        }
        jb_lit(b, "]");
    }
    jb_lit(b, "]");
}

static int secondary_factors(md_request *r, md_state *st, const int *emitted, jbuf *b)
{
    int64_t slot_count = jis(r->slots, J_ARRAY) && r->slots->count ? r->slots->count : 1;
    jb_lit(b, "[");
    int count = 0;
    for (int i = 0; i < r->nspecs; i++) {
        const md_spec *s = r->specs + i;
        if (!emitted[i] || s->primary || s->old_id <= 0) continue;
        if (s->gid_state != MD_SET || s->gid <= 0) return md_fail(MDE_VALUE, "rendered static secondary lost its Source irrep identity");
        if (count) jb_lit(b, ", ");
        jb_lit(b, "{\"slot\": ");
        jb_int(b, slot_count + count + 1);
        jb_lit(b, ", \"gid\": ");
        jb_int(b, s->gid);
        jb_lit(b, ", \"label\": ");
        jb_cstring(b, s->display_label ? s->display_label : s->label);
        jb_lit(b, s->magnetic_kind ? ", \"magnetic\": true" : ", \"magnetic\": false");
        jb_lit(b, ", \"k_parameters\": [], \"role\": \"secondary\"");
        if (s->free_state == MD_SET && s->free_count > 0) {
            jb_lit(b, ", \"parameter_count\": ");
            jb_int(b, s->free_count);
        }
        jb_lit(b, "}");
        count++;
    }
    jb_lit(b, "]");
    (void)st;
    return 0;
}

static char *parameter_name(int index)
{
    char parts[16];
    int n = 0;
    int value = index + 1;
    while (value) {
        int remainder = (value - 1) % 26;
        value = (value - 1) / 26;
        parts[n++] = (char)('a' + remainder);
    }
    char *out = malloc((size_t)n + 1);
    if (!out) return 0;
    for (int i = 0; i < n; i++) out[i] = parts[n - 1 - i];
    out[n] = 0;
    return out;
}

/* _attach_definition_factor_metadata for one definition */
int md_factor_metadata(const md_request *r, int spec_order, jbuf *b)
{
    if (spec_order < 0 || spec_order >= r->nspecs) return 0;
    const md_spec *s = r->specs + spec_order;
    if (!s->primary) {
        jb_lit(b, ", \"role\": \"secondary\"");
        return 0;
    }
    jb_lit(b, ", \"role\": \"primary\"");
    int64_t factors = jis(r->factors, J_ARRAY) ? r->factors->count : 0;
    int64_t order = -1;
    if (s->slot_order_state == MD_SET) order = s->slot_order;
    else if (factors == 1) order = 0;
    if (order < 0 || order >= factors) return 0;
    int64_t offset = 0;
    for (int64_t i = 0; i < order; i++) {
        int64_t count;
        if (py_int_of(jget(jat(r->factors, i), "parameter_count"), &count)) return md_fail(MDE_KEY, "'parameter_count'");
        offset += count;
    }
    const jvalue *factor = jat(r->factors, order);
    int64_t count;
    if (py_int_of(jget(factor, "parameter_count"), &count)) return md_fail(MDE_KEY, "'parameter_count'");
    jb_lit(b, ", \"factor\": {\"slot\": ");
    put_value_or_null(b, jget(factor, "slot"));
    jb_lit(b, ", \"gid\": ");
    put_value_or_null(b, jget(factor, "gid"));
    jb_lit(b, ", \"parameters\": [");
    for (int64_t i = 0; i < count; i++) {
        char *name = parameter_name((int)(offset + i));
        if (!name) return md_fail(MDE_MEMORY, "out of memory");
        if (i) jb_lit(b, ", ");
        jb_cstring(b, name);
        free(name);
    }
    jb_lit(b, "]}");
    return 0;
}

static void put_identity(jbuf *b, const md_state *st, const md_mode *m)
{
    const md_block *blk = st->blocks + m->block;
    if (!m->identity.present) {
        jb_lit(b, "null");
        return;
    }
    const md_identity *id = &m->identity;
    jb_lit(b, "{\"gid\": ");
    jb_int(b, id->gid);
    jb_lit(b, ", \"site_pg\": ");
    jb_int(b, id->site_pg);
    jb_lit(b, ", \"pg_irrep\": ");
    jb_int(b, id->pg_irrep);
    jb_lit(b, ", \"family\": ");
    jb_int(b, id->family);
    jb_lit(b, ", \"component\": ");
    jb_int(b, id->component);
    jb_lit(b, ", \"print_component\": ");
    jb_int(b, id->print_component);
    jb_lit(b, ", \"family_width\": ");
    jb_int(b, id->family_width);
    jb_lit(b, ", \"source_row_count\": ");
    jb_int(b, id->source_row_count);
    jb_lit(b, ", \"little_type\": ");
    jb_int(b, id->little_type);
    jb_lit(b, ", \"opd_groups\": ");
    summary_group_json(b, blk->summaries + id->summary);
    jb_lit(b, blk->dynamic ? ", \"dynamic\": true" : ", \"dynamic\": false");
    jb_lit(b, blk->phase ? ", \"phase\": true" : ", \"phase\": false");
    jb_lit(b, ", \"source_kparam\": [");
    for (int k = 0; k < blk->nkp; k++) {
        if (k) jb_lit(b, ", ");
        jb_int(b, blk->kp[k]);
    }
    jb_lit(b, "]}");
}

static int public_definitions(md_request *r, md_state *st, const int *list, int n, jbuf *b)
{
    jb_lit(b, "[");
    for (int i = 0; i < n; i++) {
        const md_mode *m = st->modes + list[i];
        const md_block *blk = st->blocks + m->block;
        const md_plan *plan = st->plans + m->plan;
        const char *public = m->display ? m->display : m->site_irrep;
        if (i) jb_lit(b, ", ");
        jb_lit(b, "{\"label\": ");
        char *label = md_format("%s%s", blk->label_prefix, public);
        if (!label) return md_fail(MDE_MEMORY, "out of memory");
        jb_cstring(b, label);
        free(label);
        jb_lit(b, ", \"mode\": {\"kind\": ");
        jb_cstring(b, blk->magnetic ? "magnetic" : "displacive");
        jb_lit(b, ", \"k\": ");
        jb_cstring(b, blk->k);
        jb_lit(b, ", \"kvector\": ");
        jb_cstring(b, blk->kvector);
        jb_lit(b, ", \"irrep\": ");
        jb_cstring(b, blk->irrep);
        jb_lit(b, ", \"gid\": ");
        jb_int(b, blk->gid);
        jb_lit(b, ", \"direction\": ");
        jb_cstring(b, blk->direction);
        const jvalue *type = jget(blk->site->site, "type"), *slabel = jget(blk->site->site, "label"), *wyck = jget(blk->site->site, "wyckoff");
        char *site_type = jtruthy(type) ? md_py_str(type) : (char *)md_strdup(blk->site->atom_label, (int64_t)strlen(blk->site->atom_label));
        char *site_label = jtruthy(slabel) ? md_py_str(slabel) : (char *)md_strdup(blk->site->atom_label, (int64_t)strlen(blk->site->atom_label));
        char *wyckoff = jtruthy(wyck) ? md_py_str(wyck) : (char *)md_strdup("", 0);
        if (!site_type || !site_label || !wyckoff) return md_fail(MDE_MEMORY, "out of memory");
        jb_lit(b, ", \"site_type\": ");
        jb_cstring(b, site_type);
        jb_lit(b, ", \"site_label\": ");
        jb_cstring(b, site_label);
        jb_lit(b, ", \"wyckoff\": ");
        jb_cstring(b, wyckoff);
        free(site_type);
        free(site_label);
        free(wyckoff);
        jb_lit(b, ", \"site_irrep\": ");
        jb_cstring(b, public);
        jb_lit(b, "}, \"normfactor\": ");
        if (m->has_normfactor) jb_float(b, m->normfactor);
        else jb_lit(b, "null");
        jb_lit(b, ", \"_spec_order\": ");
        jb_int(b, blk->spec_order);
        jb_lit(b, ", \"_site_order\": ");
        jb_int(b, blk->site_order);
        jb_lit(b, ", \"_mode_kind\": ");
        jb_cstring(b, blk->magnetic ? "mag" : "dsp");
        if (blk->stokes) jb_lit(b, ", \"_selected_stokes_columns\": true");
        jb_lit(b, ", \"_source_print_identity\": ");
        put_identity(b, st, m);
        jb_lit(b, ", \"_display_site_irrep_label\": ");
        if (m->display) jb_cstring(b, m->display);
        else jb_lit(b, "null");
        jb_lit(b, ", \"atom_block\": ");
        jb_int(b, plan->atom_block);
        jb_lit(b, ", \"vectors\": {\"$f64\": [");
        jb_int(b, (int64_t)md_vectors.used);
        jb_lit(b, ", ");
        jb_int(b, plan->natoms);
        jb_lit(b, "]}");
        MD_TRY(md_put_doubles(m->vectors, 3 * plan->natoms));
        MD_TRY(md_factor_metadata(r, blk->spec_order, b));
        jb_lit(b, "}");
    }
    jb_lit(b, "]");
    return 0;
}

md_doubles md_vectors;

int md_put_doubles(const double *values, int count)
{
    if (md_vectors.used + (size_t)count > md_vectors.cap) {
        size_t cap = md_vectors.cap ? md_vectors.cap : 4096;
        while (md_vectors.used + (size_t)count > cap) cap *= 2;
        double *grown = realloc(md_vectors.data, sizeof(double) * cap);
        if (!grown) return md_fail(MDE_MEMORY, "out of memory");
        md_vectors.data = grown;
        md_vectors.cap = cap;
    }
    memcpy(md_vectors.data + md_vectors.used, values, sizeof(double) * (size_t)count);
    md_vectors.used += (size_t)count;
    return 0;
}

static void put_undistorted(jbuf *b, const md_request *r);
static void put_lattice(jbuf *b, const double *lattice)
{
    static const char *keys[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
    jb_lit(b, "{");
    for (int i = 0; i < 6; i++) {
        if (i) jb_lit(b, ", ");
        jb_cstring(b, keys[i]);
        jb_lit(b, ": ");
        jb_float(b, lattice[i]);
    }
    jb_lit(b, "}");
}

static void put_sites(jbuf *b, const md_request *r, const md_state *st)
{
    jb_lit(b, "[");
    for (int i = 0; i < r->nsites; i++) {
        if (i) jb_lit(b, ", ");
        jb_lit(b, "{\"site\": ");
        put_value_or_null(b, jget(r->sites[i].site, "label"));
        jb_lit(b, ", \"status\": \"ok\", \"definition_count\": ");
        jb_int(b, st->site_definitions[i]);
        jb_lit(b, ", \"atom_count\": ");
        jb_int(b, st->site_atoms[i]);
        jb_lit(b, ", \"missing\": []}");
    }
    jb_lit(b, "]");
}

static int subgroup_details(jbuf *b, const md_request *r, const md_state *st)
{
    if (!st->site_isotropy && !r->iso) {
        jb_lit(b, "null");
        return 0;
    }
    jb_lit(b, "{\"row_id\": ");
    if (st->site_isotropy) {
        jb_int(b, st->iso_row_id);
        jb_lit(b, ", \"subgroup\": ");
        jb_int(b, st->iso_subgroup);
        jb_lit(b, ", \"arms\": ");
        jb_int(b, st->iso_arms);
        jb_lit(b, ", \"direction\": ");
        jb_int(b, st->iso_direction);
        jb_lit(b, ", \"basis\": [");
        for (int i = 0; i < 9; i++) {
            if (i) jb_lit(b, ", ");
            jb_int(b, st->iso_basis[i]);
        }
        jb_lit(b, "], \"basis_pml_to_cinter\": [");
        rat pml[3][3], origin[3];
        if (geo_setting_to_cinter(r->sg, 1, 0, pml, origin)) return md_fail(MDE_VALUE, "parent setting Source tables are inconsistent");
        for (int row = 0; row < 3; row++) {
            rat values[3], out[3];
            for (int k = 0; k < 3; k++) values[k] = rat_int(st->iso_basis[3 * row + k]);
            rat_row_mul(values, (const rat (*)[3])pml, out);
            if (row) jb_lit(b, ", ");
            jb_lit(b, "[");
            for (int col = 0; col < 3; col++) {
                if (col) jb_lit(b, ", ");
                put_frac_text(b, out[col]);
            }
            jb_lit(b, "]");
        }
        jb_lit(b, "]");
    } else {
        put_value_or_null(b, jget(r->iso, "row_id"));
        const char *keys[5] = {"subgroup", "arms", "direction", "basis", "basis_pml_to_cinter"};
        for (int i = 0; i < 5; i++) {
            jb_lit(b, ", ");
            jb_cstring(b, keys[i]);
            jb_lit(b, ": ");
            put_value_or_null(b, jget(r->iso, keys[i]));
        }
    }
    jb_lit(b, ", \"display_basis\": ");
    put_display_basis(b, r);
    jb_lit(b, ", \"presentation_basis\": ");
    put_display_basis(b, r);
    jb_lit(b, ", \"presentation_rule\": ");
    jb_cstring(b, !r->has_display_basis ? "none" : r->presentation_identity ? "identity" : "child_basis");
    jb_lit(b, ", \"origin\": ");
    if (r->origin_kind == 1) {
        jb_string(b, r->origin_text.data, (size_t)r->origin_text.length);
    } else if (r->origin_kind == 2) {
        jb_lit(b, "[");
        for (int i = 0; i < 4; i++) {
            if (i) jb_lit(b, ", ");
            jb_int(b, r->origin_record[i]);
        }
        jb_lit(b, "]");
    } else {
        jb_lit(b, "null");
    }
    jb_lit(b, ", \"magnetic_subgroup\": ");
    put_value_or_null(b, r->iso ? jget(r->iso, "subgroup") : 0);
    jb_lit(b, ", \"source_basis_values\": ");
    put_value_or_null(b, r->iso ? jget(r->iso, "source_basis_values") : 0);
    jb_lit(b, ", \"source_origin_values\": ");
    put_value_or_null(b, r->iso ? jget(r->iso, "source_origin_values") : 0);
    jb_lit(b, ", \"parent_inter_setting_id\": ");
    if (r->has_setting) jb_int(b, r->setting_id);
    else jb_lit(b, "null");
    jb_lit(b, "}");
    return 0;
}

static void put_undistorted(jbuf *b, const md_request *r)
{
    jbuf child = {0};
    jb_reset(&child);
    md_child_structure_json(r, &child);
    /* {"undistorted_atoms": [...], ...}: copy the first member's value */
    const char *start = strchr(child.data, '[');
    const char *end = strstr(child.data, "], \"site_params\"");
    if (start && end) jb_raw(b, start, (size_t)(end - start + 1));
    else jb_lit(b, "[]");
    jb_free(&child);
}

static struct {
    section_view irrep_pointer, irrep_label, irr_old, irr_sg;
    int ready;
} S;

/* modes.engine.source_catalog.little_record(sg, label).gid; 0 when it raises KeyError */
static int64_t little_gid_of_label(int sg, const char *label)
{
    if (!S.ready) {
        if (section_find("irreps.irrep_space_group_pointer", &S.irrep_pointer) || section_find("irreps.irrep_label", &S.irrep_label) ||
            section_find("little.little_irr_old", &S.irr_old) || section_find("little.little_irr_space_group", &S.irr_sg))
            return 0;
        S.ready = 1;
    }
    const char *a = label, *e = label + strlen(label);
    while (a < e && (*a == ' ' || (*a >= '\t' && *a <= '\r'))) a++;
    while (e > a && (e[-1] == ' ' || (e[-1] >= '\t' && e[-1] <= '\r'))) e--;
    int64_t start = section_int(&S.irrep_pointer, sg - 1), end = section_int(&S.irrep_pointer, sg), old_id = 0;
    for (int64_t oid = start; oid < end && !old_id; oid++) {
        const char *s;
        int64_t n;
        if (section_string(&S.irrep_label, oid - 1, &s, &n)) continue;
        while (n && (*s == ' ' || *s == '\t')) {
            s++;
            n--;
        }
        while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
        if (n == e - a && !memcmp(s, a, (size_t)n)) old_id = oid;
    }
    if (!old_id) return 0;
    int64_t found = 0, matches = 0, total = section_count(&S.irr_old);
    for (int64_t i = 0; i < total; i++)
        if (section_int(&S.irr_old, i) == old_id && section_int(&S.irr_sg, i) == sg) {
            found = i + 1;
            matches++;
        }
    return matches == 1 ? found : 0;
}

/* build_mode_details' strain_mode_specs: dsp spec copies with their order */
static int strain_mode_specs(md_request *r, md_spec **owned, const md_spec ***list, int **orders, int *count)
{
    int n = r->nspecs;
    *owned = calloc((size_t)n + 1, sizeof(md_spec));
    *list = calloc((size_t)n + 1, sizeof(md_spec *));
    *orders = calloc((size_t)n + 1, sizeof(int));
    if (!*owned || !*list || !*orders) return md_fail(MDE_MEMORY, "out of memory");
    *count = 0;
    for (int i = 0; i < n; i++) {
        if (r->specs[i].magnetic_kind) continue;
        md_spec *copy = *owned + *count;
        *copy = r->specs[i];
        copy->dir.values = 0;
        copy->rows.values = 0;
        copy->recomputed.values = 0;
        copy->occ = 0;
        copy->nocc = 0;
        copy->dir_state = r->specs[i].dir_state;
        if (r->specs[i].dir_state == MD_SET) {
            copy->dir.rows = r->specs[i].dir.rows;
            copy->dir.cols = r->specs[i].dir.cols;
            copy->dir.values = malloc(sizeof(double) * (size_t)(copy->dir.rows * copy->dir.cols + 1));
            if (!copy->dir.values) return md_fail(MDE_MEMORY, "out of memory");
            memcpy(copy->dir.values, r->specs[i].dir.values, sizeof(double) * (size_t)(copy->dir.rows * copy->dir.cols));
        }
        (*list)[*count] = copy;
        (*orders)[*count] = i;
        (*count)++;
    }
    const jvalue *rows = 0;
    if (r->iso) {
        rows = jget(r->iso, "source_numeric_rows");
        if (!jtruthy(rows)) rows = jget(r->iso, "source_matrix");
        if (!jtruthy(rows)) rows = 0;
    }
    if (rows)
        for (int i = 0; i < *count; i++) {
            md_spec *s = *owned + i;
            if (!s->primary || (s->dir_state == MD_SET && s->dir.rows)) continue;
            md_matrix m;
            int ndim;
            if (md_float_rows(rows, &m, &ndim)) {
                md_err.kind = 0;
                continue;
            }
            if (ndim != 2 || !m.cols) {
                free(m.values);
                continue;
            }
            free(s->dir.values);
            s->dir_state = MD_SET;
            s->dir.rows = m.cols;
            s->dir.cols = m.rows;
            s->dir.values = malloc(sizeof(double) * (size_t)(m.rows * m.cols + 1));
            if (!s->dir.values) return md_fail(MDE_MEMORY, "out of memory");
            for (int a = 0; a < m.rows; a++)
                for (int c = 0; c < m.cols; c++) s->dir.values[c * m.rows + a] = m.values[a * m.cols + c];
            free(m.values);
        }
    for (int i = 0; i < *count; i++) {
        md_spec *s = *owned + i;
        if (s->gid_state == MD_SET) continue;
        int64_t gid = little_gid_of_label(r->sg, s->label ? s->label : "");
        if (gid) {
            s->gid_state = MD_SET;
            s->gid = gid;
        }
    }
    if (!*count && jtruthy(jget(r->request, "include_strain")) && r->coupled) {
        int first = r->nspecs;
        MD_TRY(md_coupled_specs(r, 1, 0));
        for (int i = first; i < r->nspecs; i++) {
            if (r->specs[i].magnetic_kind) continue;
            md_spec *copy = realloc(*owned, sizeof(md_spec) * (size_t)(*count + 1));
            if (!copy) return md_fail(MDE_MEMORY, "out of memory");
            *owned = copy;
            (*owned)[*count] = r->specs[i];
            (*owned)[*count].dir.values = 0;
            (*owned)[*count].rows.values = 0;
            (*owned)[*count].recomputed.values = 0;
            (*owned)[*count].occ = 0;
            (*owned)[*count].nocc = 0;
            if (r->specs[i].dir_state == MD_SET) {
                md_matrix *d = &(*owned)[*count].dir;
                d->values = malloc(sizeof(double) * (size_t)(d->rows * d->cols + 1));
                if (!d->values) return md_fail(MDE_MEMORY, "out of memory");
                memcpy(d->values, r->specs[i].dir.values, sizeof(double) * (size_t)(d->rows * d->cols));
            }
            int *grown = realloc(*orders, sizeof(int) * (size_t)(*count + 1));
            if (!grown) return md_fail(MDE_MEMORY, "out of memory");
            *orders = grown;
            (*orders)[*count] = -1;
            (*count)++;
        }
        const md_spec **grown_list = realloc(*list, sizeof(md_spec *) * (size_t)(*count + 1));
        if (!grown_list) return md_fail(MDE_MEMORY, "out of memory");
        *list = grown_list;
        for (int i = 0; i < *count; i++) (*list)[i] = *owned + i;
        /* the copies own the strings of these specs; drop the rest */
        for (int i = first; i < r->nspecs; i++) {
            md_spec *s = r->specs + i;
            free(s->dir.values);
            free(s->rows.values);
            free(s->recomputed.values);
            for (int k = 0; k < s->nocc; k++) free(s->occ[k].k_label);
            free(s->occ);
        }
        r->nspecs = first;
    }
    return 0;
}

typedef struct {
    md_state st;
    md_aliases al;
    int ntext;
    md_spec *strain_specs;
    const md_spec **strain_list;
    int *strain_orders, nstrain;
    int *admitted, *lists[2], *emitted, *site_offset, *atom_named, *block_site;
    jbuf factors, strain;
} md_owned;

static void state_free(md_state *st, int ntext)
{
    for (int i = 0; i < st->nblocks; i++) {
        free(st->blocks[i].kvector);
        free(st->blocks[i].label_prefix);
        summaries_free(st->blocks[i].summaries, st->blocks[i].nsummaries);
    }
    free(st->blocks);
    for (int i = 0; i < st->nmodes; i++) {
        free(st->modes[i].site_irrep);
        free(st->modes[i].display);
        free(st->modes[i].vectors);
    }
    free(st->modes);
    for (int i = 0; i < st->nplans; i++) {
        free(st->plans[i].records);
        free(st->plans[i].raw_indices);
        free(st->plans[i].offsets);
        free(st->plans[i].equivalents);
        free(st->plans[i].negated);
    }
    free(st->plans);
    for (int i = 0; i < st->ntopologies; i++) free(st->topologies[i].records);
    free(st->topologies);
    for (int i = 0; i < st->nemissions; i++) {
        free(st->emissions[i].signature);
        free(st->emissions[i].records);
        free(st->emissions[i].fractionals);
    }
    free(st->emissions);
    for (int i = 0; i < st->npending; i++) {
        free(st->pending[i].keys);
        free(st->pending[i].rows);
    }
    free(st->pending);
    free(st->site_definitions);
    free(st->site_atoms);
    free(st->pending_counts);
    free(st->eligible);
    for (int i = 0; st->text_kvector && i < ntext; i++) free(st->text_kvector[i]);
    for (int i = 0; st->text_direction && i < ntext; i++) free(st->text_direction[i]);
    free(st->text_done);
    free(st->text_k_count);
    free(st->text_k);
    free(st->text_kvector);
    free(st->text_direction);
}

static void owned_free(md_owned *o)
{
    state_free(&o->st, o->ntext);
    free(o->al.candidate);
    free(o->al.anchor_of);
    free(o->al.observation);
    free(o->al.observed);
    free(o->al.distinct);
    for (int i = 0; o->strain_specs && i < o->nstrain; i++) free(o->strain_specs[i].dir.values);
    free(o->strain_specs);
    free(o->strain_list);
    free(o->strain_orders);
    free(o->admitted);
    free(o->lists[0]);
    free(o->lists[1]);
    free(o->emitted);
    free(o->site_offset);
    free(o->atom_named);
    free(o->block_site);
    jb_free(&o->factors);
    jb_free(&o->strain);
}

static int build(md_request *r, jbuf *b, md_owned *o)
{
    if (tables()) return md_fail(MDE_VALUE, "mode Source tables are inconsistent");
    MD_TRY(alias_orders(r, &o->al));
    MD_TRY(strain_mode_specs(r, &o->strain_specs, &o->strain_list, &o->strain_orders, &o->nstrain));
    MD_TRY(md_front(r));
    if (!r->has_child_sg || !r->has_display_basis) return md_fail(MDE_VALUE, "canonical child atom layout requires an exact presentation basis");
    MD_TRY(md_parent_setting_bridge(r->sg, r->has_setting, r->setting_id, &r->has_bridge, r->bridge_matrix, r->bridge_shift));
    MD_TRY(md_child_structure(r));
    r->has_source_basis_override = md_integer_source_basis(r, r->source_basis_override);
    int n = r->nspecs;
    o->ntext = n;
    o->st.text_done = calloc((size_t)n + 1, sizeof(int));
    o->st.text_k_count = calloc((size_t)n + 1, sizeof(int));
    o->st.text_k = calloc((size_t)n + 1, sizeof(rat[8]));
    o->st.text_kvector = calloc((size_t)n + 1, sizeof(char *));
    o->st.text_direction = calloc((size_t)n + 1, sizeof(char *));
    o->st.site_definitions = calloc((size_t)r->nsites + 1, sizeof(int));
    o->st.site_atoms = calloc((size_t)r->nsites + 1, sizeof(int));
    o->st.pending_counts = calloc((size_t)(n * r->nsites) + 1, sizeof(int64_t));
    if (!o->st.text_done || !o->st.text_k_count || !o->st.text_k || !o->st.text_kvector || !o->st.text_direction || !o->st.site_definitions || !o->st.site_atoms ||
        !o->st.pending_counts)
        return md_fail(MDE_MEMORY, "out of memory");
    MD_TRY(site_loop(r, &o->st, &o->al));
    o->admitted = calloc((size_t)n + 1, sizeof(int));
    if (!o->admitted) return md_fail(MDE_MEMORY, "out of memory");
    MD_TRY(admitted_orders(r, &o->st, &o->al, o->admitted));
    for (int spec = 0; spec < n; spec++) {
        if (!o->admitted[spec]) continue;
        for (int p = 0; p < o->st.npending; p++)
            if (o->st.pending[p].spec == spec) MD_TRY(md_present(r, &o->st, o->st.pending[p].first, o->st.pending[p].count, o->st.pending[p].keys, o->st.pending[p].rows));
        int *sites = malloc(sizeof(int) * (size_t)(o->st.neligible + 1));
        if (!sites) return md_fail(MDE_MEMORY, "out of memory");
        int count = sites_of(&o->st, r, spec, sites);
        for (int i = 0; i < count; i++) {
            int64_t value = o->st.pending_counts[(size_t)spec * r->nsites + sites[i]];
            if (value > 0 && sites[i] >= 0 && sites[i] < r->nsites) o->st.site_definitions[sites[i]] += (int)value;
        }
        free(sites);
    }
    int counts[2] = {0, 0};
    o->lists[0] = malloc(sizeof(int) * (size_t)(o->st.nmodes + 1));
    o->lists[1] = malloc(sizeof(int) * (size_t)(o->st.nmodes + 1));
    if (!o->lists[0] || !o->lists[1]) return md_fail(MDE_MEMORY, "out of memory");
    for (int i = 0; i < o->st.nmodes; i++) {
        const md_block *blk = o->st.blocks + o->st.modes[i].block;
        int pending = 0, admitted_spec = o->admitted[blk->spec_order];
        for (int p = 0; p < o->st.npending && !pending; p++)
            pending = o->st.pending[p].block == o->st.modes[i].block;
        if (pending && !admitted_spec) continue;
        o->lists[blk->magnetic][counts[blk->magnetic]++] = i;
    }
    /* Python appends o->admitted pending blocks after the loop, then sorts stably. */
    for (int kind = 0; kind < 2; kind++) {
        int *ordered = malloc(sizeof(int) * (size_t)(counts[kind] + 1));
        if (!ordered) return md_fail(MDE_MEMORY, "out of memory");
        int k = 0;
        for (int i = 0; i < counts[kind]; i++) {
            int pending = 0;
            for (int p = 0; p < o->st.npending && !pending; p++) pending = o->st.pending[p].block == o->st.modes[o->lists[kind][i]].block;
            if (!pending) ordered[k++] = o->lists[kind][i];
        }
        for (int spec = 0; spec < n; spec++) {
            if (!o->admitted[spec]) continue;
            for (int p = 0; p < o->st.npending; p++)
                if (o->st.pending[p].spec == spec && o->st.blocks[o->st.pending[p].block].magnetic == kind)
                    for (int m = 0; m < o->st.pending[p].count; m++) ordered[k++] = o->st.pending[p].first + m;
        }
        memcpy(o->lists[kind], ordered, sizeof(int) * (size_t)k);
        counts[kind] = k;
        free(ordered);
        stable_sort(&o->st, o->lists[kind], counts[kind]);
    }
    o->emitted = calloc((size_t)n + 1, sizeof(int));
    if (!o->emitted) return md_fail(MDE_MEMORY, "out of memory");
    for (int kind = 0; kind < 2; kind++)
        for (int i = 0; i < counts[kind]; i++) o->emitted[o->st.blocks[o->st.modes[o->lists[kind][i]].block].spec_order] = 1;
    jb_reset(&o->factors);
    MD_TRY(secondary_factors(r, &o->st, o->emitted, &o->factors));
    if (r->coupled) {
        MD_TRY(independent_union(r, &o->st, o->lists[0], &counts[0]));
        MD_TRY(independent_union(r, &o->st, o->lists[1], &counts[1]));
    }
    jb_reset(&o->strain);
    if (jtruthy(jget(r->request, "include_strain"))) MD_TRY(md_strain_definitions(r, o->strain_list, o->strain_orders, o->nstrain, &o->strain));
    else jb_lit(&o->strain, "[]");
    /* atom table and blocks */
    o->site_offset = malloc(sizeof(int) * (size_t)(r->nsites + 1));
    if (!o->site_offset) return md_fail(MDE_MEMORY, "out of memory");
    int total = 0;
    for (int i = 0; i < r->nsites; i++) {
        o->site_offset[i] = total;
        total += r->sites[i].natoms;
    }
    o->atom_named = calloc((size_t)total + 1, sizeof(int));
    int nblocks = 0;
    o->block_site = malloc(sizeof(int) * (size_t)(o->st.nplans + 1));
    if (!o->atom_named || !o->block_site) return md_fail(MDE_MEMORY, "out of memory");
    for (int kind = 0; kind < 2; kind++)
        for (int i = 0; i < counts[kind]; i++) {
            md_plan *plan = o->st.plans + o->st.modes[o->lists[kind][i]].plan;
            if (plan->atom_block >= 0) continue;
            for (int a = 0; a < plan->natoms; a++) o->atom_named[o->site_offset[plan->site] + a] = 1;
            int found = -1;
            for (int k = 0; k < nblocks && found < 0; k++)
                if (o->block_site[k] == plan->site) found = k;
            if (found < 0) {
                o->block_site[nblocks] = plan->site;
                found = nblocks++;
            }
            plan->atom_block = found;
        }
    int strain_empty = o->strain.used <= 2;
    if (!counts[0] && !counts[1] && strain_empty) {
        jb_reset(b);
        jb_lit(b, "{\"status\": \"missing\", \"reason\": ");
        const jvalue *symbol = jget(r->irrep, "symbol");
        if (!jtruthy(symbol)) symbol = jget(r->irrep, "label");
        char *target = jtruthy(symbol) ? md_py_str(symbol) : (char *)md_strdup("", 0);
        jbuf t = {0};
        jb_reset(&t);
        py_repr_text(&t, target, (int64_t)strlen(target));
        char *reason = md_format("selected irrep %s produced no local site definitions", t.data);
        jb_cstring(b, reason);
        free(reason);
        free(target);
        jb_free(&t);
        jb_lit(b, ", \"sites\": ");
        put_sites(b, r, &o->st);
        jb_lit(b, "}");
        jb_free(&o->factors);
        return 0;
    }
    jb_reset(b);
    jb_lit(b, "{\"status\": \"ok\", \"sites\": ");
    put_sites(b, r, &o->st);
    jb_lit(b, ", \"lattice\": ");
    put_lattice(b, r->child_lattice);
    jb_lit(b, ", \"viewer_parent_basis\": [");
    for (int row = 0; row < 3; row++) {
        if (row) jb_lit(b, ", ");
        jb_lit(b, "[");
        for (int col = 0; col < 3; col++) {
            if (col) jb_lit(b, ", ");
            jb_float(b, r->viewer_basis[row][col]);
        }
        jb_lit(b, "]");
    }
    jb_lit(b, "], \"viewer_parent_origin\": [");
    for (int axis = 0; axis < 3; axis++) {
        if (axis) jb_lit(b, ", ");
        jb_float(b, r->viewer_origin[axis]);
    }
    jb_lit(b, "], \"subgroup_details\": ");
    MD_TRY(subgroup_details(b, r, &o->st));
    jb_lit(b, ", \"undistorted_atoms\": ");
    put_undistorted(b, r);
    jb_lit(b, ", \"atoms\": [");
    int first_atom = 1;
    for (int i = 0; i < r->nsites; i++) {
        const md_site *s = r->sites + i;
        const jvalue *type = jget(s->site, "type"), *label = jget(s->site, "label");
        char *element = jtruthy(type) ? md_py_str(type) : jtruthy(label) ? md_py_str(label) : (char *)md_strdup("X", 1);
        if (!element) return md_fail(MDE_MEMORY, "out of memory");
        for (int a = 0; a < s->natoms; a++) {
            const md_child_site *c = s->children + s->atom_site[a];
            if (!first_atom) jb_lit(b, ", ");
            first_atom = 0;
            jb_lit(b, "{\"atom_id\": ");
            jb_cstring(b, s->atom_ids[a]);
            jb_lit(b, ", \"atom\": ");
            if (o->atom_named[o->site_offset[i] + a] && s->member_order[a] == 0) jb_cstring(b, c->id);
            else jb_lit(b, "null");
            jb_lit(b, ", \"label\": ");
            jb_cstring(b, c->id);
            jb_lit(b, ", \"child_site\": ");
            jb_cstring(b, c->id);
            jb_lit(b, ", \"wyckoff_site\": ");
            jb_cstring(b, c->wyckoff_site);
            jb_lit(b, ", \"element\": ");
            jb_cstring(b, element);
            jb_lit(b, ", \"site_order\": ");
            jb_int(b, i);
            jb_lit(b, ", \"atom_index\": ");
            jb_int(b, a);
            jb_lit(b, ", \"xyz\": [");
            for (int axis = 0; axis < 3; axis++) {
                if (axis) jb_lit(b, ", ");
                jb_float(b, s->positions[3 * a + axis]);
            }
            jb_lit(b, "]}");
        }
        free(element);
    }
    jb_lit(b, "], \"atom_blocks\": [");
    for (int k = 0; k < nblocks; k++) {
        if (k) jb_lit(b, ", ");
        jb_lit(b, "[");
        int site = o->block_site[k];
        for (int a = 0; a < r->sites[site].natoms; a++) {
            if (a) jb_lit(b, ", ");
            jb_int(b, o->site_offset[site] + a);
        }
        jb_lit(b, "]");
    }
    jb_lit(b, "], \"displacive_definitions\": ");
    MD_TRY(public_definitions(r, &o->st, o->lists[0], counts[0], b));
    jb_lit(b, ", \"strain_definitions\": ");
    jb_raw(b, o->strain.data, o->strain.used);
    jb_lit(b, ", \"magnetic_definitions\": ");
    MD_TRY(public_definitions(r, &o->st, o->lists[1], counts[1], b));
    jb_lit(b, ", \"secondary_opd_factors\": ");
    jb_raw(b, o->factors.data, o->factors.used);
    jb_lit(b, ", \"displacive_amplitudes\": []}");
    jb_free(&o->factors);
    jb_free(&o->strain);
    return b->failed ? md_fail(MDE_MEMORY, "out of memory") : 0;
}

int md_build(md_request *r, jbuf *b)
{
    md_owned o = {0};
    int status = build(r, b, &o);
    owned_free(&o);
    return status;
}

/* build_mode_details up to its render specs: the unsupported answers and
 * the parent cell; *done when out already holds the answer. */
int md_preflight(md_request *r, jbuf *b, int *done)
{
    *done = 0;
    int include_strain = jtruthy(jget(r->request, "include_strain"));
    if (!r->nsites && !include_strain) {
        jb_reset(b);
        jb_lit(b, "{\"status\": \"unsupported\", \"reason\": \"mode kernel requires at least one mapped displacive site\"}");
        *done = 1;
        return 0;
    }
    const jvalue *coupled = r->iso ? jget(r->iso, "coupled") : 0;
    int is_coupled = jtruthy(coupled) && jtruthy(r->slots);
    int64_t old_id = 0;
    const jvalue *old = jget(r->irrep, "old_id");
    if (jtruthy(old) && py_int_of(old, &old_id)) return md_fail(MDE_VALUE, "invalid literal for int()");
    int has_gid;
    int64_t gid;
    MD_TRY(md_dynamic_gid(r, &has_gid, &gid));
    if (!is_coupled && old_id <= 0 && !has_gid) {
        jb_reset(b);
        jb_lit(b, "{\"status\": \"unsupported\", \"reason\": \"mode kernel requires either a static old_id or a Source dynamic_gid\", \"irrep\": ");
        jb_value(b, r->irrep);
        jb_lit(b, "}");
        *done = 1;
        return 0;
    }
    MD_TRY(md_cell_params(r));
    const jvalue *symbol = jget(r->irrep, "symbol");
    if (!jtruthy(symbol)) symbol = jget(r->irrep, "label");
    if (!jtruthy(symbol)) {
        jb_reset(b);
        jb_lit(b, "{\"status\": \"unsupported\", \"reason\": \"mode kernel requires a selected irrep label\", \"irrep\": ");
        jb_value(b, r->irrep);
        jb_lit(b, "}");
        *done = 1;
    }
    return 0;
}
