/* The mode-details stage (modes.assembly.build_mode_details) in C: shared
 * request state and helpers of the md_*.c files. */

#ifndef APOSTRUCT_MD_H
#define APOSTRUCT_MD_H

#include <stdint.h>

#include "json.h"
#include "pyconv.h"
#include "rational.h"

/* Python exception classes the Python entry raises for a failed request. */
enum {
    MDE_VALUE = 1,
    MDE_KEY,
    MDE_INDEX,
    MDE_TYPE,
    MDE_ZERO,
    MDE_ARITH,
    MDE_MEMORY,
    MDE_RUNTIME,
    MDE_LINALG,
    /* a kernel status the Python entry maps through its own error table */
    MDE_CHILD_LAYOUT = 100,
    MDE_CHILD_SETTING,
    MDE_ORDERPARAM,
    MDE_MAGNETIC,
};

typedef struct {
    int kind;
    int status;
    char message[1 << 16];
} md_error;

extern md_error md_err;
int md_fail(int kind, const char *format, ...);
#define MD_TRY(call)                \
    do {                            \
        int md_status_ = (call);    \
        if (md_status_) return md_status_; \
    } while (0)

typedef struct {
    const char *data;
    int64_t length;
} md_str;

/* One child Wyckoff site of a parent site's layout. */
typedef struct {
    char *id;           /* label_prefix_N */
    char *wyckoff_site; /* multiplicity + label */
    rat representative[3];
    int first, count;   /* atom range */
} md_child_site;

/* One parent site and its canonical child layout. */
typedef struct {
    const jvalue *site;
    char *atom_label, *label_prefix;
    int has_parent_params;
    int nparams;
    char param_names[3][8];
    double param_values[3];
    const jvalue *param_raw[3]; /* the JSON value read, NULL when solved */
    int site_param_count; /* -1 for None */
    double site_params[3];
    int nchild;
    md_child_site *children;
    int natoms;
    char **atom_ids;
    int *atom_site;      /* child site order of each atom */
    int64_t *member_order;
    int64_t *topology;   /* natoms x 10: kernel[4], coset[5], centering */
    double *positions;   /* natoms x 3 */
    int presentation_rows;
} md_site;

/* Presence of a spec dict key: absent, None, or a value. */
enum { MD_ABSENT = 0, MD_NONE = 1, MD_SET = 2 };

typedef struct {
    int64_t mesh;
    char *k_label;
    int64_t kp[4];
    rat pml[3];
    int64_t bucket[6];
} md_occurrence;

typedef struct {
    int rows, cols;
    double *values; /* rows x cols */
} md_matrix;

/* One render spec: the spec dict of subduction_specs._specs (or the
 * fallback dict) with its mode kind, vector setting and display label. */
typedef struct {
    int kind;
    int64_t old_id;
    int gid_state;
    int64_t gid;
    char *label, *display_label, *k_label;
    int magnetic_state; /* MD_ABSENT, 1 False, 2 True */
    int row_id_state;
    int64_t row_id;
    int kp_state;
    int64_t kp[4];
    int carrier_state;
    int64_t carrier[4];
    int pml_state, pml_n;
    rat pml[3];
    int req_state, req_n;
    rat req[3];
    int occ_state, nocc;
    md_occurrence *occ;
    int dir_state;
    md_matrix dir;
    int rows_state;
    md_matrix rows;
    int primary;
    int free_state;
    int64_t free_count;
    int little_type_state;
    int64_t little_type;
    int coupled;
    int case_state, case_n;
    rat case_k[3];
    int dcase_state, dcase_n;
    rat dcase_k[3];
    int slot_order_state;
    int64_t slot_order;
    int replaced_state, replaced;
    int stokes;
    md_matrix recomputed;
    int obs;
    int has_anchor;
    int64_t anchor_gid, anchor_kp[4];
    int has_cand;
    int64_t cand_kp[4];
    int has_heading;
    int64_t heading_gid, heading_kp[4];
    char *heading_label;
    /* render tuple */
    int magnetic_kind; /* mode_kind == "mag" */
    int vector_setting;
    char *irrep_label;
} md_spec;

typedef struct {
    const jvalue *request, *cif, *selected_k, *irrep, *opd, *iso;
    const jvalue *slots, *factors, *k_params;
    int sg;
    int has_child_sg, child_sg;
    double cell[6];
    int has_setting;
    int64_t setting_id;
    int point_operation_count;
    /* selected presentation */
    int has_display_basis;
    double display_basis[3][3];
    int has_exact_basis;
    rat exact_basis[3][3];
    int origin_kind; /* 0 none, 1 text, 2 record */
    md_str origin_text;
    int64_t origin_record[4];
    int has_presentation_pml;
    int64_t presentation_pml[9];
    int64_t source_basis[9], source_origin[4];
    int magnetic_group;
    char centering_letter;
    int has_bridge;
    rat bridge_matrix[3][3], bridge_shift[3];
    rat parent_pml_inverse[3][3], parent_pml_origin[3];
    int setting_count;
    int setting_ids[32];
    rat child_matrix[3][3], child_origin[3];
    int nsites;
    md_site *sites;
    int64_t *label_memo;
    int coupled;
    int nspecs, spec_cap;
    md_spec *specs;
    /* front (build_mode_details before its site loop) */
    char *parent_symbol, *display_kvector, *opd_direction;
    int has_selected_gid;
    int64_t selected_gid;
    double child_lattice[6];
    double presentation_basis[3][3];
    int has_centering;
    double centering_factor;
    int has_presentation_inverse;
    double presentation_inverse[3][3];
    int has_freeparam;
    int64_t freeparam;
    int full_parameter_opd;
    int has_source_basis_override;
    int64_t source_basis_override[9];
    double viewer_origin[3];
} md_request;

int md_tables(void);
const char *md_strdup(const char *text, int64_t length);
char *md_format(const char *format, ...);
/* str(value) as an allocated C string. */
char *md_py_str(const jvalue *value);

/* selected_opd["isotropy"] or selected_opd itself, as _isotropy_from_opd_row. */
const jvalue *md_isotropy(const jvalue *opd);
int md_parent_setting_bridge(int sg, int has_setting, int64_t setting_id, int *has_bridge, rat (*matrix)[3], rat *shift);
void md_selected(md_request *r);
int md_child_structure(md_request *r);
int md_render_specs(md_request *r);
int md_front(md_request *r);
int md_dynamic_gid(const md_request *r, int *has, int64_t *gid);
/* Mode vectors leave the JSON as {"$f64": [offset, atoms]} into this buffer. */
typedef struct {
    double *data;
    size_t used, cap;
} md_doubles;
extern md_doubles md_vectors;
int md_put_doubles(const double *values, int count);
int md_build(md_request *r, jbuf *out);
int md_cell_params(md_request *r);
char *md_k_label_from_irrep_label(const char *label);
int md_coupled_specs(md_request *r, int include_displacive, int include_magnetic);
int md_factor_metadata(const md_request *r, int spec_order, jbuf *b);
int md_strain_definitions(md_request *r, const md_spec *const *specs, const int *orders, int nspecs, jbuf *b);
int md_preflight(md_request *r, jbuf *out, int *done);
int md_integer_source_basis(const md_request *r, int64_t *out);
int md_float_rows(const jvalue *rows, md_matrix *out, int *ndim);
char *md_default_symbol(int sg);
void md_origin_vector(const md_request *r, double *out);
int md_origin_rational(const md_request *r, rat *out);
void md_spec_free(md_spec *spec);
void md_child_atoms_json(const md_request *r, jbuf *b);
void md_release(md_request *r);

#endif
