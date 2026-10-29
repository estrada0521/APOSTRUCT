/* The OPD selection stage (isotropy.catalog / coupled / selection) in C:
 * helpers shared by the sel_*.c files. */

#ifndef APOSTRUCT_SEL_H
#define APOSTRUCT_SEL_H

#include <stdint.h>

#include "json.h"
#include "md.h"

/* md_err.status carries an isotropy_rows_for status. */
#define MDE_ISOTROPY_ROWS 120
/* md_err.status carries an id_subgroup_magnetic_identify status. */
#define MDE_ID_SUBGROUP 121
/* md_err.status carries a subgroup operation / coupled domain orbit / bridge weight status. */
#define MDE_SUBGROUP_OPS 122
#define MDE_DOMAIN_ORBIT 123
#define MDE_BRIDGE 124
/* OverflowError("Python int too large to convert to C long") */
#define MDE_OVERFLOW 125
#define MDE_OPD_NOT_FOUND 126
/* md_err.status carries a vector_mode_total status. */
#define MDE_PROJECTION_ENTRY 111

/* The request k parameters as normalize_k_params left them (str -> str). */
typedef struct {
    int count;
    const char *names[16];
    const char *values[16];
} sel_params;

void sel_params_of(const jvalue *object, sel_params *out);
const char *sel_param(const sel_params *p, const char *name);

int sel_tables(void);
void sel_opd_from_rows(jbuf *b, char ***rows, int nrows, const int *widths);
void sel_display_opd(jbuf *b, const char *opd);
typedef struct {
    int64_t basis[3][3], denominator, origin[4], size;
} sel_setting;
int sel_change_setting(int parent_sg, int subgroup, const int64_t *basis, const int64_t *origin, int has_parent_setting,
                       int64_t parent_setting, int has_subgroup_setting, int64_t subgroup_setting, sel_setting *out);
char *sel_default_symbol(int sg);
int64_t sel_point_group_order(int sg);

typedef struct {
    int count;
    char **items;
} strings;
void strings_free(strings *s);
void sel_put_strings(jbuf *b, const strings *s);
void sel_put_joined(jbuf *b, const strings *s);
int sel_is_space(char c);
int sel_active_k_vectors(int parent_sg, int64_t gid, int64_t full_dim, char ***rows, int nrows, const int *widths, const sel_params *p,
                         strings *out);
int sel_active_k_vectors_by(int parent_sg, int64_t gid, int64_t full_dim, int nrows, const int *widths,
                            int (*nonzero)(const void *ctx, int row, int col), const void *ctx, const sel_params *p, strings *out);
int64_t sel_subgroup_index(int parent_sg, int subgroup, int64_t size, int *has);
int sel_display(int parent_sg, int subgroup, const int64_t *basis, const int64_t *origin, int has_setting, int64_t setting_id, int has_cell,
                const double *cell, int parametric, int coupled, int magnetic_group, sel_setting *displayed, rat (*rows)[3]);
void sel_put_subgroup(jbuf *b, int parent_sg, int subgroup, const int64_t *basis, const int64_t *origin, const sel_setting *displayed,
                      const rat (*rows)[3], int64_t det, const jbuf *subgroup_json, const jbuf *before_det);
int sel_mag_tables(void);
int sel_coupled_rows(const jvalue *request, jbuf *out);
typedef struct {
    const char *key;
    const char *json;
    size_t length;
} kv;
/* {**base, key: value, ...} */
void sel_merge(jbuf *b, const jvalue *base, const kv *items, int n);
/* the k parameters normalized, as JSON text */
void sel_normalize_params(jbuf *b, const jvalue *v);
int sel_k_points(const jvalue *state, const jvalue *opts, jbuf *out);
int sel_irreps(const jvalue *state, jbuf *out);
int sel_opd_factors(jbuf *b, int sg, const jvalue *slots, const jvalue *selected);
void sel_put_basis(jbuf *b, const rat (*rows)[3]);
int sel_order_parameter_direction(const jvalue *state, const jvalue *opts, jbuf *out);
void sel_put_basis_text(jbuf *b, const rat (*rows)[3]);
void sel_put_origin_text(jbuf *b, int parent_sg, int subgroup, const rat (*rows)[3], const int64_t *origin);
/* 1 with out, 0 when absent, -kind on error. */
int sel_dynamic_kparam(int64_t gid, const sel_params *p, int64_t *out);
int sel_mag_ordinary(int group);
void sel_mag_subgroup_json(jbuf *b, int group);
int64_t sel_mag_index(int parent_group, int subgroup_group, int64_t size, int *has);
int sel_mag_ferroic(jbuf *b, int group);
void sel_put_records(jbuf *b, const int64_t *records, int64_t count);
/* 1 with words (n, d) x3, 0 when absent, -kind on error. */
int sel_visible_k_vector(int sg, int kslot, const sel_params *p, int64_t *words);

/* the OPD rows of an irrep, as a JSON array */
int sel_dynamic_rows(int magnetic, int parent_sg, int64_t gid, const sel_params *k_params, const sel_params *display, const jvalue *selected_orderparam,
                     int has_setting, int64_t setting_id, int has_cell, const double *cell, jbuf *out);

#endif
