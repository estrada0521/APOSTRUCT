/* Space-group, setting and Wyckoff geometry decoded from the Source blob. */

#ifndef APOSTRUCT_GEOMETRY_H
#define APOSTRUCT_GEOMETRY_H

#include <stdint.h>

#include "rational.h"

#define GEO_MAX_OPS 48

enum {
    GEO_OK = 0,
    GEO_TABLES = 1,
    GEO_RANGE = 2,
    GEO_COSETS = 3,
};

int geo_load(void);
/* The first error since the last call (GEO_RANGE for rational overflow), then clear. */
int geo_take_error(void);
/* The Source operation records (x, y, z, den, point op) of a space group. */
int geo_sg_records(int sg, int64_t (*out)[5]);
int64_t geo_point_group_order(int point_group);
int geo_sg_point_group(int sg);
/* iwyckoff_pg_elements settings of a site point group (zeros dropped). */
int geo_site_settings(int site_pg, int capacity, int (*out)[48], int *lengths);
int64_t geo_site_settings_pointer(int site_pg);
/* The pml point matrix of an operation (ipoint_op_psettings); 1 when Source has none. */
int geo_pml_point_matrix(int sg, int op, int64_t *out);
int geo_sg_lattice(int sg);
int64_t geo_point_op_code(int op);
void geo_point_op_matrix(int op, int64_t *out);
/* The 36-value lattice_ml record of a lattice; returns its first denominator. */
int64_t geo_ml_setting(int lattice, int64_t *record);
/* get_generators_ in the pml setting; -1 when a generator op has no record. */
int geo_pml_generators(int sg, int64_t (*out)[5]);
int64_t geo_point_product(int left, int right);
/* vrot_: a point operation of sg acting on fractional coordinates. */
void geo_vrot(int sg, int op, const rat *v, rat *out);
int geo_wyckoff_row_range(int sg, int *start, int *count);
int geo_wyckoff_site_pg(int row_id);
void geo_wyckoff_vectors(int row_id, rat (*out)[3]);
/* The PML formula vectors of a Wyckoff row offset in one inter setting. */
void geo_inter_wyckoff_vectors(int setting_id, int offset0, rat (*out)[3]);
int geo_default_setting(int sg);
/* The inter setting ids of a space group in Source order. */
int geo_inter_setting_ids(int sg, int capacity, int *out);
void geo_cinter_base(int sg, int half, int64_t *out);
void geo_inter_block(int setting_id, int half, int64_t *out);
/* wyc_pg_elements_: the site point-group records of one Wyckoff row. */
int geo_site_records(int sg, int row_id, int site_pg, int64_t (*out)[5]);
/* wyc_pg_cosets_: left coset representatives of the site group. */
int geo_cosets(int sg, const int64_t (*site)[5], int site_count, int64_t (*out)[5]);
int64_t geo_multiplicity(int sg, int site_pg);
/* The affine row transform of a Source setting (pml when pml, else cml) to cinter. */
int geo_setting_to_cinter(int sg, int pml, int setting_id, rat (*matrix)[3], rat *origin);
/* get_new_fractionals_: integer images of the quotient of basis, x-major. */
int64_t geo_kernel_fractions(const int64_t *basis, int64_t **out);
/* Lattice centering translations of sg's default setting, zero first. */
int geo_centering(int sg, rat (*out)[3]);
const int64_t *geo_point_op_codes(int64_t *count);
const int64_t *geo_point_group_codes(int64_t *count);
/* cos, sin of an angle in degrees and acos in degrees (degrees.c): exact at the
 * angles lattice symmetry fixes (0, 60, 90, 120, 180), elsewhere the same bits
 * on every platform (no libm) */
double geo_cos_deg(double degrees);
double geo_sin_deg(double degrees);
double geo_acos_deg(double cosine);

#endif
