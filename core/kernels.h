/* Entry points one kernel calls in another translation unit. */

#ifndef APOSTRUCT_KERNELS_H
#define APOSTRUCT_KERNELS_H

#include <stdint.h>

int subgroup_rows_in_context(int dim, int arms, const int64_t *k, int subtract, int element_count,
                             const int *elements, const double *sparse, const int *available, int identity_op,
                             const double *values, int64_t value_count, int row_count, int candidate_count,
                             const int64_t *candidates, const int *candidate_matrix_ops, const int *candidate_signs,
                             int max_iter, int64_t *basis_out, int64_t *kept_out, int *kept_count);
int same_subgroup(const int64_t *left_basis, int left_count, const int64_t *left, const int64_t *right_basis,
                  int right_count, const int64_t *right);
int irrep_matrix_stack(int gid, int param_count, const int64_t *kparam, int record_count, const int *records,
                       int capacity, double *out);
int id_subgroup_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag, int64_t *out);
int id_subgroup_magnetic_identify(int parent, const int64_t *basis, int op_count, const int64_t *ops, int flag,
                                  int64_t capacity, int64_t *out);
int wyckoff_child_sites(int parent_sg, int child_sg, int parent_row_id, const int64_t *basis, const int64_t *origin,
                        const int64_t *parent_params, int site_capacity, int atom_capacity, int64_t *sites_out,
                        int64_t *atoms_out, int *site_count_out, int *atom_count_out);
int child_public_wyckoff(int child_sg, int setting_count, const int *setting_ids, const int64_t *matrix_words,
                         const int64_t *origin_words, int point_count, const int64_t *xyz_words, int site_pg);
int64_t site_atom_table(int sg, int row_id, int site_pg, const int64_t *param_words, const int64_t *basis,
                        const int64_t *representative_words, int64_t capacity, int64_t *points, int64_t *records,
                        int64_t *order);
int present_rows(int count, const double *parent_xyz, const double *parent_dxyz, const double *basis,
                 const double *origin, int centering_letter, double tol, int capacity, int *out_source,
                 int *out_centering, double *out_xyz, double *out_dxyz, int *out_count);
int row_assignments(int child_sg, const int64_t *basis, int point_operation_count, int atom_count,
                    const int64_t *atoms, int row_count, const int64_t *rows, int *out);

int magnetic_site_groups(int group, int child_sg, int atom_count, const int64_t *parent_points,
                         const int64_t *atom_order, int site_count, const int64_t *site_offsets,
                         const int64_t *site_atoms, const int64_t *displayed, const int64_t *basis,
                         const int64_t *origin, int capacity, int64_t *components_out, int64_t *sites_out,
                         int *component_count, int *labels_ok);

int projection_entry(int sg, int row_id, int site_pg, int gid, int pg_irrep, int param_count, const int64_t *params,
                     const double *inverse, int64_t capacity, int64_t *out, int64_t *active_index,
                     double *active_value);
int display_families(int gid, int full_dim, int site_pg, int pg_irrep, int vector_basis_id, int vector_dim,
                     int vector_setting, int project_count, int param_count, const int64_t *params, int active_count,
                     const int64_t *active_index, const double *active_value, int direction_rows,
                     const double *direction, int lanes, const int *records, double *modes_out);
int recordwise_modes(int gid, int irrep_type, int full_dim, int old_id, int k_params, int param_count,
                     const int64_t *params, int pg_irrep, int project_count, int vector_dim, int vector_basis_id,
                     int site_pg, int vector_setting, int source_count, const int *sources, int active_count,
                     const int64_t *active_index, const double *active_value, int row_total, int row_stride,
                     const int *row_lengths, const double *rows, int source_row_id, int lanes, const int *records,
                     int transform_count, const double *transforms, int negate, int emit_freeparam, int freeparam,
                     int emit_all_rows, int emit_all_families, int64_t capacity, double *modes, int *summed,
                     int64_t summary_capacity, int64_t *summary);
void finalize_vectors(int64_t count, double *vectors, int transform_count, const double *transforms, int negate);

/* _vadd_fraction_operation and get_new_fractionals_ (orderparam statuses). */
int op_vadd(const int64_t *fraction, const int64_t *operation, int64_t *out);
int op_new_fractionals(const int64_t *matrix, int64_t *out, int capacity, int *count);
int traversal_tables(int space_count, const int64_t *translations, const int *ops, const int64_t *rotations,
                     int64_t rotation_den, int magnetic, const int64_t *point_matrices, const int *mlt,
                     int fraction_count, const int64_t *fractions, const int64_t *kernel_adjugate, int64_t kernel_det,
                     int64_t *residues, int64_t *modulus_out, int *space_product, int *space_fraction);

#endif
