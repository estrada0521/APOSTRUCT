/* Source magnetic point operations and groups (data_magnetic) for other kernels. */

#ifndef APOSTRUCT_MAGNETIC_H
#define APOSTRUCT_MAGNETIC_H

#include <stdint.h>

int mag_tables(void);
/* The first table error since the last call, then clear. */
int mag_take_error(void);
int64_t mag_group_total(void);
int mag_ordinary_group(int group);
int64_t mag_group_code(int group);
int64_t mag_op_code(int op);
int mag_nonmagnetic_op(int op);
int mag_time_reversed(int op);
/* The matrix point operation and time-reversal sign of a point operation:
 * mode 0 ordinary, 1 magnetic irrep, 2 magnetic ignoring time reversal. */
int mag_point_action(int mode, int op, int *sign);
/* The product of two magnetic point operations, from their ordinary parts and time reversals. */
int mag_point_product(int left, int right);
int64_t mag_point_product_count(void);
/* The magnetic group of orderparam_to_subgroup_magnetic_ for a parent space group. */
int mag_orderparam_group(int sg);
int mag_generators(int group, int64_t (*out)[5]);
/* The generated records of a magnetic group in a setting (0 closure, 1 binary/pml, 2 cinter). */
int magnetic_records(int group, int setting, int capacity, int64_t *records_out);
/* SourceData.generate_magnetic_space_group_records. */
int magnetic_records_iso(int group, int capacity, int64_t *records_out);

#endif
