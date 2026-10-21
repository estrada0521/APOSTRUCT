/* The get_irreps_ context of a little irrep at a k parameter (irreps.c). */

#ifndef APOSTRUCT_IRREPS_H
#define APOSTRUCT_IRREPS_H

#include <stdint.h>

enum {
    IR_OK = 0,
    IR_TABLES = 1,
    IR_NO_TABLES = 2,
    IR_RANGE = 3,
};

typedef struct {
    int dim, arms, subtract, element_count, identity_op;
    const int64_t *k;
    const int *elements;
    const double *sparse;
    const int *available;
} irrep_context;

/* kparam: the little irrep's k parameters (x, y, z, den), any shorter tuple
 * meaning none; the context stays valid until the next call for this gid. */
int irrep_context_for(int gid, int param_count, const int64_t *kparam, irrep_context *out);
/* The space group of an old (data_irreps) irrep, 0 when out of range, and the number of old irreps. */
int old_irrep_space_group(int old_id);
int old_irrep_count(void);
/* 1 when the tables list a little irrep of sg at k slot kslot. */
int little_k_has_irreps(int sg, int kslot);
/* little_irr_space_group of a little irrep, 0 when out of range. */
int little_irrep_sg(int gid);

int case_irrep_stack(int gid, int param_count, const int64_t *params, int record_count, const int *records,
                     int64_t capacity, double *out);

#endif
