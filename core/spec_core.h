/* Shared state of the mode-subduction spec stage (spec_rows.c, spec_build.c). */

#ifndef APOSTRUCT_SPEC_CORE_H
#define APOSTRUCT_SPEC_CORE_H

#include <stddef.h>
#include <stdint.h>

#include "rational.h"

/* Python exception classes carried by a failure. */
enum {
    PY_NONE = 0,
    PY_VALUE,
    PY_KEY,
    PY_INDEX,
    PY_TYPE,
    PY_ZERO,
    PY_ARITH,
    PY_MEMORY,
    PY_RUNTIME,
    PY_LINALG,
};

typedef struct {
    int kind;
    char message[160];
} sp_error;

extern sp_error sp_err;
int sp_fail(int kind, const char *message);
#define SP_TRY(expr)              \
    do {                          \
        if ((expr)) return 1;     \
    } while (0)

void *sp_alloc(size_t size);
void sp_reset(void);

typedef struct {
    const char *text;
    int64_t length;
} sp_str;

int sp_str_eq(sp_str a, sp_str b);

typedef struct {
    int64_t mesh;
    int64_t kp[4];
    rat pml[3];
    int64_t bucket[6];
    sp_str k_label;
} sp_occ;

typedef struct {
    int gid;
    sp_str irrep_name, k_label;
    int64_t kp[4];
    rat pml[3];
    int dim, count;
    const double *dir; /* dim x count, row-major */
    const int *occ;    /* indices into the occurrence pool */
    int nocc;
    int has_carrier;
    int64_t carrier[4];
} sp_row;

typedef struct {
    sp_row *items;
    int count, cap;
} sp_rows;

extern sp_occ *sp_occ_pool;
extern int sp_occ_count;

typedef struct {
    sp_str label;
    int64_t kp[4];
} sp_preferred;

typedef struct {
    int sg;
    int64_t basis[9];
    const int64_t (*ops)[5];
    int nops;
    int magnetic;
} sp_frame;

int sp_rows_push(sp_rows *rows, const sp_row *row);
int sp_frame_rows(const sp_frame *frame, const sp_preferred *preferred, int npreferred, sp_rows *out);
/* Alias rows hidden by orbit dedupe; requested may be null. */
int sp_aliases(const sp_frame *frame, const sp_row *row, const int64_t *requested, sp_rows *out);
/* kvec_standard_provenance(...).standard_kparam: *ok = 1 with out, else
 * *raised_kind is the Python class it raised (the caller decides whether it
 * catches that class). */
int sp_standard_kparam(int sg, int gid, const int64_t *kp, int64_t *out, int *ok, int *raised_kind);
/* The message of a kvec_standard class raised by k_family_standard_selection. */
const char *sp_standard_message(void);

int sp_little_old_id(int gid);
int sp_little_type(int gid);
int sp_little_full_dim(int gid);
sp_str sp_irrep_full_label(int gid);

/* source_kparam identity helpers */
int sp_kp_exact_eq(const int64_t *a, const int64_t *b);
int sp_kp_identity_eq(const int64_t *a, const int64_t *b);
int sp_kp_periodic_eq(const int64_t *a, const int64_t *b);
int sp_pml_eq(const rat *a, const rat *b);
int sp_pml_mod1_eq(const rat *a, const rat *b);

/* numpy.linalg.matrix_rank(m, tol) of a rows x cols row-major matrix. */
int sp_matrix_rank(const double *m, int rows, int cols, double tol, int *rank);
/* numpy.linalg.lstsq(a, b, rcond=None): x (acols x bcols) and its rank. */
int sp_lstsq(const double *a, int rows, int acols, const double *b, int bcols, double *x, int *rank);
/* a @ b as numpy computes it. */
void sp_matmul(const double *a, int rows, int inner, const double *b, int cols, double *out);

int case_kparams_of_pml(int gid, const rat *pml, rat *out);

#endif
