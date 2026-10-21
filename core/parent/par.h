/* The parent stage (parent.*) in C: shared helpers of the par_*.c files. */

#ifndef APOSTRUCT_PAR_H
#define APOSTRUCT_PAR_H

#include <stdint.h>

#include "json.h"
#include "md.h"
#include "rational.h"
#include "tables.h"

int par_tables(void);
void par_put_stripped(jbuf *b, const section_view *view, int64_t index);
void par_preference_text(jbuf *b, int64_t id);
extern const char *PAR_DEFAULT_PREFERENCES;
int64_t par_wyckoff_multiplicity(int sg, int row_id);
void par_wyckoff_label(jbuf *b, int row_id);
void par_site_pg_label(jbuf *b, int site_pg);
/* SourceTables._wyckoff_formula_from_vectors: the formula as a JSON string; used[i] for x, y, z */
void par_formula(const rat (*v)[3], jbuf *formula, int *used);
int par_inter_vectors(int sg, int offset0, int setting_id, rat (*out)[3]);
double par_mod1(double x);
int par_distortion(const jvalue *state, const jvalue *mode_set, const jvalue *requested, int allow_empty, jbuf *out);
/* build_parent_state_from_cif_info */
int par_parent_state(const jvalue *cif_info, jbuf *out);

#endif
