/* Python and gemmi value conversions the mode stage reproduces (io/pyconv.c). */

#ifndef APOSTRUCT_PYCONV_H
#define APOSTRUCT_PYCONV_H

#include <stdint.h>

#include "json.h"
#include "rational.h"

/* float(text): 0 with *out, 1 for ValueError. */
int py_float_of_text(const char *text, int64_t length, double *out);
/* Fraction(text): 0 with *out, 1 for ValueError or ZeroDivisionError, 2 when
 * the exact value leaves int128. */
int py_fraction_of_text(const char *text, int64_t length, rat *out);
/* Fraction(str(value)) of a float (exact decimal of its repr); 1 out of range. */
int py_fraction_of_float(double value, rat *out);
/* str(int) of an int128 / str(Fraction) ("n" or "n/d") into out (at least 96 bytes); returns the length. */
int py_i128_text(i128 value, char *out);
int py_rat_text(rat value, char *out);
/* float(Fraction): correctly rounded; 1 when n or d leave 2**53. */
int py_float_of_fraction(rat value, double *out);
/* Fraction.limit_denominator */
rat py_limit_denominator(rat value, i128 max_denominator);
/* float(value) / int(value) of a decoded JSON scalar; 1 for TypeError or ValueError. */
int py_float_of(const jvalue *value, double *out);
int py_int_of(const jvalue *value, int64_t *out);
/* str(value) of a decoded JSON scalar into b (strings verbatim). */
void py_str_of(jbuf *b, const jvalue *value);
/* a CIF number as a float: 0 with *out, 1 when absent */
int float_cif_number(const jvalue *value, double *out);

/* sum() of floats (Neumaier compensation, started from int 0). */
double py_sum(const double *values, int count);
/* repr(value) of a decoded JSON value / of a str. */
void py_repr_of(jbuf *b, const jvalue *value);
void py_repr_text(jbuf *b, const char *text, int64_t length);

/* gemmi.Op(text): rot and tran scaled by 24; 1 when gemmi would raise. */
int gemmi_parse_triplet(const char *text, int64_t length, int rot[3][3], int tran[3], char *notation);
/* Op.apply_to_xyz as the arm64 gemmi wheel evaluates it. */
void gemmi_apply_to_xyz(const int rot[3][3], const int tran[3], const double *xyz, double *out);

#endif
