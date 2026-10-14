/* The displayed K-family form shared by case_kparams.c and k_catalog.c. */

#ifndef APOSTRUCT_K_DISPLAY_H
#define APOSTRUCT_K_DISPLAY_H

#include <stdint.h>

#include "rational.h"

typedef struct {
    int count;
    rat d[4][3];
    char names[4];
    const char *fixed;
} kc_form;

int kc_display_records(int sg, int lattice, const char *symbol, int64_t symbol_length, const char *label, int64_t label_length,
                       const rat (*records)[3], int count, rat (*d)[3]);
int kc_display_form(int sg, int lattice, const char *symbol, int64_t symbol_length, const char *label, int64_t label_length,
                    const rat (*records)[3], int count, kc_form *form);

#endif
