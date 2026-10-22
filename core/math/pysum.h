/* CPython >= 3.12 sum() of three floats (Neumaier, starting from int 0). */

#ifndef APOSTRUCT_PYSUM_H
#define APOSTRUCT_PYSUM_H

#include <math.h>

static inline double apostruct_builtin_sum3(double a, double b, double c)
{
    double total = a + 0.0, compensation = 0.0;
    double values[2] = {b, c};
    for (int i = 0; i < 2; i++) {
        double value = values[i];
        double step = total + value;
        compensation = compensation + (fabs(total) >= fabs(value) ? (total - step) + value : (value - step) + total);
        total = step;
    }
    return compensation != 0.0 && isfinite(compensation) ? total + compensation : total;
}

#endif
