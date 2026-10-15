/* subgroup_to_orderparam_ for a little irrep over a list of operations: the
 * invariant directions from the reduced D - I rows (columns reversed),
 * computed exactly (exact_irrep.c). */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "irreps.h"
#include "kernels.h"
#include "magnetic.h"

int exact_subduction_directions(int gid, int param_count, const int64_t *params, int record_count, const int *mapped,
                                const int *signs, double *directions_out);

enum {
    SD_MEMORY = 1,
    SD_STACK = 2,
    SD_EXACT = 3,
};

/* mode_case must be 0 (params: the k parameter record).  point_mode: the records' point action
 * (mag_point_action).  directions_out: up to dim directions of dim
 * coordinates.  Returns the direction count or a negative status. */
int subduction_directions(int gid, int mode_case, int point_mode, int param_count, const int64_t *params,
                          int record_count, const int *records, double *directions_out)
{
    if (mode_case) return -SD_STACK;
    int *mapped = malloc(sizeof(int) * (5 * (size_t)record_count + 1)), *signs = malloc(sizeof(int) * ((size_t)record_count + 1));
    if (!mapped || !signs) {
        free(mapped);
        free(signs);
        return -SD_MEMORY;
    }
    int status = point_mode && mag_tables();
    for (int r = 0; r < record_count && !status; r++) {
        memcpy(mapped + 5 * r, records + 5 * r, sizeof(int) * 4);
        mapped[5 * r + 4] = mag_point_action(point_mode, records[5 * r + 4], signs + r);
    }
    if (!status && point_mode && mag_take_error()) status = 1;
    int count = status ? -SD_STACK : exact_subduction_directions(gid, param_count, params, record_count, mapped, signs, directions_out);
    if (count < 0 && !status) count = count == -1 ? -SD_STACK : -SD_EXACT;
    free(mapped);
    free(signs);
    return count;
}
