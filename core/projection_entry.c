/* The Source inputs of project_ for one little irrep and one site irrep of a
 * Wyckoff row: the project count (local408), the vector preparation, the
 * selected real basis blocks and the bridge source buffers. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "irreps.h"
#include "tables.h"

int project_selection(int type, int pd, int sd, int record_count, const double *parent, const double *site_rows,
                      const double *inverse, int block_count, int project_count, int basis_length, int capacity,
                      int *selected, int64_t *active_index, double *active_value, int *active_count);
int little_sparse_matrix(int gid, int dim, int position, double *out);
static int row_frequencies(int gid, int row_id, int64_t *frequency, int *site_irreps);

enum {
    PE_TABLES = 1,
    PE_KEY = 2,
    PE_VALUE = 3,
    PE_CAPACITY = 4,
};

enum { FAMILY_STRIDE = 144, SOURCE_LENGTH = 768 };

static struct {
    section_view sg, dim, type, old, kslot, star_count, pg_irrep, irrep_image, image_type,
        image_dimension, vector_reps;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (geo_load()) return PE_TABLES;
    struct {
        const char *name;
        section_view *view;
    } items[] = {
        {"little.little_irr_space_group", &T.sg},
        {"little.little_irr_full_dim", &T.dim},
        {"little.little_irr_type", &T.type},
        {"little.little_irr_old", &T.old},
        {"little.little_irr_k", &T.kslot},
        {"little.little_k_star_count", &T.star_count},
        {"wyckoff.iwyckoff_pg_irrep", &T.pg_irrep},
        {"irreps.irrep_image", &T.irrep_image},
        {"images.image_type", &T.image_type},
        {"images.image_dimension", &T.image_dimension},
        {"wyckoff.iwyckoff_pg_vector_reps", &T.vector_reps},
    };
    for (size_t i = 0; i < sizeof items / sizeof items[0]; i++)
        if (section_find(items[i].name, items[i].view) || !section_is_int(items[i].view)) return PE_TABLES;
    T.ready = 1;
    return 0;
}

static int bad;

static int64_t at(const section_view *view, int64_t index)
{
    if (index < 0 || index >= section_count(view)) {
        bad = 1;
        return 0;
    }
    return section_int(view, index);
}

static int site_old_id(int site_pg, int pg_irrep)
{
    int64_t old = at(&T.pg_irrep, (int64_t)(site_pg - 1) * 12 + pg_irrep - 1);
    return old > 0 ? (int)old : 0;
}

static int64_t image_of(int old_id) { return at(&T.irrep_image, old_id - 1); }

/* project_counts_for_little: the count of pg_irrep, after checking every
 * pair of the row's subduction. */
static int project_count_of(int gid, int row_id, int site_pg, int pg_irrep, int64_t *count)
{
    int type = (int)at(&T.type, gid - 1), n = 0;
    int64_t frequency[12];
    int status = row_frequencies(gid, row_id, frequency, &n);
    if (status) return status;
    int divisor = type == 1 ? 1 : type == 2 ? 4 : type == 3 ? 2 : 0;
    if (!divisor) return PE_KEY;
    int found = 0;
    for (int pg = 1; pg <= n; pg++) {
        if (!frequency[pg - 1]) continue;
        int site_old = site_old_id(site_pg, pg);
        if (!site_old) return PE_KEY;
        int64_t image = image_of(site_old);
        int64_t factor = (image - 1 < section_count(&T.image_type) && at(&T.image_type, image - 1) == 1) ? 1 : 2;
        int64_t numerator = frequency[pg - 1] * factor;
        if (numerator % divisor) return PE_VALUE;
        if (pg == pg_irrep) {
            *count = numerator / divisor;
            found = 1;
        }
    }
    if (bad) return PE_TABLES;
    return found ? 0 : PE_KEY;
}

static int vector_basis_id(int site_pg, const int64_t (*records)[5], int record_count, int64_t *out)
{
    static int settings[64][48];
    int lengths[64];
    int count = geo_site_settings(site_pg, 64, settings, lengths);
    int64_t pointer = geo_site_settings_pointer(site_pg);
    for (int pass = 0; pass < 2; pass++)
        for (int s = 0; s < count; s++) {
            int same = 1;
            if (!pass) {
                same = lengths[s] == record_count;
                for (int i = 0; i < record_count && same; i++) same = settings[s][i] == records[i][4];
            } else {
                for (int i = 0; i < lengths[s] && same; i++) {
                    int hit = 0;
                    for (int j = 0; j < record_count && !hit; j++) hit = settings[s][i] == records[j][4];
                    same = hit;
                }
                for (int j = 0; j < record_count && same; j++) {
                    int hit = 0;
                    for (int i = 0; i < lengths[s] && !hit; i++) hit = settings[s][i] == records[j][4];
                    same = hit;
                }
            }
            if (same) {
                *out = pointer + s;
                return 0;
            }
        }
    return PE_VALUE;
}

/* map_parent_ops_to_site_project_records: canonical site ops in parent order. */
static int site_ops_for(int site_pg, const int64_t (*records)[5], int record_count, int *out, int *out_count)
{
    static int settings[64][48];
    int lengths[64];
    int count = geo_site_settings(site_pg, 64, settings, lengths);
    if (!count) return PE_VALUE;
    for (int s = 0; s < count; s++) {
        if (lengths[s] != lengths[0]) continue;
        int mapped[GEO_MAX_OPS] = {0}, ok = 1;
        for (int i = 0; i < lengths[s] && ok; i++) {
            int index = -1;
            for (int j = 0; j < record_count && index < 0; j++)
                if (records[j][4] == settings[s][i]) index = j;
            if (index < 0) ok = 0;
            else mapped[index] = settings[0][i];
        }
        if (!ok) continue;
        *out_count = 0;
        for (int j = 0; j < record_count; j++)
            if (mapped[j]) out[(*out_count)++] = mapped[j];
        return 0;
    }
    return PE_VALUE;
}

static int little_gid_for_old(int old_id)
{
    static int *gid_by_old, old_limit;
    if (!gid_by_old) {
        int64_t total = section_count(&T.old);
        int limit = 0;
        for (int64_t i = 0; i < total; i++)
            if (section_int(&T.old, i) > limit) limit = (int)section_int(&T.old, i);
        int *index = calloc((size_t)limit + 1, sizeof *index);
        if (!index) return 0;
        for (int64_t i = 0; i < total; i++) {
            int64_t old = section_int(&T.old, i);
            if (old > 0) index[old] = index[old] ? -1 : (int)i + 1;
        }
        gid_by_old = index;
        old_limit = limit;
    }
    if (old_id <= 0 || old_id > old_limit) return 0;
    return gid_by_old[old_id] > 0 ? gid_by_old[old_id] : 0;
}

static int little_trace(int gid, int dim, int position, double *out)
{
    static double **traces;
    if (!traces && !(traces = calloc((size_t)section_count(&T.dim) + 1, sizeof *traces))) return PE_TABLES;
    if (gid < 1 || gid > section_count(&T.dim) || position < 1 || position > GEO_MAX_OPS) return PE_KEY;
    if (!traces[gid]) {
        if (!(traces[gid] = malloc(sizeof(double) * GEO_MAX_OPS))) return PE_TABLES;
        for (int i = 0; i < GEO_MAX_OPS; i++) traces[gid][i] = NAN;
    }
    double *t = &traces[gid][position - 1];
    if (isnan(*t)) {
        double *m = malloc(sizeof(double) * (size_t)dim * dim);
        if (!m) return PE_TABLES;
        int status = little_sparse_matrix(gid, dim, position, m);
        double sum = 0;
        for (int a = 0; a < dim && !status; a++) sum += m[a * dim + a];
        free(m);
        if (status) return PE_KEY;
        *t = sum;
    }
    *out = *t;
    return 0;
}

/* the subduction of a little irrep to a Wyckoff row's site group by characters: frequency[j - 1] of site irrep j is
 * sum_k tr W(k) tr S_j(k) / sum_k tr S_j(k)^2 over the site records k */
static int row_frequencies_uncached(int gid, int row_id, int64_t *frequency, int *site_irreps)
{
    int sg = (int)at(&T.sg, gid - 1), dim = (int)at(&T.dim, gid - 1), old = (int)at(&T.old, gid - 1);
    int site_pg = geo_wyckoff_site_pg(row_id);
    int64_t records[GEO_MAX_OPS][5];
    int record_count = geo_site_records(sg, row_id, site_pg, records);
    if (geo_take_error() || record_count <= 0 || dim <= 0 || dim > 48) return PE_TABLES;
    int records32[GEO_MAX_OPS * 5], site_ops[GEO_MAX_OPS], site_count = 0;
    for (int i = 0; i < record_count; i++)
        for (int j = 0; j < 5; j++) records32[5 * i + j] = (int)records[i][j];
    if (site_ops_for(site_pg, (const int64_t (*)[5])records, record_count, site_ops, &site_count) || site_count != record_count)
        return PE_VALUE;
    double *parent = malloc(sizeof(double) * (size_t)record_count * dim * dim), trace[GEO_MAX_OPS];
    const int64_t generic[6] = {1, 5, 2, 5, 3, 5};
    if (!parent) return PE_TABLES;
    if (case_irrep_stack(gid, old > 0 ? 0 : 3, generic, record_count, records32, (int64_t)record_count * dim * dim, parent)) {
        free(parent);
        return PE_VALUE;
    }
    for (int i = 0; i < record_count; i++) {
        trace[i] = 0;
        for (int a = 0; a < dim; a++) trace[i] += parent[(size_t)i * dim * dim + a * dim + a];
    }
    free(parent);
    int n = 0;
    for (int j = 1; j <= 12; j++) {
        int site_old = site_old_id(site_pg, j), sgid = site_old ? little_gid_for_old(site_old) : 0;
        if (!sgid) break;
        int sdim = (int)at(&T.dim, sgid - 1), ssg = (int)at(&T.sg, sgid - 1);
        int64_t generated[GEO_MAX_OPS][5];
        int generated_count = geo_sg_records(ssg, generated);
        double sum = 0, norm = 0;
        for (int i = 0; i < record_count; i++) {
            int position = 0;
            double ts;
            for (int g = 0; g < generated_count; g++)
                if (generated[g][4] == site_ops[i]) position = g + 1;
            if (!position || little_trace(sgid, sdim, position, &ts)) return PE_KEY;
            sum += trace[i] * ts;
            norm += ts * ts;
        }
        double value = norm > 0 ? sum / norm : -1;
        if (value < -1e-6 || fabs(value - nearbyint(value)) > 1e-6) return PE_VALUE;
        frequency[n++] = (int64_t)nearbyint(value);
    }
    *site_irreps = n;
    return n ? 0 : PE_KEY;
}

static int row_frequencies(int gid, int row_id, int64_t *frequency, int *site_irreps)
{
    typedef struct {
        int gid, row, n, status;
        int64_t frequency[12];
    } entry;
    static entry *table;
    static size_t size, used;
    if (used * 2 >= size) {
        size_t grown = size ? size * 2 : 1024;
        entry *next = calloc(grown, sizeof *next);
        if (!next) return PE_TABLES;
        for (size_t i = 0; i < size; i++)
            if (table[i].gid) {
                size_t h = ((size_t)table[i].gid * 7919u + (size_t)table[i].row) & (grown - 1);
                while (next[h].gid) h = (h + 1) & (grown - 1);
                next[h] = table[i];
            }
        free(table);
        table = next;
        size = grown;
    }
    size_t h = ((size_t)gid * 7919u + (size_t)row_id) & (size - 1);
    while (table[h].gid && (table[h].gid != gid || table[h].row != row_id)) h = (h + 1) & (size - 1);
    entry *e = &table[h];
    if (!e->gid) {
        e->status = row_frequencies_uncached(gid, row_id, e->frequency, &e->n);
        e->gid = gid;
        e->row = row_id;
        used++;
    }
    memcpy(frequency, e->frequency, sizeof e->frequency);
    *site_irreps = e->n;
    return e->status;
}

/* the site irreps of a Wyckoff row present in the subduction of a little irrep, ascending; their count, or < 0 */
int row_site_irreps(int gid, int row_id, int *out)
{
    if (load()) return -PE_TABLES;
    int64_t frequency[12];
    int n = 0, count = 0, status = row_frequencies(gid, row_id, frequency, &n);
    if (status) return -status;
    for (int j = 0; j < n; j++)
        if (frequency[j]) out[count++] = j + 1;
    return count;
}

/* Row 0 of the site irrep's Source matrix for each site op. */
static int site_rows_for(int site_old, const int *ops, int op_count, double *rows, int *dim_out)
{
    int gid = little_gid_for_old(site_old);
    if (!gid) return PE_KEY;
    int dim = (int)at(&T.dim, gid - 1), sg = (int)at(&T.sg, gid - 1);
    if (dim <= 0 || dim > 48) return PE_TABLES;
    int64_t generated[GEO_MAX_OPS][5];
    int generated_count = geo_sg_records(sg, generated);
    double *matrix = malloc(sizeof(double) * (size_t)dim * (size_t)dim);
    if (!matrix) return PE_TABLES;
    for (int i = 0; i < op_count; i++) {
        int position = 0;
        for (int g = 0; g < generated_count; g++)
            if (generated[g][4] == ops[i]) position = g + 1;
        if (!position || little_sparse_matrix(gid, dim, position, matrix)) {
            free(matrix);
            return PE_KEY;
        }
        memcpy(rows + (size_t)i * dim, matrix, sizeof(double) * (size_t)dim);
    }
    free(matrix);
    *dim_out = dim;
    return 0;
}

/* out: project count, site old id, vector dim, prep status, vector basis id,
 * parent op count n, n parent ops, basis status (0, -1 skipped, else an
 * error), project count of selected block indices, active count a, then
 * per bridge family: family, bridge project count, vector loop count,
 * source active count m; active_index/active_value hold the a basis values
 * followed by every family's m source values.  Returns 0, or the status of
 * the project count (raised by the caller). */
int projection_entry(int sg, int row_id, int site_pg, int gid, int pg_irrep, int param_count, const int64_t *params,
                     const double *inverse, int64_t capacity, int64_t *out, int64_t *active_index,
                     double *active_value)
{
    bad = 0;
    if (load()) return PE_TABLES;
    geo_take_error();
    int start, row_count;
    geo_wyckoff_row_range(sg, &start, &row_count);
    int64_t project_count = 0;
    int status = project_count_of(gid, row_id, site_pg, pg_irrep, &project_count);
    if (status) return status;
    int site_old = site_old_id(site_pg, pg_irrep);
    if (!site_old) return PE_KEY;
    int64_t *p = out;
    *p++ = project_count;
    *p++ = site_old;
    *p++ = at(&T.image_dimension, image_of(site_old) - 1);
    if (bad) return PE_TABLES;
    static int last_sg, last_row, last_pg, last_count;
    static int64_t last_records[GEO_MAX_OPS][5];
    int64_t records[GEO_MAX_OPS][5];
    int record_count, site_error = 0;
    if (last_sg == sg && last_row == row_id && last_pg == site_pg) {
        record_count = last_count;
        memcpy(records, last_records, sizeof(int64_t) * 5 * (size_t)record_count);
    } else {
        record_count = geo_site_records(sg, row_id, site_pg, records);
        site_error = geo_take_error();
        if (!site_error) {
            last_sg = sg;
            last_row = row_id;
            last_pg = site_pg;
            last_count = record_count;
            memcpy(last_records, records, sizeof(int64_t) * 5 * (size_t)record_count);
        }
    }
    int64_t basis_id = 0;
    int prep_status = site_error ? PE_TABLES : vector_basis_id(site_pg, (const int64_t (*)[5])records, record_count, &basis_id);
    *p++ = prep_status;
    *p++ = basis_id;
    *p++ = prep_status ? 0 : record_count;
    for (int i = 0; i < record_count && !prep_status; i++) *p++ = records[i][4];
    int type = (int)at(&T.type, gid - 1), dim = (int)at(&T.dim, gid - 1);
    int64_t *basis_status = p++;
    *basis_status = 0;
    if (type < 1 || type > 3) {
        *basis_status = -1;
        return 0;
    }
    int records32[GEO_MAX_OPS * 5], site_ops[GEO_MAX_OPS], site_count = 0, site_dim = 0;
    for (int i = 0; i < record_count; i++)
        for (int j = 0; j < 5; j++) records32[5 * i + j] = (int)records[i][j];
    double *parent = malloc(sizeof(double) * (size_t)(record_count + 1) * (size_t)dim * (size_t)dim);
    double *rows = malloc(sizeof(double) * (size_t)(record_count + 1) * 48);
    int64_t basis_length = project_count * FAMILY_STRIDE > 720 ? project_count * FAMILY_STRIDE : 720;
    int pd = dim, active_capacity = (int)(project_count * pd * 48 + 1);
    int *selected = malloc(sizeof(int) * (size_t)(project_count + 1));
    int64_t *indices = malloc(sizeof(int64_t) * (size_t)active_capacity);
    double *values = malloc(sizeof(double) * (size_t)active_capacity);
    int active_count = 0;
    if (!parent || !rows || !selected || !indices || !values) status = PE_TABLES;
    if (!status && case_irrep_stack(gid, param_count, params, record_count, records32,
                                    (int64_t)record_count * dim * dim, parent))
        status = PE_KEY;
    if (!status) status = site_ops_for(site_pg, (const int64_t (*)[5])records, record_count, site_ops, &site_count);
    if (!status) status = site_rows_for(site_old, site_ops, site_count, rows, &site_dim);
    if (!status && site_count != record_count) status = PE_VALUE;
    if (!status) {
        int block_count = 0;
        if (type == 3) {
            int k = (int)at(&T.kslot, gid - 1);
            block_count = (int)at(&T.star_count, (int64_t)(sg - 1) * 27 + k - 1);
        }
        int chosen = project_selection(type, pd, site_dim, record_count, parent, rows, type == 1 ? 0 : inverse,
                                       block_count, (int)project_count, (int)basis_length, active_capacity, selected,
                                       indices, values, &active_count);
        if (chosen) status = chosen == 1 ? PE_VALUE : PE_TABLES;
    }
    if (!status) {
        if (p + project_count + 1 + 4 * (project_count + 1) > out + capacity || active_count > capacity) status = PE_CAPACITY;
    }
    if (status) {
        *basis_status = status;
    } else {
        for (int64_t i = 0; i < project_count; i++) *p++ = selected[i];
        *p++ = active_count;
        int64_t written = 0;
        for (int i = 0; i < active_count; i++) {
            active_index[written] = indices[i];
            active_value[written++] = values[i];
        }
        int64_t *family_count = p++;
        *family_count = 0;
        int64_t families = 0;
        for (int i = 0; i < active_count; i++)
            if (indices[i] / FAMILY_STRIDE + 1 > families) families = indices[i] / FAMILY_STRIDE + 1;
        for (int64_t family = 0; family < families; family++) {
            int64_t top = -1, loops = 0, members = 0;
            for (int i = 0; i < active_count; i++) {
                if (indices[i] / FAMILY_STRIDE != family) continue;
                int64_t local = indices[i] % FAMILY_STRIDE;
                members++;
                if (local / 3 > top) top = local / 3;
                if (local % 3 + 1 > loops) loops = local % 3 + 1;
            }
            if (!members) continue;
            double source[SOURCE_LENGTH] = {0};
            int64_t offset = 64 - 3 * (top + 1);
            for (int i = 0; i < active_count; i++) {
                if (indices[i] / FAMILY_STRIDE != family) continue;
                int64_t target = offset + indices[i] % FAMILY_STRIDE;
                if (target >= 0 && target < SOURCE_LENGTH) source[target] = values[i];
            }
            *p++ = family;
            *p++ = top;
            *p++ = loops;
            int64_t *member_count = p++;
            *member_count = 0;
            for (int t = 0; t < SOURCE_LENGTH; t++)
                if (source[t] > 1e-12 || source[t] < -1e-12) {
                    if (written >= capacity) {
                        status = PE_CAPACITY;
                        break;
                    }
                    active_index[written] = t;
                    active_value[written++] = source[t];
                    (*member_count)++;
                }
            (*family_count)++;
        }
    }
    free(parent);
    free(rows);
    free(selected);
    free(indices);
    free(values);
    return status == PE_CAPACITY ? status : 0;
}

/* The vector modes (polar displacements, or axial moments when axial) a
 * little irrep carries over the given Wyckoff rows: per row, the project
 * counts of the site's vector irreps over their dimensions, times the full
 * dimension, or nothing unless that is whole. */
int vector_mode_total(int sg, int gid, int row_count, const int *row_ids, int axial, int64_t *total)
{
    bad = 0;
    if (load()) return PE_TABLES;
    geo_take_error();
    int start, count;
    if (geo_wyckoff_row_range(sg, &start, &count)) return PE_TABLES;
    int old = (int)at(&T.old, gid - 1), type = (int)at(&T.type, gid - 1), dim = (int)at(&T.dim, gid - 1);
    int divisor = type == 1 ? 1 : type == 2 ? 4 : type == 3 ? 2 : 0;
    *total = 0;
    for (int r = 0; r < row_count; r++) {
        int offset = row_ids[r] - start;
        if (offset < 0 || offset >= count) return PE_KEY;
        if (!divisor) continue;
        int site_pg = geo_wyckoff_site_pg(row_ids[r]);
        int64_t frequency[12];
        int pairs = 0;
        if (row_frequencies(gid, row_ids[r], frequency, &pairs)) continue;
        int64_t counts[13] = {0}, dims[13] = {0};
        int present[13] = {0}, whole = 1;
        for (int pg = 1; pg <= pairs && whole; pg++) {
            if (!frequency[pg - 1]) continue;
            int site_old = site_old_id(site_pg, pg);
            if (!site_old) return PE_KEY;
            int64_t image = image_of(site_old);
            dims[pg] = at(&T.image_dimension, image - 1);
            int64_t image_type = image - 1 < section_count(&T.image_type) ? at(&T.image_type, image - 1) : 0;
            int64_t numerator = frequency[pg - 1] * (image_type == 1 || image_type == 0 ? 1 : 2);
            if (numerator % divisor) whole = 0;
            counts[pg] = numerator / divisor;
            present[pg] = 1;
        }
        if (!whole) continue;
        rat sum = rat_int(0);
        rat_overflow = 0;
        for (int i = 0; i < 3; i++) {
            int pg = (int)at(&T.vector_reps, (int64_t)(site_pg - 1) * 6 + (axial ? 3 : 0) + i);
            if (pg >= 1 && pg <= 12 && present[pg] && dims[pg]) sum = rat_add(sum, rat_make(counts[pg], dims[pg]));
        }
        sum = rat_mul(sum, rat_int(dim));
        if (rat_overflow) return PE_VALUE;
        if (sum.d == 1) *total += (int64_t)sum.n;
    }
    if (bad || geo_take_error()) return PE_TABLES;
    return 0;
}
