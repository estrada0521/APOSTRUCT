/* The complete ordinary child layout of one parent site: the Formula15
 * partition of its child atoms, each orbit's Wyckoff label in the displayed
 * child setting, and the binding of every canonical atom to one row of the
 * Source display atom table, in presentation order. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "kernels.h"
#include "pynum.h"

enum {
    CL_MEMORY = 1,
    CL_RANGE = 2,
    CL_FORMULA = 3,
    CL_ORDINALS = 4,
    CL_PRODUCT = 5,
    CL_FLAT_ORDER = 6,
    CL_COSET = 7,
    CL_MULTIPLICITY = 8,
    CL_DUPLICATE_ATOM = 9,
    CL_BIJECTION = 10,
    CL_NONFINITE = 11,
    CL_RECORD = 12,
    CL_CAPACITY = 13,
    CL_EMPTY_TABLE = 14,
    CL_SHARED_POSITION = 15,
    CL_REGROUP = 16,
    CL_REPRESENTATIVE = 17,
    CL_MAGNETIC = 6000,
    CL_CHILD_SITES = 1000,
    CL_PUBLIC = 2000,
    CL_SITE_TABLE = 3000,
    CL_PRESENT = 4000,
    CL_ASSIGN = 5000,
};

#define FORMULA_STRIDE 37

typedef struct {
    int64_t raw, centering, flat, kernel[4], coset[5];
    int64_t parent[12], child[12];
    int site, row;
    double position[3];
} atom_t;

typedef struct {
    int64_t offset, multiplicity, representative, first, count;
} site_t;

typedef struct {
    int64_t label, multiplicity, order;
    int representative, count;
    int *members;
} group_t;

static const int *row_of_atom;

static int by_row(const void *a, const void *b)
{
    int left = row_of_atom[*(const int *)a], right = row_of_atom[*(const int *)b];
    return (left > right) - (left < right);
}

static int by_order(const void *a, const void *b)
{
    const group_t *left = a, *right = b;
    return (left->order > right->order) - (left->order < right->order);
}

/* child_atom_layout_in_presentation_order: groups by their first display
 * row, members by display row. */
static void presentation_order(group_t *groups, int count, const int *row_by_atom)
{
    row_of_atom = row_by_atom;
    for (int g = 0; g < count; g++) {
        qsort(groups[g].members, (size_t)groups[g].count, sizeof(int), by_row);
        groups[g].order = groups[g].count ? row_by_atom[groups[g].members[0]] : INT64_MAX;
    }
    qsort(groups, (size_t)count, sizeof *groups, by_order);
}

typedef struct {
    int capacity;
    site_t *sites;
    atom_t *atoms;
    int64_t *formula_sites, *formula_atoms, *points, *records, *order, *assign_atoms, *assign_rows, *xyz_words;
    double *xyz, *dxyz, *out_xyz, *out_dxyz;
    int *source, *centering, *owners, *atom_index;
} buffers;

static void release(buffers *b)
{
    free(b->sites);
    free(b->atoms);
    free(b->formula_sites);
    free(b->formula_atoms);
    free(b->points);
    free(b->records);
    free(b->order);
    free(b->assign_atoms);
    free(b->assign_rows);
    free(b->xyz_words);
    free(b->xyz);
    free(b->dxyz);
    free(b->out_xyz);
    free(b->out_dxyz);
    free(b->source);
    free(b->centering);
    free(b->owners);
    free(b->atom_index);
}

/* label_memo: [count, then (child row offset, multiplicity, public offset)
 * triples] shared by the sites of one child structure.  A positive
 * magnetic_group regroups the ordinary sites by its exact action.  Outputs
 * in final order: sites_out (public child row offset, or for a magnetic
 * group the 1-based reference row ordinal; multiplicity; representative
 * atom; first atom; atom count), atoms_out (raw index, centering, flat index, kernel[4],
 * coset[5], member order, parent xyz words[12], structure xyz words[12]),
 * positions_out (3 doubles), rows_out (atom, raw index, centering,
 * record[5]); counts: sites, atoms, rows. */
int child_site_layout(int parent_sg, int child_sg, int parent_row_id, const int64_t *split_basis,
                      const int64_t *split_origin, const int64_t *parent_param_words, int setting_count,
                      const int *setting_ids, const int64_t *matrix_words, const int64_t *origin_words,
                      int64_t *label_memo, int memo_capacity, int source_row_id, int source_site_pg,
                      const int64_t *site_param_words, const int64_t *presentation_pml,
                      const int64_t *representative_words, const int64_t *bridge_words,
                      const double *presentation_basis, const double *presentation_origin, int centering_letter,
                      int point_operation_count, int magnetic_group, int capacity, int64_t *sites_out,
                      int64_t *atoms_out,
                      double *positions_out, int64_t *rows_out, int *counts)
{
    if (geo_load()) return CL_CHILD_SITES + 1;
    rat_overflow = 0;
    buffers b = {0};
    int status = 0;
    b.capacity = capacity;
    b.sites = calloc((size_t)capacity, sizeof *b.sites);
    b.atoms = calloc((size_t)capacity, sizeof *b.atoms);
    b.formula_sites = calloc((size_t)capacity * 5, sizeof(int64_t));
    b.formula_atoms = calloc((size_t)capacity * FORMULA_STRIDE, sizeof(int64_t));
    if (!b.sites || !b.atoms || !b.formula_sites || !b.formula_atoms) {
        release(&b);
        return CL_MEMORY;
    }
    int site_count = 0, atom_count = 0;
    status = wyckoff_child_sites(parent_sg, child_sg, parent_row_id, split_basis, split_origin, parent_param_words,
                                 capacity, capacity, b.formula_sites, b.formula_atoms, &site_count, &atom_count);
    if (status) {
        release(&b);
        return CL_CHILD_SITES + status;
    }
    if (!site_count) status = CL_FORMULA;
    int kernel_count = 0, parent_count = 0;
    for (int a = 0; a < atom_count && !status; a++) {
        const int64_t *f = b.formula_atoms + (size_t)a * FORMULA_STRIDE;
        atom_t *t = &b.atoms[a];
        t->flat = f[0];
        memcpy(t->kernel, f + 2, sizeof t->kernel);
        memcpy(t->coset, f + 7, sizeof t->coset);
        t->centering = f[12];
        memcpy(t->parent, f + 13, sizeof t->parent);
        memcpy(t->child, f + 25, sizeof t->child);
        if (f[0] < 0 || f[1] < 0 || f[6] < 0 || f[12] < 0) status = CL_FORMULA;
        if (t->kernel[3] <= 0 || t->coset[3] <= 0 || t->coset[4] <= 0) status = CL_FORMULA;
        if (!f[12]) {
            if (f[1] + 1 > kernel_count) kernel_count = (int)f[1] + 1;
            if (f[6] + 1 > parent_count) parent_count = (int)f[6] + 1;
        }
    }
    int pairs = kernel_count * parent_count;
    char *covered = calloc((size_t)pairs + 1, 1);
    int64_t (*coset_of)[5] = calloc((size_t)parent_count + 1, sizeof *coset_of);
    char *coset_seen = calloc((size_t)parent_count + 1, 1);
    if (!covered || !coset_of || !coset_seen) status = CL_MEMORY;
    if (!status && !pairs) status = CL_ORDINALS;
    for (int a = 0; a < atom_count && !status; a++) {
        const int64_t *f = b.formula_atoms + (size_t)a * FORMULA_STRIDE;
        if (f[12]) continue;
        int k = (int)f[1], p = (int)f[6];
        covered[k * parent_count + p] = 1;
        if (f[0] != (int64_t)k * parent_count + p) status = CL_FLAT_ORDER;
        if (!coset_seen[p]) {
            coset_seen[p] = 1;
            memcpy(coset_of[p], f + 7, sizeof coset_of[p]);
        } else if (memcmp(coset_of[p], f + 7, sizeof coset_of[p]))
            status = CL_COSET;
    }
    for (int i = 0; i < pairs && !status; i++)
        if (!covered[i]) status = CL_PRODUCT;
    free(covered);
    free(coset_of);
    free(coset_seen);
    int start, row_count;
    geo_wyckoff_row_range(child_sg, &start, &row_count);
    for (int s = 0; s < site_count && !status; s++) {
        const int64_t *fs = b.formula_sites + (size_t)s * 5;
        site_t *site = &b.sites[s];
        site->offset = fs[0];
        site->multiplicity = fs[1];
        site->first = fs[2];
        site->count = fs[3];
        site->representative = fs[2] + fs[4];
        if (fs[4] < 0 || fs[4] >= fs[3] || fs[1] != fs[3] || fs[0] < 1 || fs[0] > row_count) {
            status = fs[1] != fs[3] ? CL_MULTIPLICITY : CL_FORMULA;
            break;
        }
        for (int64_t a = fs[2]; a < fs[2] + fs[3]; a++) {
            const int64_t *f = b.formula_atoms + (size_t)a * FORMULA_STRIDE;
            b.atoms[a].raw = f[6] * kernel_count + f[1];
            b.atoms[a].site = s;
            for (int64_t other = 0; other < a; other++)
                if (b.atoms[other].raw == b.atoms[a].raw && b.atoms[other].centering == b.atoms[a].centering)
                    status = CL_DUPLICATE_ATOM;
        }
        if (status) break;
        int64_t public_offset = -1;
        for (int64_t m = 0; m < label_memo[0]; m++)
            if (label_memo[1 + 3 * m] == fs[0] && label_memo[2 + 3 * m] == fs[1]) public_offset = label_memo[3 + 3 * m];
        if (public_offset < 0) {
            int64_t *words = malloc(sizeof(int64_t) * 12 * (size_t)fs[3]);
            if (!words) {
                status = CL_MEMORY;
                break;
            }
            for (int64_t a = 0; a < fs[3]; a++) memcpy(words + 12 * a, b.atoms[fs[2] + a].child, sizeof(int64_t) * 12);
            int found = child_public_wyckoff(child_sg, setting_count, setting_ids, matrix_words, origin_words,
                                             (int)fs[3], words, geo_wyckoff_site_pg(start + (int)fs[0] - 1));
            free(words);
            if (found < 0) {
                status = CL_PUBLIC - found;
                break;
            }
            if (label_memo[0] == memo_capacity) {
                status = CL_CAPACITY;
                break;
            }
            int64_t m = label_memo[0]++;
            label_memo[1 + 3 * m] = fs[0];
            label_memo[2 + 3 * m] = fs[1];
            label_memo[3 + 3 * m] = found;
            public_offset = found;
        }
        site->offset = public_offset;
    }

    int64_t table_capacity = capacity;
    b.points = malloc(sizeof(int64_t) * 12 * (size_t)table_capacity);
    b.records = malloc(sizeof(int64_t) * 5 * (size_t)table_capacity);
    b.order = malloc(sizeof(int64_t) * (size_t)table_capacity);
    if (!status && (!b.points || !b.records || !b.order)) status = CL_MEMORY;
    int64_t table_count = 0;
    if (!status) {
        table_count = site_atom_table(parent_sg, source_row_id, source_site_pg, site_param_words, presentation_pml,
                                      representative_words, table_capacity, b.points, b.records, b.order);
        if (table_count < 0) status = CL_SITE_TABLE - (int)table_count;
        else if (!table_count) status = CL_EMPTY_TABLE;
    }
    int presented = 0;
    if (!status) {
        b.xyz = malloc(sizeof(double) * 3 * (size_t)table_count);
        b.dxyz = calloc((size_t)table_count * 3, sizeof(double));
        b.out_xyz = malloc(sizeof(double) * 12 * (size_t)table_count);
        b.out_dxyz = malloc(sizeof(double) * 12 * (size_t)table_count);
        b.source = malloc(sizeof(int) * 4 * (size_t)table_count);
        b.centering = malloc(sizeof(int) * 4 * (size_t)table_count);
        if (!b.xyz || !b.dxyz || !b.out_xyz || !b.out_dxyz || !b.source || !b.centering) status = CL_MEMORY;
    }
    if (!status) {
        rat bridge[3][3], shift[3];
        if (bridge_words) {
            for (int row = 0; row < 3; row++) read_rat3(bridge_words + 12 * row, bridge[row]);
            read_rat3(bridge_words + 36, shift);
        }
        for (int64_t i = 0; i < table_count; i++) {
            rat point[3], shown[3];
            read_rat3(b.points + 12 * b.order[i], point);
            if (bridge_words) {
                rat_row_mul(point, (const rat (*)[3])bridge, shown);
                for (int axis = 0; axis < 3; axis++) shown[axis] = rat_add(shown[axis], shift[axis]);
            } else
                memcpy(shown, point, sizeof shown);
            for (int axis = 0; axis < 3; axis++) b.xyz[3 * i + axis] = rat_to_double(shown[axis]);
        }
        int present_status = present_rows((int)table_count, b.xyz, b.dxyz, presentation_basis, presentation_origin,
                                          centering_letter, 1e-8, 4 * (int)table_count, b.source, b.centering,
                                          b.out_xyz, b.out_dxyz, &presented);
        if (present_status) status = CL_PRESENT + present_status;
    }
    if (!status) {
        b.assign_atoms = malloc(sizeof(int64_t) * 10 * (size_t)(atom_count + 1));
        b.assign_rows = malloc(sizeof(int64_t) * 6 * (size_t)(presented + 1));
        b.owners = malloc(sizeof(int) * (size_t)(presented + 1));
        b.atom_index = malloc(sizeof(int) * (size_t)(atom_count + 1));
        if (!b.assign_atoms || !b.assign_rows || !b.owners || !b.atom_index) status = CL_MEMORY;
    }
    if (!status) {
        for (int a = 0; a < atom_count; a++) {
            memcpy(b.assign_atoms + 10 * a, b.atoms[a].kernel, sizeof(int64_t) * 4);
            memcpy(b.assign_atoms + 10 * a + 4, b.atoms[a].coset, sizeof(int64_t) * 5);
            b.assign_atoms[10 * a + 9] = b.atoms[a].centering;
            b.atoms[a].row = -1;
        }
        for (int r = 0; r < presented; r++) {
            memcpy(b.assign_rows + 6 * r, b.records + 5 * b.order[b.source[r]], sizeof(int64_t) * 5);
            b.assign_rows[6 * r + 5] = b.centering[r];
        }
        int assign_status = row_assignments(child_sg, presentation_pml, point_operation_count, atom_count,
                                            b.assign_atoms, presented, b.assign_rows, b.owners);
        if (assign_status) status = CL_ASSIGN + assign_status;
    }
    if (!status) {
        int *owned = calloc((size_t)atom_count + 1, sizeof(int));
        if (!owned) status = CL_MEMORY;
        for (int r = 0; r < presented && !status; r++) {
            int owner = b.owners[r];
            if (owner < 0 || owner >= atom_count) status = CL_BIJECTION;
            else if (owned[owner]++) status = CL_BIJECTION;
            else b.atoms[owner].row = r;
        }
        for (int a = 0; a < atom_count && !status; a++)
            if (b.atoms[a].row < 0) status = CL_BIJECTION;
        if (!status && presented != atom_count) status = CL_BIJECTION;
        free(owned);
    }
    for (int a = 0; a < atom_count && !status; a++) {
        int r = b.atoms[a].row;
        const int64_t *record = b.assign_rows + 6 * r;
        if (record[3] <= 0 || record[4] < 1 || record[4] > point_operation_count) status = CL_RECORD;
        for (int axis = 0; axis < 3; axis++) {
            double value = pymod1(b.out_xyz[3 * r + axis]);
            if (!isfinite(value)) status = CL_NONFINITE;
            b.atoms[a].position[axis] = value;
        }
    }
    rat (*structure)[3] = status ? 0 : malloc(sizeof *structure * (size_t)(atom_count + 1));
    group_t *groups = status ? 0 : calloc((size_t)site_count + 1, sizeof *groups);
    int *members = status ? 0 : malloc(sizeof(int) * (size_t)(atom_count + 1));
    int *row_by_atom = status ? 0 : malloc(sizeof(int) * (size_t)(atom_count + 1));
    if (!status && (!structure || !groups || !members || !row_by_atom)) status = CL_MEMORY;
    int group_count = site_count;
    if (!status) {
        for (int a = 0; a < atom_count; a++) {
            row_by_atom[a] = b.atoms[a].row;
            for (int axis = 0; axis < 3; axis++) structure[a][axis] = input_fraction_mod1(b.atoms[a].position[axis]);
        }
        int placed = 0;
        for (int s = 0; s < site_count; s++) {
            group_t *g = &groups[s];
            g->label = b.sites[s].offset;
            g->multiplicity = b.sites[s].count;
            g->representative = (int)b.sites[s].representative;
            g->members = members + placed;
            g->count = (int)b.sites[s].count;
            for (int i = 0; i < g->count; i++) g->members[i] = (int)b.sites[s].first + i;
            placed += g->count;
        }
        presentation_order(groups, group_count, row_by_atom);
    }
    if (!status && magnetic_group > 0) {
        int *regrouped = malloc(sizeof(int) * (size_t)(atom_count + 1));
        group_t *magnetic = calloc((size_t)atom_count + 1, sizeof *magnetic);
        int64_t *parent_words = malloc(sizeof(int64_t) * 12 * (size_t)(atom_count + 1));
        int64_t *displayed = malloc(sizeof(int64_t) * 12 * (size_t)(atom_count + 1));
        int64_t *atom_order = malloc(sizeof(int64_t) * 2 * (size_t)(atom_count + 1));
        int64_t *offsets = malloc(sizeof(int64_t) * (size_t)(group_count + 2));
        int64_t *site_atoms = malloc(sizeof(int64_t) * (size_t)(atom_count + 1));
        int64_t *components = malloc(sizeof(int64_t) * 6 * (size_t)(atom_count + 1));
        int64_t *component_sites = malloc(sizeof(int64_t) * (size_t)(atom_count + 1));
        int *ordinary_index = malloc(sizeof(int) * (size_t)(atom_count + 1));
        int *position_of = malloc(sizeof(int) * (size_t)(atom_count + 1));
        char *used_site = calloc((size_t)group_count + 1, 1);
        if (!regrouped || !magnetic || !parent_words || !displayed || !atom_order || !offsets || !site_atoms ||
            !components || !component_sites || !ordinary_index || !position_of || !used_site)
            status = CL_MEMORY;
        int n = 0;
        for (int s = 0; s < group_count && !status; s++) {
            offsets[s] = n;
            for (int i = 0; i < groups[s].count; i++) {
                int a = groups[s].members[i];
                ordinary_index[n] = a;
                position_of[a] = n;
                site_atoms[n] = n;
                memcpy(parent_words + 12 * n, b.atoms[a].parent, sizeof(int64_t) * 12);
                write_rat3(displayed + 12 * n, structure[a]);
                atom_order[2 * n] = b.atoms[a].raw;
                atom_order[2 * n + 1] = b.atoms[a].centering;
                for (int m = 0; m < n && !status; m++)
                    if (!memcmp(displayed + 12 * m, displayed + 12 * n, sizeof(int64_t) * 12)) status = CL_SHARED_POSITION;
                n++;
            }
        }
        if (!status) offsets[group_count] = n;
        int component_count = 0, labels_ok = 0;
        if (!status) {
            int mag_status = magnetic_site_groups(magnetic_group, child_sg, n, parent_words, atom_order, group_count,
                                                  offsets, site_atoms, displayed, split_basis, split_origin, n,
                                                  components, component_sites, &component_count, &labels_ok);
            if (mag_status) status = CL_MAGNETIC + mag_status;
        }
        int placed = 0, used_atoms = 0;
        for (int c = 0; c < component_count && !status; c++) {
            const int64_t *component = components + 6 * c;
            group_t *g = &magnetic[c];
            g->label = labels_ok ? component[3] : component[2];
            g->multiplicity = component[1];
            g->members = regrouped + placed;
            g->count = 0;
            if (component[5] <= 0) status = CL_REGROUP;
            for (int64_t k = component[4]; k < component[4] + component[5] && !status; k++) {
                int64_t s = component_sites[k];
                if (s < 0 || s >= group_count || used_site[s]) {
                    status = CL_REGROUP;
                    break;
                }
                used_site[s] = 1;
                for (int i = 0; i < groups[s].count; i++) g->members[g->count++] = groups[s].members[i];
            }
            if (status) break;
            int64_t representative = component[0];
            if (representative < 0 || representative >= n) {
                status = CL_REPRESENTATIVE;
                break;
            }
            g->representative = ordinary_index[representative];
            int found = 0;
            for (int i = 0; i < g->count; i++) found += g->members[i] == g->representative;
            if (found != 1) status = CL_REPRESENTATIVE;
            placed += g->count;
            used_atoms += g->count;
        }
        for (int s = 0; s < group_count && !status; s++)
            if (!used_site[s]) status = CL_REGROUP;
        if (!status && used_atoms != atom_count) status = CL_REGROUP;
        if (!status) {
            memcpy(members, regrouped, sizeof(int) * (size_t)atom_count);
            for (int c = 0; c < component_count; c++) {
                magnetic[c].members = members + (magnetic[c].members - regrouped);
                groups[c] = magnetic[c];
            }
            group_count = component_count;
            presentation_order(groups, group_count, row_by_atom);
        }
        free(regrouped);
        free(magnetic);
        free(parent_words);
        free(displayed);
        free(atom_order);
        free(offsets);
        free(site_atoms);
        free(components);
        free(component_sites);
        free(ordinary_index);
        free(position_of);
        free(used_site);
    }
    if (!status) {
        int *final_index = malloc(sizeof(int) * (size_t)(atom_count + 1));
        if (!final_index) status = CL_MEMORY;
        int placed = 0;
        for (int s = 0; s < group_count && !status; s++) {
            const group_t *g = &groups[s];
            int64_t *so = sites_out + 5 * s;
            so[0] = g->label;
            so[1] = g->multiplicity;
            so[2] = -1;
            so[3] = placed;
            so[4] = g->count;
            for (int i = 0; i < g->count; i++) {
                int a = g->members[i];
                const atom_t *t = &b.atoms[a];
                if (a == g->representative) so[2] = placed + i;
                final_index[a] = placed + i;
                int64_t *ao = atoms_out + (size_t)(placed + i) * 38;
                ao[0] = t->raw;
                ao[1] = t->centering;
                ao[2] = t->flat;
                memcpy(ao + 3, t->kernel, sizeof t->kernel);
                memcpy(ao + 7, t->coset, sizeof t->coset);
                ao[12] = i;
                memcpy(ao + 13, t->parent, sizeof t->parent);
                write_rat3(ao + 25, structure[a]);
                ao[37] = s;
                memcpy(positions_out + 3 * (placed + i), t->position, sizeof t->position);
            }
            placed += g->count;
        }
        for (int r = 0; r < presented && !status; r++) {
            int64_t *ro = rows_out + 8 * r;
            ro[0] = final_index[b.owners[r]];
            ro[1] = b.order[b.source[r]];
            ro[2] = b.centering[r];
            memcpy(ro + 3, b.assign_rows + 6 * r, sizeof(int64_t) * 5);
        }
        free(final_index);
        counts[0] = group_count;
        counts[1] = atom_count;
        counts[2] = presented;
    }
    free(structure);
    free(groups);
    free(members);
    free(row_by_atom);
    if (!status && rat_overflow) status = CL_RANGE;
    release(&b);
    return status;
}
