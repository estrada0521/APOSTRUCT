/* modes.subduction_specs: the primary/secondary mode-subduction specs of a
 * selected OPD (ordinary, static, magnetic-embedding, magnetic, coupled).
 * Selected-OPD inputs arrive decoded; specs leave as records the Python
 * boundary turns into the mode schema dicts. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "spec_core.h"
#include "tables.h"

int parent_operation_records(int parent_sg, int child_sg, const int64_t *basis, const int64_t *origin, int64_t (*out)[5]);

enum {
    SK_DYN = 1,
    SK_STATIC,
    SK_MEMB,
    SK_MAG,
    SK_CPL,
    SK_FALLBACK,
};

typedef struct {
    int kind, mode_kind;
    int old_id, gid, little_type, magnetic;
    sp_str label, k_label;
    int has_kp;
    int64_t kp[4];
    int has_carrier;
    int64_t carrier[4];
    int has_pml;
    rat pml[3];
    int has_req, req_none, req_n;
    rat req[3];
    int has_occ, nocc;
    const int *occ;
    int has_dir;
    const double *dir;
    int dim, count;
    int primary;
    int has_free, free_count;
    int obs;
    int has_anchor, anchor_gid;
    int64_t anchor_kp[4];
    int has_cand;
    int64_t cand_kp[4];
    int has_heading, heading_gid;
    sp_str heading_label;
    int64_t heading_kp[4];
    int stokes;
    const double *recomputed;
    int rdim, rcount;
    int case_none, case_n;
    rat case_k[3];
    int dcase_none, dcase_n;
    rat dcase_k[3];
    int slot_order;
    int replaced;
    int row_id, domain, domain_old, frequency;
    sp_str opd;
    int use_selected_rows;
    int static_updated;
} sp_spec;

typedef struct {
    sp_spec **items;
    int count, cap;
} sp_specs;

static int specs_push(sp_specs *specs, sp_spec *spec)
{
    if (specs->count == specs->cap) {
        int cap = specs->cap ? specs->cap * 2 : 16;
        sp_spec **grown = sp_alloc(sizeof(sp_spec *) * (size_t)cap);
        if (!grown) return 1;
        if (specs->count) memcpy(grown, specs->items, sizeof(sp_spec *) * (size_t)specs->count);
        specs->items = grown;
        specs->cap = cap;
    }
    specs->items[specs->count++] = spec;
    return 0;
}

static sp_spec *spec_new(void)
{
    sp_spec *spec = sp_alloc(sizeof(sp_spec));
    if (spec) {
        memset(spec, 0, sizeof *spec);
        spec->slot_order = -1;
        spec->case_none = spec->dcase_none = 1;
    }
    return spec;
}

static sp_spec *spec_copy(const sp_spec *spec)
{
    sp_spec *out = sp_alloc(sizeof(sp_spec));
    if (out) *out = *spec;
    return out;
}

static int caught(int mask)
{
    if (!sp_err.kind || !(mask & (1 << sp_err.kind))) return 0;
    sp_err.kind = PY_NONE;
    sp_err.message[0] = 0;
    return 1;
}

#define BIT(kind) (1 << (kind))

static struct {
    section_view mag2nonmag, sub_pointer, sub_irrep, sub_frequency, sub_subgroup, sub_domain, sub_domain_old, irrep_sg,
        irrep_label, op_label, freeparam, irr_sg, old;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (geo_load() || section_find("magnetic.mag_point_op_mag2nonmag", &T.mag2nonmag) ||
        section_find("isotropy.isotropy_subduce_pointer", &T.sub_pointer) ||
        section_find("isotropy.isotropy_subduce_irrep", &T.sub_irrep) ||
        section_find("isotropy.isotropy_subduce_frequency", &T.sub_frequency) ||
        section_find("isotropy.isotropy_subduce_subgroup", &T.sub_subgroup) ||
        section_find("isotropy.isotropy_subduce_domain", &T.sub_domain) ||
        section_find("isotropy.isotropy_subduce_domain_old", &T.sub_domain_old) ||
        section_find("irreps.irrep_space_group", &T.irrep_sg) || section_find("irreps.irrep_label", &T.irrep_label) ||
        section_find("isotropy.isotropy_orderparam_label", &T.op_label) ||
        section_find("isotropy.isotropy_orderparam_freeparam", &T.freeparam) ||
        section_find("little.little_irr_space_group", &T.irr_sg) || section_find("little.little_irr_old", &T.old))
        return sp_fail(PY_VALUE, "spec stage Source tables are missing");
    T.ready = 1;
    return 0;
}

static sp_str stripped(const section_view *view, int64_t index)
{
    sp_str out;
    if (section_string(view, index, &out.text, &out.length)) return (sp_str){"", 0};
    while (out.length && (out.text[0] == ' ' || out.text[0] == '\t')) {
        out.text++;
        out.length--;
    }
    return out;
}

/* ------------------------------------------------------------ helpers */

static int occ_eq(int a, int b)
{
    if (a == b) return 1;
    const sp_occ *x = &sp_occ_pool[a], *y = &sp_occ_pool[b];
    return x->mesh == y->mesh && sp_str_eq(x->k_label, y->k_label) && sp_kp_exact_eq(x->kp, y->kp) && sp_pml_eq(x->pml, y->pml) &&
           !memcmp(x->bucket, y->bucket, sizeof x->bucket);
}

static int occs_eq(const int *a, int na, const int *b, int nb)
{
    if (na != nb) return 0;
    for (int i = 0; i < na; i++)
        if (!occ_eq(a[i], b[i])) return 0;
    return 1;
}

static int occ_in(int o, const int *set, int n)
{
    for (int i = 0; i < n; i++)
        if (occ_eq(o, set[i])) return 1;
    return 0;
}

static int same_source_kparam(const int64_t *a, const int64_t *b) { return sp_kp_periodic_eq(a, b); }

static int dir_eq(const double *a, int adim, int acount, const double *b, int bdim, int bcount)
{
    if (adim != bdim) return 0;
    if (adim && acount != bcount) return 0;
    for (int i = 0; i < adim * acount; i++)
        if (!(a[i] == b[i])) return 0;
    return 1;
}

static int pml_cmp(const rat *a, const rat *b)
{
    for (int axis = 0; axis < 3; axis++) {
        int c = rat_cmp(a[axis], b[axis]);
        if (c) return c;
    }
    return 0;
}

static int params_cmp(const rat *a, int na, const rat *b, int nb)
{
    int n = na < nb ? na : nb;
    for (int i = 0; i < n; i++) {
        int c = rat_cmp(a[i], b[i]);
        if (c) return c;
    }
    return na < nb ? -1 : na > nb;
}

static int request_params(int gid, const rat *pml, int *none, int *n, rat *out)
{
    int count = case_kparams_of_pml(gid, pml, out);
    if (count < -1) return sp_fail(count == -2 ? PY_VALUE : PY_ARITH, "k-vector Source tables are inconsistent");
    *none = count == -1;
    *n = count < 0 ? 0 : count;
    return 0;
}

static int listed(const sp_rows *rows, int gid, const int64_t *kp)
{
    for (int i = 0; i < rows->count; i++)
        if (rows->items[i].gid == gid && sp_kp_identity_eq(rows->items[i].kp, kp)) return 1;
    return 0;
}

static int rows_copy(const sp_rows *in, sp_rows *out)
{
    out->items = 0;
    out->count = out->cap = 0;
    for (int i = 0; i < in->count; i++) SP_TRY(sp_rows_push(out, &in->items[i]));
    return 0;
}

/* ------------------------------------------------------------ spec builders */

typedef struct builder builder;
typedef int (*build_fn)(builder *b, const sp_row *row, sp_spec **out);

struct builder {
    build_fn fn;
    int primary_gid, has_primary_kp;
    int64_t primary_kp[4];
    /* magnetic */
    int selected_gid;
    int64_t selected_kp[4];
    const double *sel_dir;
    int sel_dim, sel_count;
    const sp_rows *rows;
    /* coupled */
    struct cpl *cpl;
    int mode_kind;
};

static void spec_from_row(sp_spec *spec, const sp_row *row)
{
    spec->gid = row->gid;
    spec->label = row->irrep_label;
    spec->k_label = row->k_label;
    spec->has_kp = 1;
    memcpy(spec->kp, row->kp, sizeof spec->kp);
    spec->has_pml = 1;
    memcpy(spec->pml, row->pml, sizeof spec->pml);
    spec->has_occ = 1;
    spec->occ = row->occ;
    spec->nocc = row->nocc;
    spec->has_dir = 1;
    spec->dir = row->dir;
    spec->dim = row->dim;
    spec->count = row->count;
    spec->frequency = 1;
}

static int dyn_spec(const sp_row *row, int old_id, int primary, sp_spec **out)
{
    sp_spec *spec = spec_new();
    if (!spec) return 1;
    spec->kind = SK_DYN;
    spec->mode_kind = 1;
    spec_from_row(spec, row);
    spec->old_id = old_id;
    spec->has_carrier = row->has_carrier;
    memcpy(spec->carrier, row->carrier, sizeof spec->carrier);
    spec->has_req = 1;
    SP_TRY(request_params(row->gid, row->pml, &spec->req_none, &spec->req_n, spec->req));
    spec->primary = primary;
    *out = spec;
    return 0;
}

static int build_dyn(builder *b, const sp_row *row, sp_spec **out)
{
    int primary = row->gid == b->primary_gid && (!b->has_primary_kp || same_source_kparam(row->kp, b->primary_kp));
    return dyn_spec(row, 0, primary, out);
}

static int build_memb(builder *b, const sp_row *row, sp_spec **out)
{
    (void)b;
    sp_spec *spec = spec_new();
    if (!spec) return 1;
    spec->kind = SK_MEMB;
    spec->mode_kind = 1;
    spec_from_row(spec, row);
    spec->old_id = sp_little_old_id(row->gid);
    spec->has_req = 1;
    SP_TRY(request_params(row->gid, row->pml, &spec->req_none, &spec->req_n, spec->req));
    spec->primary = 0;
    *out = spec;
    return 0;
}

static int duplicates_sibling(builder *b, const sp_row *row, int index, int *out)
{
    *out = 0;
    if (!row->nocc) return 0;
    int selected_rank;
    SP_TRY(sp_matrix_rank(b->sel_dir, b->sel_dim, b->sel_count, 1e-10, &selected_rank));
    for (int i = 0; i < b->rows->count; i++) {
        const sp_row *other = &b->rows->items[i];
        if (i == index || !occs_eq(other->occ, other->nocc, row->occ, row->nocc)) continue;
        if (!other->dim || other->dim != b->sel_dim || other->count != b->sel_count) continue;
        int width = b->sel_count + other->count;
        double *stacked = sp_alloc(sizeof(double) * (size_t)(b->sel_dim * width + 1));
        if (!stacked) return 1;
        for (int r = 0; r < b->sel_dim; r++) {
            for (int c = 0; c < b->sel_count; c++) stacked[r * width + c] = b->sel_dir[r * b->sel_count + c];
            for (int c = 0; c < other->count; c++) stacked[r * width + b->sel_count + c] = other->dir[r * other->count + c];
        }
        int rank;
        SP_TRY(sp_matrix_rank(stacked, b->sel_dim, width, 1e-10, &rank));
        if (rank == selected_rank) {
            *out = 1;
            return 0;
        }
    }
    return 0;
}

static int mag_spec(builder *b, const sp_row *row, int index, sp_spec **out)
{
    sp_spec *spec = spec_new();
    if (!spec) return 1;
    spec->kind = SK_MAG;
    spec->mode_kind = 2;
    spec_from_row(spec, row);
    spec->old_id = sp_little_old_id(row->gid);
    spec->magnetic = 1;
    spec->has_carrier = row->has_carrier;
    memcpy(spec->carrier, row->carrier, sizeof spec->carrier);
    spec->has_req = 1;
    SP_TRY(request_params(row->gid, row->pml, &spec->req_none, &spec->req_n, spec->req));
    int primary = row->gid == b->selected_gid && sp_kp_identity_eq(row->kp, b->selected_kp);
    int duplicate = 0;
    if (primary) SP_TRY(duplicates_sibling(b, row, index, &duplicate));
    if (primary && !duplicate) {
        if (b->sel_dim != row->dim) return sp_fail(PY_VALUE, "selected Stokes direction has the wrong dimension");
        spec->dir = b->sel_dir;
        spec->dim = b->sel_dim;
        spec->count = b->sel_count;
        if (!dir_eq(b->sel_dir, b->sel_dim, b->sel_count, row->dir, row->dim, row->count)) {
            spec->stokes = 1;
            spec->recomputed = row->dir;
            spec->rdim = row->dim;
            spec->rcount = row->count;
        }
    }
    spec->has_free = 1;
    spec->free_count = spec->dim ? spec->count : 0;
    spec->primary = primary;
    *out = spec;
    return 0;
}

static int build_mag(builder *b, const sp_row *row, sp_spec **out)
{
    int index = -1;
    for (int i = 0; i < b->rows->count && index < 0; i++)
        if (&b->rows->items[i] == row) index = i;
    return mag_spec(b, row, index, out);
}

/* ------------------------------------------------------------ stored occurrence aliases */

typedef struct {
    int64_t kp[4];
    rat pml[3];
} kp_pml;

typedef struct {
    int gid;
    sp_str irrep_label, k_label;
    const int *occ;
    int nocc;
    kp_pml item;
} represented_entry;

typedef struct {
    int representative_gid;
    int64_t representative_kp[4], candidate_kp[4];
    int heading_gid;
    sp_str heading_label;
    int64_t heading_kp[4];
    sp_row row;
} candidate_t;

static int kp_pml_in(const kp_pml *set, int n, const kp_pml *item)
{
    for (int i = 0; i < n; i++)
        if (sp_kp_exact_eq(set[i].kp, item->kp) && sp_pml_eq(set[i].pml, item->pml)) return 1;
    return 0;
}

static int kp_pml_add(kp_pml *set, int *n, const kp_pml *item)
{
    if (kp_pml_in(set, *n, item)) return 0;
    set[(*n)++] = *item;
    return 1;
}

static int occ_subset(const int *a, int na, const int *b, int nb)
{
    for (int i = 0; i < na; i++)
        if (!occ_in(a[i], b, nb)) return 0;
    return 1;
}

typedef struct {
    const sp_frame *frame;
    int filter_exact;
} alias_mode;

static int aliases_for_class(const alias_mode *mode, const sp_row *row, const int64_t *requested, sp_rows *out)
{
    sp_rows all = {0};
    SP_TRY(sp_aliases(mode->frame, row, requested, &all));
    out->items = 0;
    out->count = out->cap = 0;
    for (int i = 0; i < all.count; i++)
        if (!mode->filter_exact || sp_kp_exact_eq(all.items[i].kp, requested)) SP_TRY(sp_rows_push(out, &all.items[i]));
    return 0;
}

static int heading_eq(int ga, sp_str la, const int64_t *ka, int gb, sp_str lb, const int64_t *kb)
{
    return ga == gb && sp_str_eq(la, lb) && sp_kp_exact_eq(ka, kb);
}

static int stored_alias_specs(const alias_mode *mode, const sp_rows *raw, const sp_rows *represented, builder *b, sp_specs *out)
{
    int total = 0;
    for (int i = 0; i < represented->count; i++) total++;
    represented_entry *entries = sp_alloc(sizeof(represented_entry) * (size_t)(total + 1));
    if (!entries) return 1;
    int nentries = 0;
    for (int i = 0; i < represented->count; i++) {
        const sp_row *row = &represented->items[i];
        if (!row->kp[3]) continue;
        represented_entry *e = &entries[nentries++];
        e->gid = row->gid;
        e->irrep_label = row->irrep_label;
        e->k_label = row->k_label;
        e->occ = row->occ;
        e->nocc = row->nocc;
        memcpy(e->item.kp, row->kp, sizeof e->item.kp);
        memcpy(e->item.pml, row->pml, sizeof e->item.pml);
    }
    candidate_t *candidates = 0;
    int ncandidates = 0, candidate_cap = 0;
    for (int r = 0; r < raw->count; r++) {
        const sp_row *row = &raw->items[r];
        if (row->gid <= 0 || sp_little_old_id(row->gid) != 0) continue;
        kp_pml *found = sp_alloc(sizeof(kp_pml) * (size_t)(nentries + 1));
        if (!found) return 1;
        int nfound = 0;
        for (int i = 0; i < nentries; i++) {
            const represented_entry *e = &entries[i];
            if (e->gid == row->gid && sp_str_eq(e->irrep_label, row->irrep_label) && sp_str_eq(e->k_label, row->k_label) &&
                occs_eq(e->occ, e->nocc, row->occ, row->nocc))
                kp_pml_add(found, &nfound, &e->item);
        }
        if (!nfound && row->nocc)
            for (int i = 0; i < nentries; i++) {
                const represented_entry *e = &entries[i];
                if (e->gid == row->gid && sp_str_eq(e->irrep_label, row->irrep_label) && sp_str_eq(e->k_label, row->k_label) && e->nocc &&
                    occ_subset(e->occ, e->nocc, row->occ, row->nocc))
                    kp_pml_add(found, &nfound, &e->item);
            }
        kp_pml exact;
        memcpy(exact.kp, row->kp, sizeof exact.kp);
        memcpy(exact.pml, row->pml, sizeof exact.pml);
        kp_pml rep;
        int promoted;
        if (row->kp[3] && kp_pml_in(found, nfound, &exact)) {
            rep = exact;
            promoted = 0;
        } else if (nfound == 1) {
            rep = found[0];
            promoted = 1;
        } else if (!nfound) {
            return sp_fail(PY_VALUE, "Source occurrence family has no represented row");
        } else {
            return sp_fail(PY_VALUE, "Source occurrence family has multiple represented rows");
        }
        if (!rep.kp[3]) return sp_fail(PY_VALUE, "Source occurrence k-parameter denominator must be nonzero");
        int64_t (*class_kp)[4] = sp_alloc(sizeof(int64_t) * 4 * (size_t)(row->nocc + 1));
        int *class_first = sp_alloc(sizeof(int) * (size_t)(row->nocc + 1));
        int *class_of = sp_alloc(sizeof(int) * (size_t)(row->nocc + 1));
        if (!class_kp || !class_first || !class_of) return 1;
        int nclasses = 0;
        for (int i = 0; i < row->nocc; i++) {
            const int64_t *kp = sp_occ_pool[row->occ[i]].kp;
            int c = -1;
            for (int k = 0; k < nclasses && c < 0; k++)
                if (sp_kp_exact_eq(class_kp[k], kp)) c = k;
            if (c < 0) {
                c = nclasses++;
                memcpy(class_kp[c], kp, 4 * sizeof(int64_t));
                class_first[c] = row->occ[i];
            }
            class_of[i] = c;
        }
        if (nclasses <= 1) continue;
        unsigned char *represented_class = sp_alloc((size_t)nclasses + 1);
        if (!represented_class) return 1;
        memset(represented_class, 0, (size_t)nclasses + 1);
        int nsource = 0;
        for (int c = 0; c < nclasses; c++) {
            if (!class_kp[c][3]) return sp_fail(PY_VALUE, "Source occurrence k-parameter denominator must be nonzero");
            if (sp_kp_periodic_eq(class_kp[c], rep.kp)) {
                represented_class[c] = 1;
                nsource++;
            }
        }
        if (nsource > 1 || (!promoted && !nsource))
            return sp_fail(PY_VALUE, "represented k parameters must identify one Source occurrence class");
        if (promoted) {
            unsigned char *pml_class = sp_alloc((size_t)nclasses + 1);
            if (!pml_class) return 1;
            memset(pml_class, 0, (size_t)nclasses + 1);
            int npml = 0;
            for (int i = 0; i < row->nocc; i++)
                if (sp_pml_eq(sp_occ_pool[row->occ[i]].pml, rep.pml) && !pml_class[class_of[i]]) {
                    pml_class[class_of[i]] = 1;
                    npml++;
                }
            if (!npml)
                for (int i = 0; i < row->nocc; i++)
                    if (sp_pml_mod1_eq(sp_occ_pool[row->occ[i]].pml, rep.pml) && !pml_class[class_of[i]]) {
                        pml_class[class_of[i]] = 1;
                        npml++;
                    }
            if (npml != 1) return sp_fail(PY_VALUE, "represented reciprocal vector must identify one Source mesh class");
            for (int c = 0; c < nclasses; c++) represented_class[c] |= pml_class[c];
        }
        for (int c = 0; c < nclasses; c++) {
            if (represented_class[c]) continue;
            sp_rows aliases;
            if (aliases_for_class(mode, row, class_kp[c], &aliases)) {
                if (caught(BIT(PY_INDEX) | BIT(PY_KEY) | BIT(PY_TYPE) | BIT(PY_VALUE))) continue;
                return 1;
            }
            int nexact = 0, which = -1;
            for (int a = 0; a < aliases.count; a++)
                if (aliases.items[a].gid == row->gid && sp_kp_exact_eq(aliases.items[a].kp, class_kp[c])) {
                    nexact++;
                    which = a;
                }
            if (nexact != 1) continue;
            sp_row alias = aliases.items[which];
            memcpy(alias.kp, class_kp[c], sizeof alias.kp);
            memcpy(alias.pml, sp_occ_pool[class_first[c]].pml, sizeof alias.pml);
            int n = 0;
            for (int i = 0; i < row->nocc; i++) n += class_of[i] == c;
            int *occ = sp_alloc(sizeof(int) * (size_t)(n + 1));
            if (!occ) return 1;
            n = 0;
            for (int i = 0; i < row->nocc; i++)
                if (class_of[i] == c) occ[n++] = row->occ[i];
            alias.occ = occ;
            alias.nocc = n;
            if (!alias.dim || !alias.count) continue;
            int finite = 1;
            for (int i = 0; i < alias.dim * alias.count; i++) finite &= isfinite(alias.dir[i]) != 0;
            if (!finite) continue;
            int rank;
            if (sp_matrix_rank(alias.dir, alias.dim, alias.count, 1e-10, &rank)) {
                if (caught(BIT(PY_LINALG) | BIT(PY_TYPE) | BIT(PY_VALUE))) continue;
                return 1;
            }
            if (rank <= 0) continue;
            if (alias.gid <= 0 || !alias.kp[3]) continue;
            int seen = 0;
            for (int i = 0; i < represented->count && !seen; i++) {
                const sp_row *rr = &represented->items[i];
                if (rr->gid > 0 && rr->kp[3] && heading_eq(rr->gid, rr->k_label, rr->kp, alias.gid, alias.k_label, alias.kp)) seen = 1;
            }
            if (seen) continue;
            if (ncandidates == candidate_cap) {
                int cap = candidate_cap ? candidate_cap * 2 : 8;
                candidate_t *grown = sp_alloc(sizeof(candidate_t) * (size_t)cap);
                if (!grown) return 1;
                if (ncandidates) memcpy(grown, candidates, sizeof(candidate_t) * (size_t)ncandidates);
                candidates = grown;
                candidate_cap = cap;
            }
            candidate_t *cand = &candidates[ncandidates++];
            cand->representative_gid = row->gid;
            memcpy(cand->representative_kp, rep.kp, sizeof cand->representative_kp);
            memcpy(cand->candidate_kp, class_kp[c], sizeof cand->candidate_kp);
            cand->heading_gid = alias.gid;
            cand->heading_label = alias.k_label;
            memcpy(cand->heading_kp, alias.kp, sizeof cand->heading_kp);
            cand->row = alias;
        }
    }
    unsigned char *keep = sp_alloc((size_t)ncandidates + 1);
    if (!keep) return 1;
    for (int i = 0; i < ncandidates; i++) {
        int n = 0;
        for (int j = 0; j < ncandidates; j++)
            n += heading_eq(candidates[i].heading_gid, candidates[i].heading_label, candidates[i].heading_kp, candidates[j].heading_gid,
                            candidates[j].heading_label, candidates[j].heading_kp);
        keep[i] = n == 1;
    }
    for (int r = 0; r < raw->count; r++) {
        const sp_row *row = &raw->items[r];
        if (row->gid <= 0) continue;
        int required = 0;
        for (int i = 0; i < ncandidates && !required; i++)
            required = keep[i] && candidates[i].representative_gid == row->gid && sp_kp_exact_eq(candidates[i].representative_kp, row->kp);
        if (!required) continue;
        int represented_exact = 0;
        for (int i = 0; i < represented->count && !represented_exact; i++)
            represented_exact = represented->items[i].gid > 0 && represented->items[i].gid == row->gid &&
                                sp_kp_exact_eq(represented->items[i].kp, row->kp);
        if (represented_exact) continue;
        sp_spec *spec;
        SP_TRY(b->fn(b, row, &spec));
        spec->primary = 0;
        spec->obs = 1;
        SP_TRY(specs_push(out, spec));
    }
    for (int i = 0; i < ncandidates; i++) {
        if (!keep[i]) continue;
        sp_spec *spec;
        SP_TRY(b->fn(b, &candidates[i].row, &spec));
        spec->primary = 0;
        spec->has_anchor = 1;
        spec->anchor_gid = candidates[i].representative_gid;
        memcpy(spec->anchor_kp, candidates[i].representative_kp, sizeof spec->anchor_kp);
        spec->has_cand = 1;
        memcpy(spec->cand_kp, candidates[i].candidate_kp, sizeof spec->cand_kp);
        spec->has_heading = 1;
        spec->heading_gid = candidates[i].heading_gid;
        spec->heading_label = candidates[i].heading_label;
        memcpy(spec->heading_kp, candidates[i].heading_kp, sizeof spec->heading_kp);
        SP_TRY(specs_push(out, spec));
    }
    return 0;
}

static int spec_gid_kp(const sp_spec *spec, int *gid, const int64_t **kp)
{
    *gid = spec->gid > 0 ? spec->gid : 0;
    *kp = spec->has_kp && spec->kp[3] ? spec->kp : 0;
    return *gid && *kp;
}

static int merge_alias_specs(const sp_specs *represented, const sp_specs *aliases, sp_specs *out)
{
    int *placed = sp_alloc(sizeof(int) * (size_t)(aliases->count + 1));
    if (!placed) return 1;
    for (int a = 0; a < aliases->count; a++) {
        const sp_spec *alias = aliases->items[a];
        placed[a] = -1;
        if (!alias->has_anchor || alias->anchor_gid <= 0 || !alias->anchor_kp[3]) continue;
        int match = -1, n = 0;
        for (int i = 0; i < represented->count; i++) {
            int gid;
            const int64_t *kp;
            if (!spec_gid_kp(represented->items[i], &gid, &kp)) continue;
            if (gid == alias->anchor_gid && sp_kp_exact_eq(kp, alias->anchor_kp)) {
                if (!n) match = i;
                n++;
            }
        }
        if (n == 1) placed[a] = match;
    }
    sp_spec **cluster = sp_alloc(sizeof(sp_spec *) * (size_t)(aliases->count + 2));
    if (!cluster) return 1;
    for (int i = 0; i < represented->count; i++) {
        int n = 0;
        cluster[n++] = represented->items[i];
        for (int a = 0; a < aliases->count; a++)
            if (placed[a] == i) cluster[n++] = aliases->items[a];
        for (int x = 1; x < n; x++) {
            sp_spec *item = cluster[x];
            int y = x - 1;
            while (y >= 0) {
                int item_none = !item->has_pml, prev_none = !cluster[y]->has_pml;
                int greater = prev_none != item_none ? prev_none > item_none : (!item_none && pml_cmp(cluster[y]->pml, item->pml) > 0);
                if (!greater) break;
                cluster[y + 1] = cluster[y];
                y--;
            }
            cluster[y + 1] = item;
        }
        for (int k = 0; k < n; k++) SP_TRY(specs_push(out, cluster[k]));
    }
    for (int a = 0; a < aliases->count; a++)
        if (placed[a] < 0) SP_TRY(specs_push(out, aliases->items[a]));
    return 0;
}

static int canonicalize(const sp_specs *specs, sp_specs *out)
{
    int n = specs->count;
    int *replacement_of = sp_alloc(sizeof(int) * (size_t)(n + 1));
    sp_spec **replacement = sp_alloc(sizeof(sp_spec *) * (size_t)(n + 1));
    unsigned char *removed = sp_alloc((size_t)n + 1), *grouped = sp_alloc((size_t)n + 1);
    if (!replacement_of || !replacement || !removed || !grouped) return 1;
    memset(removed, 0, (size_t)n + 1);
    memset(grouped, 0, (size_t)n + 1);
    for (int i = 0; i < n; i++) replacement_of[i] = -1;
    for (int g = 0; g < n; g++) {
        const sp_spec *head = specs->items[g];
        if (grouped[g] || head->gid <= 0 || !head->has_pml) continue;
        int group_gid = head->gid;
        int cand = -1, ncand = 0, ord = -1, nord = 0;
        for (int i = g; i < n; i++) {
            const sp_spec *s = specs->items[i];
            if (grouped[i] || s->gid != group_gid || !s->has_pml || !sp_pml_eq(s->pml, head->pml)) continue;
            grouped[i] = 1;
            if (s->has_cand) {
                cand = i;
                ncand++;
            } else if (!s->obs) {
                ord = i;
                nord++;
            }
        }
        if (ncand != 1 || nord != 1) continue;
        const sp_spec *candidate = specs->items[cand], *row = specs->items[ord];
        if (!candidate->cand_kp[3] || !row->has_kp || !row->kp[3]) continue;
        const int64_t *source_kparam = candidate->cand_kp;
        int generic = candidate->has_anchor && candidate->anchor_gid > 0 && candidate->anchor_kp[3] && candidate->anchor_gid == group_gid &&
                      candidate->has_heading && candidate->heading_gid > 0 && candidate->heading_kp[3] &&
                      heading_eq(candidate->heading_gid, candidate->heading_label, candidate->heading_kp, group_gid, candidate->k_label,
                                 source_kparam);
        if (!generic || candidate->obs || candidate->primary || row->primary || !sp_str_eq(candidate->label, row->label) ||
            !sp_str_eq(candidate->k_label, row->k_label) ||
            !dir_eq(candidate->dir, candidate->dim, candidate->count, row->dir, row->dim, row->count) || !candidate->has_kp ||
            !candidate->kp[3] || !sp_kp_exact_eq(candidate->kp, source_kparam) || sp_kp_exact_eq(row->kp, source_kparam) || !row->nocc)
            continue;
        int bad = 0;
        for (int i = 0; i < row->nocc && !bad; i++) {
            const sp_occ *o = &sp_occ_pool[row->occ[i]];
            bad = !o->kp[3] || !sp_kp_exact_eq(o->kp, source_kparam) || !sp_pml_eq(o->pml, candidate->pml);
        }
        if (bad) continue;
        int anchors = 0, anchor = -1;
        for (int i = 0; i < n; i++) {
            const sp_spec *s = specs->items[i];
            int gid;
            const int64_t *kp;
            spec_gid_kp(s, &gid, &kp);
            if (s->has_cand || !gid || !kp || gid != candidate->anchor_gid || !sp_kp_exact_eq(kp, candidate->anchor_kp) || s->primary ||
                !sp_str_eq(s->label, candidate->label) || !sp_str_eq(s->k_label, candidate->k_label) ||
                !dir_eq(s->dir, s->dim, s->count, candidate->dir, candidate->dim, candidate->count))
                continue;
            anchors++;
            anchor = i;
        }
        if (anchors != 1) continue;
        const sp_spec *anchor_row = specs->items[anchor];
        int hit = 0;
        for (int i = 0; i < anchor_row->nocc && !hit; i++) {
            const sp_occ *o = &sp_occ_pool[anchor_row->occ[i]];
            hit = o->kp[3] && sp_kp_exact_eq(o->kp, source_kparam) && sp_pml_eq(o->pml, candidate->pml);
        }
        if (!hit) continue;
        sp_spec *canonical = spec_copy(candidate);
        if (!canonical) return 1;
        canonical->has_anchor = canonical->has_cand = canonical->has_heading = 0;
        canonical->obs = 0;
        canonical->occ = row->occ;
        canonical->nocc = row->nocc;
        canonical->has_occ = row->has_occ;
        replacement_of[ord] = cand;
        replacement[ord] = canonical;
    }
    for (int i = 0; i < n; i++)
        if (replacement_of[i] >= 0) removed[replacement_of[i]] = 1;
    for (int i = 0; i < n; i++) {
        const sp_spec *s = specs->items[i];
        if (replacement_of[i] >= 0) {
            SP_TRY(specs_push(out, replacement[i]));
            continue;
        }
        if (removed[i]) continue;
        if (s->obs) {
            int gid;
            const int64_t *kp;
            if (spec_gid_kp(s, &gid, &kp)) {
                int active = 0;
                for (int k = 0; k < n && !active; k++) {
                    const sp_spec *a = specs->items[k];
                    if (removed[k] || !a->has_anchor || a->anchor_gid <= 0 || !a->anchor_kp[3]) continue;
                    active = a->anchor_gid == gid && sp_kp_exact_eq(a->anchor_kp, kp);
                }
                if (!active) continue;
            }
        }
        SP_TRY(specs_push(out, specs->items[i]));
    }
    return 0;
}

/* ------------------------------------------------------------ promotion */

static int promote_dynamic(const sp_frame *frame, const sp_rows *rows, const int *gids, const int64_t (*kps)[4], int ncases,
                           int preserve, sp_rows *out)
{
    SP_TRY(rows_copy(rows, out));
    unsigned char *claimed = sp_alloc((size_t)rows->count + 1);
    int *promo_index = sp_alloc(sizeof(int) * (size_t)(ncases + 1));
    sp_row *promo_alias = sp_alloc(sizeof(sp_row) * (size_t)(ncases + 1));
    if (!claimed || !promo_index || !promo_alias) return 1;
    memset(claimed, 0, (size_t)rows->count + 1);
    int npromo = 0;
    for (int s = 0; s < ncases; s++) {
        if (listed(out, gids[s], kps[s])) continue;
        int snapshot = out->count;
        for (int r = 0; r < snapshot; r++) {
            const sp_row row = out->items[r];
            if (claimed[r] || row.gid != gids[s]) continue;
            sp_rows aliases = {0};
            SP_TRY(sp_aliases(frame, &row, 0, &aliases));
            int found = -1;
            for (int a = 0; a < aliases.count && found < 0; a++)
                if (sp_kp_identity_eq(aliases.items[a].kp, kps[s])) found = a;
            sp_row selected;
            if (found >= 0) {
                selected = aliases.items[found];
            } else {
                int64_t standard[4];
                int ok, raised;
                SP_TRY(sp_standard_kparam(frame->sg, gids[s], kps[s], standard, &ok, &raised));
                if (!ok && raised && !(BIT(raised) & (BIT(PY_KEY) | BIT(PY_TYPE) | BIT(PY_VALUE))))
                    return sp_fail(raised, sp_standard_message());
                if (!ok || !sp_kp_identity_eq(row.kp, standard)) continue;
                selected = row;
                memcpy(selected.kp, kps[s], sizeof selected.kp);
            }
            promo_index[npromo] = r;
            promo_alias[npromo++] = selected;
            out->items[r] = selected;
            claimed[r] = 1;
            break;
        }
    }
    int *family_count = sp_alloc(sizeof(int) * (size_t)(npromo + 1));
    if (!family_count) return 1;
    for (int p = 0; p < npromo; p++) {
        const sp_row *row = &rows->items[promo_index[p]];
        family_count[p] = 0;
        if (!preserve || !row->nocc) continue;
        int members = 0;
        for (int r = 0; r < rows->count; r++) {
            const sp_row *other = &rows->items[r];
            int selected_identity = 0;
            for (int s = 0; s < ncases && !selected_identity; s++)
                selected_identity = other->gid == gids[s] && sp_kp_identity_eq(other->kp, kps[s]);
            int in = (selected_identity && other->nocc && occs_eq(other->occ, other->nocc, row->occ, row->nocc));
            if (!in)
                for (int q = 0; q < npromo && !in; q++)
                    in = promo_index[q] == r && rows->items[r].nocc && occs_eq(rows->items[r].occ, rows->items[r].nocc, row->occ, row->nocc);
            members += in;
        }
        family_count[p] = members;
    }
    for (int p = 0; p < npromo; p++) {
        if (family_count[p] < 2) continue;
        const sp_row *row = &rows->items[promo_index[p]];
        sp_row replaced = promo_alias[p];
        replaced.dir = row->dir;
        replaced.dim = row->dim;
        replaced.count = row->count;
        replaced.occ = row->occ;
        replaced.nocc = row->nocc;
        replaced.has_carrier = 1;
        memcpy(replaced.carrier, row->kp, sizeof replaced.carrier);
        out->items[promo_index[p]] = replaced;
    }
    return 0;
}

static int projects_case(int gid, const sp_row *row, const rat *case_k, int ncase, int *out)
{
    *out = 0;
    if (row->gid != gid || !row->kp[3]) return 0;
    int none, n;
    rat params[3];
    if (request_params(gid, row->pml, &none, &n, params)) {
        caught(~0);
        return 0;
    }
    *out = !none && n == ncase && !params_cmp(params, n, case_k, ncase);
    return 0;
}

typedef struct {
    int gid;
    int64_t kp[4];
    int has_case, ncase;
    rat case_k[3];
} mag_case;

static int promote_magnetic(const sp_frame *frame, const sp_rows *rows, const mag_case *cases, int ncases, sp_rows *out)
{
    SP_TRY(rows_copy(rows, out));
    unsigned char *claimed = sp_alloc((size_t)rows->count + 1);
    if (!claimed) return 1;
    memset(claimed, 0, (size_t)rows->count + 1);
    for (int s = 0; s < ncases; s++) {
        const mag_case *c = &cases[s];
        if (c->gid <= 0 || !c->kp[3]) continue;
        int skip = 0;
        for (int r = 0; r < out->count && !skip; r++)
            skip = out->items[r].gid == c->gid && out->items[r].kp[3] && sp_kp_identity_eq(out->items[r].kp, c->kp);
        if (skip) continue;
        if (c->has_case)
            for (int r = 0; r < out->count && !skip; r++) {
                SP_TRY(projects_case(c->gid, &out->items[r], c->case_k, c->ncase, &skip));
                if (skip && !sp_kp_periodic_eq(out->items[r].kp, c->kp)) skip = 0;
            }
        if (skip) continue;
        int snapshot = out->count;
        for (int r = 0; r < snapshot; r++) {
            const sp_row row = out->items[r];
            if (claimed[r] || row.gid != c->gid || !row.kp[3]) continue;
            sp_rows aliases = {0};
            if (sp_aliases(frame, &row, c->kp, &aliases)) {
                if (!caught(BIT(PY_INDEX) | BIT(PY_KEY) | BIT(PY_TYPE) | BIT(PY_VALUE) | BIT(PY_ZERO))) return 1;
                aliases.count = 0;
            }
            int found = -1;
            for (int a = 0; a < aliases.count && found < 0; a++)
                if (aliases.items[a].gid == c->gid && aliases.items[a].kp[3] && sp_kp_identity_eq(aliases.items[a].kp, c->kp)) found = a;
            if (found < 0) continue;
            sp_row replaced = aliases.items[found];
            replaced.dir = row.dir;
            replaced.dim = row.dim;
            replaced.count = row.count;
            replaced.occ = row.occ;
            replaced.nocc = row.nocc;
            replaced.has_carrier = 1;
            memcpy(replaced.carrier, row.kp, sizeof replaced.carrier);
            out->items[r] = replaced;
            claimed[r] = 1;
            break;
        }
    }
    return 0;
}

static int replay_type3(const sp_frame *frame, sp_rows *rows, int gid, const int64_t *kp)
{
    if (listed(rows, gid, kp)) return 0;
    if (sp_little_type(gid) != 3) return 0;
    int cap = rows->count + 64, nseen = 0;
    rat (*seen)[3] = malloc(sizeof(rat) * 3 * (size_t)cap);
    if (!seen) return sp_fail(PY_MEMORY, "spec stage ran out of memory");
#define ADD_SEEN(v)                                                                   \
    do {                                                                              \
        int dup = 0;                                                                  \
        for (int q = 0; q < nseen && !dup; q++) dup = sp_pml_eq(seen[q], (v));        \
        if (!dup) {                                                                   \
            if (nseen == cap) {                                                       \
                cap *= 2;                                                             \
                seen = realloc(seen, sizeof(rat) * 3 * (size_t)cap);                  \
                if (!seen) return sp_fail(PY_MEMORY, "spec stage ran out of memory"); \
            }                                                                         \
            memcpy(seen[nseen++], (v), sizeof(rat) * 3);                              \
        }                                                                             \
    } while (0)
    for (int r = 0; r < rows->count; r++) {
        if (rows->items[r].gid != gid) continue;
        rat physical[3];
        for (int axis = 0; axis < 3; axis++) physical[axis] = rat_mod1(rat_make(rows->items[r].kp[axis], rows->items[r].kp[3]));
        ADD_SEEN(physical);
    }
    rat selected[3];
    for (int axis = 0; axis < 3; axis++) selected[axis] = rat_mod1(rat_make(kp[axis], kp[3]));
    ADD_SEEN(selected);
    int snapshot = rows->count;
    for (int r = 0; r < snapshot; r++) {
        if (rows->items[r].gid != gid) continue;
        sp_row row = rows->items[r];
        sp_rows aliases = {0};
        if (sp_aliases(frame, &row, 0, &aliases)) {
            free(seen);
            return 1;
        }
        for (int a = 0; a < aliases.count; a++) {
            rat physical[3];
            for (int axis = 0; axis < 3; axis++) physical[axis] = rat_mod1(rat_make(aliases.items[a].kp[axis], aliases.items[a].kp[3]));
            int dup = 0;
            for (int q = 0; q < nseen && !dup; q++) dup = sp_pml_eq(seen[q], physical);
            if (dup) continue;
            ADD_SEEN(physical);
            if (sp_rows_push(rows, &aliases.items[a])) {
                free(seen);
                return 1;
            }
        }
    }
#undef ADD_SEEN
    free(seen);
    return 0;
}

static int order_promoted_families(sp_rows *rows, const int *promoted_gids, int npromoted)
{
    if (!npromoted) return 0;
    int n = rows->count;
    unsigned char *done = sp_alloc((size_t)n + 1);
    if (!done) return 1;
    memset(done, 0, (size_t)n + 1);
    for (int start = 0; start < n; start++) {
        int is_promoted_label = 0;
        sp_str label = rows->items[start].k_label;
        if (done[start]) continue;
        for (int r = 0; r < n && !is_promoted_label; r++)
            for (int g = 0; g < npromoted && !is_promoted_label; g++)
                is_promoted_label = rows->items[r].gid == promoted_gids[g] && sp_str_eq(rows->items[r].k_label, label);
        int *indices = sp_alloc(sizeof(int) * (size_t)(n + 1));
        if (!indices) return 1;
        int m = 0;
        for (int r = 0; r < n; r++)
            if (sp_str_eq(rows->items[r].k_label, label)) {
                indices[m++] = r;
                done[r] = 1;
            }
        if (!is_promoted_label || m <= 1) continue;
        rat (*params)[3] = sp_alloc(sizeof(rat) * 3 * (size_t)(m + 1));
        int *nparams = sp_alloc(sizeof(int) * (size_t)(m + 1)), *order = sp_alloc(sizeof(int) * (size_t)(m + 1));
        sp_row *saved = sp_alloc(sizeof(sp_row) * (size_t)(m + 1));
        if (!params || !nparams || !order || !saved) return 1;
        int ok = 1;
        for (int k = 0; k < m && ok; k++) {
            const sp_row *row = &rows->items[indices[k]];
            int have = 0;
            int nvec = row->nocc ? row->nocc : 1;
            for (int v = 0; v < nvec && ok; v++) {
                const rat *vector = row->nocc ? sp_occ_pool[row->occ[v]].pml : row->pml;
                int none, np;
                rat p[3];
                SP_TRY(request_params(row->gid, vector, &none, &np, p));
                if (none) {
                    ok = 0;
                    break;
                }
                if (!have || params_cmp(p, np, params[k], nparams[k]) < 0) {
                    memcpy(params[k], p, sizeof p);
                    nparams[k] = np;
                    have = 1;
                }
            }
        }
        if (!ok) continue;
        for (int k = 0; k < m; k++) {
            order[k] = k;
            saved[k] = rows->items[indices[k]];
        }
        for (int x = 1; x < m; x++) {
            int item = order[x], y = x - 1;
            while (y >= 0 && params_cmp(params[order[y]], nparams[order[y]], params[item], nparams[item]) > 0) {
                order[y + 1] = order[y];
                y--;
            }
            order[y + 1] = item;
        }
        for (int k = 0; k < m; k++) rows->items[indices[k]] = saved[order[k]];
    }
    return 0;
}

/* ------------------------------------------------------------ results */

static sp_specs result;
static int64_t *result_ints;
static double *result_doubles;
static size_t ints_used, ints_cap, doubles_used, doubles_cap;

static int put_int(int64_t value)
{
    if (ints_used == ints_cap) {
        size_t cap = ints_cap ? ints_cap * 2 : 4096;
        int64_t *grown = realloc(result_ints, sizeof(int64_t) * cap);
        if (!grown) return sp_fail(PY_MEMORY, "spec stage ran out of memory");
        result_ints = grown;
        ints_cap = cap;
    }
    result_ints[ints_used++] = value;
    return 0;
}

static int put_double(double value)
{
    if (doubles_used == doubles_cap) {
        size_t cap = doubles_cap ? doubles_cap * 2 : 4096;
        double *grown = realloc(result_doubles, sizeof(double) * cap);
        if (!grown) return sp_fail(PY_MEMORY, "spec stage ran out of memory");
        result_doubles = grown;
        doubles_cap = cap;
    }
    result_doubles[doubles_used++] = value;
    return 0;
}

static int put_str(sp_str s)
{
    SP_TRY(put_int(s.length));
    for (int64_t i = 0; i < s.length; i += 8) {
        int64_t word = 0;
        memcpy(&word, s.text + i, (size_t)(s.length - i < 8 ? s.length - i : 8));
        SP_TRY(put_int(word));
    }
    return 0;
}

static int put_rats(const rat *values, int n)
{
    SP_TRY(put_int(n));
    for (int i = 0; i < n; i++) {
        SP_TRY(put_int((int64_t)values[i].n));
        SP_TRY(put_int((int64_t)values[i].d));
    }
    return 0;
}

static int put_matrix(const double *m, int dim, int count)
{
    SP_TRY(put_int(dim));
    SP_TRY(put_int(count));
    SP_TRY(put_int((int64_t)doubles_used));
    for (int i = 0; i < dim * count; i++) SP_TRY(put_double(m[i]));
    return 0;
}

/* Serialized record per spec: fixed header then optional parts, in the
 * order the Python boundary reads them. */
static int emit(const sp_specs *specs)
{
    ints_used = doubles_used = 0;
    SP_TRY(put_int(specs->count));
    for (int i = 0; i < specs->count; i++) {
        const sp_spec *s = specs->items[i];
        int64_t header[] = {s->kind, s->mode_kind, s->old_id, s->gid, s->little_type, s->magnetic, s->primary, s->has_kp, s->has_carrier,
                            s->has_pml, s->has_req, s->req_none, s->has_occ, s->has_dir, s->has_free, s->free_count, s->obs,
                            s->has_anchor, s->anchor_gid, s->has_cand, s->has_heading, s->heading_gid, s->stokes, s->case_none,
                            s->dcase_none, s->slot_order, s->replaced, s->row_id, s->domain, s->domain_old, s->frequency,
                            s->use_selected_rows, s->static_updated};
        SP_TRY(put_int((int64_t)(sizeof header / sizeof header[0])));
        for (size_t k = 0; k < sizeof header / sizeof header[0]; k++) SP_TRY(put_int(header[k]));
        SP_TRY(put_str(s->label));
        SP_TRY(put_str(s->k_label));
        SP_TRY(put_str(s->opd));
        SP_TRY(put_str(s->heading_label));
        for (int k = 0; k < 4; k++) SP_TRY(put_int(s->kp[k]));
        for (int k = 0; k < 4; k++) SP_TRY(put_int(s->carrier[k]));
        for (int k = 0; k < 4; k++) SP_TRY(put_int(s->anchor_kp[k]));
        for (int k = 0; k < 4; k++) SP_TRY(put_int(s->cand_kp[k]));
        for (int k = 0; k < 4; k++) SP_TRY(put_int(s->heading_kp[k]));
        SP_TRY(put_rats(s->pml, s->has_pml ? 3 : 0));
        SP_TRY(put_rats(s->req, s->req_n));
        SP_TRY(put_rats(s->case_k, s->case_none ? 0 : s->case_n));
        SP_TRY(put_rats(s->dcase_k, s->dcase_none ? 0 : s->dcase_n));
        SP_TRY(put_int(s->nocc));
        for (int k = 0; k < s->nocc; k++) SP_TRY(put_int(s->occ[k]));
        SP_TRY(put_matrix(s->dir, s->has_dir ? s->dim : 0, s->has_dir ? s->count : 0));
        SP_TRY(put_matrix(s->recomputed, s->stokes ? s->rdim : 0, s->stokes ? s->rcount : 0));
    }
    return 0;
}

int64_t spec_result_ints(const int64_t **ints, const double **doubles)
{
    *ints = result_ints;
    *doubles = result_doubles;
    return (int64_t)ints_used;
}

/* The occurrence pool: mesh, kp[4], pml (n, d) x3, bucket[6], label. */
int64_t spec_result_occurrences(int64_t *out, int capacity)
{
    if (out)
        for (int i = 0; i < sp_occ_count && i < capacity; i++) {
            const sp_occ *o = &sp_occ_pool[i];
            int64_t *w = out + 18 * i;
            w[0] = o->mesh;
            memcpy(w + 1, o->kp, 4 * sizeof(int64_t));
            for (int axis = 0; axis < 3; axis++) {
                w[5 + 2 * axis] = (int64_t)o->pml[axis].n;
                w[6 + 2 * axis] = (int64_t)o->pml[axis].d;
            }
            memcpy(w + 11, o->bucket, 6 * sizeof(int64_t));
            w[17] = 0;
        }
    return sp_occ_count;
}

int spec_result_occurrence_label(int index, const char **text, int64_t *length)
{
    if (index < 0 || index >= sp_occ_count) return 1;
    *text = sp_occ_pool[index].k_label.text;
    *length = sp_occ_pool[index].k_label.length;
    return 0;
}

const char *spec_error(int *kind)
{
    *kind = sp_err.kind;
    return sp_err.message;
}

static int finish(int status)
{
    if (status) return sp_err.kind ? sp_err.kind : PY_RUNTIME;
    if (emit(&result)) return sp_err.kind;
    return 0;
}

static int frame_operations(int sg, int child_sg, const int64_t *basis, const int64_t *origin, sp_frame *frame)
{
    int64_t (*ops)[5] = sp_alloc(sizeof(int64_t) * 5 * GEO_MAX_OPS);
    if (!ops) return 1;
    int n = parent_operation_records(sg, child_sg, basis, origin, ops);
    if (n <= -10) {
        char message[160];
        snprintf(message, sizeof message, "child SG%d operation %d does not map to a point operation in SG%d", child_sg, -10 - n, sg);
        return sp_fail(PY_VALUE, message);
    }
    if (n < 0) return sp_fail(PY_VALUE, "operation records failed");
    frame->sg = sg;
    memcpy(frame->basis, basis, sizeof frame->basis);
    frame->ops = (const int64_t (*)[5])ops;
    frame->nops = n;
    frame->magnetic = 0;
    return 0;
}

/* ------------------------------------------------------------ entry A: ordinary dynamic */

/* _subduced_mode_specs for old_id <= 0 after the fallback checks: specs
 * in final order; a SK_FALLBACK record marks where the fallback goes. */
int spec_ordinary_dynamic(int sg, int child_sg, const int64_t *basis, const int64_t *origin, int selected_gid, int has_kp,
                          const int64_t *selected_kp, int has_fallback_direction)
{
    sp_reset();
    result.items = 0;
    result.count = result.cap = 0;
    if (load()) return finish(1);
    sp_frame frame;
    if (frame_operations(sg, child_sg, basis, origin, &frame)) return finish(1);
    sp_rows raw = {0}, rows = {0};
    if (sp_frame_rows(&frame, 0, 0, &raw)) return finish(1);
    if (selected_gid && has_kp) {
        int64_t kps[1][4];
        memcpy(kps[0], selected_kp, sizeof kps[0]);
        if (promote_dynamic(&frame, &raw, &selected_gid, (const int64_t (*)[4])kps, 1, 0, &rows)) return finish(1);
        if (replay_type3(&frame, &rows, selected_gid, selected_kp)) return finish(1);
    } else if (rows_copy(&raw, &rows)) {
        return finish(1);
    }
    builder b = {build_dyn, selected_gid, has_kp, {0}, 0, {0}, 0, 0, 0, 0, 0, 1};
    if (has_kp) memcpy(b.primary_kp, selected_kp, sizeof b.primary_kp);
    sp_specs specs = {0}, canonical = {0};
    for (int r = 0; r < rows.count; r++) {
        sp_spec *spec;
        if (build_dyn(&b, &rows.items[r], &spec) || specs_push(&specs, spec)) return finish(1);
    }
    alias_mode mode = {&frame, 1};
    if (stored_alias_specs(&mode, &raw, &rows, &b, &specs)) return finish(1);
    if (canonicalize(&specs, &canonical)) return finish(1);
    int any_primary = 0;
    for (int i = 0; i < canonical.count; i++) {
        sp_spec *spec = canonical.items[i];
        if (spec->primary && has_kp && spec->has_kp && sp_kp_exact_eq(spec->kp, selected_kp)) {
            spec->has_free = 1;
            spec->free_count = -1;
            spec->use_selected_rows = has_fallback_direction;
        }
        any_primary |= spec->primary;
    }
    if (!canonical.count || !any_primary) {
        sp_spec *fallback = spec_new();
        if (!fallback) return finish(1);
        fallback->kind = SK_FALLBACK;
        fallback->primary = 1;
        sp_specs with = {0};
        if (specs_push(&with, fallback)) return finish(1);
        if (canonical.count)
            for (int i = 0; i < canonical.count; i++)
                if (specs_push(&with, canonical.items[i])) return finish(1);
        canonical = with;
    }
    for (int i = 0; i < canonical.count; i++)
        if (canonical.items[i]->primary && specs_push(&result, canonical.items[i])) return finish(1);
    for (int i = 0; i < canonical.count; i++)
        if (!canonical.items[i]->primary && specs_push(&result, canonical.items[i])) return finish(1);
    return finish(0);
}

/* ------------------------------------------------------------ entry A: static */

static sp_str k_label_from_irrep_label(sp_str label)
{
    sp_str text = label;
    while (text.length && (text.text[text.length - 1] == ' ')) text.length--;
    if (text.length > 1 && text.text[0] == 'm' && text.text[1] >= 'A' && text.text[1] <= 'Z') {
        text.text++;
        text.length--;
    }
    int64_t n = 0;
    while (n < text.length && ((text.text[n] >= 'A' && text.text[n] <= 'Z') || (text.text[n] >= 'a' && text.text[n] <= 'z'))) n++;
    if (!n) return text;
    return (sp_str){text.text, n};
}

static int upper_prefix(sp_str label, sp_str *out)
{
    int64_t n = 0;
    while (n < label.length && label.text[n] >= 'A' && label.text[n] <= 'Z') n++;
    *out = (sp_str){label.text, n};
    return n > 0;
}

/* The canonical little gid of an old irrep among little_records_for_k(sg, k_label), 0 for none. */
static int canonical_little(int sg, sp_str k_label, int old_id, int *gid_out, sp_str *label_out)
{
    int64_t total = section_count(&T.irr_sg);
    int *seen_old = sp_alloc(sizeof(int) * 2048);
    if (!seen_old) return 1;
    int nseen = 0;
    *gid_out = 0;
    for (int64_t index = 0; index < total; index++) {
        if (section_int(&T.irr_sg, index) != sg) continue;
        int gid = (int)index + 1;
        int old = (int)section_int(&T.old, index);
        int key = old > 0 ? old : -gid;
        int dup = 0;
        for (int i = 0; i < nseen && !dup; i++) dup = seen_old[i] == key;
        if (dup) continue;
        sp_str label = sp_irrep_full_label(gid), prefix;
        if (!upper_prefix(label, &prefix)) {
            char message[160];
            snprintf(message, sizeof message, "cannot infer k label from irrep label '%.*s'", (int)label.length, label.text);
            return sp_fail(PY_VALUE, message);
        }
        if (!sp_str_eq(prefix, k_label)) continue;
        if (nseen < 2048) seen_old[nseen++] = key;
        if (old == old_id && !*gid_out) {
            *gid_out = gid;
            *label_out = label;
        }
    }
    return 0;
}

int spec_static(int sg, int selected_old_id, int selected_row_id, int has_frame, int child_sg, const int64_t *basis,
                const int64_t *origin)
{
    sp_reset();
    result.items = 0;
    result.count = result.cap = 0;
    if (load()) return finish(1);
    sp_specs specs = {0};
    int64_t pointers = section_count(&T.sub_pointer);
    if (selected_row_id >= 1 && selected_row_id <= pointers) {
        int64_t start = section_int(&T.sub_pointer, selected_row_id - 1);
        int64_t end = selected_row_id < pointers ? section_int(&T.sub_pointer, selected_row_id) : section_count(&T.sub_irrep) + 1;
        int64_t (*seen)[2] = sp_alloc(sizeof(int64_t) * 2 * (size_t)(end - start + 1));
        if (!seen) return finish(1);
        int nseen = 0;
        for (int64_t index = start - 1; index < end - 1; index++) {
            int old_id = (int)section_int(&T.sub_irrep, index);
            int entry_sg = (int)section_int(&T.irrep_sg, old_id - 1);
            if (entry_sg != sg) continue;
            int row_id = (int)section_int(&T.sub_subgroup, index);
            int dup = 0;
            for (int i = 0; i < nseen && !dup; i++) dup = seen[i][0] == old_id && seen[i][1] == row_id;
            if (dup) continue;
            seen[nseen][0] = old_id;
            seen[nseen++][1] = row_id;
            sp_str entry_label = stripped(&T.irrep_label, old_id - 1);
            int gid;
            sp_str little_label = entry_label;
            if (canonical_little(sg, k_label_from_irrep_label(entry_label), old_id, &gid, &little_label)) return finish(1);
            sp_spec *spec = spec_new();
            if (!spec) return finish(1);
            spec->kind = SK_STATIC;
            spec->mode_kind = 1;
            spec->old_id = old_id;
            spec->gid = gid;
            spec->label = gid ? little_label : entry_label;
            spec->k_label = k_label_from_irrep_label(entry_label);
            spec->row_id = row_id;
            spec->opd = stripped(&T.op_label, row_id - 1);
            spec->domain = (int)section_int(&T.sub_domain, index);
            spec->domain_old = (int)section_int(&T.sub_domain_old, index);
            spec->frequency = (int)section_int(&T.sub_frequency, index);
            spec->primary = old_id == selected_old_id && row_id == selected_row_id;
            if (specs_push(&specs, spec)) return finish(1);
        }
    }
    for (int i = 0; i < specs.count; i++) {
        sp_spec *spec = specs.items[i];
        spec->has_free = 1;
        spec->free_count = (int)section_int(&T.freeparam, spec->row_id - 1);
    }
    if (has_frame) {
        sp_frame frame;
        if (frame_operations(sg, child_sg, basis, origin, &frame)) return finish(1);
        sp_rows rows = {0};
        if (sp_frame_rows(&frame, 0, 0, &rows)) return finish(1);
        for (int i = 0; i < specs.count; i++) {
            sp_spec *spec = specs.items[i];
            if (spec->primary) continue;
            const sp_row *selected = 0;
            for (int r = 0; r < rows.count; r++) {
                int old = sp_little_old_id(rows.items[r].gid);
                if (old > 0 && old == spec->old_id) selected = &rows.items[r];
            }
            if (!selected) continue;
            spec->static_updated = 1;
            spec->gid = selected->gid;
            spec->has_kp = 1;
            memcpy(spec->kp, selected->kp, sizeof spec->kp);
            spec->has_pml = 1;
            memcpy(spec->pml, selected->pml, sizeof spec->pml);
            spec->has_dir = 1;
            spec->dir = selected->dir;
            spec->dim = selected->dim;
            spec->count = selected->count;
        }
        for (int r = 0; r < rows.count; r++) {
            if (sp_little_old_id(rows.items[r].gid) != 0) continue;
            sp_spec *spec;
            if (dyn_spec(&rows.items[r], 0, 0, &spec) || specs_push(&specs, spec)) return finish(1);
        }
    }
    for (int i = 0; i < specs.count; i++)
        if (specs.items[i]->primary && specs_push(&result, specs.items[i])) return finish(1);
    for (int i = 0; i < specs.count; i++)
        if (!specs.items[i]->primary && specs_push(&result, specs.items[i])) return finish(1);
    return finish(0);
}

/* ------------------------------------------------------------ entry C: ordinary from a magnetic embedding */

static int ordinary_ops(const int64_t *ops, int nops, int64_t (**out)[5])
{
    int64_t (*o)[5] = sp_alloc(sizeof(int64_t) * 5 * (size_t)(nops + 1));
    if (!o) return 1;
    for (int i = 0; i < nops; i++) {
        memcpy(o[i], ops + 5 * i, 4 * sizeof(int64_t));
        o[i][4] = section_int(&T.mag2nonmag, ops[5 * i + 4] - 1);
    }
    *out = o;
    return 0;
}

int spec_magnetic_embedding_ordinary(int sg, const int64_t *basis, int nops, const int64_t *magnetic_ops)
{
    sp_reset();
    result.items = 0;
    result.count = result.cap = 0;
    if (load()) return finish(1);
    int64_t (*ops)[5];
    if (ordinary_ops(magnetic_ops, nops, &ops)) return finish(1);
    sp_frame frame = {sg, {0}, (const int64_t (*)[5])ops, nops, 0};
    memcpy(frame.basis, basis, sizeof frame.basis);
    sp_rows rows = {0};
    if (sp_frame_rows(&frame, 0, 0, &rows)) return finish(1);
    builder b = {build_memb, 0, 0, {0}, 0, {0}, 0, 0, 0, 0, 0, 1};
    sp_specs specs = {0}, aliases = {0}, merged = {0};
    for (int r = 0; r < rows.count; r++) {
        sp_spec *spec;
        if (build_memb(&b, &rows.items[r], &spec) || specs_push(&specs, spec)) return finish(1);
    }
    alias_mode mode = {&frame, 1};
    if (stored_alias_specs(&mode, &rows, &rows, &b, &aliases)) return finish(1);
    if (merge_alias_specs(&specs, &aliases, &merged)) return finish(1);
    if (canonicalize(&merged, &result)) return finish(1);
    return finish(0);
}

/* ------------------------------------------------------------ entry B: magnetic */

int spec_magnetic(int sg, const int64_t *basis, int nops, const int64_t *ops, int selected_gid, const int64_t *selected_kp,
                  const char *k_label, int64_t k_label_length, int sel_dim, int sel_count, const double *selected_direction)
{
    sp_reset();
    result.items = 0;
    result.count = result.cap = 0;
    if (load()) return finish(1);
    sp_frame frame = {sg, {0}, (const int64_t (*)[5])ops, nops, 1};
    memcpy(frame.basis, basis, sizeof frame.basis);
    sp_preferred preferred = {{k_label, k_label_length}, {0}};
    memcpy(preferred.kp, selected_kp, sizeof preferred.kp);
    sp_rows raw = {0}, rows = {0};
    if (sp_frame_rows(&frame, &preferred, 1, &raw)) return finish(1);
    mag_case c = {selected_gid, {0}, 0, 0, {{0}}};
    memcpy(c.kp, selected_kp, sizeof c.kp);
    if (promote_magnetic(&frame, &raw, &c, 1, &rows)) return finish(1);
    double *sel = sp_alloc(sizeof(double) * (size_t)(sel_dim * sel_count + 1));
    if (!sel) return finish(1);
    memcpy(sel, selected_direction, sizeof(double) * (size_t)(sel_dim * sel_count));
    builder b = {build_mag, 0, 0, {0}, selected_gid, {0}, sel, sel_dim, sel_count, &rows, 0, 2};
    memcpy(b.selected_kp, selected_kp, sizeof b.selected_kp);
    sp_specs specs = {0}, aliases = {0}, merged = {0}, canonical = {0};
    for (int r = 0; r < rows.count; r++) {
        sp_spec *spec;
        if (mag_spec(&b, &rows.items[r], r, &spec) || specs_push(&specs, spec)) return finish(1);
    }
    alias_mode mode = {&frame, 0};
    if (stored_alias_specs(&mode, &raw, &rows, &b, &aliases)) return finish(1);
    if (merge_alias_specs(&specs, &aliases, &merged)) return finish(1);
    if (canonicalize(&merged, &canonical)) return finish(1);
    for (int i = 0; i < canonical.count; i++)
        if (canonical.items[i]->primary && specs_push(&result, canonical.items[i])) return finish(1);
    for (int i = 0; i < canonical.count; i++)
        if (!canonical.items[i]->primary && specs_push(&result, canonical.items[i])) return finish(1);
    return finish(0);
}

/* ------------------------------------------------------------ entry D: coupled */

typedef struct {
    int order, magnetic, gid;
    sp_str k_label;
    int64_t kp[4];
    int ncase, ndcase;
    rat case_k[3], dcase_k[3];
    int dim, count;
    const double *dir; /* selected direction: transpose of the slot's source rows */
} cpl_slot;

struct cpl {
    int sg;
    cpl_slot *slots;
    int nslots;
    unsigned char *listed_mag;
};

static int primary_slot_for_row(struct cpl *c, const sp_row *row, int magnetic, const cpl_slot **out)
{
    *out = 0;
    for (int i = 0; i < c->nslots; i++) {
        const cpl_slot *slot = &c->slots[i];
        if (slot->magnetic != magnetic || slot->gid != row->gid) continue;
        if (same_source_kparam(slot->kp, row->kp)) {
            *out = slot;
            return 0;
        }
        int64_t standard[4];
        int ok, raised;
        SP_TRY(sp_standard_kparam(c->sg, row->gid, slot->kp, standard, &ok, &raised));
        if (!ok) {
            if (raised && !(BIT(raised) & (BIT(PY_INDEX) | BIT(PY_KEY) | BIT(PY_TYPE) | BIT(PY_VALUE) | BIT(PY_ZERO))))
                return sp_fail(raised, sp_standard_message());
            continue;
        }
        if (same_source_kparam(standard, row->kp)) {
            *out = slot;
            return 0;
        }
    }
    if (magnetic)
        for (int i = 0; i < c->nslots; i++) {
            const cpl_slot *slot = &c->slots[i];
            if (slot->magnetic != magnetic || slot->gid != row->gid || c->listed_mag[i]) continue;
            int projects;
            SP_TRY(projects_case(slot->gid, row, slot->case_k, slot->ncase, &projects));
            if (projects && !sp_kp_periodic_eq(row->kp, slot->kp)) projects = 0;
            if (projects) {
                *out = slot;
                return 0;
            }
        }
    return 0;
}

static int same_subspace(const sp_row *row, const cpl_slot *slot, int *out)
{
    int occ_rank, sel_rank, transform_rank;
    SP_TRY(sp_matrix_rank(row->dir, row->dim, row->count, 1e-10, &occ_rank));
    SP_TRY(sp_matrix_rank(slot->dir, slot->dim, slot->count, 1e-10, &sel_rank));
    double *transform = sp_alloc(sizeof(double) * (size_t)(row->count * slot->count + 1));
    if (!transform) return 1;
    SP_TRY(sp_lstsq(row->dir, row->dim, row->count, slot->dir, slot->count, transform, &transform_rank));
    double limit = row->has_carrier ? 1e-10 : nextafter(5e-7, INFINITY);
    *out = 0;
    if (occ_rank != sel_rank || transform_rank != sel_rank) return 0;
    int rank;
    SP_TRY(sp_matrix_rank(transform, row->count, slot->count, 1e-10, &rank));
    if (rank != sel_rank) return 0;
    double *product = sp_alloc(sizeof(double) * (size_t)(row->dim * slot->count + 1));
    if (!product) return 1;
    sp_matmul(row->dir, row->dim, row->count, transform, slot->count, product);
    if (row->dim * slot->count == 0) return sp_fail(PY_VALUE, "zero-size array to reduction operation maximum which has no identity");
    double worst = -INFINITY;
    int nan = 0;
    for (int i = 0; i < row->dim * slot->count; i++) {
        double value = fabs(product[i] - slot->dir[i]);
        if (isnan(value)) nan = 1;
        if (value > worst) worst = value;
    }
    *out = !nan && worst <= limit;
    return 0;
}

static int cpl_spec(builder *b, const sp_row *row, sp_spec **out)
{
    struct cpl *c = b->cpl;
    int magnetic = b->mode_kind == 2;
    const cpl_slot *slot;
    SP_TRY(primary_slot_for_row(c, row, magnetic, &slot));
    sp_spec *spec = spec_new();
    if (!spec) return 1;
    spec->kind = SK_CPL;
    spec->mode_kind = b->mode_kind;
    spec_from_row(spec, row);
    spec->old_id = sp_little_old_id(row->gid);
    spec->little_type = sp_little_type(row->gid);
    spec->magnetic = magnetic;
    spec->has_carrier = row->has_carrier;
    memcpy(spec->carrier, row->carrier, sizeof spec->carrier);
    if (slot) {
        if (slot->dim != row->dim || (row->dim && slot->count != row->count))
            return sp_fail(PY_VALUE, "coupled Source direction has the wrong shape");
        int same;
        SP_TRY(same_subspace(row, slot, &same));
        if (!same && (!row->has_carrier || magnetic)) {
            spec->dir = slot->dir;
            spec->dim = slot->dim;
            spec->count = slot->count;
            spec->replaced = 1;
        }
        spec->case_none = 0;
        spec->case_n = slot->ncase;
        memcpy(spec->case_k, slot->case_k, sizeof spec->case_k);
        spec->dcase_none = 0;
        spec->dcase_n = slot->ndcase;
        memcpy(spec->dcase_k, slot->dcase_k, sizeof spec->dcase_k);
        spec->slot_order = slot->order;
    }
    spec->has_free = 1;
    spec->free_count = spec->dim ? spec->count : 0;
    spec->primary = slot != 0;
    *out = spec;
    return 0;
}

typedef struct {
    int64_t rank;
    sp_spec *spec;
} ranked;

static ranked *ranked_items;
static int nranked, ranked_cap, source_order;

static int append_specs(struct cpl *c, const sp_specs *specs, int magnetic)
{
    for (int i = 0; i < specs->count; i++) {
        sp_spec *spec = specs->items[i];
        const cpl_slot *primary = 0;
        if (spec->primary)
            for (int s = 0; s < c->nslots && !primary; s++) {
                const cpl_slot *slot = &c->slots[s];
                int spec_n = spec->case_none ? 0 : spec->case_n;
                if (slot->magnetic == magnetic && slot->gid == spec->gid && !params_cmp(slot->case_k, slot->ncase, spec->case_k, spec_n))
                    primary = slot;
            }
        spec->primary = spec->primary && primary;
        spec->slot_order = primary ? primary->order : -1;
        int64_t rank = primary ? primary->order : c->nslots + source_order;
        if (nranked == ranked_cap) {
            int cap = ranked_cap ? ranked_cap * 2 : 32;
            ranked *grown = sp_alloc(sizeof(ranked) * (size_t)cap);
            if (!grown) return 1;
            if (nranked) memcpy(grown, ranked_items, sizeof(ranked) * (size_t)nranked);
            ranked_items = grown;
            ranked_cap = cap;
        }
        ranked_items[nranked].rank = (spec->primary ? 0 : 1) * 1000000 + rank;
        ranked_items[nranked++].spec = spec;
        source_order++;
    }
    return 0;
}

/* slots: per slot (order, magnetic, gid, kp[4], k_label index into labels,
 * ncase, case (n, d) x3, ndcase, dcase (n, d) x3, dim, count, direction
 * offset into directions); preferred: ordinary entries (label index, kp[4]). */
int spec_coupled(int sg, const int64_t *basis, int nops, const int64_t *ops, int magnetic_embedding, int include_displacive,
                 int include_magnetic, int nslots, const int64_t *slot_words, const double *directions, const char *const *labels,
                 const int64_t *label_lengths, int npreferred, const int64_t *preferred_words, int nmag_cases,
                 const int64_t *mag_case_words)
{
    sp_reset();
    result.items = 0;
    result.count = result.cap = 0;
    ranked_items = 0;
    nranked = ranked_cap = source_order = 0;
    if (load()) return finish(1);
    struct cpl c = {sg, sp_alloc(sizeof(cpl_slot) * (size_t)(nslots + 1)), nslots, sp_alloc((size_t)nslots + 1)};
    if (!c.slots || !c.listed_mag) return finish(1);
    memset(c.listed_mag, 0, (size_t)nslots + 1);
    for (int s = 0; s < nslots; s++) {
        const int64_t *w = slot_words + 25 * s;
        cpl_slot *slot = &c.slots[s];
        slot->order = (int)w[0];
        slot->magnetic = (int)w[1];
        slot->gid = (int)w[2];
        memcpy(slot->kp, w + 3, sizeof slot->kp);
        slot->k_label = (sp_str){labels[w[7]], label_lengths[w[7]]};
        slot->ncase = (int)w[8];
        for (int i = 0; i < 3; i++) slot->case_k[i] = i < slot->ncase ? rat_make(w[9 + 2 * i], w[10 + 2 * i]) : rat_int(0);
        slot->ndcase = (int)w[15];
        for (int i = 0; i < 3; i++) slot->dcase_k[i] = i < slot->ndcase ? rat_make(w[16 + 2 * i], w[17 + 2 * i]) : rat_int(0);
        slot->dim = (int)w[22];
        slot->count = (int)w[23];
        slot->dir = directions + w[24];
    }
    builder b = {cpl_spec, 0, 0, {0}, 0, {0}, 0, 0, 0, 0, &c, 1};
    if (include_displacive) {
        int64_t (*ord_ops)[5] = (int64_t (*)[5])ops;
        if (magnetic_embedding && ordinary_ops(ops, nops, &ord_ops)) return finish(1);
        sp_frame frame = {sg, {0}, (const int64_t (*)[5])ord_ops, nops, 0};
        memcpy(frame.basis, basis, sizeof frame.basis);
        sp_preferred *preferred = sp_alloc(sizeof(sp_preferred) * (size_t)(npreferred + 1));
        if (!preferred) return finish(1);
        for (int i = 0; i < npreferred; i++) {
            preferred[i].label = (sp_str){labels[preferred_words[5 * i]], label_lengths[preferred_words[5 * i]]};
            memcpy(preferred[i].kp, preferred_words + 5 * i + 1, 4 * sizeof(int64_t));
        }
        sp_rows raw = {0}, rows = {0};
        if (sp_frame_rows(&frame, preferred, npreferred, &raw)) return finish(1);
        int *gids = sp_alloc(sizeof(int) * (size_t)(nslots + 1)), *promoted = sp_alloc(sizeof(int) * (size_t)(nslots + 1));
        int64_t (*kps)[4] = sp_alloc(sizeof(int64_t) * 4 * (size_t)(nslots + 1));
        if (!gids || !promoted || !kps) return finish(1);
        int ncases = 0, npromoted = 0;
        for (int s = 0; s < nslots; s++) {
            if (c.slots[s].magnetic || c.slots[s].gid <= 0) continue;
            gids[ncases] = c.slots[s].gid;
            memcpy(kps[ncases], c.slots[s].kp, sizeof kps[0]);
            ncases++;
        }
        for (int i = 0; i < ncases; i++)
            if (!listed(&raw, gids[i], kps[i])) {
                int dup = 0;
                for (int k = 0; k < npromoted && !dup; k++) dup = promoted[k] == gids[i];
                if (!dup) promoted[npromoted++] = gids[i];
            }
        if (promote_dynamic(&frame, &raw, gids, (const int64_t (*)[4])kps, ncases, 1, &rows)) return finish(1);
        if (order_promoted_families(&rows, promoted, npromoted)) return finish(1);
        b.mode_kind = 1;
        sp_specs specs = {0}, canonical = {0};
        for (int r = 0; r < rows.count; r++) {
            sp_spec *spec;
            if (cpl_spec(&b, &rows.items[r], &spec) || specs_push(&specs, spec)) return finish(1);
        }
        for (int i = 0; i < specs.count; i++) {
            sp_spec *spec = specs.items[i];
            if (!spec->primary || !spec->nocc) continue;
            int first_index = -1;
            for (int k = 0; k < i; k++) {
                sp_spec *other = specs.items[k];
                if (other->primary && occs_eq(other->occ, other->nocc, spec->occ, spec->nocc)) {
                    first_index = k;
                    break;
                }
            }
            if (first_index >= 0) continue;
            sp_spec *first = spec;
            for (int k = i + 1; k < specs.count; k++) {
                sp_spec *other = specs.items[k];
                if (other->primary && occs_eq(other->occ, other->nocc, spec->occ, spec->nocc) && other->slot_order < first->slot_order) first = other;
            }
            sp_spec displayed = *first;
            for (int k = i; k < specs.count; k++) {
                sp_spec *other = specs.items[k];
                if (other->primary && occs_eq(other->occ, other->nocc, spec->occ, spec->nocc)) {
                    other->dcase_none = displayed.dcase_none;
                    other->dcase_n = displayed.dcase_n;
                    memcpy(other->dcase_k, displayed.dcase_k, sizeof other->dcase_k);
                }
            }
        }
        alias_mode mode = {&frame, 1};
        if (stored_alias_specs(&mode, &raw, &rows, &b, &specs)) return finish(1);
        if (canonicalize(&specs, &canonical)) return finish(1);
        if (append_specs(&c, &canonical, 0)) return finish(1);
    }
    if (include_magnetic && magnetic_embedding) {
        sp_frame frame = {sg, {0}, (const int64_t (*)[5])ops, nops, 1};
        memcpy(frame.basis, basis, sizeof frame.basis);
        sp_preferred *preferred = sp_alloc(sizeof(sp_preferred) * (size_t)(nslots + 1));
        if (!preferred) return finish(1);
        int npref = 0;
        for (int s = 0; s < nslots; s++) {
            if (!c.slots[s].magnetic) continue;
            int seen = 0;
            for (int k = 0; k < s && !seen; k++) seen = c.slots[k].magnetic && sp_str_eq(c.slots[k].k_label, c.slots[s].k_label);
            if (seen) continue;
            for (int k = s; k < nslots; k++)
                if (c.slots[k].magnetic && sp_str_eq(c.slots[k].k_label, c.slots[s].k_label)) {
                    preferred[npref].label = c.slots[k].k_label;
                    memcpy(preferred[npref++].kp, c.slots[k].kp, 4 * sizeof(int64_t));
                }
        }
        sp_rows raw = {0}, rows = {0};
        if (sp_frame_rows(&frame, preferred, npref, &raw)) return finish(1);
        mag_case *cases = sp_alloc(sizeof(mag_case) * (size_t)(nmag_cases + 1));
        if (!cases) return finish(1);
        for (int i = 0; i < nmag_cases; i++) {
            const int64_t *w = mag_case_words + 12 * i;
            cases[i].gid = (int)w[0];
            memcpy(cases[i].kp, w + 1, 4 * sizeof(int64_t));
            cases[i].has_case = 1;
            cases[i].ncase = (int)w[5];
            for (int k = 0; k < 3; k++) cases[i].case_k[k] = k < cases[i].ncase ? rat_make(w[6 + 2 * k], w[7 + 2 * k]) : rat_int(0);
        }
        if (promote_magnetic(&frame, &raw, cases, nmag_cases, &rows)) return finish(1);
        for (int s = 0; s < nslots; s++) {
            if (!c.slots[s].magnetic) continue;
            for (int r = 0; r < raw.count && !c.listed_mag[s]; r++)
                c.listed_mag[s] = raw.items[r].gid == c.slots[s].gid && same_source_kparam(raw.items[r].kp, c.slots[s].kp);
        }
        b.mode_kind = 2;
        sp_specs specs = {0}, canonical = {0};
        for (int r = 0; r < rows.count; r++) {
            sp_spec *spec;
            if (cpl_spec(&b, &rows.items[r], &spec) || specs_push(&specs, spec)) return finish(1);
        }
        alias_mode mode = {&frame, 0};
        if (stored_alias_specs(&mode, &raw, &rows, &b, &specs)) return finish(1);
        if (canonicalize(&specs, &canonical)) return finish(1);
        if (append_specs(&c, &canonical, 1)) return finish(1);
    }
    for (int x = 1; x < nranked; x++) {
        ranked item = ranked_items[x];
        int y = x - 1;
        while (y >= 0 && ranked_items[y].rank > item.rank) {
            ranked_items[y + 1] = ranked_items[y];
            y--;
        }
        ranked_items[y + 1] = item;
    }
    for (int i = 0; i < nranked; i++)
        if (specs_push(&result, ranked_items[i].spec)) return finish(1);
    return finish(0);
}
