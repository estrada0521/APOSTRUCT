/* The mode-details entry: one JSON request in, one JSON result out. */

#include <stdlib.h>
#include <string.h>

#include "md.h"
#include "tables.h"

static jbuf out;

static void put_error(void)
{
    jb_reset(&out);
    jb_lit(&out, "{\"error\": {\"kind\": ");
    jb_int(&out, md_err.kind);
    jb_lit(&out, ", \"status\": ");
    jb_int(&out, md_err.status);
    jb_lit(&out, ", \"message\": ");
    jb_cstring(&out, md_err.message);
    jb_lit(&out, "}}");
}

static int load_request(md_request *r, const jvalue *request)
{
    memset(r, 0, sizeof *r);
    r->request = request;
    r->cif = jget(request, "cif_info");
    r->selected_k = jget(request, "selected_k");
    r->irrep = jget(request, "selected_irrep");
    r->opd = jget(request, "selected_opd");
    r->slots = jget(request, "selected_slots");
    r->factors = jget(request, "selected_factors");
    r->k_params = jget(request, "k_params");
    r->iso = md_isotropy(r->opd);
    const jvalue *setting = jget(request, "parent_inter_setting_id");
    if (jis(setting, J_INT)) {
        r->has_setting = 1;
        r->setting_id = setting->u.i;
    }
    const jvalue *parent = jget(r->cif, "parent");
    int64_t sg;
    if (!parent || !jget(parent, "number")) return md_fail(MDE_KEY, "'number'");
    if (py_int_of(jget(parent, "number"), &sg)) return md_fail(MDE_VALUE, "invalid parent space group");
    r->sg = (int)sg;
    if (md_tables()) return md_fail(MDE_VALUE, "mode Source tables are inconsistent");
    md_selected(r);
    const jvalue *sites = jget(r->cif, "atom_sites");
    int64_t total = jtruthy(sites) && jis(sites, J_ARRAY) ? sites->count : 0;
    r->sites = calloc((size_t)total + 1, sizeof(md_site));
    if (!r->sites) return md_fail(MDE_MEMORY, "out of memory");
    for (int64_t i = 0; i < total; i++) {
        const jvalue *site = jat(sites, i);
        if (!jis(site, J_OBJECT) || !jtruthy(jget(site, "wyckoff"))) continue;
        md_site *s = r->sites + r->nsites++;
        s->site = site;
        const jvalue *type = jget(site, "type"), *label = jget(site, "label");
        s->atom_label = jtruthy(type) ? md_py_str(type) : jtruthy(label) ? md_py_str(label) : (char *)md_strdup("X", 1);
        s->label_prefix = jtruthy(label) ? md_py_str(label) : (char *)md_strdup(s->atom_label, (int64_t)strlen(s->atom_label));
        s->site_param_count = -1;
        if (!s->atom_label || !s->label_prefix) return md_fail(MDE_MEMORY, "out of memory");
    }
    return 0;
}

int md_details(const jvalue *request, jbuf *b)
{
    md_request r;
    md_vectors.used = 0;
    int status = load_request(&r, request);
    int done = 0;
    if (!status) status = md_preflight(&r, b, &done);
    if (!status && !done) status = md_render_specs(&r);
    if (!status && !done) status = md_build(&r, b);
    md_release(&r);
    return status;
}

const char *mode_details_json(const char *text, int64_t length)
{
    jdoc *doc = jdoc_parse(text, (size_t)length);
    if (!doc) {
        md_fail(MDE_VALUE, "mode request is not JSON");
        put_error();
        return out.data;
    }
    jb_reset(&out);
    if (md_details(jdoc_root(doc), &out)) put_error();
    jdoc_free(doc);
    return out.failed ? 0 : out.data;
}

/* The mode vectors of the last mode_details_json result. */
int64_t mode_details_vectors(const double **data)
{
    *data = md_vectors.data;
    return (int64_t)md_vectors.used;
}
