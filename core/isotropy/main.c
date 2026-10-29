/* The selection-stage entry: one JSON request in, one JSON result out. */

#include <string.h>

#include "par.h"
#include "sel.h"

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

/* pipeline.build_state_from_parent_state's selection steps */
int sel_stage(const jvalue *state, const jvalue *opts, jbuf *b)
{
    jbuf k = {0}, r = {0}, p = {0};
    jb_reset(&k);
    jb_reset(&r);
    jb_reset(&p);
    int status = par_distortion(state, jget(opts, "mode_set"), jget(opts, "distortion_selection"), jtruthy(jget(opts, "allow_empty_distortion")), &p);
    jdoc *pdoc = status ? 0 : jdoc_parse(p.data, p.used);
    if (!status) status = sel_k_points(jdoc_root(pdoc), opts, &k);
    jdoc *kdoc = status ? 0 : jdoc_parse(k.data, k.used), *rdoc = 0;
    if (!status) status = sel_irreps(jdoc_root(kdoc), &r);
    if (!status && !jtruthy(jget(opts, "include_opd"))) {
        jb_raw(b, r.data, r.used);
    } else if (!status) {
        rdoc = jdoc_parse(r.data, r.used);
        status = sel_order_parameter_direction(jdoc_root(rdoc), opts, b);
    }
    jdoc_free(pdoc);
    jdoc_free(kdoc);
    jdoc_free(rdoc);
    jb_free(&p);
    jb_free(&k);
    jb_free(&r);
    return status;
}

const char *selection_json(const char *text, int64_t length)
{
    jdoc *doc = jdoc_parse(text, (size_t)length);
    if (!doc) {
        md_fail(MDE_VALUE, "selection request is not JSON");
        put_error();
        return out.data;
    }
    const jvalue *request = jdoc_root(doc);
    jb_reset(&out);
    int status;
    if (jstr_is(jget(request, "call"), "coupled_opd_rows")) {
        jb_lit(&out, "{\"rows\": ");
        status = sel_coupled_rows(request, &out);
        jb_lit(&out, "}");
    } else if (jstr_is(jget(request, "call"), "stage")) {
        status = sel_stage(jget(request, "state"), jget(request, "options"), &out);
    } else if (jstr_is(jget(request, "call"), "k_points")) {
        status = sel_k_points(jget(request, "state"), jget(request, "options"), &out);
    } else if (jstr_is(jget(request, "call"), "irreps")) {
        status = sel_irreps(jget(request, "state"), &out);
    } else if (jstr_is(jget(request, "call"), "order_parameter_direction")) {
        status = sel_order_parameter_direction(jget(request, "state"), jget(request, "options"), &out);
    } else {
        status = md_fail(MDE_VALUE, "unknown selection call");
    }
    if (status) put_error();
    jdoc_free(doc);
    return out.failed ? 0 : out.data;
}
