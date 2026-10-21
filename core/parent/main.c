/* The parent-stage entry: one JSON request in, one JSON result out. */

#include "par.h"
#include "sel.h"

static jbuf out;

const char *parent_json(const char *text, int64_t length)
{
    jdoc *doc = jdoc_parse(text, (size_t)length);
    int status;
    jb_reset(&out);
    if (!doc) {
        status = md_fail(MDE_VALUE, "parent request is not JSON");
    } else {
        const jvalue *request = jdoc_root(doc);
        const jvalue *call = jget(request, "call");
        if (jstr_is(call, "parent_state")) status = par_parent_state(jget(request, "cif_info"), &out);
        else if (jstr_is(call, "distortion"))
            status = par_distortion(jget(request, "state"), jget(request, "mode_set"), jget(request, "distortion_selection"),
                                    jtruthy(jget(request, "allow_empty_distortion")), &out);
        else status = md_fail(MDE_VALUE, "unknown parent call");
    }
    if (status) {
        jb_reset(&out);
        jb_lit(&out, "{\"error\": {\"kind\": ");
        jb_int(&out, md_err.kind);
        jb_lit(&out, ", \"status\": ");
        jb_int(&out, md_err.status);
        jb_lit(&out, ", \"message\": ");
        jb_cstring(&out, md_err.message);
        jb_lit(&out, "}}");
    }
    jdoc_free(doc);
    return out.failed ? 0 : out.data;
}
