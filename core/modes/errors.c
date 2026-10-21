/* The error kind and text a stage error code stands for, as the CLI and the GUI report it. */

#include <stdio.h>
#include <string.h>

#include "json.h"
#include "md_reason_tables.h"
#include "pyconv.h"

typedef struct {
    const char *type;
    jbuf *text;
} reason;

static const char *const generic_types[] = {0,           "ValueError",      "KeyError",    "IndexError",   "TypeError",
                                            "ZeroDivisionError", "ArithmeticError", "MemoryError", "RuntimeError", "LinAlgError"};

static void set(reason *r, const char *type, const char *text)
{
    r->type = type;
    jb_lit(r->text, text);
}

static void setf(reason *r, const char *type, const char *format, int value)
{
    char text[96];
    snprintf(text, sizeof text, format, value);
    set(r, type, text);
}

static int lookup(const md_reason_entry *table, int status, reason *r)
{
    for (; table->type; table++)
        if (table->status == status) {
            set(r, table->type, table->text);
            return 1;
        }
    return 0;
}

static void table_or(reason *r, const md_reason_entry *table, int status, const char *fallback, int shown)
{
    if (!lookup(table, status, r)) setf(r, "RuntimeError", fallback, shown);
}

/* Python divmod(status, 1000) */
static void divmod1000(int status, int *family, int *code)
{
    *family = status / 1000;
    *code = status % 1000;
    if (*code < 0) {
        *family -= 1;
        *code += 1000;
    }
}

static void child_layout(int status, reason *r)
{
    int family, code;
    divmod1000(status, &family, &code);
    switch (family) {
    case 1:
        if (code >= 100) table_or(r, wyckoff_split_errors, code - 100, "Wyckoff split failed (%d)", code - 100);
        else table_or(r, child_site_errors, code, "child sites failed (%d)", code);
        return;
    case 2: table_or(r, child_setting_errors, code, "child setting failed (%d)", code); return;
    case 3: table_or(r, site_atom_errors, code, "site atom table failed (%d)", code); return;
    case 4:
        if (code == 1) set(r, "ValueError", "mode presentation basis is singular");
        else setf(r, "ValueError", "mode presentation failed (%d)", code);
        return;
    case 5: table_or(r, row_assignment_errors, code, "row assignment failed (%d)", code); return;
    case 6: table_or(r, magnetic_errors, code, "magnetic kernel failed (%d)", code); return;
    default: table_or(r, child_layout_errors, code, "child layout failed (%d)", status); return;
    }
}

static void orderparam(int status, reason *r) { table_or(r, orderparam_errors, status, "orderparam_to_subgroup_ failed (%d)", status); }

static void mode_detail(int kind, int status, const char *message, reason *r)
{
    switch (kind) {
    case 100: child_layout(status, r); return;
    case 101: table_or(r, child_setting_errors, status, "child setting failed (%d)", status); return;
    case 102: orderparam(status, r); return;
    case 103: table_or(r, magnetic_errors, status, "magnetic kernel failed (%d)", status); return;
    case 110: table_or(r, recordwise_errors, status, "mode block failed (%d)", status); return;
    case 111: table_or(r, projection_entry_errors, status, "mode block failed (%d)", status); return;
    case 112: table_or(r, mode_block_errors, status, "mode block failed (%d)", status); return;
    case 113:
        if (status >= 100) table_or(r, row_assignment_errors, status - 100, "row assignment failed (%d)", status - 100);
        else if (status == 6) set(r, "MemoryError", "block plan ran out of memory");
        else if (status >= 0 && status < (int)(sizeof block_plan_errors / sizeof block_plan_errors[0]) && block_plan_errors[status])
            set(r, "ValueError", block_plan_errors[status]);
        else setf(r, "ValueError", "block plan failed (%d)", status);
        return;
    default: break;
    }
    r->type = kind >= 1 && kind <= 9 ? generic_types[kind] : "RuntimeError";
    size_t n = strlen(message);
    if (kind == 2 && n >= 2 && message[0] == '\'' && message[n - 1] == '\'') jb_raw(r->text, message + 1, n - 2);
    else jb_lit(r->text, message);
}

static void isotropy_rows(int status, reason *r)
{
    if (status >= 3000) set(r, "ValueError", "magnetic Source tables are inconsistent");
    else if (status >= 2000) {
        if (status - 2000 >= 100) orderparam(status - 2100, r);
        else table_or(r, bridge_weight_errors, status - 2000, "bridge weights failed (%d)", status - 2000);
    } else if (status >= 1000) table_or(r, table_errors, status - 1000, "traversal tables failed (%d)", status - 1000);
    else if (status > 200) table_or(r, id_subgroup_errors, 200 - status, "id_subgroup_ status %d", status);
    else if (status > 100) orderparam(status, r);
    else table_or(r, kernel_errors, status, "native isotropy traversal failed (%d)", status);
}

static void coupled_domain_orbit(int status, reason *r)
{
    if (status == 20) {
        set(r, "KeyError", "domain action leaves the retained domain subspaces");
        return;
    }
    if (status == 21) {
        set(r, "ValueError", "magnetic point-operation tables are inconsistent");
        return;
    }
    if (status > 2000) {
        table_or(r, subgroup_ops_errors, status - 2000, "subgroup operation failed (%d)", status - 2000);
        return;
    }
    if (status >= 1100) {
        orderparam(status - 1100, r);
        return;
    }
    if (status >= 1000) {
        table_or(r, bridge_weight_errors, status - 1000, "bridge weights failed (%d)", status - 1000);
        return;
    }
    if (status > 100) {
        orderparam(status, r);
        return;
    }
    table_or(r, kernel_errors, status, "coupled domain orbit failed (%d)", status);
}

static void selection(int kind, int status, const char *message, reason *r)
{
    char fallback[64];
    switch (kind) {
    case 120: isotropy_rows(status, r); return;
    case 102: orderparam(status, r); return;
    case 103: table_or(r, magnetic_errors, status, "magnetic kernel failed (%d)", status); return;
    case 111:
    case 121:
    case 122:
    case 124: {
        const md_reason_entry *table = kind == 111 ? projection_entry_errors : kind == 121 ? id_subgroup_errors
                                                                           : kind == 122 ? subgroup_ops_errors
                                                                                         : bridge_weight_errors;
        if (!lookup(table, status, r)) {
            snprintf(fallback, sizeof fallback, "selection kernel failed (%d, %d)", kind, status);
            set(r, "RuntimeError", fallback);
        }
        return;
    }
    case 123: coupled_domain_orbit(status, r); return;
    case 125: r->type = "OverflowError"; jb_lit(r->text, message); return;
    case 126: r->type = "OpdNotFoundError"; jb_lit(r->text, message); return;
    default: mode_detail(kind, status, message, r); return;
    }
}

static void compose(const reason *r, jbuf *out)
{
    jb_lit(out, r->type);
    jb_lit(out, ": ");
    if (strcmp(r->type, "KeyError")) jb_raw(out, r->text->data ? r->text->data : "", r->text->used);
    else py_repr_text(out, r->text->data ? r->text->data : "", (int64_t)r->text->used);
}

/* "Type: text" of a failed mode-detail request */
void md_reason(int kind, int status, const char *message, jbuf *out)
{
    jbuf text = {0};
    jb_reset(&text);
    reason r = {0, &text};
    mode_detail(kind, status, message, &r);
    compose(&r, out);
    jb_free(&text);
}

/* The exception a failed selection or CLI request raises: its type name, and
 * str(exception) into text (a KeyError's repr of its message). */
const char *apo_error(int kind, int status, const char *message, jbuf *text)
{
    jbuf raw = {0};
    jb_reset(&raw);
    reason r = {0, &raw};
    selection(kind, status, message, &r);
    if (strcmp(r.type, "KeyError")) jb_raw(text, raw.data ? raw.data : "", raw.used);
    else py_repr_text(text, raw.data ? raw.data : "", (int64_t)raw.used);
    jb_free(&raw);
    return r.type;
}
