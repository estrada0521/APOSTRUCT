/* reciprocal.catalog: the k-point catalog of a space group (kpoints) as JSON,
 * and the displayed k-vector texts of Source k records. */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "k_display.h"
#include "tables.h"

static struct {
    section_view k, k_dim, k_label, k_kov, k_count, label_bc, label_ml, label, lattice_label, sg_point_group, pg_order, star_count,
        star_minusk, ops_count, star_pointer, star, star_special, irr_sg, irr_k, irr_old, full_label, full_kov_label, full_dim,
        irr_type, irr_lif, real_pointer, real2_pointer, irrep_image, image_label, image_dimension, image_order, image_type,
        image_landau, image_subgroup_count, image_subgroup_pointer;
    int ready;
} T;

static int load(void)
{
    if (T.ready) return 0;
    if (geo_load() || section_find("little.little_k", &T.k) || section_find("little.little_k_dim", &T.k_dim) ||
        section_find("little.little_k_label", &T.k_label) || section_find("little.little_k_kov", &T.k_kov) ||
        section_find("little.little_k_count", &T.k_count) || section_find("space.space_label_bc", &T.label_bc) ||
        section_find("space.space_label_ml", &T.label_ml) || section_find("space.space_label", &T.label) ||
        section_find("space.lattice_label", &T.lattice_label) || section_find("space.ispace_point_group", &T.sg_point_group) ||
        section_find("space.ipoint_group_order", &T.pg_order) || section_find("little.little_k_star_count", &T.star_count) ||
        section_find("little.little_k_star_minusk", &T.star_minusk) || section_find("little.little_ops_count", &T.ops_count) ||
        section_find("little.little_k_star_ml_pointer", &T.star_pointer) || section_find("little.little_k_star", &T.star) ||
        section_find("little.little_k_star_special", &T.star_special) ||
        section_find("little.little_irr_space_group", &T.irr_sg) || section_find("little.little_irr_k", &T.irr_k) ||
        section_find("little.little_irr_old", &T.irr_old) || section_find("little.little_irr_full_label", &T.full_label) ||
        section_find("little.little_irr_full_kov_label", &T.full_kov_label) ||
        section_find("little.little_irr_full_dim", &T.full_dim) || section_find("little.little_irr_type", &T.irr_type) ||
        section_find("little.little_irr_lif", &T.irr_lif) || section_find("little.little_irr_real_pointer", &T.real_pointer) ||
        section_find("little.little_irr_real2_pointer", &T.real2_pointer) || section_find("irreps.irrep_image", &T.irrep_image) ||
        section_find("images.image_label", &T.image_label) || section_find("images.image_dimension", &T.image_dimension) ||
        section_find("images.image_order", &T.image_order) || section_find("images.image_type", &T.image_type) ||
        section_find("images.image_landau", &T.image_landau) ||
        section_find("images.image_subgroup_count", &T.image_subgroup_count) ||
        section_find("images.image_subgroup_pointer", &T.image_subgroup_pointer))
        return 1;
    T.ready = 1;
    return 0;
}

typedef struct {
    const char *text;
    int64_t length;
} str;

/* Python str(value).strip() of a string section entry. */
static str strip(const section_view *view, int64_t index)
{
    str out;
    if (section_string(view, index, &out.text, &out.length)) return (str){"", 0};
    while (out.length && (out.text[0] == ' ' || out.text[0] == '\t')) {
        out.text++;
        out.length--;
    }
    while (out.length && (out.text[out.length - 1] == ' ' || out.text[out.length - 1] == '\t')) out.length--;
    return out;
}

static int str_is(str a, const char *b) { return a.length == (int64_t)strlen(b) && !memcmp(a.text, b, (size_t)a.length); }

/* ------------------------------------------------------------ text buffer */

typedef struct {
    char *data;
    size_t used, cap;
    int failed;
} buffer;

static void put(buffer *b, const char *format, ...) __attribute__((format(printf, 2, 3)));
static void put(buffer *b, const char *format, ...)
{
    if (b->failed) return;
    for (;;) {
        va_list args;
        va_start(args, format);
        int n = vsnprintf(b->data ? b->data + b->used : 0, b->data ? b->cap - b->used : 0, format, args);
        va_end(args);
        if (n < 0) {
            b->failed = 1;
            return;
        }
        if (b->data && b->used + (size_t)n < b->cap) {
            b->used += (size_t)n;
            return;
        }
        size_t cap = b->cap ? b->cap * 2 : 4096;
        while (cap < b->used + (size_t)n + 1) cap *= 2;
        char *grown = realloc(b->data, cap);
        if (!grown) {
            b->failed = 1;
            return;
        }
        b->data = grown;
        b->cap = cap;
    }
}

static void put_json_str(buffer *b, const char *text, int64_t length)
{
    put(b, "\"");
    for (int64_t i = 0; i < length; i++) {
        unsigned char c = (unsigned char)text[i];
        if (c == '"' || c == '\\')
            put(b, "\\%c", c);
        else if (c < 0x20)
            put(b, "\\u%04x", c);
        else
            put(b, "%c", c);
    }
    put(b, "\"");
}

static void put_frac(buffer *b, rat v)
{
    if (v.d == 1)
        put(b, "%lld", (long long)v.n);
    else
        put(b, "%lld/%lld", (long long)v.n, (long long)v.d);
}

/* _term(value, symbol) */
static void put_term(buffer *b, rat v, const char *symbol)
{
    if (!symbol[0]) {
        put_frac(b, v);
        return;
    }
    if (v.n == 1 && v.d == 1) {
        put(b, "%s", symbol);
        return;
    }
    if (v.n == -1 && v.d == 1) {
        put(b, "-%s", symbol);
        return;
    }
    put_frac(b, v);
    put(b, "%s", symbol);
}

/* _fmt_vector over Source records with the α β γ parameter names. */
static void put_vector(buffer *b, const rat (*records)[3], int count)
{
    static const char *const params[4] = {"", "\xce\xb1", "\xce\xb2", "\xce\xb3"};
    put(b, "(");
    for (int axis = 0; axis < 3; axis++) {
        if (axis) put(b, ", ");
        int parts = 0;
        for (int index = 0; index < count; index++) {
            rat v = records[index][axis];
            if (rat_zero(v)) continue;
            buffer term = {0};
            put_term(&term, v, params[index]);
            if (term.failed || !term.data) {
                b->failed = 1;
                free(term.data);
                return;
            }
            if (parts && term.data[0] != '-') put(b, "+");
            put(b, "%s", term.data);
            free(term.data);
            parts++;
        }
        if (!parts) put(b, "0");
    }
    put(b, ")");
}

/* _fmt_display_k_component / _fmt_display_k_vector from a display form. */
static void put_display_vector(buffer *b, const kc_form *form)
{
    if (form->fixed) {
        put(b, "(");
        for (int axis = 0; axis < 3; axis++) {
            if (axis) put(b, ",");
            char c = form->fixed[axis];
            if (c == 'h')
                put(b, "1/2");
            else
                put(b, "%c", c);
        }
        put(b, ")");
        return;
    }
    static const char names[3] = {'a', 'b', 'g'};
    put(b, "(");
    for (int axis = 0; axis < 3; axis++) {
        if (axis) put(b, ",");
        rat constant = form->d[0][axis];
        char *parts[3];
        int nparts = 0;
        for (int index = 1; index < form->count; index++) {
            rat v = form->d[index][axis];
            if (rat_zero(v)) continue;
            char symbol[2] = {form->names[index] ? form->names[index] : names[index - 1 < 2 ? index - 1 : 2], 0};
            buffer term = {0};
            put_term(&term, v, symbol);
            parts[nparts++] = term.data;
        }
        if (nparts) {
            if (constant.d == 1 && rat_cmp(constant, rat_int(1)) > 0) constant = rat_mod1(constant);
            if (constant.d == 1 && rat_cmp(constant, rat_int(0)) > 0 && parts[0][0] == '-') {
                for (int i = 0; i < nparts; i++) {
                    if (i && parts[i][0] != '-') put(b, "+");
                    put(b, "%s", parts[i]);
                }
                put(b, "+");
                put_frac(b, constant);
                for (int i = 0; i < nparts; i++) free(parts[i]);
                continue;
            }
        }
        int written = 0;
        if (!rat_zero(constant)) {
            put_frac(b, constant);
            written = 1;
        }
        for (int i = 0; i < nparts; i++) {
            if (written && parts[i][0] != '-') put(b, "+");
            put(b, "%s", parts[i]);
            written = 1;
            free(parts[i]);
        }
        if (!written) put(b, "0");
    }
    put(b, ")");
}

/* ------------------------------------------------------------ catalog */

static void slot_records(int64_t lattice_slot, rat (*out)[3])
{
    for (int index = 0; index < 4; index++) {
        int64_t den = section_int(&T.k, lattice_slot * 16 + 4 * index + 3);
        for (int axis = 0; axis < 3; axis++)
            out[index][axis] = den ? rat_make(section_int(&T.k, lattice_slot * 16 + 4 * index + axis), den) : rat_int(0);
    }
}

static int display_form_of(int sg, int lattice, str label, const rat (*records)[3], int count, kc_form *form)
{
    str symbol = strip(&T.label_bc, sg - 1);
    return kc_display_form(sg, lattice, symbol.text, symbol.length, label.text, label.length, records, count, form);
}

static int gids_of(int sg, int kslot, int *out, int capacity)
{
    int64_t total = section_count(&T.irr_sg);
    int n = 0;
    for (int64_t index = 0; index < total; index++)
        if (section_int(&T.irr_sg, index) == sg && section_int(&T.irr_k, index) == kslot && n < capacity) out[n++] = (int)index + 1;
    return n;
}

static int has_extra_slot(int sg, int kslot, int lattice)
{
    str label = strip(&T.k_label, (int64_t)(lattice - 1) * 27 + kslot - 1);
    if (!label.length) return 0;
    if (section_int(&T.star_count, (int64_t)(sg - 1) * 27 + kslot - 1) > 0) return 1;
    int gids[1];
    return gids_of(sg, kslot, gids, 1) > 0;
}

static const char *greek(str label)
{
    if (str_is(label, "GM")) return "\xce\x93";
    if (str_is(label, "DT")) return "\xce\x94";
    if (str_is(label, "LD")) return "\xce\x9b";
    if (str_is(label, "SM")) return "\xce\xa3";
    return 0;
}

static int put_slot(buffer *b, int sg, int kslot, int lattice, int extra)
{
    int64_t lattice_slot = (int64_t)(lattice - 1) * 27 + kslot - 1, sg_slot = (int64_t)(sg - 1) * 27 + kslot - 1;
    str label = strip(&T.k_label, lattice_slot);
    int dim = (int)section_int(&T.k_dim, lattice_slot);
    int64_t kov = section_int(&T.k_kov, sg_slot);
    rat records[4][3];
    slot_records(lattice_slot, records);
    put(b, "{\"kslot\": %d, \"kid\": ", kslot);
    if (kov == 0 && !str_is(label, "GP"))
        put(b, "null");
    else
        put(b, "%lld", (long long)kov);
    put(b, ", \"label\": ");
    put_json_str(b, label.text, label.length);
    put(b, ", \"display\": ");
    const char *g = greek(label);
    if (g)
        put_json_str(b, g, (int64_t)strlen(g));
    else
        put_json_str(b, label.text, label.length);
    put(b, ", \"dimension\": %d, \"basis_records\": [", dim);
    for (int r = 0; r < 4; r++) {
        put(b, r ? ", [" : "[");
        for (int c = 0; c < 3; c++) {
            if (c) put(b, ", ");
            put(b, "\"");
            put_frac(b, records[r][c]);
            put(b, "\"");
        }
        put(b, "]");
    }
    put(b, "], \"kvector\": ");
    {
        buffer text = {0};
        put_vector(&text, (const rat (*)[3])records, dim + 1);
        if (text.failed) return 1;
        put_json_str(b, text.data, (int64_t)text.used);
        free(text.data);
    }
    put(b, ", \"display_kvector\": ");
    {
        kc_form form;
        if (display_form_of(sg, lattice, label, (const rat (*)[3])records, dim + 1, &form)) return 1;
        buffer text = {0};
        put_display_vector(&text, &form);
        if (text.failed) return 1;
        put_json_str(b, text.data, (int64_t)text.used);
        free(text.data);
    }
    int star_count = (int)section_int(&T.star_count, sg_slot);
    put(b, ", \"star_size\": %d, \"minus_k\": %lld, \"little_order\": %lld", star_count,
        (long long)section_int(&T.star_minusk, sg_slot), (long long)section_int(&T.ops_count, sg_slot));
    put(b, ", \"sg_specific\": %s, \"star\": [", extra ? "true" : "false");
    int64_t pointer = section_int(&T.star_pointer, sg_slot);
    if (star_count > 0 && pointer > 0) {
        for (int arm = 0; arm < star_count; arm++) {
            rat arm_records[4][3];
            int count;
            if (dim == 0) {
                int64_t base = 4 * (pointer - 1 + arm);
                if (base + 4 > section_count(&T.star_special)) break;
                int64_t den = section_int(&T.star_special, base + 3);
                if (!den) den = 1;
                for (int axis = 0; axis < 3; axis++) arm_records[0][axis] = rat_make(section_int(&T.star_special, base + axis), den);
                count = 1;
            } else {
                int64_t base = 16 * (pointer - 1 + arm);
                if (base + 16 > section_count(&T.star)) break;
                for (int piece = 0; piece < 4; piece++) {
                    int64_t den = section_int(&T.star, base + 4 * piece + 3);
                    for (int axis = 0; axis < 3; axis++)
                        arm_records[piece][axis] = den ? rat_make(section_int(&T.star, base + 4 * piece + axis), den) : rat_int(0);
                }
                count = dim + 1;
            }
            buffer text = {0};
            put_vector(&text, (const rat (*)[3])arm_records, count);
            if (text.failed) return 1;
            if (arm) put(b, ", ");
            put_json_str(b, text.data, (int64_t)text.used);
            free(text.data);
        }
    }
    int gids[512];
    int ngids = gids_of(sg, kslot, gids, 512);
    put(b, "], \"n_irreps\": %d, \"irreps\": [", ngids);
    for (int i = 0; i < ngids; i++) {
        int gid = gids[i];
        int64_t idx = gid - 1;
        int full_dim = (int)section_int(&T.full_dim, idx), type = (int)section_int(&T.irr_type, idx);
        int old_id = (int)section_int(&T.irr_old, idx);
        str symbol = strip(&T.full_label, idx), kov_label = strip(&T.full_kov_label, idx);
        int64_t k_length = 0;
        while (k_length < symbol.length &&
               ((symbol.text[k_length] >= 'A' && symbol.text[k_length] <= 'Z') || (symbol.text[k_length] >= 'a' && symbol.text[k_length] <= 'z')))
            k_length++;
        if (i) put(b, ", ");
        put(b, "{\"gid\": %d, \"old_id\": %d, \"symbol\": ", gid, old_id);
        put_json_str(b, symbol.text, symbol.length);
        put(b, ", \"kov_label\": ");
        put_json_str(b, kov_label.text, kov_label.length);
        put(b, ", \"k_label\": ");
        put_json_str(b, symbol.text, k_length ? k_length : symbol.length);
        put(b, ", \"full_dim\": %d, \"little_dim\": ", full_dim);
        if (star_count && full_dim % star_count == 0)
            put(b, "%d", full_dim / star_count);
        else
            put(b, "null");
        put(b, ", \"type\": %d, \"type_label\": ", type);
        if (type == 1)
            put(b, "\"real\"");
        else
            put(b, "\"type %d\"", type);
        put(b, ", \"lif\": %lld, \"real_pointer\": %lld, \"real2_pointer\": %lld, \"image\": ", (long long)section_int(&T.irr_lif, idx),
            (long long)section_int(&T.real_pointer, idx), (long long)section_int(&T.real2_pointer, idx));
        if (old_id > 0) {
            int64_t image = section_int(&T.irrep_image, old_id - 1), im = image - 1;
            str image_label = strip(&T.image_label, im);
            put(b, "{\"id\": %lld, \"label\": ", (long long)image);
            put_json_str(b, image_label.text, image_label.length);
            put(b, ", \"dimension\": %lld, \"order\": %lld, \"type\": ", (long long)section_int(&T.image_dimension, im),
                (long long)section_int(&T.image_order, im));
            if (im < section_count(&T.image_type))
                put(b, "%lld", (long long)section_int(&T.image_type, im));
            else
                put(b, "null");
            put(b, ", \"landau\": %lld, \"direction_count\": %lld, \"subgroup_pointer\": %lld}", (long long)section_int(&T.image_landau, im),
                (long long)section_int(&T.image_subgroup_count, im), (long long)section_int(&T.image_subgroup_pointer, im));
        } else {
            put(b, "null");
        }
        put(b, "}");
    }
    put(b, "]}");
    return b->failed;
}

static buffer catalog_text;

/* kpoints(sg) as JSON in the catalog's key order; NULL on failure. */
const char *k_catalog_json(int sg)
{
    free(catalog_text.data);
    memset(&catalog_text, 0, sizeof catalog_text);
    if (load()) return 0;
    buffer *b = &catalog_text;
    int lattice = geo_sg_lattice(sg);
    int count = (int)section_int(&T.k_count, lattice - 1);
    int point_group = (int)section_int(&T.sg_point_group, sg - 1);
    /* rows in display order: Source slots, then each sg-specific slot before the first higher-dimensional row */
    int order[27], dims[27], extra[27], n = 0;
    for (int kslot = 1; kslot <= 27; kslot++) {
        int is_extra = kslot > count;
        if (is_extra && !has_extra_slot(sg, kslot, lattice)) continue;
        int dim = (int)section_int(&T.k_dim, (int64_t)(lattice - 1) * 27 + kslot - 1);
        int at = n;
        if (is_extra)
            for (int i = 0; i < n; i++)
                if (dims[i] > dim) {
                    at = i;
                    break;
                }
        for (int i = n; i > at; i--) {
            order[i] = order[i - 1];
            dims[i] = dims[i - 1];
            extra[i] = extra[i - 1];
        }
        order[at] = kslot;
        dims[at] = dim;
        extra[at] = is_extra;
        n++;
    }
    str symbol = strip(&T.label_bc, sg - 1), ml = strip(&T.label_ml, sg - 1), crystal = strip(&T.label, sg - 1),
        lattice_type = strip(&T.lattice_label, lattice - 1);
    put(b, "{\"source\": \"Source/data_*\", \"space_group\": {\"number\": %d, \"symbol\": ", sg);
    put_json_str(b, symbol.text, symbol.length);
    put(b, ", \"ml_symbol\": ");
    put_json_str(b, ml.text, ml.length);
    put(b, ", \"crystal_class\": ");
    put_json_str(b, crystal.text, crystal.length);
    put(b, ", \"lattice\": %d, \"lattice_type\": ", lattice);
    put_json_str(b, lattice_type.text, lattice_type.length);
    put(b, ", \"point_group\": %d, \"point_group_order\": %lld}, \"kpoints\": [", point_group,
        (long long)section_int(&T.pg_order, point_group - 1));
    for (int i = 0; i < n; i++) {
        if (i) put(b, ", ");
        if (put_slot(b, sg, order[i], lattice, extra[i])) return 0;
    }
    put(b, "]}");
    return b->failed ? 0 : b->data;
}

/* ------------------------------------------------------------ Source <-> visible K parameters */

static section_view conv2ml;
static int conv2ml_ready;

/* source_kparam_from_visible on an evaluated visible vector: 1 with out[4],
 * 0 for none. */
int k_source_kparam_of_visible(int sg, int kslot, const int64_t *visible_words, int64_t *out)
{
    if (load()) return -1;
    if (!conv2ml_ready) {
        if (section_find("little.little_k_star_conv2ml_pointer", &conv2ml)) return -1;
        conv2ml_ready = 1;
    }
    int lattice = geo_sg_lattice(sg);
    if (section_int(&T.k_dim, (int64_t)(lattice - 1) * 27 + kslot - 1) <= 0) return 0;
    int64_t pointer = section_int(&conv2ml, (int64_t)(sg - 1) * 27 + kslot - 1);
    if (pointer <= 0) return 0;
    rat_overflow = 0;
    rat visible[3];
    for (int axis = 0; axis < 3; axis++) visible[axis] = rat_make(visible_words[2 * axis], visible_words[2 * axis + 1]);
    int64_t raw[16];
    for (int i = 0; i < 16; i++) raw[i] = section_int(&T.star, 16 * (pointer - 1) + i);
    int64_t den = raw[15] ? raw[15] : 1;
    rat converted[3];
    i128 common = 1;
    for (int col = 0; col < 3; col++) {
        rat value = rat_make(raw[12 + col], den);
        for (int row = 0; row < 3; row++) value = rat_add(value, rat_mul(visible[row], rat_make(raw[4 * row + col], den)));
        converted[col] = value;
        common = common / rat_gcd(common, value.d) * value.d;
    }
    for (int col = 0; col < 3; col++) out[col] = (int64_t)(converted[col].n * (common / converted[col].d));
    out[3] = (int64_t)common;
    return rat_overflow ? -2 : 1;
}

/* source_k_little_order: the little-group order at one Source K-family point. */
int k_source_little_order(int sg, int kslot, int count, const int64_t *parameter_words, int *order_out)
{
    if (load()) return -1;
    int lattice = geo_sg_lattice(sg);
    int64_t lattice_slot = (int64_t)(lattice - 1) * 27 + kslot - 1;
    int dim = (int)section_int(&T.k_dim, lattice_slot);
    if (count != dim) return -3;
    rat_overflow = 0;
    rat basis[4][3];
    slot_records(lattice_slot, basis);
    rat vector[3];
    for (int axis = 0; axis < 3; axis++) {
        vector[axis] = basis[0][axis];
        for (int i = 0; i < dim; i++)
            vector[axis] = rat_add(vector[axis], rat_mul(rat_make(parameter_words[2 * i], parameter_words[2 * i + 1]), basis[i + 1][axis]));
    }
    int64_t records[GEO_MAX_OPS][5];
    int nrecords = geo_sg_records(sg, records), ops[GEO_MAX_OPS], nops = 0;
    for (int i = 0; i < nrecords; i++) {
        int op = (int)records[i][4], seen = 0;
        for (int k = 0; k < nops && !seen; k++) seen = ops[k] == op;
        if (!seen) ops[nops++] = op;
    }
    int order = 0;
    for (int i = 0; i < nops; i++) {
        rat columns[3][3], direct[3][3], inverse[3][3];
        for (int unit = 0; unit < 3; unit++) {
            rat e[3] = {rat_int(unit == 0), rat_int(unit == 1), rat_int(unit == 2)};
            geo_vrot(sg, ops[i], e, columns[unit]);
        }
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) direct[row][col] = columns[col][row];
        if (rat_inverse3((const rat (*)[3])direct, inverse)) return -2;
        int integral = 1;
        for (int row = 0; row < 3 && integral; row++) {
            rat rotated = rat_int(0);
            for (int col = 0; col < 3; col++) rotated = rat_add(rotated, rat_mul(inverse[col][row], vector[col]));
            integral = rat_sub(rotated, vector[row]).d == 1;
        }
        order += integral;
    }
    *order_out = order;
    return rat_overflow ? -2 : 0;
}

/* ------------------------------------------------------------ display texts for the isotropy catalog */

static buffer display_text;

/* _k_slot_info(...)["display_kvector"] */
const char *k_display_kvector(int sg, int kslot)
{
    free(display_text.data);
    memset(&display_text, 0, sizeof display_text);
    if (load()) return 0;
    int lattice = geo_sg_lattice(sg);
    int64_t lattice_slot = (int64_t)(lattice - 1) * 27 + kslot - 1;
    int dim = (int)section_int(&T.k_dim, lattice_slot);
    rat records[4][3];
    slot_records(lattice_slot, records);
    kc_form form;
    if (display_form_of(sg, lattice, strip(&T.k_label, lattice_slot), (const rat (*)[3])records, dim + 1, &form)) return 0;
    put_display_vector(&display_text, &form);
    put(&display_text, "%s", "");
    return display_text.failed ? 0 : display_text.data;
}

/* _display_star_vectors(sg, kslot, dim) as a JSON list of strings. */
const char *k_display_star_json(int sg, int kslot, int dim)
{
    free(display_text.data);
    memset(&display_text, 0, sizeof display_text);
    if (load()) return 0;
    buffer *b = &display_text;
    int lattice = geo_sg_lattice(sg);
    int64_t lattice_slot = (int64_t)(lattice - 1) * 27 + kslot - 1, sg_slot = (int64_t)(sg - 1) * 27 + kslot - 1;
    int count = (int)section_int(&T.star_count, sg_slot);
    int64_t pointer = section_int(&T.star_pointer, sg_slot);
    str label = strip(&T.k_label, lattice_slot);
    put(b, "[");
    int written = 0;
    if (count > 0 && pointer > 0) {
        for (int arm = 0; arm < count; arm++) {
            rat records[4][3];
            int n;
            if (dim == 0) {
                int64_t base = 4 * (pointer - 1 + arm);
                if (base + 4 > section_count(&T.star_special)) break;
                int64_t den = section_int(&T.star_special, base + 3);
                if (!den) den = 1;
                for (int axis = 0; axis < 3; axis++) records[0][axis] = rat_make(section_int(&T.star_special, base + axis), den);
                n = 1;
            } else {
                int64_t base = 16 * (pointer - 1 + arm);
                if (base + 16 > section_count(&T.star)) break;
                for (int piece = 0; piece < 4; piece++) {
                    int64_t den = section_int(&T.star, base + 4 * piece + 3);
                    for (int axis = 0; axis < 3; axis++)
                        records[piece][axis] = den ? rat_make(section_int(&T.star, base + 4 * piece + axis), den) : rat_int(0);
                }
                n = dim + 1;
            }
            kc_form form;
            if (display_form_of(sg, lattice, label, (const rat (*)[3])records, n, &form)) return 0;
            buffer text = {0};
            if (dim == 0 && lattice == 3 && str_is(label, "L") && count == 2 && arm == 1)
                put(&text, "(-1/2,1/2,1/2)");
            else
                put_display_vector(&text, &form);
            if (text.failed) return 0;
            if (written) put(b, ", ");
            put_json_str(b, text.data, (int64_t)text.used);
            free(text.data);
            written++;
        }
    }
    put(b, "]");
    return b->failed ? 0 : b->data;
}

/* ------------------------------------------------------------ numeric display stars */

static void put_point(buffer *b, const rat *v)
{
    put(b, "(");
    for (int axis = 0; axis < 3; axis++) {
        if (axis) put(b, ",");
        put_frac(b, v[axis]);
    }
    put(b, ")");
}

static void star_arm_records(int64_t pointer, int arm, int dim, rat (*out)[3])
{
    for (int piece = 0; piece < 4; piece++) {
        int64_t base = 16 * (pointer - 1 + arm) + 4 * piece, den = section_int(&T.star, base + 3);
        for (int axis = 0; axis < 3; axis++) out[piece][axis] = den ? rat_make(section_int(&T.star, base + axis), den) : rat_int(0);
    }
    (void)dim;
}

/* isotropy.catalog._solve_star_parameter_values: 1 with values, 0 for None. */
static int solve_star_parameters(const rat (*records)[3], int count, const rat *visible, rat *values)
{
    int dim = count - 1;
    if (dim <= 0) return 1;
    rat eq[3][4];
    for (int axis = 0; axis < 3; axis++) {
        for (int p = 0; p < dim; p++) eq[axis][p] = records[p + 1][axis];
        eq[axis][dim] = rat_sub(visible[axis], records[0][axis]);
    }
    int row = 0, pivots = 0;
    for (int col = 0; col < dim; col++) {
        int pivot = -1;
        for (int c = row; c < 3 && pivot < 0; c++)
            if (!rat_zero(eq[c][col])) pivot = c;
        if (pivot < 0) continue;
        rat swap[4];
        memcpy(swap, eq[row], sizeof swap);
        memcpy(eq[row], eq[pivot], sizeof swap);
        memcpy(eq[pivot], swap, sizeof swap);
        rat scale = eq[row][col];
        for (int i = 0; i <= dim; i++) eq[row][i] = rat_div(eq[row][i], scale);
        for (int other = 0; other < 3; other++) {
            if (other == row) continue;
            rat factor = eq[other][col];
            if (rat_zero(factor)) continue;
            for (int i = 0; i <= dim; i++) eq[other][i] = rat_sub(eq[other][i], rat_mul(factor, eq[row][i]));
        }
        pivots++;
        row++;
        if (row == 3) break;
    }
    for (int r = 0; r < 3; r++) {
        int any = 0;
        for (int i = 0; i < dim; i++) any |= !rat_zero(eq[r][i]);
        if (!any && !rat_zero(eq[r][dim])) return 0;
    }
    if (pivots != dim) return 0;
    for (int i = 0; i < dim; i++) values[i] = rat_int(0);
    for (int r = 0; r < 3; r++)
        for (int i = 0; i < dim; i++)
            if (!rat_zero(eq[r][i])) {
                values[i] = eq[r][dim];
                break;
            }
    return 1;
}

int k_source_kparam_of_visible(int sg, int kslot, const int64_t *visible_words, int64_t *out);

/* isotropy.catalog._numeric_display_star_vectors as a JSON list (empty for
 * no star).  visible: the request's visible k under parent_sg's display map;
 * source_visible: the same under the gid's own space group, for the Source
 * kparam fallback. */
const char *k_numeric_star_json(int parent_sg, int gid, int has_visible, const int64_t *visible_words, int has_source_visible,
                                const int64_t *source_visible_words)
{
    free(display_text.data);
    memset(&display_text, 0, sizeof display_text);
    if (load()) return 0;
    buffer *b = &display_text;
    rat_overflow = 0;
    int kslot = (int)section_int(&T.irr_k, gid - 1), lattice = geo_sg_lattice(parent_sg);
    int64_t lattice_slot = (int64_t)(lattice - 1) * 27 + kslot - 1, sg_slot = (int64_t)(parent_sg - 1) * 27 + kslot - 1;
    int dim = (int)section_int(&T.k_dim, lattice_slot);
    int count = (int)section_int(&T.star_count, sg_slot);
    int64_t pointer = section_int(&T.star_pointer, sg_slot);
    put(b, "[");
    if (dim <= 0 || count <= 0 || pointer <= 0 || 16 * pointer > section_count(&T.star)) {
        put(b, "]");
        return b->data;
    }
    str label = strip(&T.k_label, lattice_slot), symbol = strip(&T.label_bc, parent_sg - 1);
    rat source[4][3], first[4][3];
    star_arm_records(pointer, 0, dim, source);
    if (kc_display_records(parent_sg, lattice, symbol.text, symbol.length, label.text, label.length, (const rat (*)[3])source, dim + 1,
                           first))
        return 0;
    rat visible[3];
    for (int axis = 0; axis < 3 && has_visible; axis++) visible[axis] = rat_make(visible_words[2 * axis], visible_words[2 * axis + 1]);
    if (lattice == 3 && str_is(label, "GP") && has_visible && count == 2) {
        rat sign = rat_int(-1);
        if (16 * (pointer + 1) <= section_count(&T.star)) {
            rat second_source[4][3];
            star_arm_records(pointer, 1, dim, second_source);
            if (dim + 1 > 1) sign = second_source[1][0];
        }
        rat second[3] = {visible[0], rat_neg(visible[1]), visible[2]};
        if (!(rat_cmp(sign, rat_int(0)) > 0)) {
            second[0] = rat_neg(visible[0]);
            second[1] = visible[1];
            second[2] = rat_neg(visible[2]);
        }
        put(b, "\"");
        put_point(b, visible);
        put(b, "\", \"");
        put_point(b, second);
        put(b, "\"]");
        return b->failed ? 0 : b->data;
    }
    rat values[3];
    int have = has_visible && solve_star_parameters((const rat (*)[3])first, dim + 1, visible, values);
    if (!have && has_visible && has_source_visible) {
        int64_t kp[4];
        int sg_gid = (int)section_int(&T.irr_sg, gid - 1);
        int found = k_source_kparam_of_visible(sg_gid, kslot, source_visible_words, kp);
        if (found < 0) return 0;
        if (found) {
            rat candidate[3];
            for (int i = 0; i < dim; i++) candidate[i] = rat_make(kp[i], kp[3]);
            int integral = 1;
            for (int axis = 0; axis < 3; axis++) {
                rat represented = first[0][axis];
                for (int i = 0; i < dim; i++) represented = rat_add(represented, rat_mul(first[i + 1][axis], candidate[i]));
                integral &= rat_sub(represented, visible[axis]).d == 1;
            }
            if (integral) {
                memcpy(values, candidate, sizeof candidate);
                have = 1;
            }
        }
    }
    if (have) {
        for (int arm = 0; arm < count; arm++) {
            if (16 * (pointer + arm) > section_count(&T.star)) break;
            rat arm_source[4][3], records[4][3];
            star_arm_records(pointer, arm, dim, arm_source);
            if (kc_display_records(parent_sg, lattice, symbol.text, symbol.length, label.text, label.length,
                                   (const rat (*)[3])arm_source, dim + 1, records))
                return 0;
            rat point[3];
            for (int axis = 0; axis < 3; axis++) {
                point[axis] = records[0][axis];
                for (int i = 0; i < dim; i++) point[axis] = rat_add(point[axis], rat_mul(records[i + 1][axis], values[i]));
            }
            if (arm) put(b, ", ");
            put(b, "\"");
            put_point(b, point);
            put(b, "\"");
        }
    }
    put(b, "]");
    return b->failed || rat_overflow ? 0 : b->data;
}
