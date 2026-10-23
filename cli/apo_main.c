/* apo: the command line in one executable. Arguments and case files are read
 * here; the result comes from cli_run_json. */

#include <errno.h>
#include <limits.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

#include "platform.h"
#include "json.h"
#include "pyconv.h"
#include "rational.h"

int tables_open(const char *path);
const char *cli_run_json(const char *text, int64_t length);
int cli_run_status(void);
int64_t cli_run_length(void);
const char *apo_error(int kind, int status, const char *message, jbuf *text);
int apo_serve(const char *host, int port, int open_browser, const char *frontend_dir, char *error, size_t error_size);
int apo_post(const char *server, const char *path, const jbuf *body, jbuf *reply, char *base, size_t base_size, char *error, size_t error_size);

static const char *command = "";

/* ------------------------------------------------------------ errors */

static void usage_error(const char *message)
{
    fprintf(stderr, "apo %s: error: %s\n", command, message);
    exit(2);
}

static void usage_errorf(const char *format, const char *a, const char *b)
{
    char text[8192];
    snprintf(text, sizeof text, format, a, b);
    usage_error(text);
}

static void fail_buf(jbuf *message)
{
    fprintf(stderr, "apo %s: error: %.*s\n", command, (int)message->used, message->data ? message->data : "");
    exit(2);
}

/* ------------------------------------------------------------ Python texts */

static void put_repr(jbuf *b, const char *s, size_t n) { py_repr_text(b, s, (int64_t)n); }

static int is_space(char c) { return c == ' ' || (c >= '\t' && c <= '\r') || (c >= 0x1c && c <= 0x1f); }

typedef struct {
    const char *s;
    size_t n;
} text;

static text strip(const char *s, size_t n)
{
    while (n && is_space(*s)) {
        s++;
        n--;
    }
    while (n && is_space(s[n - 1])) n--;
    return (text){s, n};
}

/* case_input._text: the stripped string or the field error */
static text text_field(const jvalue *v, const char *field)
{
    if (jis(v, J_STRING)) {
        text t = strip(v->u.s, (size_t)v->count);
        if (t.n) return t;
    }
    usage_errorf("%s must be a non-empty string", field, 0);
    return (text){0, 0};
}

static text text_arg(const char *s, const char *field)
{
    text t = s ? strip(s, strlen(s)) : (text){0, 0};
    if (!t.n) usage_errorf("%s must be a non-empty string", field, 0);
    return t;
}

/* case_input._rational: str(Fraction(str(value).strip())) */
static void rational_text(const char *s, size_t n, const char *field, jbuf *out)
{
    text t = strip(s, n);
    rat r;
    rat_overflow = 0;
    if (!t.n || py_fraction_of_text(t.s, (int64_t)t.n, &r) || rat_overflow) usage_errorf("%s must be rational", field, 0);
    char value[96];
    py_rat_text(r, value);
    jb_lit(out, value);
}

static void rational_value(const jvalue *v, const char *field, jbuf *out)
{
    if (jis(v, J_TRUE) || jis(v, J_FALSE) || !v || jis(v, J_NULL) || jis(v, J_ARRAY) || jis(v, J_OBJECT))
        usage_errorf("%s must be rational", field, 0);
    jbuf s = {0};
    jb_reset(&s);
    py_str_of(&s, v);
    rational_text(s.data ? s.data : "", s.used, field, out);
    jb_free(&s);
}

/* ------------------------------------------------------------ paths */

static int is_absolute(const char *path)
{
#ifdef _WIN32
    if (((path[0] >= 'A' && path[0] <= 'Z') || (path[0] >= 'a' && path[0] <= 'z')) && path[1] == ':') return 1;
    if (path[0] == '\\') return 1;
#endif
    return path[0] == '/';
}

static void join(char *out, size_t size, const char *a, const char *b)
{
    if (is_absolute(b)) snprintf(out, size, "%s", b);
    else snprintf(out, size, "%s/%s", a, b);
}

/* Path.expanduser().resolve() without requiring the path to exist */
static void resolve(const char *raw, const char *base, char *out)
{
    char path[PATH_MAX * 2];
    const char *home = getenv("HOME");
#ifdef _WIN32
    if (!home) home = getenv("USERPROFILE");
#endif
    if (raw[0] == '~' && (raw[1] == '/' || !raw[1])) snprintf(path, sizeof path, "%s%s", home ? home : "", raw + 1);
    else join(path, sizeof path, base, raw);
#ifdef _WIN32
    if (!_fullpath(out, path, PATH_MAX)) snprintf(out, PATH_MAX, "%s", path);
    for (char *c = out; *c; c++)
        if (*c == '\\') *c = '/';
    return;
#else
    if (realpath(path, out)) return;
    char parts[PATH_MAX * 2], *stack[512];
    int depth = 0;
    snprintf(parts, sizeof parts, "%s", path);
    for (char *p = strtok(parts, "/"); p; p = strtok(0, "/")) {
        if (!strcmp(p, ".")) continue;
        if (!strcmp(p, "..")) {
            if (depth) depth--;
            continue;
        }
        if (depth < 512) stack[depth++] = p;
    }
    char lexical[PATH_MAX * 2] = "";
    for (int i = 0; i < depth; i++) {
        strcat(lexical, "/");
        strcat(lexical, stack[i]);
    }
    if (!depth) strcpy(lexical, "/");
    for (int keep = depth; keep > 0; keep--) {
        char prefix[PATH_MAX * 2] = "", real[PATH_MAX];
        for (int i = 0; i < keep; i++) {
            strcat(prefix, "/");
            strcat(prefix, stack[i]);
        }
        if (plat_realpath(prefix, real)) {
            snprintf(out, PATH_MAX, "%s", real);
            for (int i = keep; i < depth; i++) {
                strncat(out, "/", PATH_MAX - strlen(out) - 1);
                strncat(out, stack[i], PATH_MAX - strlen(out) - 1);
            }
            return;
        }
    }
    snprintf(out, PATH_MAX, "%s", lexical);
#endif
}

static int is_file(const char *path)
{
    struct stat st;
    return !stat(path, &st) && S_ISREG(st.st_mode);
}

static int parent_of(char *path)
{
    char *slash = strrchr(path, '/');
    if (!slash || !strcmp(path, "/")) return 0;
    if (slash == path) path[1] = 0;
    else *slash = 0;
    return 1;
}

static char exe_dir[PATH_MAX];
static char cwd[PATH_MAX];

/* case_input._cif_asset_path */
static void cif_asset(const jvalue *v, const char *base, char *out)
{
    text t = text_field(v, "cif");
    char id[16];
    int ok = t.n == 12;
    for (size_t i = 0; ok && i < 12; i++) {
        char c = t.s[i];
        if (c >= 'A' && c <= 'F') c = (char)(c - 'A' + 'a');
        id[i] = c;
        ok = (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
    }
    if (!ok) usage_error("cif must be a 12-character content ID");
    id[12] = 0;
    const char *starts[3] = {base, cwd, exe_dir};
    for (int s = 0; s < 3; s++) {
        char root[PATH_MAX];
        resolve(starts[s], "/", root);
        for (;;) {
            static const char *const homes[] = {"Assets/cif", "assets/cif", "Assets/magndata/cif", "assets/magndata/cif"};
            for (int h = 0; h < 4; h++) {
                char candidate[PATH_MAX * 2];
                snprintf(candidate, sizeof candidate, "%s/%s/%s.cif", strcmp(root, "/") ? root : "", homes[h], id);
                if (is_file(candidate)) {
                    resolve(candidate, "/", out);
                    return;
                }
            }
            if (!parent_of(root)) break;
        }
    }
    char message[128];
    snprintf(message, sizeof message, "CIF asset does not exist: Assets/cif/%s.cif", id);
    usage_error(message);
}

/* ------------------------------------------------------------ the case */

typedef struct {
    char structure[PATH_MAX];
    jbuf k;     /* JSON array of {"label", "values", "params", "ir"} */
    jbuf sites; /* JSON array of [name, displacive, magnetic] */
    int nslots, strain, nsites, any_site_mode;
    char **slot_labels;
    jbuf opd; /* JSON value */
} apo_case;

static void case_init(apo_case *c)
{
    memset(c, 0, sizeof *c);
    jb_reset(&c->k);
    jb_reset(&c->sites);
    jb_reset(&c->opd);
    jb_lit(&c->opd, "null");
}

static void push_label(apo_case *c, text label)
{
    c->slot_labels = realloc(c->slot_labels, sizeof(char *) * (size_t)(c->nslots + 1));
    char *copy = malloc(label.n + 1);
    memcpy(copy, label.s, label.n);
    copy[label.n] = 0;
    c->slot_labels[c->nslots++] = copy;
}

static int slot_is_gamma(const char *label) { return (label[0] == 'g' || label[0] == 'G') && (label[1] == 'm' || label[1] == 'M') && !label[2]; }

static int key_order(const void *a, const void *b)
{
    const jvalue *x = *(const jvalue *const *)a, *y = *(const jvalue *const *)b;
    int64_t n = x->count < y->count ? x->count : y->count;
    int c = memcmp(x->u.s, y->u.s, (size_t)n);
    return c ? c : (x->count > y->count) - (x->count < y->count);
}

static void put_text_json(jbuf *b, text t) { jb_string(b, t.s, t.n); }

/* case_input._json_slot */
static void json_slot(const jvalue *raw, int ordinal, apo_case *c)
{
    char field[64];
    if (!jis(raw, J_OBJECT)) {
        snprintf(field, sizeof field, "k[%d] must be an object", ordinal);
        usage_error(field);
    }
    int has_label = 0, other = 0;
    for (int64_t i = 0; i < raw->count; i++) {
        const jvalue *key = raw->u.members + 2 * i;
        if (jstr_is(key, "label")) has_label = 1;
        else if (!jstr_is(key, "params") && !jstr_is(key, "ir")) other = 1;
    }
    if (!has_label || other) {
        snprintf(field, sizeof field, "k[%d] fields must be drawn from ['ir', 'label', 'params']", ordinal);
        usage_error(field);
    }
    const jvalue *params = jget(raw, "params");
    if (params && !jis(params, J_NULL) && !jis(params, J_OBJECT)) {
        snprintf(field, sizeof field, "k[%d].params must be an object", ordinal);
        usage_error(field);
    }
    jbuf pairs = {0};
    jb_reset(&pairs);
    jb_lit(&pairs, "[");
    if (jis(params, J_OBJECT) && params->count) {
        const jvalue **keys = malloc(sizeof *keys * (size_t)params->count);
        for (int64_t i = 0; i < params->count; i++) keys[i] = params->u.members + 2 * i;
        qsort(keys, (size_t)params->count, sizeof *keys, key_order);
        for (int64_t i = 0; i < params->count; i++) {
            char name_field[160];
            snprintf(name_field, sizeof name_field, "k[%d].params.%.*s", ordinal, (int)keys[i]->count, keys[i]->u.s);
            if (i) jb_lit(&pairs, ", ");
            jb_lit(&pairs, "[");
            jb_value(&pairs, keys[i]);
            jb_lit(&pairs, ", \"");
            rational_value(keys[i] + 1, name_field, &pairs);
            jb_lit(&pairs, "\"]");
        }
        free(keys);
    }
    jb_lit(&pairs, "]");
    snprintf(field, sizeof field, "k[%d].label", ordinal);
    text label = text_field(jget(raw, "label"), field);
    const jvalue *ir = jget(raw, "ir");
    text irrep = {0, 0};
    if (ir && !jis(ir, J_NULL)) {
        snprintf(field, sizeof field, "k[%d].ir", ordinal);
        irrep = text_field(ir, field);
    }
    if (c->nslots) jb_lit(&c->k, ", ");
    jb_lit(&c->k, "{\"label\": ");
    put_text_json(&c->k, label);
    jb_lit(&c->k, ", \"values\": [], \"params\": ");
    jb_raw(&c->k, pairs.data, pairs.used);
    jb_lit(&c->k, ", \"ir\": ");
    if (irrep.s) put_text_json(&c->k, irrep);
    else jb_lit(&c->k, "null");
    jb_lit(&c->k, "}");
    jb_free(&pairs);
    push_label(c, label);
}

/* case_input.case_from_json */
static void case_from_json(const jvalue *payload, const char *base, apo_case *c)
{
    case_init(c);
    if (!jis(payload, J_OBJECT)) usage_error("case JSON must be an object");
    static const char *const allowed[] = {"cif", "structure", "k", "sites", "strain", "opd"};
    int structure_fields = 0;
    for (int64_t i = 0; i < payload->count; i++) {
        const jvalue *key = payload->u.members + 2 * i;
        int ok = 0;
        for (int a = 0; a < 6 && !ok; a++) ok = jstr_is(key, allowed[a]);
        if (!ok) usage_error("case fields must be drawn from ['cif', 'k', 'opd', 'sites', 'strain', 'structure']");
        structure_fields += jstr_is(key, "cif") || jstr_is(key, "structure");
    }
    if (structure_fields != 1) usage_error("case must contain exactly one of cif or structure");
    const jvalue *raw_k = jget(payload, "k");
    if (raw_k && (!jis(raw_k, J_ARRAY) || raw_k->count > 4)) usage_error("k must contain zero to four ordered slots");
    const jvalue *raw_sites = jget(payload, "sites");
    if (raw_sites && !jis(raw_sites, J_NULL) && !jis(raw_sites, J_OBJECT)) usage_error("sites must be an object");
    jb_lit(&c->sites, "[");
    for (int64_t i = 0; jis(raw_sites, J_OBJECT) && i < raw_sites->count; i++) {
        const jvalue *name = raw_sites->u.members + 2 * i, *flags = name + 1;
        text element = text_field(name, "sites element");
        jbuf message = {0};
        jb_reset(&message);
        const jvalue *o = jis(flags, J_OBJECT) ? jget(flags, "occupancy") : 0;
        int exact = jis(flags, J_OBJECT) && flags->count == 2 + (o != 0) && jget(flags, "displacive") && jget(flags, "magnetic");
        if (!exact) {
            jb_lit(&message, "sites[");
            put_repr(&message, element.s, element.n);
            jb_lit(&message, "] must contain displacive and magnetic, and optionally occupancy");
            fail_buf(&message);
        }
        const jvalue *d = jget(flags, "displacive"), *m = jget(flags, "magnetic");
        int bools = (jis(d, J_TRUE) || jis(d, J_FALSE)) && (jis(m, J_TRUE) || jis(m, J_FALSE)) && (!o || jis(o, J_TRUE) || jis(o, J_FALSE));
        if (!bools) {
            jb_lit(&message, "sites[");
            put_repr(&message, element.s, element.n);
            jb_lit(&message, "] flags must be boolean");
            fail_buf(&message);
        }
        jb_free(&message);
        if (jis(d, J_TRUE) || jis(m, J_TRUE) || jis(o, J_TRUE)) {
            if (c->nsites++) jb_lit(&c->sites, ", ");
            jb_lit(&c->sites, "[");
            put_text_json(&c->sites, element);
            jb_lit(&c->sites, jis(d, J_TRUE) ? ", true" : ", false");
            jb_lit(&c->sites, jis(m, J_TRUE) ? ", true" : ", false");
            jb_lit(&c->sites, jis(o, J_TRUE) ? ", true]" : ", false]");
            c->any_site_mode = 1;
        }
    }
    jb_lit(&c->sites, "]");
    const jvalue *strain = jget(payload, "strain");
    if (strain && !jis(strain, J_TRUE) && !jis(strain, J_FALSE)) usage_error("strain must be boolean");
    c->strain = jis(strain, J_TRUE);
    if (jget(payload, "cif")) {
        cif_asset(jget(payload, "cif"), base, c->structure);
    } else {
        text raw = text_field(jget(payload, "structure"), "structure");
        char path[PATH_MAX * 2];
        snprintf(path, sizeof path, "%.*s", (int)raw.n, raw.s);
        resolve(path, base, c->structure);
        if (!is_file(c->structure)) usage_errorf("structure CIF does not exist: %s", c->structure, 0);
    }
    jb_lit(&c->k, "[");
    for (int64_t i = 0; jis(raw_k, J_ARRAY) && i < raw_k->count; i++) json_slot(jat(raw_k, i), (int)i + 1, c);
    jb_lit(&c->k, "]");
    const jvalue *opd = jget(payload, "opd");
    if (opd && !jis(opd, J_NULL)) {
        text t = text_field(opd, "opd");
        jb_reset(&c->opd);
        put_text_json(&c->opd, t);
    }
}

static char *read_all(FILE *f, size_t *size)
{
    size_t cap = 1 << 16, used = 0;
    char *data = malloc(cap);
    for (;;) {
        size_t got = fread(data + used, 1, cap - used, f);
        used += got;
        if (used < cap) break;
        cap *= 2;
        data = realloc(data, cap);
    }
    *size = used;
    return data;
}

/* cli._load_cli_case */
static void load_case(const char *raw, apo_case *c)
{
    size_t size;
    char *data;
    char base[PATH_MAX];
    int from_stdin = !strcmp(raw, "-");
    if (from_stdin) {
        data = read_all(stdin, &size);
        snprintf(base, sizeof base, "%s", cwd);
    } else {
        char path[PATH_MAX];
        resolve(raw, cwd, path);
        if (!is_file(path)) usage_errorf("case file does not exist: %s", path, 0);
        FILE *f = fopen(path, "rb");
        if (!f) usage_errorf("case file does not exist: %s", path, 0);
        data = read_all(f, &size);
        fclose(f);
        size_t n = strlen(path);
        if (n < 5 || strcasecmp(path + n - 5, ".json")) usage_error("case input must be JSON");
        snprintf(base, sizeof base, "%s", path);
        parent_of(base);
    }
    jdoc *doc = jdoc_parse(data, size);
    if (!doc) usage_error(from_stdin ? "invalid case JSON from stdin" : "invalid case JSON");
    case_from_json(jdoc_root(doc), base, c);
    jdoc_free(doc);
    free(data);
}

/* ------------------------------------------------------------ arguments */

extern int exact_certify;
extern long exact_certified[4];
extern long exact_walk_certified[8];

typedef struct {
    const char *structure, *case_path, *opd, *output, *format, *upto, *host, *server;
    int strain, full_state, indent, port, open_browser, certify;
    char **amplitudes;
    int namplitudes;
    int nk;
    char ***k;
    int *k_counts;
    char **irreps, **displacive, **magnetic, **occupancy;
    int nirreps, ndisplacive, nmagnetic, noccupancy;
    int has_k_option, has_irrep_option, has_case_option, has_opd_option, has_format, has_full_state, has_distortion;
    const char *sg;
} args_t;

static int option_like(const char *s) { return s[0] == '-' && s[1] && !(s[1] >= '0' && s[1] <= '9') && s[1] != '.'; }

static void push(char ***list, int *count, char *value)
{
    *list = realloc(*list, sizeof(char *) * (size_t)(*count + 1));
    (*list)[(*count)++] = value;
}

static void parse_args(int argc, char **argv, int start, args_t *a)
{
    memset(a, 0, sizeof *a);
    a->indent = 2;
    a->upto = "modes";
    a->host = "127.0.0.1";
    a->port = 8300;
    a->server = "http://127.0.0.1:8300";
    int serve = !strcmp(command, "serve"), show = !strcmp(command, "show");
    int positional = 0;
    for (int i = start; i < argc; i++) {
        char *arg = argv[i], *value = 0, *eq = 0;
        if (arg[0] == '-' && arg[1] == '-' && (eq = strchr(arg, '='))) {
            *eq = 0;
            value = eq + 1;
        }
#define TAKE()                                                                                                                      \
    do {                                                                                                                            \
        if (!value) {                                                                                                               \
            if (i + 1 >= argc || option_like(argv[i + 1])) usage_errorf("argument %s: expected one argument", arg, 0);              \
            value = argv[++i];                                                                                                      \
        }                                                                                                                           \
    } while (0)
#define TAKE_MANY(list, count)                                                                                                      \
    do {                                                                                                                            \
        int taken = 0;                                                                                                              \
        if (value) {                                                                                                                \
            push(list, count, value);                                                                                               \
            taken++;                                                                                                                \
        }                                                                                                                           \
        while (i + 1 < argc && !option_like(argv[i + 1])) {                                                                         \
            push(list, count, argv[++i]);                                                                                           \
            taken++;                                                                                                                \
        }                                                                                                                           \
        if (!taken) usage_errorf("argument %s: expected at least one argument", arg, 0);                                            \
    } while (0)
        int accepts_selection = !strcmp(command, "irreps") || !strcmp(command, "opds") || !strcmp(command, "modes");
        if (serve && !strcmp(arg, "--host")) {
            TAKE();
            a->host = value;
        } else if (serve && !strcmp(arg, "--port")) {
            TAKE();
            char *end;
            long n = strtol(value, &end, 10);
            while (*end == ' ' || (*end >= '\t' && *end <= '\r')) end++;
            if (*end || !*value || n < INT_MIN || n > INT_MAX) usage_errorf("argument --port: invalid int value: '%s'", value, 0);
            a->port = (int)n;
        } else if (serve && !strcmp(arg, "--open-browser")) {
            a->open_browser = 1;
        } else if (show && !strcmp(arg, "--server")) {
            TAKE();
            a->server = value;
        } else if (show && !strcmp(arg, "--amplitude")) {
            TAKE();
            push(&a->amplitudes, &a->namplitudes, value);
        } else if (!serve && (!strcmp(arg, "-o") || !strcmp(arg, "--output"))) {
            TAKE();
            a->output = value;
        } else if (!serve && !strcmp(arg, "--indent")) {
            TAKE();
            char *end;
            long n = strtol(value, &end, 10);
            if (*end || !*value) usage_errorf("argument --indent: invalid int value: '%s'", value, 0);
            a->indent = (int)n;
        } else if (!strcmp(arg, "--certify") && !show && !serve) {
            a->certify = 1;
        } else if (!strcmp(arg, "--full-state") && !show && !serve) {
            a->full_state = 1;
            a->has_full_state = 1;
        } else if (!strcmp(arg, "--format") && (!strcmp(command, "modes") || !strcmp(command, "run"))) {
            TAKE();
            if (strcmp(value, "json") && strcmp(value, "text"))
                usage_errorf("argument --format: invalid choice: '%s' (choose from json, text)", value, 0);
            a->format = value;
        } else if (!strcmp(arg, "--case") && (!strcmp(command, "modes") || !strcmp(command, "run") || show)) {
            TAKE();
            a->case_path = value;
            a->has_case_option = 1;
        } else if (!strcmp(arg, "--upto") && !strcmp(command, "run")) {
            TAKE();
            static const char *const stages[] = {"kpoints", "irreps", "opds", "modes", "spec"};
            int ok = 0;
            for (int s = 0; s < 5 && !ok; s++) ok = !strcmp(value, stages[s]);
            if (!ok) usage_errorf("argument --upto: invalid choice: '%s' (choose from kpoints, irreps, opds, modes)", value, 0);
            a->upto = value;
        } else if (!strcmp(arg, "--k") && accepts_selection) {
            a->k = realloc(a->k, sizeof(char **) * (size_t)(a->nk + 1));
            a->k_counts = realloc(a->k_counts, sizeof(int) * (size_t)(a->nk + 1));
            a->k[a->nk] = 0;
            a->k_counts[a->nk] = 0;
            TAKE_MANY(&a->k[a->nk], &a->k_counts[a->nk]);
            a->nk++;
            a->has_k_option = 1;
        } else if (!strcmp(arg, "--sg") && (!strcmp(command, "kpoints") || !strcmp(command, "irreps") || !strcmp(command, "opds"))) {
            TAKE();
            a->sg = value;
        } else if (!strcmp(arg, "--irrep") && (!strcmp(command, "opds") || !strcmp(command, "modes"))) {
            TAKE();
            push(&a->irreps, &a->nirreps, value);
            a->has_irrep_option = 1;
        } else if (!strcmp(arg, "--displacive") && accepts_selection) {
            TAKE_MANY(&a->displacive, &a->ndisplacive);
        } else if (!strcmp(arg, "--magnetic") && accepts_selection) {
            TAKE_MANY(&a->magnetic, &a->nmagnetic);
        } else if (!strcmp(arg, "--occupancy") && accepts_selection) {
            TAKE_MANY(&a->occupancy, &a->noccupancy);
        } else if (!strcmp(arg, "--strain") && accepts_selection) {
            a->strain = 1;
        } else if (!strcmp(arg, "--opd") && !strcmp(command, "modes")) {
            TAKE();
            a->opd = value;
            a->has_opd_option = 1;
        } else if (option_like(arg)) {
            if (eq) *eq = '=';
            fprintf(stderr, "apo: error: unrecognized arguments: %s\n", arg);
            exit(2);
        } else if (!positional && strcmp(command, "run") && !show && !serve) {
            a->structure = arg;
            positional = 1;
        } else {
            fprintf(stderr, "apo: error: unrecognized arguments: %s\n", arg);
            exit(2);
        }
#undef TAKE
#undef TAKE_MANY
    }
    if (!strcmp(command, "opds") && (!a->has_k_option || !a->has_irrep_option)) {
        if (!a->has_k_option && !a->has_irrep_option) usage_error("the following arguments are required: --k, --irrep");
        usage_errorf("the following arguments are required: %s", a->has_k_option ? "--irrep" : "--k", 0);
    }
    if ((!strcmp(command, "run") || show) && !a->case_path) usage_error("the following arguments are required: --case");
}

static int in_list(char **list, int n, const char *s)
{
    for (int i = 0; i < n; i++)
        if (!strcmp(list[i], s)) return 1;
    return 0;
}

static char sg_cif_path[PATH_MAX];

static void remove_sg_cif(void)
{
    if (sg_cif_path[0]) remove(sg_cif_path);
}

/* --sg N: the parent is space group N in its default setting, given to the pipeline as a CIF with one atom (label A) at
 * a general position, whose displacements and moments carry every irrep; k points, irreps and OPDs do not depend on it. */
static void write_sg_cif(args_t *a)
{
    char *end;
    long sg = strtol(a->sg, &end, 10);
    if (*end || !*a->sg || sg < 1 || sg > 230) usage_errorf("argument --sg: invalid space group number: '%s'", a->sg, 0);
    if (a->structure) usage_error("give either a structure or --sg, not both");
    const char *cell = sg <= 2     ? "5.1 5.7 6.3 81 86 79"
                       : sg <= 15  ? "5.1 5.7 6.3 90 97 90"
                       : sg <= 74  ? "5.1 5.7 6.3 90 90 90"
                       : sg <= 142 ? "5.1 5.1 6.3 90 90 90"
                       : sg <= 194 ? "5.1 5.1 6.3 90 90 120"
                                   : "5.1 5.1 5.1 90 90 90";
    double v[6];
    sscanf(cell, "%lf %lf %lf %lf %lf %lf", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]);
    const char *dir = getenv("TMPDIR");
    if (!dir || !*dir) dir = getenv("TEMP");
    if (!dir || !*dir) dir = "/tmp";
    snprintf(sg_cif_path, sizeof sg_cif_path, "%s/apo-sg%ld-%ld.cif", dir, sg, (long)getpid());
    FILE *f = fopen(sg_cif_path, "w");
    if (!f) usage_errorf("cannot write a temporary CIF in %s", dir, 0);
    fprintf(f,
            "data_sg%ld\n_symmetry_Int_Tables_number %ld\n_cell_length_a %g\n_cell_length_b %g\n_cell_length_c %g\n"
            "_cell_angle_alpha %g\n_cell_angle_beta %g\n_cell_angle_gamma %g\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
            "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nA C 0.1234 0.2345 0.3456\n",
            sg, sg, v[0], v[1], v[2], v[3], v[4], v[5]);
    fclose(f);
    atexit(remove_sg_cif);
    a->structure = sg_cif_path;
    static char *site = "A";
    if (!a->ndisplacive && !a->nmagnetic && !a->noccupancy && !a->strain) push(&a->displacive, &a->ndisplacive, site);
}

/* case_input.case_from_cli */
static void case_from_cli(const args_t *a, apo_case *c)
{
    case_init(c);
    resolve(a->structure, cwd, c->structure);
    if (a->nk > 4) usage_error("--k must contain at most four ordered slots");
    if (a->nirreps && a->nirreps != a->nk) usage_error("each --k must have one corresponding --irrep");
    jb_lit(&c->k, "[");
    for (int s = 0; s < a->nk; s++) {
        int ordinal = s + 1;
        char field[96];
        snprintf(field, sizeof field, "k[%d].label", ordinal);
        text label = text_arg(a->k_counts[s] ? a->k[s][0] : 0, field);
        int nvalues = a->k_counts[s] - 1, named = 0;
        for (int v = 0; v < nvalues; v++) named += strchr(a->k[s][v + 1], '=') != 0;
        if (named && named != nvalues) {
            snprintf(field, sizeof field, "k[%d] parameters must be entirely positional or NAME=VALUE", ordinal);
            usage_error(field);
        }
        jbuf values = {0}, params = {0};
        jb_reset(&values);
        jb_reset(&params);
        jb_lit(&values, "[");
        jb_lit(&params, "[");
        if (nvalues && named == nvalues) {
            for (int v = 0; v < nvalues; v++) {
                char *raw = a->k[s][v + 1], *eqs = strchr(raw, '=');
                snprintf(field, sizeof field, "k[%d] parameter name", ordinal);
                char name_copy[256];
                snprintf(name_copy, sizeof name_copy, "%.*s", (int)(eqs - raw), raw);
                text name = text_arg(name_copy, field);
                for (int w = 0; w < v; w++) {
                    char *other = a->k[s][w + 1], *oeq = strchr(other, '=');
                    char other_name[256];
                    snprintf(other_name, sizeof other_name, "%.*s", (int)(oeq - other), other);
                    text o = strip(other_name, strlen(other_name));
                    if (o.n == name.n && !memcmp(o.s, name.s, name.n)) {
                        jbuf message = {0};
                        jb_reset(&message);
                        char head[32];
                        snprintf(head, sizeof head, "k[%d] parameter ", ordinal);
                        jb_lit(&message, head);
                        put_repr(&message, name.s, name.n);
                        jb_lit(&message, " is duplicated");
                        fail_buf(&message);
                    }
                }
                char value_field[320];
                snprintf(value_field, sizeof value_field, "k[%d].params.%.*s", ordinal, (int)name.n, name.s);
                if (v) jb_lit(&params, ", ");
                jb_lit(&params, "[");
                put_text_json(&params, name);
                jb_lit(&params, ", \"");
                rational_text(eqs + 1, strlen(eqs + 1), value_field, &params);
                jb_lit(&params, "\"]");
            }
        } else {
            snprintf(field, sizeof field, "k[%d].values", ordinal);
            for (int v = 0; v < nvalues; v++) {
                if (v) jb_lit(&values, ", ");
                jb_lit(&values, "\"");
                rational_text(a->k[s][v + 1], strlen(a->k[s][v + 1]), field, &values);
                jb_lit(&values, "\"");
            }
        }
        jb_lit(&values, "]");
        jb_lit(&params, "]");
        if (s) jb_lit(&c->k, ", ");
        jb_lit(&c->k, "{\"label\": ");
        put_text_json(&c->k, label);
        jb_lit(&c->k, ", \"values\": ");
        jb_raw(&c->k, values.data, values.used);
        jb_lit(&c->k, ", \"params\": ");
        jb_raw(&c->k, params.data, params.used);
        jb_lit(&c->k, ", \"ir\": ");
        if (a->nirreps) jb_cstring(&c->k, a->irreps[s]);
        else jb_lit(&c->k, "null");
        jb_lit(&c->k, "}");
        jb_free(&values);
        jb_free(&params);
        push_label(c, label);
    }
    jb_lit(&c->k, "]");
    jb_lit(&c->sites, "[");
    char **seen = 0;
    int nseen = 0;
    for (int pass = 0; pass < 3; pass++) {
        char **list = pass == 2 ? a->occupancy : pass ? a->magnetic : a->displacive;
        int n = pass == 2 ? a->noccupancy : pass ? a->nmagnetic : a->ndisplacive;
        for (int i = 0; i < n; i++) {
            if (in_list(seen, nseen, list[i])) continue;
            push(&seen, &nseen, list[i]);
            if (c->nsites++) jb_lit(&c->sites, ", ");
            jb_lit(&c->sites, "[");
            jb_cstring(&c->sites, list[i]);
            jb_lit(&c->sites, in_list(a->displacive, a->ndisplacive, list[i]) ? ", true" : ", false");
            jb_lit(&c->sites, in_list(a->magnetic, a->nmagnetic, list[i]) ? ", true" : ", false");
            jb_lit(&c->sites, in_list(a->occupancy, a->noccupancy, list[i]) ? ", true]" : ", false]");
            c->any_site_mode = 1;
        }
    }
    free(seen);
    jb_lit(&c->sites, "]");
    c->strain = a->strain;
    if (a->opd) {
        jb_reset(&c->opd);
        jb_cstring(&c->opd, a->opd);
    }
}

/* cli._cli_case */
static void cli_case(const args_t *a, const char *stage, apo_case *c)
{
    if ((!strcmp(command, "modes") || !strcmp(command, "run")) && a->case_path) {
        if (a->structure || a->nk || a->nirreps || a->ndisplacive || a->nmagnetic || a->noccupancy || a->strain || a->opd)
            usage_error("--case cannot be combined with direct mode-selection arguments");
        load_case(a->case_path, c);
        int pure_strain_irreps = !strcmp(stage, "irreps") && c->strain && !c->any_site_mode;
        if (pure_strain_irreps)
            for (int s = 0; s < c->nslots; s++)
                if (!slot_is_gamma(c->slot_labels[s])) usage_error("strain-only irreps require k label GM or no k");
        if (!c->nslots && strcmp(stage, "kpoints") && !pure_strain_irreps) usage_errorf("%s case requires one to four ordered k slots", stage, 0);
        if (!strcmp(stage, "modes") && !strcmp(c->opd.data, "null")) usage_errorf("%s case requires an OPD", stage, 0);
        return;
    }
    if (!a->structure) {
        if (!strcmp(command, "modes")) usage_error("modes requires structure or --case");
        usage_errorf("%s requires structure", command, 0);
    }
    if (!strcmp(command, "modes") && (!a->nk || !a->nirreps || !a->opd || !*a->opd)) usage_error("direct modes requires --k, --irrep, and --opd");
    int pure_strain = a->strain && !a->ndisplacive && !a->nmagnetic && !a->noccupancy;
    if (!strcmp(command, "irreps") && !a->nk && !pure_strain) usage_error("irreps requires --k unless --strain is the only distortion");
    if (!strcmp(command, "irreps") && pure_strain && a->nk)
        for (int s = 0; s < a->nk; s++)
            if (a->k_counts[s] && !slot_is_gamma(a->k[s][0])) usage_error("strain-only irreps require --k GM or no --k");
    if ((!strcmp(command, "irreps") || !strcmp(command, "opds") || !strcmp(command, "modes")) && !a->ndisplacive && !a->nmagnetic && !a->noccupancy && !a->strain)
        usage_error("select at least one distortion with --displacive SITE, --magnetic SITE, --occupancy SITE, or --strain; SITE is an info.sites[].label or type "
                    "value");
    case_from_cli(a, c);
}

/* ------------------------------------------------------------ output */

static void write_result(const char *output, const char *data, size_t n)
{
    if (!output) {
        fwrite(data, 1, n, stdout);
        fflush(stdout);
        return;
    }
    char path[PATH_MAX];
    resolve(output, cwd, path);
    char dir[PATH_MAX];
    snprintf(dir, sizeof dir, "%s", path);
    parent_of(dir);
    char partial[PATH_MAX];
    snprintf(partial, sizeof partial, "%s", dir);
    for (char *p = partial + 1; *p; p++)
        if (*p == '/') {
            *p = 0;
            plat_mkdir(partial);
            *p = '/';
        }
    plat_mkdir(partial);
    FILE *f = fopen(path, "wb");
    if (!f) usage_errorf("[Errno %s] could not write output: '%s'", strerror(errno), path);
    fwrite(data, 1, n, f);
    fclose(f);
}

static int report_error(const char *result)
{
    {
        jdoc *doc = jdoc_parse(result, (size_t)cli_run_length());
        const jvalue *error = jget(jdoc_root(doc), "error");
        int64_t kind = 0, code = 0;
        py_int_of(jget(error, "kind"), &kind);
        py_int_of(jget(error, "status"), &code);
        const jvalue *message = jget(error, "message");
        char *text = plat_strndup(message ? message->u.s : "", message ? (size_t)message->count : 0);
        jbuf what = {0};
        jb_reset(&what);
        const char *type = apo_error((int)kind, (int)code, text, &what);
        static const char *const caught[] = {"ValueError", "KeyError", "IndexError", "OpdNotFoundError", "LinAlgError"};
        int is_caught = 0;
        for (int i = 0; i < 5 && !is_caught; i++) is_caught = !strcmp(type, caught[i]);
        if (is_caught) fprintf(stderr, "apo %s: error: %.*s\n", command, (int)what.used, what.data ? what.data : "");
        else fprintf(stderr, "%s: %.*s\n", type, (int)what.used, what.data ? what.data : "");
        jb_free(&what);
        free(text);
        jdoc_free(doc);
        return is_caught ? 2 : 1;
    }
}

static int run_request(jbuf *request, const char *output)
{
    const char *result = cli_run_json(request->data, (int64_t)request->used);
    if (!result) {
        fprintf(stderr, "MemoryError: CLI run ran out of memory\n");
        return 1;
    }
    int status = cli_run_status();
    if (status < 0) return report_error(result);
    write_result(output, result, (size_t)cli_run_length());
    return status;
}

/* ------------------------------------------------------------ GUI */

static void stop_serving(int signal_number)
{
    (void)signal_number;
    _exit(0);
}

static int serve_command(const args_t *a)
{
    char frontend[PATH_MAX * 2], error[512];
    const char *env = getenv("APO_FRONTEND");
    if (env) snprintf(frontend, sizeof frontend, "%s", env);
    else snprintf(frontend, sizeof frontend, "%s/../gui/web", exe_dir);
    signal(SIGINT, stop_serving);
    if (apo_serve(a->host, a->port, a->open_browser, frontend, error, sizeof error)) usage_error(error);
    return 0;
}

static text trimmed(const char *s, size_t n) { return strip(s, n); }

/* gui_launch.parse_definition_weights, as launch amplitude rows */
static void put_amplitudes(const args_t *a, jbuf *b)
{
    jb_lit(b, "[");
    for (int i = 0; i < a->namplitudes; i++) {
        const char *raw = a->amplitudes[i], *eq = strchr(raw, '=');
        text id = trimmed(raw, eq ? (size_t)(eq - raw) : strlen(raw));
        text value = eq ? trimmed(eq + 1, strlen(eq + 1)) : (text){"", 0};
        if (!eq || !id.n || !value.n) usage_error("weights must use DEFINITION_ID=VALUE");
        for (int k = 0; k < i; k++) {
            const char *other = a->amplitudes[k], *other_eq = strchr(other, '=');
            text other_id = trimmed(other, (size_t)(other_eq - other));
            if (other_id.n == id.n && !memcmp(other_id.s, id.s, id.n)) {
                jbuf message = {0};
                jb_reset(&message);
                jb_lit(&message, "duplicate weight for definition ");
                jb_raw(&message, id.s, id.n);
                fail_buf(&message);
            }
        }
        rat r;
        if (py_fraction_of_text(value.s, (int64_t)value.n, &r)) {
            jbuf message = {0};
            jb_reset(&message);
            jb_lit(&message, "weight for definition ");
            jb_raw(&message, id.s, id.n);
            jb_lit(&message, " must be an exact number");
            fail_buf(&message);
        }
        char exact[96];
        py_rat_text(r, exact);
        if (i) jb_lit(b, ", ");
        jb_lit(b, "{\"definition_id\": ");
        jb_string(b, id.s, id.n);
        jb_lit(b, ", \"value\": ");
        jb_cstring(b, exact);
        jb_lit(b, "}");
    }
    jb_lit(b, "]");
}

static void put_case(jbuf *b, const apo_case *c);

/* gui_launch.gui_launch_from_case and send_gui_launch */
static int show_command(const args_t *a)
{
    apo_case c;
    load_case(a->case_path, &c);
    jdoc *kdoc = jdoc_parse(c.k.data, c.k.used);
    const jvalue *slots = jdoc_root(kdoc);
    int nslots = jis(slots, J_ARRAY) ? (int)slots->count : 0, all_irreps = nslots > 0;
    for (int i = 0; i < nslots; i++) all_irreps = all_irreps && jtruthy(jget(jat(slots, i), "ir"));
    jdoc_free(kdoc);
    const char *stage = strcmp(c.opd.data, "null") ? "modes" : all_irreps ? "opds" : nslots ? "irreps" : "kpoints";
    int deep = !strcmp(stage, "opds") || !strcmp(stage, "modes");
    jbuf request = {0};
    jb_reset(&request);
    jb_lit(&request, "{\"call\": \"request\", ");
    put_case(&request, &c);
    jb_lit(&request, deep ? ", \"require_irreps\": true" : ", \"require_irreps\": false");
    jb_lit(&request, !strcmp(stage, "kpoints") ? ", \"allow_empty_distortion\": true}" : ", \"allow_empty_distortion\": false}");
    const char *result = cli_run_json(request.data, (int64_t)request.used);
    jb_free(&request);
    if (!result) {
        fprintf(stderr, "MemoryError: CLI run ran out of memory\n");
        return 1;
    }
    if (cli_run_status() < 0) return report_error(result);
    jdoc *rdoc = jdoc_parse(result, (size_t)cli_run_length());
    const jvalue *options = jdoc_root(rdoc);
    jbuf launch = {0};
    jb_reset(&launch);
    jb_lit(&launch, "{\"schema\": \"APOSTRUCT.gui.launch.v1\", \"stage\": ");
    jb_cstring(&launch, stage);
    const char *name = strrchr(c.structure, '/');
    jb_lit(&launch, ", \"source_name\": ");
    jb_cstring(&launch, name ? name + 1 : c.structure);
    jb_lit(&launch, ", \"state_request\": {");
    static const char *const fields[] = {"mode_set", "distortion_selection", "k_label", "k_index", "k_params",
                                         "irrep", "irrep_index", "orderparam", "orderparam_index", "selections"};
    for (int i = 0; i < 10; i++) {
        const jvalue *v = jget(options, fields[i]);
        if (!v || jis(v, J_NULL)) continue;
        jb_lit(&launch, "\"");
        jb_lit(&launch, fields[i]);
        jb_lit(&launch, "\": ");
        jb_value(&launch, v);
        jb_lit(&launch, ", ");
    }
    jdoc_free(rdoc);
    FILE *f = fopen(c.structure, "rb");
    if (!f) usage_errorf("could not read structure CIF: %s", c.structure, 0);
    size_t size;
    char *cif = read_all(f, &size);
    fclose(f);
    jb_lit(&launch, "\"cif\": ");
    jb_string(&launch, cif, size);
    free(cif);
    jb_lit(&launch, deep ? ", \"include_opd\": true" : ", \"include_opd\": false");
    jb_lit(&launch, !strcmp(stage, "modes") ? ", \"include_mode_details\": true}" : ", \"include_mode_details\": false}");
    if (a->namplitudes && strcmp(stage, "modes")) usage_error("--amplitude requires a case with opd");
    jb_lit(&launch, ", \"amplitudes\": ");
    put_amplitudes(a, &launch);
    jb_lit(&launch, "}");
    jbuf reply = {0};
    jb_reset(&reply);
    char base[1024], error[32768];
    if (apo_post(a->server, "/api/APOSTRUCT/gui_launch", &launch, &reply, base, sizeof base, error, sizeof error)) usage_error(error);
    jb_free(&launch);
    jdoc *doc = jdoc_parse(reply.data ? reply.data : "", reply.used);
    const jvalue *root = doc ? jdoc_root(doc) : 0;
    if (!jis(root, J_OBJECT) || !jis(jget(root, "revision"), J_INT)) usage_error("GUI server returned an invalid launch response");
    int n = (int)root->count + 1;
    const jvalue **keys = calloc((size_t)n, sizeof *keys);
    static const jvalue server_key = {J_STRING, 6, {.s = "server"}};
    int nkeys = 0;
    for (int64_t i = 0; i < root->count; i++)
        if (!jstr_is(root->u.members + 2 * i, "server")) keys[nkeys++] = root->u.members + 2 * i;
    keys[nkeys++] = &server_key;
    qsort(keys, (size_t)nkeys, sizeof *keys, key_order);
    jbuf out = {0};
    jb_reset(&out);
    int indent = a->indent > 0 ? a->indent : 0;
    jb_lit(&out, "{");
    for (int i = 0; i < nkeys; i++) {
        jb_lit(&out, i ? ",\n" : "\n");
        for (int k = 0; k < indent; k++) jb_lit(&out, " ");
        jb_value(&out, keys[i]);
        jb_lit(&out, ": ");
        if (keys[i] == &server_key) jb_cstring(&out, base);
        else jb_value(&out, keys[i] + 1);
    }
    jb_lit(&out, "\n}\n");
    write_result(a->output, out.data, out.used);
    free(keys);
    jdoc_free(doc);
    jb_free(&reply);
    jb_free(&out);
    return 0;
}

static void put_case(jbuf *b, const apo_case *c)
{
    jb_lit(b, "\"structure\": ");
    jb_cstring(b, c->structure);
    jb_lit(b, ", \"k\": ");
    jb_raw(b, c->k.data, c->k.used);
    jb_lit(b, ", \"sites\": ");
    jb_raw(b, c->sites.data, c->sites.used);
    jb_lit(b, c->strain ? ", \"strain\": true" : ", \"strain\": false");
    jb_lit(b, ", \"opd\": ");
    jb_raw(b, c->opd.data, c->opd.used);
}

static const char *usage_text =
    "usage: apo [-h] [--version] {info,kpoints,irreps,opds,modes,run,serve,show} ...\n\n"
    "Compute symmetry-adapted distortion modes from a parent CIF.\n\n"
    "commands:\n"
    "  info       show the normalized parent structure\n"
    "  kpoints    enumerate parent k points\n"
    "  irreps     enumerate irreps at selected k points\n"
    "  opds       enumerate order-parameter directions\n"
    "  modes      compute complete mode details\n"
    "  run        execute a JSON case through one pipeline stage\n"
    "  serve      serve the standalone frontend\n"
    "  show       show a JSON case in a running graphical interface\n";

static void locate(const char *argv0)
{
    char raw[PATH_MAX] = "";
#if defined(_WIN32)
    if (plat_exe_path(raw, sizeof raw)) raw[0] = 0;
#elif defined(__APPLE__)
    uint32_t size = sizeof raw;
    if (_NSGetExecutablePath(raw, &size)) raw[0] = 0;
#else
    ssize_t n = readlink("/proc/self/exe", raw, sizeof raw - 1);
    raw[n > 0 ? n : 0] = 0;
#endif
    if (!raw[0]) snprintf(raw, sizeof raw, "%s", argv0);
    if (!getcwd(cwd, sizeof cwd)) strcpy(cwd, "/");
#ifdef _WIN32
    for (char *c = cwd; *c; c++)
        if (*c == '\\') *c = '/';
#endif
    resolve(raw, cwd, exe_dir);
    parent_of(exe_dir);
}

int main(int argc, char **argv)
{
    locate(argv[0]);
    if (argc < 2 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        fputs(usage_text, argc < 2 ? stderr : stdout);
        return argc < 2 ? 2 : 0;
    }
    if (!strcmp(argv[1], "--version")) {
        puts("APOSTRUCT 0.1.0");
        return 0;
    }
    command = argv[1];
    static const char *const commands[] = {"info", "kpoints", "irreps", "opds", "modes", "run", "serve", "show"};
    int known = 0;
    for (int i = 0; i < 8 && !known; i++) known = !strcmp(command, commands[i]);
    if (!known) {
        fprintf(stderr, "apo: error: argument command: invalid choice: '%s' (choose from info, kpoints, irreps, opds, modes, run, serve, show)\n", command);
        return 2;
    }
    char blob[PATH_MAX * 2];
    const char *env = getenv("APO_SOURCE_BLOB");
    if (env) snprintf(blob, sizeof blob, "%s", env);
    else snprintf(blob, sizeof blob, "%s/source.blob", exe_dir);
    if (tables_open(blob)) {
        fprintf(stderr, "apo: error: the Source blob is missing or unreadable: %s\n", blob);
        return 2;
    }
    args_t a;
    parse_args(argc, argv, 2, &a);
    if (a.sg) write_sg_cif(&a);
    if (!strcmp(command, "serve")) return serve_command(&a);
    if (!strcmp(command, "show")) return show_command(&a);
    const char *stage = !strcmp(command, "run") ? a.upto : command;
    jbuf request = {0};
    jb_reset(&request);
    if (!strcmp(command, "info")) {
        if (!a.structure) usage_error("info requires structure");
        char path[PATH_MAX];
        resolve(a.structure, cwd, path);
        jb_lit(&request, "{\"structure\": ");
        jb_cstring(&request, path);
        jb_lit(&request, ", \"stage\": \"info\"");
    } else {
        apo_case c;
        cli_case(&a, stage, &c);
        const char *format = a.format ? a.format : "json";
        if (a.format && !strcmp(a.format, "text") && strcmp(stage, "modes")) usage_error("text output is available only for modes");
        jb_lit(&request, "{");
        put_case(&request, &c);
        jb_lit(&request, ", \"stage\": ");
        jb_cstring(&request, stage);
        jb_lit(&request, ", \"format\": ");
        jb_cstring(&request, format);
    }
    jb_lit(&request, a.full_state ? ", \"full_state\": true" : ", \"full_state\": false");
    char indent[48];
    snprintf(indent, sizeof indent, ", \"indent\": %d}", a.indent);
    jb_lit(&request, indent);
    exact_certify = a.certify;
    int status = run_request(&request, a.output);
    jb_free(&request);
    if (a.certify)
        fprintf(stderr,
                "{\"certify\": {\"checked\": %ld, \"failed\": %ld, \"outside_field\": %ld, \"undecided\": %ld, \"walks\": %ld, "
                "\"walk_nodes\": %ld, \"rank_failed\": %ld, \"intersections\": %ld, \"missing\": %ld, \"walk_exact\": %ld, "
                "\"walk_exact_failed\": %ld, \"walk_undecided\": %ld}}\n",
                exact_certified[0], exact_certified[1], exact_certified[2], exact_certified[3], exact_walk_certified[0],
                exact_walk_certified[1], exact_walk_certified[2], exact_walk_certified[3], exact_walk_certified[4], exact_walk_certified[5], exact_walk_certified[6],
                exact_walk_certified[7]);
    return status;
}
