/* A small JSON document model (io/json.c): parsed values live in one arena per
 * document; the writer appends to a growable buffer and prints floats the
 * way Python's repr (and so json.dumps) does. */

#ifndef APOSTRUCT_JSON_H
#define APOSTRUCT_JSON_H

#include <stddef.h>
#include <stdint.h>

typedef enum { J_NULL, J_FALSE, J_TRUE, J_INT, J_FLOAT, J_STRING, J_ARRAY, J_OBJECT } jkind;

typedef struct jvalue jvalue;

struct jvalue {
    jkind kind;
    int64_t count; /* string bytes, array items or object members */
    union {
        int64_t i;
        double f;
        const char *s;
        jvalue *items;   /* array */
        jvalue *members; /* object: key (J_STRING), value pairs */
    } u;
};

typedef struct jdoc jdoc;

/* NULL on a syntax error or memory failure. */
jdoc *jdoc_parse(const char *text, size_t length);
const jvalue *jdoc_root(const jdoc *doc);
void jdoc_free(jdoc *doc);

/* Object member by key, NULL when absent; array item, NULL out of range. */
const jvalue *jget(const jvalue *object, const char *key);
const jvalue *jat(const jvalue *array, int64_t index);
int jis(const jvalue *value, jkind kind);
int jtruthy(const jvalue *value);
int jstr_is(const jvalue *value, const char *text);

typedef struct {
    char *data;
    size_t used, cap;
    int failed;
} jbuf;

void jb_reset(jbuf *b);
void jb_free(jbuf *b);
void jb_raw(jbuf *b, const char *text, size_t length);
void jb_lit(jbuf *b, const char *text);
void jb_int(jbuf *b, int64_t value);
void jb_float(jbuf *b, double value);
void jb_string(jbuf *b, const char *text, size_t length);
void jb_cstring(jbuf *b, const char *text);
/* Re-serialize a parsed value with json.dumps' default separators. */
void jb_value(jbuf *b, const jvalue *value);
/* Python repr(float) into out (at least 32 bytes); returns its length. */
int py_float_repr(double value, char *out);

#endif
