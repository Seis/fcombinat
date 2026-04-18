#ifndef SOLVER_OBJECT_H
#define SOLVER_OBJECT_H

#include <stdlib.h>

// --- Object AST ---
typedef enum {
  OBJ_ATOM,    // Z
  OBJ_EPSILON, // Eps
  OBJ_STRUCT   // Set, Cycle, Seq, Prod, etc.
} ObjectType;

typedef struct Object_s {
  ObjectType type;
  char *name; // "Z", "Set", "Cycle", etc.
  int label;  // For Z(k), label=k. Otherwise 0.

  struct Object_s **children;
  int num_children;
} Object;

Object *new_object(ObjectType type, char *name);
void free_object(Object *obj);
void print_object(Object *obj);

Object *parse_object(char *str);

// Utilities that operate on Objects
int get_obj_size(Object *obj);
void collect_labels(Object *obj, int **list, int *size, int *cap);
int cmp_int(const void *a, const void *b);

// Canonicalization: sort Set children by min-label (recursive)
void canonicalize_sets(Object *obj);

// Serialize object tree to string (caller must free)
char *object_to_string(Object *obj);

#endif
