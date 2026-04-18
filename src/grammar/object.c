#include "grammar/object.h"

#include <ctype.h>
#include <stdio.h>
#include <string.h>

Object *new_object(ObjectType type, char *name) {
  Object *obj = malloc(sizeof(Object));
  obj->type = type;
  obj->name = name ? strdup(name) : NULL;
  obj->label = 0;
  obj->children = NULL;
  obj->num_children = 0;
  return obj;
}

void free_object(Object *obj) {
  if (!obj)
    return;
  if (obj->name)
    free(obj->name);
  if (obj->children) {
    for (int i = 0; i < obj->num_children; i++)
      free_object(obj->children[i]);
    free(obj->children);
  }
  free(obj);
}

void print_object(Object *obj) {
  if (!obj)
    return;
  if (obj->type == OBJ_ATOM) {
    if (obj->label == 0)
      printf("Z()");
    else
      printf("%s(%d)", obj->name, obj->label);
  } else if (obj->type == OBJ_EPSILON) {
    printf("Eps");
  } else {
    printf("%s(", obj->name);
    for (int i = 0; i < obj->num_children; i++) {
      print_object(obj->children[i]);
      if (i < obj->num_children - 1)
        printf(", ");
    }
    printf(")");
  }
}

static char *skip_ws(char *s) {
  while (isspace(*s))
    s++;
  return s;
}

static char *parse_name(char **ptr) {
  char *s = *ptr;
  s = skip_ws(s);
  if (!isalpha(*s))
    return NULL;
  char *start = s;
  while (isalnum(*s) || *s == '_')
    s++;
  int len = s - start;
  char *name = malloc(len + 1);
  strncpy(name, start, len);
  name[len] = '\0';
  *ptr = s;
  return name;
}

static Object *parse_obj_recursive(char **ptr) {
  char *s = *ptr;
  s = skip_ws(s);

  if (*s == '\0')
    return NULL;

  char *name = parse_name(&s);
  if (!name)
    return NULL;

  s = skip_ws(s);

  if (*s == '(') {
    s++;
    s = skip_ws(s);

    if (strcmp(name, "Z") == 0) {
      int label = 0;
      if (*s == ')') {
        /* Z() — unlabeled atom with label=0 */
      } else if (isdigit(*s)) {
        label = strtol(s, &s, 10);
        s = skip_ws(s);
        if (*s != ')') {
          free(name);
          return NULL;
        }
      } else {
        free(name);
        return NULL;
      }
      s++;

      Object *obj = new_object(OBJ_ATOM, name);
      obj->label = label;
      free(name);
      *ptr = s;
      return obj;
    }

    Object *obj = new_object(OBJ_STRUCT, name);
    free(name);

    Object **kids = malloc(sizeof(Object *) * 16);
    int cap = 16;
    int count = 0;

    while (*s != ')' && *s != '\0') {
      Object *child = parse_obj_recursive(&s);
      if (!child)
        break;

      if (count >= cap) {
        cap *= 2;
        kids = realloc(kids, sizeof(Object *) * cap);
      }
      kids[count++] = child;

      s = skip_ws(s);
      if (*s == ',') {
        s++;
        s = skip_ws(s);
      }
    }

    if (*s == ')')
      s++;

    obj->children = kids;
    obj->num_children = count;
    *ptr = s;
    return obj;
  } else {
    if (strcmp(name, "Eps") == 0) {
      Object *obj = new_object(OBJ_EPSILON, "Eps");
      free(name);
      *ptr = s;
      return obj;
    }
    free(name);
    return NULL;
  }
}

Object *parse_object(char *str) {
  char *s = str;
  return parse_obj_recursive(&s);
}

int get_obj_size(Object *obj) {
  if (obj->type == OBJ_ATOM)
    return 1;
  int s = 0;
  for (int i = 0; i < obj->num_children; i++)
    s += get_obj_size(obj->children[i]);
  return s;
}

void collect_labels(Object *obj, int **list, int *size, int *cap) {
  if (obj->type == OBJ_ATOM) {
    if (*size >= *cap) {
      *cap *= 2;
      *list = realloc(*list, sizeof(int) * (*cap));
    }
    (*list)[(*size)++] = obj->label;
  } else {
    for (int i = 0; i < obj->num_children; i++) {
      collect_labels(obj->children[i], list, size, cap);
    }
  }
}

// --- Canonicalization ---

static int obj_get_min_label(Object *obj) {
  if (obj->type == OBJ_ATOM)
    return obj->label;
  int min = 0x7fffffff;
  for (int i = 0; i < obj->num_children; i++) {
    int m = obj_get_min_label(obj->children[i]);
    if (m < min)
      min = m;
  }
  return min;
}

static int cmp_obj_by_min_label(const void *a, const void *b) {
  Object *oa = *(Object **)a;
  Object *ob = *(Object **)b;
  return obj_get_min_label(oa) - obj_get_min_label(ob);
}

void canonicalize_sets(Object *obj) {
  if (!obj)
    return;
  // Recurse into children first
  for (int i = 0; i < obj->num_children; i++)
    canonicalize_sets(obj->children[i]);

  if (obj->type != OBJ_STRUCT || !obj->name || obj->num_children < 2)
    return;

  // Sort Set children by min-label
  if (strcmp(obj->name, "Set") == 0) {
    qsort(obj->children, obj->num_children, sizeof(Object *),
          cmp_obj_by_min_label);
  }
  // Rotate Cycle children so min-label child is first
  else if (strcmp(obj->name, "Cycle") == 0) {
    int min_pos = 0;
    int min_val = obj_get_min_label(obj->children[0]);
    for (int i = 1; i < obj->num_children; i++) {
      int v = obj_get_min_label(obj->children[i]);
      if (v < min_val) {
        min_val = v;
        min_pos = i;
      }
    }
    if (min_pos != 0) {
      int nc = obj->num_children;
      Object **rotated = malloc(sizeof(Object *) * nc);
      for (int i = 0; i < nc; i++)
        rotated[i] = obj->children[(i + min_pos) % nc];
      memcpy(obj->children, rotated, sizeof(Object *) * nc);
      free(rotated);
    }
  }
}

// --- Serialization ---

static int obj_to_string_len(Object *obj) {
  if (!obj)
    return 0;
  if (obj->type == OBJ_ATOM) {
    // "Z(label)" - label can be up to ~10 digits
    return 16;
  }
  if (obj->type == OBJ_EPSILON) {
    return 5; // "Eps()"
  }
  // "Name(" + children with ", " separators + ")"
  int len = (obj->name ? strlen(obj->name) : 0) + 2;
  for (int i = 0; i < obj->num_children; i++) {
    len += obj_to_string_len(obj->children[i]);
    if (i > 0)
      len += 2; // ", "
  }
  return len;
}

static char *obj_to_string_write(Object *obj, char *out) {
  if (!obj)
    return out;
  if (obj->type == OBJ_ATOM) {
    if (obj->label == 0)
      out += sprintf(out, "Z()");
    else
      out += sprintf(out, "Z(%d)", obj->label);
    return out;
  }
  if (obj->type == OBJ_EPSILON) {
    out += sprintf(out, "Eps()");
    return out;
  }
  out += sprintf(out, "%s(", obj->name ? obj->name : "?");
  for (int i = 0; i < obj->num_children; i++) {
    if (i > 0) {
      *out++ = ',';
      *out++ = ' ';
    }
    out = obj_to_string_write(obj->children[i], out);
  }
  *out++ = ')';
  *out = '\0';
  return out;
}

char *object_to_string(Object *obj) {
  int len = obj_to_string_len(obj) + 1;
  char *buf = malloc(len);
  obj_to_string_write(obj, buf);
  return buf;
}
