#include "solver/rank/union.h"


void rank_union_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res, int depth) {
  ExprList *el = (ExprList *)expr->component;
  int n = get_obj_size(obj);

  int branch_tag = -1;
  Object *inner_obj = obj;

  if (strncmp(obj->name, "Union_", 6) == 0) {
    branch_tag = atoi(obj->name + 6);
    if (obj->num_children == 1) {
      inner_obj = obj->children[0];
    }
  }

  fmpz_t accum;
  fmpz_init(accum);
  fmpz_zero(accum);
  fmpz_t count;
  fmpz_init(count);
  fmpz_t sub_rank;
  fmpz_init(sub_rank);

  int found = 0;
  for (int i = 0; i < el->size; i++) {
    if (branch_tag != -1 && i != branch_tag) {
      get_expr_count(count, ctx, el->components[i], n);
      fmpz_add(accum, accum, count);
      continue;
    }

    int match = 0;
    Expr *child = el->components[i];
    switch (child->type) {
    case SET:
      match = (strcmp(inner_obj->name, "Set") == 0);
      break;
    case POWERSET:
      match = (strcmp(inner_obj->name, "Set") == 0) ||
              (strcmp(inner_obj->name, "PowerSet") == 0);
      match = (strcmp(inner_obj->name, "PowerSet") == 0);
      break;
    case CYCLE:
      match = (strcmp(inner_obj->name, "Cycle") == 0);
      break;
    case SEQUENCE:
      match = (strcmp(inner_obj->name, "Seq") == 0);
      break;
    case PROD:
      match = (strcmp(inner_obj->name, "Prod") == 0) ||
              (n > 1 && strcmp(inner_obj->name, "Seq") != 0);
      break;
    case ATOM:
    case Z:
      match = (inner_obj->type == OBJ_ATOM);
      break;
    case EPSILON:
      match = (inner_obj->type == OBJ_EPSILON);
      break;
    case ID:
      match = 1;
      break;
    default:
      match = 1;
    }

    if (match || branch_tag == i) {
      rank_e(ctx, child, inner_obj, sub_rank, depth);
      if (fmpz_sgn(sub_rank) >= 0) {
        fmpz_add(res, accum, sub_rank);
        found = 1;
        break;
      }
    }

    get_expr_count(count, ctx, child, n);
    fmpz_add(accum, accum, count);
  }

  fmpz_clear(accum);
  fmpz_clear(count);
  fmpz_clear(sub_rank);
  if (!found)
    fmpz_set_si(res, -1);
}
