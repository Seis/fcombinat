#include "solver/rank/cycle.h"

#include <limits.h>
#include <string.h>

#include "solver/math.h"

void rank_cycle_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res, int depth) {
  if (strcmp(obj->name, "Cycle") != 0) {
    fmpz_set_si(res, -1);
    return;
  }

  int n = get_obj_size(obj);
  if (n == 0) {
    fmpz_zero(res);
    return;
  }
  if (obj->num_children < 1) {
    fmpz_set_si(res, -1);
    return;
  }
  // Canonicalize: rotate children so min-label child is first
  if (obj->num_children > 1) {
    int min_pos = 0;
    int min_val = INT_MAX;
    for (int i = 0; i < obj->num_children; i++) {
      int *lbl = malloc(sizeof(int) * 128);
      int ls = 0, lc = 128;
      collect_labels(obj->children[i], &lbl, &ls, &lc);
      for (int j = 0; j < ls; j++) {
        if (lbl[j] < min_val) {
          min_val = lbl[j];
          min_pos = i;
        }
      }
      free(lbl);
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

  Object *Comp_obj = obj->children[0];
  Object *Rest_obj;
  if (obj->num_children == 1) {
    Rest_obj = new_object(OBJ_STRUCT, "Seq");
  } else {
    Rest_obj = new_object(OBJ_STRUCT, "Seq");
    Rest_obj->num_children = obj->num_children - 1;
    Rest_obj->children = malloc(sizeof(Object *) * Rest_obj->num_children);
    for (int i = 0; i < Rest_obj->num_children; i++)
      Rest_obj->children[i] = obj->children[i + 1];
  }

  Expr *child_expr = (Expr *)expr->component;
  Expr seq_expr_struct;
  seq_expr_struct.type = SEQUENCE;
  seq_expr_struct.component = child_expr;
  if (expr->restriction != NONE) {
    if (expr->restriction == GREATER && expr->limit == 0) {
      seq_expr_struct.restriction = NONE;
      seq_expr_struct.limit = 0;
    } else {
      seq_expr_struct.restriction = expr->restriction;
      seq_expr_struct.limit = expr->limit - 1;
    }
  } else {
    seq_expr_struct.restriction = NONE;
    seq_expr_struct.limit = 0;
  }
  Expr *seq_expr = &seq_expr_struct;

  int k = get_obj_size(Comp_obj);
  int n_total = n;

  int *pool = malloc(sizeof(int) * 128);
  int ps = 0;
  int pc = 128;
  collect_labels(obj, &pool, &ps, &pc);
  qsort(pool, ps, sizeof(int), cmp_int);
  int min_label = pool[0];

  int *labels_C = malloc(sizeof(int) * 128);
  int cs = 0;
  int cc = 128;
  collect_labels(Comp_obj, &labels_C, &cs, &cc);
  qsort(labels_C, cs, sizeof(int), cmp_int);

  int *pool_rem = &pool[1];
  int *labels_rem = malloc(sizeof(int) * k);
  int idx = 0;
  for (int i = 0; i < cs; i++)
    if (labels_C[i] != min_label)
      labels_rem[idx++] = labels_C[i];

  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  rank_combination(subset_rank, n_total - 1, k - 1, pool_rem, labels_rem);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  rank_e(ctx, child_expr, Comp_obj, rank_A, depth);

  fmpz_t rank_Rest;
  fmpz_init(rank_Rest);
  rank_e(ctx, seq_expr, Rest_obj, rank_Rest, depth);

  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_S;
  fmpz_init(count_S);
  get_expr_count(count_A, ctx, child_expr, k);
  get_expr_count(count_S, ctx, seq_expr, n_total - k);

  fmpz_t term1;
  fmpz_init(term1);
  fmpz_mul(term1, count_A, count_S);
  fmpz_mul(term1, term1, subset_rank);
  fmpz_t term2;
  fmpz_init(term2);
  fmpz_mul(term2, rank_A, count_S);

  fmpz_add(res, term1, term2);
  fmpz_add(res, res, rank_Rest);

  fmpz_t bin;
  fmpz_init(bin);
  fmpz_t weight;
  fmpz_init(weight);
  for (int i = 1; i < k; i++) {
    get_expr_count(count_A, ctx, child_expr, i);
    get_expr_count(count_S, ctx, seq_expr, n_total - i);
    fmpz_bin_uiui(bin, n_total - 1, i - 1);
    fmpz_mul(weight, count_A, count_S);
    fmpz_mul(weight, weight, bin);
    fmpz_add(res, res, weight);
  }

  if (Rest_obj->num_children > 0)
    free(Rest_obj->children);
  free(Rest_obj);
  free(labels_C);
  free(labels_rem);
  free(pool);
  fmpz_clear(subset_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_Rest);
  fmpz_clear(count_A);
  fmpz_clear(count_S);
  fmpz_clear(term1);
  fmpz_clear(term2);
  fmpz_clear(bin);
  fmpz_clear(weight);
}
