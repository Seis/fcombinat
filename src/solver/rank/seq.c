#include "solver/rank/seq.h"


#include "solver/math.h"

void rank_seq_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res) {
  if (strcmp(obj->name, "Seq") != 0) {
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
  int k = get_obj_size(Comp_obj);
  int n_total = n;

  int *pool = malloc(sizeof(int) * 128);
  int pool_size = 0;
  int caps = 128;
  collect_labels(obj, &pool, &pool_size, &caps);
  qsort(pool, pool_size, sizeof(int), cmp_int);

  int *labels_A = malloc(sizeof(int) * 128);
  int size_A = 0;
  int capA = 128;
  collect_labels(Comp_obj, &labels_A, &size_A, &capA);
  qsort(labels_A, size_A, sizeof(int), cmp_int);

  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  rank_combination(subset_rank, n_total, k, pool, labels_A);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  rank_e(ctx, child_expr, Comp_obj, rank_A);

  Expr tail_expr_struct = *expr;
  if (expr->restriction != NONE && expr->limit > 1) {
    tail_expr_struct.limit = expr->limit - 1;
  } else if (expr->restriction != NONE && expr->limit <= 1) {
    tail_expr_struct.restriction = NONE;
    tail_expr_struct.limit = 0;
  }
  Expr *tail_expr = &tail_expr_struct;

  fmpz_t rank_Rest;
  fmpz_init(rank_Rest);
  rank_e(ctx, tail_expr, Rest_obj, rank_Rest);

  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_S;
  fmpz_init(count_S);
  get_expr_count(count_A, ctx, child_expr, k);
  get_expr_count(count_S, ctx, tail_expr, n_total - k);

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
    Expr offset_tail_expr_struct = *expr;
    if (expr->restriction != NONE) {
      if (expr->restriction == GREATER && expr->limit == 0) {
        offset_tail_expr_struct.restriction = NONE;
        offset_tail_expr_struct.limit = 0;
      } else {
        offset_tail_expr_struct.limit = expr->limit - 1;
      }
    }
    Expr *offset_tail_expr = &offset_tail_expr_struct;

    get_expr_count(count_A, ctx, child_expr, i);
    get_expr_count(count_S, ctx, offset_tail_expr, n_total - i);
    fmpz_bin_uiui(bin, n_total, i);
    fmpz_mul(weight, count_A, count_S);
    fmpz_mul(weight, weight, bin);
    fmpz_add(res, res, weight);
  }

  if (Rest_obj->num_children > 0)
    free(Rest_obj->children);
  free(Rest_obj);
  free(labels_A);
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

void rank_seq_unlabeled(Context *ctx, Expr *expr, Object *obj, fmpz_t res) {
  if (strcmp(obj->name, "Seq") != 0) {
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

  Expr *child_expr = (Expr *)expr->component;
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

  Expr tail_expr_struct = *expr;
  if (expr->restriction != NONE && expr->limit > 1) {
    tail_expr_struct.limit = expr->limit - 1;
  } else if (expr->restriction != NONE && expr->limit <= 1) {
    tail_expr_struct.restriction = NONE;
    tail_expr_struct.limit = 0;
  }
  Expr *tail_expr = &tail_expr_struct;

  int k = get_obj_size(Comp_obj);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  rank_e(ctx, child_expr, Comp_obj, rank_A);

  fmpz_t rank_Rest;
  fmpz_init(rank_Rest);
  rank_e(ctx, tail_expr, Rest_obj, rank_Rest);

  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_S;
  fmpz_init(count_S);
  get_expr_count(count_A, ctx, child_expr, k);
  get_expr_count(count_S, ctx, tail_expr, n - k);

  /* accumulated = sum_{i=1}^{k-1} f_A(i) * f_Seq(n-i) */
  fmpz_t accumulated;
  fmpz_init(accumulated);
  fmpz_zero(accumulated);
  fmpz_t ca;
  fmpz_init(ca);
  fmpz_t cs;
  fmpz_init(cs);
  fmpz_t tmp;
  fmpz_init(tmp);

  for (int i = 1; i < k; i++) {
    Expr offset_tail_expr_struct = *expr;
    if (expr->restriction != NONE) {
      if (expr->restriction == GREATER && expr->limit == 0) {
        offset_tail_expr_struct.restriction = NONE;
        offset_tail_expr_struct.limit = 0;
      } else {
        offset_tail_expr_struct.limit = expr->limit - 1;
      }
    }
    Expr *offset_tail_expr = &offset_tail_expr_struct;

    get_expr_count(ca, ctx, child_expr, i);
    get_expr_count(cs, ctx, offset_tail_expr, n - i);
    fmpz_mul(tmp, ca, cs);
    fmpz_add(accumulated, accumulated, tmp);
  }

  /* res = accumulated + rank_A * count_S + rank_Rest */
  fmpz_mul(tmp, rank_A, count_S);
  fmpz_add(res, accumulated, tmp);
  fmpz_add(res, res, rank_Rest);

  if (Rest_obj->num_children > 0)
    free(Rest_obj->children);
  free(Rest_obj);
  fmpz_clear(rank_A);
  fmpz_clear(rank_Rest);
  fmpz_clear(count_A);
  fmpz_clear(count_S);
  fmpz_clear(accumulated);
  fmpz_clear(ca);
  fmpz_clear(cs);
  fmpz_clear(tmp);
}
