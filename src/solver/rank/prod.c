#include "solver/rank/prod.h"


#include "solver/math.h"

void rank_prod_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res, int depth) {
  ExprList *el = (ExprList *)expr->component;
  if (el->size == 0) {
    fmpz_zero(res);
    return;
  }
  if (el->size == 1) {
    rank_e(ctx, el->components[0], obj, res, depth);
    return;
  }

  if (strcmp(obj->name, "Prod") != 0) {
    fmpz_set_si(res, -1);
    return;
  }
  if (obj->num_children < 2) {
    fmpz_set_si(res, -1);
    return;
  }

  Expr *A_expr = el->components[0];
  Object *A_obj = obj->children[0];
  Object *B_obj = obj->children[1];

  Expr *B_expr;
  ExprList *b_list = malloc(sizeof(ExprList));
  b_list->size = el->size - 1;
  b_list->components = &el->components[1];
  Expr b_expr_struct;
  b_expr_struct.type = PROD;
  b_expr_struct.component = b_list;
  B_expr = &b_expr_struct;

  int k = get_obj_size(A_obj);
  int n = get_obj_size(obj);

  int *labels_A = malloc(sizeof(int) * 128);
  int size_A = 0;
  int cap_A = 128;
  collect_labels(A_obj, &labels_A, &size_A, &cap_A);
  qsort(labels_A, size_A, sizeof(int), cmp_int);

  int *pool = malloc(sizeof(int) * 128);
  int pool_size = 0;
  int pool_cap = 128;
  collect_labels(obj, &pool, &pool_size, &pool_cap);
  qsort(pool, pool_size, sizeof(int), cmp_int);

  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  rank_combination(subset_rank, pool_size, k, pool, labels_A);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  rank_e(ctx, A_expr, A_obj, rank_A, depth);

  fmpz_t rank_B;
  fmpz_init(rank_B);
  rank_e(ctx, B_expr, B_obj, rank_B, depth);

  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_B;
  fmpz_init(count_B);
  get_expr_count(count_A, ctx, A_expr, k);
  get_expr_count(count_B, ctx, B_expr, n - k);

  fmpz_t term1;
  fmpz_init(term1);
  fmpz_mul(term1, count_A, count_B);
  fmpz_mul(term1, term1, subset_rank);

  fmpz_t term2;
  fmpz_init(term2);
  fmpz_mul(term2, rank_A, count_B);

  fmpz_add(res, term1, term2);
  fmpz_add(res, res, rank_B);

  fmpz_t bin;
  fmpz_init(bin);
  fmpz_t weight;
  fmpz_init(weight);

  for (int i = 0; i < k; i++) {
    get_expr_count(count_A, ctx, A_expr, i);
    get_expr_count(count_B, ctx, B_expr, n - i);
    fmpz_bin_uiui(bin, n, i);
    fmpz_mul(weight, count_A, count_B);
    fmpz_mul(weight, weight, bin);
    fmpz_add(res, res, weight);
  }

  free(b_list);
  free(labels_A);
  free(pool);
  fmpz_clear(subset_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_B);
  fmpz_clear(count_A);
  fmpz_clear(count_B);
  fmpz_clear(term1);
  fmpz_clear(term2);
  fmpz_clear(bin);
  fmpz_clear(weight);
}

void rank_prod_unlabeled(Context *ctx, Expr *expr, Object *obj, fmpz_t res, int depth) {
  ExprList *el = (ExprList *)expr->component;
  if (el->size == 0) {
    fmpz_zero(res);
    return;
  }
  if (el->size == 1) {
    rank_e(ctx, el->components[0], obj, res, depth);
    return;
  }

  if (strcmp(obj->name, "Prod") != 0) {
    fmpz_set_si(res, -1);
    return;
  }
  if (obj->num_children < 2) {
    fmpz_set_si(res, -1);
    return;
  }

  Expr *A_expr = el->components[0];
  Object *A_obj = obj->children[0];
  Object *B_obj = obj->children[1];

  ExprList *b_list = malloc(sizeof(ExprList));
  b_list->size = el->size - 1;
  b_list->components = &el->components[1];
  Expr b_expr_struct;
  b_expr_struct.type = PROD;
  b_expr_struct.component = b_list;
  b_expr_struct.restriction = expr->restriction;
  b_expr_struct.limit = expr->limit;
  Expr *B_expr = &b_expr_struct;

  int k = get_obj_size(A_obj);
  int n = get_obj_size(obj);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  rank_e(ctx, A_expr, A_obj, rank_A, depth);

  fmpz_t rank_B;
  fmpz_init(rank_B);
  rank_e(ctx, B_expr, B_obj, rank_B, depth);

  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_B;
  fmpz_init(count_B);
  get_expr_count(count_A, ctx, A_expr, k);
  get_expr_count(count_B, ctx, B_expr, n - k);

  /* accumulated = sum_{i=0}^{k-1} f_A(i) * f_B(n-i) */
  fmpz_t accumulated;
  fmpz_init(accumulated);
  fmpz_zero(accumulated);
  fmpz_t ca;
  fmpz_init(ca);
  fmpz_t cb;
  fmpz_init(cb);
  fmpz_t tmp;
  fmpz_init(tmp);

  for (int i = 0; i < k; i++) {
    get_expr_count(ca, ctx, A_expr, i);
    get_expr_count(cb, ctx, B_expr, n - i);
    fmpz_mul(tmp, ca, cb);
    fmpz_add(accumulated, accumulated, tmp);
  }

  /* res = accumulated + rank_A * count_B + rank_B */
  fmpz_mul(tmp, rank_A, count_B);
  fmpz_add(res, accumulated, tmp);
  fmpz_add(res, res, rank_B);

  free(b_list);
  fmpz_clear(rank_A);
  fmpz_clear(rank_B);
  fmpz_clear(count_A);
  fmpz_clear(count_B);
  fmpz_clear(accumulated);
  fmpz_clear(ca);
  fmpz_clear(cb);
  fmpz_clear(tmp);
}
