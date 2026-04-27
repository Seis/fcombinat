#include "solver/unrank/prod.h"

#include <string.h>
#include "solver/math.h"

char *unrank_prod_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                       int *labels, int depth) {
  ExprList *el = (ExprList *)expr->component;
  if (el->size == 0)
    return unrank_e(ctx, NULL, n, rank, labels, depth);
  if (el->size == 1)
    return unrank_e(ctx, el->components[0], n, rank, labels, depth);

  Expr *A_expr = el->components[0];

  Expr *B_expr;
  ExprList *b_list = malloc(sizeof(ExprList));
  b_list->size = el->size - 1;
  b_list->components = &el->components[1];
  Expr b_expr_struct;
  b_expr_struct.type = PROD;
  b_expr_struct.component = b_list;
  B_expr = &b_expr_struct;

  fmpz_t weight;
  fmpz_init(weight);
  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_B;
  fmpz_init(count_B);
  fmpz_t bin;
  fmpz_init(bin);
  fmpz_t current_rank;
  fmpz_init_set(current_rank, rank);

  int chosen_k = -1;

  for (int k = 0; k <= n; k++) {
    get_expr_count(count_A, ctx, A_expr, k);
    get_expr_count(count_B, ctx, B_expr, n - k);
    binom_ui(bin, n, k);

    fmpz_mul(weight, count_A, count_B);
    fmpz_mul(weight, weight, bin);

    if (fmpz_cmp(current_rank, weight) < 0) {
      chosen_k = k;
      break;
    }
    fmpz_sub(current_rank, current_rank, weight);
  }

  if (chosen_k == -1) {
    fmpz_clear(weight);
    fmpz_clear(count_A);
    fmpz_clear(count_B);
    fmpz_clear(bin);
    fmpz_clear(current_rank);
    free(b_list);
    return strdup("ErrorProdBounds");
  }

  fmpz_t total_structs;
  fmpz_init(total_structs);
  fmpz_mul(total_structs, count_A, count_B);

  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  fmpz_t structure_rank;
  fmpz_init(structure_rank);

  fmpz_fdiv_qr(subset_rank, structure_rank, current_rank, total_structs);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  fmpz_t rank_B;
  fmpz_init(rank_B);

  fmpz_fdiv_qr(rank_A, rank_B, structure_rank, count_B);

  int *selected_labels = malloc(sizeof(int) * chosen_k);
  int *other_labels = malloc(sizeof(int) * (n - chosen_k));
  unrank_combination(subset_rank, n, chosen_k, labels, selected_labels,
                     other_labels);

  char *res_A = unrank_e(ctx, A_expr, chosen_k, rank_A, selected_labels, depth);
  char *res_B = unrank_e(ctx, B_expr, n - chosen_k, rank_B, other_labels, depth);

  int len = strlen(res_A) + strlen(res_B) + 10;
  char *res = malloc(len);
  sprintf(res, "Prod(%s, %s)", res_A, res_B);

  free(res_A);
  free(res_B);
  free(selected_labels);
  free(other_labels);
  free(b_list);

  fmpz_clear(weight);
  fmpz_clear(count_A);
  fmpz_clear(count_B);
  fmpz_clear(bin);
  fmpz_clear(current_rank);
  fmpz_clear(total_structs);
  fmpz_clear(subset_rank);
  fmpz_clear(structure_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_B);

  return res;
}

char *unrank_prod_unlabeled(Context *ctx, Expr *expr, int n, fmpz_t rank,
                             int depth) {
  ExprList *el = (ExprList *)expr->component;
  if (el->size == 0)
    return unrank_e(ctx, NULL, n, rank, NULL, depth);
  if (el->size == 1)
    return unrank_e(ctx, el->components[0], n, rank, NULL, depth);

  Expr *A_expr = el->components[0];

  ExprList *b_list = malloc(sizeof(ExprList));
  b_list->size = el->size - 1;
  b_list->components = &el->components[1];
  Expr b_expr_struct;
  b_expr_struct.type = PROD;
  b_expr_struct.component = b_list;
  b_expr_struct.restriction = expr->restriction;
  b_expr_struct.limit = expr->limit;
  Expr *B_expr = &b_expr_struct;

  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_B;
  fmpz_init(count_B);
  fmpz_t block;
  fmpz_init(block);
  fmpz_t current_rank;
  fmpz_init_set(current_rank, rank);

  int chosen_k = -1;

  for (int k = 0; k <= n; k++) {
    get_expr_count(count_A, ctx, A_expr, k);
    get_expr_count(count_B, ctx, B_expr, n - k);
    fmpz_mul(block, count_A, count_B);

    if (fmpz_cmp(current_rank, block) < 0) {
      chosen_k = k;
      break;
    }
    fmpz_sub(current_rank, current_rank, block);
  }

  if (chosen_k == -1) {
    fmpz_clear(count_A);
    fmpz_clear(count_B);
    fmpz_clear(block);
    fmpz_clear(current_rank);
    free(b_list);
    return strdup("ErrorProdUnlabeledBounds");
  }

  get_expr_count(count_A, ctx, A_expr, chosen_k);
  get_expr_count(count_B, ctx, B_expr, n - chosen_k);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  fmpz_t rank_B;
  fmpz_init(rank_B);

  /* current_rank encodes: rank_A * count_B + rank_B */
  fmpz_fdiv_qr(rank_A, rank_B, current_rank, count_B);

  char *res_A = unrank_e(ctx, A_expr, chosen_k, rank_A, NULL, depth);
  char *res_B = unrank_e(ctx, B_expr, n - chosen_k, rank_B, NULL, depth);

  int len = strlen(res_A) + strlen(res_B) + 10;
  char *res = malloc(len);
  sprintf(res, "Prod(%s, %s)", res_A, res_B);

  free(res_A);
  free(res_B);
  free(b_list);
  fmpz_clear(count_A);
  fmpz_clear(count_B);
  fmpz_clear(block);
  fmpz_clear(current_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_B);

  return res;
}
