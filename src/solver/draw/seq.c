#include "solver/draw/seq.h"

#include "solver/math.h"

char *draw_seq_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                    int *labels) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return xstrdup("Seq()");

  fmpz_t weight;
  fmpz_init(weight);
  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_Seq;
  fmpz_init(count_Seq);
  fmpz_t bin;
  fmpz_init(bin);

  // Compute total weight
  fmpz_t total;
  fmpz_init(total);
  for (int k = 1; k <= n; k++) {
    Expr k_tail_expr_struct = *expr;
    if (expr->restriction != NONE) {
      if (expr->restriction == GREATER && expr->limit == 0) {
        k_tail_expr_struct.restriction = NONE;
        k_tail_expr_struct.limit = 0;
      } else {
        k_tail_expr_struct.limit = expr->limit - 1;
      }
    }
    Expr *k_tail_expr = &k_tail_expr_struct;

    get_expr_count(count_A, ctx, child, k);
    get_expr_count(count_Seq, ctx, k_tail_expr, n - k);
    binom_ui(bin, n, k);

    fmpz_mul(weight, count_A, count_Seq);
    fmpz_mul(weight, weight, bin);
    fmpz_add(total, total, weight);
  }

  if (fmpz_is_zero(total)) {
    fmpz_clear(weight);
    fmpz_clear(count_A);
    fmpz_clear(count_Seq);
    fmpz_clear(bin);
    fmpz_clear(total);
    return xstrdup("ErrorSeqEmpty");
  }

  // Draw uniform random in [0, total)
  fmpz_t current_rank;
  fmpz_init(current_rank);
  fmpz_randm(current_rank, state, total);

  int chosen_k = -1;

  for (int k = 1; k <= n; k++) {
    Expr k_tail_expr_struct = *expr;
    if (expr->restriction != NONE) {
      if (expr->restriction == GREATER && expr->limit == 0) {
        k_tail_expr_struct.restriction = NONE;
        k_tail_expr_struct.limit = 0;
      } else {
        k_tail_expr_struct.limit = expr->limit - 1;
      }
    }
    Expr *k_tail_expr = &k_tail_expr_struct;

    get_expr_count(count_A, ctx, child, k);
    get_expr_count(count_Seq, ctx, k_tail_expr, n - k);
    binom_ui(bin, n, k);

    fmpz_mul(weight, count_A, count_Seq);
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
    fmpz_clear(count_Seq);
    fmpz_clear(bin);
    fmpz_clear(total);
    fmpz_clear(current_rank);
    return xstrdup("ErrorSeqBounds");
  }

  // Draw independent random sub-ranks
  fmpz_t rank_A;
  fmpz_init(rank_A);
  fmpz_t rank_Seq;
  fmpz_init(rank_Seq);

  if (!fmpz_is_zero(count_A))
    fmpz_randm(rank_A, state, count_A);
  if (!fmpz_is_zero(count_Seq))
    fmpz_randm(rank_Seq, state, count_Seq);

  // Draw random label subset
  fmpz_t combo_total;
  fmpz_init(combo_total);
  binom_ui(combo_total, n, chosen_k);
  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  if (!fmpz_is_zero(combo_total))
    fmpz_randm(subset_rank, state, combo_total);

  int *component_labels = xmalloc(sizeof(int) * chosen_k);
  int *rest_labels = xmalloc(sizeof(int) * (n - chosen_k));
  unrank_combination(subset_rank, n, chosen_k, labels, component_labels,
                     rest_labels);

  char *res_A = draw_e(ctx, child, chosen_k, state, component_labels);

  Expr final_tail_expr_struct = *expr;
  if (expr->restriction != NONE && expr->limit > 1) {
    final_tail_expr_struct.limit = expr->limit - 1;
  } else if (expr->restriction != NONE && expr->limit <= 1) {
    final_tail_expr_struct.restriction = NONE;
    final_tail_expr_struct.limit = 0;
  }
  Expr *final_tail_expr = &final_tail_expr_struct;

  char *res_S = NULL;
  if (n - chosen_k == 0) {
    res_S = strdup("Seq()");
  } else {
    res_S = draw_e(ctx, final_tail_expr, n - chosen_k, state, rest_labels);
  }

  char *res = xmalloc(strlen(res_A) + strlen(res_S) + 20);
  if (strncmp(res_S, "Seq(", 4) == 0) {
    int lenS = strlen(res_S);
    if (lenS > 5) {
      res_S[lenS - 1] = '\0';
      sprintf(res, "Seq(%s, %s)", res_A, res_S + 4);
    } else {
      sprintf(res, "Seq(%s)", res_A);
    }
  } else {
    sprintf(res, "Seq(%s, %s)", res_A, res_S);
  }

  free(res_A);
  free(res_S);
  free(component_labels);
  free(rest_labels);
  fmpz_clear(weight);
  fmpz_clear(count_A);
  fmpz_clear(count_Seq);
  fmpz_clear(bin);
  fmpz_clear(total);
  fmpz_clear(current_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_Seq);
  fmpz_clear(combo_total);
  fmpz_clear(subset_rank);
  return res;
}
