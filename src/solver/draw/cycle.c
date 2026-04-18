#include "solver/draw/cycle.h"

#include "solver/math.h"

char *draw_cycle_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                      int *labels) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return xstrdup("Cycle()");

  int min_label = labels[0];
  int *remaining_labels = &labels[1];
  int n_rem = n - 1;

  fmpz_t weight;
  fmpz_init(weight);
  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_Seq;
  fmpz_init(count_Seq);
  fmpz_t bin;
  fmpz_init(bin);

  Expr seq_expr_struct;
  seq_expr_struct.type = SEQUENCE;
  seq_expr_struct.component = child;
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

  // Compute total weight
  fmpz_t total;
  fmpz_init(total);
  for (int k = 1; k <= n; k++) {
    get_expr_count(count_A, ctx, child, k);
    get_expr_count(count_Seq, ctx, seq_expr, n - k);
    binom_ui(bin, n_rem, k - 1);

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
    return xstrdup("ErrorCycleEmpty");
  }

  // Draw uniform random in [0, total)
  fmpz_t current_rank;
  fmpz_init(current_rank);
  fmpz_randm(current_rank, state, total);

  int chosen_k = -1;

  for (int k = 1; k <= n; k++) {
    get_expr_count(count_A, ctx, child, k);
    get_expr_count(count_Seq, ctx, seq_expr, n - k);
    binom_ui(bin, n_rem, k - 1);

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
    return xstrdup("ErrorCycleBounds");
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
  binom_ui(combo_total, n_rem, chosen_k - 1);
  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  if (!fmpz_is_zero(combo_total))
    fmpz_randm(subset_rank, state, combo_total);

  int *component_labels = xmalloc(sizeof(int) * chosen_k);
  int *rest_labels = xmalloc(sizeof(int) * (n - chosen_k));
  component_labels[0] = min_label;
  unrank_combination(subset_rank, n_rem, chosen_k - 1, remaining_labels,
                     &component_labels[1], rest_labels);

  char *res_A = draw_e(ctx, child, chosen_k, state, component_labels);

  char *res_Seq = NULL;
  if (n - chosen_k == 0) {
    res_Seq = strdup("Seq()");
  } else {
    res_Seq = draw_e(ctx, seq_expr, n - chosen_k, state, rest_labels);
  }

  char *res = xmalloc(strlen(res_A) + strlen(res_Seq) + 20);
  if (strncmp(res_Seq, "Seq(", 4) == 0) {
    int lenS = strlen(res_Seq);
    if (lenS > 5) {
      res_Seq[lenS - 1] = '\0';
      sprintf(res, "Cycle(%s, %s)", res_A, res_Seq + 4);
    } else {
      sprintf(res, "Cycle(%s)", res_A);
    }
  } else {
    sprintf(res, "Cycle(%s, %s)", res_A, res_Seq);
  }

  free(res_A);
  free(res_Seq);
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
