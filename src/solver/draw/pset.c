#include "solver/draw/pset.h"

#include "solver/math.h"

char *draw_powerset_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                         int *labels) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return xstrdup("PowerSet()");

  int min_label = labels[0];
  int *remaining_labels = &labels[1];
  int n_rem = n - 1;

  Expr tail_expr_struct = *expr;
  tail_expr_struct.restriction = NONE;
  tail_expr_struct.limit = 0;
  Expr *tail_expr = &tail_expr_struct;

  fmpz_t weight;
  fmpz_init(weight);
  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_P;
  fmpz_init(count_P);
  fmpz_t bin;
  fmpz_init(bin);

  // Compute total weight
  fmpz_t total;
  fmpz_init(total);
  for (int k = 1; k <= n; k++) {
    get_expr_count(count_A, ctx, child, k);
    get_expr_count(count_P, ctx, tail_expr, n - k);
    binom_ui(bin, n_rem, k - 1);

    fmpz_mul(weight, count_A, count_P);
    fmpz_mul(weight, weight, bin);
    fmpz_add(total, total, weight);
  }

  if (fmpz_is_zero(total)) {
    fmpz_clear(weight);
    fmpz_clear(count_A);
    fmpz_clear(count_P);
    fmpz_clear(bin);
    fmpz_clear(total);
    return xstrdup("ErrorPowerSetEmpty");
  }

  // Draw uniform random in [0, total)
  fmpz_t current_rank;
  fmpz_init(current_rank);
  fmpz_randm(current_rank, state, total);

  int chosen_k = -1;

  for (int k = 1; k <= n; k++) {
    get_expr_count(count_A, ctx, child, k);
    get_expr_count(count_P, ctx, tail_expr, n - k);
    binom_ui(bin, n_rem, k - 1);

    fmpz_mul(weight, count_A, count_P);
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
    fmpz_clear(count_P);
    fmpz_clear(bin);
    fmpz_clear(total);
    fmpz_clear(current_rank);
    return xstrdup("ErrorPowerSetBounds");
  }

  // Draw independent random sub-ranks
  fmpz_t rank_A;
  fmpz_init(rank_A);
  fmpz_t rank_P;
  fmpz_init(rank_P);

  if (!fmpz_is_zero(count_A))
    fmpz_randm(rank_A, state, count_A);
  if (!fmpz_is_zero(count_P))
    fmpz_randm(rank_P, state, count_P);

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

  char *res_P = NULL;
  if (n - chosen_k == 0) {
    res_P = strdup("");
  } else {
    res_P = draw_e(ctx, tail_expr, n - chosen_k, state, rest_labels);
  }

  char *res = xmalloc(strlen(res_A) + (res_P ? strlen(res_P) : 0) + 20);
  if (n - chosen_k == 0) {
    sprintf(res, "PowerSet(%s)", res_A);
  } else {
    if (strncmp(res_P, "PowerSet(", 9) == 0) {
      res_P[strlen(res_P) - 1] = '\0';
      sprintf(res, "PowerSet(%s, %s)", res_A, res_P + 9);
    } else {
      sprintf(res, "PowerSet(%s, %s)", res_A, res_P);
    }
  }

  free(res_A);
  if (res_P)
    free(res_P);
  free(component_labels);
  free(rest_labels);
  fmpz_clear(weight);
  fmpz_clear(count_A);
  fmpz_clear(count_P);
  fmpz_clear(bin);
  fmpz_clear(total);
  fmpz_clear(current_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_P);
  fmpz_clear(combo_total);
  fmpz_clear(subset_rank);

  return res;
}
