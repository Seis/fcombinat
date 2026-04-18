#include "solver/unrank/cycle.h"


#include "solver/math.h"

char *unrank_cycle_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                        int *labels) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return strdup("Cycle()");

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
  fmpz_t current_rank;
  fmpz_init_set(current_rank, rank);

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
    fmpz_clear(current_rank);
    return strdup("ErrorCycleBounds");
  }

  fmpz_t total_structs;
  fmpz_init(total_structs);
  fmpz_mul(total_structs, count_A, count_Seq);
  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  fmpz_t structure_rank;
  fmpz_init(structure_rank);
  fmpz_fdiv_qr(subset_rank, structure_rank, current_rank, total_structs);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  fmpz_t rank_Seq;
  fmpz_init(rank_Seq);
  fmpz_fdiv_qr(rank_A, rank_Seq, structure_rank, count_Seq);

  int *component_labels = malloc(sizeof(int) * chosen_k);
  int *rest_labels = malloc(sizeof(int) * (n - chosen_k));
  component_labels[0] = min_label;
  unrank_combination(subset_rank, n_rem, chosen_k - 1, remaining_labels,
                     &component_labels[1], rest_labels);

  char *res_A = unrank_e(ctx, child, chosen_k, rank_A, component_labels);

  char *res_Seq = NULL;
  if (n - chosen_k == 0) {
    res_Seq = strdup("Seq()"); // Empty sequence
  } else {
    res_Seq = unrank_e(ctx, seq_expr, n - chosen_k, rank_Seq, rest_labels);
  }

  char *res = malloc(strlen(res_A) + strlen(res_Seq) + 20);
  if (strncmp(res_Seq, "Seq(", 4) == 0) { // Strip
    int lenS = strlen(res_Seq);
    if (lenS > 5) { // Seq(contents)
      res_Seq[lenS - 1] = '\0';
      sprintf(res, "Cycle(%s, %s)", res_A, res_Seq + 4);
    } else { // Seq()
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
  fmpz_clear(current_rank);
  fmpz_clear(total_structs);
  fmpz_clear(subset_rank);
  fmpz_clear(structure_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_Seq);
  return res;
}
