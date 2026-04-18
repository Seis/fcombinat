#include "solver/unrank/seq.h"


#include "solver/math.h"

char *unrank_seq_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                      int *labels) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return strdup("Seq()");

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

  Expr tail_expr_struct = *expr;
  if (expr->restriction != NONE && expr->limit > 1) {
    tail_expr_struct.limit = expr->limit - 1;
  } else if (expr->restriction != NONE && expr->limit <= 1) {
    tail_expr_struct.restriction = NONE;
    tail_expr_struct.limit = 0;
  }
  Expr *tail_expr = &tail_expr_struct;

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
    fmpz_clear(current_rank);
    return strdup("ErrorSeqBounds");
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
  unrank_combination(subset_rank, n, chosen_k, labels, component_labels,
                     rest_labels);

  char *res_A = unrank_e(ctx, child, chosen_k, rank_A, component_labels);

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
    res_S = unrank_e(ctx, final_tail_expr, n - chosen_k, rank_Seq, rest_labels);
  }

  char *res = malloc(strlen(res_A) + strlen(res_S) + 20);
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
  fmpz_clear(current_rank);
  fmpz_clear(total_structs);
  fmpz_clear(subset_rank);
  fmpz_clear(structure_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_Seq);
  return res;
}

char *unrank_seq_unlabeled(Context *ctx, Expr *expr, int n, fmpz_t rank) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return strdup("Seq()");

  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_Seq;
  fmpz_init(count_Seq);
  fmpz_t block;
  fmpz_init(block);
  fmpz_t current_rank;
  fmpz_init_set(current_rank, rank);

  /* Build the tail expression (same seq, possibly with limit decremented). */
  Expr tail_expr_struct = *expr;
  if (expr->restriction != NONE && expr->limit > 1) {
    tail_expr_struct.limit = expr->limit - 1;
  } else if (expr->restriction != NONE && expr->limit <= 1) {
    tail_expr_struct.restriction = NONE;
    tail_expr_struct.limit = 0;
  }
  Expr *tail_expr = &tail_expr_struct;

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
    fmpz_mul(block, count_A, count_Seq);

    if (fmpz_cmp(current_rank, block) < 0) {
      chosen_k = k;
      break;
    }
    fmpz_sub(current_rank, current_rank, block);
  }

  if (chosen_k == -1) {
    fmpz_clear(count_A);
    fmpz_clear(count_Seq);
    fmpz_clear(block);
    fmpz_clear(current_rank);
    return strdup("ErrorSeqUnlabeledBounds");
  }

  get_expr_count(count_A, ctx, child, chosen_k);

  Expr final_tail_expr_struct = *expr;
  if (expr->restriction != NONE && expr->limit > 1) {
    final_tail_expr_struct.limit = expr->limit - 1;
  } else if (expr->restriction != NONE && expr->limit <= 1) {
    final_tail_expr_struct.restriction = NONE;
    final_tail_expr_struct.limit = 0;
  }
  Expr *final_tail_expr = &final_tail_expr_struct;
  get_expr_count(count_Seq, ctx, final_tail_expr, n - chosen_k);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  fmpz_t rank_Seq;
  fmpz_init(rank_Seq);
  /* current_rank = rank_A * count_Seq + rank_Seq */
  fmpz_fdiv_qr(rank_A, rank_Seq, current_rank, count_Seq);

  char *res_A = unrank_e(ctx, child, chosen_k, rank_A, NULL);

  char *res_S = NULL;
  if (n - chosen_k == 0) {
    res_S = strdup("Seq()");
  } else {
    res_S = unrank_e(ctx, final_tail_expr, n - chosen_k, rank_Seq, NULL);
  }

  char *res = malloc(strlen(res_A) + strlen(res_S) + 20);
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
  fmpz_clear(count_A);
  fmpz_clear(count_Seq);
  fmpz_clear(block);
  fmpz_clear(current_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_Seq);
  return res;
}
