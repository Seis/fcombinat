#include "solver/unrank/set.h"

#include <string.h>
#include "solver/math.h"

char *unrank_set_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                      int *labels, int depth) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return strdup("Set()");

  int min_label = labels[0];

  int *remaining_labels = &labels[1];
  int n_rem = n - 1;

  Expr tail_expr_struct = *expr;
  if (expr->restriction != NONE && expr->limit > 1) {
    tail_expr_struct.limit = expr->limit - 1;
  } else if (expr->restriction != NONE && expr->limit <= 1) {
    tail_expr_struct.restriction = NONE;
    tail_expr_struct.limit = 0;
  }
  Expr *tail_expr = &tail_expr_struct;

  fmpz_t weight;
  fmpz_init(weight);
  fmpz_t count_A;
  fmpz_init(count_A);
  fmpz_t count_S;
  fmpz_init(count_S);
  fmpz_t bin;
  fmpz_init(bin);
  fmpz_t current_rank;
  fmpz_init_set(current_rank, rank);

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
    get_expr_count(count_S, ctx, k_tail_expr, n - k);
    binom_ui(bin, n_rem, k - 1);

    fmpz_mul(weight, count_A, count_S);
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
    fmpz_clear(count_S);
    fmpz_clear(bin);
    fmpz_clear(current_rank);
    return strdup("ErrorSetBounds");
  }

  fmpz_t total_structs;
  fmpz_init(total_structs);
  fmpz_mul(total_structs, count_A, count_S);

  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  fmpz_t structure_rank;
  fmpz_init(structure_rank);

  fmpz_fdiv_qr(subset_rank, structure_rank, current_rank, total_structs);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  fmpz_t rank_S;
  fmpz_init(rank_S);

  fmpz_fdiv_qr(rank_A, rank_S, structure_rank, count_S);

  int *component_labels = malloc(sizeof(int) * chosen_k);
  int *rest_labels = malloc(sizeof(int) * (n - chosen_k));

  component_labels[0] = min_label;

  unrank_combination(subset_rank, n_rem, chosen_k - 1, remaining_labels,
                     &component_labels[1], rest_labels);

  char *res_A = unrank_e(ctx, child, chosen_k, rank_A, component_labels, depth);

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
    res_S = strdup("");
  } else {
    res_S = unrank_e(ctx, final_tail_expr, n - chosen_k, rank_S, rest_labels, depth);
  }

  char *res = malloc(strlen(res_A) + (res_S ? strlen(res_S) : 0) + 20);
  if (n - chosen_k == 0) {
    sprintf(res, "Set(%s)", res_A);
  } else {
    if (strncmp(res_S, "Set(", 4) == 0) {
      res_S[strlen(res_S) - 1] = '\0';
      sprintf(res, "Set(%s, %s)", res_A, res_S + 4);
    } else {
      sprintf(res, "Set(%s, %s)", res_A, res_S);
    }
  }

  free(res_A);
  if (res_S)
    free(res_S);
  free(component_labels);
  free(rest_labels);

  fmpz_clear(weight);
  fmpz_clear(count_A);
  fmpz_clear(count_S);
  fmpz_clear(bin);
  fmpz_clear(current_rank);
  fmpz_clear(total_structs);
  fmpz_clear(subset_rank);
  fmpz_clear(structure_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_S);

  return res;
}

/*
 * unrank_set_unlabeled:
 * Iteratively decode: at each step, find the smallest global rank g0 >= g_min
 * such that count_unl_set_restricted(remaining_w - w(g0), g0, restriction-1) > rank.
 * Accumulate rank offsets, then unrank each A-structure by its local rank.
 */
char *unrank_set_unlabeled(Context *ctx, Expr *expr, int n, fmpz_t rank, int depth) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return strdup("Set()");

  /* Build global rank table. */
  int total_g = 0;
  fmpz_t *A_counts_arr = malloc(sizeof(fmpz_t) * (n + 1));
  for (int k = 0; k <= n; k++) {
    fmpz_init(A_counts_arr[k]);
    get_expr_count(A_counts_arr[k], ctx, child, k);
    if (k >= 1)
      total_g += (int)fmpz_get_si(A_counts_arr[k]);
  }

  int *g_to_k = malloc(sizeof(int) * (total_g + 1));
  int *g_to_local = malloc(sizeof(int) * (total_g + 1));
  {
    int g = 0;
    for (int k = 1; k <= n; k++) {
      int cnt = (int)fmpz_get_si(A_counts_arr[k]);
      for (int r = 0; r < cnt; r++) {
        g_to_k[g] = k;
        g_to_local[g] = r;
        g++;
      }
    }
  }

  UnlSetCtx uctx = {ctx, child, n, g_to_k, g_to_local, total_g};

  fmpz_t current_rank;
  fmpz_init_set(current_rank, rank);
  fmpz_t block;
  fmpz_init(block);

  int *chosen_gs = NULL;
  int num_chosen = 0;
  int cap_chosen = 0;

  int rem_weight = n;
  int g_min = 0;
  int cur_restriction = (int)expr->restriction;
  long long cur_limit = expr->limit;

  while (rem_weight > 0) {
    /* Restrict: must be able to place at least one more element */
    if (cur_restriction == 2 && cur_limit == 0) break; /* EQUAL 0: done */
    if (cur_restriction == 1 && cur_limit <= 0) break; /* LESS 0: done */

    int placed = 0;
    for (int g0 = g_min; g0 < total_g; g0++) {
      int w = g_to_k[g0];
      if (w > rem_weight) continue;

      /* Tail restriction after placing one more element */
      int tail_restriction = cur_restriction;
      long long tail_limit = cur_limit;
      if (cur_restriction == 2) {
        tail_limit = cur_limit - 1;
      } else if (cur_restriction == 1) {
        tail_limit = cur_limit - 1;
        if (tail_limit < 0) continue;
      } else if (cur_restriction == 3) {
        tail_limit = cur_limit - 1;
        if (tail_limit < 0) tail_restriction = 0;
      }

      count_unl_set_restricted(block, &uctx, rem_weight - w, g0,
                               tail_restriction, tail_limit);

      if (fmpz_cmp(current_rank, block) < 0) {
        if (num_chosen >= cap_chosen) {
          cap_chosen = cap_chosen ? cap_chosen * 2 : 8;
          chosen_gs = realloc(chosen_gs, sizeof(int) * cap_chosen);
        }
        chosen_gs[num_chosen++] = g0;
        rem_weight -= w;
        g_min = g0; /* repetition allowed */
        cur_restriction = tail_restriction;
        cur_limit = tail_limit;
        placed = 1;
        break;
      }
      fmpz_sub(current_rank, current_rank, block);
    }
    if (!placed) break;
    if (rem_weight == 0) break;
  }

  /* Unrank each chosen A-structure. */
  char **parts = malloc(sizeof(char *) * (num_chosen + 1));
  int total_len = 5;

  for (int i = 0; i < num_chosen; i++) {
    int gi = chosen_gs[i];
    int sz = g_to_k[gi];
    int lr = g_to_local[gi];
    fmpz_t local_rank;
    fmpz_init_set_si(local_rank, lr);
    parts[i] = unrank_e(ctx, child, sz, local_rank, NULL, depth);
    fmpz_clear(local_rank);
    total_len += strlen(parts[i]) + 2;
  }

  char *res = malloc(total_len + 10);
  char *p = res;
  p += sprintf(p, "Set(");
  for (int i = 0; i < num_chosen; i++) {
    if (i > 0) { *p++ = ','; *p++ = ' '; }
    p += sprintf(p, "%s", parts[i]);
    free(parts[i]);
  }
  *p++ = ')';
  *p = '\0';

  free(parts);
  if (chosen_gs) free(chosen_gs);
  free(g_to_k);
  free(g_to_local);
  for (int k = 0; k <= n; k++)
    fmpz_clear(A_counts_arr[k]);
  free(A_counts_arr);
  fmpz_clear(current_rank);
  fmpz_clear(block);

  return res;
}
