#include "solver/unrank/pset.h"

#include <string.h>
#include "solver/math.h"

char *unrank_powerset_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                           int *labels, int depth) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return strdup("PowerSet()");

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
  fmpz_t current_rank;
  fmpz_init_set(current_rank, rank);

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
    fmpz_clear(current_rank);
    return strdup("ErrorPowerSetBounds");
  }

  fmpz_t total_structs;
  fmpz_init(total_structs);
  fmpz_mul(total_structs, count_A, count_P);

  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  fmpz_t structure_rank;
  fmpz_init(structure_rank);

  fmpz_fdiv_qr(subset_rank, structure_rank, current_rank, total_structs);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  fmpz_t rank_P;
  fmpz_init(rank_P);

  fmpz_fdiv_qr(rank_A, rank_P, structure_rank, count_P);

  int *component_labels = malloc(sizeof(int) * chosen_k);
  int *rest_labels = malloc(sizeof(int) * (n - chosen_k));
  component_labels[0] = min_label;
  unrank_combination(subset_rank, n_rem, chosen_k - 1, remaining_labels,
                     &component_labels[1], rest_labels);

  char *res_A = unrank_e(ctx, child, chosen_k, rank_A, component_labels, depth);

  char *res_P = NULL;
  if (n - chosen_k == 0) {
    res_P = strdup("");
  } else {
    res_P = unrank_e(ctx, tail_expr, n - chosen_k, rank_P, rest_labels, depth);
  }

  char *res = malloc(strlen(res_A) + (res_P ? strlen(res_P) : 0) + 20);
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
  fmpz_clear(current_rank);
  fmpz_clear(total_structs);
  fmpz_clear(subset_rank);
  fmpz_clear(structure_rank);
  fmpz_clear(rank_A);
  fmpz_clear(rank_P);

  return res;
}

/*
 * unrank_pset_unlabeled:
 * Like unrank_set_unlabeled but each A-structure may appear at most once
 * (PowerSet = set of distinct A-structures, strictly increasing global ranks).
 * Handles cardinality restrictions via count_unl_pset_restricted.
 */
char *unrank_pset_unlabeled(Context *ctx, Expr *expr, int n, fmpz_t rank, int depth) {
  Expr *child = (Expr *)expr->component;
  if (n == 0)
    return strdup("PowerSet()");

  /* Build global rank table. */
  int total_g = 0;
  fmpz_t *A_counts = malloc(sizeof(fmpz_t) * (n + 1));
  for (int k = 0; k <= n; k++) {
    fmpz_init(A_counts[k]);
    get_expr_count(A_counts[k], ctx, child, k);
    if (k >= 1)
      total_g += (int)fmpz_get_si(A_counts[k]);
  }

  int *g_to_k = malloc(sizeof(int) * (total_g + 1));
  int *g_to_local = malloc(sizeof(int) * (total_g + 1));
  {
    int g = 0;
    for (int k = 1; k <= n; k++) {
      int cnt = (int)fmpz_get_si(A_counts[k]);
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
    if (cur_restriction == 2 && cur_limit == 0) break;
    if (cur_restriction == 1 && cur_limit <= 0) break;

    int placed = 0;
    for (int g0 = g_min; g0 < total_g; g0++) {
      int w = g_to_k[g0];
      if (w > rem_weight) continue;

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

      /* Next element must have strictly greater global rank */
      count_unl_pset_restricted(block, &uctx, rem_weight - w, g0 + 1,
                                tail_restriction, tail_limit);

      if (fmpz_cmp(current_rank, block) < 0) {
        if (num_chosen >= cap_chosen) {
          cap_chosen = cap_chosen ? cap_chosen * 2 : 8;
          chosen_gs = realloc(chosen_gs, sizeof(int) * cap_chosen);
        }
        chosen_gs[num_chosen++] = g0;
        rem_weight -= w;
        g_min = g0 + 1;
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

  /* Unrank each A-structure. */
  char **parts = malloc(sizeof(char *) * (num_chosen + 1));
  int total_len = 12;

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
  p += sprintf(p, "PowerSet(");
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
    fmpz_clear(A_counts[k]);
  free(A_counts);
  fmpz_clear(current_rank);
  fmpz_clear(block);

  return res;
}
