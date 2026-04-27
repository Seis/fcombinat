#include "solver/rank/set.h"

#include <limits.h>
#include <string.h>

#include "solver/math.h"

// Get the minimum label in an object subtree
static int get_min_label(Object *obj) {
  if (obj->type == OBJ_ATOM)
    return obj->label;
  int min = INT_MAX;
  for (int i = 0; i < obj->num_children; i++) {
    int m = get_min_label(obj->children[i]);
    if (m < min)
      min = m;
  }
  return min;
}

// Comparator for sorting Set children by min-label
static int cmp_by_min_label(const void *a, const void *b) {
  Object *oa = *(Object **)a;
  Object *ob = *(Object **)b;
  return get_min_label(oa) - get_min_label(ob);
}

void rank_set_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res, int depth) {
  if (strcmp(obj->name, "Set") != 0) {
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

  // Canonicalize: sort children by min-label so children[0] has the global min
  qsort(obj->children, obj->num_children, sizeof(Object *), cmp_by_min_label);

  Object *Comp_obj = obj->children[0];
  Object *Rest_obj;
  if (obj->num_children == 1) {
    Rest_obj = new_object(OBJ_STRUCT, "Set");
  } else {
    Rest_obj = new_object(OBJ_STRUCT, "Set");
    Rest_obj->num_children = obj->num_children - 1;
    Rest_obj->children = malloc(sizeof(Object *) * Rest_obj->num_children);
    for (int i = 0; i < Rest_obj->num_children; i++) {
      Rest_obj->children[i] = obj->children[i + 1];
    }
  }

  Expr *child_expr = (Expr *)expr->component;

  int k = get_obj_size(Comp_obj);
  int n_total = n;

  int *pool = malloc(sizeof(int) * 128);
  int pool_size = 0;
  int pool_cap = 128;
  collect_labels(obj, &pool, &pool_size, &pool_cap);
  qsort(pool, pool_size, sizeof(int), cmp_int);
  int min_label = pool[0];

  int *labels_Comp = malloc(sizeof(int) * 128);
  int size_C = 0;
  int cap_C = 128;
  collect_labels(Comp_obj, &labels_Comp, &size_C, &cap_C);
  qsort(labels_Comp, size_C, sizeof(int), cmp_int);


  int *pool_rem = &pool[1];
  int *labels_rem = malloc(sizeof(int) * k);
  int idx = 0;
  for (int i = 0; i < size_C; i++)
    if (labels_Comp[i] != min_label)
      labels_rem[idx++] = labels_Comp[i];

  fmpz_t subset_rank;
  fmpz_init(subset_rank);
  rank_combination(subset_rank, n_total - 1, k - 1, pool_rem, labels_rem);

  fmpz_t rank_A;
  fmpz_init(rank_A);
  rank_e(ctx, child_expr, Comp_obj, rank_A, depth);

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
  rank_e(ctx, tail_expr, Rest_obj, rank_Rest, depth);

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
    fmpz_bin_uiui(bin, n_total - 1, i - 1);
    fmpz_mul(weight, count_A, count_S);
    fmpz_mul(weight, weight, bin);
    fmpz_add(res, res, weight);
  }

  if (Rest_obj->num_children > 0)
    free(Rest_obj->children);
  free(Rest_obj);
  free(labels_Comp);
  free(labels_rem);
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

/*
 * rank_set_unlabeled:
 * Inverse of unrank_set_unlabeled.
 * Given a Set object, compute the global rank of each child, sort by global
 * rank (non-decreasing), then compute the rank using count_unl_set_restricted.
 */
void rank_set_unlabeled(Context *ctx, Expr *expr, Object *obj, fmpz_t res, int depth) {
  if (strcmp(obj->name, "Set") != 0) {
    fmpz_set_si(res, -1);
    return;
  }

  int n = get_obj_size(obj);
  if (n == 0) {
    fmpz_zero(res);
    return;
  }

  Expr *child_expr = (Expr *)expr->component;

  /* Build global rank table. */
  int total_g = 0;
  fmpz_t *A_counts = malloc(sizeof(fmpz_t) * (n + 1));
  for (int k = 0; k <= n; k++) {
    fmpz_init(A_counts[k]);
    get_expr_count(A_counts[k], ctx, child_expr, k);
    if (k >= 1)
      total_g += (int)fmpz_get_si(A_counts[k]);
  }

  int *g_to_k = malloc(sizeof(int) * (total_g + 1));
  int *g_to_local_dummy = malloc(sizeof(int) * (total_g + 1));
  {
    int g = 0;
    for (int k = 1; k <= n; k++) {
      int cnt = (int)fmpz_get_si(A_counts[k]);
      for (int r = 0; r < cnt; r++) {
        g_to_k[g] = k;
        g_to_local_dummy[g] = r;
        g++;
      }
    }
  }

  UnlSetCtx uctx = {ctx, child_expr, n, g_to_k, g_to_local_dummy, total_g};

  /* Compute global offset for each size class. */
  int *size_to_goffset = malloc(sizeof(int) * (n + 1));
  size_to_goffset[0] = 0;
  {
    int off = 0;
    for (int k = 1; k <= n; k++) {
      size_to_goffset[k] = off;
      off += (int)fmpz_get_si(A_counts[k]);
    }
  }

  /* Rank each child A-structure and compute its global rank. */
  int nc = obj->num_children;
  int *child_grank = malloc(sizeof(int) * nc);

  for (int i = 0; i < nc; i++) {
    Object *ch = obj->children[i];
    int sz = get_obj_size(ch);
    if (sz < 1 || sz > n) { child_grank[i] = 0; continue; }
    fmpz_t lr;
    fmpz_init(lr);
    rank_e(ctx, child_expr, ch, lr, depth);
    int local_r = (int)fmpz_get_si(lr);
    child_grank[i] = size_to_goffset[sz] + local_r;
    fmpz_clear(lr);
  }

  /* Sort child_grank (non-decreasing). */
  for (int i = 1; i < nc; i++) {
    int key = child_grank[i];
    int j = i - 1;
    while (j >= 0 && child_grank[j] > key) {
      child_grank[j + 1] = child_grank[j];
      j--;
    }
    child_grank[j + 1] = key;
  }

  /* Compute rank by summing blocks skipped before each chosen element. */
  fmpz_zero(res);
  fmpz_t block;
  fmpz_init(block);

  int rem_weight = n;
  int g_prev = 0;
  int cur_restriction = (int)expr->restriction;
  long long cur_limit = expr->limit;

  for (int i = 0; i < nc; i++) {
    int gi = child_grank[i];
    int wi = g_to_k[gi];

    /* Tail restriction after placing element i */
    int tail_restriction = cur_restriction;
    long long tail_limit = cur_limit;
    if (cur_restriction == 2) {
      tail_limit = cur_limit - 1;
    } else if (cur_restriction == 1) {
      tail_limit = cur_limit - 1;
    } else if (cur_restriction == 3) {
      tail_limit = cur_limit - 1;
      if (tail_limit < 0) tail_restriction = 0;
    }

    /* Add counts for global ranks g = g_prev .. gi-1 */
    for (int g = g_prev; g < gi; g++) {
      int wg = g_to_k[g];
      if (wg > rem_weight) continue;
      count_unl_set_restricted(block, &uctx, rem_weight - wg, g,
                               tail_restriction, tail_limit);
      fmpz_add(res, res, block);
    }

    rem_weight -= wi;
    g_prev = gi;
    cur_restriction = tail_restriction;
    cur_limit = tail_limit;
  }

  free(child_grank);
  free(g_to_k);
  free(g_to_local_dummy);
  free(size_to_goffset);
  for (int k = 0; k <= n; k++)
    fmpz_clear(A_counts[k]);
  free(A_counts);
  fmpz_clear(block);
}
