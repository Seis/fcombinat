#include "solver/math.h"


// --- Helpers ---

void counts_to_egf(fmpq_poly_t res, fmpz_poly_t counts, int n) {
  fmpz_t c;
  fmpz_init(c);
  fmpq_t val;
  fmpq_init(val);
  fmpz_t fact;
  fmpz_init(fact);

  for (int i = 0; i <= n; i++) {
    fmpz_poly_get_coeff_fmpz(c, counts, i);
    if (fmpz_is_zero(c)) {
      fmpq_poly_set_coeff_ui(res, i, 0);
      continue;
    }

    fmpz_fac_ui(fact, i);
    fmpq_set_fmpz_frac(val, c, fact);
    fmpq_poly_set_coeff_fmpq(res, i, val);
  }

  fmpz_clear(c);
  fmpq_clear(val);
  fmpz_clear(fact);
}

void egf_to_counts(fmpz_poly_t res, fmpq_poly_t egf, int n) {
  fmpq_t c;
  fmpq_init(c);
  fmpz_t fact;
  fmpz_init(fact);

  fmpz_poly_zero(res);

  for (int i = 0; i <= n; i++) {
    fmpq_poly_get_coeff_fmpq(c, egf, i);
    if (fmpq_is_zero(c))
      continue;

    fmpz_fac_ui(fact, i);
    fmpq_mul_fmpz(c, c, fact); // c = c * i!

    if (fmpz_is_one(fmpq_denref(c))) {
      fmpz_poly_set_coeff_fmpz(res, i, fmpq_numref(c));
    }
  }

  fmpq_clear(c);
  fmpz_clear(fact);
}

void binomial_convolution(fmpz_poly_t res, fmpz_poly_t A, fmpz_poly_t B,
                          int n) {
  fmpq_poly_t Aq;
  fmpq_poly_init(Aq);
  fmpq_poly_t Bq;
  fmpq_poly_init(Bq);
  fmpq_poly_t Cq;
  fmpq_poly_init(Cq);

  counts_to_egf(Aq, A, n);
  counts_to_egf(Bq, B, n);

  fmpq_poly_mullow(Cq, Aq, Bq, n + 1);

  egf_to_counts(res, Cq, n);

  fmpq_poly_clear(Aq);
  fmpq_poly_clear(Bq);
  fmpq_poly_clear(Cq);
}

// --- SeriesPoly Logic ---

void sp_init(SeriesPoly *sp, int limit) {
  sp->limit = limit;
  sp->coeffs = malloc(sizeof(fmpq_poly_t) * (limit + 1));
  for (int i = 0; i <= limit; i++) {
    fmpq_poly_init(sp->coeffs[i]);
  }
}

void sp_clear(SeriesPoly *sp) {
  for (int i = 0; i <= sp->limit; i++) {
    fmpq_poly_clear(sp->coeffs[i]);
  }
  free(sp->coeffs);
}

void sp_mul(SeriesPoly *res, SeriesPoly *A, SeriesPoly *B, int n) {
  fmpq_poly_t tmp;
  fmpq_poly_init(tmp);

  // Clear Res first (Assume Res != A, Res != B)
  for (int i = 0; i <= res->limit; i++) {
    fmpq_poly_zero(res->coeffs[i]);

    for (int j = 0; j <= i; j++) {
      if (j > A->limit || (i - j) > B->limit)
        continue;

      fmpq_poly_mullow(tmp, A->coeffs[j], B->coeffs[i - j], n + 1);
      fmpq_poly_add(res->coeffs[i], res->coeffs[i], tmp);
    }
  }
  fmpq_poly_clear(tmp);
}

void sp_exp(SeriesPoly *res, SeriesPoly *A, int n) {
  for (int i = 0; i <= res->limit; i++)
    fmpq_poly_zero(res->coeffs[i]);
  fmpq_poly_set_coeff_ui(res->coeffs[0], 0, 1);

  SeriesPoly term;
  sp_init(&term, res->limit);
  fmpq_poly_set_coeff_ui(term.coeffs[0], 0, 1);

  SeriesPoly next_term;
  sp_init(&next_term, res->limit);

  fmpq_t fact;
  fmpq_init(fact);

  for (int k = 1; k <= res->limit; k++) {
    sp_mul(&next_term, &term, A, n);

    for (int i = 0; i <= res->limit; i++)
      fmpq_poly_set(term.coeffs[i], next_term.coeffs[i]);

    fmpq_set_si(fact, 1, 1);
    fmpz_t f_z;
    fmpz_init(f_z);
    for (int f = 1; f <= k; f++) {
      fmpz_set_si(f_z, f);
      fmpq_div_fmpz(fact, fact, f_z);
    }
    fmpz_clear(f_z);

    for (int i = 0; i <= res->limit; i++) {
      fmpq_poly_scalar_mul_fmpq(next_term.coeffs[i], next_term.coeffs[i], fact);
      fmpq_poly_add(res->coeffs[i], res->coeffs[i], next_term.coeffs[i]);
    }
  }

  fmpq_clear(fact);
  sp_clear(&term);
  sp_clear(&next_term);
}

void binom_ui(fmpz_t res, int n, int k) {
  if (k < 0 || k > n) {
    fmpz_zero(res);
    return;
  }
  fmpz_bin_uiui(res, n, k);
}

void rank_combination(fmpz_t res, int n, int k, int *labels, int *selected) {
  fmpz_zero(res);
  int sel_idx = 0;
  int cur_k = k;
  int rem_n = n;

  fmpz_t bin;
  fmpz_init(bin);

  for (int i = 0; i < n; i++) {
    if (cur_k == 0)
      break;

    if (labels[i] == selected[sel_idx]) {
      sel_idx++;
      cur_k--;
    } else {
      fmpz_bin_uiui(bin, rem_n - 1, cur_k - 1);
      fmpz_add(res, res, bin);
    }
    rem_n--;
  }
  fmpz_clear(bin);
}

void unrank_combination(fmpz_t rank, int n, int k, int *labels, int *selected,
                        int *complement) {
  if (k == 0) {
    for (int i = 0; i < n; i++)
      complement[i] = labels[i];
    return;
  }
  if (k == n) {
    for (int i = 0; i < n; i++)
      selected[i] = labels[i];
    return;
  }

  fmpz_t r;
  fmpz_init_set(r, rank);
  fmpz_t threshold;
  fmpz_init(threshold);

  int s_idx = 0;
  int c_idx = 0;

  int current_k = k;
  int remaining_n = n;

  for (int i = 0; i < n; i++) {
    if (current_k > 0) {
      binom_ui(threshold, remaining_n - 1, current_k - 1);

      if (fmpz_cmp(r, threshold) < 0) {
        selected[s_idx++] = labels[i];
        current_k--;
      } else {
        complement[c_idx++] = labels[i];
        fmpz_sub(r, r, threshold);
      }
    } else {
      complement[c_idx++] = labels[i];
    }
    remaining_n--;
  }

  fmpz_clear(r);
  fmpz_clear(threshold);
}

int cmp_int(const void *a, const void *b) { return (*(int *)a - *(int *)b); }

/* ---------- Unlabeled helpers ---------- */

/*
 * rank_multiset_comb: rank a sorted c-tuple (ascending, with repetition)
 * drawn from {0..M-1} using the multiset combinatorial number system.
 *
 * The bijection maps a sorted tuple (a_0 <= a_1 <= ... <= a_{c-1}) where
 * each a_i in {0..M-1} to a rank in [0, C(M+c-1, c) - 1].
 *
 * We transform: b_i = a_i + i  (so b_0 < b_1 < ... < b_{c-1}, each in
 * {0..M+c-2}) then apply the standard combinatorial number system for sets.
 */
void rank_multiset_comb(fmpz_t res, int M, int c, fmpz_t *sorted_ranks) {
  fmpz_zero(res);
  if (c == 0 || M == 0)
    return;

  fmpz_t bin;
  fmpz_init(bin);
  fmpz_t b;
  fmpz_init(b);

  /* Transform to strict sequence and apply combinatorial number system. */
  for (int i = 0; i < c; i++) {
    /* b_i = sorted_ranks[i] + i */
    fmpz_add_ui(b, sorted_ranks[i], (unsigned long)i);
    int bi = (int)fmpz_get_si(b);
    /* contribution: C(bi, i+1) */
    fmpz_bin_uiui(bin, (unsigned long)bi, (unsigned long)(i + 1));
    fmpz_add(res, res, bin);
  }

  fmpz_clear(bin);
  fmpz_clear(b);
}

/*
 * unrank_multiset_comb: fills out_ranks[0..c-1] (ascending) given rank in
 * [0, C(M+c-1,c)-1].
 */
void unrank_multiset_comb(fmpz_t rank, int M, int c, fmpz_t *out_ranks) {
  if (c == 0)
    return;

  fmpz_t r;
  fmpz_init_set(r, rank);
  fmpz_t bin;
  fmpz_init(bin);

  /* Decode combinatorial number system: find b_{c-1} >= b_{c-2} >= ... */
  int N = M + c - 1; /* values of b_i range in {0..N-1} */

  for (int i = c - 1; i >= 0; i--) {
    /* Find largest bi such that C(bi, i+1) <= r */
    int bi = i; /* minimum value */
    while (1) {
      fmpz_bin_uiui(bin, (unsigned long)(bi + 1), (unsigned long)(i + 1));
      if (fmpz_cmp(bin, r) > 0)
        break;
      bi++;
      if (bi >= N)
        break;
    }
    /* b_i = bi, contribution C(bi, i+1) */
    fmpz_bin_uiui(bin, (unsigned long)bi, (unsigned long)(i + 1));
    fmpz_sub(r, r, bin);
    /* a_i = b_i - i */
    fmpz_set_si(out_ranks[i], bi - i);
  }

  fmpz_clear(r);
  fmpz_clear(bin);
}

/*
 * count_canonical_set_dp:
 * Count canonical multisets of total weight n using A-structures
 * with global rank >= g_min.
 *
 * We iterate global ranks g from g_min upward.  For each g, the
 * A-structure has weight w = g_to_k[g] and f_A(w) = A_counts[w].
 * We DP over (remaining weight, min allowed global rank).
 *
 * DP state: dp[weight_remaining] = number of canonical multisets of that
 * remaining weight using global ranks >= current g.
 *
 * We process global ranks in increasing order.  For rank g with weight w,
 * we pick multiplicity m=0,1,2,... until m*w > remaining weight, and add.
 *
 * This is a forward DP: initialise dp[0]=1, dp[j]=0 for j>0.
 * For each global rank g (in order from g_min), let w = g_to_k[g].
 * The number of A-structures of weight w with rank exactly g within their
 * weight class is 1 (each global rank is one A-structure).
 * So we treat it like a coin of denomination w with unlimited supply.
 *
 * Wait — each global rank corresponds to exactly one A-structure.  The
 * "canonical" condition says: in the sorted multiset, ties on weight are
 * broken by global rank.  So the DP must track the *last used global rank*
 * to enforce non-decreasing order.
 *
 * Simpler recursive formulation:
 *   f(n, g) = sum over m=0,1,... : m*w <= n, A-struct at rank g contributes m
 *           = 1 * f(n - m*w, g+1)   for each valid m
 * But this is O(n * total_g * max_mult) which is fine for small n.
 *
 * We implement this with memoisation.
 */
void count_canonical_set_dp(fmpz_t res, int n, int g_min,
                             fmpz_t *A_counts, int *g_to_k, int total_g) {
  /* Base case */
  if (n == 0) {
    fmpz_one(res);
    return;
  }
  if (g_min >= total_g) {
    fmpz_zero(res);
    return;
  }

  fmpz_zero(res);
  fmpz_t sub;
  fmpz_init(sub);

  for (int g = g_min; g < total_g; g++) {
    int w = g_to_k[g];
    if (w <= 0 || w > n)
      continue;
    /* Try multiplicities m = 1, 2, ... as long as m*w <= n */
    for (int m = 1; m * w <= n; m++) {
      count_canonical_set_dp(sub, n - m * w, g, A_counts, g_to_k, total_g);
      fmpz_add(res, res, sub);
    }
  }

  fmpz_clear(sub);
}

/*
 * count_canonical_pset_dp: same but multiplicity is 0 or 1.
 */
void count_canonical_pset_dp(fmpz_t res, int n, int g_min,
                              fmpz_t *A_counts, int *g_to_k, int total_g) {
  if (n == 0) {
    fmpz_one(res);
    return;
  }
  if (g_min >= total_g) {
    fmpz_zero(res);
    return;
  }

  fmpz_zero(res);
  fmpz_t sub;
  fmpz_init(sub);

  for (int g = g_min; g < total_g; g++) {
    int w = g_to_k[g];
    if (w <= 0 || w > n)
      continue;
    /* multiplicity 1 only */
    count_canonical_pset_dp(sub, n - w, g + 1, A_counts, g_to_k, total_g);
    fmpz_add(res, res, sub);
  }

  fmpz_clear(sub);
}

void count_unl_set_restricted(fmpz_t res, UnlSetCtx *uctx,
                               int n, int g_min,
                               int restriction, long long limit) {
  /* Base cases */
  if (n == 0) {
    switch (restriction) {
    case 0: /* NONE */    fmpz_one(res); return;
    case 1: /* LESS */    fmpz_set_si(res, (limit >= 0) ? 1 : 0); return;
    case 2: /* EQUAL */   fmpz_set_si(res, (limit == 0) ? 1 : 0); return;
    case 3: /* GREATER */ fmpz_set_si(res, (limit <= 0) ? 1 : 0); return;
    default: fmpz_zero(res); return;
    }
  }
  if (restriction == 2 && limit == 0) { fmpz_zero(res); return; }
  if (restriction == 1 && limit <= 0) { fmpz_zero(res); return; }
  if (g_min >= uctx->total_g)         { fmpz_zero(res); return; }

  fmpz_zero(res);
  fmpz_t sub;
  fmpz_init(sub);

  for (int g = g_min; g < uctx->total_g; g++) {
    int w = uctx->g_to_k[g];
    if (w > n) continue;

    int tail_restriction = restriction;
    long long tail_limit = limit;
    if (restriction == 2) {
      tail_limit = limit - 1;
    } else if (restriction == 1) {
      tail_limit = limit - 1;
      if (tail_limit < 0) continue;
    } else if (restriction == 3) {
      tail_limit = limit - 1;
      if (tail_limit < 0) tail_restriction = 0;
    }

    count_unl_set_restricted(sub, uctx, n - w, g,
                             tail_restriction, tail_limit);
    fmpz_add(res, res, sub);
  }

  fmpz_clear(sub);
}

void count_unl_pset_restricted(fmpz_t res, UnlSetCtx *uctx,
                                int n, int g_min,
                                int restriction, long long limit) {
  if (n == 0) {
    switch (restriction) {
    case 0: /* NONE */    fmpz_one(res); return;
    case 1: /* LESS */    fmpz_set_si(res, (limit >= 0) ? 1 : 0); return;
    case 2: /* EQUAL */   fmpz_set_si(res, (limit == 0) ? 1 : 0); return;
    case 3: /* GREATER */ fmpz_set_si(res, (limit <= 0) ? 1 : 0); return;
    default: fmpz_zero(res); return;
    }
  }
  if (restriction == 2 && limit == 0) { fmpz_zero(res); return; }
  if (restriction == 1 && limit <= 0) { fmpz_zero(res); return; }
  if (g_min >= uctx->total_g)         { fmpz_zero(res); return; }

  fmpz_zero(res);
  fmpz_t sub;
  fmpz_init(sub);

  for (int g = g_min; g < uctx->total_g; g++) {
    int w = uctx->g_to_k[g];
    if (w > n) continue;

    int tail_restriction = restriction;
    long long tail_limit = limit;
    if (restriction == 2) {
      tail_limit = limit - 1;
    } else if (restriction == 1) {
      tail_limit = limit - 1;
      if (tail_limit < 0) continue;
    } else if (restriction == 3) {
      tail_limit = limit - 1;
      if (tail_limit < 0) tail_restriction = 0;
    }

    /* PowerSet: next element must have strictly greater global rank */
    count_unl_pset_restricted(sub, uctx, n - w, g + 1,
                              tail_restriction, tail_limit);
    fmpz_add(res, res, sub);
  }

  fmpz_clear(sub);
}
