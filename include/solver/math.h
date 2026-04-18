#ifndef SOLVER_MATH_H
#define SOLVER_MATH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

/* Abort-on-failure allocators for internal solver code. */
static inline void *xmalloc(size_t n) {
  void *p = malloc(n ? n : 1);
  if (!p) { fprintf(stderr, "fcombinat: out of memory\n"); abort(); }
  return p;
}
static inline void *xrealloc(void *ptr, size_t n) {
  void *p = realloc(ptr, n ? n : 1);
  if (!p) { fprintf(stderr, "fcombinat: out of memory\n"); abort(); }
  return p;
}
static inline char *xstrdup(const char *s) {
  if (!s) return NULL;
  char *d = strdup(s);
  if (!d) { fprintf(stderr, "fcombinat: out of memory\n"); abort(); }
  return d;
}

// Helper structs
typedef struct {
  fmpq_poly_t *coeffs; // Array of series
  int limit;           // Max degree in u
} SeriesPoly;

// Common Helpers
void sp_init(SeriesPoly *sp, int limit);
void sp_clear(SeriesPoly *sp);
void sp_mul(SeriesPoly *res, SeriesPoly *A, SeriesPoly *B, int n);
void sp_exp(SeriesPoly *res, SeriesPoly *A, int n);

void counts_to_egf(fmpq_poly_t res, fmpz_poly_t counts, int n);
void egf_to_counts(fmpz_poly_t res, fmpq_poly_t egf, int n);
void binomial_convolution(fmpz_poly_t res, fmpz_poly_t A, fmpz_poly_t B, int n);

// Combinatorics / Combination Helpers
void binom_ui(fmpz_t res, int n, int k);
void rank_combination(fmpz_t res, int n, int k, int *labels, int *selected);
void unrank_combination(fmpz_t rank, int n, int k, int *labels, int *selected,
                        int *complement);

int cmp_int(const void *a, const void *b);

// Unlabeled helpers

/* Rank a sorted c-tuple drawn (with repetition) from {0..M-1} using
   the combinatorial number system for multiset combinations. */
void rank_multiset_comb(fmpz_t res, int M, int c, fmpz_t *sorted_ranks);

/* Unrank: fills out_ranks[0..c-1] (ascending) from {0..M-1}. */
void unrank_multiset_comb(fmpz_t rank, int M, int c, fmpz_t *out_ranks);

/* Count canonical unlabeled Sets of total weight n using A-structures
   with global rank >= g_min.
   A_counts[k] = f_A(k) for k=0..max_k
   g_to_k[g] = size of A-structure with global rank g
   total_g   = total number of A-structures across all sizes */
void count_canonical_set_dp(fmpz_t res, int n, int g_min,
                             fmpz_t *A_counts, int *g_to_k, int total_g);

/* Same but each A-structure may appear at most once (for PowerSet). */
void count_canonical_pset_dp(fmpz_t res, int n, int g_min,
                              fmpz_t *A_counts, int *g_to_k, int total_g);

/*
 * Context for unlabeled Set rank/unrank.
 * Holds the global rank table (g_to_k, g_to_local) and the Context
 * pointer so we can call count functions.
 */
typedef struct {
  void    *ctx;       /* actually Context*, avoid circular include */
  void    *child_expr; /* actually Expr* */
  int      n;
  int     *g_to_k;
  int     *g_to_local;
  int      total_g;
} UnlSetCtx;

/*
 * Count canonical multisets (with repetition) of weight n,
 * using A-structures with global rank >= g_min, subject to a
 * cardinality restriction (restriction and limit match Restriction enum values:
 *   NONE=0, LESS=1, EQUAL=2, GREATER=3).
 */
void count_unl_set_restricted(fmpz_t res, UnlSetCtx *uctx,
                               int n, int g_min,
                               int restriction, long long limit);

/*
 * Same but each A-structure may appear at most once (PowerSet).
 */
void count_unl_pset_restricted(fmpz_t res, UnlSetCtx *uctx,
                                int n, int g_min,
                                int restriction, long long limit);

#endif
