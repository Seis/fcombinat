#include "solver/count/seq.h"


#include "solver/count/eval.h"
#include "solver/math.h"

void compute_seq_expr(Context *ctx, Expr *child, Expr *expr_full,
                      fmpz_poly_t res, int n, int is_labeled) {
  fmpz_poly_t A;
  fmpz_poly_init(A);
  compute_e(ctx, child, A, n, is_labeled);

  fmpz_poly_set_coeff_ui(A, 0, 0);

  if (is_labeled) {
    // Labeled Sequence
    // A (Counts) -> Aq (EGF)
    fmpq_poly_t Aq;
    fmpq_poly_init(Aq);
    counts_to_egf(Aq, A, n);

    fmpq_poly_t Rq;
    fmpq_poly_init(Rq);

    if (expr_full->restriction == NONE) {
      // 1/(1-Aq)
      fmpq_poly_t denom;
      fmpq_poly_init(denom);
      fmpq_poly_zero(denom);
      fmpq_poly_set_coeff_ui(denom, 0, 1);
      fmpq_poly_sub(denom, denom, Aq);
      fmpq_poly_inv_series(Rq, denom, n + 1);
      fmpq_poly_clear(denom);
    } else {
      // Bounded Sequence
      // Sum A^k
      long long limit = expr_full->limit;
      if (limit < 0)
        limit = 0;

      fmpq_poly_zero(Rq);
      fmpq_poly_t term;
      fmpq_poly_init(term);
      fmpq_poly_set_coeff_ui(term, 0, 1); // A^0

      int max_idx = limit;
      if (expr_full->restriction == GREATER)
        max_idx = limit - 1; // Subtract logic

      fmpq_poly_t sum;
      fmpq_poly_init(sum);
      fmpq_poly_zero(sum);

      // Optimization: iterative power
      for (int k = 0; k <= max_idx; k++) {
        if (expr_full->restriction == LESS ||
            expr_full->restriction == GREATER ||
            (expr_full->restriction == EQUAL && k == limit)) {
          fmpq_poly_add(sum, sum, term);
        }

        if (k < max_idx) {
          fmpq_poly_mullow(term, term, Aq, n + 1);
        }
      }

      if (expr_full->restriction == LESS || expr_full->restriction == EQUAL) {
        fmpq_poly_set(Rq, sum);
      } else if (expr_full->restriction == GREATER) {
        // Unrestricted - Less(k-1)
        fmpq_poly_t denom;
        fmpq_poly_init(denom);
        fmpq_poly_zero(denom);
        fmpq_poly_set_coeff_ui(denom, 0, 1);
        fmpq_poly_sub(denom, denom, Aq);
        fmpq_poly_inv_series(Rq, denom, n + 1);
        fmpq_poly_sub(Rq, Rq, sum);
        fmpq_poly_clear(denom);
      }
      fmpq_poly_clear(sum);
      fmpq_poly_clear(term);
    }

    egf_to_counts(res, Rq, n);
    fmpq_poly_clear(Aq);
    fmpq_poly_clear(Rq);
  } else {
    // Unlabeled Sequence
    if (expr_full->restriction == NONE) {
      fmpz_poly_t denom;
      fmpz_poly_init(denom);
      fmpz_poly_zero(denom);
      fmpz_poly_set_coeff_ui(denom, 0, 1);
      fmpz_poly_sub(denom, denom, A);

      fmpz_t d0;
      fmpz_init(d0);
      fmpz_poly_get_coeff_fmpz(d0, denom, 0);
      if (fmpz_is_zero(d0)) {
        fprintf(stderr,
                "DEBUG: Sequence DivZero! A(0)=1. Denom(0)=0. Labeled=%d\n",
                is_labeled);
      }
      fmpz_clear(d0);

      fmpz_poly_inv_series(res, denom, n + 1);
      fmpz_poly_clear(denom);
    } else {
      // Restricted Unlabeled Sequence
      fmpz_poly_zero(res);
      fmpz_poly_t term;
      fmpz_poly_init(term);
      fmpz_poly_set_coeff_ui(term, 0, 1); // A^0
      fmpz_poly_t sum;
      fmpz_poly_init(sum);
      fmpz_poly_zero(sum);

      long long limit = expr_full->limit;
      if (limit < 0)
        limit = 0;
      int max_idx = limit;
      if (expr_full->restriction == GREATER)
        max_idx = limit - 1;

      for (int k = 0; k <= max_idx; k++) {
        if (expr_full->restriction == LESS ||
            expr_full->restriction == GREATER ||
            (expr_full->restriction == EQUAL && k == limit)) {
          fmpz_poly_add(sum, sum, term);
        }
        if (k < max_idx) {
          fmpz_poly_mullow(term, term, A, n + 1);
        }
      }

      if (expr_full->restriction == LESS || expr_full->restriction == EQUAL) {
        fmpz_poly_set(res, sum);
      } else if (expr_full->restriction == GREATER) {
        fmpz_poly_t denom;
        fmpz_poly_init(denom);
        fmpz_poly_zero(denom);
        fmpz_poly_set_coeff_ui(denom, 0, 1);
        fmpz_poly_sub(denom, denom, A);
        fmpz_poly_inv_series(res, denom, n + 1);
        fmpz_poly_sub(res, res, sum);
        fmpz_poly_clear(denom);
      }

      fmpz_poly_clear(term);
      fmpz_poly_clear(sum);
    }
  }

  fmpz_poly_clear(A);
}
