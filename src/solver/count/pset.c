#include "solver/count/pset.h"


#include "solver/count/eval.h"
#include "solver/count/mset.h"
#include "solver/math.h"

void compute_powerset_unlabeled(fmpz_poly_t res, fmpz_poly_t A, int n) {
  fmpz_poly_t A_z;
  fmpz_poly_init(A_z);
  fmpz_poly_set(A_z, A);

  fmpz_t a0;
  fmpz_init(a0);
  fmpz_poly_get_coeff_fmpz(a0, A_z, 0);
  if (!fmpz_is_zero(a0)) {
    fmpz_poly_set_coeff_ui(A_z, 0, 0);
  }

  fmpq_poly_t sum_poly;
  fmpq_poly_init(sum_poly);
  fmpq_poly_zero(sum_poly);

  for (int k = 1; k <= n; k++) {
    fmpq_poly_t term;
    fmpq_poly_init(term);

    for (int i = 1; i <= n / k; i++) {
      fmpz_t coeff_val;
      fmpz_init(coeff_val);
      fmpz_poly_get_coeff_fmpz(coeff_val, A_z, i);
      if (!fmpz_is_zero(coeff_val)) {
        fmpq_poly_set_coeff_fmpz(term, i * k, coeff_val);
      }
      fmpz_clear(coeff_val);
    }

    fmpq_t scale;
    fmpq_init(scale);
    if ((k - 1) % 2 == 0)
      fmpq_set_si(scale, 1, k);
    else
      fmpq_set_si(scale, -1, k);

    fmpq_poly_scalar_mul_fmpq(term, term, scale);
    fmpq_poly_add(sum_poly, sum_poly, term);

    fmpq_clear(scale);
    fmpq_poly_clear(term);
  }

  fmpq_poly_t result_rational;
  fmpq_poly_init(result_rational);
  fmpq_poly_exp_series(result_rational, sum_poly, n + 1);

  fmpz_poly_zero(res);
  fmpq_t q_coeff;
  fmpq_init(q_coeff);
  fmpz_t z_coeff;
  fmpz_init(z_coeff);

  for (int i = 0; i <= n; i++) {
    fmpq_poly_get_coeff_fmpq(q_coeff, result_rational, i);
    fmpz_set(z_coeff, fmpq_numref(q_coeff));
    fmpz_poly_set_coeff_fmpz(res, i, z_coeff);
  }

  fmpq_clear(q_coeff);
  fmpz_clear(z_coeff);
  fmpz_clear(a0);
  fmpz_poly_clear(A_z);
  fmpq_poly_clear(sum_poly);
  fmpq_poly_clear(result_rational);
}

void compute_powerset_expr(Context *ctx, Expr *child, Expr *expr_full,
                           fmpz_poly_t res, int n, int is_labeled) {
  fmpz_poly_t A;
  fmpz_poly_init(A);
  compute_e(ctx, child, A, n, is_labeled);

  if (is_labeled) {
    fmpq_poly_t Aq;
    fmpq_poly_init(Aq);
    counts_to_egf(Aq, A, n);

    fmpq_poly_t Rq;
    fmpq_poly_init(Rq);

    if (expr_full->restriction == NONE) {
      fmpq_poly_exp_series(Rq, Aq, n + 1);
    } else {
      long long k_limit = expr_full->limit;
      if (k_limit < 0)
        k_limit = 0;

      fmpq_poly_t sum;
      fmpq_poly_init(sum);
      fmpq_poly_zero(sum);
      fmpq_poly_t term;
      fmpq_poly_init(term);
      fmpq_poly_set_coeff_ui(term, 0, 1);
      fmpq_t fact;
      fmpq_init(fact);
      fmpq_set_si(fact, 1, 1);
      fmpz_t f_z;
      fmpz_init(f_z);

      int max_loop = k_limit;
      if (expr_full->restriction == GREATER)
        max_loop = k_limit - 1;

      for (int j = 0; j <= max_loop; j++) {
        if (expr_full->restriction == LESS ||
            expr_full->restriction == GREATER) {
          fmpq_poly_t scaled_term;
          fmpq_poly_init(scaled_term);
          fmpq_poly_scalar_mul_fmpq(scaled_term, term, fact);
          fmpq_poly_add(sum, sum, scaled_term);
          fmpq_poly_clear(scaled_term);
        } else if (expr_full->restriction == EQUAL && j == k_limit) {
          fmpq_poly_scalar_mul_fmpq(sum, term, fact);
        }

        if (j < max_loop) {
          fmpq_poly_mullow(term, term, Aq, n + 1);
          fmpz_set_si(f_z, j + 1);
          fmpq_div_fmpz(fact, fact, f_z);
        }
      }

      if (expr_full->restriction == LESS || expr_full->restriction == EQUAL) {
        fmpq_poly_set(Rq, sum);
      } else if (expr_full->restriction == GREATER) {
        fmpq_poly_t exp_A;
        fmpq_poly_init(exp_A);
        fmpq_poly_exp_series(exp_A, Aq, n + 1);
        fmpq_poly_sub(Rq, exp_A, sum);
        fmpq_poly_clear(exp_A);
      }

      fmpq_poly_clear(sum);
      fmpq_poly_clear(term);
      fmpq_clear(fact);
      fmpz_clear(f_z);
    }

    egf_to_counts(res, Rq, n);
    fmpq_poly_clear(Aq);
    fmpq_poly_clear(Rq);
  } else {
    if (expr_full->restriction != NONE) {
      compute_multiset_restricted(res, A, n, expr_full->restriction,
                                  expr_full->limit, 1);
    } else {
      compute_powerset_unlabeled(res, A, n);
    }
  }

  fmpz_poly_clear(A);
}
