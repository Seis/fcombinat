#include "solver/count/cycle.h"

#include <flint/ulong_extras.h>

#include "solver/count/eval.h"
#include "solver/math.h"

void compute_cycle_unlabeled(fmpz_poly_t res, fmpz_poly_t A, int n) {
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

  fmpq_poly_t term;
  fmpq_poly_init(term);
  fmpq_poly_t inner_log;
  fmpq_poly_init(inner_log);
  fmpq_t scale;
  fmpq_init(scale);
  fmpz_t phi;
  fmpz_init(phi);

  for (int k = 1; k <= n; k++) {
    fmpq_poly_zero(term);

    for (int i = 1; i <= n / k; i++) {
      fmpz_t coeff_val;
      fmpz_init(coeff_val);
      fmpz_poly_get_coeff_fmpz(coeff_val, A_z, i);
      if (!fmpz_is_zero(coeff_val)) {
        fmpq_poly_set_coeff_fmpz(term, i * k, coeff_val);
      }
      fmpz_clear(coeff_val);
    }

    fmpq_poly_neg(term, term);
    fmpq_poly_set_coeff_ui(term, 0, 1);

    fmpq_poly_zero(inner_log);
    fmpq_poly_log_series(inner_log, term, n + 1);

    ulong phi_val = n_euler_phi(k);
    fmpz_set_ui(phi, phi_val);

    fmpz_t f_k;
    fmpz_init(f_k);
    fmpz_set_ui(f_k, k);
    fmpq_set_fmpz_frac(scale, phi, f_k);
    fmpz_clear(f_k);

    fmpq_neg(scale, scale);

    fmpq_poly_scalar_mul_fmpq(inner_log, inner_log, scale);
    fmpq_poly_add(sum_poly, sum_poly, inner_log);
  }

  fmpz_poly_zero(res);
  fmpq_t q_coeff;
  fmpq_init(q_coeff);
  fmpz_t z_coeff;
  fmpz_init(z_coeff);

  for (int i = 0; i <= n; i++) {
    fmpq_poly_get_coeff_fmpq(q_coeff, sum_poly, i);
    fmpz_set(z_coeff, fmpq_numref(q_coeff));
    fmpz_poly_set_coeff_fmpz(res, i, z_coeff);
  }

  fmpq_clear(q_coeff);
  fmpz_clear(z_coeff);
  fmpz_clear(a0);
  fmpz_clear(phi);
  fmpq_clear(scale);
  fmpq_poly_clear(term);
  fmpq_poly_clear(inner_log);
  fmpq_poly_clear(sum_poly);
  fmpz_poly_clear(A_z);
}

void compute_cycle_restricted_unlabeled(fmpz_poly_t res, fmpz_poly_t A, int n,
                                        int k) {
  fmpz_poly_zero(res);
  if (k == 0) {
    fmpz_poly_set_coeff_ui(res, 0, 1);
    return;
  }

  fmpz_poly_t sum;
  fmpz_poly_init(sum);
  fmpz_poly_zero(sum);
  fmpz_poly_t term;
  fmpz_poly_init(term);
  fmpz_poly_t A_pow;
  fmpz_poly_init(A_pow);

  fmpz_t phi;
  fmpz_init(phi);
  fmpz_t d_z;
  fmpz_init(d_z);

  for (int d = 1; d <= k; d++) {
    if (k % d == 0) {
      fmpz_poly_zero(term);
      for (int i = 1; i <= n / d; i++) {
        fmpz_t coeff_val;
        fmpz_init(coeff_val);
        fmpz_poly_get_coeff_fmpz(coeff_val, A, i);
        if (!fmpz_is_zero(coeff_val)) {
          fmpz_poly_set_coeff_fmpz(term, i * d, coeff_val);
        }
        fmpz_clear(coeff_val);
      }
      fmpz_t a0;
      fmpz_init(a0);
      fmpz_poly_get_coeff_fmpz(a0, A, 0);
      if (!fmpz_is_zero(a0)) {
        fmpz_poly_set_coeff_fmpz(term, 0, a0);
      }
      fmpz_clear(a0);

      fmpz_poly_pow_trunc(A_pow, term, k / d, n + 1);

      ulong phi_val = n_euler_phi(d);
      fmpz_set_ui(phi, phi_val);
      fmpz_poly_scalar_mul_fmpz(A_pow, A_pow, phi);

      fmpz_poly_add(sum, sum, A_pow);
    }
  }

  fmpz_set_ui(d_z, k);
  for (int i = 0; i <= n; i++) {
    fmpz_t c;
    fmpz_init(c);
    fmpz_poly_get_coeff_fmpz(c, sum, i);
    fmpz_divexact(c, c, d_z);
    fmpz_poly_set_coeff_fmpz(res, i, c);
    fmpz_clear(c);
  }

  fmpz_clear(phi);
  fmpz_clear(d_z);
  fmpz_poly_clear(sum);
  fmpz_poly_clear(term);
  fmpz_poly_clear(A_pow);
}

void compute_cycle_expr(Context *ctx, Expr *child, Expr *expr_full,
                        fmpz_poly_t res, int n, int is_labeled) {
  fmpz_poly_t A;
  fmpz_poly_init(A);
  compute_e(ctx, child, A, n, is_labeled);

  if (is_labeled) {
    fmpq_poly_t Aq;
    fmpq_poly_init(Aq);
    counts_to_egf(Aq, A, n);
    fmpq_poly_set_coeff_ui(Aq, 0, 0);

    fmpq_poly_t Rq;
    fmpq_poly_init(Rq);

    if (expr_full->restriction == NONE) {
      fmpq_poly_t one;
      fmpq_poly_init(one);
      fmpq_poly_set_coeff_ui(one, 0, 1);
      fmpq_poly_sub(one, one, Aq);
      fmpq_poly_log_series(Rq, one, n + 1);
      fmpq_poly_scalar_mul_si(Rq, Rq, -1);
      fmpq_poly_clear(one);
    } else {
      long long limit = expr_full->limit;
      if (limit < 0)
        limit = 0;

      fmpq_poly_zero(Rq);
      fmpq_poly_t term;
      fmpq_poly_init(term);
      fmpq_poly_t sum;
      fmpq_poly_init(sum);
      fmpq_poly_zero(sum);

      int max_idx = limit;
      if (expr_full->restriction == GREATER)
        max_idx = limit - 1;

      fmpq_poly_set(term, Aq);

      for (int k = 1; k <= max_idx; k++) {
        if ((expr_full->restriction == LESS ||
             expr_full->restriction == GREATER) ||
            (expr_full->restriction == EQUAL && k == limit)) {
          fmpq_poly_t tmp;
          fmpq_poly_init(tmp);
          fmpq_t scale;
          fmpq_init(scale);
          fmpq_set_si(scale, 1, k);
          fmpq_poly_scalar_mul_fmpq(tmp, term, scale);
          fmpq_poly_add(sum, sum, tmp);
          fmpq_poly_clear(tmp);
          fmpq_clear(scale);
        }

        if (k < max_idx) {
          fmpq_poly_mullow(term, term, Aq, n + 1);
        }
      }
      fmpq_poly_clear(term);

      if (expr_full->restriction == EQUAL || expr_full->restriction == LESS) {
        fmpq_poly_set(Rq, sum);
      } else if (expr_full->restriction == GREATER) {
        fmpq_poly_t one;
        fmpq_poly_init(one);
        fmpq_poly_set_coeff_ui(one, 0, 1);
        fmpq_poly_sub(one, one, Aq);
        fmpq_poly_log_series(Rq, one, n + 1);
        fmpq_poly_scalar_mul_si(Rq, Rq, -1);
        fmpq_poly_sub(Rq, Rq, sum);
        fmpq_poly_clear(one);
      }
      fmpq_poly_clear(sum);
    }

    egf_to_counts(res, Rq, n);
    fmpq_poly_clear(Aq);
    fmpq_poly_clear(Rq);
  } else {
    if (expr_full->restriction == NONE) {
      compute_cycle_unlabeled(res, A, n);
    } else if (expr_full->restriction == EQUAL) {
      compute_cycle_restricted_unlabeled(res, A, n, (int)expr_full->limit);
    } else {
      fmpz_poly_zero(res);
      fmpz_poly_t tmp;
      fmpz_poly_init(tmp);

      int max_idx = expr_full->limit;
      if (expr_full->restriction == GREATER) {
        compute_cycle_unlabeled(res, A, n);
        max_idx = expr_full->limit - 1;
        for (int k = 1; k <= max_idx; k++) {
          compute_cycle_restricted_unlabeled(tmp, A, n, k);
          fmpz_poly_sub(res, res, tmp);
        }
      } else {
        for (int k = 1; k <= max_idx; k++) {
          compute_cycle_restricted_unlabeled(tmp, A, n, k);
          fmpz_poly_add(res, res, tmp);
        }
      }
      fmpz_poly_clear(tmp);
    }
  }

  fmpz_poly_clear(A);
}
