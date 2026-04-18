#include "solver/count/mset.h"


#include "solver/count/pset.h"
#include "solver/math.h"

void compute_multiset_unlabeled(fmpz_poly_t res, fmpz_poly_t A, int n) {
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
    fmpq_set_si(scale, 1, k);

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

void compute_multiset_restricted(fmpz_poly_t res, fmpz_poly_t A, int n,
                                 int restriction_type, long long limit,
                                 int is_powerset) {
  if (limit < 0) {
    fmpz_poly_zero(res);
    return;
  }
  if (limit > n && restriction_type != GREATER)
    limit = n;

  int compute_limit = limit;
  if (restriction_type == GREATER)
    compute_limit = limit - 1;

  SeriesPoly H;
  sp_init(&H, compute_limit);

  fmpz_poly_t A_z;
  fmpz_poly_init(A_z);
  fmpz_poly_set(A_z, A);
  fmpz_poly_set_coeff_ui(A_z, 0, 0);

  fmpq_t scale;
  fmpq_init(scale);

  for (int k = 1; k <= compute_limit; k++) {
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

    if (is_powerset && (k - 1) % 2 != 0)
      fmpq_set_si(scale, -1, k);
    else
      fmpq_set_si(scale, 1, k);

    fmpq_poly_scalar_mul_fmpq(H.coeffs[k], term, scale);
    fmpq_poly_clear(term);
  }

  SeriesPoly G;
  sp_init(&G, compute_limit);
  sp_exp(&G, &H, n);

  fmpq_poly_t result_rational;
  fmpq_poly_init(result_rational);
  fmpq_poly_zero(result_rational);

  if (restriction_type == 2) {
    if (limit <= compute_limit)
      fmpq_poly_set(result_rational, G.coeffs[limit]);
  } else if (restriction_type == 1) {
    for (int k = 0; k <= limit; k++) {
      fmpq_poly_add(result_rational, result_rational, G.coeffs[k]);
    }
  } else if (restriction_type == 3) {
    fmpq_poly_t less_part;
    fmpq_poly_init(less_part);
    for (int k = 0; k <= compute_limit; k++) {
      fmpq_poly_add(less_part, less_part, G.coeffs[k]);
    }

    fmpz_poly_t unrest;
    fmpz_poly_init(unrest);
    if (is_powerset)
      compute_powerset_unlabeled(unrest, A, n);
    else
      compute_multiset_unlabeled(unrest, A, n);

    fmpq_poly_t unrest_q;
    fmpq_poly_init(unrest_q);
    fmpq_poly_set_fmpz_poly(unrest_q, unrest);

    fmpq_poly_sub(result_rational, unrest_q, less_part);

    fmpq_poly_clear(unrest_q);
    fmpz_poly_clear(unrest);
    fmpq_poly_clear(less_part);
  }

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

  fmpq_clear(scale);
  fmpq_poly_clear(result_rational);
  fmpq_clear(q_coeff);
  fmpz_clear(z_coeff);
  fmpz_poly_clear(A_z);
  sp_clear(&H);
  sp_clear(&G);
}
