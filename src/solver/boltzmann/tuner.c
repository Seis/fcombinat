#include "solver/boltzmann/tuner.h"

#include <math.h>
#include <stdio.h>

#include <flint/fmpz.h>

// Declared in oracle.c
extern double oracle_expected_size(Context *ctx, char *symbol, double x,
                                   int max_n, int is_labeled);

double estimate_radius(Context *ctx, char *symbol, int max_n,
                       int is_labeled) {
  // Use ratio test: radius ≈ lim |a_n / a_{n+1}|
  // For EGF coefficients (a_n/n!), radius = lim |a_n/n!| / |a_{n+1}/(n+1)!|
  //                                       = lim |a_n / a_{n+1}| * (n+1)
  // Use the highest available k for best convergence.
  // Work in log-space to avoid double overflow for large coefficients.
  fmpz_poly_struct *poly = get_poly(ctx, symbol);
  if (!poly)
    return 1.0;

  // Iterate from the highest k downward, return the first valid ratio
  for (int k = max_n - 1; k >= 1; k--) {
    fmpz_t ck, ck1;
    fmpz_init(ck);
    fmpz_init(ck1);
    fmpz_poly_get_coeff_fmpz(ck, poly, k);
    fmpz_poly_get_coeff_fmpz(ck1, poly, k + 1);

    if (!fmpz_is_zero(ck) && !fmpz_is_zero(ck1)) {
      double log_ck = fmpz_dlog(ck);
      double log_ck1 = fmpz_dlog(ck1);
      double log_ratio = log_ck - log_ck1;
      if (is_labeled) {
        log_ratio += log(k + 1);
      }
      double r = exp(log_ratio);
      fmpz_clear(ck);
      fmpz_clear(ck1);
      if (r > 0.0 && r < 1e9)
        return r;
    }

    fmpz_clear(ck);
    fmpz_clear(ck1);
  }

  return 1.0;
}

double boltzmann_tune(Context *ctx, char *symbol, int target_n, int max_n,
                      int is_labeled) {
  double radius = estimate_radius(ctx, symbol, max_n, is_labeled);

  // Binary search for x such that E[size](x) = target_n
  // E[size](x) is monotonically increasing in x on (0, radius)
  double lo = 0.0;
  double hi = radius * 0.999; // Stay below singularity

  for (int iter = 0; iter < 100; iter++) {
    double mid = (lo + hi) / 2.0;
    double e_size = oracle_expected_size(ctx, symbol, mid, max_n, is_labeled);

    if (e_size < (double)target_n) {
      lo = mid;
    } else {
      hi = mid;
    }

    if (fabs(hi - lo) < 1e-15)
      break;
  }

  return (lo + hi) / 2.0;
}
