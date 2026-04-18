#include "solver/boltzmann/oracle.h"

#include <math.h>
#include <string.h>

#include <flint/fmpz.h>

#include "solver/count/eval.h"
#include "solver/math.h"

// Evaluate a polynomial at point x: sum_k coeff[k] * x^k
// For labeled (EGF): sum_k coeff[k] * x^k / k!
// Uses log-space computation for EGF to avoid double overflow.
static double eval_poly_at(fmpz_poly_t poly, double x, int max_n,
                           int is_labeled) {
  double result = 0.0;
  double x_pow = 1.0;
  double log_x = (x > 0.0) ? log(x) : -1e30;

  for (int k = 0; k <= max_n; k++) {
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_poly_get_coeff_fmpz(coeff, poly, k);

    if (!fmpz_is_zero(coeff)) {
      if (is_labeled) {
        // Compute a_k * x^k / k! in log-space to avoid overflow
        double log_abs_coeff = fmpz_dlog(coeff);
        double log_term = log_abs_coeff + k * log_x - lgamma(k + 1);
        double term = exp(log_term);
        if (fmpz_sgn(coeff) < 0)
          term = -term;
        result += term;
      } else {
        double c = fmpz_get_d(coeff);
        result += c * x_pow;
      }
    }

    fmpz_clear(coeff);
    if (!is_labeled)
      x_pow *= x;
  }

  return result;
}

// Evaluate the derivative A'(x) of a polynomial at point x
// For labeled (EGF): sum_k k * coeff[k] * x^(k-1) / k! = sum_k coeff[k] * x^(k-1) / (k-1)!
// Uses log-space computation for EGF to avoid double overflow.
static double eval_poly_deriv_at(fmpz_poly_t poly, double x, int max_n,
                                 int is_labeled) {
  double result = 0.0;
  double x_pow = 1.0; // x^(k-1) starting at k=1
  double log_x = (x > 0.0) ? log(x) : -1e30;

  for (int k = 1; k <= max_n; k++) {
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_poly_get_coeff_fmpz(coeff, poly, k);

    if (!fmpz_is_zero(coeff)) {
      if (is_labeled) {
        // Compute a_k * x^(k-1) / (k-1)! in log-space
        double log_abs_coeff = fmpz_dlog(coeff);
        double log_term = log_abs_coeff + (k - 1) * log_x - lgamma(k);
        double term = exp(log_term);
        if (fmpz_sgn(coeff) < 0)
          term = -term;
        result += term;
      } else {
        double c = fmpz_get_d(coeff);
        result += (double)k * c * x_pow;
      }
    }

    fmpz_clear(coeff);
    if (!is_labeled)
      x_pow *= x;
  }

  return result;
}

Oracle *oracle_create(Context *ctx, int max_n, double x, int is_labeled) {
  Oracle *orc = malloc(sizeof(Oracle));
  if (!orc)
    return NULL;
  orc->num_entries = ctx->num_entries;
  orc->x = x;
  orc->max_n = max_n;
  orc->is_labeled = is_labeled;
  orc->names = xmalloc(sizeof(char *) * ctx->num_entries);
  orc->values = xmalloc(sizeof(double) * ctx->num_entries);

  for (int i = 0; i < ctx->num_entries; i++) {
    orc->names[i] = xstrdup(ctx->entries[i].name);
    orc->values[i] = eval_poly_at(ctx->entries[i].poly, x, max_n, is_labeled);
  }

  return orc;
}

double oracle_get(Oracle *orc, char *name) {
  for (int i = 0; i < orc->num_entries; i++) {
    if (strcmp(orc->names[i], name) == 0)
      return orc->values[i];
  }
  return 0.0;
}

double oracle_eval_expr(Oracle *orc, Context *ctx, Expr *expr, double x,
                        int max_n, int is_labeled) {
  // Fast path: atoms evaluate to x (both OGF and EGF, since coeff is 1)
  if (expr->type == Z || expr->type == ATOM)
    return x;

  // Fast path: epsilon evaluates to 1
  if (expr->type == EPSILON)
    return 1.0;

  // Fast path: named ID - use pre-cached oracle value
  if (expr->type == ID) {
    Id *id = (Id *)expr->component;
    double val = oracle_get(orc, id->name);
    if (val != 0.0 || strcmp(id->name, "Z") == 0)
      return (strcmp(id->name, "Z") == 0) ? x : val;
  }

  // Slow path: compute count polynomial for this expression, then evaluate
  fmpz_poly_t tmp;
  fmpz_poly_init(tmp);
  compute_e(ctx, expr, tmp, max_n, is_labeled);
  double val = eval_poly_at(tmp, x, max_n, is_labeled);
  fmpz_poly_clear(tmp);
  return val;
}

// Get the expected size E[size] = x * A'(x) / A(x) for entry 'symbol'
double oracle_expected_size(Context *ctx, char *symbol, double x, int max_n,
                            int is_labeled) {
  for (int i = 0; i < ctx->num_entries; i++) {
    if (strcmp(ctx->entries[i].name, symbol) == 0) {
      double ax = eval_poly_at(ctx->entries[i].poly, x, max_n, is_labeled);
      double axp =
          eval_poly_deriv_at(ctx->entries[i].poly, x, max_n, is_labeled);
      if (ax <= 0.0)
        return 0.0;
      return x * axp / ax;
    }
  }
  return 0.0;
}

void oracle_free(Oracle *orc) {
  if (!orc)
    return;
  for (int i = 0; i < orc->num_entries; i++)
    free(orc->names[i]);
  free(orc->names);
  free(orc->values);
  free(orc);
}
