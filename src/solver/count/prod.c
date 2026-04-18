#include "solver/count/prod.h"

#include "solver/count/eval.h"
#include "solver/math.h"

void compute_prod_expr(Context *ctx, ExprList *el, fmpz_poly_t res, int n,
                       int is_labeled) {
  fmpz_poly_zero(res);
  fmpz_poly_set_coeff_ui(res, 0, 1);

  fmpz_poly_t tmp;
  fmpz_poly_init(tmp);

  for (int i = 0; i < el->size; i++) {
    compute_e(ctx, el->components[i], tmp, n, is_labeled);
    if (is_labeled) {
      fmpz_poly_t partial;
      fmpz_poly_init(partial);
      binomial_convolution(partial, res, tmp, n);
      fmpz_poly_set(res, partial);
      fmpz_poly_clear(partial);
    } else {
      fmpz_poly_mullow(res, res, tmp, n + 1);
    }
  }
  fmpz_poly_clear(tmp);
}
