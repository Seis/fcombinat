#include "solver/count/union.h"

#include "solver/count/eval.h"

void compute_union_expr(Context *ctx, ExprList *el, fmpz_poly_t res, int n,
                        int is_labeled) {
  fmpz_poly_zero(res);
  fmpz_poly_t tmp;
  fmpz_poly_init(tmp);

  for (int i = 0; i < el->size; i++) {
    compute_e(ctx, el->components[i], tmp, n, is_labeled);
    fmpz_poly_add(res, res, tmp);
  }
  fmpz_poly_clear(tmp);
}
