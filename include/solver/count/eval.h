#ifndef SOLVER_COUNT_EVAL_H
#define SOLVER_COUNT_EVAL_H

#include <flint/fmpz_poly.h>

#include "grammar/utils.h"

void compute_e(Context *ctx, Expr *expr, fmpz_poly_t res, int n,
               int is_labeled);

#endif
