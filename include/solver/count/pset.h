#ifndef SOLVER_COUNT_PSET_H
#define SOLVER_COUNT_PSET_H

#include <flint/fmpz_poly.h>

#include "grammar/utils.h"

void compute_powerset_unlabeled(fmpz_poly_t res, fmpz_poly_t A, int n);
void compute_powerset_expr(Context *ctx, Expr *child, Expr *expr_full,
                           fmpz_poly_t res, int n, int is_labeled);

#endif
