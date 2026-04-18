#ifndef SOLVER_COUNT_CYCLE_H
#define SOLVER_COUNT_CYCLE_H

#include <flint/fmpz_poly.h>

#include "grammar/utils.h"

void compute_cycle_unlabeled(fmpz_poly_t res, fmpz_poly_t A, int n);
void compute_cycle_restricted_unlabeled(fmpz_poly_t res, fmpz_poly_t A, int n,
                                        int k);
void compute_cycle_expr(Context *ctx, Expr *child, Expr *expr_full,
                        fmpz_poly_t res, int n, int is_labeled);

#endif
