#ifndef SOLVER_COUNT_UNION_H
#define SOLVER_COUNT_UNION_H

#include <flint/fmpz_poly.h>

#include "grammar/utils.h"

void compute_union_expr(Context *ctx, ExprList *el, fmpz_poly_t res, int n,
                        int is_labeled);

#endif
