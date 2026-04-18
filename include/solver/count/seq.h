#ifndef SOLVER_COUNT_SEQ_H
#define SOLVER_COUNT_SEQ_H

#include <flint/fmpz_poly.h>

#include "grammar/utils.h"

void compute_seq_expr(Context *ctx, Expr *child, Expr *expr_full,
                      fmpz_poly_t res, int n, int is_labeled);

#endif
