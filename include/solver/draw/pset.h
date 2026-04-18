#ifndef DRAW_PSET_H
#define DRAW_PSET_H

#include "solver/draw/eval.h"

char *draw_powerset_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                         int *labels);

#endif
