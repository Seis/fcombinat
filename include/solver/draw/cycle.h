#ifndef DRAW_CYCLE_H
#define DRAW_CYCLE_H

#include "solver/draw/eval.h"

char *draw_cycle_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                      int *labels);

#endif
