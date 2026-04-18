#ifndef DRAW_PROD_H
#define DRAW_PROD_H

#include "solver/draw/eval.h"

char *draw_prod_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                     int *labels);

#endif
