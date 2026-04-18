#ifndef DRAW_SET_H
#define DRAW_SET_H

#include "solver/draw/eval.h"

char *draw_set_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                    int *labels);

#endif
