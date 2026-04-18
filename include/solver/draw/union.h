#ifndef DRAW_UNION_H
#define DRAW_UNION_H

#include "solver/draw/eval.h"

char *draw_union_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                      int *labels);

#endif
