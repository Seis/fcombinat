#ifndef DRAW_SEQ_H
#define DRAW_SEQ_H

#include "solver/draw/eval.h"

char *draw_seq_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                    int *labels);

#endif
