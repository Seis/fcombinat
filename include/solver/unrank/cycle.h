#ifndef UNRANK_CYCLE_H
#define UNRANK_CYCLE_H

#include "solver/unrank/eval.h"

char *unrank_cycle_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                        int *labels);

#endif
