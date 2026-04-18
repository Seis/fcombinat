#ifndef UNRANK_SET_H
#define UNRANK_SET_H

#include "solver/unrank/eval.h"

char *unrank_set_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                      int *labels);

char *unrank_set_unlabeled(Context *ctx, Expr *expr, int n, fmpz_t rank);

#endif
