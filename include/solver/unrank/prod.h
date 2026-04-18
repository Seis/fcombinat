#ifndef UNRANK_PROD_H
#define UNRANK_PROD_H

#include "solver/unrank/eval.h"

char *unrank_prod_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                       int *labels);

char *unrank_prod_unlabeled(Context *ctx, Expr *expr, int n, fmpz_t rank);

#endif
