#ifndef UNRANK_PSET_H
#define UNRANK_PSET_H

#include "solver/unrank/eval.h"

char *unrank_powerset_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                           int *labels);

char *unrank_pset_unlabeled(Context *ctx, Expr *expr, int n, fmpz_t rank);

#endif
