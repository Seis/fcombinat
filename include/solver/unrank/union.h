#ifndef UNRANK_UNION_H
#define UNRANK_UNION_H

#include "solver/unrank/eval.h"

char *unrank_union_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                        int *labels, int depth);

#endif
