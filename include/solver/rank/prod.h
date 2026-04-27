#ifndef RANK_PROD_H
#define RANK_PROD_H

#include "solver/rank/eval.h"

void rank_prod_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res, int depth);
void rank_prod_unlabeled(Context *ctx, Expr *expr, Object *obj, fmpz_t res, int depth);

#endif
