#ifndef RANK_PSET_H
#define RANK_PSET_H

#include "solver/rank/eval.h"

void rank_powerset_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res);
void rank_pset_unlabeled(Context *ctx, Expr *expr, Object *obj, fmpz_t res);

#endif
