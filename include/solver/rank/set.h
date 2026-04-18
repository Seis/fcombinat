#ifndef RANK_SET_H
#define RANK_SET_H

#include "solver/rank/eval.h"

void rank_set_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res);
void rank_set_unlabeled(Context *ctx, Expr *expr, Object *obj, fmpz_t res);

#endif
