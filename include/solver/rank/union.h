#ifndef RANK_UNION_H
#define RANK_UNION_H

#include "solver/rank/eval.h"

void rank_union_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res);

#endif
