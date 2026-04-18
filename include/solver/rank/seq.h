#ifndef RANK_SEQ_H
#define RANK_SEQ_H

#include "solver/rank/eval.h"

void rank_seq_expr(Context *ctx, Expr *expr, Object *obj, fmpz_t res);
void rank_seq_unlabeled(Context *ctx, Expr *expr, Object *obj, fmpz_t res);

#endif
