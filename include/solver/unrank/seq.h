#ifndef UNRANK_SEQ_H
#define UNRANK_SEQ_H

#include "solver/unrank/eval.h"

char *unrank_seq_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                      int *labels);

char *unrank_seq_unlabeled(Context *ctx, Expr *expr, int n, fmpz_t rank);

#endif
