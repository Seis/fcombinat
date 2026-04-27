#ifndef UNRANK_H
#define UNRANK_H

#include <flint/fmpz.h>

#include "grammar/utils.h"

char *unrank(Context *ctx, char *symbol, int n, fmpz_t rank);

char *unrank_e(Context *ctx, Expr *expr, int n, fmpz_t rank, int *labels, int depth);

#endif
