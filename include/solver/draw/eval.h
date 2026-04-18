#ifndef DRAW_H
#define DRAW_H

#include <flint/fmpz.h>

#include "grammar/utils.h"

char *draw(Context *ctx, char *symbol, int n, flint_rand_t state);

char *draw_e(Context *ctx, Expr *expr, int n, flint_rand_t state, int *labels);

#endif
