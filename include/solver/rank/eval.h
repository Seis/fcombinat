#ifndef RANK_H
#define RANK_H

#include <flint/fmpz.h>

#include "grammar/object.h"
#include "grammar/utils.h"

void rank(Context *ctx, char *symbol, Object *obj, fmpz_t res);

void rank_e(Context *ctx, Expr *expr, Object *obj, fmpz_t res);

int cmp_int(const void *a, const void *b);

#endif
