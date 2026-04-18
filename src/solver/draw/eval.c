#include "solver/draw/eval.h"

#include <ctype.h>
#include <flint/fmpz.h>

#include "solver/math.h"
#include "solver/draw/cycle.h"
#include "solver/draw/prod.h"
#include "solver/draw/pset.h"
#include "solver/draw/seq.h"
#include "solver/draw/set.h"
#include "solver/draw/union.h"

char *draw_e(Context *ctx, Expr *expr, int n, flint_rand_t state,
             int *labels) {
  switch (expr->type) {
  case ATOM:
  case Z: {
    if (n != 1)
      return xstrdup("ErrorAtomSize");
    char buf[32];
    sprintf(buf, "Z(%d)", labels[0]);
    return xstrdup(buf);
  }
  case EPSILON: {
    if (n != 0)
      return xstrdup("ErrorEpsilonSize");
    return xstrdup("Eps()");
  }
  case ID: {
    Id *id = (Id *)expr->component;
    for (int i = 0; i < ctx->num_entries; i++) {
      if (strcmp(ctx->entries[i].name, id->name) == 0) {
        return draw_e(ctx, ctx->entries[i].expr, n, state, labels);
      }
    }
    if (strcmp(id->name, "Z") == 0) {
      if (n != 1)
        return xstrdup("ErrorZ");
      char buf[32];
      sprintf(buf, "Z(%d)", labels[0]);
      return xstrdup(buf);
    }
    return xstrdup("ErrorIDNotFound");
  }
  case UNION:
    return draw_union_expr(ctx, expr, n, state, labels);
  case PROD:
    return draw_prod_expr(ctx, expr, n, state, labels);
  case SET:
    return draw_set_expr(ctx, expr, n, state, labels);
  case SEQUENCE:
    return draw_seq_expr(ctx, expr, n, state, labels);
  case CYCLE:
    return draw_cycle_expr(ctx, expr, n, state, labels);
  case POWERSET:
    return draw_powerset_expr(ctx, expr, n, state, labels);
  default:
    return xstrdup("ErrorUnknownType");
  }
}

// Wrapper
char *draw(Context *ctx, char *symbol, int n, flint_rand_t state) {
  // 1. Find root expression
  Expr *root_expr = NULL;
  for (int i = 0; i < ctx->num_entries; i++) {
    if (strcmp(ctx->entries[i].name, symbol) == 0) {
      root_expr = ctx->entries[i].expr;
      break;
    }
  }
  if (!root_expr)
    return xstrdup("SymbolNotFound");

  // 2. Prepare labels {1..n}
  int *labels = n > 0 ? xmalloc(sizeof(int) * n) : NULL;
  for (int i = 0; i < n; i++)
    labels[i] = i + 1;

  char *res = draw_e(ctx, root_expr, n, state, labels);
  free(labels);
  return res;
}
