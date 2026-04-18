#include "solver/unrank/eval.h"

#include <ctype.h>
#include <flint/fmpz.h>

#include "solver/unrank/cycle.h"
#include "solver/unrank/prod.h"
#include "solver/unrank/pset.h"
#include "solver/unrank/seq.h"
#include "solver/unrank/set.h"
#include "solver/unrank/union.h"

static char *my_strdup(const char *s) {
  if (!s)
    return NULL;
  char *d = malloc(strlen(s) + 1);
  if (d)
    strcpy(d, s);
  return d;
}

char *unrank_e(Context *ctx, Expr *expr, int n, fmpz_t rank, int *labels) {
  switch (expr->type) {
  case ATOM:
  case Z: {
    if (n != 1)
      return my_strdup("ErrorAtomSize");
    if (!ctx->is_labeled)
      return my_strdup("Z()");
    // Return "Z(label)".
    char buf[32];
    sprintf(buf, "Z(%d)", labels[0]);
    return my_strdup(buf);
  }
  case EPSILON: {
    if (n != 0)
      return my_strdup("ErrorEpsilonSize");
    return my_strdup("Eps()");
  }
  case ID: {
    Id *id = (Id *)expr->component;
    for (int i = 0; i < ctx->num_entries; i++) {
      if (strcmp(ctx->entries[i].name, id->name) == 0) {
        return unrank_e(ctx, ctx->entries[i].expr, n, rank, labels);
      }
    }
    if (strcmp(id->name, "Z") == 0) {
      if (n != 1)
        return my_strdup("ErrorZ");
      if (!ctx->is_labeled)
        return my_strdup("Z()");
      char buf[32];
      sprintf(buf, "Z(%d)", labels[0]);
      return my_strdup(buf);
    }
    return my_strdup("ErrorIDNotFound");
  }
  case UNION:
    if (ctx->is_labeled)
      return unrank_union_expr(ctx, expr, n, rank, labels);
    else
      return unrank_union_expr(ctx, expr, n, rank, NULL);
  case PROD:
    if (ctx->is_labeled)
      return unrank_prod_expr(ctx, expr, n, rank, labels);
    else
      return unrank_prod_unlabeled(ctx, expr, n, rank);
  case SET:
    if (ctx->is_labeled)
      return unrank_set_expr(ctx, expr, n, rank, labels);
    else
      return unrank_set_unlabeled(ctx, expr, n, rank);
  case SEQUENCE:
    if (ctx->is_labeled)
      return unrank_seq_expr(ctx, expr, n, rank, labels);
    else
      return unrank_seq_unlabeled(ctx, expr, n, rank);
  case CYCLE:
    if (ctx->is_labeled)
      return unrank_cycle_expr(ctx, expr, n, rank, labels);
    else
      return my_strdup("ErrorCycleUnlabeledNotSupported");
  case POWERSET:
    if (ctx->is_labeled)
      return unrank_powerset_expr(ctx, expr, n, rank, labels);
    else
      return unrank_pset_unlabeled(ctx, expr, n, rank);
  default:
    return my_strdup("ErrorUnknownType");
  }
}

// Wrapper
char *unrank(Context *ctx, char *symbol, int n, fmpz_t rank) {
  // 1. Find root expression
  Expr *root_expr = NULL;
  for (int i = 0; i < ctx->num_entries; i++) {
    if (strcmp(ctx->entries[i].name, symbol) == 0) {
      root_expr = ctx->entries[i].expr;
      break;
    }
  }
  if (!root_expr)
    return strdup("SymbolNotFound");

  if (ctx->is_labeled) {
    // 2. Prepare labels {1..n}
    int *labels = malloc(sizeof(int) * (n > 0 ? n : 1));
    for (int i = 0; i < n; i++)
      labels[i] = i + 1;

    char *res = unrank_e(ctx, root_expr, n, rank, labels);
    free(labels);
    return res;
  } else {
    return unrank_e(ctx, root_expr, n, rank, NULL);
  }
}