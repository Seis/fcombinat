#include "solver/rank/eval.h"

#include <flint/fmpz.h>
#include <string.h>

#include "solver/math.h"
#include "solver/rank/cycle.h"
#include "solver/rank/prod.h"
#include "solver/rank/pset.h"
#include "solver/rank/seq.h"
#include "solver/rank/set.h"
#include "solver/rank/union.h"
#include "solver/unrank/eval.h"

void rank(Context *ctx, char *symbol, Object *obj, fmpz_t res) {
  Expr *root_expr = NULL;
  for (int i = 0; i < ctx->num_entries; i++) {
    if (strcmp(ctx->entries[i].name, symbol) == 0) {
      root_expr = ctx->entries[i].expr;
      break;
    }
  }
  if (!root_expr) {
    printf("Error: Symbol %s not found\n", symbol);
    fmpz_set_si(res, -1);
    return;
  }
  rank_e(ctx, root_expr, obj, res);
}

void rank_e(Context *ctx, Expr *expr, Object *obj, fmpz_t res) {
  switch (expr->type) {
  case ATOM:
  case Z:
    if (obj->type != OBJ_ATOM) {
      fmpz_set_si(res, -1);
      return;
    }
    fmpz_zero(res);
    return;

  case EPSILON:
    fmpz_zero(res);
    return;

  case ID: {
    Id *id = (Id *)expr->component;
    for (int i = 0; i < ctx->num_entries; i++) {
      if (strcmp(ctx->entries[i].name, id->name) == 0) {
        rank_e(ctx, ctx->entries[i].expr, obj, res);
        return;
      }
    }
    if (strcmp(id->name, "Z") == 0) {
      fmpz_zero(res);
      return;
    }
    fmpz_set_si(res, -1);
    return;
  }

  case UNION:
    rank_union_expr(ctx, expr, obj, res);
    return;

  case PROD:
    if (ctx->is_labeled)
      rank_prod_expr(ctx, expr, obj, res);
    else
      rank_prod_unlabeled(ctx, expr, obj, res);
    return;

  case SET:
    if (ctx->is_labeled)
      rank_set_expr(ctx, expr, obj, res);
    else
      rank_set_unlabeled(ctx, expr, obj, res);
    return;

  case SEQUENCE:
    if (ctx->is_labeled)
      rank_seq_expr(ctx, expr, obj, res);
    else
      rank_seq_unlabeled(ctx, expr, obj, res);
    return;

  case CYCLE:
    if (ctx->is_labeled)
      rank_cycle_expr(ctx, expr, obj, res);
    else
      fmpz_set_si(res, -1); /* Cycle unlabeled not supported */
    return;

  case POWERSET:
    if (ctx->is_labeled)
      rank_powerset_expr(ctx, expr, obj, res);
    else
      rank_pset_unlabeled(ctx, expr, obj, res);
    return;

  default:
    fmpz_set_si(res, -1);
  }
}
