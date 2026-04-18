#include "solver/count/eval.h"

#include <flint/fmpz.h>

#include "solver/count/cycle.h"
#include "solver/count/prod.h"
#include "solver/count/pset.h"
#include "solver/count/seq.h"
#include "solver/count/set.h"
#include "solver/count/union.h"

void compute_e(Context *ctx, Expr *expr, fmpz_poly_t res, int n,
               int is_labeled) {
  switch (expr->type) {
  case ATOM:
  case Z:
    fmpz_poly_zero(res);
    fmpz_poly_set_coeff_ui(res, 1, 1);
    break;

  case EPSILON:
    fmpz_poly_zero(res);
    fmpz_poly_set_coeff_ui(res, 0, 1);
    break;

  case ID: {
    Id *id = (Id *)expr->component;
    fmpz_poly_struct *cached = get_poly(ctx, id->name);
    if (cached) {
      fmpz_poly_set(res, cached);
    } else {
      fmpz_poly_zero(res);
    }
    break;
  }

  case UNION: {
    ExprList *el = (ExprList *)expr->component;
    compute_union_expr(ctx, el, res, n, is_labeled);
    break;
  }

  case PROD: {
    ExprList *el = (ExprList *)expr->component;
    compute_prod_expr(ctx, el, res, n, is_labeled);
    break;
  }

  case SEQUENCE: {
    Expr *child = (Expr *)expr->component;
    compute_seq_expr(ctx, child, expr, res, n, is_labeled);
    break;
  }

  case SET: {
    Expr *child = (Expr *)expr->component;
    compute_set_expr(ctx, child, expr, res, n, is_labeled);
    break;
  }

  case POWERSET: {
    Expr *child = (Expr *)expr->component;
    compute_powerset_expr(ctx, child, expr, res, n, is_labeled);
    break;
  }
  case CYCLE: {
    Expr *child = (Expr *)expr->component;
    compute_cycle_expr(ctx, child, expr, res, n, is_labeled);
    break;
  }
  default:
    fmpz_poly_zero(res);
    break;
  }
}
