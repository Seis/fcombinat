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

char *unrank_e(Context *ctx, Expr *expr, int n, fmpz_t rank, int *labels,
               int depth) {
  /* Boustrophedonic inversion: at odd depths, map rank -> (total-1-rank). */
  fmpz_t effective_rank;
  fmpz_init_set(effective_rank, rank);

  if (ctx->order == ORDER_BOUSTROPHEDON && (depth % 2) == 1) {
    fmpz_t total;
    fmpz_init(total);
    get_expr_count(total, ctx, expr, n);
    fmpz_sub_ui(total, total, 1);
    fmpz_sub(effective_rank, total, effective_rank);
    fmpz_clear(total);
  }

  char *result;
  switch (expr->type) {
  case ATOM:
  case Z: {
    if (n != 1) {
      result = my_strdup("ErrorAtomSize");
      break;
    }
    if (!ctx->is_labeled) {
      result = my_strdup("Z()");
      break;
    }
    char buf[32];
    sprintf(buf, "Z(%d)", labels[0]);
    result = my_strdup(buf);
    break;
  }
  case EPSILON: {
    if (n != 0) {
      result = my_strdup("ErrorEpsilonSize");
      break;
    }
    result = my_strdup("Eps()");
    break;
  }
  case ID: {
    Id *id = (Id *)expr->component;
    for (int i = 0; i < ctx->num_entries; i++) {
      if (strcmp(ctx->entries[i].name, id->name) == 0) {
        fmpz_clear(effective_rank);
        return unrank_e(ctx, ctx->entries[i].expr, n, rank, labels, depth);
      }
    }
    if (strcmp(id->name, "Z") == 0) {
      if (n != 1) {
        result = my_strdup("ErrorZ");
        break;
      }
      if (!ctx->is_labeled) {
        result = my_strdup("Z()");
        break;
      }
      char buf[32];
      sprintf(buf, "Z(%d)", labels[0]);
      result = my_strdup(buf);
      break;
    }
    result = my_strdup("ErrorIDNotFound");
    break;
  }
  case UNION:
    if (ctx->is_labeled)
      result = unrank_union_expr(ctx, expr, n, effective_rank, labels, depth + 1);
    else
      result = unrank_union_expr(ctx, expr, n, effective_rank, NULL, depth + 1);
    break;
  case PROD:
    if (ctx->is_labeled)
      result = unrank_prod_expr(ctx, expr, n, effective_rank, labels, depth + 1);
    else
      result = unrank_prod_unlabeled(ctx, expr, n, effective_rank, depth + 1);
    break;
  case SET:
    if (ctx->is_labeled)
      result = unrank_set_expr(ctx, expr, n, effective_rank, labels, depth + 1);
    else
      result = unrank_set_unlabeled(ctx, expr, n, effective_rank, depth + 1);
    break;
  case SEQUENCE:
    if (ctx->is_labeled)
      result = unrank_seq_expr(ctx, expr, n, effective_rank, labels, depth + 1);
    else
      result = unrank_seq_unlabeled(ctx, expr, n, effective_rank, depth + 1);
    break;
  case CYCLE:
    if (ctx->is_labeled)
      result = unrank_cycle_expr(ctx, expr, n, effective_rank, labels, depth + 1);
    else {
      fprintf(stderr, "Error: unlabeled Cycle unranking is not supported\n");
      result = NULL;
    }
    break;
  case POWERSET:
    if (ctx->is_labeled)
      result = unrank_powerset_expr(ctx, expr, n, effective_rank, labels, depth + 1);
    else
      result = unrank_pset_unlabeled(ctx, expr, n, effective_rank, depth + 1);
    break;
  default:
    result = my_strdup("ErrorUnknownType");
    break;
  }

  fmpz_clear(effective_rank);
  return result;
}

// Wrapper
char *unrank(Context *ctx, char *symbol, int n, fmpz_t rank) {
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
    int *labels = malloc(sizeof(int) * (n > 0 ? n : 1));
    for (int i = 0; i < n; i++)
      labels[i] = i + 1;

    char *res = unrank_e(ctx, root_expr, n, rank, labels, 0);
    free(labels);
    return res;
  } else {
    return unrank_e(ctx, root_expr, n, rank, NULL, 0);
  }
}