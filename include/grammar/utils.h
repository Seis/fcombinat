#ifndef GRAMMAR_UTILS_H
#define GRAMMAR_UTILS_H

#include <flint/fmpz_poly.h>

#include "grammar/absyn.h"

typedef struct {
  char *name;
  Expr *expr;
  fmpz_poly_t poly;
} Entry;

typedef enum {
  ORDER_LEX           = 0,
  ORDER_BOUSTROPHEDON = 1
} order_t;

typedef struct {
  Entry   *entries;
  int      num_entries;
  int      is_labeled;
  order_t  order;
} Context;

void context_init(Context *ctx);
void context_clear(Context *ctx);
fmpz_poly_struct *get_poly(Context *ctx, char *name);

Context *solve_spec(Spec *g, int max_n, int is_labeled);

void get_expr_count(fmpz_t res, Context *ctx, Expr *expr, int n);

#endif