#include "grammar/utils.h"


#include "solver/count/eval.h"

// --- Context Management ---

void context_init(Context *ctx) {
  ctx->entries = NULL;
  ctx->num_entries = 0;
}

void context_clear(Context *ctx) {
  if (ctx->entries) {
    for (int i = 0; i < ctx->num_entries; i++) {
      fmpz_poly_clear(ctx->entries[i].poly);
    }
    free(ctx->entries);
  }
}

fmpz_poly_struct *get_poly(Context *ctx, char *name) {
  for (int i = 0; i < ctx->num_entries; i++) {
    if (strcmp(ctx->entries[i].name, name) == 0) {
      return ctx->entries[i].poly;
    }
  }
  return NULL;
}

void get_expr_count(fmpz_t res, Context *ctx, Expr *expr, int n) {
  fmpz_poly_t tmp;
  fmpz_poly_init(tmp);
  compute_e(ctx, expr, tmp, n, ctx->is_labeled);
  fmpz_poly_get_coeff_fmpz(res, tmp, n);
  fmpz_poly_clear(tmp);
}

Context *solve_spec(Spec *g, int max_n, int is_labeled) {
  if (g->type == ISERROR)
    return NULL;

  Rules *sl = (Rules *)g->component;

  Context *ctx = malloc(sizeof(Context));
  ctx->num_entries = sl->size;
  ctx->entries = malloc(sizeof(Entry) * ctx->num_entries);
  ctx->is_labeled = is_labeled;

  for (int i = 0; i < ctx->num_entries; i++) {
    ctx->entries[i].name = sl->components[i]->variable->name;
    ctx->entries[i].expr = sl->components[i]->expression;
    fmpz_poly_init(ctx->entries[i].poly);
    fmpz_poly_zero(ctx->entries[i].poly);
  }

  int limit_iter = max_n * (ctx->num_entries + 1);

  for (int iter = 0; iter <= limit_iter; iter++) {
    for (int i = 0; i < ctx->num_entries; i++) {
      Rule *s = sl->components[i];
      if (strcmp(s->variable->name, "Z") == 0)
        continue;

      fmpz_poly_t tmp;
      fmpz_poly_init(tmp);

      compute_e(ctx, s->expression, tmp, max_n, is_labeled);

      fmpz_poly_set(ctx->entries[i].poly, tmp);
      fmpz_poly_clear(tmp);
    }
  }

  return ctx;
}
