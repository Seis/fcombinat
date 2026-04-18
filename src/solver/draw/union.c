#include "solver/draw/union.h"

#include "solver/math.h"

char *draw_union_expr(Context *ctx, Expr *expr, int n, flint_rand_t state,
                      int *labels) {
  ExprList *el = (ExprList *)expr->component;
  fmpz_t count;
  fmpz_init(count);

  // Compute total count across all branches
  fmpz_t total;
  fmpz_init(total);
  for (int i = 0; i < el->size; i++) {
    get_expr_count(count, ctx, el->components[i], n);
    fmpz_add(total, total, count);
  }

  if (fmpz_is_zero(total)) {
    fmpz_clear(count);
    fmpz_clear(total);
    return xstrdup("ErrorUnionEmpty");
  }

  // Draw uniform random in [0, total)
  fmpz_t current_rank;
  fmpz_init(current_rank);
  fmpz_randm(current_rank, state, total);

  // Same subtract-and-compare loop as unrank
  for (int i = 0; i < el->size; i++) {
    get_expr_count(count, ctx, el->components[i], n);
    if (fmpz_cmp(current_rank, count) < 0) {
      char *child_res =
          draw_e(ctx, el->components[i], n, state, labels);

      char *tagged_res = xmalloc(strlen(child_res) + 30);
      sprintf(tagged_res, "Union_%d(%s)", i, child_res);
      free(child_res);

      fmpz_clear(count);
      fmpz_clear(total);
      fmpz_clear(current_rank);
      return tagged_res;
    }
    fmpz_sub(current_rank, current_rank, count);
  }
  fmpz_clear(count);
  fmpz_clear(total);
  fmpz_clear(current_rank);
  return xstrdup("ErrorUnionBounds");
}
