#include "solver/unrank/union.h"


char *unrank_union_expr(Context *ctx, Expr *expr, int n, fmpz_t rank,
                        int *labels, int depth) {
  ExprList *el = (ExprList *)expr->component;
  fmpz_t count;
  fmpz_init(count);
  fmpz_t current_rank;
  fmpz_init_set(current_rank, rank);

  for (int i = 0; i < el->size; i++) {
    get_expr_count(count, ctx, el->components[i], n);
    if (fmpz_cmp(current_rank, count) < 0) {
      char *child_res =
          unrank_e(ctx, el->components[i], n, current_rank, labels, depth);

      char *tagged_res = malloc(strlen(child_res) + 30);
      sprintf(tagged_res, "Union_%d(%s)", i, child_res);
      free(child_res);

      fmpz_clear(count);
      fmpz_clear(current_rank);
      return tagged_res;
    }
    fmpz_sub(current_rank, current_rank, count);
  }
  fmpz_clear(count);
  fmpz_clear(current_rank);
  return strdup("ErrorUnionBounds");
}
