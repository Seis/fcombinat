#ifndef BOLTZMANN_ORACLE_H
#define BOLTZMANN_ORACLE_H

#include "grammar/utils.h"

// Oracle: evaluates generating functions at a real point x
typedef struct {
  double *values;   // values[i] = A_i(x) for each entry
  int num_entries;
  char **names;
  double x;         // The evaluation point
  int max_n;
  int is_labeled;
} Oracle;

Oracle *oracle_create(Context *ctx, int max_n, double x, int is_labeled);
double oracle_get(Oracle *orc, char *name);
double oracle_eval_expr(Oracle *orc, Context *ctx, Expr *expr, double x,
                        int max_n, int is_labeled);
void oracle_free(Oracle *orc);

#endif
