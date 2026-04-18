#ifndef BOLTZMANN_EVAL_H
#define BOLTZMANN_EVAL_H

#include <flint/fmpz.h>

#include "grammar/utils.h"
#include "solver/boltzmann/oracle.h"

// Generate a random object from the Boltzmann distribution.
// Returns NULL if generation fails.
// actual_size is set to the size of the generated object.
char *boltzmann_draw_e(Oracle *orc, Context *ctx, Expr *expr, double x,
                       flint_rand_t state, int *actual_size, int is_labeled);

// Main entry: generate a random object of exact size target_n.
// Uses rejection sampling: generates objects and rejects those not of target_n.
// tolerance: accept objects of size target_n ± tolerance (0 for exact).
// max_attempts: maximum number of rejection attempts.
char *boltzmann_generate(Context *ctx, char *symbol, int target_n,
                         int tolerance, flint_rand_t state, int max_attempts,
                         int *actual_size, int max_n, int is_labeled);

#endif
