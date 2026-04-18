#ifndef BOLTZMANN_TUNER_H
#define BOLTZMANN_TUNER_H

#include "grammar/utils.h"
#include "solver/boltzmann/oracle.h"

// Find the parameter x such that the expected size of generated
// objects from the Boltzmann distribution equals target_n.
// Returns x in (0, radius_of_convergence).
double boltzmann_tune(Context *ctx, char *symbol, int target_n, int max_n,
                      int is_labeled);

// Estimate the radius of convergence from polynomial coefficients.
double estimate_radius(Context *ctx, char *symbol, int max_n, int is_labeled);

#endif
