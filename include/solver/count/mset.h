#ifndef SOLVER_COUNT_MSET_H
#define SOLVER_COUNT_MSET_H

#include <flint/fmpz_poly.h>

void compute_multiset_unlabeled(fmpz_poly_t res, fmpz_poly_t A, int n);
void compute_multiset_restricted(fmpz_poly_t res, fmpz_poly_t A, int n,
                                 int restriction_type, long long limit,
                                 int is_powerset);

#endif
