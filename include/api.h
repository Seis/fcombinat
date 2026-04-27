#ifndef FCOMB_API_H
#define FCOMB_API_H

#include <flint/fmpz.h>

#include "grammar/absyn.h"
#include "grammar/utils.h"

struct fcomb_ctx_s {
  Spec     *spec;
  Context  *solve_ctx;
  char     *symbol;   /* owned (strdup'd); NULL means first rule */
  int       n;
  int       labeled;
  ulong     seed;     /* 0 = auto (time ^ pid) */
  order_t   order;
};
typedef struct fcomb_ctx_s *fcomb_ctx_t;

/* Lifecycle */
fcomb_ctx_t fcomb_ctx_init(void);
void        fcomb_ctx_destroy(fcomb_ctx_t ctx);

/* Setters - return 0 on success, non-zero on bad input.
   Any setter that changes spec/n/labeled invalidates solve_ctx. */
int fcomb_ctx_set_spec_str (const char *spec,       fcomb_ctx_t ctx);
int fcomb_ctx_set_spec_json(const char *path,       fcomb_ctx_t ctx);
int fcomb_ctx_set_symbol   (const char *sym,        fcomb_ctx_t ctx);
int fcomb_ctx_set_size     (int n,                  fcomb_ctx_t ctx);
int  fcomb_ctx_set_labeled (int labeled,   fcomb_ctx_t ctx);
void fcomb_ctx_set_seed    (ulong seed,    fcomb_ctx_t ctx);
int  fcomb_ctx_set_order   (order_t order, fcomb_ctx_t ctx);

/* Solve: compute GF coefficients up to n.  Must be called before any
   operation below (except fcomb_boltzmann, which handles its own context). */
int fcomb_solve(fcomb_ctx_t ctx);

/* Operations - require fcomb_solve() to have succeeded first,
   unless noted otherwise. */
void  fcomb_count    (fmpz_t res,                     fcomb_ctx_t ctx);
void  fcomb_rank     (fmpz_t res, const char *obj,    fcomb_ctx_t ctx);
char *fcomb_unrank   (fmpz_t rank,                    fcomb_ctx_t ctx);
char *fcomb_draw     (fcomb_ctx_t ctx);

/* fcomb_boltzmann builds its own internal solve context with larger headroom;
   fcomb_solve() does not need to be called beforehand.
   tolerance: accept objects of size n ± tolerance.
   max_attempts: rejection-sampling budget (0 → sensible default). */
char *fcomb_boltzmann(int tolerance, int max_attempts, fcomb_ctx_t ctx);

#endif
