#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "api.h"
#include "grammar/object.h"
#include "solver/boltzmann/eval.h"
#include "solver/draw/eval.h"
#include "solver/rank/eval.h"
#include "solver/unrank/eval.h"

extern Spec *readSpec(char *filename);
extern Spec *readSpecStr(const char *input);

/* ---- internal helpers ---- */

static void drop_solve_ctx(fcomb_ctx_t ctx) {
  if (ctx->solve_ctx) {
    context_clear(ctx->solve_ctx);
    free(ctx->solve_ctx);
    ctx->solve_ctx = NULL;
  }
}

static char *resolve_symbol(fcomb_ctx_t ctx) {
  if (ctx->symbol)
    return ctx->symbol;
  Rules *sl = (Rules *)ctx->spec->component;
  return sl->components[0]->variable->name;
}

static Entry *find_entry(fcomb_ctx_t ctx, const char *symbol) {
  for (int i = 0; i < ctx->solve_ctx->num_entries; i++) {
    if (strcmp(ctx->solve_ctx->entries[i].name, symbol) == 0)
      return &ctx->solve_ctx->entries[i];
  }
  return NULL;
}

static void seed_rng(flint_rand_t rng, ulong seed) {
  flint_rand_init(rng);
  if (seed == 0) {
    srand((unsigned)time(NULL) ^ (unsigned)getpid());
    seed = (ulong)rand();
  }
  flint_rand_set_seed(rng, seed, seed ^ 0xdeadbeefUL);
}

/* ---- lifecycle ---- */

fcomb_ctx_t fcomb_ctx_init(void) {
  fcomb_ctx_t ctx = malloc(sizeof(struct fcomb_ctx_s));
  if (!ctx)
    return NULL;
  ctx->spec      = NULL;
  ctx->solve_ctx = NULL;
  ctx->symbol    = NULL;
  ctx->n         = -1;
  ctx->labeled   = 0;
  ctx->seed      = 0;
  ctx->order     = ORDER_LEX;
  return ctx;
}

void fcomb_ctx_destroy(fcomb_ctx_t ctx) {
  if (!ctx)
    return;
  drop_solve_ctx(ctx);
  free(ctx->symbol);
  free(ctx);
}

/* ---- setters ---- */

int fcomb_ctx_set_spec_str(const char *spec, fcomb_ctx_t ctx) {
  Spec *g = readSpecStr(spec);
  if (!g || g->type == ISERROR)
    return 1;
  ctx->spec = g;
  drop_solve_ctx(ctx);
  return 0;
}

int fcomb_ctx_set_spec_json(const char *path, fcomb_ctx_t ctx) {
  Spec *g = readSpec((char *)path);
  if (!g || g->type == ISERROR)
    return 1;
  ctx->spec = g;
  drop_solve_ctx(ctx);
  return 0;
}

int fcomb_ctx_set_symbol(const char *sym, fcomb_ctx_t ctx) {
  free(ctx->symbol);
  if (sym) {
    ctx->symbol = strdup(sym);
    if (!ctx->symbol)
      return -1;
  } else {
    ctx->symbol = NULL;
  }
  return 0;
}

int fcomb_ctx_set_size(int n, fcomb_ctx_t ctx) {
  ctx->n = n;
  drop_solve_ctx(ctx);
  return 0;
}

int fcomb_ctx_set_labeled(int labeled, fcomb_ctx_t ctx) {
  ctx->labeled = labeled;
  drop_solve_ctx(ctx);
  return 0;
}

void fcomb_ctx_set_seed(ulong seed, fcomb_ctx_t ctx) {
  ctx->seed = seed;
}

int fcomb_ctx_set_order(order_t order, fcomb_ctx_t ctx) {
  ctx->order = order;
  return 0;
}

/* ---- solve ---- */

int fcomb_solve(fcomb_ctx_t ctx) {
  if (!ctx->spec || ctx->n < 0)
    return 1;
  drop_solve_ctx(ctx);
  /* n+1 so the coefficient at index n is always accessible */
  ctx->solve_ctx = solve_spec(ctx->spec, ctx->n + 1, ctx->labeled);
  return ctx->solve_ctx ? 0 : 1;
}

/* ---- operations ---- */

void fcomb_count(fmpz_t res, fcomb_ctx_t ctx) {
  Entry *e = find_entry(ctx, resolve_symbol(ctx));
  if (!e) {
    fmpz_set_si(res, -1);
    return;
  }
  fmpz_poly_get_coeff_fmpz(res, e->poly, ctx->n);
}

void fcomb_rank(fmpz_t res, const char *obj_str, fcomb_ctx_t ctx) {
  Object *obj = parse_object((char *)obj_str);
  if (!obj) {
    fmpz_set_si(res, -1);
    return;
  }
  ctx->solve_ctx->order = ctx->order;
  rank(ctx->solve_ctx, resolve_symbol(ctx), obj, res);
  free_object(obj);
}

char *fcomb_unrank(fmpz_t r, fcomb_ctx_t ctx) {
  ctx->solve_ctx->order = ctx->order;
  return unrank(ctx->solve_ctx, resolve_symbol(ctx), ctx->n, r);
}

char *fcomb_draw(fcomb_ctx_t ctx) {
  flint_rand_t rng;
  seed_rng(rng, ctx->seed);
  char *res = draw(ctx->solve_ctx, resolve_symbol(ctx), ctx->n, rng);
  flint_rand_clear(rng);
  return res;
}

char *fcomb_boltzmann(int tolerance, int max_attempts, fcomb_ctx_t ctx) {
  if (!ctx->spec || ctx->n < 0)
    return NULL;

  int target     = ctx->n;
  int log_hdroom = (int)(target * log((double)(target + 1)) * 1.5) + 100;
  int max_n      = target * 2 + 50;
  if (log_hdroom > max_n)
    max_n = log_hdroom;

  if (max_attempts <= 0) {
    max_attempts = target * 5000;
    if (max_attempts < 100000)
      max_attempts = 100000;
  }

  Context *bctx = solve_spec(ctx->spec, max_n, ctx->labeled);
  if (!bctx)
    return NULL;

  flint_rand_t rng;
  seed_rng(rng, ctx->seed);

  int actual_size = 0;
  char *res =
      boltzmann_generate(bctx, resolve_symbol(ctx), target, tolerance, rng,
                         max_attempts, &actual_size, max_n, ctx->labeled);

  flint_rand_clear(rng);
  context_clear(bctx);
  free(bctx);
  return res;
}
