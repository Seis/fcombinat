#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "api.h"

typedef struct {
  char *input_spec;
  char *symbol_arg;
  char *unrank_str;
  char *rank_str;
  int   max_n;
  int   is_labeled;
  int   is_poly;
  int   is_verbose;
  int   is_unrank;
  int   is_rank;
  int   is_bijection;
  int   is_draw;
  int   is_boltzmann;
  int   is_enum_all;
  int   is_terms;
  ulong seed;
} Config;

static void print_usage(char *prog) {
  printf("Usage: %s [OPTIONS] [ACTION]\n", prog);
  printf("  -j, -g, --parse-json <FILE|SPEC>  JSON spec file or inline spec string\n");
  printf("  -n, --size <N>                    Size parameter (required)\n");
  printf("  -l, --labeled                     Labeled combinatorial class\n");
  printf("  -P, --poly                        Output as polynomial\n");
  printf("  -v, --verbose                     Verbose output\n");
  printf("  --symbol <S>                      Target symbol (default: first rule)\n");
  printf("  --seed <N>                        RNG seed for --draw and --boltzmann (0 = random)\n");
  printf("\nActions:\n");
  printf("  --terms                           Print first N coefficients (default)\n");
  printf("  --rank <OBJ>                      Compute rank of object string\n");
  printf("  --unrank <INT>                    Generate object at rank\n");
  printf("  --draw                            Uniformly random object of size N\n");
  printf("  --boltzmann                       Approximate-size Boltzmann sample\n");
  printf("  --enum-all                        Enumerate all objects of size N\n");
  printf("  --bijection                       Bijection roundtrip test\n");
}

static fcomb_ctx_t build_ctx(Config *cfg) {
  fcomb_ctx_t ctx = fcomb_ctx_init();
  if (!ctx)
    return NULL;

  FILE *f = fopen(cfg->input_spec, "r");
  int is_file = (f != NULL);
  if (f)
    fclose(f);

  int err = is_file ? fcomb_ctx_set_spec_json(cfg->input_spec, ctx)
                    : fcomb_ctx_set_spec_str(cfg->input_spec, ctx);
  if (err) {
    fprintf(stderr, "Error: failed to parse spec '%s'\n", cfg->input_spec);
    fcomb_ctx_destroy(ctx);
    return NULL;
  }

  fcomb_ctx_set_size(cfg->max_n, ctx);
  fcomb_ctx_set_labeled(cfg->is_labeled, ctx);

  if (cfg->symbol_arg)
    fcomb_ctx_set_symbol(cfg->symbol_arg, ctx);

  if (cfg->seed)
    fcomb_ctx_set_seed(cfg->seed, ctx);

  return ctx;
}

/* ---- actions ---- */

static int run_terms(fcomb_ctx_t ctx, Config *cfg) {
  if (fcomb_solve(ctx))
    return 1;

  Context *sc = ctx->solve_ctx;
  for (int i = 0; i < sc->num_entries; i++) {
    Entry *e = &sc->entries[i];

    if (cfg->symbol_arg) {
      if (strcmp(e->name, cfg->symbol_arg) != 0)
        continue;
    } else {
      if (strcmp(e->name, "Z") == 0)
        continue;
    }

    if (cfg->is_poly) {
      printf("%s: ", e->name);
      int started = 0;
      for (int k = 0; k < cfg->max_n; k++) {
        fmpz_t c;
        fmpz_init(c);
        fmpz_poly_get_coeff_fmpz(c, e->poly, k);
        if (!fmpz_is_zero(c)) {
          if (started && fmpz_sgn(c) > 0)
            printf("+");
          if (k > 0) {
            if (!fmpz_is_one(c)) {
              if (fmpz_equal_si(c, -1))
                printf("-");
              else {
                fmpz_print(c);
                if (fmpz_sgn(c) > 0)
                  printf("* ");
              }
            }
            printf("z");
            if (k > 1)
              printf("^%d", k);
          } else {
            fmpz_print(c);
          }
          started = 1;
        }
        fmpz_clear(c);
      }
      if (!started)
        printf("0");
      printf("\n");
    } else {
      for (int k = 0; k < cfg->max_n; k++) {
        fmpz_t c;
        fmpz_init(c);
        fmpz_poly_get_coeff_fmpz(c, e->poly, k);
        fmpz_print(c);
        printf("\n");
        fmpz_clear(c);
      }
    }
  }
  return 0;
}

static int run_rank(fcomb_ctx_t ctx, Config *cfg) {
  if (fcomb_solve(ctx))
    return 1;

  fmpz_t r;
  fmpz_init(r);
  fcomb_rank(r, cfg->rank_str, ctx);

  if (fmpz_sgn(r) < 0) {
    fprintf(stderr, "Error: ranking failed\n");
    fmpz_clear(r);
    return 1;
  }
  fmpz_print(r);
  printf("\n");
  fmpz_clear(r);
  return 0;
}

static int run_unrank(fcomb_ctx_t ctx, Config *cfg) {
  if (fcomb_solve(ctx))
    return 1;

  fmpz_t r;
  fmpz_init(r);
  if (fmpz_set_str(r, cfg->unrank_str, 10) != 0) {
    fprintf(stderr, "Error: invalid rank integer '%s'\n", cfg->unrank_str);
    fmpz_clear(r);
    return 1;
  }

  char *obj = fcomb_unrank(r, ctx);
  fmpz_clear(r);
  if (!obj)
    return 1;
  printf("%s\n", obj);
  free(obj);
  return 0;
}

static int run_draw(fcomb_ctx_t ctx) {
  if (fcomb_solve(ctx))
    return 1;

  char *obj = fcomb_draw(ctx);
  if (!obj)
    return 1;
  printf("%s\n", obj);
  free(obj);
  return 0;
}

static int run_boltzmann(fcomb_ctx_t ctx, Config *cfg) {
  int tolerance = cfg->max_n / 5;
  if (tolerance < 1)
    tolerance = 1;

  char *obj = fcomb_boltzmann(tolerance, 0, ctx);
  if (!obj) {
    fprintf(stderr, "Boltzmann: failed to generate object of size %d±%d\n",
            cfg->max_n, tolerance);
    return 1;
  }
  printf("%s\n", obj);
  if (cfg->is_verbose)
    fprintf(stderr, "Boltzmann: target=%d tolerance=%d\n", cfg->max_n, tolerance);
  free(obj);
  return 0;
}

static int run_enum_all(fcomb_ctx_t ctx, Config *cfg) {
  if (fcomb_solve(ctx))
    return 1;

  fmpz_t count, r;
  fmpz_init(count);
  fmpz_init(r);
  fcomb_count(count, ctx);

  if (!fmpz_is_zero(count) && fmpz_sgn(count) > 0) {
    fmpz_set_ui(r, 0);
    while (fmpz_cmp(r, count) < 0) {
      char *obj = fcomb_unrank(r, ctx);
      if (obj) {
        printf("%s\n", obj);
        free(obj);
      }
      fmpz_add_ui(r, r, 1);
    }
  }

  fmpz_clear(count);
  fmpz_clear(r);
  return 0;
}

static int run_bijection(fcomb_ctx_t ctx, Config *cfg) {
  if (!cfg->symbol_arg) {
    fprintf(stderr, "Error: --bijection requires --symbol\n");
    return 1;
  }
  if (fcomb_solve(ctx))
    return 1;

  int failed = 0;
  fmpz_t count, r_in, r_out;
  fmpz_init(count);
  fmpz_init(r_in);
  fmpz_init(r_out);

  /* Test all sizes 1..max_n-1 */
  for (int n = 1; n < cfg->max_n; n++) {
    /* Temporarily adjust ctx->n so fcomb_count works at size n.
       solve_ctx covers 0..max_n so coefficients are available. */
    int saved_n = ctx->n;
    ctx->n = n;
    fcomb_count(count, ctx);
    ctx->n = saved_n;

    if (fmpz_is_zero(count))
      continue;

    int small    = fmpz_cmp_si(count, 50) <= 0;
    int num_tests = small ? (int)fmpz_get_si(count) : 10;

    for (int t = 0; t < num_tests; t++) {
      if (small) {
        fmpz_set_si(r_in, t);
      } else {
        if      (t == 0) fmpz_set_si(r_in, 0);
        else if (t == 1) fmpz_set_si(r_in, 1);
        else if (t == 2) fmpz_sub_si(r_in, count, 1);
        else if (t == 3) fmpz_fdiv_q_si(r_in, count, 2);
        else {
          fmpz_set_si(r_in, rand() % 10000);
          fmpz_mod(r_in, r_in, count);
        }
      }

      int saved_n2 = ctx->n;
      ctx->n = n;
      char *str = fcomb_unrank(r_in, ctx);
      ctx->n = saved_n2;

      if (!str) {
        fprintf(stderr, "Unrank failed n=%d rank=", n);
        fmpz_fprint(stderr, r_in);
        fprintf(stderr, "\n");
        failed = 1;
        continue;
      }

      fcomb_rank(r_out, str, ctx);

      if (fmpz_sgn(r_out) < 0) {
        fprintf(stderr, "Rank error n=%d obj=%s\n", n, str);
        failed = 1;
      } else if (!fmpz_equal(r_in, r_out)) {
        fprintf(stderr, "Bijection FAIL n=%d: expected ");
        fmpz_fprint(stderr, r_in);
        fprintf(stderr, " got ");
        fmpz_fprint(stderr, r_out);
        fprintf(stderr, " (obj=%s)\n", str);
        failed = 1;
      }
      free(str);
    }
  }

  fmpz_clear(count);
  fmpz_clear(r_in);
  fmpz_clear(r_out);
  return failed ? 2 : 0;
}

/* ---- main ---- */

int main(int argc, char *argv[]) {
  Config cfg = {.max_n = -1};

  static struct option long_opts[] = {
      {"grammar",    required_argument, 0, 'j'},
      {"parse-json", required_argument, 0, 'j'},
      {"size",       required_argument, 0, 'n'},
      {"labeled",    no_argument,       0, 'l'},
      {"labelled",   no_argument,       0, 'l'},
      {"poly",       no_argument,       0, 'P'},
      {"verbose",    no_argument,       0, 'v'},
      {"symbol",     required_argument, 0, 'S'},
      {"unrank",     required_argument, 0, 'u'},
      {"rank",       required_argument, 0, 'r'},
      {"draw",       no_argument,       0, 'd'},
      {"boltzmann",  no_argument,       0, 'B'},
      {"enum-all",   no_argument,       0, 'e'},
      {"terms",      no_argument,       0, 't'},
      {"bijection",  no_argument,       0, 'b'},
      {"seed",       required_argument, 0, 's'},
      {0, 0, 0, 0},
  };

  int opt, opt_idx = 0, action = 0;
  while ((opt = getopt_long(argc, argv, "g:j:n:lLPvS:u:r:dBetbs:",
                            long_opts, &opt_idx)) != -1) {
    switch (opt) {
    case 'g': case 'j': cfg.input_spec = optarg; break;
    case 'n': {
      char *end;
      long v = strtol(optarg, &end, 10);
      if (*end != '\0' || v <= 0 || v > INT_MAX) {
        fprintf(stderr, "Error: -n requires a positive integer\n");
        return 1;
      }
      cfg.max_n = (int)v;
      break;
    }
    case 'l': case 'L': cfg.is_labeled = 1;  break;
    case 'P':  cfg.is_poly    = 1;             break;
    case 'v':  cfg.is_verbose = 1;             break;
    case 'S':  cfg.symbol_arg = optarg;        break;
    case 'u':  cfg.is_unrank  = 1; cfg.unrank_str = optarg; action = 1; break;
    case 'r':  cfg.is_rank    = 1; cfg.rank_str   = optarg; action = 1; break;
    case 'd':  cfg.is_draw      = 1; action = 1; break;
    case 'B':  cfg.is_boltzmann = 1; action = 1; break;
    case 'e':  cfg.is_enum_all  = 1; action = 1; break;
    case 't':  cfg.is_terms     = 1; action = 1; break;
    case 'b':  cfg.is_bijection = 1; action = 1; break;
    case 's': {
      char *end;
      cfg.seed = strtoul(optarg, &end, 10);
      if (*end != '\0') {
        fprintf(stderr, "Error: --seed requires a non-negative integer\n");
        return 1;
      }
      break;
    }
    default:   print_usage(argv[0]); return 1;
    }
  }

  if (optind < argc) {
    fprintf(stderr, "Error: unexpected argument '%s'\n", argv[optind]);
    print_usage(argv[0]);
    return 1;
  }
  if (!cfg.input_spec) {
    fprintf(stderr, "Error: no spec given (-j / -g)\n");
    print_usage(argv[0]);
    return 1;
  }
  if (cfg.max_n < 0) {
    fprintf(stderr, "Error: -n <N> required\n");
    print_usage(argv[0]);
    return 1;
  }
  if (!action) {
    fprintf(stderr, "Error: no action specified\n");
    print_usage(argv[0]);
    return 1;
  }

  if (cfg.is_verbose) {
    fprintf(stderr, "Mode: %s\n", cfg.is_labeled ? "labeled" : "unlabeled");
    fprintf(stderr, "N:    %d\n", cfg.max_n);
  }

  fcomb_ctx_t ctx = build_ctx(&cfg);
  if (!ctx)
    return 1;

  int ret = 0;
  if      (cfg.is_unrank)    ret = run_unrank(ctx, &cfg);
  else if (cfg.is_rank)      ret = run_rank(ctx, &cfg);
  else if (cfg.is_draw)      ret = run_draw(ctx);
  else if (cfg.is_boltzmann) ret = run_boltzmann(ctx, &cfg);
  else if (cfg.is_enum_all)  ret = run_enum_all(ctx, &cfg);
  else if (cfg.is_bijection) ret = run_bijection(ctx, &cfg);
  else                       ret = run_terms(ctx, &cfg);

  fcomb_ctx_destroy(ctx);
  return ret;
}
