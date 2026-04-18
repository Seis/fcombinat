#include "solver/boltzmann/eval.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grammar/object.h"
#include "solver/math.h"
#include "solver/boltzmann/oracle.h"
#include "solver/boltzmann/tuner.h"
#include "solver/count/eval.h"

// Generate a uniform random double in [0, 1)
static double rand_uniform(flint_rand_t state) {
  fmpz_t big;
  fmpz_init(big);
  fmpz_t mod;
  fmpz_init(mod);
  fmpz_set_ui(mod, 1000000000UL);
  fmpz_randm(big, state, mod);
  double r = fmpz_get_d(big) / 1e9;
  fmpz_clear(big);
  fmpz_clear(mod);
  return r;
}

// Poisson random variable with parameter lambda
static int rand_poisson(double lambda, flint_rand_t state) {
  if (lambda <= 0.0)
    return 0;
  // Knuth's algorithm
  double L = exp(-lambda);
  int k = 0;
  double p = 1.0;
  do {
    k++;
    p *= rand_uniform(state);
  } while (p > L);
  return k - 1;
}

// Geometric random variable: P(X=k) = (1-p)*p^k for k=0,1,2,...
// Expected value = p/(1-p)
static int rand_geometric(double p, flint_rand_t state) {
  if (p <= 0.0)
    return 0;
  if (p >= 1.0)
    return 1000; // Cap
  // Using inverse CDF: k = floor(ln(U) / ln(p))
  double u = rand_uniform(state);
  if (u <= 0.0)
    return 0;
  return (int)(log(u) / log(p));
}

// Returns the minimum cardinality if expr is Sequence(Z) (possibly with geq
// restriction), so that the Pólya-Boltzmann decomposition can be applied.
// Returns -1 if the child is not a plain atom sequence.
static int atom_seq_min_card(Expr *expr) {
  if (expr->type != SEQUENCE)
    return -1;
  Expr *inner = (Expr *)expr->component;
  if (inner->type != Z && inner->type != ATOM)
    return -1;
  // GREATER means card >= limit (enforced in the sampling loop)
  if (expr->restriction == GREATER)
    return (int)expr->limit;
  if (expr->restriction == EQUAL)
    return (int)expr->limit;
  return 0; // NONE: sequences of any length >= 0
}

// Build "Seq(Z(0), ..., Z(0))" with k atoms, for unlabeled structure output.
static char *make_atom_sequence_str(int k) {
  if (k == 0)
    return strdup("Seq()");
  // "Seq(" + k*"Z(0)" + (k-1)*", " + ")" + "\0"  =>  4 + 6k chars
  char *buf = malloc(4 + 6 * k + 2);
  strcpy(buf, "Seq(");
  for (int i = 0; i < k; i++) {
    if (i > 0)
      strcat(buf, ", ");
    strcat(buf, "Z(0)");
  }
  strcat(buf, ")");
  return buf;
}

// Build string result for atom
static char *make_atom_str(int label) {
  char buf[32];
  sprintf(buf, "Z(%d)", label);
  return strdup(buf);
}

// Internal recursive Boltzmann draw
// Returns an object string, sets *size to the actual size of the object
char *boltzmann_draw_e(Oracle *orc, Context *ctx, Expr *expr, double x,
                       flint_rand_t state, int *actual_size, int is_labeled) {
  switch (expr->type) {
  case ATOM:
  case Z: {
    *actual_size = 1;
    return make_atom_str(0); // Labels are assigned post-hoc
  }

  case EPSILON: {
    *actual_size = 0;
    return strdup("Eps()");
  }

  case ID: {
    Id *id = (Id *)expr->component;
    // Look up the expression for this ID
    for (int i = 0; i < ctx->num_entries; i++) {
      if (strcmp(ctx->entries[i].name, id->name) == 0) {
        return boltzmann_draw_e(orc, ctx, ctx->entries[i].expr, x, state,
                                actual_size, is_labeled);
      }
    }
    if (strcmp(id->name, "Z") == 0) {
      *actual_size = 1;
      return make_atom_str(0);
    }
    *actual_size = 0;
    return strdup("ErrorIDNotFound");
  }

  case UNION: {
    // Choose branch i with probability A_i(x) / (A_0(x) + A_1(x) + ...)
    ExprList *el = (ExprList *)expr->component;
    double total = 0.0;
    double *weights = malloc(sizeof(double) * el->size);
    for (int i = 0; i < el->size; i++) {
      weights[i] = oracle_eval_expr(orc, ctx, el->components[i], x,
                                    orc->max_n, is_labeled);
      if (weights[i] < 0.0)
        weights[i] = 0.0;
      total += weights[i];
    }

    if (total <= 0.0) {
      free(weights);
      *actual_size = 0;
      return strdup("ErrorUnionEmpty");
    }

    double r = rand_uniform(state) * total;
    int chosen = 0;
    double cum = 0.0;
    for (int i = 0; i < el->size; i++) {
      cum += weights[i];
      if (r < cum) {
        chosen = i;
        break;
      }
    }
    free(weights);

    char *child =
        boltzmann_draw_e(orc, ctx, el->components[chosen], x, state,
                         actual_size, is_labeled);
    char *res = malloc(strlen(child) + 30);
    sprintf(res, "Union_%d(%s)", chosen, child);
    free(child);
    return res;
  }

  case PROD: {
    // Generate from A independently, then from B independently
    // Sizes add up
    ExprList *el = (ExprList *)expr->component;
    if (el->size == 0) {
      *actual_size = 0;
      return strdup("Eps()");
    }
    if (el->size == 1)
      return boltzmann_draw_e(orc, ctx, el->components[0], x, state,
                              actual_size, is_labeled);

    int total_size = 0;
    char **parts = malloc(sizeof(char *) * el->size);
    int *sizes = malloc(sizeof(int) * el->size);

    for (int i = 0; i < el->size; i++) {
      int sz = 0;
      parts[i] = boltzmann_draw_e(orc, ctx, el->components[i], x, state, &sz,
                                  is_labeled);
      sizes[i] = sz;
      total_size += sz;
    }

    // Build result string
    int len = 10;
    for (int i = 0; i < el->size; i++)
      len += strlen(parts[i]) + 2;
    char *res = malloc(len);
    strcpy(res, "Prod(");
    for (int i = 0; i < el->size; i++) {
      if (i > 0)
        strcat(res, ", ");
      strcat(res, parts[i]);
      free(parts[i]);
    }
    strcat(res, ")");
    free(parts);
    free(sizes);

    *actual_size = total_size;
    return res;
  }

  case SET: {
    Expr *child = (Expr *)expr->component;
    double child_val = oracle_eval_expr(orc, ctx, child, x, orc->max_n,
                                        is_labeled);

    // Unlabeled SET = MSet: Pólya-Boltzmann via Fristedt's algorithm.
    // When child is Seq(Z, card >= k_min), draw m_k ~ Geometric_0(x^k) for
    // each part size k; each drawn part is output as a sequence of k atoms.
    if (!is_labeled) {
      int min_k = atom_seq_min_card(child);
      if (min_k >= 0) {
        int k_start = (min_k > 0) ? min_k : 1;
        int cap = 64;
        char **polya_parts = malloc(sizeof(char *) * cap);
        int nparts = 0, total_sz = 0, overflow = 0;
        double xk = 1.0;
        for (int k = 1; k < k_start; k++)
          xk *= x; // advance to x^(k_start-1)
        for (int k = k_start; k <= orc->max_n && !overflow; k++) {
          xk *= x; // xk = x^k
          if (xk < 1e-300)
            break;
          int mk = rand_geometric(xk, state);
          for (int j = 0; j < mk && !overflow; j++) {
            if (nparts >= cap) {
              cap *= 2;
              polya_parts = xrealloc(polya_parts, sizeof(char *) * cap);
            }
            polya_parts[nparts++] = make_atom_sequence_str(k);
            total_sz += k;
            if (total_sz > orc->max_n * 4)
              overflow = 1;
          }
        }
        int len = 10;
        for (int i = 0; i < nparts; i++)
          len += strlen(polya_parts[i]) + 2;
        char *res = malloc(len);
        strcpy(res, "Set(");
        for (int i = 0; i < nparts; i++) {
          if (i > 0)
            strcat(res, ", ");
          strcat(res, polya_parts[i]);
          free(polya_parts[i]);
        }
        strcat(res, ")");
        free(polya_parts);
        *actual_size = total_sz;
        return res;
      }
    }

    // Labeled SET (EGF): number of components ~ Poisson(A(x)).
    // Unlabeled SET with non-atom-sequence child: Poisson approximation.
    int num_components;

    // Apply cardinality restriction via rejection sampling
    if (expr->restriction == GREATER) {
      // Rejection: redraw Poisson until >= limit
      int limit = (int)expr->limit;
      int max_redraw = 10000;
      do {
        num_components = rand_poisson(child_val, state);
      } while (num_components < limit && --max_redraw > 0);
      if (max_redraw <= 0)
        num_components = limit; // Fallback
    } else if (expr->restriction == LESS) {
      // Rejection: redraw Poisson until < limit
      int limit = (int)expr->limit;
      int max_redraw = 10000;
      do {
        num_components = rand_poisson(child_val, state);
      } while (num_components >= limit && --max_redraw > 0);
      if (max_redraw <= 0)
        num_components = limit - 1;
    } else if (expr->restriction == EQUAL) {
      num_components = (int)expr->limit;
    } else {
      num_components = rand_poisson(child_val, state);
    }

    if (num_components <= 0) {
      *actual_size = 0;
      return strdup("Set()");
    }

    int total_size = 0;
    char **parts = malloc(sizeof(char *) * num_components);
    for (int i = 0; i < num_components; i++) {
      int sz = 0;
      parts[i] = boltzmann_draw_e(orc, ctx, child, x, state, &sz, is_labeled);
      total_size += sz;
    }

    int len = 10;
    for (int i = 0; i < num_components; i++)
      len += strlen(parts[i]) + 2;
    char *res = malloc(len);
    strcpy(res, "Set(");
    for (int i = 0; i < num_components; i++) {
      if (i > 0)
        strcat(res, ", ");
      strcat(res, parts[i]);
      free(parts[i]);
    }
    strcat(res, ")");
    free(parts);

    *actual_size = total_size;
    return res;
  }

  case SEQUENCE: {
    // Boltzmann for Seq(B): k components with P(k) = (1-B(x)) * B(x)^k
    // This is a geometric distribution with parameter p = B(x)
    Expr *child = (Expr *)expr->component;
    double child_val = oracle_eval_expr(orc, ctx, child, x, orc->max_n,
                                        is_labeled);
    // Guard: need B(x) < 1 for the geometric series to converge
    if (child_val >= 1.0)
      child_val = 0.999;
    double p = child_val;
    int length = rand_geometric(p, state);

    // Apply restriction
    if (expr->restriction == GREATER && length < (int)expr->limit) {
      length = (int)expr->limit;
    } else if (expr->restriction == LESS && length >= (int)expr->limit) {
      length = (int)expr->limit - 1;
    }

    if (length <= 0) {
      *actual_size = 0;
      return strdup("Seq()");
    }

    int total_size = 0;
    char **parts = malloc(sizeof(char *) * length);
    for (int i = 0; i < length; i++) {
      int sz = 0;
      parts[i] = boltzmann_draw_e(orc, ctx, child, x, state, &sz, is_labeled);
      total_size += sz;
    }

    int len = 10;
    for (int i = 0; i < length; i++)
      len += strlen(parts[i]) + 2;
    char *res = malloc(len);
    strcpy(res, "Seq(");
    for (int i = 0; i < length; i++) {
      if (i > 0)
        strcat(res, ", ");
      strcat(res, parts[i]);
      free(parts[i]);
    }
    strcat(res, ")");
    free(parts);

    *actual_size = total_size;
    return res;
  }

  case CYCLE: {
    // For labeled: similar to Sequence but divided by cycle symmetry
    // Use Logarithmic distribution: P(X=k) = -p^k / (k*ln(1-p))
    // where p = A(x)
    Expr *child = (Expr *)expr->component;
    double child_val = oracle_eval_expr(orc, ctx, child, x, orc->max_n,
                                        is_labeled);

    // Logarithmic distribution via Kemp's algorithm
    int length = 1; // Minimum cycle length
    if (child_val > 0.0 && child_val < 1.0) {
      double u = rand_uniform(state);
      double p = child_val;
      double q = 1.0 - p;
      double s = p;
      double prob = -p / log(q);
      double cum = prob;
      length = 1;
      while (u > cum && length < 1000) {
        length++;
        s *= p;
        prob = s / (length * (-log(q)));
        cum += prob;
      }
    }

    if (length <= 0)
      length = 1;

    int total_size = 0;
    char **parts = malloc(sizeof(char *) * length);
    for (int i = 0; i < length; i++) {
      int sz = 0;
      parts[i] = boltzmann_draw_e(orc, ctx, child, x, state, &sz, is_labeled);
      total_size += sz;
    }

    int len = 12;
    for (int i = 0; i < length; i++)
      len += strlen(parts[i]) + 2;
    char *res = malloc(len);
    strcpy(res, "Cycle(");
    for (int i = 0; i < length; i++) {
      if (i > 0)
        strcat(res, ", ");
      strcat(res, parts[i]);
      free(parts[i]);
    }
    strcat(res, ")");
    free(parts);

    *actual_size = total_size;
    return res;
  }

  case POWERSET: {
    // Pólya-Boltzmann for PowerSet(Seq(Z, card >= k_min)):
    // include part of size k independently with probability x^k / (1 + x^k)
    // (Bernoulli sampling - each part size appears at most once).
    Expr *child = (Expr *)expr->component;
    double child_val = oracle_eval_expr(orc, ctx, child, x, orc->max_n,
                                        is_labeled);
    int min_k = atom_seq_min_card(child);
    if (min_k >= 0) {
      int k_start = (min_k > 0) ? min_k : 1;
      int cap = 64;
      char **bern_parts = malloc(sizeof(char *) * cap);
      int nparts = 0, total_sz = 0;
      double xk = 1.0;
      for (int k = 1; k < k_start; k++)
        xk *= x; // advance to x^(k_start-1)
      for (int k = k_start; k <= orc->max_n; k++) {
        xk *= x; // xk = x^k
        if (xk < 1e-300)
          break;
        double prob = xk / (1.0 + xk);
        if (rand_uniform(state) < prob) {
          if (nparts >= cap) {
            cap *= 2;
            bern_parts = xrealloc(bern_parts, sizeof(char *) * cap);
          }
          bern_parts[nparts++] = make_atom_sequence_str(k);
          total_sz += k;
          if (total_sz > orc->max_n * 4)
            break;
        }
      }
      int blen = 15;
      for (int i = 0; i < nparts; i++)
        blen += strlen(bern_parts[i]) + 2;
      char *bres = malloc(blen);
      strcpy(bres, "PowerSet(");
      for (int i = 0; i < nparts; i++) {
        if (i > 0)
          strcat(bres, ", ");
        strcat(bres, bern_parts[i]);
        free(bern_parts[i]);
      }
      strcat(bres, ")");
      free(bern_parts);
      *actual_size = total_sz;
      return bres;
    }

    // Fallback: Poisson approximation for non-atom-sequence child
    int num_components = rand_poisson(child_val, state);
    if (num_components <= 0) {
      *actual_size = 0;
      return strdup("PowerSet()");
    }

    int total_size = 0;
    char **fb_parts = malloc(sizeof(char *) * num_components);
    for (int i = 0; i < num_components; i++) {
      int sz = 0;
      fb_parts[i] =
          boltzmann_draw_e(orc, ctx, child, x, state, &sz, is_labeled);
      total_size += sz;
    }

    int fb_len = 15;
    for (int i = 0; i < num_components; i++)
      fb_len += strlen(fb_parts[i]) + 2;
    char *fb_res = malloc(fb_len);
    strcpy(fb_res, "PowerSet(");
    for (int i = 0; i < num_components; i++) {
      if (i > 0)
        strcat(fb_res, ", ");
      strcat(fb_res, fb_parts[i]);
      free(fb_parts[i]);
    }
    strcat(fb_res, ")");
    free(fb_parts);

    *actual_size = total_size;
    return fb_res;
  }

  default:
    *actual_size = 0;
    return strdup("ErrorUnknownType");
  }
}

// Relabel: assign random labels {1..n} to atoms in a Boltzmann-generated object
static void collect_atom_count(const char *str, int *count) {
  // Count occurrences of "Z(" in the string
  *count = 0;
  const char *p = str;
  while ((p = strstr(p, "Z(")) != NULL) {
    (*count)++;
    p += 2;
  }
}

static char *relabel_object(const char *str, int n, flint_rand_t state) {
  // Assign labels 1..n to atoms in order (for labeled case, this is
  // already a valid random labeling since the Boltzmann process generates
  // structures uniformly)
  int count = 0;
  collect_atom_count(str, &count);
  if (count != n)
    return strdup(str); // Size mismatch, return as-is

  // Create a random permutation of {1..n}
  int *perm = xmalloc(sizeof(int) * n);
  for (int i = 0; i < n; i++)
    perm[i] = i + 1;
  // Fisher-Yates shuffle
  for (int i = n - 1; i > 0; i--) {
    fmpz_t idx;
    fmpz_init(idx);
    fmpz_t mod;
    fmpz_init(mod);
    fmpz_set_ui(mod, i + 1);
    fmpz_randm(idx, state, mod);
    int j = fmpz_get_ui(idx);
    int tmp = perm[i];
    perm[i] = perm[j];
    perm[j] = tmp;
    fmpz_clear(idx);
    fmpz_clear(mod);
  }

  // Replace Z(0) occurrences with Z(perm[k])
  char *result = xmalloc(strlen(str) + n * 10);
  char *out = result;
  const char *p = str;
  int label_idx = 0;

  while (*p) {
    if (p[0] == 'Z' && p[1] == '(') {
      out += sprintf(out, "Z(%d)", perm[label_idx++]);
      // Skip past "Z(0)"
      p += 2;
      while (*p && *p != ')')
        p++;
      if (*p == ')')
        p++;
    } else {
      *out++ = *p++;
    }
  }
  *out = '\0';

  free(perm);
  return result;
}

char *boltzmann_generate(Context *ctx, char *symbol, int target_n,
                         int tolerance, flint_rand_t state, int max_attempts,
                         int *actual_size, int max_n, int is_labeled) {
  // 1. Tune the parameter x
  double x = boltzmann_tune(ctx, symbol, target_n, max_n, is_labeled);

  if (x <= 0.0 || isnan(x) || isinf(x)) {
    fprintf(stderr, "Boltzmann: tuning failed (x=%f)\n", x);
    return NULL;
  }

  // 2. Create oracle
  Oracle *orc = oracle_create(ctx, max_n, x, is_labeled);

  // 3. Find root expression
  Expr *root_expr = NULL;
  for (int i = 0; i < ctx->num_entries; i++) {
    if (strcmp(ctx->entries[i].name, symbol) == 0) {
      root_expr = ctx->entries[i].expr;
      break;
    }
  }
  if (!root_expr) {
    oracle_free(orc);
    return NULL;
  }

  // 4. Rejection loop
  for (int attempt = 0; attempt < max_attempts; attempt++) {
    int size = 0;
    char *obj = boltzmann_draw_e(orc, ctx, root_expr, x, state, &size,
                                 is_labeled);

    if (abs(size - target_n) <= tolerance) {
      *actual_size = size;
      // Relabel for labeled case
      if (is_labeled && size > 0) {
        char *relabeled = relabel_object(obj, size, state);
        free(obj);
        // Canonicalize: sort Set children by min-label
        Object *parsed = parse_object(relabeled);
        free(relabeled);
        if (parsed) {
          canonicalize_sets(parsed);
          char *canonical = object_to_string(parsed);
          free_object(parsed);
          oracle_free(orc);
          return canonical;
        }
        oracle_free(orc);
        return strdup("ErrorCanonicalParse");
      }
      oracle_free(orc);
      return obj;
    }

    free(obj);
  }

  oracle_free(orc);
  return NULL; // Failed after max_attempts
}
