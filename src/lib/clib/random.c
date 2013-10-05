/*
 * random.c: useful functions for random numbers.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include <string.h>
#include <gsl/gsl_rng.h>

local gsl_rng *rng = NULL;
local unsigned long rng_min, rng_max;

//  __________________________________________________________________
//  init_random: select and initialize random number generator.  For
//  backward compatibility, default reproduces unix random() function;
//  set GSL_RNG_TYPE environment variable to select better generators.

void init_random(unsigned long seed)
{
  if (rng == NULL) {				// select on first call
    if (getenv("GSL_RNG_TYPE") != NULL && !strnull(getenv("GSL_RNG_TYPE"))) {
      rng = gsl_rng_alloc(gsl_rng_env_setup());
      eprintf("[%s.init_random: generator %s, range %lu:%lu]\n", getprog(),
	      gsl_rng_name(rng), gsl_rng_min(rng), gsl_rng_max(rng));
    } else
      rng = gsl_rng_alloc(gsl_rng_random128_glibc2);
  }
  rng_min = gsl_rng_min(rng);
  rng_max = gsl_rng_max(rng);
  gsl_rng_set(rng, seed);
}

//  __________________________________________________________________
//  xrandom: floating-point random number in range [xl,xh], inclusive;
//  computed by scaling unsigned long, so at most 32 bits are random.

double xrandom(double xl, double xh)
{
  if (rng == NULL) {
    eprintf("[%s.xrandom: WARNING: using default seed]\n", getprog());
    init_random(0);
  }
  return (xl + (xh - xl) *
	    (gsl_rng_get(rng) - rng_min) / ((double) (rng_max - rng_min)));
}

//  ___________________________________________________________
//  grandom: normally distributed random number (polar method).
//  Reference: Knuth, vol. 2, p. 104.

double grandom(double mean, double sdev)
{
  double v1, v2, s;

  do {
    v1 = xrandom(-1.0, 1.0);
    v2 = xrandom(-1.0, 1.0);
    s = v1*v1 + v2*v2;
  } while (s >= 1.0);
  return (mean + sdev * v1 * sqrt(-2.0 * log(s) / s));
}

//  _________________________________________________________________
//  get_random_state: get snapshot of current random generator state.
//  State pointer *st must be initialized to NULL before first call.

void get_random_state(int *nb, void **st)
{
  if (*st == NULL) {
    *nb = (int) gsl_rng_size(rng);
    *st = (void *) allocate(*nb);
  } else if (*nb != (int) gsl_rng_size(rng))
    error("%s.get_random_state: inconsistent state size\n", getprog());
  memcpy(*st, gsl_rng_state(rng), *nb);
}

//  _______________________________________________________________
//  set_random_state: copy state snapshot back to random generator.

void set_random_state(int *nb, void **st)
{
  if (rng == NULL)
    init_random(0);				// make sure rng exists
  if (*nb != (int) gsl_rng_size(rng))
    error("%s.set_random_state: inconsistent state size\n", getprog());
  memcpy(gsl_rng_state(rng), *st, *nb);
}

#if defined(TESTBED)

string defv[] = {
  "seed=123",
  "count=10",
  "state=false",
  NULL,
};

int main(int argc, string argv[])
{
  int n, nb;
  double x;
  void *st = NULL;

  initparam(argv, defv);
  init_random(getiparam("seed"));
  for (n = getiparam("count"); n > 0; n--) {
    x = xrandom(0.0, 1.0);
    printf("xrandom(0,1) -> %.15f (%a)\n", x, x);
  }
  if (getbparam("state")) {
    get_random_state(&nb, &st);
    printf("state size: %d bytes\n", nb);
    printf("state: ");
    for (n = 0; n < nb; n++)
      printf("%02x%c", ((byte *) st)[n], n < nb-1 ? ' ' : '\n');
    for (n = getiparam("count"); n > 0; n--) {
      x = xrandom(0.0, 1.0);
      printf("xrandom(0,1) -> %.15f (%a)\n", x, x);
    }
    printf("resetting random state\n");
    set_random_state(&nb, &st);
    for (n = getiparam("count"); n > 0; n--) {
      x = xrandom(0.0, 1.0);
      printf("xrandom(0,1) -> %.15f (%a)\n", x, x);
    }
  }
  return (0);
}

#endif
