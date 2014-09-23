/*
 * XRANDOM.C: useful functions for random numbers.
 */

#include "stdinc.h"
#include "mathfns.h"

#if (!defined(LINUX) && !defined(MACOSX))
long random(void);
#endif

/*
 * XRANDOM: floating-point random number routine.
 */

double xrandom(double xl, double xh)
{

    return (xl + (xh - xl) * ((double) random()) / 2147483647.0);
}

/*
 * GRANDOM: normally distributed random number (polar method).
 * Reference: Knuth, vol. 2, p. 104.
 */

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

#if defined(TESTBED)

#include "getparam.h"

string defv[] = {
  "seed=123",
  "count=10",
  NULL,
};

int main(int argc, string argv[])
{
  int n;
  double x;

  initparam(argv, defv);
  srandom(getiparam("seed"));
  for (n = getiparam("count"); n > 0; n--) {
    x = xrandom(0.0, 1.0);
    printf("xrandom(0,1) -> %.15f (%a)\n", x, x);
  }
  return (0);
}

#endif
