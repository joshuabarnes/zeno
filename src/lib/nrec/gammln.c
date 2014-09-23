/*
 * GAMMLN.C: log of gamma function [Numerical Recipes, Ch. 6.1].
 */

#include <stdio.h>
#include <math.h>

#ifdef TESTBED

#include <gsl/gsl_sf_gamma.h>

int main(int argc, char *argv[])
{
  double atof(char *);
  float gammln(float xx);

  if (argc == 2)
    printf("gammln(%f) = %f (%f)\n", atof(argv[1]),
	   gammln(atof(argv[1])), gsl_sf_lngamma(atof(argv[1])));
  else
    printf("Usage: %s <arg>\n", argv[0]);
  return (0);
}

#endif

float gammln(float xx)
{
  double x, tmp, ser;
  static double cof[6] =
    {  76.18009173, -86.50532033,     24.01409822,
       -1.231739516,  0.120858003e-2, -0.536382e-5 };
  int j;

  x = xx - 1.0;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.0;
  for (j = 0; j <= 5; j++) {
    x = x + 1.0;
    ser += cof[j]/x;
  }
  return (log(2.50662827465 * ser) - tmp);
}
