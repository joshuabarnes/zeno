/*
 * GAMMP.C: incomplete gamma function P(a,x) [Numerical Recipes, Ch. 6.1].
 * Equivalent to Mathematica's Gamma[a, 0, x] .
 */

#include <stdio.h>
#include <math.h>

#ifdef TESTBED

#include <gsl/gsl_sf_gamma.h>

int main(int argc, char *argv[])
{
  double atof(char *);
  float gammp(float, float);

  if (argc == 3)
    printf("gammp(%f, %f) = %f (%f)\n", atof(argv[1]), atof(argv[2]),
	   gammp(atof(argv[1]), atof(argv[2])),
	   gsl_sf_gamma_inc_P(atof(argv[1]), atof(argv[2])));
  else
    printf("Usage: %s a x\n", argv[0]);
  return (0);
}

#endif

float gammp(float a, float x)
{
  float gamser, gamcnf, gln;
  void gser(float *, float, float, float *);
  void gcnf(float *, float, float, float *);

  if (x < 0.0 || a <= 0.0)
    nrerror("gammp: invalid arguments");
  if (x < (a + 1)) {
    gser(&gamser, a, x, &gln);
    return (gamser);
  } else {
    gcnf(&gamcnf, a, x, &gln);
    return (1.0 - gamcnf);
  }
}

#define ITMAX  100
#define EPS    3.0e-7

void gser(float *gamser, float a, float x, float *gln)
{
  float sum, del, ap, gammln(float);
  int n;

  *gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0)
      nrerror("gser: x less than 0");
    *gamser = 0.0;
  } else {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; fabs(del) > EPS * fabs(sum); n++) {
      if (n > ITMAX)
	nrerror("gser: a too large for convergence");
      ap += 1.0;
      del *= x / ap;
      sum += del;
    }
    *gamser = sum * exp(a * log(x) - *gln - x);
  }
}

void gcnf(float *gamcnf, float a, float x, float *gln)
{
  float g, anf, ana, an, a1, gammln(float);
  float gold = 0.0, a0 = 1.0, b0 = 0.0, b1 = 1.0, fac = 1.0;
  int n;

  *gln = gammln(a);
  a1 = x;
  for (n = 1; ; n++) {
    if (n > ITMAX)
      nrerror("gcnf: a too large for convergence");
    an = (float) n;
    ana = an - a;
    a0 = (a1 + a0 * ana) * fac;
    b0 = (b1 + b0 * ana) * fac;
    anf = an * fac;
    a1 = x * a0 + anf * a1;
    b1 = x * b0 + anf * b1;
    if (a1 != 0.0) {
      fac = 1.0 / a1;
      g = b1 * fac;
      if (fabs((g - gold) / g) < EPS)
	break;
      gold = g;
    }
  }
  *gamcnf = g * exp(a * log(x)- *gln - x);
}
