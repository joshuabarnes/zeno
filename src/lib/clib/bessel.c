/*
 * bessel.c: interfaces to gsl bessel functions provide error checking.
 */

#include "stdinc.h"
#include "getparam.h"

#include <gsl/gsl_sf_bessel.h>

//  bessel_I0: regular modified cylindrical Bessel function I, order zero.
//  ______________________________________________________________________

double bessel_I0(double x)
{
  gsl_sf_result res;
  int stat;

  stat = gsl_sf_bessel_I0_e(x, &res);
  if (stat != 0)
    error("%s.bessel_I0: error status: %s\n", getprog(), gsl_strerror(stat));
  return (res.val);
}

//  bessel_I1: regular modified cylindrical Bessel function I, order one.
//  _____________________________________________________________________

double bessel_I1(double x)
{
  gsl_sf_result res;
  int stat;

  stat = gsl_sf_bessel_I1_e(x, &res);
  if (stat != 0)
    error("%s.bessel_I1: error status: %s\n", getprog(), gsl_strerror(stat));
  return (res.val);
}

//  bessel_I: regular modified cylindrical Bessel function I, order nu.
//  ___________________________________________________________________

double bessel_I(double nu, double x)
{
  gsl_sf_result res;
  int stat;

  stat = gsl_sf_bessel_Inu_e(nu, x, &res);
  if (stat != 0)
    error("%s.bessel_I: error status: %s\n", getprog(), gsl_strerror(stat));
  return (res.val);
}

//  bessel_K0: irregular modified cylindrical Bessel function K, order zero.
//  ________________________________________________________________________

double bessel_K0(double x)
{
  gsl_sf_result res;
  int stat;

  stat = gsl_sf_bessel_K0_e(x, &res);
  if (stat != 0)
    error("%s.bessel_K0: error status: %s\n", getprog(), gsl_strerror(stat));
  return (res.val);
}

//  bessel_K1: irregular modified cylindrical Bessel function K, order one.
//  _______________________________________________________________________

double bessel_K1(double x)
{
  gsl_sf_result res;
  int stat;

  stat = gsl_sf_bessel_K1_e(x, &res);
  if (stat != 0)
    error("%s.bessel_K1: error status: %s\n", getprog(), gsl_strerror(stat));
  return (res.val);
}

//  bessel_K: irregular modified cylindrical Bessel function K, order nu.
//  _____________________________________________________________________

double bessel_K(double nu, double x)
{
  gsl_sf_result res;
  int stat;

  stat = gsl_sf_bessel_Knu_e(nu, x, &res);
  if (stat != 0)
    error("%s.bessel_k: error status: %s\n", getprog(), gsl_strerror(stat));
  return (res.val);
}

#if defined(TESTBED)

#include "mathfns.h"

string defv[] = {		";Test bessel functions",
  "x=4.5",			";Argument for functions",
  "nu=0.5",			";Order for functions",
  "VERSION=1.0",		";Joshua Barnes  11 May 2012",    
  NULL,
};

int main(int argc, string argv[])
{
  double x, nu, ex;

  initparam(argv, defv);
  x = getdparam("x");
  nu = getdparam("nu");
  ex = exp(x);
  printf("bessel_I0(%f) = %f\t%f\n", x, bessel_I0(x), bessel_I0(x) / ex);
  printf("bessel_I1(%f) = %f\t%f\n", x, bessel_I1(x), bessel_I1(x) / ex);
  printf("bessel_K0(%f) = %f\t%f\n", x, bessel_K0(x), ex * bessel_K0(x));
  printf("bessel_K1(%f) = %f\t%f\n", x, bessel_K1(x), ex * bessel_K1(x));
  printf("bessel_I(%f,%f) = %f\t%f\n",
	 nu, x, bessel_I(nu, x), bessel_I(nu, x) / ex);
  printf("bessel_K(%f,%f) = %f\t%f\n",
	 nu, x, bessel_K(nu, x), bessel_K(nu, x) / ex);
  return (0);
}

#endif
