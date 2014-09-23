/*
 * expdgsp.c: generate profile tables for spherical configuration with
 * same cumulative mass profile as an exponential disk.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include "assert.h"

#ifdef UTILITY

string defv[] = {		";Generate \"spherical\" exp. disk",
    "out=???",			";Output file for profile tables",
    "mtot=1.0",			";Total mass of disk",
    "alpha=1.0",		";Inverse length scale of disk",
    "zdisk=0.0",		";Vertical scale height of disk",
    "npoint=257",		";Number of points in tables",
    "rrange=1/4096:16",		";Range of radii tabulated",
    "VERSION=1.1",		";Josh Barnes  26 January 2013",
    NULL,
};

int main(int argc, string argv[])
{
  real rrange[2];
  gsprof *gsp;
  stream ostr;

  initparam(argv, defv);
  setrange(rrange, getparam("rrange"));
  gsp = expdgsp(getdparam("mtot"), getdparam("alpha"), getdparam("zdisk"),
		getiparam("npoint"), rrange[0], rrange[1]);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  put_gsprof(ostr, gsp);
  strclose(ostr);
  return (0);
}

#endif

#define NWKSP  1000				// integration workspace size
#define ERRTOL 1.0e-7				// integration relative tol

local real msphr(real, real, real, real, gsl_integration_workspace *);
local real dsphr(real, real, real, real, gsl_integration_workspace *);
local double dsint(double z, void *params);
local double msint(double z, void *params);

//  ____________________________________________________________
//  expdgsp: initialize tables for "spherical exponential disk".

gsprof *expdgsp(real mtot, real alpha, real zdisk,
		int np, real rmin, real rmax)
{
  gsprof *gsp;
  real *rtab, *dtab, *mtab, lgrs, x;
  int i;
  gsl_integration_workspace *wksp;

  assert(mtot > 0.0 && alpha > 0.0 && zdisk >= 0);
  gsp = (gsprof *) allocate(sizeof(gsprof));
  rtab = (real *) allocate(np * sizeof(real));
  dtab = (real *) allocate(np * sizeof(real));
  mtab = (real *) allocate(np * sizeof(real));
  lgrs = rlog2(rmax / rmin) / (np - 1);
  if (zdisk == 0) {
    for (i = 0; i < np; i++) {
      rtab[i] = rmin * rexp2(lgrs * i);
      x = - alpha * rtab[i];
      dtab[i] = mtot * (rsqr(alpha) / (4 * PI * rtab[i])) * rexp(x);
      mtab[i] = mtot * (x * rexp(x) - gsl_expm1(x));	// use stable formula
    }
  } else {
    wksp = gsl_integration_workspace_alloc(NWKSP);
    for (i = 0; i < np; i++) {
      rtab[i] = rmin * rexp2(lgrs * i);
      dtab[i] = dsphr(rtab[i], mtot, alpha, zdisk, wksp);
      mtab[i] = msphr(rtab[i], mtot, alpha, zdisk, wksp);
    }
    gsl_integration_workspace_free(wksp);
  }
  for (i = 1; i < np; i++) {				// check monotonicity
    if (dtab[i-1] < dtab[i] || mtab[i-1] > mtab[i])
      error("%s: %s not monotonic near r = %f\n", 
	    getprog(), mtab[i-1] > mtab[i] ? "mass" : "density", rtab[i]);
    if (dtab[i-1] <= dtab[i] || mtab[i-1] >= mtab[i]) {
      eprintf("[%s.expdgsp: warning: %s not strictly monotonic near r = %f]\n",
	      getprog(), mtab[i-1] >= mtab[i] ? "mass" : "density", rtab[i]);
      break;
    }
  }
  gsp->npoint = np;
  gsp->radius = rtab;
  gsp->density = dtab;
  gsp->mass = mtab;
  gsp->alpha = (zdisk == 0 ? -1.0 : 0.0);
  gsp->beta = - (1 + alpha * rmax);
  gsp->mtot = mtot;
  return (gsp);
}

//  _____________________________________________________
//  dsphr: compute average density at radius r.  Value is
//  rhoDiskAvg = (4 Pi r^2)^-1 D[massDisk[r], r], 

local real dsphr(real r, real mtot, real alpha, real zdisk,
		 gsl_integration_workspace *wksp)
{
  double params[3], result, abserr;
  gsl_function dsint_func;
  int stat;

  params[0] = alpha;
  params[1] = zdisk;
  params[2] = r;
  dsint_func.function = &dsint;
  dsint_func.params = params;
  stat = gsl_integration_qag(&dsint_func, 0.0, (double) r,
			     0.0, ERRTOL, NWKSP, GSL_INTEG_GAUSS61,
			     wksp, &result, &abserr);
  assert(stat == 0);
  return (rsqr(alpha) * mtot * result / (4 * M_PI * zdisk * r));
}

//  ________________________________________________
//  dsint: integrand for dsphr() function.  Value is
//  Sech[z/z0]^2 Exp[-alpha Sqrt[r^2 - z^2]]

local double dsint(double z, void *params)
{
  double alpha = ((double *) params)[0];
  double zdisk = ((double *) params)[1];
  double r = ((double *) params)[2];

  assert(r*r - z*z >= 0);
  return (exp(- alpha * sqrt(r*r - z*z)) / gsl_pow_2(cosh(z / zdisk)));
}

//  __________________________________________________________________________
//  msphr: compute mass interior to radius r by integrating in z.  Value is
//  massDisk = Integrate[2 Pi R rhoDisk[R,z], {z,-r,r}, {R,0,Sqrt[r^2 - z^2]}]

local real msphr(real r, real mtot, real alpha, real zdisk,
		 gsl_integration_workspace *wksp)
{
  double params[3], result, abserr;
  gsl_function msint_func;
  int stat;

  params[0] = alpha;
  params[1] = zdisk;
  params[2] = r;
  msint_func.function = &msint;
  msint_func.params = params;
  stat = gsl_integration_qag(&msint_func, 0.0, (double) r,
			     0.0, ERRTOL, NWKSP, GSL_INTEG_GAUSS61,
			     wksp, &result, &abserr);
  assert(stat == 0);
  return (mtot * result / zdisk);
}

//  __________________________________________________________________________
//  msint: integrand for msphr() function.  Value is
//  Sech[z/z0]^2 (1 - (1 + alpha Sqrt[r^2 - z^2]) Exp[-alpha Sqrt[r^2 - z^2]])

local double msint(double z, void *params)
{
  double alpha = ((double *) params)[0];
  double zdisk = ((double *) params)[1];
  double r = ((double *) params)[2];

  assert(r*r - z*z >= 0);
  return ((1 - (1 + alpha * sqrt(r*r - z*z)) / exp(alpha * sqrt(r*r - z*z))) /
	    gsl_pow_2(cosh(z / zdisk)));
}
