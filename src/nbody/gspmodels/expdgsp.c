/*
 * expdgsp.c: generate profile tables for spherical configuration with
 * same cumulative mass profile as an exponential disk.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <assert.h>

#define NWKSP  1000				// integration workspace size
#define ERRTOL 1.0e-7				// integration relative tol

local double msphr(double, double, double, double,
		   gsl_integration_workspace *);
local double dsphr(double, double, double, double,
		   gsl_integration_workspace *);
local double dsint(double z, void *params);
local double msint(double z, void *params);

//  gsp_expd: initialize tables for "spherical exponential disk".
//  _____________________________________________________________

gsprof *gsp_expd(double mtot, double alpha, double zdisk,
		 int np, double rmin, double rmax)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  double lgrs = log2(rmax / rmin) / (np - 1), x;
  gsl_integration_workspace *wksp;

  assert(np > 0 && mtot > 0.0 && alpha > 0.0 && zdisk >= 0);
  gsp->npoint = np;
  gsp->radius  = (double *) allocate(np * sizeof(double));
  gsp->density = (double *) allocate(np * sizeof(double));
  gsp->mass    = (double *) allocate(np * sizeof(double));
  if (zdisk > 0)				// will need integration?
    wksp = gsl_integration_workspace_alloc(NWKSP);
  for (int i = 0; i < np; i++) {
    gsp->radius[i] = rmin * exp2(lgrs * i);
    x = - alpha * gsp->radius[i];
    gsp->density[i] = (zdisk > 0 ?
		       dsphr(gsp->radius[i], mtot, alpha, zdisk, wksp) :
		       mtot * exp(x) * alpha*alpha / (4*M_PI*gsp->radius[i]));
    gsp->mass[i]    = (zdisk > 0 ?
		       msphr(gsp->radius[i], mtot, alpha, zdisk, wksp) :
		       mtot * (x * exp(x) - gsl_expm1(x)));
  }
  if (zdisk > 0)
    gsl_integration_workspace_free(wksp);
  for (int i = 1; i < np; i++) {		// check strict monotonicity
    assert(gsp->mass[i-1] < gsp->mass[i]);	// OK if alpha*rmax <= ~32
    assert(gsp->density[i-1] > gsp->density[i]);
  }
  gsp->alpha = (zdisk == 0 ? -1.0 : 0.0);
  gsp->beta = - (1 + alpha * rmax);
  gsp->mtot = mtot;
  return gsp;
}

//  dsphr: compute average density at radius r.  Value is
//  rhoDiskAvg = (4 Pi r^2)^-1 D[massDisk[r], r], 
//  _____________________________________________________

local double dsphr(double r, double mtot, double alpha, double zdisk,
		 gsl_integration_workspace *wksp)
{
  double params[3] = { alpha, zdisk, r }, result, abserr;
  gsl_function dsint_func = { &dsint, params };

  assert(gsl_integration_qag(&dsint_func, 0.0, r,
			     0.0, ERRTOL, NWKSP, GSL_INTEG_GAUSS61,
			     wksp, &result, &abserr) == GSL_SUCCESS);
  return (alpha*alpha * mtot * result / (4 * M_PI * zdisk * r));
}

//  dsint: integrand for dsphr() function.  Value is
//  Sech[z/z0]^2 Exp[-alpha Sqrt[r^2 - z^2]]
//  ________________________________________________

local double dsint(double z, void *params)
{
  double alpha = ((double *) params)[0];
  double zdisk = ((double *) params)[1];
  double r = ((double *) params)[2];

  assert(r*r - z*z >= 0);
  return (exp(- alpha * sqrt(r*r - z*z)) / gsl_pow_2(cosh(z / zdisk)));
}

//  msphr: compute mass interior to radius r by integrating in z.  Value is
//  massDisk = Integrate[2 Pi R rhoDisk[R,z], {z,-r,r}, {R,0,Sqrt[r^2 - z^2]}]
//  __________________________________________________________________________

local double msphr(double r, double mtot, double alpha, double zdisk,
		 gsl_integration_workspace *wksp)
{
  double params[3] = { alpha, zdisk, r }, result, abserr;
  gsl_function msint_func = { &msint, params };

  assert(gsl_integration_qag(&msint_func, 0.0, (double) r,
			     0.0, ERRTOL, NWKSP, GSL_INTEG_GAUSS61,
			     wksp, &result, &abserr) == GSL_SUCCESS);
  return (mtot * result / zdisk);
}

//  msint: integrand for msphr() function.  Value is
//  Sech[z/z0]^2 (1 - (1 + alpha Sqrt[r^2 - z^2]) Exp[-alpha Sqrt[r^2 - z^2]])
//  __________________________________________________________________________

local double msint(double z, void *params)
{
  double alpha = ((double *) params)[0];
  double zdisk = ((double *) params)[1];
  double r = ((double *) params)[2];

  assert(r*r - z*z >= 0);
  return ((1 - (1 + alpha * sqrt(r*r - z*z)) / exp(alpha * sqrt(r*r - z*z))) /
	    gsl_pow_2(cosh(z / zdisk)));
}

#ifdef UTILITY

string defv[] = {		";Generate \"spherical\" exponential disk",
  "out=???",			";Output file for profile tables",
  "mtot=1.0",			";Total mass of disk model",
  "alpha=1.0",			";Inverse length scale of disk",
  "zdisk=0.0",			";Vertical scale height of disk",
  "npoint=1025",		";Number of points in tables",
  "rrange=1/4096:16",		";Range of radii tabulated",
  "smartrange=true",		";Insure alpha*rrange[1] < ~32",
  "VERSION=2.0",		";Josh Barnes  9 July 2017",
  NULL,
};

int main(int argc, string argv[])
{
  real rrange[2], rmax;
  int np;
  double lgrs;
  gsprof *gsp;
  stream ostr;

  initparam(argv, defv);
  setrange(rrange, getparam("rrange"));
  np = getiparam("npoint");
  rmax = pow(2.0, floor(log2(32.0 / getdparam("alpha"))));
  if (rmax < rrange[1] && getbparam("smartrange")) {
    lgrs = log2(rrange[1] / rrange[0]) / (np - 1);
    np = 1 + log2(rmax / rrange[0]) / lgrs;
    eprintf("[%s: warning: npoint = %d -> %d  rrange[1] = %f -> %f]\n",
	    getprog(), getiparam("npoint"), np, rrange[1], rmax);
    rrange[1] = rmax;
  }
  gsp = gsp_expd(getdparam("mtot"), getdparam("alpha"), getdparam("zdisk"),
		 np, rrange[0], rrange[1]);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  gsp_write(ostr, gsp);
  fflush(NULL);
  return 0;
}

#endif
