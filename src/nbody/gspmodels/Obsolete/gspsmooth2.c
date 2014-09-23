/*
 * GSPSMOOTH.C: smooth mass distribution of general spherical profile.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "assert.h"
#include "gsp.h"

/*
 * GSPSMOOTH1: old-style approximation, apparently based on power-law
 * model for mass profile, M(r) \propto r^(3 + alpha) .
 */

gsprof *gspsmooth1(gsprof *gsp, real eps, real kappa, string trace)
{
  gsprof *sgsp;
  int i;
  real c0, r, f, d;

  assert(eps > 0.0 && kappa > 0.0 && gsp->alpha < 0.0);
  if (trace != NULL)
    eprintf("[%s:  rho(eps) = %f  r[0]/eps = %f]\n", trace,
	    rho_gsp(gsp, eps), gsp->radius[0] / eps);
  sgsp = (gsprof *) allocate(sizeof(gsprof));
  sgsp->npoint = gsp->npoint;
  sgsp->radius = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->density = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->mass = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->alpha = 0.0;
  sgsp->beta = gsp->beta;
  sgsp->mtot = gsp->mtot;
  c0 = rpow(2.0/3.0, kappa/gsp->alpha);
  for (i = 0; i < gsp->npoint; i++) {
    r = gsp->radius[i];
    f = rpow(1 + c0 * rpow(eps/r, kappa), gsp->alpha/kappa);
    d = - gsp->alpha * c0 * (rpow(eps/r, kappa) / r) *
          rpow(1 + c0 * rpow(eps/r, kappa), gsp->alpha/kappa - 1);
    sgsp->radius[i] = gsp->radius[i];
    sgsp->density[i] = f * gsp->density[i] +
                        d * gsp->mass[i] / (4 * PI * rsqr(r));
    sgsp->mass[i] = f * gsp->mass[i];
  }
  return (sgsp);
}

/*
 * GSPSMOOTH2: new version based on mass-interpolation formula,
 * from Barnes (2011).
 */

gsprof *gspsmooth2(gsprof *gsp, real rho_0, real eta, string trace)
{
  gsprof *sgsp;
  int i;
  real r_i, mass_0, mass_r, rho_r;

  assert(rho_0 > 0.0 && eta > 0.0 && gsp->alpha <= 0.0);
  if (trace != NULL)
    eprintf("[%s:  rho_0 = %f  eta = %f]\n", trace, rho_0, eta);
  sgsp = (gsprof *) allocate(sizeof(gsprof));
  sgsp->npoint = gsp->npoint;
  sgsp->radius = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->density = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->mass = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->alpha = 0.0;
  sgsp->beta = gsp->beta;
  sgsp->mtot = gsp->mtot;
  for (i = 0; i < gsp->npoint; i++) {
    r_i = gsp->radius[i];
    mass_0 = (4 * PI / 3.0) * rqbe(r_i) * rho_0;
    mass_r = rpow(rpow(gsp->mass[i], -eta) + rpow(mass_0, -eta), -1/eta);
    rho_r = rpow(mass_r / gsp->mass[i], eta+1) * gsp->density[i] +
            rpow(mass_r / mass_0, eta+1) * rho_0;
    sgsp->radius[i] = r_i;
    sgsp->density[i] = rho_r;
    sgsp->mass[i] = mass_r;
  }
  return (sgsp);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Smooth mass distribution of GSP",
    "in=???",			";Input file with GSP",
    "out=",			";Output file for GSP",
    "eps=0.025",		";Plummer smoothing length",
    "kappa=",			";Old interpolation parameter",
    "eta=",			";New interpolation parameter",
    "VERSION=2",		";Josh Barnes  11 June 2011",
    NULL,
};

int main(int argc, string argv[])
{
  stream istr, ostr;
  gsprof *gsp, *sgsp;
  real rho_0;

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  gsp = get_gsprof(istr);
  assert(gsp->alpha < 0);
  if (! strnull(getparam("eta"))) {
    rho_0 = - gsp->alpha * rho_gsp(gsp, getdparam("eps"));
					/* use semi-bogus expr for rho(0) */
    sgsp = gspsmooth2(gsp, rho_0, getdparam("eta"), getargv0());
  } else if (! strnull(getparam("kappa"))) {
    sgsp = gspsmooth1(gsp, getdparam("eps"), getdparam("kappa"), getargv0());
  } else
    error("%s: one of eta or kappa must be specified\n", getargv0());
  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_gsprof(ostr, sgsp);
    strclose(ostr);
  }
  return (0);
}

#endif
