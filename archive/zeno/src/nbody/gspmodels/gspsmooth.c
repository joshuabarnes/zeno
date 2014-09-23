/*
 * GSPSMOOTH.C: smooth mass distribution of general spherical profile.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "assert.h"
#include "gsp.h"

gsprof *gspsmooth(gsprof *gsp, real eps, real kappa, string trace)
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
#ifdef DEBUG
  printf("#%11s %12s %12s %12s %12s %12s\n",
	 "r", "f", "d", "rho1", "rho2", "mass");
#endif
  for (i = 0; i < gsp->npoint; i++) {
    r = gsp->radius[i];
    f = rpow(1 + c0 * rpow(eps/r, kappa), gsp->alpha/kappa);
    d = - gsp->alpha * c0 * (rpow(eps/r, kappa) / r) *
          rpow(1 + c0 * rpow(eps/r, kappa), gsp->alpha/kappa - 1);
    sgsp->radius[i] = gsp->radius[i];
    sgsp->density[i] = f * gsp->density[i] +
                        d * gsp->mass[i] / (4 * PI * rsqr(r));
    sgsp->mass[i] = f * gsp->mass[i];
#ifdef DEBUG
    printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	   r, f, d, f * gsp->density[i],
	   d * gsp->mass[i] / (4 * PI * rsqr(r)), f * gsp->mass[i]);
#endif
  }
  return (sgsp);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Smooth mass distribution of GSP",
    "in=???",			";Input file with GSP",
    "out=???",			";Output file for GSP",
    "eps=0.025",		";Plummer smoothing length",
    "kappa=1.5",		";Interpolation parameter",
    "VERSION=1.2",		";Josh Barnes  11 July 2010",
    NULL,
};

int main(int argc, string argv[])
{
  stream istr, ostr;
  gsprof *gsp, *sgsp;

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  gsp = get_gsprof(istr);
  sgsp = gspsmooth(gsp, getdparam("eps"), getdparam("kappa"), getargv0());
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  put_gsprof(ostr, sgsp);
  strclose(ostr);
  return (0);
}

#endif
