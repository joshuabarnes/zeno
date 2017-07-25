/*
 * gspsmooth.c: smooth mass distribution of general spherical profile.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "assert.h"
#include "gsp.h"

gsprof *gsp_smooth(gsprof *gsp, double eps, double kappa, string trace)
{
  gsprof *sgsp = (gsprof *) allocate(sizeof(gsprof));
  double c0, r, f, d;

  assert(eps > 0.0 && kappa > 0.0 && gsp->alpha < 0.0);
  if (trace != NULL)
    eprintf("[%s: rho(eps) = %f  r[0]/eps = %f]\n", trace,
	    gsp_rho(gsp, eps), gsp->radius[0] / eps);
  sgsp->npoint = gsp->npoint;
  sgsp->radius = (double *) allocate(sgsp->npoint * sizeof(double));
  sgsp->density = (double *) allocate(sgsp->npoint * sizeof(double));
  sgsp->mass = (double *) allocate(sgsp->npoint * sizeof(double));
  sgsp->alpha = 0.0;
  sgsp->beta = gsp->beta;
  sgsp->mtot = gsp->mtot;
  c0 = pow(2.0/3.0, kappa/gsp->alpha);
#ifdef DEBUG
  printf("#%11s %12s %12s %12s %12s %12s\n",
	 "r", "f", "d", "rho1", "rho2", "mass");
#endif
  for (int i = 0; i < gsp->npoint; i++) {
    r = gsp->radius[i];
    f = pow(1 + c0 * pow(eps/r, kappa), gsp->alpha/kappa);
    d = - (gsp->alpha * c0 * (pow(eps/r, kappa) / r) *
	   pow(1 + c0 * pow(eps/r, kappa), gsp->alpha/kappa - 1));
    sgsp->radius[i] = gsp->radius[i];
    sgsp->density[i] = (f * gsp->density[i] +
                        d * gsp->mass[i] / (4 * M_PI * r * r));
    sgsp->mass[i] = f * gsp->mass[i];
#ifdef DEBUG
    printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	   r, f, d, f * gsp->density[i],
	   d * gsp->mass[i] / (4 * M_PI * r * r), f * gsp->mass[i]);
#endif
  }
  return sgsp;
}

#ifdef UTILITY

#include "getparam.h"

string defv[] = {		";Smooth mass distribution of GSP",
  "in=???",			";Input file with GSP",
  "out=???",			";Output file for GSP",
  "eps=0.025",			";Plummer smoothing length",
  "kappa=1.5",			";Interpolation parameter",
  "VERSION=2.0",		";Josh Barnes  21 June 2017",
  NULL,
};

int main(int argc, string argv[])
{
  stream istr, ostr;
  gsprof *gsp, *sgsp;

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  gsp = gsp_read(istr);
  sgsp = gsp_smooth(gsp, getdparam("eps"), getdparam("kappa"), getprog());
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  gsp_write(ostr, sgsp);
  fflush(NULL);
  return 0;
}

#endif
