/*
 * gspgravity.c: calculate gravitational potential of GSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"
#include <assert.h>

#define NITER  3				// iterations for gsp_phi_rad

//  gsp_phi: evaluate potential at given radius.
//  ____________________________________________

double gsp_phi(gsprof *gsp, double r)
{
  int N = gsp->npoint - 1, xind;
  double r0 = gsp->radius[0], m0 = gsp->mass[0], phi0;
  double rN = gsp->radius[N], mN = gsp->mass[N], lgr;

  if (gsp->phi == NULL)
    gsp_calc_phi(gsp);
  if (r < 0)
    error("%s.gsp_phi: undefined for r = %e < 0\n", getprog(), r);
  phi0 = gsp->phi[0];
  if (r < r0)
    return (phi0 - (m0 / r0) *
	    (gsp->alpha == -2 ? rlog(r0/r) :
	     (1 - pow(r/r0, 2 + gsp->alpha)) / (2 + gsp->alpha)));
  else if (r > rN)
    return (- (gsp->mtot + (gsp->mtot - mN) *
	       pow(r/rN, 3 + gsp->beta) / (2 + gsp->beta)) / r);
  if (gsp->pr_spline == NULL) {
    gsp->pr_spline = gsl_interp_alloc(gsl_interp_akima, gsp->npoint);
#ifdef USELOG2R
    gsl_interp_init(gsp->pr_spline, gsp->lg2rad, gsp->phi, gsp->npoint);
#else
    gsl_interp_init(gsp->pr_spline, gsp->radius, gsp->phi, gsp->npoint);
#endif
  }
#ifdef USELOG2R
  lgr = log2(r);
  xind = N * (lgr - gsp->lg2rad[0]) / (gsp->lg2rad[N] - gsp->lg2rad[0]);
  gsp->r_acc->cache = MAX(0, MIN(N - 1, xind));
  return gsl_interp_eval(gsp->pr_spline, gsp->lg2rad, gsp->phi,
			 lgr, gsp->r_acc);
#else
  return gsl_interp_eval(gsp->pr_spline, gsp->radius, gsp->phi,
			 r, gsp->r_acc);
#endif  
}

//  gsp_phi_rad: evaluate radius corresponding to given potential.
//  ______________________________________________________________

double gsp_phi_rad(gsprof *gsp, double p)
{
  int N = gsp->npoint - 1;
  double r0 = gsp->radius[0], m0 = gsp->mass[0], phi0;
  double rN = gsp->radius[N], mN = gsp->mass[N], phiN, x, r;

  if (gsp->phi == NULL)
    gsp_calc_phi(gsp);
  if (p > 0)
    error("%s.gsp_phi_rad: undefined for phi = %e > 0\n", getprog(), p);
  phi0 = gsp->phi[0];
  phiN = gsp->phi[N];
  if (p < phi0) {
    if (gsp->alpha > -2 && p < phi0 - m0 / (r0*(2 + gsp->alpha)))
      error("%s.gsp_phi_rad: phi = %e < phi(0)\n", getprog(), p);
    x = (r0 / m0) * (p - phi0);
    return (gsp->alpha == -2 ? r0 * exp(x) :
	    r0 * pow(1 + (2 + gsp->alpha) * x, 1 / (2 + gsp->alpha)));
  } else if (p > phiN) {
    r = - gsp->mtot / p;			// make initial guess
    for (int i = 0; i < NITER; i++)
      r = - (gsp->mtot + pow(r / rN, 3 + gsp->beta) *
	     (gsp->mtot - mN) / (2 + gsp->beta)) / p;
    return r;					// return refined estimate
  }
  if (gsp->rp_spline == NULL) {
    gsp->rp_spline = gsl_interp_alloc(gsl_interp_akima, gsp->npoint);
    gsl_interp_init(gsp->rp_spline, gsp->phi, gsp->radius, gsp->npoint);
  }
  return gsl_interp_eval(gsp->rp_spline, gsp->phi, gsp->radius, p, NULL);
}

//  gsp_calc_phi: calculate gravitational potential of the given profile.
//  _____________________________________________________________________

void gsp_calc_phi(gsprof *gsp)
{
  double *acc = (double *) allocate(gsp->npoint * sizeof(double));
  double r, phi, delphi;
  gsl_interp *ar_spline = gsl_interp_alloc(gsl_interp_akima, gsp->npoint);

  for (int i = 0; i < gsp->npoint; i++) {	// tabulate acceleration
    r = gsp->radius[i];
    acc[i] = - gsp_mass(gsp, r) / (r > 0 ? r*r : 1.0);
  }
  gsl_interp_init(ar_spline, gsp->radius, acc, gsp->npoint);
  gsp->phi = (double *) allocate(gsp->npoint * sizeof(double));
  phi = - (gsp->mtot + (gsp->mtot - gsp->mass[gsp->npoint - 1]) /
	     (2 + gsp->beta)) / gsp->radius[gsp->npoint - 1];
  gsp->phi[gsp->npoint - 1] = phi;		// store outermost potential
  for (int i = gsp->npoint - 2; i >= 0; i--) {	// integrate acceleration
    delphi = gsl_interp_eval_integ(ar_spline, gsp->radius, acc,
				   gsp->radius[i], gsp->radius[i+1], NULL);
    phi += delphi;				// sum potentail differences
    gsp->phi[i] = phi;				// store potential at radius
    if (i > 0 && gsp->phi[i] == gsp->phi[i-1])
      error("%s.gsp_calc_phi: potential degenerate (i = %d)\n", getprog(), i);
  }
  free(acc);
  gsl_interp_free(ar_spline);
}

#ifdef UTILITY

string defv[] = {		";Tabulate gravitational potential",
  "gsp=???",			";GSP for mass model",
  "npoint=161",			";Number of potential values to list",
  "rrange=1/1024:1024",		";Range of radii tabulated",
  "VERSION=2.0",		";Josh Barnes  7 June 2017",
  NULL,
};

int main(int argc, string argv[])
{
  stream str;
  gsprof *gsp;
  int npnt;
  real rrange[2];
  double lgrs, r;

  initparam(argv, defv);
  str = stropen(getparam("gsp"), "r");
  get_history(str);
  gsp = gsp_read(str);
  gsp_calc_phi(gsp);				// perform calculation
  npnt = getiparam("npoint");
  setrange(rrange, getparam("rrange"));
  lgrs = log2(rrange[1] / rrange[0]) / (npnt - 1);
  eprintf("[%s: lgrs = %f (%a)]\n", getprog(), lgrs, lgrs);
  printf("#%5s  %16s  %16s  %16s\n", "i", "r", "phi", "acc");
  for (int i = 0; i < npnt; i++) {
    r = rrange[0] * pow(2.0, lgrs * i);
    printf("%6d  %16.8e  %16.8e  %16.8e\n", i, r,
	   gsp_phi(gsp, r), - gsp_mass(gsp, r) / (r*r));
  }
  gsp_free(gsp);
  return 0;
}

#endif

#ifdef TESTBED

string defv[] = {		";Calculate potential of GSP",
  "gsp=???",			";Input file with GSP",
  "npoint=27",			";Number of points tabulated",
  "rrange=1/128:64",		";Range of radii tabulated",
  "VERSION=2.0",		";Josh Barnes  15 May 2017",
  NULL,
};

int main(int argc, string argv[])
{
  stream istr;
  gsprof *gsp;
  int np;
  real rrange[2];
  double lgrs, r;

  initparam(argv, defv);
  istr = stropen(getparam("gsp"), "r");
  get_history(istr);
  gsp = gsp_read(istr);
  (void) gsp_phi(gsp, 1.0);			// precalculate phi table
  np = getiparam("npoint");
  setrange(rrange, getparam("rrange"));
  lgrs = log2(rrange[1] / rrange[0]) / (np - 1);
  printf("#%15s  %16s  %16s\n", "radius", "phi(r)", "1-r(phi)/r");
  for (int i = 0; i < np; i++) {
    r = rrange[0] * pow(2.0, lgrs * i);
    printf("%16.8e  %16.8e  %16.8e\n",
	   r, gsp_phi(gsp, r), 1 - gsp_phi_rad(gsp, gsp_phi(gsp, r)) / r);
  }
  return (0);
}

#endif
