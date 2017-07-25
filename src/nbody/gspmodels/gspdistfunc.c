/*
 * gspdistfunc.c: calculate distribution function of GSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

//  gsp_dist: evaluate distribution function f = f(E) = dF/dE.
//  __________________________________________________________

double gsp_dist(gsprof *dgsp, double E)
{
  double E0, EN, FN, f = 0.0;
  static bool warn = FALSE;
  
  if (dgsp->df_spline == NULL)
    error("%s.gsp_dist: call gsp_calc_dist first\n", getprog());
  E0 = dgsp->energy[0];
  EN = dgsp->energy[dgsp->npoint-1];
  FN = dgsp->dfint[dgsp->npoint-1];
  if (E0 <= E && E <= EN) {			// in range: use spline
    f = gsl_interp_eval_deriv(dgsp->df_spline,
			      dgsp->energy, dgsp->dfint, E, NULL);
  } else if (EN < E && E < 0) {			// use asymptotic form
    f = - (dgsp->beta + (dgsp->raniso <= 0 ? 0.5 : 2.5)) * FN *
          pow(EN / E, dgsp->beta + (dgsp->raniso <= 0 ? 1.5 : 3.5)) / EN;
  } else if (E < E0) {				// danger: fake it!
    f = gsl_interp_eval_deriv(dgsp->df_spline,
			      dgsp->energy, dgsp->dfint, E0, NULL);
    eprintf("[%s.gsp_dist: WARNING: E = %.8e < E0 = %.8e;\n"
	    " returning endpoint value f(E0) = %.8e]\n", getprog(), E, E0, f);
  }
  if ((f < 0.0 && !warn) || E >= 0) {
    eprintf("[%s.gsp_dist: WARNING: f = %e%s for E = %e]\n",
	    getprog(), f, f < 0 ? " < 0" : "", E);
    warn = TRUE;
  }
  return f;
}

//  gsp_dist_integ: evaluate integral of distribution function F(E).
//  ________________________________________________________________

double gsp_dist_integ(gsprof *dgsp, double E)
{
  if (dgsp->df_spline == NULL)
    error("%s.gsp_dist_integ: call gsp_calc_dist first\n", getprog());
  if (dgsp->energy[dgsp->npoint-1] < E && E < 0)
    return (dgsp->dfint[dgsp->npoint-1] *
	    pow(dgsp->energy[dgsp->npoint-1] / E,
		dgsp->beta + (dgsp->raniso <= 0 ? 0.5 : 2.5)));
  return gsl_interp_eval(dgsp->df_spline, dgsp->energy, dgsp->dfint, E, NULL);
}

//  dFpars: parameter block used by integrand for F(E).
//  ___________________________________________________

typedef struct {
  gsprof *dgsp;				// GSP specifying density profile
  gsprof *ggsp;				// GSP specifying gravational field
  double E;				// energy value
  double ra;				// Osipkov-Merritt anisotropy radius
} dFpars;

//  dFinteg: integrand for F(E) or F(Q).
//  ____________________________________

local double dFinteg(double x, void *pars)
{
  gsprof *dgsp = ((dFpars *) pars)->dgsp;
  gsprof *ggsp = ((dFpars *) pars)->ggsp;
  double E =     ((dFpars *) pars)->E;
  double ra =    ((dFpars *) pars)->ra;
  double r, grad;

  r = gsp_phi_rad(ggsp, MIN(x*x + E, 0.0));	// radius for Phi = x*x+E
  grad = gsp_grad(dgsp, r);			// grad. for isotropic model
  if (ra > 0)					// got Osipkov-Merritt model?
    grad = (1 + (r/ra)*(r/ra)) * grad +  (2 * r / (ra*ra)) * gsp_rho(dgsp, r);
  return (grad * r*r / (M_SQRT2*M_PI*M_PI * gsp_mass(ggsp, r)));
}

// Constants and parameters for numerical integration.

#define NWKSP  1000			// integration workspace size
local double epsabs = 0.0;		// integration absolute tolerance
local double epsrel = 1.0e-4;		// integration relative tolerance
local bool usephi = TRUE;		// use phi values for energy points

//  gsp_calc_dist: calculate tables for distribution function.
//  __________________________________________________________

void gsp_calc_dist(gsprof *dgsp, gsprof *ggsp, double ra)
{
  int npnt = dgsp->npoint, stat;
  double Emin, Emax, abserr, avgerr, maxerr;
  gsl_integration_workspace *wksp = gsl_integration_workspace_alloc(NWKSP);
  gsl_function Finteg_func;
  dFpars params = { dgsp, ggsp, 0.0, ra };

  if (dgsp->energy != NULL && dgsp->dfint == NULL) {
    free(dgsp->energy);				// free old arrays, if any
    free(dgsp->dfint);
  }
  dgsp->energy = (double *) allocate(npnt * sizeof(double));
  dgsp->dfint = (double *) allocate(npnt * sizeof(double));
  (void) gsp_phi(ggsp, 1.0);			// precompute potential table
  if (usephi) {
    for (int i = 0; i < npnt; i++)
      dgsp->energy[i] = gsp_phi(ggsp, dgsp->radius[i]);
  } else {
    Emin = ggsp->phi[0];			// ABSTRACTION VIOLATION!!!
    Emax = ggsp->phi[ggsp->npoint-1];		// ggsp should be black box
    for (int i = 0; i < npnt; i++)
      dgsp->energy[i] =	Emin + i * (Emax - Emin) / (npnt-1.0);
  }
  Finteg_func.function = &dFinteg;
  Finteg_func.params = &params;
  avgerr = maxerr = 0.0;
  for (int i = 0; i < npnt; i++) {
    params.E = dgsp->energy[i];			// integrate 0 to sqrt(-E)
    assert(params.E < 0);
    stat = gsl_integration_qag(&Finteg_func, 0.0, sqrt(- params.E),
    			       epsabs, epsrel, NWKSP, GSL_INTEG_GAUSS61,
    			       wksp, &dgsp->dfint[i], &abserr);
    assert(stat == 0);
    avgerr += ABS(abserr / dgsp->dfint[i]) / npnt;
    maxerr = MAX(maxerr, ABS(abserr / dgsp->dfint[i]));
  }
  eprintf("[%s.gsp_calc_dist: relerr = %.8e (avg), %.8e (max)]\n",
	  getprog(), avgerr, maxerr);
  dgsp->df_spline = gsl_interp_alloc(gsl_interp_akima, npnt);
  gsl_interp_init(dgsp->df_spline, dgsp->energy, dgsp->dfint, npnt);
  dgsp->raniso = ra;				// save aniso. radius 
}

//  gsp_calc_dist_pars: set parameters for numerical integration.
//  _____________________________________________________________

void gsp_calc_dist_pars(double *eap, double *erp, bool *upp)
{
  epsabs = (eap != NULL ? *eap : epsabs);	// set non-null parameters
  epsrel = (erp != NULL ? *erp : epsrel);
  usephi = (upp != NULL ? *upp : usephi);
}

#ifdef UTILITY

string defv[] = {		";Tabulate distribution function",
  "gsp=???",			";GSP for density profile",
  "grav=",			";GSP for gravitational field.",
				";If blank, use density for field.",
  "raniso=-1",			";Osipkov-Merritt anisotropy radius.",
				";If raniso <= 0, find isotropic DF.",
  "epsabs=0.0",			";Integration absolute tolerance",
  "epsrel=1.0e-4",		";Integration relative tolerance.",
				";Need to relax for anisotropic models.",
  "npoint=161",			";Number of energy values to list",
  "rrange=1/1024:1024",		";Range of radii (mapped to energy)",
  "VERSION=2.0",		";Josh Barnes  1y July 2017",
  NULL,
};

int main(int argc, string argv[])
{
  stream dstr, gstr;
  gsprof *dgsp, *ggsp;
  int npnt;
  real rrange[2];
  double ea, er, lgrs, E;

  initparam(argv, defv);
  dstr = stropen(getparam("gsp"), "r");
  get_history(dstr);
  dgsp = gsp_read(dstr);
  if (! strnull(getparam("grav"))) {
    gstr = stropen(getparam("grav"), "r");
    get_history(gstr);
    ggsp = gsp_read(gstr);
  } else
    ggsp = dgsp;
  ea = getdparam("epsabs");
  er = getdparam("epsrel");
  gsp_calc_dist_pars(&ea, &er, NULL);	// set all parameters
  gsp_calc_dist(dgsp, ggsp, getdparam("raniso"));
  npnt = getiparam("npoint");
  setrange(rrange, getparam("rrange"));
  lgrs = log2(rrange[1] / rrange[0]) / (npnt - 1);
  eprintf("[%s: lgrs = %f (%a)]\n", getprog(), lgrs, lgrs);
  if (getdparam("raniso") <= 0.0)
    printf("#%5s  %16s  %16s  %16s\n", "i", "E", "F(E)", "f(E)");
  else
    printf("#%5s  %16s  %16s  %16s\n", "i", "Q", "F(Q)", "f(Q)");
  for (int i = 0; i < npnt; i++) {
    E = gsp_phi(ggsp, rrange[0] * pow(2.0, lgrs * i));
    printf("%6d  %16.8e  %16.8e  %16.8e\n", i, E,
	   gsp_dist_integ(dgsp, E), gsp_dist(dgsp, E));
  }
  if (ggsp != dgsp)
    gsp_free(ggsp);
  gsp_free(dgsp);
  return 0;
}

#endif

#ifdef TESTBED

string defv[] = {		";Test distribution function calculation",
  "model=hernquist",		";Select hernquist or plummer",
  "raniso=-1",			";Osipkov-Merritt anisotropy radius.",
				";If raniso <= 0, find isotropic DF.",
  "npoint=1281",		";Number of radial points in tables",
  "rrange=1/1024:1024",		";Range of radii used in tables",
  "epsabs=0.0",			";Integration absolute tolerance",
  "epsrel=1.0e-4",		";Integration relative tolerance",
				";Need to relax for anisotropic models.",
  "usephi=true",		";Use phi values for energy points",
  "ntest=-1",			";Compare DF at values in energy array.",
				";If positive, use equally-spaced values.",
  "VERSION=2.0",		";Josh Barnes  14 July 2017",
  NULL,
};

#define MTEST  1.0
#define ATEST  1.0

//  f_hernquist: compute DF for Hernquist model.
//  ____________________________________________

double f_hernquist(double E, double ra)
{
  double v = sqrt(MTEST / ATEST);
  double q = sqrt(- E * ATEST / MTEST);
  double fiso;

  fiso = (MTEST / (8*M_SQRT2 * gsl_pow_3(M_PI*ATEST*v))) * pow(1-q*q, -2.5) *
         (3*asin(q) + q*sqrt(1-q*q) * (1-2*q*q) * (8*q*q*q*q - 8*q*q - 3));
  if (ra <= 0.0)
    return fiso;
  else
    return (fiso + (MTEST / (M_SQRT2 * gsl_pow_3(M_PI * ATEST * v))) *
	           gsl_pow_2(ATEST / ra) * q * (1 - 2 * q*q));
}

//  f_plummer: compute DF for Plummer model.
//  ________________________________________

double f_plummer(double E, double ra)
{
  double F = (24.0 * M_SQRT2 * ATEST*ATEST /
	      (7.0 * gsl_pow_3(M_PI) * gsl_pow_4(MTEST)));
  
  if (ra > 0.0)
    error("%s.f_plummer: anisotropic version not available\n", getprog());
  return (F * pow(-E, 3.5));
}

int main(int argc, string argv[])
{
  real rrange[2];
  gsprof *gsp;
  double (*f_func)(double, double), ra, ea, er, E, f_exact, f_numer;
  double relerr, avgerr = 0.0;
  bool up;
  int npnt, ntst;

  initparam(argv, defv);
  setrange(rrange, getparam("rrange"));
  if (streq(getparam("model"), "hernquist")) {
    gsp = gsp_gamma(1.0, MTEST, ATEST, getiparam("npoint"),
		    rrange[0], rrange[1]);
    f_func = &f_hernquist;
  } else if (streq(getparam("model"), "plummer")) {
    gsp = gsp_plum(MTEST, ATEST, getiparam("npoint"),
		   rrange[0], rrange[1]);
    f_func = &f_plummer;
  } else
    error("%s: model=%s not recognized\n", getprog(), getparam("model"));
  ra = getdparam("raniso");
  ea = getdparam("epsabs");
  er = getdparam("epsrel");
  up = getbparam("usephi");
  gsp_calc_dist_pars(&ea, &er, &up);	// set all parameters
  gsp_calc_dist(gsp, gsp, ra);		// perform calculation
  npnt = gsp->npoint;
  ntst = getiparam("ntest");
  if (ntst != 0) {
    printf("#%5s  %16s  %16s  %16s  %16s  %16s\n",
	   "i", "E", "F(E)", "f(E)", "f_exact(E)", "rel. error");
    for (int i = 0; i < (ntst > 0 ? ntst : npnt) ; i++) {
      if (ntst > 0)
	E = gsp->energy[0] + (gsp->energy[npnt-1] - gsp->energy[0]) *
	                       (ntst > 1 ? ((double) i) / (ntst - 1) : 1);
      else
	E = gsp->energy[i];
      f_numer = gsp_dist(gsp, E);
      f_exact = (*f_func)(E, ra);
      relerr = (f_numer - f_exact) / f_exact;
      printf("%6d  %16.8e  %16.8e  %16.8e  %16.8e  %16.8e\n",
	     i, E, gsp_dist_integ(gsp, E), f_numer, f_exact, relerr);
      if (1 < i && i < (ntst > 0 ? ntst : npnt) - 2)
	avgerr += ABS(relerr) / ((ntst > 0 ? ntst : npnt) - 4);
    }
    eprintf((avgerr < 1.0e-6 ? "[%s: average error = %e (less endpoints)]\n" :
	     "[%s: warning: average error = %e (less endpoints)]\n"),
	    getprog(), avgerr);
  }
  gsp_free(gsp);
  return 0;
}

#endif
