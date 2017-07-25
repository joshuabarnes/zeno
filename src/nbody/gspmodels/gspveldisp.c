/*
 * gspveldisp.c: calculate radial velocity dispersion for GSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

//  gsp_sig2: evaluate square of velocity dispersion at given radius.
//  _________________________________________________________________

double gsp_sig2(gsprof *dgsp, double r)
{
  int N = dgsp->npoint - 1, xind;
  double r0 = dgsp->radius[0], rN = dgsp->radius[N], gamma, m0, c0, lgr;

  if (dgsp->sig2 == NULL)
    error("%s.gsp_sig2: call gsp_calc_sig2 first\n", getprog());
  if (r < r0) {
    gamma = 2 * dgsp->beta_a + dgsp->ggsp->alpha + dgsp->alpha + 2;
    m0 = gsp_mass(dgsp->ggsp, r0);
    if (gamma != 0.0) {
      c0 = dgsp->sig2[0] + m0 / (gamma * r0);
      return (- (m0 / r0) * pow(r / r0, 2 + dgsp->ggsp->alpha) / gamma +
	      c0 * pow(r0 / r, 2 * dgsp->beta_a + dgsp->alpha));
    } else {
      c0 = dgsp->sig2[0] + log(r0) * m0 / r0;
      return (- log(r) * m0 / pow(r0, 3 + dgsp->ggsp->alpha) /
	      pow(r, 2 * dgsp->beta_a + dgsp->alpha) +
	      c0 * pow(r0 / r, 2 * dgsp->beta_a + dgsp->alpha));
    }
  }
  if (r > rN) {
    if (dgsp->density[N] > 0) {
      gamma = 2 * dgsp->beta_a + dgsp->ggsp->beta + dgsp->beta + 2;
      return (- dgsp->ggsp->mtot / ((2*dgsp->beta_a + dgsp->beta - 1) * r) +
	      (dgsp->ggsp->mtot - gsp_mass(dgsp->ggsp, r)) / (gamma * r));
    } else
      return 0.0;
  }
#ifdef USELOG2R
  lgr = log2(r);
  xind = N * (lgr - dgsp->lg2rad[0]) / (dgsp->lg2rad[N] - dgsp->lg2rad[0]);
  dgsp->r_acc->cache = MAX(0, MIN(N - 1, xind));
  return gsl_interp_eval(dgsp->sr_spline, dgsp->lg2rad, dgsp->sig2,
			 lgr, dgsp->r_acc);
#else
  return gsl_interp_eval(dgsp->sr_spline, dgsp->radius, dgsp->sig2,
			 r, dgsp->r_acc);
#endif
}

//  dsig2func: compute -d(sig2)/dr.  Sign accounts for integration direction.
//  _________________________________________________________________________

local int dsig2func(double rm, const double sig2[], double dsig2[],
		    void *params)
{
  gsprof *dgsp = ((gsprof *) params);

  dsig2[0] = 0.0;
  if (- rm > 0)
    dsig2[0] += gsp_mass(dgsp->ggsp, - rm) / (rm*rm);
  if (- rm > 0 && sig2[0] > 0)
    dsig2[0] += sig2[0] * (2 * dgsp->beta_a / (- rm) +
			   gsp_grad(dgsp, - rm) / gsp_rho(dgsp, - rm));
  return GSL_SUCCESS;
}

//  gsp_calc_sig2: calculate radial dispersion needed to support the
//  density profile dgsp in the potential of profile ggsp.
//  ________________________________________________________________

void gsp_calc_sig2(gsprof *dgsp, gsprof *ggsp, double beta_a)
{
  gsl_odeiv2_system sys = { dsig2func, NULL, 1, dgsp };
  gsl_odeiv2_driver *drv;
  double rm, sig2[1], gamma;
  int N = dgsp->npoint - 1, stat;

  dgsp->ggsp = ggsp;
  dgsp->beta_a = beta_a;
  if (dgsp->sig2 == NULL)
    dgsp->sig2 = (double *) allocate(dgsp->npoint * sizeof(double));
  drv = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
				      1.0e-6, 1.0e-6, 0.0);
  rm = - dgsp->radius[N];			// start neg to integ forward
  if (dgsp->density[N] > 0.0) {
    gamma = 2 * beta_a + ggsp->beta + dgsp->beta + 2;
    sig2[0] = - ggsp->mtot / ((2 * beta_a + dgsp->beta - 1) * (- rm)) +
      (ggsp->mtot - gsp_mass(ggsp, - rm)) / (gamma * (- rm));
  } else
    sig2[0] = 0.0;
  dgsp->sig2[N] = sig2[0];
  for (int i = N - 1; i >= 0; i--) {
    stat = gsl_odeiv2_driver_apply(drv, &rm, - dgsp->radius[i], sig2);
    if (stat != GSL_SUCCESS)
      error("%s.gsp_calc_sig2: GSL error, status = %d\n", getprog(), stat);
    dgsp->sig2[i] = sig2[0];
  }
  gsl_odeiv2_driver_free(drv);
  if (dgsp->sr_spline == NULL)
    dgsp->sr_spline = gsl_interp_alloc(gsl_interp_akima, dgsp->npoint);
#ifdef USELOG2R
  gsl_interp_init(dgsp->sr_spline, dgsp->lg2rad, dgsp->sig2, dgsp->npoint);
#else
  gsl_interp_init(dgsp->sr_spline, dgsp->radius, dgsp->sig2, dgsp->npoint);
#endif
}

#ifdef UTILITY

string defv[] = {		";Tabulate velocity dispersions",
  "gsp=???",			";Input GSP for density distribution",
  "grav=",			";Input GSP for potential computation",
  "beta=0.0",			";Anisotropy parameter: beta <= 1",
  "npoint=161",			";Number of potential values to list",
  "rrange=1/1024:1024",		";Range of radii tabulated",
  "VERSION=2.0",		";Josh Barnes  16 June 2017",
  NULL,
};

int main(int argc, string argv[])
{
  stream istr;
  gsprof *dgsp, *ggsp;
  double beta, lgrs, r, sig_r, sig_t;
  int npnt;
  real rrange[2];

  initparam(argv, defv);
  istr = stropen(getparam("gsp"), "r");
  get_history(istr);
  dgsp = gsp_read(istr);
  if (! strnull(getparam("grav"))) {
    istr = stropen(getparam("grav"), "r");
    get_history(istr);
    ggsp = gsp_read(istr);
  } else
    ggsp = dgsp;
  beta = getdparam("beta");
  gsp_calc_sig2(dgsp, ggsp, beta);
  npnt = getiparam("npoint");
  setrange(rrange, getparam("rrange"));
  lgrs = log2(rrange[1] / rrange[0]) / (npnt - 1);
  eprintf("[%s: lgrs = %f (%a)]\n", getprog(), lgrs, lgrs);
  printf("#%5s  %16s  %16s  %16s\n", "i", "r", "sig_r", "sig_t");
  for (int i = 0; i < npnt; i++) {
    r = rrange[0] * pow(2.0, lgrs * i);
    sig_r = sqrt(gsp_sig2(dgsp, r));
    sig_t = sqrt(1 - beta) * sig_r;
    printf("%6d  %16.8e  %16.8e  %16.8e\n", i, r,
	   sig_r, sig_t);
  }
  if (ggsp != dgsp)
    gsp_free(ggsp);
  gsp_free(dgsp);
  return 0;
}
#endif

#ifdef TESTBED

string defv[] = {		";Test velocity dispersion calculation",
  "model=hernquist",		";Select hernquist or ...",
  "npoint=769",			";Number of radial points in tables",
  "rrange=1/64:64",		";Range of radii used in tables",
  "VERSION=2.0",		";Josh Barnes  16 June 2017",
  NULL,
};

//  sig2_hernquist: compute sig^2 for isotropic Hernquist (1990) model.
//  ___________________________________________________________________

local double sig2_hernquist(void *p, double r)
{
  double M = ((double *) p)[0], a = ((double *) p)[1];

  return (M / (12 * a)) *
    (log((r + a) / r) * (12 * r * gsl_pow_3(r + a)) / gsl_pow_4(a)  -
     (25 + 52 * r/a + 42 * gsl_pow_2(r/a)  + 12 * gsl_pow_3(r/a)) * r/(r+a));
}

#define MTEST  1.0
#define ATEST  1.0

int main(int argc, string argv[])
{
  real rrange[2];
  gsprof *gsp;
  double (*sig2_func)(void *, double), pars[2] = { MTEST, ATEST };
  double lgrs, r, sig2_numer, sig2_exact, relerr;
  int np;

  initparam(argv, defv);
  setrange(rrange, getparam("rrange"));
  if (streq(getparam("model"), "hernquist")) {
    gsp = gsp_gamma(1.0, MTEST, ATEST, getiparam("npoint"),
		    rrange[0], rrange[1]);
    sig2_func = &sig2_hernquist;
  } else
    error("%s: model=%s not recognized\n", getprog(), getparam("model"));
  gsp_calc_sig2(gsp, gsp, 0.0);

  gsp_test_rad(gsp, gsp_sig2, sig2_func, pars, "sigma^2");

  np = getiparam("npoint");
  setrange(rrange, getparam("rrange"));
  lgrs = log2(rrange[1] / rrange[0]) / (np - 1);
  printf("#%5s  %16s  %16s  %16s  %16s\n",
	 "i", "radius", "sig2", "sig2_exact", "rel. error");
  for (int i = 0; i < np; i++) {
    r = rrange[0] * pow(2.0, lgrs * i);
    sig2_numer = gsp_sig2(gsp, r);
    sig2_exact = (*sig2_func)(pars, r);
    relerr = (sig2_numer - sig2_exact) / sig2_exact;
    printf("%6d  %16.8e  %16.8e  %16.8e  %16.8e\n",
	   i, r, sig2_numer, sig2_exact, relerr);
  }
  gsp_free(gsp);
  return 0;
}

#endif
