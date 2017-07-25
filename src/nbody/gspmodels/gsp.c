/*
 * gsp.c: routines for operating on general spherical profiles.
 */

#include <assert.h>
#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "filestruct.h"
#include "gsp.h"

local void init_lg2rad(gsprof *gsp);

//  gsp_rho: evaluate density at given radius.
//  __________________________________________

double gsp_rho(gsprof *gsp, double r)
{
  int N = gsp->npoint - 1, xind;
  double lgr;

  if (r < 0)
    error("%s.gsp_rho: undefined for r = %g < 0\n", getprog(), r);
  if (r <= gsp->radius[0])
    return gsp->density[0] * pow(r / gsp->radius[0], gsp->alpha);
  if (r >= gsp->radius[N])
    return gsp->density[N] * pow(r / gsp->radius[N], gsp->beta);
  if (gsp->dr_spline == NULL) {
    gsp->dr_spline = gsl_interp_alloc(gsl_interp_akima, gsp->npoint);
#ifdef USELOG2R
    init_lg2rad(gsp);
    gsl_interp_init(gsp->dr_spline, gsp->lg2rad, gsp->density, gsp->npoint);
#else
    gsl_interp_init(gsp->dr_spline, gsp->radius, gsp->density, gsp->npoint);
#endif
    if (gsp->r_acc == NULL)
      gsp->r_acc = gsl_interp_accel_alloc();
  }
#ifdef USELOG2R
  lgr = log2(r);
  xind = N * (lgr - gsp->lg2rad[0]) / (gsp->lg2rad[N] - gsp->lg2rad[0]);
  gsp->r_acc->cache = MAX(0, MIN(N - 1, xind));
  return gsl_interp_eval(gsp->dr_spline, gsp->lg2rad, gsp->density,
			 lgr, gsp->r_acc);
#else
  return gsl_interp_eval(gsp->dr_spline, gsp->radius, gsp->density,
			 r, gsp->r_acc);
#endif
}

//  gsp_grad: evaluate radial derivative of density at given radius.
//  ________________________________________________________________

double gsp_grad(gsprof *gsp, double r)
{
  int N = gsp->npoint - 1, xind;
  double lgr;

  if (r < 0)
    error("%s.gsp_grad: undefined for r = %g < 0\n", getprog(), r);
  if (r <= gsp->radius[0])
    return (gsp->density[0] * pow(r / gsp->radius[0], gsp->alpha) *
	    gsp->alpha / r);
  if (r >= gsp->radius[N])			// less glitches at large r
    return (gsp->density[N] * pow(r / gsp->radius[N], gsp->beta) *
	    gsp->beta / r);
  if (gsp->dr_spline == NULL) {
    gsp->dr_spline = gsl_interp_alloc(gsl_interp_akima, gsp->npoint);
#ifdef USELOG2R
    init_lg2rad(gsp);
    gsl_interp_init(gsp->dr_spline, gsp->lg2rad, gsp->density, gsp->npoint);
#else
    gsl_interp_init(gsp->dr_spline, gsp->radius, gsp->density, gsp->npoint);
#endif
    if (gsp->r_acc == NULL)
      gsp->r_acc = gsl_interp_accel_alloc();
  }
#ifdef USELOG2R
  lgr = log2(r);
  xind = N * (lgr - gsp->lg2rad[0]) / (gsp->lg2rad[N] - gsp->lg2rad[0]);
  gsp->r_acc->cache = MAX(0, MIN(N - 1, xind));
  return gsl_interp_eval_deriv(gsp->dr_spline,	gsp->lg2rad, gsp->density,
			       lgr, gsp->r_acc) / (r * M_LN2);
#else
  return gsl_interp_eval_deriv(gsp->dr_spline,	gsp->radius, gsp->density,
			       r, gsp->r_acc);
#endif
}

//  gsp_mass: evaluate enclosed mass at given radius.
//  _________________________________________________

double gsp_mass(gsprof *gsp, double r)
{
  int N = gsp->npoint - 1, xind;
  double lgr;

  if (r < 0)
    error("%s.gsp_mass: undefined for r = %g < 0\n", getprog(), r);
  if (r <= gsp->radius[0])
    return (gsp->mass[0] * pow(r / gsp->radius[0], 3 + gsp->alpha));
  else if (r >= gsp->radius[N])
    return (gsp->mtot - (gsp->mtot - gsp->mass[N]) *
	    pow(r / gsp->radius[N], 3 + gsp->beta));
  if (gsp->mr_spline == NULL) {
    gsp->mr_spline = gsl_interp_alloc(gsl_interp_akima, gsp->npoint);
#ifdef USELOG2R
    init_lg2rad(gsp);
    gsl_interp_init(gsp->mr_spline, gsp->lg2rad, gsp->mass, gsp->npoint);
#else
    gsl_interp_init(gsp->mr_spline, gsp->radius, gsp->mass, gsp->npoint);
#endif
    if (gsp->r_acc == NULL)
      gsp->r_acc = gsl_interp_accel_alloc();
  }
#ifdef USELOG2R
  lgr  = log2(r);
  xind = N * (lgr - gsp->lg2rad[0]) / (gsp->lg2rad[N] - gsp->lg2rad[0]);
  gsp->r_acc->cache = MAX(0, MIN(N - 1, xind));
  return gsl_interp_eval(gsp->mr_spline, gsp->lg2rad, gsp->mass,
			  lgr, gsp->r_acc);
#else
  return gsl_interp_eval(gsp->mr_spline, gsp->radius, gsp->mass,
			  r, gsp->r_acc);
#endif
}

//  gsp_mass_rad: evaluate radius enclosing given mass.
//  ___________________________________________________

double gsp_mass_rad(gsprof *gsp, double m)
{
  int N = gsp->npoint - 1;

  if (m < 0 || m > gsp->mtot)
    error("%s.gsp_mass_rad: undefined for m = %g\n", getprog(), m);
  if (m < gsp->mass[0])
    return (gsp->radius[0] * pow(m / gsp->mass[0], 1 / (3 + gsp->alpha)));
  if (m > gsp->mass[N])
    return (gsp->radius[N] *
	    pow((gsp->mtot - m) / (gsp->mtot - gsp->mass[N]),
		1 / (3 + gsp->beta)));
  if (gsp->rm_spline == NULL) {
    for (int n = 0; n < N; n++)
      if (gsp->mass[n] == gsp->mass[n+1]) {
	eprintf("[%s.gsp_mass_radius: warning: truncating mass table;"
		" N = %d -> %d]\n", getprog(), N, n);
	N = n;
	break;
      }
    gsp->rm_spline = gsl_interp_alloc(gsl_interp_akima, N + 1);
    gsl_interp_init(gsp->rm_spline, gsp->mass, gsp->radius, N + 1);
  }
  return gsl_interp_eval(gsp->rm_spline, gsp->mass, gsp->radius,
			  m, NULL);
}

//  init_lg2rad: if needed, construct array of log2(radius).
//  ________________________________________________________

local void init_lg2rad(gsprof *gsp)
{
  if (gsp->lg2rad == NULL) {
    gsp->lg2rad = (double *) allocate(gsp->npoint * sizeof(double));
    for (int i = 0; i < gsp->npoint; i++)
      gsp->lg2rad[i] = log2(gsp->radius[i]);
  }
}

//  gsp_read: read profile tables from input stream.
//  ________________________________________________

gsprof *gsp_read(stream istr)
{
  string gsptag = "GeneralSphericalProfile";
  gsprof *gsp;
    
  if (get_tag_ok(istr, "FiniteSphericalProfile"))	// old-style file?
    gsptag = "FiniteSphericalProfile";			// use old tag
  gsp = (gsprof *) allocate(sizeof(gsprof));
  get_set(istr, gsptag);
  get_data(istr, "Npoint", IntType, &gsp->npoint, 0);
  gsp->radius = (double *) allocate(gsp->npoint * sizeof(double));
  gsp->density = (double *) allocate(gsp->npoint * sizeof(double));
  gsp->mass = (double *) allocate(gsp->npoint * sizeof(double));
  get_data(istr, "Radius", DoubleType, gsp->radius, gsp->npoint, 0);
  get_data(istr, "Density", DoubleType, gsp->density, gsp->npoint, 0);
  get_data(istr, "Mass", DoubleType, gsp->mass, gsp->npoint, 0);
  get_data(istr, "Alpha", DoubleType, &gsp->alpha, 0);
  get_data(istr, "Beta", DoubleType, &gsp->beta, 0);
  get_data(istr, "Mtot", DoubleType, &gsp->mtot, 0);
  get_tes(istr, gsptag);
  return (gsp);
}

//  gsp_write: write profile tables to output stream.
//  _________________________________________________

void gsp_write(stream ostr, gsprof *gsp)
{
  put_set(ostr, "GeneralSphericalProfile");
  put_data(ostr, "Npoint", IntType, &gsp->npoint, 0);
  put_data(ostr, "Radius", DoubleType, gsp->radius, gsp->npoint, 0);
  put_data(ostr, "Density", DoubleType, gsp->density, gsp->npoint, 0);
  put_data(ostr, "Mass", DoubleType, gsp->mass, gsp->npoint, 0);
  put_data(ostr, "Alpha", DoubleType, &gsp->alpha, 0);
  put_data(ostr, "Beta", DoubleType, &gsp->beta, 0);
  put_data(ostr, "Mtot", DoubleType, &gsp->mtot, 0);
  put_tes(ostr, "GeneralSphericalProfile");
}

//  gsp_free: deallocate a gsp and associated tables.
//  _________________________________________________

void gsp_free(gsprof *gsp)
{
  if (gsp->radius != NULL)
    free((void *) gsp->radius);
  if (gsp->lg2rad != NULL)
    free((void *) gsp->lg2rad);
  if (gsp->density != NULL)
    free((void *) gsp->density);
  if (gsp->mass != NULL)
    free((void *) gsp->mass);
  if (gsp->phi != NULL)
    free((void *) gsp->phi);
  if (gsp->sig2 != NULL)
    free((void *) gsp->sig2);
  if (gsp->dfint != NULL)
    free((void *) gsp->dfint);
  if (gsp->energy != NULL)			// WARNING: don't want to do
    free((void *) gsp->energy);			// this if energy == phi
  if (gsp->dr_spline != NULL)
    gsl_interp_free(gsp->dr_spline);
  if (gsp->mr_spline != NULL)
    gsl_interp_free(gsp->mr_spline);
  if (gsp->rm_spline != NULL)
    gsl_interp_free(gsp->rm_spline);
  if (gsp->pr_spline != NULL)
    gsl_interp_free(gsp->pr_spline);
  if (gsp->rp_spline != NULL)
    gsl_interp_free(gsp->rp_spline);
  if (gsp->sr_spline != NULL)
    gsl_interp_free(gsp->sr_spline);
  if (gsp->df_spline != NULL)
    gsl_interp_free(gsp->df_spline);
  if (gsp->r_acc != NULL) {
    eprintf("[%s.gsp_free: radial evals = %d  hit_rate = %.4f]\n",
	    getprog(), gsp->r_acc->miss_count + gsp->r_acc->hit_count,
	    ((real) gsp->r_acc->hit_count) /
	    (gsp->r_acc->miss_count + gsp->r_acc->hit_count));
    gsl_interp_accel_free(gsp->r_acc);
  }
  free((void *) gsp);
}

//  gsp_test_rad: check function of r for accuracy.
//  Extrapolation is tested and reported separately.
//  ________________________________________________

void gsp_test_rad(gsprof *gsp, double (*gsp_func)(gsprof *, double),
		  double (*ref_func)(void *, double), void *pars, string fun)
{
  int n = gsp->npoint, k;
  double eta = pow(gsp->radius[n-1] / gsp->radius[0], 1.0 / (n - 1));
  double err_avg[3] = {0, 0, 0}, err_max[3] = {0, 0, 0}, r, gf, rf, err;

  for (int i = 1 - n/4; i < n + n/4; i++) {
    r = i < 1 ? gsp->radius[0] * pow(eta, i - 0.5) :
        i < n ? sqrt(gsp->radius[i-1] * gsp->radius[i]) :
                gsp->radius[n-1] * pow(eta, i - n + 0.5);
    gf = (*gsp_func)(gsp, r);
    rf = (*ref_func)(pars, r);
    // if (streq(fun, "gradient"))
    //   printf("%8d  %16.8f  %16.8e  %16.8e  %16.8e\n", i, log2(r),
    //	        rf, gf, (gf - rf) / rf);
    err = ABS(gf - rf) / ABS(rf);
    k = i < 1 ? 0 : i < n ? 1 : 2;
    err_avg[k] += err / (k == 1 ? n - 1 : n / 4);
    err_max[k] = MAX(err_max[k], err);
  }
  eprintf("[%s.gsp_test_rad: %s  rmin, rmax = %.6e, %.6e\n"
	  "  average error:  %.6e  %.6e (int)  %.6e (ext)\n"
	  "  maximum error:  %.6e  %.6e (int)  %.6e (ext)]\n",
	  getprog(), fun, gsp->radius[0], gsp->radius[n-1],
	  err_avg[1], err_avg[0], err_avg[2],
	  err_max[1], err_max[0], err_max[2]);
}

//  gsp_test_mass: check function of m for accuracy on internal points.
//  ___________________________________________________________________

void gsp_test_mass(gsprof *gsp, double (*gsp_func)(gsprof *, double),
		   double (*ref_func)(void *, double), void *pars, string fun)
{
  double err_sum = 0, err_max = 0, m, gf, rf, err;

  for (int i = 1; i < gsp->npoint; i++) {
    m = sqrt(gsp->mass[i-1] * gsp->mass[i]);
    gf = (*gsp_func)(gsp, m);
    rf = (*ref_func)(pars, m);
    err = ABS(gf - rf) / ABS(rf);
    err_sum += err;
    err_max = MAX(err_max, err);
  }
  eprintf("[%s.gsp_test_mass: %s  mmin, mmax = %.6e, %.6e\n"
	  "  average error:  %.6e  maximum error:  %6e]\n",
	  getprog(), fun, gsp->mass[0], gsp->mass[gsp->npoint - 1],
	  err_sum / (gsp->npoint - 1), err_max);
}

#ifdef TESTBED

//  main: test profile evaluation and input/output routines.
//  ________________________________________________________

#include "getparam.h"

string defv[] = {		";Test general spherical profile",
  "in=???",			";Input file with GSP",
  "out=",			";Output file for GSP",
  "npoint=17",			";Number of points to list",
  "r0=1/256.0",			";radius of first point",
  "lgrstep=1.0",		";Log2 increment in radius",
  "VERSION=2.0",		";Josh Barnes  4 May 2017",
    NULL,
};

int main(int argc, string argv[])
{
  stream istr, ostr;
  gsprof *gsp;
  int np;
  double r0, lgrs, r;

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  gsp = gsp_read(istr);
  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    gsp_write(ostr, gsp);
    fflush(NULL);
  }
  np = getiparam("npoint");
  r0 = getdparam("r0");
  lgrs = getdparam("lgrstep");
  printf("#%11s%12s%12s%12s%12s%12s\n", "radius",
	 "log rho", "gradient", "mass", "mtot-mass", "radius(m)");
  for (int i = 0; i < np; i++) {
    r = r0 * pow(2.0, lgrs * i);
    printf("%12.5f%12.7f%12.3e%12.8f%12.8f%12.5f\n", r,
	   rlog10(gsp_rho(gsp, r)), gsp_grad(gsp, r),
	   gsp_mass(gsp, r), gsp->mtot - gsp_mass(gsp, r),
	   gsp_mass_rad(gsp, gsp_mass(gsp, r)));
  }
  gsp_free(gsp);
  return 0;
}

#endif
