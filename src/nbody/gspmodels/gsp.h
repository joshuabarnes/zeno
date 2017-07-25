/*
 * gsp.h: definitions for general spherical profile model.
 */

#include <gsl/gsl_interp.h>

//  gsprof: structure representing general spherical profile.
//  _________________________________________________________

typedef struct _gsprof {
  int npoint;			// number of points in radial arrays
  double *radius;		// radii used to define model profile
  double *lg2rad;		// log2 of radii, for USELOG2R option
  double *density;		// density tabulated at each radius
  double *mass;			// enclosed mass within each radius
  double *phi;			// grav. potential at each radius
  double *sig2;			// square of radial velocity dispersion
  double *dfint;		// integral of distribution function
  double *energy;		// energy for distribution function
  double alpha;			// density power-law index for small r
  double beta;			// density power-law index for large r
  double mtot;			// total mass of model (inside r=infinity)
  double raniso;		// anisotropy radius for O.-M. models
  struct _gsprof *ggsp;		// gravitational GSP for sig2(radius)
  double beta_a;		// anisotropy parameter for sig2(radius)
  gsl_interp *dr_spline;	// spline-fit for density(radius)
  gsl_interp *mr_spline;	// spline-fit for mass(radius)
  gsl_interp *rm_spline;	// spline-fit for radius(mass)
  gsl_interp *pr_spline;	// spline-fit for phi(radius)
  gsl_interp *rp_spline;	// spline-fit for radius(phi)
  gsl_interp *sr_spline;	// spline-fit for sig2(radius)
  gsl_interp *df_spline;	// spline-fit for distibution function
  gsl_interp_accel *r_acc;	// interpolation accelerator for radius
} gsprof;

//  Basic GSP access macros and functions; see gsp.c for code.
//  __________________________________________________________

#define gsp_mtot(gsp)   ((gsp)->mtot)
#define gsp_alpha(gsp)  ((gsp)->alpha)
#define gsp_beta(gsp)   ((gsp)->beta)

double gsp_rho(gsprof *gsp, double r);
double gsp_grad(gsprof *gsp, double r);
double gsp_mass(gsprof *gsp, double r);
double gsp_mass_rad(gsprof *gsp, double m);

gsprof *gsp_read(stream istr);
void gsp_write(stream ostr, gsprof *);

void gsp_free(gsprof *gsp);

//  Test routines; see gsp.c for code.
//  __________________________________

void gsp_test_rad(gsprof *gsp, double (*gsp_func)(gsprof *, double),
		  double (*ref_func)(void *, double), void *pars, string);
void gsp_test_mass(gsprof *gsp, double (*gsp_func)(gsprof *, double),
		   double (*ref_func)(void *, double), void *pars, string);

//  High-level GSP functions.
//  _________________________

double gsp_phi(gsprof *gsp, double r);
double gsp_phi_rad(gsprof *gsp, double phi);
void gsp_calc_phi(gsprof *gsp);

double gsp_dist(gsprof *dgsp, double E);
double gsp_dist_integ(gsprof *dgsp, double E);
void gsp_calc_dist(gsprof *dgsp, gsprof *ggsp, double ra);
void gsp_calc_dist_pars(double *eap, double *erp, bool *upp);

double gsp_sig2(gsprof *dgsp, double r);
void gsp_calc_sig2(gsprof *dgsp, gsprof *ggsp, double beta_a);

//  Constructor functions for GSP models.
//  _____________________________________

gsprof *gsp_expd(double mtot, double alpha, double zdisk,
		 int np, double rmin, double rmax);

gsprof *gsp_gamma(double gam, double mtot, double ascale,
		  int np, double rmin, double rmax);

gsprof *gsp_halo_e(double, double, double, int, double, double);

gsprof *gsp_halo_g(double, double, double, int, double, double);

gsprof *gsp_halo_sw(double, double, double, int, double, double);

gsprof *gsp_isoth(double, double, double, int, double, double);

gsprof *gsp_plum(double mtot, double ascale,
		 int np, double rmin, double rmax);

gsprof *gsp_poly(double, double, double, int);

//  Transformation functions for GSP models.
//  ________________________________________

gsprof *gsp_smooth(gsprof *gsp, double eps, double kappa, string trace);

//  USELOG2R: if defined, interpolate using log2(radius) instead of radius.
//  Appears to yield more accurate results, at cost of log2() evaluation.
//  _______________________________________________________________________

#define USELOG2R
