/*
 * halogsp.c: generate profile tables for Navarro, Frenk, White model.
 */

#include "stdinc.h"
#include "getparam.h"
#include "mathfns.h"
#include "assert.h"
#include "gsp.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

#if (!defined(LINUX) && !defined(MACOSX))
#include <ieeefp.h>
#endif

//  gsp_halo_e: initialize Navarro et al model with exponential taper.
//  __________________________________________________________________

gsprof *gsp_halo_e(double m_a, double a, double b,
		   int np, double rmin, double rmax)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  double lgrs, mu, rho_b, gam, m_b, r_i;

  assert(np > 0 && m_a > 0.0 && a > 0.0 && a < b);
  gsp->npoint = np;
  gsp->radius  = (double *) allocate(np * sizeof(double));
  gsp->density = (double *) allocate(np * sizeof(double));
  gsp->mass    = (double *) allocate(np * sizeof(double));
  lgrs = log2(rmax / rmin) / (np - 1);
  mu = m_a / (log(2.0) - 0.5);
  rho_b = (mu / (4*M_PI)) / (b * gsl_pow_2(a + b));
  gam = b / (b + a) - 0.5;
  m_b = mu * (log((b + a) / a) - b / (a + b));
  eprintf("[%s.gsp_halo_e: mu = %f  rho_b = %f]\n", getprog(), mu, rho_b);
  eprintf("[%s.gsp_halo_e: gamma = %f  m_b = %f]\n", getprog(), gam, m_b);
  for (int i = 0; i < np; i++) {
    gsp->radius[i] = r_i = rmin * exp2(lgrs * i);
    if (r_i <= b) {
      gsp->density[i] = (mu / (4*M_PI)) / (r_i * gsl_pow_2(a + r_i));
      gsp->mass[i] = mu * (log1p(r_i / a) - r_i / (a + r_i));
    } else {
      gsp->density[i] = rho_b * gsl_pow_2(b/r_i) *
	exp(-2 * gam * (r_i/b - 1));
      gsp->mass[i] = m_b + (2*M_PI / gam) * gsl_pow_3(b) * rho_b *
	(1 - exp(-2 * gam * (r_i/b - 1)));
    }
  }
  gsp->alpha = -1.0;
  gsp->beta = -2 * gam * rmax / b - 2;
  gsp->mtot = m_b + (2*M_PI / gam) * gsl_pow_3(b) * rho_b;
  eprintf("[%s.gsp_halo_e: beta = %f  mtot = %f]\n",
	  getprog(), gsp->beta, gsp->mtot);
  return gsp;
}

//  gsp_halo_g: initialize Navarro et al model with gaussian taper.
//  _______________________________________________________________

gsprof *gsp_halo_g(double m_a, double a, double b,
		   int np, double rmin, double rmax)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  double lgrs, mu, rho_b, gam, m_b, r_i, ghalf, pi3half;
  
  assert(np > 0 && m_a > 0.0 && a > 0.0 && a < b);
  gsp->npoint = np;
  gsp->radius  = (double *) allocate(np * sizeof(double));
  gsp->density = (double *) allocate(np * sizeof(double));
  gsp->mass    = (double *) allocate(np * sizeof(double));
  lgrs = log2(rmax / rmin) / (np - 1);
  mu = m_a / (log(2.0) - 0.5);
  rho_b = (mu / (4*M_PI)) / (b * gsl_pow_2(a + b));
  gam = b / (b + a) - 0.5;
  m_b = mu * (log((b + a) / a) - b / (a + b));
  eprintf("[%s.gsp_halo_g: mu = %f  rho_b = %f]\n", getprog(), mu, rho_b);
  eprintf("[%s.gsp_halo_g: gamma = %f  m_b = %f]\n", getprog(), gam, m_b);
  ghalf = sqrt(gam);
  pi3half = sqrt(gsl_pow_3(PI));
  for (int i = 0; i < np; i++) {
    r_i = gsp->radius[i] = rmin * exp2(lgrs * i);
    if (r_i <= b) {
      gsp->density[i] = (mu / (4*M_PI)) / (r_i * gsl_pow_2(a + r_i));
      gsp->mass[i] = mu * (log1p(r_i / a) - r_i / (a + r_i));
    } else {
      gsp->density[i] = rho_b * gsl_pow_2(b/r_i) *
	exp(- gam * (gsl_pow_2(r_i/b) - 1));
      gsp->mass[i] = m_b + (2 * pi3half / ghalf) * gsl_pow_3(b) * rho_b *
	exp(gam) * (erf(ghalf * r_i/b) - erf(ghalf));
    }
  }
  gsp->alpha = -1.0;
  gsp->beta = -2 * gam * gsl_pow_2(rmax / b) - 2;
  gsp->mtot = m_b + (2 * pi3half / ghalf) * gsl_pow_3(b) * rho_b *
    exp(gam) * (1.0 - erf(ghalf));
  eprintf("[%s.gsp_halo_g: beta = %f  mtot = %f]\n",
	  getprog(), gsp->beta, gsp->mtot);
  return gsp;
}

//  gsp_halo_sw: initialize Navarro et al model with fast taper
//  (Springel & White 1998, astro-ph/9807320).
//  ___________________________________________________________

gsprof *gsp_halo_sw(double m_a, double a, double b,
		    int np, double rmin, double rmax)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  double lgrs, mu, rho_b, gam, m_b, gscale, r_i;
  
  assert(np > 0 && m_a > 0.0 && a > 0.0 && a < b);
  gsp->npoint = np;
  gsp->radius  = (double *) allocate(np * sizeof(double));
  gsp->density = (double *) allocate(np * sizeof(double));
  gsp->mass    = (double *) allocate(np * sizeof(double));
  lgrs = log2(rmax / rmin) / (np - 1);
  mu = m_a / (log(2.0) - 0.5);
  rho_b = (mu / (4*M_PI)) / (b * gsl_pow_2(a + b));
  gam = (b / a) - (a + 3 * b) / (a + b);
  m_b = mu * (log((b + a) / a) - b / (a + b));
  gscale = exp(gam * log(a/b) + b/a + gsl_sf_lngamma(3 + gam));
  eprintf("[%s.gsp_halo_sw: mu = %f  rho_b = %f]\n", getprog(), mu, rho_b);
  eprintf("[%s.gsp_halo_sw: gamma = %f  m_b = %f]\n", getprog(), gam, m_b);
  for (int i = 0; i < np; i++) {
    gsp->radius[i] = r_i = rmin * exp2(lgrs * i);
    if (r_i <= b) {
      gsp->density[i] = (mu / (4*M_PI)) / (r_i * gsl_pow_2(a + r_i));
      gsp->mass[i] = mu * (log1p(r_i / a) - r_i / (a + r_i));
    } else {
      gsp->density[i] = rho_b * pow(r_i / b, gam) * exp(- (r_i - b) / a);
      if (isnan((double) gsp->density[i])) {
	gsp->density[i] = 0.0;
	eprintf("[%s.gsp_halo_sw: warning: density is nan]\n", getprog());
      }
      gsp->mass[i] = m_b + (4*M_PI * rho_b * gsl_pow_3(a) * gscale *
			    (gsl_sf_gamma_inc_P(3 + gam, r_i/a) -
			     gsl_sf_gamma_inc_P(3 + gam, b/a)));
    }
  }
  if (gsp->density[np-1] == 0.0)
    eprintf("[%s.gsp_halo_sw: WARNING: density vanishes]\n", getprog());
  if (gsp->mass[np-1] == gsp->mass[np-2])
    eprintf("[%s.gsp_halo_sw: WARNING: mass converges]\n", getprog());
  gsp->alpha = -1.0;
  gsp->beta = gam - rmax / a;
  gsp->mtot = m_b + (4*M_PI * rho_b * gsl_pow_3(a) * gscale *
		     (1 - gsl_sf_gamma_inc_P(3 + gam, b/a)));
  eprintf("[%s.gsp_halo_sw: beta = %f  mtot = %f  m_b/mtot = %f]\n",
	  getprog(), gsp->beta, gsp->mtot, m_b / gsp->mtot);
  return gsp;
}

#ifdef UTILITY

string defv[] = {		";Generate profile for halo model",
  "out=???",			";Output file for profile tables",
  "m_a=1.0",			";Mass within radius a",
  "a=1.0",			";Radial scale of model",
  "b=4.0",			";Radius to begin taper",
  "taper=exp",			";Tapering: exp, gauss, or sw",
  "npoint=257",			";Number of points in tables",
  "rrange=1/4096:16",		";Range of radii tabulated",
  "VERSION=2.0",		";Josh Barnes  17 June 2017",
  NULL,
};

int main(int argc, string argv[])
{
  real rrange[2];
  gsprof *gsp;
  stream ostr;

  initparam(argv, defv);
  setrange(rrange, getparam("rrange"));
  if (getdparam("a") < rrange[0] || rrange[1] < getdparam("b"))
    error("%s: rrange does not include both a and b\n", getprog());
  if (streq(getparam("taper"), "exp"))
    gsp = gsp_halo_e(getdparam("m_a"), getdparam("a"), getdparam("b"),
		     getiparam("npoint"), rrange[0], rrange[1]);
  else if (streq(getparam("taper"), "gauss"))
    gsp = gsp_halo_g(getdparam("m_a"), getdparam("a"), getdparam("b"),
		     getiparam("npoint"), rrange[0], rrange[1]);
  else if (streq(getparam("taper"), "sw"))
    gsp = gsp_halo_sw(getdparam("m_a"), getdparam("a"), getdparam("b"),
		      getiparam("npoint"), rrange[0], rrange[1]);
  else
    error("%s: unknown taper option %s\n", getprog(), getparam("taper"));
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  gsp_write(ostr, gsp);
  fflush(NULL);
  return 0;
}

#endif
