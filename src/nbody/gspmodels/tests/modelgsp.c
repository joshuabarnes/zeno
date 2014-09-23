/*
 * MODELGSP.C: convert smoothed density profile to GSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_exp.h>

string defv[] = {			";Convert smoothed profile to GSP",
  "in=???",				";Input file from smoothing calc",
  "out=",				";Output file with GSP",
  "rtaper=4.0",				";Radius where taper begins.",
					";If past last radii, no taper",
  "mtot=",				";If given, rescale to tot mass",
  "nmax=1026",				";Maximum table length",
  "VERSION=1.0",			";Josh Barnes  23 June 2011",
  NULL,
};

double *radi_tab, *dens_tab, *mass_tab;
int nmax, npts;

void read_model(string in);
void integ_mass(void);
void taper_model(real rtaper, real *mass_tot, real *beta_out);
void write_gsp(string out, real mass_scale, real mass_tot, real beta_out);

int main(int argc, string argv[])
{
  real mass_tot, beta_out, mass_scale;

  initparam(argv, defv);
  nmax = getiparam("nmax");
  radi_tab = (double *) allocate(nmax * sizeof(double));
  dens_tab = (double *) allocate(nmax * sizeof(double));
  mass_tab = (double *) allocate(nmax * sizeof(double));
  read_model(getparam("in"));
  integ_mass();
  taper_model(getdparam("rtaper"), &mass_tot, &beta_out);
  mass_scale = (strnull(getparam("mtot")) ? 1 : getdparam("mtot") / mass_tot);
  if (!strnull(getparam("out")))
    write_gsp(getparam("out"), mass_scale, mass_tot, beta_out);
  return (0);
}

void read_model(string in)
{
  stream instr = stropen(in, "r");
  double radi, dens;

  if (fscanf(instr, "#%*[^\n]\n") != 0 || fscanf(instr, "#%*[^\n]\n") != 0)
    error("%s: can't read header\n", getargv0());
  npts = 0;
  while (fscanf(instr, "%lf %*s %*s %lf %*s\n", &radi, &dens) == 2) {
    if (npts == nmax)
      error("%s: input table overflow at radius = %f\n", getargv0(), radi);
    radi_tab[npts] = radi;
    dens_tab[npts] = dens;
    npts++;
  }
  if (radi_tab[0] != 0.0)
    error("%s: expect data to start at radii = 0.0\n", getargv0());
}

void integ_mass(void)
{
  double *dmdr_tab = (double *) allocate(npts * sizeof(double));
  gsl_spline *spl_dmdr = gsl_spline_alloc(gsl_interp_cspline, npts);
  gsl_interp_accel *acc_dmdr = gsl_interp_accel_alloc();
  int i, stat;
  double del_mass;
  
  for (i = 0; i < npts; i++)
    dmdr_tab[i] = 4 * M_PI * gsl_pow_2(radi_tab[i]) * dens_tab[i];
  stat = gsl_spline_init(spl_dmdr, radi_tab, dmdr_tab, npts);
  if (stat != 0)
    error("%s: spl_dmdr init status: %s\n", getargv0(), gsl_strerror(stat));
  mass_tab[0] = 0.0;
  for (i = 1; i < npts; i++) {
    stat = gsl_spline_eval_integ_e(spl_dmdr, radi_tab[i-1], radi_tab[i],
				   acc_dmdr, &del_mass);
    if (stat != 0)
      error("%s: spl_dmdr eval status: %s\n", getargv0(), gsl_strerror(stat));
    mass_tab[i] = mass_tab[i-1] + del_mass;
  }
  gsl_interp_accel_free(acc_dmdr);
  gsl_spline_free(spl_dmdr);
  free(dmdr_tab);
}

void taper_model(real rtaper, real *mass_tot, real *beta_out)
{
  gsl_spline *spl_dens = gsl_spline_alloc(gsl_interp_cspline, npts);
  gsl_interp_accel *acc_dens = gsl_interp_accel_alloc();
  double x, dddr, beta, rstar, kstar;
  int stat, i, itaper;

  stat = gsl_spline_init(spl_dens, radi_tab, dens_tab, npts);
  if (stat != 0)
    error("%s: spl_dens init status: %s\n", getargv0(), gsl_strerror(stat));
  if (rtaper < radi_tab[npts-1]) {
    x = ABS(rtaper - radi_tab[0]);
    for (i = 1; i < npts; i++)
      if (ABS(rtaper - radi_tab[i]) < x) {
	x = ABS(rtaper - radi_tab[i]);
	itaper = i;
      }
    stat = gsl_spline_eval_deriv_e(spl_dens, radi_tab[itaper], acc_dens,
				   &dddr);
    if (stat != 0)
      error("%s: spl_dens eval status: %s\n", getargv0(), gsl_strerror(stat));
    beta = (radi_tab[itaper] / dens_tab[itaper]) * dddr;
    rstar = - radi_tab[itaper] / (2 + beta);
    kstar = gsl_pow_2(radi_tab[itaper]) * dens_tab[itaper] / gsl_sf_exp(2+beta);
    eprintf("[%s: beta = %f  rstar = %f  kstar = %e]\n",
	    getargv0(), beta, rstar, kstar);
    for (i = itaper + 1; i < npts; i++) {
      dens_tab[i] = kstar * gsl_sf_exp(- radi_tab[i]/rstar) /
		       gsl_pow_2(radi_tab[i]);
      mass_tab[i] = mass_tab[itaper] + 4 * M_PI * kstar * rstar *
		      (gsl_sf_exp(2+beta) - gsl_sf_exp(- radi_tab[i]/rstar));
    }
    *mass_tot = mass_tab[itaper] + 4 * M_PI * kstar * rstar * gsl_sf_exp(2+beta);
    *beta_out = - (2 + radi_tab[npts-1]/rstar);
  } else {
    if (npts < 5)
      error("%s: table too small\n", getargv0());
    beta = 0.0;
    eprintf("[%s: beta est =", getargv0());
    for (i = npts - 5; i < npts - 1; i++) {
      stat = gsl_spline_eval_deriv_e(spl_dens, radi_tab[i], acc_dens, &dddr);
      if (stat != 0)
	error("%s: spl_dens eval status: %s\n", getargv0(), gsl_strerror(stat));
      beta += (radi_tab[i] / dens_tab[i]) * dddr / 4.0;
      eprintf(" %g", (radi_tab[i] / dens_tab[i]) * dddr);
    }
    eprintf("]\n");
    *mass_tot = mass_tab[npts-1] -
                (4 * M_PI / (beta + 3)) *
                  gsl_pow_3(radi_tab[npts-1]) * dens_tab[npts-1];
    *beta_out = beta;
  }
  eprintf("[%s: mass_tot = %e  beta_out = %e]\n",
	  getargv0(), *mass_tot, *beta_out);
  gsl_interp_accel_free(acc_dens);
  gsl_spline_free(spl_dens);
}

void write_gsp(string out, real mass_scale, real mass_tot, real beta_out)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  stream outstr = stropen(getparam("out"), "w");
  int i;

  gsp->npoint = npts - 1;
  gsp->radius = (real *) allocate((npts - 1) * sizeof(real));
  gsp->density = (real *) allocate((npts - 1) * sizeof(real));
  gsp->mass = (real *) allocate((npts - 1) * sizeof(real));
  for (i = 1; i < npts; i++) {
    gsp->radius[i-1] = radi_tab[i];
    gsp->density[i-1] = mass_scale * dens_tab[i];
    gsp->mass[i-1] = mass_scale * mass_tab[i];
  }
  gsp->alpha = 0;
  gsp->beta = beta_out;
  gsp->mtot = mass_scale * mass_tot;
  put_history(outstr);
  put_gsprof(outstr, gsp);
  strclose(outstr);
}
