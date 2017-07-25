/*
 * tablegsp.c: convert text profile table to GSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_exp.h>

string defv[] = {			";Convert profile table to GSP",
  "in=???",				";Input file from smoothing calc",
  "out=",				";Output file with GSP",
  "smooth=true",			";Read smoothed densities",
  "mtot=",				";If given, rescale mass",
  "nmax=4097",				";Maximum table length",
  "VERSION=2.0",			";Josh Barnes  30 May 2017",
  NULL,
};

double *radius, *density, *mass;
int nmax, ntab;
double alpha, beta, mtot;

void read_table(string in, bool smooth);
void integrate_mass(bool update);
void fit_params(void);
void save_gsp(string out, double mass_scale);

int main(int argc, string argv[])
{
  double mass_scale;

  initparam(argv, defv);
  nmax = getiparam("nmax");
  radius = (double *) allocate(nmax * sizeof(double));
  density = (double *) allocate(nmax * sizeof(double));
  mass = (double *) allocate(nmax * sizeof(double));
  read_table(getparam("in"), getbparam("smooth"));
  integrate_mass(getbparam("smooth"));
  fit_params();
  mass_scale = (strnull(getparam("mtot")) ? 1.0 : getdparam("mtot") / mtot);
  if (!strnull(getparam("out")))
    save_gsp(getparam("out"), mass_scale);
  return (0);
}

void read_table(string in, bool smooth)
{
  stream instr = stropen(in, "r");
  char inbuf[257];
  double rad, rho0, mass0, rho1;

  if (fscanf(instr, "#%256[^\n]\n", inbuf) != 1)
    error("%s: can't read first line\n", getprog());
  eprintf("[%s: \"#%s\"]\n", getprog(), inbuf);
  if (fscanf(instr, "#%256[^\n]\n", inbuf) != 1)
    error("%s: can't read second line\n", getprog());
  ntab = 0;
  while (fscanf(instr, "%256[^\n]\n", inbuf) == 1) {
    if (sscanf(inbuf, "%lf Infinity %lf %*lf %*s", &rad, &mass0) == 2)
      eprintf("[%s: skipping infinite density at radius = %.8e]\n",
	      getprog(), rad);
    else if (sscanf(inbuf, "Infinity %lf %lf %*lf %*s", &rho0, &mass0) == 2)
      eprintf("[%s: skipping infinite radius at mass = %.8e]\n",
	      getprog(), mass0);
    else if (sscanf(inbuf, "%lf %lf %lf %lf %*s",
		    &rad, &rho0, &mass0, &rho1) != 4)
      error("%s: can't parse data line \"%s\"\n", getprog(), inbuf);
    else {
      if (ntab < nmax) {
	radius[ntab] = rad;
	density[ntab] = (smooth ? rho1 : rho0);
	mass[ntab] = mass0;
      }
      ntab++;
    }
  }
  if (ntab > nmax) {
    eprintf("[%s: warning: truncating table to %d values]\n", getprog(), nmax);
    ntab = nmax;
  }
  eprintf("[%s: radius[] = %.8e:%.8e  (%d values)]\n", getprog(),
	  radius[0], radius[ntab-1], ntab);
}

void integrate_mass(bool update)
{
  double *dmdr_tab = (double *) allocate(ntab * sizeof(double));
  gsl_spline *spl_dmdr = gsl_spline_alloc(gsl_interp_cspline, ntab);
  gsl_interp_accel *acc_dmdr = gsl_interp_accel_alloc();
  int stat;
  double alpha, mass_int, delta_mass, mass_end = mass[ntab - 1];
  
  for (int i = 0; i < ntab; i++)
    dmdr_tab[i] = 4*M_PI * gsl_pow_2(radius[i]) * density[i];
  gsl_spline_init(spl_dmdr, radius, dmdr_tab, ntab);
  if (radius[0] > 0.0) {			// estimate internal mass
    alpha = log(density[1] / density[0]) / log(radius[1] / radius[0]);
    mass_int = (4*M_PI / (alpha+3)) * gsl_pow_3(radius[0]) * density[0];
  } else
    mass_int = 0.0;
  if (update)
    mass[0] = mass_int;
  for (int i = 1; i < ntab; i++) {
    stat = gsl_spline_eval_integ_e(spl_dmdr, radius[i-1], radius[i],
				   acc_dmdr, &delta_mass);
    if (stat != 0)
      error("%s: spline error: %s\n", getprog(), gsl_strerror(stat));
    mass_int = mass_int + delta_mass;
    if (update)
      mass[i] = mass_int;
  }
  gsl_interp_accel_free(acc_dmdr);
  gsl_spline_free(spl_dmdr);
  free(dmdr_tab);
  eprintf("[%s: mass[] = %.8e:%.8e]\n", getprog(), mass[0], mass[ntab - 1]);
  if (mass_int < 0.999 * mass_end || mass_int > 1.001 * mass_end)
    eprintf("[%s: WARNING: final mass = %e  integ mass = %e]\n",
	    getprog(), mass_end, mass_int);
}

void fit_params(void)
{
  int N = ntab - 1;

  if (radius[0] != 0.0)
    alpha = log(density[1] / density[0]) / log(radius[1] / radius[0]);
  else
    alpha = 0.0;
  beta = log(density[N-2] / density[N]) / log(radius[N-2] / radius[N]);
  if (beta >= -3.0)
    error("%s: total mass diverges (beta = %.8f)\n", getprog(), beta);
  mtot = mass[N] - (4*M_PI / (3+beta)) * gsl_pow_3(radius[N]) * density[N];
  eprintf("[%s: alpha = %.8f  beta = %.8f  mtot = %.8f]\n",
	  getprog(), alpha, beta, mtot);
}

void save_gsp(string out, double mass_scale)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  stream outstr = stropen(getparam("out"), "w");
  int k;

  k = (radius[0] > 0.0 ? 0 : 1);
  gsp->npoint = ntab - k;
  gsp->radius = &radius[k];
  gsp->density = &density[k];
  gsp->mass = (double *) allocate(gsp->npoint * sizeof(double));
  for (int i = k; i < ntab; i++) {
    gsp->mass[i - k] = mass_scale * mass[i];
  }
  gsp->alpha = alpha;
  gsp->beta = beta;
  gsp->mtot = mass_scale * mtot;
  put_history(outstr);
  gsp_write(outstr, gsp);
  fflush(NULL);
}
