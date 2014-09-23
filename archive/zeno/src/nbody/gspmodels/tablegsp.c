/*
 * tablegsp.c: convert ascii density profile table to GSP.
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
  "nmax=1026",				";Maximum table length",
  "VERSION=1.1",			";Josh Barnes  19 June 2012",
  NULL,
};

double *radius_tab, *density_tab, *mass_tab;
int nmax, ntab;
real alpha, beta, mtot;

void read_table(string in, bool smooth);
void integ_mass(bool update);
void fit_params(void);
void write_gsp(string out, real mass_scale);

int main(int argc, string argv[])
{
  real mass_scale;

  initparam(argv, defv);
  nmax = getiparam("nmax");
  radius_tab = (double *) allocate(nmax * sizeof(double));
  density_tab = (double *) allocate(nmax * sizeof(double));
  mass_tab = (double *) allocate(nmax * sizeof(double));
  read_table(getparam("in"), getbparam("smooth"));
  integ_mass(getbparam("smooth"));
  fit_params();
  mass_scale = (strnull(getparam("mtot")) ? 1 : getdparam("mtot") / mtot);
  if (!strnull(getparam("out")))
    write_gsp(getparam("out"), mass_scale);
  return (0);
}

void read_table(string in, bool smooth)
{
  stream instr = stropen(in, "r");
  char inbuf[129];
  double radius, rho0, rho1, mass;

  if (fscanf(instr, "#%128[^\n]\n", inbuf) != 1)
    error("%s: can't read first line\n", getprog());
  eprintf("[%s: \"#%s\"]\n", getprog(), inbuf);
  if (fscanf(instr, "#%128[^\n]\n", inbuf) != 1)
    error("%s: can't read second line\n", getprog());
  if (fscanf(instr, "%128[^\n]\n", inbuf) != 1)
    error("%s: can't read third line\n", getprog());
  if (sscanf(inbuf, "%lf Infinity %lf %lf %*s", &radius, &mass, &rho1) == 3) {
    eprintf("[%s: skipping third line; radius = %f]\n", getprog(), radius);    
    ntab = 0;
  } else if (sscanf(inbuf, "%lf %lf %lf %lf %*s",
		    &radius, &rho0, &mass, &rho1) == 4) {
    radius_tab[0] = radius;
    density_tab[0] = (smooth ? rho1 : rho0);
    mass_tab[0] = mass;
    ntab = 1;
  } else
    error("%s: can't interpret third line\n", getprog());
  while (fscanf(instr, "%lf %lf %lf %lf %*s\n",
		&radius, &rho0, &mass, &rho1) == 4) {
    if (ntab < nmax) {
      radius_tab[ntab] = radius;
      density_tab[ntab] = (smooth ? rho1 : rho0);
      mass_tab[ntab] = mass;
    }
    ntab++;
  }
  if (ntab > nmax) {
    eprintf("[%s: warning: truncating table to %d values]\n", getprog(), nmax);
    ntab = nmax;
  }
  eprintf("[%s: radius[] = %e:%e  (%d values)\n", getprog(),
	  radius_tab[0], radius_tab[ntab-1], ntab);
}

void integ_mass(bool update)
{
  double *dmdr_tab = (double *) allocate(ntab * sizeof(double));
  gsl_spline *spl_dmdr = gsl_spline_alloc(gsl_interp_cspline, ntab);
  gsl_interp_accel *acc_dmdr = gsl_interp_accel_alloc();
  int i, stat;
  double alpha, mass_int, del_mass, mass_end = mass_tab[ntab - 1];
  
  for (i = 0; i < ntab; i++)
    dmdr_tab[i] = 4 * M_PI * gsl_pow_2(radius_tab[i]) * density_tab[i];
  gsl_spline_init(spl_dmdr, radius_tab, dmdr_tab, ntab);
  if (radius_tab[0] > 0.0) {
    alpha = rlog10(density_tab[1] / density_tab[0]) /
              rlog10(radius_tab[1] / radius_tab[0]);
    mass_int = (4 * M_PI / (alpha + 3)) *
                    gsl_pow_3(radius_tab[0]) * density_tab[0];
  } else
    mass_int = 0.0;
  if (update)
    mass_tab[0] = mass_int;
  for (i = 1; i < ntab; i++) {
    stat = gsl_spline_eval_integ_e(spl_dmdr, radius_tab[i-1], radius_tab[i],
				   acc_dmdr, &del_mass);
    if (stat != 0)
      error("%s: spline error: %s\n", getprog(), gsl_strerror(stat));
    mass_int = mass_int + del_mass;
    if (update)
      mass_tab[i] = mass_int;
  }
  gsl_interp_accel_free(acc_dmdr);
  gsl_spline_free(spl_dmdr);
  free(dmdr_tab);
  eprintf("[%s: mass[] = %e:%e]\n",
	  getprog(), mass_tab[0], mass_tab[ntab - 1]);
  if (mass_int < 0.99 * mass_end || mass_int > 1.01 * mass_end)
    eprintf("[%s: WARNING: final mass = %e  integ mass = %e]\n",
	    getprog(), mass_end, mass_int);
}

void fit_params(void)
{
  if (radius_tab[0] != 0.0)
    alpha = rlog10(density_tab[1] / density_tab[0]) /
              rlog10(radius_tab[1] / radius_tab[0]);
  else
    alpha = 0.0;
  beta = rlog10(density_tab[ntab-2] / density_tab[ntab-1]) /
              rlog10(radius_tab[ntab-2] / radius_tab[ntab-1]);
  if (beta > -3.0)
    error("%s: total mass diverges (beta = %f)\n", getprog(), beta);
  mtot = mass_tab[ntab-1] -
         (4 * M_PI / (3 + beta)) *
           gsl_pow_3(radius_tab[ntab-1]) * density_tab[ntab-1];
  eprintf("[%s: alpha = %f  beta = %f  mtot = %f]\n", getprog(),
	  alpha, beta, mtot);
}

void write_gsp(string out, real mass_scale)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  stream outstr = stropen(getparam("out"), "w");
  int k, i;

  k = (radius_tab[0] > 0.0 ? 0 : 1);
  gsp->npoint = ntab - k;
  gsp->radius = (real *) allocate((ntab - k) * sizeof(real));
  gsp->density = (real *) allocate((ntab - k) * sizeof(real));
  gsp->mass = (real *) allocate((ntab - k) * sizeof(real));
  for (i = k; i < ntab; i++) {
    gsp->radius[i - k] = radius_tab[i];
    gsp->density[i - k] = mass_scale * density_tab[i];
    gsp->mass[i - k] = mass_scale * mass_tab[i];
  }
  gsp->alpha = alpha;
  gsp->beta = beta;
  gsp->mtot = mass_scale * mtot;
  put_history(outstr);
  put_gsprof(outstr, gsp);
  strclose(outstr);
}
