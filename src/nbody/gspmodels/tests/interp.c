#include <stdinc.h>
#include <mathfns.h>
#include <getparam.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

string defv[] = {
  "in=???",
  "step=1",
  "nmax=1026",
  "list=false",
  "VERSION=0",
  NULL,
};

int main(int argc, string argv[])
{
  stream instr;
  int step, nmax, npts, nsub, i, stat;
  double *rad0, *rho0, *rad1, *rho1, radi, rhoi, sum1, sum2;
  bool list;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  initparam(argv, defv);
  instr = stropen(getparam("in"), "r");
  step = getiparam("step");
  nmax = getiparam("nmax");
  list = getbparam("list");
  rad0 = (double *) allocate(nmax * sizeof(double));
  rho0 = (double *) allocate(nmax * sizeof(double));
  rad1 = (double *) allocate(nmax * sizeof(double));
  rho1 = (double *) allocate(nmax * sizeof(double));
  if (fscanf(instr, "%*[^\n]\n") != 0 || fscanf(instr, "%*[^\n]\n") != 0)
    error("%s: can't read header\n", getargv0());
  npts = 0;
  while (fscanf(instr, "%lf %*s %*s %lf %*s\n", &radi, &rhoi) == 2) {
    if (npts == nmax)
      error("%s: input table overflow at rad = %f\n", getargv0(), radi);
    rad0[npts] = radi;
    rho0[npts] = rhoi;
    npts++;
  }
  nsub = 0;
  for (i = 0; i < npts; i++)
    if (i == 0 || (i - 1) % step == 0) {
      rad1[nsub] = rad0[i];
      rho1[nsub] = rho0[i];
      nsub++;
    }
  spline = gsl_spline_alloc(gsl_interp_cspline, nsub);
  stat = gsl_spline_init(spline, rad1, rho1, nsub);
  if (stat != 0)
    error("%s: spline init status: %s\n", getargv0(), gsl_strerror(stat));
  acc = gsl_interp_accel_alloc();
  sum1 = sum2 = 0.0;
  for (i = 1; i < npts; i++) {
    radi = rad0[i];
    stat = gsl_spline_eval_e(spline, radi, acc, &rhoi);
    if (stat != 0)
      error("%s: spline eval status: %s\n", getargv0(), gsl_strerror(stat));
    if (radi <= rad1[nsub - 2]) {		/* don't use last interval */
      sum1 += (rhoi - rho0[i]) / rho0[i];
      sum2 += gsl_pow_2((rhoi - rho0[i]) / rho0[i]);
    }
    if (list)
      printf("%12e  %12e  %12e  %12e\n",
	     radi, rho0[i], rhoi, (rhoi - rho0[i]) / rho0[i]);
  }
  eprintf("[%s: sum1 = %e  sum2 = %e]\n", getargv0(), sum1, sum2);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return 0;
}
