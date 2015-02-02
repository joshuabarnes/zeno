/*
 * GSPLIST.C: list properties of GSP as a function of radius.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

string defv[] = {		";List general spherical profile",
  "gsp=???",			";Input GSP for density distribution",
  "grav=",		        ";Input GSP for potential computation",
  "beta_a=0.0",			";Anisotropy parameter: beta_a <= 1",
  "npoint=17",			";Number of points to list",
  "rrange=1/256:256",		";Range of radii to list",
  "VERSION=1.1",		";Josh Barnes  31 January 2015",
  NULL,
};

int main(int argc, string argv[])
{
  stream istr;
  gsprof *gsp, *ggsp;
  real beta_a, *sig2, rrange[2], lgrs, r;
  int np, i;

  initparam(argv, defv);
  istr = stropen(getparam("gsp"), "r");
  get_history(istr);
  gsp = get_gsprof(istr);
  strclose(istr);
  if (! strnull(getparam("grav"))) {
    istr = stropen(getparam("grav"), "r");
    get_history(istr);
    ggsp = get_gsprof(istr);
    strclose(istr);
  } else
    ggsp = gsp;
  beta_a = getdparam("beta_a");
  np = getiparam("npoint");
  setrange(rrange, getparam("rrange"));
  lgrs = np > 1 ? rlog2(rrange[1] / rrange[0]) / (np - 1) : 0.0;
  (void) phi_gsp(ggsp, rrange[0]);		// precalculate phi table
  sig2 = calc_sig2_gsp(gsp, ggsp, beta_a);
  printf("#%12s%13s%13s%13s%13s%13s\n",
	 "radius", "rho", "drho/dr", "mass", "phi", "sig^2");
  for (i = 0; i < np; i++) {
    r = rrange[0] * rexp2(lgrs * i);
    printf("%13.4e%13.4e%13.4e%13.4e%13.7f%13.4e\n", r,
	   rho_gsp(gsp, r), drho_gsp(gsp, r), mass_gsp(gsp, r),
	   phi_gsp(ggsp, r), sig2_gsp(gsp, ggsp, beta_a, sig2, r));
  }
  return (0);
}
