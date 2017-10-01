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
  "VERSION=2.0",		";Josh Barnes  30 September 2017",
  NULL,
};

int main(int argc, string argv[])
{
  stream istr;
  gsprof *gsp, *ggsp;
  double beta_a, lgrs, r;
  real rrange[2];
  int np, i;

  initparam(argv, defv);
  istr = stropen(getparam("gsp"), "r");
  get_history(istr);
  gsp = gsp_read(istr);
  strclose(istr);
  if (! strnull(getparam("grav"))) {
    istr = stropen(getparam("grav"), "r");
    get_history(istr);
    ggsp = gsp_read(istr);
    strclose(istr);
  } else
    ggsp = gsp;
  beta_a = getdparam("beta_a");
  np = getiparam("npoint");
  setrange(rrange, getparam("rrange"));
  lgrs = np > 1 ? rlog2(rrange[1] / rrange[0]) / (np - 1) : 0.0;
  gsp_calc_phi(ggsp);
  gsp_calc_sig2(gsp, ggsp, beta_a);
  printf("#%12s%13s%13s%13s%13s%13s\n",
	 "radius", "rho", "drho/dr", "mass", "phi", "sig^2");
  for (i = 0; i < np; i++) {
    r = rrange[0] * rexp2(lgrs * i);
    printf("%13.4e%13.4e%13.4e%13.4e%13.7f%13.4e\n", r,
	   gsp_rho(gsp, r), gsp_grad(gsp, r), gsp_mass(gsp, r),
	   gsp_phi(ggsp, r), gsp_sig2(gsp, r));
  }
  return (0);
}
