/*
 * FSPLIST.C: list properties of FSP as a function of radius.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "fsp.h"

string defv[] = {		";List finite spherical profile",
    "fsp=???",			";Input FSP for density distribution",
    "grav=",		        ";Input FSP for potential computation",
    "beta_a=0.0",		";Anisotropy parameter: beta_a <= 1",
    "npoint=17",		";Number of points to list",
    "rrange=1/256:256",		";Range of radii to list",
    "VERSION=1.1",		";Josh Barnes  9 November 2000",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr;
    fsprof *fsp, *gfsp;
    real beta_a, *sig2, rrange[2], lgrs, r;
    int np, i;

    initparam(argv, defv);
    istr = stropen(getparam("fsp"), "r");
    get_history(istr);
    fsp = get_fsprof(istr);
    strclose(istr);
    if (! strnull(getparam("grav"))) {
	istr = stropen(getparam("grav"), "r");
	get_history(istr);
	gfsp = get_fsprof(istr);
	strclose(istr);
    } else
	gfsp = fsp;
    beta_a = getdparam("beta_a");
    np = getiparam("npoint");
    setrange(rrange, getparam("rrange"));
    lgrs = rlog2(rrange[1] / rrange[0]) / (np - 1);
    (void) phi_fsp(gfsp, rrange[0]);		/* precalculate phi table   */
    sig2 = calc_sig2_fsp(fsp, gfsp, beta_a);
    printf("#%11s%12s%12s%12s%12s%12s\n",
	   "radius", "rho", "drho/dr", "mass", "phi", "sig^2");
    for (i = 0; i < np; i++) {
	r = rrange[0] * rexp2(lgrs * i);
	printf("%12.5f%12.3e%12.3e%12.7f%12.7f%12.7f\n", r,
	       rho_fsp(fsp, r), drho_fsp(fsp, r), mass_fsp(fsp, r),
	       phi_fsp(gfsp, r), sig2_fsp(fsp, gfsp, beta_a, sig2, r));
    }
    return (0);
}
