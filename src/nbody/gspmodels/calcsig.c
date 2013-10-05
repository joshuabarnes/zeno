/*
 * CALCSIG.C: calculate radial velocity dispersion for GSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

/*
 * SIG2_GSP: evaluate velocity dispersion at given radius.
 */

real sig2_gsp(gsprof *tgsp, gsprof *mgsp, real beta_a, real *sig2, real r)
{
    int n = tgsp->npoint - 1;
    real r0, rn, gamma, m0, c0;

    r0 = tgsp->radius[0];
    rn = tgsp->radius[n];
    if (r < r0) {
	gamma = 2 * beta_a + mgsp->alpha + tgsp->alpha + 2;
	m0 = mass_gsp(mgsp, r0);
	if (gamma != 0.0) {
	    c0 = sig2[0] + m0 / (gamma * r0);
	    return (- (m0 / r0) * rpow(r / r0, 2 + mgsp->alpha) / gamma +
		      c0 * rpow(r0 / r, 2 * beta_a + tgsp->alpha));
	} else {
	    c0 = sig2[0] + rlog(r0) * m0 / r0;
	    return (- rlog(r) * m0 / rpow(r0, 3 + mgsp->alpha) /
		      rpow(r, 2 * beta_a + tgsp->alpha) +
		        c0 * rpow(r0 / r, 2 * beta_a + tgsp->alpha));
	}
    } else if (r > rn) {
	if (tgsp->density[n] > 0) {
	    gamma = 2 * beta_a + mgsp->beta + tgsp->beta + 2;
	    return (- mgsp->mtot / ((2*beta_a + tgsp->beta - 1) * r) +
		      (mgsp->mtot - mass_gsp(mgsp, r)) / (gamma * r));
	} else
	    return (0.0);
    } else
	return (seval(r, tgsp->radius, sig2, sig2 + n + 1, n + 1));
}

/*
 * CALC_SIG2_GSP: calculate radial dispersion needed to support the
 * density profile testprof in the potential of profile massprof, and
 * return the resulting table (testprof->npoint elements) and spline
 * coefficients (3 * testprof->npoint elements).
 */

real diffstep(real *, real *, int, void (*)(real *, real *), real);

local void sigdiff(real *, real *);

local gsprof *tgsp, *mgsp;		/* density, total mass profiles     */

local real beta_a;			/* anisotropy parameter             */

real *calc_sig2_gsp(gsprof *testprof, gsprof *massprof, real beta_aniso)
{
    real *sig2, rsig[2], avgerr, maxerr, locerr;
    int n, i;

    tgsp = testprof;
    mgsp = massprof;
    beta_a = beta_aniso;
    sig2 = (real *) allocate(4 * tgsp->npoint * sizeof(real));
    n = tgsp->npoint - 1;
    rsig[0] = tgsp->radius[n];
    if (tgsp->density[n] > 0.0)
	rsig[1] = - mgsp->mtot / ((2*beta_a + tgsp->beta - 1) * rsig[0]) +
		    (mgsp->mtot - mass_gsp(mgsp, rsig[0])) /
		      ((2 * beta_a + mgsp->beta + tgsp->beta + 2) * rsig[0]);
    else
	rsig[1] = 0.0;
    sig2[n] = rsig[1];
    avgerr = maxerr = 0.0;
    for (i = n - 1; i >= 0; i--) {
	locerr = diffstep(rsig, rsig, 2, sigdiff, tgsp->radius[i] - rsig[0]);
	sig2[i] = rsig[1];
	avgerr = avgerr + locerr / n;
	maxerr = MAX(maxerr, locerr);
    }
    eprintf("[%s.calc_sig2_gsp:  avgerr = %f  maxerr = %f]\n",
	    getargv0(), avgerr, maxerr);
    spline(sig2 + n + 1, tgsp->radius, sig2, n + 1);
    return (sig2);
}

local void sigdiff(real *drsig, real *rsig)
{
    drsig[0] = 1.0;
    drsig[1] = 0.0;
    if (rsig[0] > 0)
	drsig[1] -= mass_gsp(mgsp, rsig[0]) / rsqr(rsig[0]);
    if (rsig[1] > 0 && rsig[0] > 0)
	drsig[1] -= rsig[1] * (2*beta_a / rsig[0] + drho_gsp(tgsp, rsig[0]) /
						      rho_gsp(tgsp, rsig[0]));
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Calculate radial dispersion for GSP",
    "gsp=???",			";Input GSP for density distribution",
    "grav=",			";Input GSP for potential computation",
    "beta_a=0.0",		";Anisotropy parameter: beta_a <= 1",
    "npoint=27",		";Number of points tabulated",
    "rrange=1/128:64",		";Range of radii tabulated",
    "VERSION=1.3",		";Josh Barnes  19 December 1998",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr;
    gsprof *tgsp, *mgsp;
    real beta_a, *sig2, rrange[2], lgrs, r;
    int np, i;

    initparam(argv, defv);
    istr = stropen(getparam("gsp"), "r");
    get_history(istr);
    tgsp = get_gsprof(istr);
    strclose(istr);
    if (! strnull(getparam("grav"))) {
	istr = stropen(getparam("grav"), "r");
	get_history(istr);
	mgsp = get_gsprof(istr);
	strclose(istr);
    } else
	mgsp = tgsp;
    beta_a = getdparam("beta_a");
    sig2 = calc_sig2_gsp(tgsp, mgsp, beta_a);
    np = getiparam("npoint");
    setrange(rrange, getparam("rrange"));
    lgrs = rlog2(rrange[1] / rrange[0]) / (np - 1);
    printf("%12s  %12s\n", "radius", "sig_r^2");
    for (i = 0; i < np; i++) {
	r = rrange[0] * rpow(2.0, lgrs * i);
	printf("%12.5f  %12.7f\n",
	       r, sig2_gsp(tgsp, mgsp, beta_a, sig2, r));
    }
    return (0);
}

#endif
