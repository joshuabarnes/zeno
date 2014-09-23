/*
 * CALCSIG.C: calculate radial velocity dispersion for FSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "fsp.h"

/*
 * SIG2_FSP: evaluate velocity dispersion at given radius.
 */

real sig2_fsp(fsprof *tfsp, fsprof *mfsp, real beta_a, real *sig2, real r)
{
    int n = tfsp->npoint - 1;
    real r0, rn, gamma, m0, c0;

    r0 = tfsp->radius[0];
    rn = tfsp->radius[n];
    if (r < r0) {
	gamma = 2 * beta_a + mfsp->alpha + tfsp->alpha + 2;
	m0 = mass_fsp(mfsp, r0);
	if (gamma != 0.0) {
	    c0 = sig2[0] + m0 / (gamma * r0);
	    return (- (m0 / r0) * rpow(r / r0, 2 + mfsp->alpha) / gamma +
		      c0 * rpow(r0 / r, 2 * beta_a + tfsp->alpha));
	} else {
	    c0 = sig2[0] + rlog(r0) * m0 / r0;
	    return (- rlog(r) * m0 / rpow(r0, 3 + mfsp->alpha) /
		      rpow(r, 2 * beta_a + tfsp->alpha) +
		        c0 * rpow(r0 / r, 2 * beta_a + tfsp->alpha));
	}
    } else if (r > rn) {
	if (tfsp->density[n] > 0) {
	    gamma = 2 * beta_a + mfsp->beta + tfsp->beta + 2;
	    return (- mfsp->mtot / ((2*beta_a + tfsp->beta - 1) * r) +
		      (mfsp->mtot - mass_fsp(mfsp, r)) / (gamma * r));
	} else
	    return (0.0);
    } else
	return (seval(r, tfsp->radius, sig2, sig2 + n + 1, n + 1));
}

/*
 * CALC_SIG2_FSP: calculate radial dispersion needed to support the
 * density profile testprof in the potential of profile massprof, and
 * return the resulting table (testprof->npoint elements) and spline
 * coefficients (3 * testprof->npoint elements).
 */

real diffstep(real *, real *, int, void (*)(real *, real *), real);

local void sigdiff(real *, real *);

local fsprof *tfsp, *mfsp;		/* density, total mass profiles     */

local real beta_a;			/* anisotropy parameter             */

real *calc_sig2_fsp(fsprof *testprof, fsprof *massprof, real beta_aniso)
{
    real *sig2, rsig[2], avgerr, maxerr, locerr;
    int n, i;

    tfsp = testprof;
    mfsp = massprof;
    beta_a = beta_aniso;
    sig2 = (real *) allocate(4 * tfsp->npoint * sizeof(real));
    n = tfsp->npoint - 1;
    rsig[0] = tfsp->radius[n];
    if (tfsp->density[n] > 0.0)
	rsig[1] = - mfsp->mtot / ((2*beta_a + tfsp->beta - 1) * rsig[0]) +
		    (mfsp->mtot - mass_fsp(mfsp, rsig[0])) /
		      ((2 * beta_a + mfsp->beta + tfsp->beta + 2) * rsig[0]);
    else
	rsig[1] = 0.0;
    sig2[n] = rsig[1];
    avgerr = maxerr = 0.0;
    for (i = n - 1; i >= 0; i--) {
	locerr = diffstep(rsig, rsig, 2, sigdiff, tfsp->radius[i] - rsig[0]);
	sig2[i] = rsig[1];
	avgerr = avgerr + locerr / n;
	maxerr = MAX(maxerr, locerr);
    }
    eprintf("[%s.calc_sig2_fsp:  avgerr = %f  maxerr = %f]\n",
	    getargv0(), avgerr, maxerr);
    spline(sig2 + n + 1, tfsp->radius, sig2, n + 1);
    return (sig2);
}

local void sigdiff(real *drsig, real *rsig)
{
    drsig[0] = 1.0;
    drsig[1] = 0.0;
    if (rsig[0] > 0)
	drsig[1] -= mass_fsp(mfsp, rsig[0]) / rsqr(rsig[0]);
    if (rsig[1] > 0 && rsig[0] > 0)
	drsig[1] -= rsig[1] * (2*beta_a / rsig[0] + drho_fsp(tfsp, rsig[0]) /
						      rho_fsp(tfsp, rsig[0]));
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Calculate radial dispersion for FSP",
    "fsp=???",			";Input FSP for density distribution",
    "grav=",			";Input FSP for potential computation",
    "beta_a=0.0",		";Anisotropy parameter: beta_a <= 1",
    "npoint=27",		";Number of points tabulated",
    "rrange=1/128:64",		";Range of radii tabulated",
    "VERSION=1.3",		";Josh Barnes  19 December 1998",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr;
    fsprof *tfsp, *mfsp;
    real beta_a, *sig2, rrange[2], lgrs, r;
    int np, i;

    initparam(argv, defv);
    istr = stropen(getparam("fsp"), "r");
    get_history(istr);
    tfsp = get_fsprof(istr);
    strclose(istr);
    if (! strnull(getparam("grav"))) {
	istr = stropen(getparam("grav"), "r");
	get_history(istr);
	mfsp = get_fsprof(istr);
	strclose(istr);
    } else
	mfsp = tfsp;
    beta_a = getdparam("beta_a");
    sig2 = calc_sig2_fsp(tfsp, mfsp, beta_a);
    np = getiparam("npoint");
    setrange(rrange, getparam("rrange"));
    lgrs = rlog2(rrange[1] / rrange[0]) / (np - 1);
    printf("%12s  %12s\n", "radius", "sig_r^2");
    for (i = 0; i < np; i++) {
	r = rrange[0] * rpow(2.0, lgrs * i);
	printf("%12.5f  %12.7f\n",
	       r, sig2_fsp(tfsp, mfsp, beta_a, sig2, r));
    }
    return (0);
}

#endif
