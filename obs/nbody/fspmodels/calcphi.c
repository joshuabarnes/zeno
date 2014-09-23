/*
 * CALCPHI.C: calculate gravitational potential of FSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "fsp.h"

local void calc_phi_fsp(fsprof *);		/* solve for potential      */

#define NITER  3				/* iterations for r_phi_fsp */

/*
 * PHI_FSP: evaluate potential at given radius.
 */

real phi_fsp(fsprof *fsp, real r)
{
    int n = fsp->npoint - 1;

    if (fsp->phi == NULL)
	calc_phi_fsp(fsp);
    if (r < fsp->radius[0])
	return (fsp->phi[0] - (fsp->mass[0] / fsp->radius[0]) *
		  (fsp->alpha == -2 ? rlog(fsp->radius[0] / r) :
		     (1 - rpow(r / fsp->radius[0], 2 + fsp->alpha)) /
		       (2 + fsp->alpha)));
    else if (r > fsp->radius[n])
	return (- (fsp->mtot + (fsp->mtot - fsp->mass[n]) *
		     rpow(r / fsp->radius[n], 3 + fsp->beta) /
		       (2 + fsp->beta)) / r);
    return (seval(r, fsp->radius, fsp->phi, fsp->pr_coef, n + 1));
}

/*
 * R_PHI_FSP: evaluate radius corresponding to given potential.
 */

real r_phi_fsp(fsprof *fsp, real p)
{
    int n = fsp->npoint - 1, k, i;
    real x, r;

    if (fsp->phi == NULL)
	calc_phi_fsp(fsp);
    if (p < fsp->phi[0]) {
        if (fsp->alpha > -2 &&
	      p < fsp->phi[0] - fsp->mass[0] / (fsp->radius[0]*(2+fsp->alpha)))
	    error("%s.r_phi_fsp: phi < phi(0)\n", getargv0());
	x = (fsp->radius[0] / fsp->mass[0]) * (p - fsp->phi[0]);
	return ((fsp->alpha == -2 ? rexp(x) :
		   rpow(1 + (2 + fsp->alpha) * x, 1 / (2 + fsp->alpha))) *
		  fsp->radius[0]);
    } else if (p > fsp->phi[n]) {
        r = - fsp->mtot / p;			/* make initial guess       */
	for (i = 0; i < NITER; i++)
	    r = - (fsp->mtot + rpow(r / fsp->radius[n], 3 + fsp->beta) *
		     (fsp->mtot - fsp->mass[n]) / (2 + fsp->beta)) / p;
	return (r);				/* return refined estimate  */
    }
    for (k = 0; k <= n && fsp->phi[k] == fsp->phi[k+1]; k++);
    return (seval(p, fsp->phi + k, fsp->radius + k, fsp->rp_coef, n + 1 - k));
}

/*
 * Declarations for potential calculation.
 */

real diffstep(real *, real *, int, void (*)(real *, real *), real);

local void phidiff(real *, real *);

local fsprof *fsp;

/*
 * CALC_PHI_FSP: calculate gravitational potential of the given profile.
 */

void calc_phi_fsp(fsprof *fsp1)
{
    int n = fsp1->npoint - 1, i, k;
    real rphi[2], avgerr, maxerr, locerr;

    fsp = fsp1;
    fsp->phi = (real *) allocate((n + 1) * sizeof(real));
    rphi[0] = fsp->radius[n];
    rphi[1] = - (fsp->mtot + (fsp->mtot - fsp->mass[n]) /
		               (2 + fsp->beta)) / rphi[0];
    fsp->phi[n] = rphi[1];
    avgerr = maxerr = 0.0;
    for (i = n - 1; i >= 0; i--) {
	locerr = diffstep(rphi, rphi, 2, phidiff, fsp->radius[i] - rphi[0]);
	fsp->phi[i] = rphi[1];
	avgerr = avgerr + locerr / n;
	maxerr = MAX(maxerr, locerr);
    }
    eprintf("[%s.calc_phi_fsp:  avgerr = %f  maxerr = %f]\n",
	    getargv0(), avgerr, maxerr);
    fsp->pr_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
    fsp->rp_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
    spline(fsp->pr_coef, fsp->radius, fsp->phi, n + 1);
    for (k = 0; k <= n && fsp->phi[k] == fsp->phi[k+1]; k++);
    if (k > 0)
        eprintf("[%s.calc_phi_fsp: skipping %d values for r_phi_fsp]\n",
		getargv0(), k);
    spline(fsp->rp_coef, fsp->phi + k, fsp->radius + k, n + 1 - k);
}

local void phidiff(real *drphi, real *rphi)
{
    drphi[0] = 1.0;
    drphi[1] = (rphi[0] > 0 ? mass_fsp(fsp, rphi[0]) / rsqr(rphi[0]) : 0.0);
}

#ifdef TESTBED

string defv[] = {		";Calculate potential of FSP",
    "fsp=???",			";Input file with FSP",
    "npoint=27",		";Number of points tabulated",
    "rrange=1/128:64",		";Range of radii tabulated",
    "VERSION=1.4",		";Josh Barnes  20 November 2000",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr;
    fsprof *fsp;
    int np, i;
    real rrange[2], lgrs, r;

    initparam(argv, defv);
    istr = stropen(getparam("fsp"), "r");
    get_history(istr);
    fsp = get_fsprof(istr);
    strclose(istr);
    (void) phi_fsp(fsp, rrange[0]);		/* precalculate phi table   */
    np = getiparam("npoint");
    setrange(rrange, getparam("rrange"));
    lgrs = rlog2(rrange[1] / rrange[0]) / (np - 1);
    printf("%12s  %12s  %12s\n", "radius", "phi(r)", "r(phi)/r");
    for (i = 0; i < np; i++) {
	r = rrange[0] * rpow(2.0, lgrs * i);
	printf("%12.7f  %12.7f  %12.7f\n",
	       r, phi_fsp(fsp, r), r_phi_fsp(fsp, phi_fsp(fsp, r)) / r);
    }
    return (0);
}

#endif
