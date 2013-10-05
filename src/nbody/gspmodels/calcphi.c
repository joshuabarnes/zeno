/*
 * CALCPHI.C: calculate gravitational potential of GSP.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

#define NITER  3				/* iterations for r_phi_gsp */

/*
 * PHI_GSP: evaluate potential at given radius.
 */

real phi_gsp(gsprof *gsp, real r)
{
    int n = gsp->npoint - 1;

    if (gsp->phi == NULL)
        calc_phi_gsp(gsp, getargv0());
    if (r < gsp->radius[0])
	return (gsp->phi[0] - (gsp->mass[0] / gsp->radius[0]) *
		  (gsp->alpha == -2 ? rlog(gsp->radius[0] / r) :
		     (1 - rpow(r / gsp->radius[0], 2 + gsp->alpha)) /
		       (2 + gsp->alpha)));
    else if (r > gsp->radius[n])
	return (- (gsp->mtot + (gsp->mtot - gsp->mass[n]) *
		     rpow(r / gsp->radius[n], 3 + gsp->beta) /
		       (2 + gsp->beta)) / r);
    return (seval(r, gsp->radius, gsp->phi, gsp->pr_coef, n + 1));
}

/*
 * R_PHI_GSP: evaluate radius corresponding to given potential.
 */

real r_phi_gsp(gsprof *gsp, real p)
{
    int n = gsp->npoint - 1, k, i;
    real x, r;

    if (gsp->phi == NULL)
        calc_phi_gsp(gsp, getargv0());
    if (p < gsp->phi[0]) {
        if (gsp->alpha > -2 &&
	      p < gsp->phi[0] - gsp->mass[0] / (gsp->radius[0]*(2+gsp->alpha)))
	    error("%s.r_phi_gsp: phi < phi(0)\n", getargv0());
	x = (gsp->radius[0] / gsp->mass[0]) * (p - gsp->phi[0]);
	return ((gsp->alpha == -2 ? rexp(x) :
		   rpow(1 + (2 + gsp->alpha) * x, 1 / (2 + gsp->alpha))) *
		  gsp->radius[0]);
    } else if (p > gsp->phi[n]) {
        r = - gsp->mtot / p;			/* make initial guess       */
	for (i = 0; i < NITER; i++)
	    r = - (gsp->mtot + rpow(r / gsp->radius[n], 3 + gsp->beta) *
		     (gsp->mtot - gsp->mass[n]) / (2 + gsp->beta)) / p;
	return (r);				/* return refined estimate  */
    }
    for (k = 0; k <= n && gsp->phi[k] == gsp->phi[k+1]; k++);
    return (seval(p, gsp->phi + k, gsp->radius + k, gsp->rp_coef, n + 1 - k));
}

/*
 * Declarations for potential calculation.
 */

real diffstep(real *, real *, int, void (*)(real *, real *), real);

local void phidiff(real *, real *);

local gsprof *gsp;

/*
 * CALC_PHI_GSP: calculate gravitational potential of the given profile.
 */

void calc_phi_gsp(gsprof *gsp1, string trace)
{
    int n = gsp1->npoint - 1, i, k;
    real rphi[2], avgerr, maxerr, locerr;

    gsp = gsp1;
    gsp->phi = (real *) allocate((n + 1) * sizeof(real));
    rphi[0] = gsp->radius[n];
    rphi[1] = - (gsp->mtot + (gsp->mtot - gsp->mass[n]) /
		               (2 + gsp->beta)) / rphi[0];
    gsp->phi[n] = rphi[1];
    avgerr = maxerr = 0.0;
    for (i = n - 1; i >= 0; i--) {
	locerr = diffstep(rphi, rphi, 2, phidiff, gsp->radius[i] - rphi[0]);
	gsp->phi[i] = rphi[1];
	avgerr = avgerr + locerr / n;
	maxerr = MAX(maxerr, locerr);
    }
    if (trace != NULL)
        eprintf("[%s.calc_phi_gsp:  avgerr = %f  maxerr = %f]\n",
		trace, avgerr, maxerr);
    gsp->pr_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
    gsp->rp_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
    spline(gsp->pr_coef, gsp->radius, gsp->phi, n + 1);
    for (k = 0; k <= n && gsp->phi[k] == gsp->phi[k+1]; k++);
    if (k > 0 && trace != NULL)
        eprintf("[%s.calc_phi_gsp: skipping %d values for r_phi_gsp]\n",
		trace, k);
    spline(gsp->rp_coef, gsp->phi + k, gsp->radius + k, n + 1 - k);
}

local void phidiff(real *drphi, real *rphi)
{
    drphi[0] = 1.0;
    drphi[1] = (rphi[0] > 0 ? mass_gsp(gsp, rphi[0]) / rsqr(rphi[0]) : 0.0);
}

#ifdef TESTBED

string defv[] = {		";Calculate potential of GSP",
    "gsp=???",			";Input file with GSP",
    "npoint=27",		";Number of points tabulated",
    "rrange=1/128:64",		";Range of radii tabulated",
    "VERSION=1.4",		";Josh Barnes  20 November 2000",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr;
    gsprof *gsp;
    int np, i;
    real rrange[2], lgrs, r;

    initparam(argv, defv);
    istr = stropen(getparam("gsp"), "r");
    get_history(istr);
    gsp = get_gsprof(istr);
    strclose(istr);
    (void) phi_gsp(gsp, rrange[0]);		/* precalculate phi table   */
    np = getiparam("npoint");
    setrange(rrange, getparam("rrange"));
    lgrs = rlog2(rrange[1] / rrange[0]) / (np - 1);
    printf("%12s  %12s  %12s\n", "radius", "phi(r)", "r(phi)/r");
    for (i = 0; i < np; i++) {
	r = rrange[0] * rpow(2.0, lgrs * i);
	printf("%12.7f  %12.7f  %12.7f\n",
	       r, phi_gsp(gsp, r), r_phi_gsp(gsp, phi_gsp(gsp, r)) / r);
    }
    return (0);
}

#endif
