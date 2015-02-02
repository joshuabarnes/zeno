/*
 * POLYGSP.C: generate profile tables for polytropic spheres.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

real diffstep(real *, real *, int, void (*)(real *, real *), real);

local real indx;			/* index of polytropic gas sphere   */

local void polydiff(real *, real *);	/* compute derivatives for model    */

#define MSTEP  2048			/* maximum integration steps	    */

/*
 * POLYDIFF: compute derivatives for polytropic sphere.
 */

local void polydiff(real *drdm, real *rdm)
{
    drdm[0] = 1.0;
    if (rdm[0] > 0.0 && rdm[1] > 0.0) {
	drdm[1] = - (indx / FOUR_PI) * rpow(rdm[1], 1.0 - 1.0/indx) *
		  rdm[2] / rsqr(rdm[0]);
	drdm[2] = FOUR_PI * rsqr(rdm[0]) * rdm[1];
    } else {
	drdm[1] = 0.0;
	drdm[2] = 0.0;
    }
}

/*
 * POLYGSP: initialize tables for polytropic spheres.
 */

gsprof *polygsp(real indx0, real mtot, real rsur, int np)
{
    gsprof *gsp;
    real dr, rdm[3], rtab[MSTEP], dtab[MSTEP], mtab[MSTEP];
    real r0, m0, dcoef[3*MSTEP], mcoef[3*MSTEP], *rptr, *dptr, *mptr;
    int n, i;

    if (! (0 < indx0 && indx0 <= 4.5))
	error("%s: must have 0 < index <= 4.5\n", getargv0());
    indx = indx0;
    rdm[0] = rtab[0] = 0.0;
    rdm[1] = dtab[0] = 1.0;
    rdm[2] = mtab[0] = 0.0;
    dr = (indx <= 4 ? 16.0 : 32.0) / MSTEP;
    for (n = 1; n < MSTEP && rdm[1] > 0; n++) {
	(void) diffstep(rdm, rdm, 3, polydiff, dr);
	rtab[n] = rdm[0];
	dtab[n] = rdm[1];
	mtab[n] = rdm[2];
    }
    if (rdm[1] > 0)
	error("%s: integration didn't reach surface\n", getargv0());
    r0 = rtab[n-1] = rtab[n-2] + dr * dtab[n-2] / (dtab[n-2] - dtab[n-1]);
    m0 = mtab[n-1];
    dtab[n-1] = 0.0;
    eprintf("[%s: reached surface at xi_n = %g]\n", getargv0(), rtab[n-1]);
    for (i = 0; i < n; i++) {
	rtab[i] = rtab[i] * (rsur / r0);
	mtab[i] = mtab[i] * (mtot / m0);
	dtab[i] = dtab[i] * (mtot / m0) / rqbe(rsur / r0);
    }
    spline(dcoef, rtab, dtab, n);
    spline(mcoef, rtab, mtab, n);
    gsp = (gsprof *) allocate(sizeof(gsprof));
    gsp->npoint = np;
    gsp->radius = rptr = (real *) allocate(np * sizeof(real));
    gsp->density = dptr = (real *) allocate(np * sizeof(real));
    gsp->mass = mptr = (real *) allocate(np * sizeof(real));
    for (i = 0; i < np; i++) {
	rptr[i] = rsur * i / (np - 1.0);
	dptr[i] = seval(rptr[i], rtab, dtab, dcoef, n);
	mptr[i] = seval(rptr[i], rtab, mtab, mcoef, n);
    }
    gsp->alpha = 0;
    gsp->beta = 0;
    gsp->mtot = mtot;
    return (gsp);
}

#ifdef TESTBED

string defv[] = {		";Generate profile for polytropic sphere",
    "out=",			";Output file for profile tables",
    "index=1.5",		";Structure parameter: 0 < index <= 4.5",
    "mtot=1.0",			";Total mass of model",
    "rsur=1.0",			";Surface radius of model",
    "npoint=65",		";Number of points in tables",
    "VERSION=1.0",		";Josh Barnes  18 December 1998",
    NULL,
};

int main(int argc, string argv[])
{
    gsprof *gsp;
    stream ostr;

    initparam(argv, defv);
    gsp = polygsp(getdparam("index"), getdparam("mtot"),
		  getdparam("rsur"), getiparam("npoint"));
    if (! strnull(getparam("out"))) {
	ostr = stropen(getparam("out"), "w");
	put_history(ostr);
	put_gsprof(ostr, gsp);
	strclose(ostr);
    }
    return (0);
}

#endif
