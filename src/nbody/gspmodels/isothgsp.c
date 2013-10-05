/*
 * ISOTHGSP.C: generate tables for non-singular truncated isothermal.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "assert.h"
#include "gsp.h"

real diffstep(real *, real *, int, void (*)(real *, real *), real);

local void mdiff(real *drm, real *rm);	/* compute dm/dr for mass profile   */

local gsprof *gsp;			/* for use in mdiff() function      */

/*
 * ISOTHGSP: non-singular truncated isothermal.
 */

gsprof *isothgsp(real m_h, real a, real b, int np, real rmin, real rmax)
{
    real *rtab, *dtab, *mtab, lgrs;
    real q, fnorm, r_i, rm[2], maxerr, locerr;
    int i;

    assert(m_h > 0.0 && a > 0.0 && a < b);
    gsp = (gsprof *) allocate(sizeof(gsprof));
    gsp->npoint = np;
    gsp->radius =  rtab = (real *) allocate(np * sizeof(real));
    gsp->density = dtab = (real *) allocate(np * sizeof(real));
    gsp->mass =    mtab = (real *) allocate(np * sizeof(real));
    gsp->alpha = 0.0;
    gsp->beta = -2 * rsqr(rmax) * (1 / rsqr(b) + 1 / (rsqr(rmax) + rsqr(a)));
    gsp->mtot = m_h;
    lgrs = rlog2(rmax / rmin) / (np - 1);
    q = a / b;
    fnorm = 1 / (1 - rsqrt(PI) * q * rexp(q * q) * erfc(q));
    eprintf("[q = %f  fnorm = %f  beta = %f]\n", q, fnorm, gsp->beta);
    for (i = 0; i < np; i++) {
	rtab[i] = r_i = rmin * rexp2(lgrs * i);
	dtab[i] = (fnorm * m_h / (2 * PI * rsqrt(PI) * b)) *
		    rexp(- rsqr(r_i / b)) / (rsqr(r_i) + rsqr(a));
    }
    mtab[0] = FRTHRD_PI * rqbe(rtab[0]) * dtab[0];
    rm[0] = rtab[0];
    rm[1] = mtab[0];
    for (i = 1; i < np; i++) {
        locerr = diffstep(rm, rm, 2, mdiff, rtab[i] - rm[0]);
	mtab[i] = rm[1];
	maxerr = (i == 1 ? locerr : MAX(locerr, maxerr));
    }
    eprintf("[maxerr = %f  mass(rmax) = %f]\n", maxerr, mtab[np - 1]);
    return (gsp);
}

/*
 * MDIFF: compute dm/dr, assuming density given by existing gsp.
 */

local void mdiff(real *drm, real *rm)
{
    drm[0] = 1.0;
    drm[1] = 4 * PI * rsqr(rm[0]) * rho_gsp(gsp, rm[0]);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Generate isothermal halo profile",
    "out=???",			";Output file for profile tables",
    "mtot=1.0",			";Total mass of halo",
    "a=1.0",			";Radius of core",
    "b=10.0",			";Radius of taper",
    "npoint=257",		";Number of points in tables",
    "rrange=1/4096:16",		";Range of radii tabulated",
    "VERSION=1.0",		";Josh Barnes  9 November 2000",
    NULL,
};

int main(int argc, string argv[])
{
    real rrange[2];
    gsprof *gsp;
    stream ostr;

    initparam(argv, defv);
    setrange(rrange, getparam("rrange"));
    gsp = isothgsp(getdparam("mtot"), getdparam("a"), getdparam("b"),
		   getiparam("npoint"), rrange[0], rrange[1]);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_gsprof(ostr, gsp);
    strclose(ostr);
    return (0);
}

#endif
