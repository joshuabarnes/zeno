/*
 * PLUMGSP.C: generate profile tables for Plummer model.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "assert.h"
#include "gsp.h"

/*
 * PLUMGSP: initialize tables for Plummer model.
 */

gsprof *plumgsp(real mtot, real a, int np, real rmin, real rmax)
{
    gsprof *gsp;
    real *rtab, *dtab, *mtab, lgrs;
    int i;

    assert(mtot > 0.0 && a > 0.0);
    gsp = (gsprof *) allocate(sizeof(gsprof));
    rtab = (real *) allocate(np * sizeof(real));
    dtab = (real *) allocate(np * sizeof(real));
    mtab = (real *) allocate(np * sizeof(real));
    lgrs = rlog2(rmax / rmin) / (np - 1);
    for (i = 0; i < np; i++) {
	rtab[i] = rmin * rexp2(lgrs * i);
	dtab[i] = (3 / FOUR_PI) * mtot * a*a /
		    rpow(rsqr(rtab[i]) + a*a, 2.5);
	mtab[i] = mtot * rqbe(rtab[i]) /
		    rpow(rsqr(rtab[i]) + a*a, 1.5);
    }
    gsp->npoint = np;
    gsp->radius = rtab;
    gsp->density = dtab;
    gsp->mass = mtab;
    gsp->alpha = 0.0;
    gsp->beta = -5.0;
    gsp->mtot = mtot;
    return (gsp);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Generate profile for Plummer model",
    "out=???",			";Output file for profile tables",
    "mtot=1.0",			";Total mass of model",
    "a=1.0",			";Length scale of model",
    "npoint=257",		";Number of points in tables",
    "rrange=1/4096:16",		";Range of radii tabulated",
    "VERSION=1.2",		";Josh Barnes  4 February 1995",
    NULL,
};

int main(int argc, string argv[])
{
    real rrange[2];
    gsprof *gsp;
    stream ostr;

    initparam(argv, defv);
    setrange(rrange, getparam("rrange"));
    gsp = plumgsp(getdparam("mtot"), getdparam("a"),
		  getiparam("npoint"), rrange[0], rrange[1]);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_gsprof(ostr, gsp);
    strclose(ostr);
    return (0);
}

#endif
