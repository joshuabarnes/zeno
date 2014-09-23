/*
 * EXPDFSP.C: generate profile tables for spherical configuration with
 * same cumulative mass profile as an exponential disk.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "assert.h"
#include "fsp.h"

/*
 * EXPDFSP: initialize tables for "spherical exponential disk".
 */

fsprof *expdfsp(real mtot, real alpha, int np, real rmin, real rmax)
{
    fsprof *fsp;
    real *rtab, *dtab, *mtab, lgrs;
    int i;

    assert(mtot > 0.0 && alpha > 0.0);
    fsp = (fsprof *) allocate(sizeof(fsprof));
    rtab = (real *) allocate(np * sizeof(real));
    dtab = (real *) allocate(np * sizeof(real));
    mtab = (real *) allocate(np * sizeof(real));
    lgrs = rlog2(rmax / rmin) / (np - 1);
    for (i = 0; i < np; i++) {
	rtab[i] = rmin * rexp2(lgrs * i);
	dtab[i] = mtot * (rsqr(alpha) / (4 * PI * rtab[i])) *
		    rexp(- alpha * rtab[i]);
	mtab[i] = mtot * (1 - rexp(- alpha * rtab[i]) *
			    (1 + alpha * rtab[i]));
    }
    fsp->npoint = np;
    fsp->radius = rtab;
    fsp->density = dtab;
    fsp->mass = mtab;
    fsp->alpha = -1.0;
    fsp->beta = - (1 + alpha * rmax);
    fsp->mtot = mtot;
    return (fsp);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Generate \"spherical\" exp. disk",
    "out=???",			";Output file for profile tables",
    "mtot=1.0",			";Total mass of disk",
    "alpha=1.0",		";Inverse length scale of disk",
    "npoint=257",		";Number of points in tables",
    "rrange=1/4096:16",		";Range of radii tabulated",
    "VERSION=1.0",		";Josh Barnes  4 February 1995",
    NULL,
};

int main(int argc, string argv[])
{
    real rrange[2];
    fsprof *fsp;
    stream ostr;

    initparam(argv, defv);
    setrange(rrange, getparam("rrange"));
    fsp = expdfsp(getdparam("mtot"), getdparam("alpha"),
		  getiparam("npoint"), rrange[0], rrange[1]);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_fsprof(ostr, fsp);
    strclose(ostr);
    return (0);
}

#endif
