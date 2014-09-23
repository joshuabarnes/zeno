/*
 * GAMMAFSP.C: generate profile tables for Dehnen models.
 * See Dehnen, W. 1993, MNRAS, 265, 250.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "fsp.h"

/*
 * GAMMAFSP: initialize tables for Gamma model.
 */

fsprof *gammafsp(real gam, real mtot, real a, int np, real rmin, real rmax)
{
    fsprof *fsp;
    real *rtab, *dtab, *mtab, lgrs;
    int i;

    fsp = (fsprof *) allocate(sizeof(fsprof));
    rtab = (real *) allocate(np * sizeof(real));
    dtab = (real *) allocate(np * sizeof(real));
    mtab = (real *) allocate(np * sizeof(real));
    lgrs = rlog2(rmax / rmin) / (np - 1);
    for (i = 0; i < np; i++) {
	rtab[i] = rmin * rexp2(lgrs * i);
	dtab[i] = ((3 - gam) / FOUR_PI) * mtot * a /
		    (rpow(rtab[i], gam) * rpow(rtab[i] + a, 4 - gam));
	mtab[i] = mtot * rpow(rtab[i] / (rtab[i] + a), 3 - gam);
    }
    fsp->npoint = np;
    fsp->radius = rtab;
    fsp->density = dtab;
    fsp->mass = mtab;
    fsp->alpha = -gam;
    fsp->beta = -4;
    fsp->mtot = mtot;
    return (fsp);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Generate profile for gamma model",
    "out=???",			";Output file for profile tables",
    "gamma=1.0",		";Structure parameter: 0 < gamma < 3",
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
    fsprof *fsp;
    stream ostr;

    initparam(argv, defv);
    setrange(rrange, getparam("rrange"));
    fsp = gammafsp(getdparam("gamma"), getdparam("mtot"), getdparam("a"),
		   getiparam("npoint"), rrange[0], rrange[1]);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_fsprof(ostr, fsp);
    strclose(ostr);
    return (0);
}

#endif
