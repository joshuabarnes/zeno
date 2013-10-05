/*
 * SNAPCENTER.C: routines to compute center-of-mass coordinates.
 */

#include "stdinc.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"

/*
 * WEIGHT: macro to select weight field for c. of m. calculations.
 */

#define Weight(bp, woff)  SelectReal(bp, woff)

/*
 * SNAPCMPOS, SNAPCMVEL, SNAPCMACC: compute c. of m. position, velocity, ...
 */

void snapcmpos(vector cmpos, bodyptr btab, int nbody, int woff)
{
    double wtot, cmtmp[NDIM];
    bodyptr bp;

    wtot = 0.0;
    CLRV(cmtmp);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	wtot = wtot + Weight(bp, woff);
	ADDMULVS(cmtmp, Pos(bp), Weight(bp, woff));
    }
    DIVVS(cmpos, cmtmp, wtot);
}

void snapcmvel(vector cmvel, bodyptr btab, int nbody, int woff)
{
    double wtot, cmtmp[NDIM];
    bodyptr bp;

    wtot = 0.0;
    CLRV(cmtmp);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	wtot = wtot + Weight(bp, woff);
	ADDMULVS(cmtmp, Vel(bp), Weight(bp, woff));
    }
    DIVVS(cmvel, cmtmp, wtot);
}

void snapcmacc(vector cmacc, bodyptr btab, int nbody, int woff)
{
    double wtot, cmtmp[NDIM];
    bodyptr bp;

    wtot = 0.0;
    CLRV(cmtmp);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	wtot = wtot + Weight(bp, woff);
	ADDMULVS(cmtmp, Acc(bp), Weight(bp, woff));
    }
    DIVVS(cmacc, cmtmp, wtot);
}

/*
 * SNAPCENTER: transform to the weighted center of mass.
 */

void snapcenter(bodyptr btab, int nbody, int woff)
{
    vector cmpos, cmvel;
    bodyptr bp;

    snapcmpos(cmpos, btab, nbody, woff);
    snapcmvel(cmvel, btab, nbody, woff);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	SUBV(Pos(bp), Pos(bp), cmpos);
	SUBV(Vel(bp), Vel(bp), cmvel);
    }
}
