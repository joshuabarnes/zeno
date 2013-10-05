/*
 * SNAPKINEM.C: routines to compute kinematic parameters.
 */

#include "stdinc.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapkinem.h"

/*
 * WEIGHT: macro to select weight field for kinematic calculations.
 */

#define Weight(bp,woff)  SelectReal(bp, woff)

/*
 * SNAPKE: return total kinetic energy.
 */

real snapke(bodyptr btab, int nbody, int woff)
{
    real ke;
    bodyptr bp;

    ke = 0.0;
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
        ke += 0.5 * Weight(bp, woff) * dotvp(Vel(bp), Vel(bp));
    return (ke);
}

/*
 * SNAPPE: return total potential energy.
 */

real snappe(bodyptr btab, int nbody, int woff)
{
    real pe;
    bodyptr bp;

    pe = 0.0;
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
        pe += Weight(bp, woff) * Phi(bp);
    return (pe);
}

/*
 * SNAPAMVEC: compute angular momentum vector.
 */

void snapamvec(vector amvec, bodyptr btab, int nbody, int woff)
{
    bodyptr bp;
    vector tmpv;

    CLRV(amvec);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	CROSSVP(tmpv, Vel(bp), Pos(bp));
	ADDMULVS(amvec, tmpv, Weight(bp, woff));
    }
}

/*
 * SNAPKETEN: compute kinetic energy tensor.
 */

void snapketen(matrix keten, bodyptr btab, int nbody, int woff)
{
    bodyptr bp;
    vector vtmp;
    matrix mtmp;

    CLRM(keten);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
        MULVS(vtmp, Vel(bp), 0.5 * Weight(bp, woff));
	OUTVP(mtmp, Vel(bp), vtmp);
	ADDM(keten, keten, mtmp);
    }
}

/*
 * SNAPPETEN: compute potential energy tensor.
 */

void snappeten(matrix peten, bodyptr btab, int nbody, int woff)
{
    bodyptr bp;
    vector vtmp;
    matrix mtmp;

    CLRM(peten);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
        MULVS(vtmp, Pos(bp), Weight(bp, woff));
	OUTVP(mtmp, Acc(bp), vtmp);
	ADDM(peten, peten, mtmp);
    }
}

/*
 * SNAPMITEN: compute moment of inertia tensor.
 */

void snapmiten(matrix miten, bodyptr btab, int nbody, int woff)
{
    bodyptr bp;
    vector vtmp;
    matrix mtmp;

    CLRM(miten);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
        MULVS(vtmp, Pos(bp), Weight(bp, woff));
	OUTVP(mtmp, Pos(bp), vtmp);
	ADDM(miten, miten, mtmp);
    }
}
