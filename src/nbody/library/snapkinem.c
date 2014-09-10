/*
 * snapkinem.c: routines to compute kinematic parameters.
 */

#include "stdinc.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapkinem.h"

//  Weight: macro to select weight field for kinematic calculations.
//  ________________________________________________________________

#define Weight(bp,woff)  SelectReal(bp, woff)

//  snapke: return total kinetic energy.
//  ____________________________________

real snapke(bodyptr btab, int nbody, int woff)
{
  real ke;
  bodyptr bp;

  ke = 0.0;
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
    ke += 0.5 * Weight(bp, woff) * dotvp(Vel(bp), Vel(bp));
  return (ke);
}

//  snappe: return total potential energy.
//  ______________________________________

real snappe(bodyptr btab, int nbody, int woff)
{
  real pe;
  bodyptr bp;

  pe = 0.0;
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
    pe += Weight(bp, woff) * Phi(bp);
  return (pe);
}

//  snapamvec: compute angular momentum vector.
//  ___________________________________________

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

//  snapketen: compute kinetic energy tensor.
//  _________________________________________

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

//  snappeten: compute potential energy tensor.
//  ___________________________________________

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

//  snapmiten: compute moment of inertia tensor.
//  ____________________________________________

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
