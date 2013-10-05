/*
 * PLUMSNAP: create Plummer model.
 */
 
#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
 
string defv[] = {		";Construct Plummer model",
    "out=???",			";Snapshot output file",
    "nbody=4096",		";Number of bodies to generate",
    "seed=54321",		";Random number seed",
    "nmodel=1",			";Number of copies to generate",
    "besort=true",		";Sort particles by binding energy",
    "zerocm=false",		";Transform to center of mass coords",
    "VERSION=1.0",		";Josh Barnes	12 May 2012",
    NULL,
};

string bodyfields[] = { PosTag, VelTag, MassTag, PhiTag, NULL };

bodyptr btab;			/* body array generated below */
int nbody;			/* number of bodies in array */

void plummodel();
int berank(const void *p1, const void *p2);

int main(int argc, string argv[])
{
  stream outstr;
  int nmodel;
  real tzero = 0.0;

  initparam(argv, defv);
  layout_body(bodyfields, Precision, NDIM);
  nbody = getiparam("nbody");
  btab = (bodyptr) allocate(nbody * SizeofBody);
  init_random(getiparam("seed"));
  nmodel = getiparam("nmodel");
  outstr = stropen(getparam("out"), "w");
  put_history(outstr);
  while (--nmodel >= 0) {
    plummodel();
    if (getbparam("besort"))
      qsort(btab, nbody, SizeofBody, berank);
    if (getbparam("zerocm"))
      snapcenter(btab, nbody, MassField.offset);
    put_snap(outstr, &btab, &nbody, &tzero, bodyfields);
  }
  strclose(outstr);
  return (0);
}

/*
 * PLUMMODEL: generate N-body realization of Plummer model,
 * scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 */

#define MFRAC  0.999				/* cut off 1-MFRAC of mass  */

void plummodel(void)
{
  real rsc, vsc, r, v, x, y;
  bodyptr p;

  if (nbody < 1)				/* check for silly values   */
    error("%s: absurd value for nbody\n", getargv0());
  rsc = (3 * PI) / 16;				/* set length scale factor  */
  vsc = rsqrt(1.0 / rsc);			/* find speed scale factor  */
  for (p = btab; p < NthBody(btab,nbody); p = NextBody(p)) {
						/* loop over particles      */
    Mass(p) = 1.0 / nbody;			/* set masses equal         */
    x = xrandom(0.0, MFRAC);			/* pick enclosed mass       */
    r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1);	/* find enclosing radius    */
    pickshell(Pos(p), NDIM, rsc * r);		/* pick position vector     */
    Phi(p) = -1.0 / (rsc * rsqrt(rsqr(r) + 1));	/* compute model potential  */
    do {					/* select from fn g(x)      */
      x = xrandom(0.0, 1.0);			/* for x in range 0:1       */
      y = xrandom(0.0, 0.1);			/* max of g(x) is 0.092     */
    } while (y > x*x * rpow(1 - x*x, 3.5));	/* using von Neumann tech   */
    v = x * rsqrt(2.0 / rsqrt(1 + r*r));	/* find resulting speed     */
    pickshell(Vel(p), NDIM, vsc * v);		/* pick velocity vector     */
  }
}

int berank(const void *p1, const void *p2)
{
  real e1, e2;

  e1 = 0.5 * rsqr(absv(Vel((bodyptr) p1))) + Phi((bodyptr) p1);
  e2 = 0.5 * rsqr(absv(Vel((bodyptr) p2))) + Phi((bodyptr) p2);
  return (e1 < e2 ? -1 : e1 > e2 ? 1 : 0);
}
