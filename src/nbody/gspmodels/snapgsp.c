/*
 * SNAPGSP.C: Construct GSP from N-body snapshot.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "gsp.h"

string defv[] = {		";Construct GSP from N-body snapshot",
    "in=???",			";Input snapshot file",
    "out=",			";Output GSP file",
    "npoint=256",		";Number of data points in GSP.",
				";Must divide nbody evenly.",
    "alpha=0.0",		";Power-law index at small r",
    "beta=0.0",			";Power-law index at large r",
    "VERSION=1.0",		";Josh Barnes  18 October 2004",
    NULL,
};

string bodyfields[] = { PosTag, MassTag, NULL };

gsprof *snapgsp(bodyptr btab, int nbody, int npoint, real alpha, real beta);
int radrank(const void *, const void *);

int main(int argc, string argv[])
{
  stream istr, ostr;
  bodyptr btab = NULL;
  int nbody;
  real tnow;
  string intags[MaxBodyFields];
  gsprof *gsp;

  initparam(argv, defv);
  layout_body(bodyfields, Precision, NDIM);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab, &nbody, &tnow, intags, FALSE))
    error("%s: snapshot input failed\n", getargv0());
  if (! set_member(intags, PosTag))
    error("%s: position data missing\n", getargv0());
  if (! set_member(intags, MassTag))
    error("%s: mass data missing\n", getargv0());
  gsp = snapgsp(btab, nbody, getiparam("npoint"),
		getdparam("alpha"), getdparam("beta"));
  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_gsprof(ostr, gsp);
    strclose(ostr);
  }
  return (0);
}

gsprof *snapgsp(bodyptr btab, int nbody, int npoint, real alpha, real beta)
{
  gsprof *gsp;
  real *rtab, *dtab, *mtab, *coef, mtot;
  int nsamp, i, j;

  if (nbody % npoint != 0)
    error("%s: npoint must divide nbody\n", getargv0());
  gsp = (gsprof *) allocate(sizeof(gsprof));
  rtab = (real *) allocate(npoint * sizeof(real));
  dtab = (real *) allocate(npoint * sizeof(real));
  mtab = (real *) allocate(npoint * sizeof(real));
  coef = (real *) allocate(3 * npoint * sizeof(real));
  qsort(btab, nbody, SizeofBody, radrank);
  nsamp = nbody / npoint;
  mtot = 0.0;
  j = 0;
  for (i = 0; i < nbody; i++) {
    mtot += Mass(NthBody(btab, i));
    if (i % nsamp == nsamp - 1) {
      rtab[j] = absv(Pos(NthBody(btab, i)));
      mtab[j] = mtot;
      j++;
    }
  }
  spline(coef, rtab, mtab, npoint);
  for (j = 0; j < npoint; j++)
    dtab[j] =
      spldif(rtab[j], rtab, mtab, coef, npoint) / (FOUR_PI * rsqr(rtab[j]));
  gsp->npoint = npoint;
  gsp->radius = rtab;
  gsp->density = dtab;
  gsp->alpha = alpha;
  gsp->beta = beta;
  gsp->mass = mtab;
  gsp->mtot = mtot;
  return (gsp);
}

int radrank(const void *a, const void *b)
{
  real Ra, Rb;

  Ra = absv(Pos((bodyptr) a));
  Rb = absv(Pos((bodyptr) b));
  return (Ra < Rb ? -1 : Ra > Rb ? 1 : 0);
}
