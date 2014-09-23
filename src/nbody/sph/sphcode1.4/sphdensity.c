/*
 * sphdensity.c: compute densities and smoothing lengths.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "datatypes.h"
#include "strset.h"
#include "gsp.h"

#define global					// prevent default to extern
#include "sphcode.h"
#include "kdtree.h"
#include "smooth.h"
#include "fixbody.h"

//  Command-line parameters and defaults.

string defv[] = {		";Compute densities and smoothing lengths",
  "in=???",			";Input, including " PosTag " data.",
				";" MassTag " data is optional.",
  "out=???",			";Output, with " RhoTag " and " SmoothTag
                                  " data",
  "rhomass=false",		";If TRUE, output particle mass density.",
				";Otherwise, output number density.",
  "nsmooth=40",			";Number of bodies in smoothing volume",
  "nbucket=16",			";Number of bodies in leaves of KD tree",
  "slope=0.0",			";Smoothing kernel slope at origin",
  "VERSION=1.0",		";Joshua Barnes  18 December 2012",
  NULL,
};

//  Local procedure prototypes and variables.

local void setupbody(void);			// set offsets for dynbody
local void sphdensity(bool, real, int, int);	// perform smoothing
local void sph_numb_density(smxptr, int, int);	// compute number density
local void sph_mass_density(smxptr, int, int);	// compute mass density

//  _______________________
//  main: toplevel routine.

int main(int argc, string argv[])
{
  stream istr, ostr = NULL;
  string intags[MaxBodyFields], *newtags, *outtags;
  bool massdata = FALSE;
  bodyptr p;

  initparam(argv, defv);			// initialize param access
  setupbody();					// interface to dynbody
  newtags = set_cons(SmoothTag, RhoTag, NULL);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  while (get_snap(istr, &btab, &nbody, &tnow, intags, FALSE)) {
    if (! set_member(intags, PosTag))
      error("%s: %s input data missing\n", getprog(), PosTag);
    massdata = massdata || set_member(intags, MassTag);
    if (getbparam("rhomass") && ! massdata)
      error("%s: %s input data required for mass density\n",
	    getprog(), MassTag);
    for (p = btab; p < btab+nbody; p++)
      Type(p) = BODY | GAS;
    sphdensity(getbparam("rhomass"), getdparam("slope"),
	       getiparam("nsmooth"), getiparam("nbucket"));
    if (ostr == NULL) {
      ostr = stropen(getparam("out"), "w");
      put_history(ostr);
    }
    outtags = set_union(intags, newtags);
    put_snap(ostr, &btab, &nbody, &tnow, outtags);
    free(outtags);
  }
  return (0);
}

//  ________________________________________________________________
//  setupbody: inform dynamic body routines of relevant body fields.

local void setupbody(void)
{
  define_body(sizeof(body), Precision, NDIM);
  define_body_offset(TypeTag, BodyOffset(Type));
  define_body_offset(PosTag, BodyOffset(Pos));
  define_body_offset(VelTag, BodyOffset(Vel));
  define_body_offset(MassTag, BodyOffset(Mass));
  define_body_offset(SmoothTag, BodyOffset(Smooth));
  define_body_offset(PhiTag, BodyOffset(Phi));
  define_body_offset(AccTag, BodyOffset(Acc));
  define_body_offset(RhoTag, BodyOffset(Rho));
}

//  _________________________________________________
//  sphdensity: compute smoothing length and density.

local void sphdensity(bool rhomass, real slope, int nsmooth, int nbucket)
{
  kdxptr kd;
  smxptr sm;
  bodyptr p;

  kd = init_kdtree(btab, nbody, nbody);		// prepare for kd tree
  build_kdtree(kd, nbucket);			// do tree construction
  sm = init_smooth(kd, nsmooth, slope);		// prepare for smoothing
  for (p = btab; p < btab+nbody; p++)		// loop over bodies
    Rho(p) = 0.0;				// prepare to sum density
  if (rhomass)					// mass density reqested?
    smooth(sm, sph_mass_density);		// compute mass density
  else
    smooth(sm, sph_numb_density);		// compute number density
  finish_smooth(sm);				// deallocate smooth data
  finish_kdtree(kd);				// deallocate kdtree data
}

//  ____________________________________________________________________
//  sph_numb_density: compute number density for body and its neighbors.

local void sph_numb_density(smxptr sm, int pi, int nball)
{
  bodyptr bi = sm->kd->bptr[pi], bj;
  real hinv2, wsc, rhinv2, wsm;
  int j;

  hinv2 = 4 / sm->r2ball[pi];
  wsc = 0.5 * rsqrt(hinv2) * hinv2 / PI;
  for (j = 0; j < nball; ++j) {
    bj = sm->kd->bptr[sm->inlist[j]];
    rhinv2 = sm->r2list[j] * hinv2;
    WSmooth(wsm, wsc, rhinv2, sm->coefs);
    Rho(bi) += wsm;
    Rho(bj) += wsm;
  }
}

//  __________________________________________________________________
//  sph_mass_density: compute mass density for body and its neighbors.

local void sph_mass_density(smxptr sm, int pi, int nball)
{
  bodyptr bi = sm->kd->bptr[pi], bj;
  real hinv2, wsc, rhinv2, wsm;
  int j;

  hinv2 = 4 / sm->r2ball[pi];
  wsc = 0.5 * rsqrt(hinv2) * hinv2 / PI;
  for (j = 0; j < nball; ++j) {
    bj = sm->kd->bptr[sm->inlist[j]];
    rhinv2 = sm->r2list[j] * hinv2;
    WSmooth(wsm, wsc, rhinv2, sm->coefs);
    Rho(bi) += wsm * Mass(bj);
    Rho(bj) += wsm * Mass(bi);
  }
}
