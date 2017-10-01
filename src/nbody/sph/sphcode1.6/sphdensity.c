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
//  _____________________________________

string defv[] = {		";Compute densities and smoothing lengths",
  "in=???",			";Input, including " PosTag " data.",
				";" MassTag " data is optional.",
  "out=???",			";Output, with " RhoTag " and " SmoothTag
                                  " data",
  "rhomass=false",		";If TRUE, output particle mass density.",
				";Otherwise, output number density.",
  "outgrad=false",		";If TRUE, output density gradient.",
				";Uses " AuxVecTag " for gradient data.",
  "nsmooth=40",			";Number of bodies in smoothing volume",
  "nbucket=16",			";Number of bodies in leaves of KD tree",
  "slope=0.0",			";Smoothing kernel slope at origin",
  "VERSION=1.1",		";Joshua Barnes  12 August 2017",
  NULL,
};

//  Local procedure prototypes.
//  ___________________________

local void setupbody(bool);			// set offsets for dynbody
local void sphdensity(bool, bool, real, int, int);	// perform smoothing
local void sph_mass_density(smxptr, int, int);	// compute mass density
local void sph_mass_gradient(smxptr, int, int);	// compute gradient
local void sph_numb_density(smxptr, int, int);	// compute number density
local void sph_numb_gradient(smxptr, int, int);	// compute gradient

//  main: toplevel routine.
//  _______________________

int main(int argc, string argv[])
{
  bool outgrad, massdata = FALSE;
  stream istr, ostr = NULL;
  string intags[MaxBodyFields], *newtags, *outtags;
  bodyptr p;

  initparam(argv, defv);			// initialize param access
  outgrad = getbparam("outgrad");
  setupbody(outgrad);				// interface to dynbody
  newtags = set_cons(SmoothTag, RhoTag, (outgrad ? AuxVecTag : NULL), NULL);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  while (get_snap(istr, &btab, &nbody, &tnow, intags, FALSE, NULL)) {
    if (! set_member(intags, PosTag))
      error("%s: %s input data missing\n", getprog(), PosTag);
    massdata = massdata || set_member(intags, MassTag);
    if (getbparam("rhomass") && ! massdata)
      error("%s: %s input data required for mass density\n",
	    getprog(), MassTag);
    for (p = btab; p < btab+nbody; p++)
      Type(p) = BODY | GAS;
    sphdensity(getbparam("rhomass"), outgrad, getdparam("slope"),
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

//  setupbody: inform dynamic body routines of relevant body fields.
//  ________________________________________________________________

local void setupbody(bool outgrad)
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
  if (outgrad)					// use vmid for gradient
    define_body_offset(AuxVecTag, BodyOffset(Vmid));
}

//  sphdensity: compute smoothing length and density.
//  _________________________________________________

local void sphdensity(bool usemass, bool outgrad,
		      real slope, int nsmooth, int nbucket)
{
  kdxptr kd;
  smxptr sm;
  bodyptr p;

  kd = init_kdtree(btab, nbody, nbody);		// prepare for kd tree
  build_kdtree(kd, nbucket);			// do tree construction
  sm = init_smooth(kd, nsmooth, slope);		// prepare for smoothing
  for (p = btab; p < btab+nbody; p++)		// loop over bodies
    Rho(p) = 0.0;				// prepare to sum density
  if (usemass) {				// mass density reqested?
    smooth(sm, sph_mass_density);		// compute mass density
    if (outgrad)
      smooth(sm, sph_mass_gradient);		// compute density gradient
  } else {
    smooth(sm, sph_numb_density);		// compute number density
    if (outgrad)
      smooth(sm, sph_numb_gradient);		// compute density gradient
  }
  finish_smooth(sm);				// deallocate smooth data
  finish_kdtree(kd);				// deallocate kdtree data
}

//  sph_mass_density: compute density for body and its neighbors.
//  _____________________________________________________________

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

//  sph_mass_gradient: compute density gradient for body.
//  _____________________________________________________

local void sph_mass_gradient(smxptr sm, int pi, int nball)
{
  bodyptr bi = sm->kd->bptr[pi], bj;
  real hinv2, dwsc, rhinv2, dwsm, dwsmrinv, fij;
  vector gradrho, rij;
  int j;

  hinv2 = 4 / sm->r2ball[pi];
  dwsc = hinv2 * hinv2 / PI;
  CLRV(gradrho);
  for (j = 0; j < nball; ++j) {
    bj = sm->kd->bptr[sm->inlist[j]];
    if (bi != bj) {
      rhinv2 = sm->r2list[j] * hinv2;
      dWSmooth(dwsm, dwsc, rhinv2, sm->coefs);
      dwsmrinv = dwsm / sqrt(sm->r2list[j]);
      SUBV(rij, Pos(bi), Pos(bj));
      fij = Rho(bj) / Rho(bi) - 1.0;
      ADDMULVS(gradrho, rij, fij * dwsmrinv * Mass(bj));
    }
  }
  SETV(Vmid(bi), gradrho);			// store gradient in vmid
}

//  sph_numb_density: compute density for body and its neighbors.
//  _____________________________________________________________

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

//  sph_numb_gradient: compute density gradient for body.
//  _____________________________________________________

local void sph_numb_gradient(smxptr sm, int pi, int nball)
{
  bodyptr bi = sm->kd->bptr[pi], bj;
  real hinv2, dwsc, rhinv2, dwsm, dwsmrinv, fij;
  vector gradrho, rij;
  int j;

  hinv2 = 4 / sm->r2ball[pi];
  dwsc = hinv2 * hinv2 / PI;
  CLRV(gradrho);
  for (j = 0; j < nball; ++j) {
    bj = sm->kd->bptr[sm->inlist[j]];
    if (bi != bj) {
      rhinv2 = sm->r2list[j] * hinv2;
      dWSmooth(dwsm, dwsc, rhinv2, sm->coefs);
      dwsmrinv = dwsm / sqrt(sm->r2list[j]);
      SUBV(rij, Pos(bi), Pos(bj));
      fij = Rho(bj) / Rho(bi) - 1.0;
      ADDMULVS(gradrho, rij, fij * dwsmrinv);
    }
  }
  SETV(Vmid(bi), gradrho);			// store gradient in vmid
}
