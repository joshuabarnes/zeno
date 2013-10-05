/*
 * gsprealize.c: make N-body realization of GSP.  This version uses
 * the range of radii listed in the density-profile gsp to define the
 * range of particle radii; this should be fixed!
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"

#if (!defined(LINUX) && !defined(MACOSX))
#include <ieeefp.h>
#endif

string defv[] = {		";Construct N-body realization of GSP.",
				";Uses Abel integral to compute DF.",
    "gsp=???",			";Input GSP for density profile",
    "out=",			";Output snapshot file with bodies",
    "grav=",		        ";Input GSP for gravitational potential.",
				";If blank, use density profile GSP.",
    "nstep=128",		";Integration steps for DF calculation",
    "dflist=false",		";Print out distribution function",
    "copies=1",			";Number of realizations to produce",
    "nbody=4096",		";Number of bodies per realization",
    "seed=54321",		";Seed for random number generator",
    "randrad=true",		";Pick radii randomly from M'(r).",
				";If false, sample radii uniformly.",
    "besort=true",		";Sort particles by binding energy",
    "zerocm=false",		";Transform to center of mass coords",
    "hmaxpar=1024,1.25",	";Parameters for hmax function",
    "VERSION=2.2",		";Josh Barnes  14 July 2012",
    NULL,
};

// Prototypes for model construction.

local void computedf(void);			// tabulate dist. func.
local void dfderiv(real *, real *);		// derivatives for f(E)
local void inithmax(void);			// init table for hmax(phi)
local void gsprealize(void);			// construct realization
local bool pickspeed(real *, real);		// pick speed from h(v)
local real h_v(real, real);			// speed distribution func
local real f_E(real);				// phase-space dist. func
local int berank(const void *, const void *);	// compare binding energies

real diffstep(real *, real *, int, void (*)(real *, real *), real);

// Global data for communication between routines.

local gsprof *gsp, *ggsp;		// profiles for mass and gravity
local real *Etab, *Ftab, *Fcoef;	// tables for f(E) = dF/dE
local real *htab, *hcoef, hsafe;	// tables for hmax(phi)
local int ntab;				// number of values in tables
local bodyptr btab = NULL;		// pointer to array of bodies
local int nbody;			// number of bodies in array

// Miscellaneous parameters and constants.

#define NTRIAL 1000			// rejection cycles before note

local string bodyfields[] = { PosTag, VelTag, MassTag, PhiTag, AuxTag, NULL };

//  ________________________________________________
//  main: handle I/O and call construction routines.

int main(int argc, string argv[])
{
  stream fstr, gstr, ostr = NULL;
  int ncopy;
  real tsnap = 0.0;

  initparam(argv, defv);
  fstr = stropen(getparam("gsp"), "r");
  get_history(fstr);
  gsp = get_gsprof(fstr);
  if (! strnull(getparam("grav"))) {
    gstr = stropen(getparam("grav"), "r");
    get_history(gstr);
    ggsp = get_gsprof(gstr);
  } else
    ggsp = gsp;
  computedf();
  inithmax();
  layout_body(bodyfields, Precision, NDIM);
  nbody = getiparam("nbody");
  if (nbody <= 0)
    error("%s: absurd value for nbody\n", getprog());
  init_random(getiparam("seed"));
  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
  }
  ncopy = getiparam("copies");
  while (--ncopy >= 0) {
    gsprealize();
    if (ostr != NULL)
      put_snap(ostr, &btab, &nbody, &tsnap, bodyfields);
  }
  if (ostr != NULL)
    strclose(ostr);
  return (0);
}

//  __________________________________________________________
//  computedf: numerically integrate F(E) and tabulate result.

local void computedf(void)
{
  real xmax, xFE[3], err, avgerr, maxerr;
  int k, nstep, i, j;

  (void) phi_gsp(ggsp, 1.0);			// compute potential table
  for (k = 0; ggsp->phi[k] == ggsp->phi[k+1]; k++);
  if (k > 0)
    eprintf("[%s.computedf: skipping first %d values of phi]\n",
	    getprog(), k);
  ntab = ggsp->npoint - k;
  Etab = ggsp->phi + k;
  Ftab = (real *) allocate(ntab * sizeof(real));
  Fcoef = (real *) allocate(3 * ntab * sizeof(real));
  nstep = getiparam("nstep");
  if (nstep <= 0)
    error("%s: absurd value for nstep\n", getprog());
  avgerr = maxerr = 0.0;
  for (i = 0; i < ntab; i++) {
    xmax = rsqrt(- (Etab[i] - Etab[ntab-1]));
    xFE[0] = 0.0;				// init independent var
    xFE[1] = 0.0;				// init integration result
    xFE[2] = Etab[i];				// pass E on to dfderiv()
    for (j = 0; j < nstep; j++) {
      err = diffstep(xFE, xFE, 3, dfderiv, xmax / nstep);
      avgerr += err / (ntab * nstep);
      maxerr = MAX(maxerr, err);
    }
    Ftab[i] = xFE[1];
  }
  eprintf("[%s.computedf: avgerr = %f  maxerr = %f]\n",
	  getprog(), avgerr, maxerr);
  spline(Fcoef, Etab, Ftab, ntab);
  if (getbparam("dflist")) {
    printf("#%11s  %12s  %12s\n", "E", "-F(E)", "f(E)");
    for (i = 0; i < ntab; i++)
      printf("%12.6f  %12.5e  %12.5e\n", Etab[i],
	     - seval(Etab[i], Etab, Ftab, Fcoef, ntab),
	     spldif(Etab[i], Etab, Ftab, Fcoef, ntab));
  }
}    

//  ________________________________________
//  dfderiv: evaluate dF/dx for computedf().

local void dfderiv(real *dxFE, real *xFE)
{
  real r, C = 1.0 / (M_SQRT2 * M_PI * M_PI);

  r = r_phi_gsp(ggsp, MIN(rsqr(xFE[0]) + xFE[2], 0.0));
  dxFE[0] = 1.0;
  dxFE[1] = C * (rsqr(r) / mass_gsp(ggsp, r)) * drho_gsp(gsp, r);
  dxFE[2] = 0.0;
  if (isnan((double) dxFE[1]))
    error("%s: distfunc undefined for x = %f, E = %f, r = %f\n",
	  getprog(), xFE[0], xFE[2], r);
}

//  _____________________________________________________
//  inithmax: initalize the table used to find hmax(phi).

local void inithmax(void)
{
  int nsamp, i, j;
  real vtop, h, v, hmax, vmax;

  if (sscanf(getparam("hmaxpar"), "%i,%f", &nsamp, &hsafe) != 2)
    error("%s: error scanning hmaxpar\n", getprog());
  htab = (real *) allocate(ntab * sizeof(real));
  hcoef = (real *) allocate(3 * ntab * sizeof(real));
  for (i = 0; i < ntab-1; i++) {
    vtop = rsqrt(-2.0 * (Etab[i] - Etab[ntab-1]));
    for (j = 0; j <= nsamp; j++) {
      v = vtop * j / ((real) nsamp);
      h = h_v(v, Etab[i]);
      if (j == 0 || h > hmax) {
	hmax = h;
	vmax = v;
      }
    }
    htab[i] = hmax;
  }
  htab[ntab-1] = 0.0;
  spline(hcoef, Etab, htab, ntab);
}

//  _____________________________________________________________
//  gsprealize: construct realization from distribution function.

local void gsprealize(void)
{
  bool randrad, speedfault;
  real mbody, mmin, mmax, x, rx, vx;
  int i;
  bodyptr bp;

  randrad = getbparam("randrad");
  mbody = gsp->mtot / nbody;
  mmin = gsp->mass[0];
  mmax = gsp->mass[gsp->npoint-1];
  if (mbody < mmin || mbody < (gsp->mtot - mmax))
    eprintf("[%s: WARNING: gsp limits range of radii\n"
	    " mbody = %g  mmin = %g  mtot-mmax = %g]\n",
	    getprog(), mbody, mmin, gsp->mtot - mmax);
  if (btab == NULL)
    btab = (bodyptr) allocate(nbody * SizeofBody);
  do {
    speedfault = FALSE;
    for (i = 0; i < nbody; i++) {
      bp = NthBody(btab, i);
      Mass(bp) = mbody;
      if (randrad)				// use random sampling
	x = xrandom(mmin, mmax);
      else					// use uniform sampling
	x = mmin + ((i + 0.5) / nbody) * (mmax - mmin);
      rx = r_mass_gsp(gsp, x);
      if (isnan((double) rx))
	error("%s: rx = NAN for x = %f\n", getprog(), x);
      pickshell(Pos(bp), NDIM, rx);
      Phi(bp) = phi_gsp(ggsp, rx);
      if (! pickspeed(&vx, Phi(bp))) {		// VN rejection failed?
	speedfault = TRUE;
	hsafe = 2 * hsafe;
	eprintf("[%s: warning: increasing hsafe to %g]\n", getprog(), hsafe);
	break;
      }
      if (isnan((double) vx))
	error("%s: vx = NAN for x = %f\n", getprog(), x);
      pickshell(Vel(bp), NDIM, vx);
      Aux(bp) = f_E(0.5 * rsqr(vx) + Phi(bp));
    }
  } while (speedfault);
  if (getbparam("besort"))
    qsort(btab, nbody, SizeofBody, berank);
  if (getbparam("zerocm"))
    snapcenter(btab, nbody, MassField.offset);
}

//  _____________________________________________________
//  pickspeed: chose speed distributed according to h(v).

local bool pickspeed(real *vp, real phi)
{
  real vmax, hmax, v0, h0, hv;
  int nloop = 0, nwarn = NTRIAL, i;

  if (phi - Etab[ntab-1] > 0)
    error("%s.pickspeed: imaginary vmax;  phi = %.8g  Etab[%d] = %.8g\n",
	  getprog(), phi, ntab-1, Etab[ntab-1]);
  vmax = rsqrt(-2.0 * (phi - Etab[ntab-1]));
  hmax = (Etab[0] < phi ? seval(phi, Etab, htab, hcoef, ntab) : htab[0]);
  if (hmax < 0)
    error("%s.pickspeed: hmax = %.8g < 0;  phi = %.8g  Etab[0] = %.8g\n",
	  getprog(), hmax, phi, Etab[0]);
  do {
    if (nloop >= nwarn) {
      eprintf("[%s.pickspeed: %d iterations]\n", getprog(), nloop);
      nwarn += NTRIAL;
    }
    nloop++;
    v0 = xrandom(0.0, vmax);
    h0 = xrandom(0.0, hsafe * hmax);
    hv = h_v(v0, phi);
    if (hv > hsafe * hmax) {			// fail if guess is OOB!
      for (i = 0; i < ntab-1; i++)
	if (Etab[i] <= phi && phi < Etab[i+1])
	  eprintf("[%s.pickspeed: warning: guess out of bounds\n"
		  "  hsafe*hmax = %g < hv = %g  v0 = %g  phi = %.8g\n"
		  "  Etab[%d:%d] = %.8g:%.8g  htab[%d:%d] = %g:%g]\n",
		  getprog(), hsafe * hmax, hv, v0, phi,
		  i, i+1, Etab[i], Etab[i+1], i, i+1, htab[i], htab[i+1]);
      if (phi < Etab[0])
	eprintf("[%s.pickspeed: warning: guess out of bounds\n"
		"  hsafe*hmax = %g < hv = %g  v0 = %g  phi = %.8g\n"
		"  Etab[0] = %.8g  htab[0] = %g]\n",
		getprog(), hsafe * hmax, hv, v0, phi, Etab[0], htab[0]);
      return (FALSE);
    }
  } while (hv < h0);
  *vp = v0;
  return (TRUE);
}

//  _________________________________________________________
//  h_v: compute speed distribution (up to a factor of 4 PI).

local real h_v(real v, real phi)
{
  real E;

  E = 0.5 * v * v + phi;
  if (E > 0)
    error("%s.h_v: E > 0 for v = %f, phi = %f\n", getprog(), v, phi);
  return (v * v * f_E(E));
}

//  ___________________________________
//  f_e: compute distribution function.

local real f_E(real E)
{
  real f;
  static bool warn = FALSE;

  f = spldif(E, Etab, Ftab, Fcoef, ntab);
  if (f < 0.0 && !warn) {
    eprintf("[%s.f_E: WARNING: f < 0 for E = %f]\n", getprog(), E);
    warn = TRUE;
  }
  return (f);
}

//  ______________________________________
//  berank: rank bodies by binding energy.

local int berank(const void *a, const void *b)
{
  real Ea, Eb;

  Ea = 0.5 * dotvp(Vel((bodyptr) a), Vel((bodyptr) a)) + Phi((bodyptr) a);
  Eb = 0.5 * dotvp(Vel((bodyptr) b), Vel((bodyptr) b)) + Phi((bodyptr) b);
  return (Ea < Eb ? -1 : Ea > Eb ? 1 : 0);
}
