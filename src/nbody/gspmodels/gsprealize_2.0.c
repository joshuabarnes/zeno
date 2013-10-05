/*
 * GSPREALIZE.C: make N-body realization of GSP.  This version still
 * uses the range of radii listed in the density-profile gsp to set
 * the range of particle radii; this should be fixed!
 */

#include "stdinc.h"
#include "assert.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"

#if (!defined(LINUX) && !defined(MACOSX))
#include <ieeefp.h>
#endif

#include <gsl/gsl_rng.h>

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
    "hmaxpar=256,1.125",	";Parameters for hmax function",
    "VERSION=2.0",		";Josh Barnes  12 May 2012",
    NULL,
};

/* Prototypes for model construction.					    */

local void computedf(void);			/* tabulate dist. func.     */
local void dfderiv(real *, real *);		/* derivatives for f(E)     */
local void inithmax(void);			/* init table for hmax(phi) */
local void gsprealize(void);			/* construct realization    */
local real pickspeed(real);			/* pick speed from h(v)     */
local real h_v(real, real);			/* speed distribution func  */
local real f_E(real);				/* phase-space dist. func   */
local int berank(const void *, const void *);	/* compare binding energies */

real diffstep(real *, real *, int, void (*)(real *, real *), real);

/* Global data for communication between routines.			    */

local gsprof *gsp, *ggsp;		/* profiles for mass and gravity    */
local real *Etab, *Ftab, *Fcoef;	/* tables for f(E) = dF/dE          */
local int ntab;				/* number of values in tables	    */
local real *htab, *hcoef;		/* tables for hmax(phi)		    */
local bodyptr btab = NULL;		/* pointer to array of bodies	    */
local int nbody;			/* number of bodies in array	    */

local gsl_rng *rng;

/* Miscellaneous parameters and constants.				    */

#define NTRIAL 1000			/* rejection cycles before warning  */

local string bodyfields[] = { PosTag, VelTag, MassTag, PhiTag, AuxTag, NULL };

/*
 * MAIN: handle I/O and call construction routines.
 */

int main(int argc, string argv[])
{
  stream fstr, gstr, ostr;
  int n;
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
  if (! strnull(getparam("out"))) {
    layout_body(bodyfields, Precision, NDIM);
    nbody = getiparam("nbody");
    assert(nbody > 0);
    rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng, getiparam("seed"));
    init_random(getiparam("seed"));		/* set for pickshell calls */
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    n = getiparam("copies");
    while (--n >= 0) {
      gsprealize();
      put_snap(ostr, &btab, &nbody, &tsnap, bodyfields);
    }
    strclose(ostr);
  }
  return (0);
}

/*
 * COMPUTEDF: numerically integrate F(E) and tabulate result.
 */

local void computedf(void)
{
  real xmax, xFE[3], err, avgerr, maxerr;
  int k, nstep, i, j;

  (void) phi_gsp(ggsp, 1.0);			/* compute potential table  */
  for (k = 0; ggsp->phi[k] == ggsp->phi[k+1]; k++);
  if (k > 0)
    eprintf("[%s.computedf: skipping first %d values of phi]\n",
	    getargv0(), k);
  ntab = ggsp->npoint - k;
  Etab = ggsp->phi + k;
  Ftab = (real *) allocate(ntab * sizeof(real));
  Fcoef = (real *) allocate(3 * ntab * sizeof(real));
  nstep = getiparam("nstep");
  assert(nstep > 0);
  avgerr = maxerr = 0.0;
  for (i = 0; i < ntab; i++) {
    xmax = rsqrt(- (Etab[i] - Etab[ntab-1]));
    xFE[0] = 0.0;				/* init independent var     */
    xFE[1] = 0.0;				/* init integration result  */
    xFE[2] = Etab[i];				/* pass E on to dfderiv()   */
    for (j = 0; j < nstep; j++) {
      err = diffstep(xFE, xFE, 3, dfderiv, xmax / nstep);
      avgerr += err / (ntab * nstep);
      maxerr = MAX(maxerr, err);
    }
    Ftab[i] = xFE[1];
  }
  eprintf("[%s.computedf: avgerr = %f  maxerr = %f]\n",
	  getargv0(), avgerr, maxerr);
  spline(Fcoef, Etab, Ftab, ntab);
  if (getbparam("dflist")) {
    printf("#%11s  %12s  %12s\n", "E", "-F(E)", "f(E)");
    for (i = 0; i < ntab; i++)
      printf("%12.6f  %12.5e  %12.5e\n", Etab[i],
	     - seval(Etab[i], Etab, Ftab, Fcoef, ntab),
	     spldif(Etab[i], Etab, Ftab, Fcoef, ntab));
  }
}    

/*
 * DFDERIV: evaluate dF/dx for computedf().
 */

local void dfderiv(real *dxFE, real *xFE)
{
  real r, C = 1.0 / (M_SQRT2 * M_PI * M_PI);

  r = r_phi_gsp(ggsp, MIN(rsqr(xFE[0]) + xFE[2], 0.0));
  dxFE[0] = 1.0;
  dxFE[1] = C * (rsqr(r) / mass_gsp(ggsp, r)) * drho_gsp(gsp, r);
  dxFE[2] = 0.0;
  if (isnan((double) dxFE[1]))
    error("%s: distfunc undefined for x = %f, E = %f, r = %f\n",
	  getargv0(), xFE[0], xFE[2], r);
}

/*
 * INITHMAX: initalize the table used to find hmax(phi).
 */

local void inithmax(void)
{
  int nsamp, i, j;
  real hsafe, vtop, h, v, hmax, vmax;

  if (sscanf(getparam("hmaxpar"), "%i,%f", &nsamp, &hsafe) != 2)
    error("%s: error scanning hmaxpar\n", getargv0());
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
    htab[i] = hsafe * hmax;
  }
  htab[ntab-1] = 0.0;
  spline(hcoef, Etab, htab, ntab);
}

/*
 * GSPREALIZE: construct realization from distribution function.
 */

local void gsprealize(void)
{
  bool randrad;
  real mmin, mmax, x, rx, vx;
  int i;
  bodyptr bp;

  randrad = getbparam("randrad");
  mmin = gsp->mass[0];
  mmax = gsp->mass[gsp->npoint-1];
  if (btab == NULL)
    btab = (bodyptr) allocate(nbody * SizeofBody);
  for (i = 0; i < nbody; i++) {
    bp = NthBody(btab, i);
    Mass(bp) = gsp->mtot / nbody;
    if (randrad)				/* use random sampling	    */
      x = mmin + (mmax - mmin) * gsl_rng_uniform(rng);
    else					/* use uniform sampling	    */
      x = mmin + ((i + 0.5) / nbody) * (mmax - mmin);
    rx = r_mass_gsp(gsp, x);
    if (isnan((double) rx))
      error("%s: rx = NAN for x = %f\n", getargv0(), x);
    pickshell(Pos(bp), NDIM, rx);
    Phi(bp) = phi_gsp(ggsp, rx);
    vx = pickspeed(Phi(bp));
    if (isnan((double) vx))
      error("%s: vx = NAN for x = %f\n", getargv0(), x);
    pickshell(Vel(bp), NDIM, vx);
    Aux(bp) = f_E(0.5 * rsqr(vx) + Phi(bp));
  }
  if (getbparam("besort"))
    qsort(btab, nbody, SizeofBody, berank);
  if (getbparam("zerocm"))
    snapcenter(btab, nbody, MassField.offset);
}

/*
 * PICKSPEED: chose speed distributed according to h(v).
 */

local real pickspeed(real phi)
{
  real vmax, hmax, v0, h0, hv;
  int n, nwarn;

  vmax = rsqrt(-2.0 * (phi - Etab[ntab-1]));
  hmax = (Etab[0] < phi ? seval(phi, Etab, htab, hcoef, ntab) : htab[0]);
  if (hmax <= 0)
    error("%s.pickspeed: hmax = %g < 0;  phi = %g  Etab[0] = %g\n",
	  getargv0(), hmax, phi, Etab[0]);
  nwarn = NTRIAL;
  n = 0;
  do {
    if (n > nwarn) {
      eprintf("[%s.pickspeed: %d iterations of pickspeed]\n",
	      getargv0(), n);
      nwarn += NTRIAL;
    }
    v0 = vmax * gsl_rng_uniform(rng);
    h0 = hmax * gsl_rng_uniform(rng);
    hv = h_v(v0, phi);
    if (hv > hmax)
      error("%s.pickspeed: guess out of bounds\n"
	    "  hmax = %g < hv = %g for v0 = %g, phi = %g\n",
	    getargv0(), hmax, hv, v0, phi);
    n++;
  } while (hv < h0);
  return (v0);
}

/*
 * H_V: compute speed distribution (up to a factor of 4 PI).
 */

local real h_v(real v, real phi)
{
  real E;

  E = 0.5 * v * v + phi;
  if (E > 0)
    error("%s.h_v: E > 0 for v = %f, phi = %f\n", getargv0(), v, phi);
  return (v * v * f_E(E));
}

/*
 * F_E: compute distribution function.
 */

local real f_E(real E)
{
  real f;
  static bool warn = FALSE;

  f = spldif(E, Etab, Ftab, Fcoef, ntab);
  if (f < 0.0 && !warn) {
    eprintf("[%s.f_E: WARNING: f < 0 for E = %f]\n", getargv0(), E);
    warn = TRUE;
  }
  return (f);
}

/*
 * BERANK: rank bodies by binding energy.
 */

local int berank(const void *a, const void *b)
{
  real Ea, Eb;

  Ea = 0.5 * dotvp(Vel((bodyptr) a), Vel((bodyptr) a)) + Phi((bodyptr) a);
  Eb = 0.5 * dotvp(Vel((bodyptr) b), Vel((bodyptr) b)) + Phi((bodyptr) b);
  return (Ea < Eb ? -1 : Ea > Eb ? 1 : 0);
}
