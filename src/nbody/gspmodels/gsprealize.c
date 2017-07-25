/*
 * gsprealize.c: make N-body realization of GSP with isotropic
 * or Osipkov-Merritt distribution function
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"
#include <gsl/gsl_integration.h>
#include <assert.h>

string defv[] = {		";Construct N-body realization of GSP.",
				";Uses Abel integral to compute DF.",
  "gsp=???",			";Input GSP for density profile",
  "out=",			";Output snapshot file with bodies",
  "grav=",		        ";Input GSP for gravitational potential.",
				";If none, use density profile GSP.",
  "raniso=-1",			";Osipkov-Merritt anisotropy radius.",
				";If raniso <= 0, find isotropic DF.",
  "epsrel=1.0e-4",		";Integration relative tolerance.",
				";Need to relax for anisotropic models.",
  "copies=1",			";Number of realizations to produce",
  "nbody=16384",		";Number of bodies per realization",
  "seed=54321",			";Seed for random number generator",
  "randrad=true",		";Pick radii randomly from M(r).",
				";If FALSE, sample M(r) uniformly.",
  "extmass=true",		";Extend model beyond limit of gsp",
  "besort=true",		";Sort particles by binding energy",
  "zerocm=false",		";Transform to center of mass coords",
  "hmaxpar=1024,1.25",		";Parameters for hmax function",
  "VERSION=2.5",		";Josh Barnes  14 July 2017",
  NULL,
};

// Prototypes for model construction.

void gsp_realize(bool extmass, bool randrad);	// construct realization
void init_hmax(void);				// init table for hmax(phi)
bool pick_speed(double *, double);		// pick speed from h(v)
double hfunc(double, double);			// speed distribution func
int berank(const void *, const void *);		// compare binding energies

// Global data for communication between routines.

gsprof *gsp, *ggsp;			// profiles for mass and gravity
double *ptab, *htab, hsafe;		// table for hmax(phi), safe factor
gsl_interp *hp_spline;			// spline fit to hmax(phi)

string bodyfields[] = { PosTag, VelTag, MassTag, PhiTag, AuxTag, NULL };
bodyptr btab = NULL;			// pointer to array of bodies
int nbody;				// number of bodies in array

#define NTRIAL 1000

//  main: handle I/O and call construction routines.
//  ________________________________________________

int main(int argc, string argv[])
{
  stream fstr, gstr, ostr = NULL;
  int ncopy;
  real tsnap = 0.0;
  double raniso, epsrel;

  initparam(argv, defv);
  fstr = stropen(getparam("gsp"), "r");
  get_history(fstr);
  gsp = gsp_read(fstr);
  if (! strnull(getparam("grav"))) {
    gstr = stropen(getparam("grav"), "r");
    get_history(gstr);
    ggsp = gsp_read(gstr);
  } else
    ggsp = gsp;
  raniso = getdparam("raniso");
  epsrel = getdparam("epsrel");
  gsp_calc_dist_pars(NULL, &epsrel, NULL);
  gsp_calc_dist(gsp, ggsp, raniso);		// compute distribution func
  init_hmax();
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
    gsp_realize(getbparam("extmass"), getbparam("randrad"));
    if (getbparam("besort"))
      qsort(btab, nbody, SizeofBody, berank);
    if (getbparam("zerocm"))
      snapcenter(btab, nbody, MassField.offset);
    if (ostr != NULL) {
      put_snap(ostr, &btab, &nbody, &tsnap, bodyfields);
      fflush(ostr);
    }
  }
  fflush(NULL);
  if (ggsp != gsp)
    gsp_free(ggsp);
  gsp_free(gsp);
  return 0;
}

//  gsp_realize: construct realization from distribution function.
//  ______________________________________________________________

void gsp_realize(bool extmass, bool randrad)
{
  double mbody = gsp_mtot(gsp) / nbody, mmin, mmax, m, r, v, p, q;
  bodyptr bp;
  int nint = 0, next = 0;
  bool speedfault;
  vector vrad, vtan;

  mmin = (extmass ? 0.0           : gsp->mass[0]);
  mmax = (extmass ? gsp_mtot(gsp) : gsp->mass[gsp->npoint-1]);
  if (mbody < mmin || mbody < (gsp_mtot(gsp) - mmax))
    eprintf("[%s.gsp_realize: warning: gsp limits range of radii\n"
	    " mbody = %g  mmin = %g  mtot-mmax = %g]\n",
	    getprog(), mbody, mmin, gsp_mtot(gsp) - mmax);
  if (btab == NULL)
    btab = (bodyptr) allocate(nbody * SizeofBody);
  for (int i = 0; i < nbody; i++) {
    bp = NthBody(btab, i);
    Mass(bp) = mbody;
    m = (randrad ? xrandom(mmin, mmax) :
	           mmin + ((i + 0.5) / nbody) * (mmax - mmin));
    if (m < gsp->mass[0])
      nint++;
    if (m > gsp->mass[gsp->npoint-1])
      next++;
    r = gsp_mass_rad(gsp, m);
    pickshell(Pos(bp), NDIM, r);
    Phi(bp) = gsp_phi(ggsp, r);
  }
  if (nint > 0 || next > 0)
    eprintf("[%s.gsp_realize: nint, next = %d, %d]\n", getprog(), nint, next);
  speedfault = TRUE;				// always take 1st loop
  while (speedfault) {
    speedfault = FALSE;
    for (int i = 0; i < nbody; i++) {
      bp = NthBody(btab, i);
      if (! pick_speed(&v, Phi(bp))) {		// VN rejection failed?
	hsafe = 2 * hsafe;
	speedfault = TRUE;			// try again w/ bigger margin
	eprintf("[%s.gsp_realize: warning: increasing hsafe to %g]\n",
		getprog(), hsafe);
	break;					// break out of for loop
      }
      Aux(bp) = gsp_dist(gsp, 0.5 * v*v + Phi(bp));
      pickshell(Vel(bp), NDIM, v);
      if (gsp->raniso > 0.0) {			// making aniso. model?
	r = absv(Pos(bp));
	p = 1 / sqrt(1 + gsl_pow_2(r / gsp->raniso));
	q = dotvp(Vel(bp), Pos(bp)) / (r * r);
	MULVS(vrad, Pos(bp), q);		// get radial vel. vector
	SUBV(vtan, Vel(bp), vrad);		// and tangent. vel. vector
	MULVS(vtan, vtan, p);			// squeeze tangent. component
	ADDV(Vel(bp), vrad, vtan);		// put components together
      }
    }
  }
}

//  init_hmax: initalize the table used to find hmax(phi).
//  ______________________________________________________

void init_hmax(void)
{
  int nsamp, ntab;
  double vesc, h, v, hmax, vmax;

  if (sscanf(getparam("hmaxpar"), "%i,%lf", &nsamp, &hsafe) != 2)
    error("%s: error scanning hmaxpar\n", getprog());
  ntab = gsp->npoint;				// ABSTRACTION VIOLATION
  ptab = (double *) allocate(ntab * sizeof(double));
  htab = (double *) allocate(ntab * sizeof(double));
  for (int i = 0; i < ntab; i++) {
    ptab[i] = gsp->energy[i];			// ABSTRACTION VIOLATION
    vesc = sqrt(-2.0 * ptab[i]);
    vmax = 2 * sqrt(ptab[i] / (1 + 2 * gsp_beta(gsp)));
    hmax = hfunc(vmax, ptab[i]);		// set asymptotic values
    for (int j = 0; j < nsamp; j++) {
      v = vesc * j / ((double) nsamp);		// don't include v = vesc
      h = hfunc(v, ptab[i]);
      if (h > hmax) {
	hmax = h;
	vmax = v;
      }
    }
    htab[i] = hmax;
    // printf("%16.8e  %16.8e  %16.8e  %16.8e  %16.8e\n", ptab[i], hmax, vmax,
    //	      hfunc(2 * sqrt(ptab[i] / (1 + 2 * gsp_beta(gsp))), ptab[i]),
    //        2 * sqrt(ptab[i] / (1 + 2 * gsp_beta(gsp))));
  }
  hp_spline = gsl_interp_alloc(gsl_interp_akima, ntab);
  gsl_interp_init(hp_spline, ptab, htab, ntab);  
}

//  pick_speed: chose speed distributed according to h(v).
//  ______________________________________________________

bool pick_speed(double *vp, double phi)
{
  double vesc, phi0 = ptab[0], phiN = ptab[gsp->npoint-1], hmax, v0, h0, hv;
  int nloop = 0, nwarn = NTRIAL;

  if (phi > 0)
    error("%s.pick_speed: phi = %.8g > 0\n", getprog(), phi);
  vesc = sqrt(-2.0 * phi);
  hmax = (phi < phi0 ? htab[0] :
	  phi > phiN ? hfunc(2 * sqrt(phi / (1 + 2*gsp_beta(gsp))), phi) :
	               gsl_interp_eval(hp_spline, ptab, htab, phi, NULL));
  if (hmax <= 0)
    error("%s.pick_speed: hmax = %.8g <= 0  phi = %.8g  ptab[0] = %.8g\n",
	  getprog(), hmax, phi, phi0);
  do {
    if (nloop >= nwarn) {
      eprintf("[%s.pick_speed: warning: %d iterations]\n", getprog(), nloop);
      nwarn += NTRIAL;
    }
    nloop++;
    v0 = xrandom(0.0, vesc);
    h0 = xrandom(0.0, hsafe * hmax);
    hv = hfunc(v0, phi);
    if (hv > hsafe * hmax) {			// if guess is out of bounds
      eprintf("[%s.pick_speed: warning: guess out of bounds\n"
	      " hsafe*hmax = %g < hv = %g  v0 = %g  phi = %.8g]\n",
	      getprog(), hsafe * hmax, hv, v0, phi);
      return (FALSE);
    }
  } while (hv < h0);
  *vp = v0;
  return (TRUE);
}

//  hfunc: compute speed distribution (up to a factor of 4 PI).
//  ___________________________________________________________

double hfunc(double v, double phi)
{
  double E = 0.5 * v * v + phi;

  if (E > 0) {
    eprintf("[%s.hfunc: warning: E = %e > 0 for v = %e, phi = %e]\n",
	    getprog(), E, v, phi);
    return 0;
  }
  return (v * v * gsp_dist(gsp, E));
}

//  berank: rank bodies by binding energy.
//  ______________________________________

int berank(const void *a, const void *b)
{
  real Ea, Eb;

  Ea = 0.5 * dotvp(Vel((bodyptr) a), Vel((bodyptr) a)) + Phi((bodyptr) a);
  Eb = 0.5 * dotvp(Vel((bodyptr) b), Vel((bodyptr) b)) + Phi((bodyptr) b);
  return (Ea < Eb ? -1 : Ea > Eb ? 1 : 0);
}
