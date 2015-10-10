/*
 * polysnap.c: Construct generalized polytrope.
 */
 
#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "snapcenter.h"
#include <assert.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
 
string defv[] = {		";Construct N-body generalized polytrope",
  "out=",			";Snapshot output file",
  "n=1.5",			";Binding energy index; must be > 1/2",
  "m=1.5",			";Angular momentum index; must be > -1.",
				";Exceptions: n,m = (1,-1) or (0.5,-0.5).",
  "nbody=16384",		";Number of bodies per realization",
  "nmodel=1",			";Number of realizations to generate",
  "seed=54321",			";Random number seed",
  "hstep=0.01",			";Initial integration step",
  "besort=true",		";Sort particles by binding energy",
  "zerocm=false",		";Transform to center of mass coords",
  "listmodel=false",		";Tabulate model solution on stdout",
  "VERSION=1.5",		";Josh Barnes	18 June 2015",
  NULL,
};

//  Procedure prototypes.
//  _____________________

void polysolve(double, bool);		// model solution (double precision)
void asmpstep(double [], double);
int diffrpmw(double t, const double rpmw[], double drpmw[], void *params);
double rho(double, double);
void storestep(double []);
void fixsurface(void);
void polyscale(void);

void polymodel(void);			// model realization (real precision)
real rad_m(real);
real pick_psi(void);
real hfunc(real);
real pick_v(real);
real gfunc1(real);
real gfunc2(real);
void polymodel1(void);
void polymodel2(void);

real vnpick(real (*)(real), real, real, real, string);
int berank(const void *, const void *);

extern double lgamma(double);

//  Global data and parameters.
//  ___________________________

string bodyfields[] = { PosTag, VelTag, MassTag, PhiTag, AuxTag, NULL };

bodyptr btab;			// body array generated below
int nbody;			// number of bodies in array

#define K	1.0		// constant in rho(r, phi)
#define PHI0	-1.0		// phi(r=0) before rescaling
#define MAXSTEP 1024		// maximum integration steps

double npol, mpol;		// polytropic indicies
double Kprime;			// constant in f(E, J)
double rtab[MAXSTEP];		// radial integration points
double ptab[MAXSTEP];		// potential as a fn of radius
double mtab[MAXSTEP];		// enclosed mass as a fn of radius
int nstep;			// number of integration steps
double rad1, phi1, mtot;	// surface radius, potential, mass
gsl_spline *pmspline;		// spline for phi as function of mass

int main(int argc, string argv[])
{
  stream outstr = NULL;
  int nmodel;
  real tzero = 0.0;

  initparam(argv, defv);
  layout_body(bodyfields, Precision, NDIM);
  npol = getdparam("n");
  mpol = getdparam("m");
  nbody = getiparam("nbody");
  btab = (bodyptr) allocate(nbody * SizeofBody);
  init_random(getiparam("seed"));
  nmodel = getiparam("nmodel");
  if (! ((npol == 1.0 && mpol == -1.0) ||
	 (npol == 0.5 && mpol == -0.5)))
    polysolve(getdparam("hstep"), getbparam("listmodel"));
  if (! strnull(getparam("out"))) {
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
  }
  while (--nmodel >= 0) {
    if (npol == 1.0 && mpol == -1.0)
      polymodel1();
    else if (npol == 0.5 && mpol == -0.5)
      polymodel2();
    else
      polymodel();
    if (getbparam("besort"))
      qsort(btab, nbody, SizeofBody, berank);
    if (getbparam("zerocm"))
      snapcenter(btab, nbody, MassField.offset);
    if (outstr != NULL)
      put_snap(outstr, &btab, &nbody, &tzero, bodyfields);
    fflush(outstr);
  }
  return (0);
}

//  polysolve: solve structure equation for generalized polytrope.
//  ______________________________________________________________
 
void polysolve(double hstep, bool listmodel)
{
  double h, rpmw[4], w1, w2, t;
  gsl_odeiv2_system sys = { diffrpmw, NULL, 4, NULL };
  gsl_odeiv2_driver *drv;
  int stat;
 
  if (npol <= 0.5 || mpol <= -1.0)
    error("%s: illegal value for n or m\n", getprog());
  if (npol >= 3*mpol + 5)
    error("%s: model would have infinite radius\n", getprog());
  Kprime = K * (mpol+1) * pow(PI, -1.5) * pow(2.0, - (mpol+1.5)) *
           exp(lgamma(mpol+npol+1) - lgamma(mpol+2) - lgamma(npol-0.5));
  h = hstep;
  do {
    drv = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4,
					h, 1.0, 1.0);
    rpmw[0] = 0.0;				// initialize radius
    rpmw[1] = PHI0;				// initialize potential
    rpmw[2] = 0.0;				// initialize enclosed mass
    rpmw[3] = 0.0;				// initialize binding energy
    nstep = 0;
    storestep(rpmw);
    asmpstep(rpmw, h);				// use asymp. approximation
    storestep(rpmw);
    while (rpmw[1] < 0.0) {			// while potential is neg.
      stat = gsl_odeiv2_driver_apply_fixed_step(drv, &t, h, 1, rpmw);
      if (stat)
	eprintf("[%s.polysolve: stat = %d]\n", getprog(), stat);
      storestep(rpmw);
    }
    fixsurface();
    w1 = rpmw[3] - 0.5 * rsqr(mtot) / rad1;
    w2 = - (2*mpol + 3) / (3*mpol - npol + 5) * rsqr(mtot) / rad1;
    eprintf("[%s.polysolve: nstep = %3d  W = %f  error = %f]\n",
	    getprog(), nstep, w2, (w1 - w2)/w2);
    gsl_odeiv2_driver_free(drv);
    h = 0.5 * h;				// refine step-size by 2
  } while (nstep < MAXSTEP/4);
  polyscale();
  pmspline = gsl_spline_alloc(gsl_interp_cspline, nstep);
  gsl_spline_init(pmspline, rtab, ptab, nstep);
  if (listmodel) {
    printf("#%11s %11s %11s\n", "radius", "mass", "potential");
    for (int i = 0; i < nstep - 1; i++)
      printf(" %11.6f %11.6f %11.6f\n", rtab[i], mtab[i], ptab[i]);
    printf(" %11.6f %11.6f %11.6f\n", rad1, mtot, phi1);
  }
}

//  asmpstep: asymptotic approx. for generalized polytrope near r=0.
//  ________________________________________________________________

void asmpstep(double rpmw[], double h)
{
  rpmw[0] = h;
  rpmw[1] = (FOUR_PI * K / ((2*mpol + 2) * (2*mpol + 3))) *
              pow(- PHI0, npol + mpol) * pow(h, 2*mpol + 2) + PHI0;
  rpmw[2] = (FOUR_PI * K / (2*mpol + 3)) *
              pow(- PHI0, npol + mpol) * pow(h, 2*mpol + 3);
  rpmw[3] = 0.5 * PHI0 * rpmw[2];
}

//  diffrpmw: differentials of radius, potential, mass, binding energy.
//  ___________________________________________________________________

int diffrpmw(double t, const double rpmw[], double drpmw[], void *params)
{
  drpmw[0] = 1.0;
  drpmw[1] = rpmw[2] / (rpmw[0]*rpmw[0]);
  drpmw[2] = FOUR_PI * (rpmw[0]*rpmw[0]) * rho(rpmw[0], rpmw[1]);
  drpmw[3] = 0.5 * rpmw[1] * drpmw[2];
  return (GSL_SUCCESS);
}

double rho(double r, double phi)
{
  return (phi < 0.0 ? K * pow(- phi, npol + mpol) * pow(r, 2 * mpol) : 0.0);
}

//  storestep: tabulate results of integration.
//  ___________________________________________

void storestep(double rpmw[])
{
  if (nstep >= MAXSTEP)
    error("%s.storestep: too many steps\n", getprog());
  rtab[nstep] = rpmw[0];
  ptab[nstep] = rpmw[1];
  mtab[nstep] = rpmw[2];
  nstep++;
}

//  fixsurface: locate radius of surface by linear interpolation.
//  _____________________________________________________________

void fixsurface(void)
{
  double f;

  while (mtab[nstep-2] == mtab[nstep-1]) {
    eprintf("[%s.fixsurface: reducing nstep to %d]\n", getprog(), nstep-1);
    ptab[nstep-2] = ptab[nstep-1];
    rtab[nstep-2] = rtab[nstep-1];
    nstep--;
  }
  f = -ptab[nstep-2] / (ptab[nstep-1] - ptab[nstep-2]);
  rad1 = f * rtab[nstep-1] + (1 - f) * rtab[nstep-2];
  mtot = f * mtab[nstep-1] + (1 - f) * mtab[nstep-2];
  phi1 = 0.0;					// fixed by definition
}

//  polyscale: scale model to Henon's units.
//  ________________________________________

void polyscale(void)
{
  double r_henon, m_henon, p_henon, scl_r, scl_m, scl_v;
  int i;

  r_henon = (4*mpol + 6) / (3*mpol - npol + 5);
  m_henon = 1.0;
  p_henon = m_henon / r_henon;
  scl_r = r_henon / rad1;
  scl_m = m_henon / mtot;
  scl_v = rsqrt(scl_m / scl_r);
  for (i = 0; i < nstep; i++) {
    rtab[i] = scl_r * rtab[i];
    mtab[i] = scl_m * mtab[i];
    ptab[i] = rsqr(scl_v) * ptab[i] - p_henon;
    if (i > 0 && mtab[i-1] >= mtab[i])
      error("%s: mass not monotonic!  mtab[%d:%d] = %f,%f\n", getprog(),
	    i-1, i, mtab[i-1], mtab[i]);
  }
  rad1 = scl_r * rad1;
  mtot = scl_m * mtot;
  phi1 = rsqr(scl_v) * phi1 - p_henon;		// == - p_henon by def
  Kprime = Kprime * scl_m / rqbe(scl_r * scl_v);
  eprintf("[%s.polyscale: Kprime = %f]\n", getprog(), Kprime);
}

//  polymodel: construct N-body realization of polytrope.
//  _____________________________________________________

void polymodel(void)
{
  gsl_interp_accel *pmsplacc = gsl_interp_accel_alloc();
  bodyptr p;
  real rad, phi, vel, psi, vr, vp, a, E, J;
  vector rhat, vtmp, vper;

  for (p = btab; p < NthBody(btab, nbody); p = NextBody(p)) {
    rad = rad_m(xrandom(0.0, mtot));
    phi = gsl_spline_eval(pmspline, (double) rad, pmsplacc);
    vel = pick_v(phi);
    psi = pick_psi();
    vr = vel * rcos(psi);
    vp = vel * rsin(psi);
    Mass(p) = mtot / nbody;
    pickshell(rhat, NDIM, 1.0);
    MULVS(Pos(p), rhat, rad);
    pickshell(vtmp, NDIM, 1.0);
    a = dotvp(vtmp, rhat);
    MULVS(vper, rhat, - a);
    ADDV(vper, vper, vtmp);
    a = absv(vper);
    MULVS(vper, vper, vp / a);
    MULVS(Vel(p), rhat, vr);
    ADDV(Vel(p), Vel(p), vper);
    Phi(p) = phi;
    E = phi + 0.5 * rsqr(vel);
    J = rad * ABS(vp);
    Aux(p) = Kprime * rpow(phi1 - E, npol - 1.5) * rpow(J, 2 * mpol);
  }
  gsl_interp_accel_free(pmsplacc);
}

//  rad_m: compute radius as a function of enclosed mass.  Cannot use
//  a spline since mass increments can become very small near r_1.
//  __________________________________________________________________

real rad_m(real x)
{
  int i, j, k;
  real f;

  i = 0;
  k = nstep - 1;
  assert(mtab[i] <= x && x <= mtab[k]);
  while (k - i > 1) {
    j = (i + k) / 2;
    if (mtab[j] <= x)
      i = j;
    else
      k = j;
  }
  if (! (mtab[i] <= x && x <= mtab[k]))
    error("%s.rad_m: x=%f not between mtab[%d]=%f and mtab[%d]=%f\n",
	  getprog(), x, i, mtab[i], k, mtab[k]);
  f = (x - mtab[i]) / (mtab[k] - mtab[i]);
  return (f * rtab[k] + (1 - f) * rtab[i]);
}

//  pick_psi: pick angle between velocity and radial vectors.
//  _________________________________________________________

real pick_psi(void)
{
  static real hmax = -1.0;
  real x, psi;

  if (hmax < 0.0)
    hmax = (mpol > 0 ? rpow(2.0, mpol) : 1);
  x = vnpick(hfunc, 0.0, 1.0, hmax, "hfunc");
  psi = racos(1 - rpow(x, 1 / (mpol + 1)));
  return (xrandom(-1.0, 1.0) < 0.0 ? PI - psi : psi);
}

real hfunc(real x)
{
  return (rpow(2 - rpow(x, 1 / (mpol + 1)), mpol));
}

//  pick_v: pick magnitude of velocity.
//  ___________________________________

real pick_v(real phi)
{
  real vmax, x;
  static real gmax = -1.0;

  vmax = rsqrt(2 * phi1 - 2 * phi);
  if (npol > 1.5) {
    x = vnpick(gfunc1, 0.0, 1.0, 1.0, "gfunc1");
    return (vmax * x);
  } else {
    if (gmax < 0.0) {
      x = (npol + 4*mpol + 2.5) / (npol + 2*mpol + 0.5);
      if (0.0 <= x && x < 1.0)
	gmax = rpow(2 - x, npol - 1.5) * rpow(1 - x, 2*mpol + 2);
      else
	gmax = MAX(gfunc2(0.0), gfunc2(0.99));
    }
    x = vnpick(gfunc2, 0.0, 1.0, gmax, "gfunc2");
    return (vmax * (1 - rpow(x, 1 / (npol - 0.5))));
  }
}

real gfunc1(real x)
{
  if (rlog10(x) * (2*mpol + 2) < -36.0)
    return (0.0);
  else
    return (rpow(1 - rsqr(x), npol - 1.5) * rpow(x, 2*mpol + 2));
}

real gfunc2(real x)
{
  real y;
  
  y = rpow(x, 1 / (npol - 0.5));
  return (rpow(2 - y, npol - 1.5) * rpow(1 - y, 2*mpol + 2));
}

//  polymodel1: generate polytrope model for n = 1, m = -1.
//  _______________________________________________________

#define SQRT2 1.414214

void polymodel1(void)
{
  bodyptr p;
  real r, x, v;

  for (p = btab; p < NthBody(btab, nbody); p = NextBody(p)) {
    Mass(p) = 1.0 / nbody;
    r = xrandom(0.0, 2.0);
    pickshell(Pos(p), 3, r);
    x = SQRT2 * rcos(PI * (2.0 - xrandom(0.0, 1.0)) / 4.0);
    v = (xrandom(-1.0, 1.0) < 0.0 ? -1.0 : 1.0) *
      (1 - x*x) * rsqrt(rlog(2 / r));
    MULVS(Vel(p), Pos(p), v/r);
    Phi(p) = 0.5 * rlog(r / 2.0) - 0.5;
  }
  bodyfields[4] = NULL;				// don't output Aux data
}

//  polymodel2: generate polytrope model for n = 1/2, m = -1/2.
//  ___________________________________________________________

void polymodel2(void)
{
  error("%s: special case n = 1/2, m = -1/2 not implemented\n", getprog());
  bodyfields[4] = NULL;				// don't output Aux data
}

//  vnpick: chose from a distribution by von Neumann technique.
//  Ref: von Neumann, J. 1963. Collected Works, Vol. 5.
//  ___________________________________________________________

#define NCYC  1024

real vnpick(real (*fun)(real), real xmin, real xmax, real fmax, string name)
{
  int ncyc;
  bool warn;
  real fr, fx, x;

  ncyc = 0;
  warn = FALSE;
  fr = 1.0;
  fx = 0.0;
  while (fr > fx) {
    x = xrandom(xmin, xmax);
    fx = (*fun)(x);
    fr = xrandom(0.0, 1.1 * fmax);
    if (fx > 1.01 * fmax && ! warn) {
      eprintf("[%s.vnpick: %s(%g) = %g > %s_max = %g]\n",
	      getprog(), name, x, fx, name, fmax);
      warn = TRUE;
    }
    if (fx > 1.1 * fmax)
      error("%s.vnpick: %s(x) out of bounds\n", getprog(), name);
    ncyc++;
    if (ncyc % NCYC == 0)
      eprintf("[%s.vnpick: %d cycles picking %s(x)]\n",getprog(), ncyc, name);
  }
  return (x);
}

//  berank: compare binding energies of particles, for sorting.
//  ___________________________________________________________

int berank(const void *p1, const void *p2)
{
  real e1 = 0.5 * rsqr(absv(Vel((bodyptr) p1))) + Phi((bodyptr) p1);
  real e2 = 0.5 * rsqr(absv(Vel((bodyptr) p2))) + Phi((bodyptr) p2);
  return (e1 < e2 ? -1 : e1 > e2 ? 1 : 0);
}
