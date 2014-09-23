/*
 * GSPDISK.C: set up an exponential disk in a gsp model.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"

#include <gsl/gsl_sf_bessel.h>

string defv[] = {		";Make exponential disk in a gsp model",
    "grav=",			";Input gsp file sets spheroid gravity.",
				";If blank, just use disk self-gravity.",
    "out=",			";Output N-body file with disk model.",
				";Contains " MassTag ", " PosTag ", "
				  VelTag " data.",
    "mdisk=0.1875",		";Exponential disk mass (or masses).",
				";Can take two comma-separated values;",
				";uses 1st for N-body model, adds 2nd to",
				";grav. field.  If so, must also specify",
				";two values for alpha and epsilon.",
    "alpha=12.0",		";Inverse radial scale length(s)",
    "epsilon=-1",		";Plummer smoothing parameter(s).",
				";No smoothing is used if epsilon < 0.",
    "zdisk=0.01",		";Vertical scale height; sets dispersions",
    "mu=2.0",			";Ratio of radial to vertical dispersion.",
				";Constant unless r_mu > 0.",
    "r_mu=-1",			";Scale radius for mu(R) function.",
				";If r_mu > 0 then mu -> 1 as R -> 0.",
    "eta=-1",			";Velocity distribution parameter.",
				";Use Gaussian distribution if eta < 0.",
    "rcut=1.0",			";Outer disk cutoff radius",
    "nlist=256",		";Number of radii in disk listing",
    "ndisk=12288",		";Number of disk particles",
    "seed=54321",		";Seed for random number generator",
    "VERSION=1.4",		";Josh Barnes  02 September 2010",
    NULL,
};

/* Global parameters. */

real mdisk1, mdisk2;			/* total mass of expo. disk(s)      */
real alpha1, alpha2;			/* inverse radial scale length(s)   */
real epsilon1, epsilon2;		/* value for smoothing length(s)    */
real zdisk;				/* vertical scale height            */
real mu;				/* ratio of sig_r to sig_z          */
real r_mu;				/* scale radius for mu(R) function  */
real eta;				/* shape parameter for vel. dist.   */
real rcut;				/* cut-off radius for disk model    */

/* Global tables and data structures. */

#define NTAB  (256 + 1)

real mdtab[NTAB];			/* use disk mass as indp var        */
real rdtab[4*NTAB];			/* radius as fcn of mass            */
real vctab[4*NTAB];			/* circ. velocity as fcn of radius  */

gsprof *sphr = NULL;			/* sphroid mass as fcn of radius    */

string bodyfields[] = { PosTag, VelTag, MassTag, NULL };

bodyptr disk = NULL;			/* array of disk particles          */
int ndisk;				/* number of particles in the disk  */

/* Function prototypes. */

void readgsp(void);
void writemodel(void);
void setprof(void);
real gdisk(real);
real dgdisk(real, real, real, real);
real simpson(real (*)(real, real, real, real),
	     real, real, real, real, real, real);
void listdisk(int);
void makedisk(void);
real ratanh(real);
real pickdist(real, real);
real bessel_k(double nu, double x);

int main(int argc, string argv[])
{
  float tmp1, tmp2;

  initparam(argv, defv);
  if (sscanf(getparam("mdisk"), "%f,%f", &tmp1, &tmp2) == 2) {
    mdisk1 = tmp1;
    mdisk2 = tmp2;
    if (sscanf(getparam("alpha"), "%f,%f", &tmp1, &tmp2) != 2)
      error("%s: must specify two alpha values\n", getargv0());
    alpha1 = tmp1;
    alpha2 = tmp2;
    if (sscanf(getparam("epsilon"), "%f,%f", &tmp1, &tmp2) != 2)
      error("%s: must specify two epsilon values\n", getargv0());
    epsilon1 = tmp1;
    epsilon2 = tmp2;
  } else {
    mdisk1 = getdparam("mdisk");
    mdisk2 = 0.0;
    alpha1 = getdparam("alpha");
    alpha2 = 12.0;			 /* nonzero value stops div by zero */
    epsilon1 = getdparam("epsilon");
    epsilon2 = -1.0;
  }
  zdisk = getdparam("zdisk");
  mu = getdparam("mu");
  r_mu = getdparam("r_mu");
  eta = getdparam("eta");
  rcut = getdparam("rcut");
  readgsp();
  setprof();
  layout_body(bodyfields, Precision, NDIM);
  ndisk = getiparam("ndisk");
  disk = (bodyptr) allocate(ndisk * SizeofBody);
  srandom(getiparam("seed"));
  if (getiparam("nlist") > 0)
    listdisk(getiparam("nlist"));
  makedisk();
  writemodel();
  return (0);
}

/*
 * READGSP: read sphr GSP from input file.
 */

void readgsp(void)
{
  stream istr;

  if (! strnull(getparam("grav"))) {
    istr = stropen(getparam("grav"), "r");
    get_history(istr);
    sphr = get_gsprof(istr);
    strclose(istr);
  }
}

/*
 * WRITEMODEL: write N-body model to output file.
 */

void writemodel(void)
{
  stream ostr;
  real tsnap = 0.0;

  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_snap(ostr, &disk, &ndisk, &tsnap, bodyfields);
    strclose(ostr);
  }
}

/*
 * SETPROF: initialize disk tables for radius and circular velocity.
 */

void setprof(void)
{
  int j;
  real r, msphr;

  rdtab[0] = mdtab[0] = vctab[0] = 0.0;
  for (j = 1; j < NTAB; j++) {
    r = rcut * rpow(((real) j) / (NTAB - 1), 2.0);
    rdtab[j] = r;
    mdtab[j] = 1 - rexp(- alpha1 * r) - alpha1 * r * rexp(- alpha1 * r);
    msphr = (sphr != NULL ? mass_gsp(sphr, r) : 0.0);
    vctab[j] = rsqrt(msphr / r - gdisk(r) * r);
  }
  eprintf("[%s: rcut = %8.4f/alpha  M(rcut) = %8.6f mdisk]\n",
	  getargv0(), rdtab[NTAB-1] * alpha1, mdtab[NTAB-1]);
  if ((mdtab[0] == mdtab[1]) || (mdtab[NTAB-2] == mdtab[NTAB-1]))
      error("%s: disk mass table is degenerate\n", getargv0());
  spline(&rdtab[NTAB], &mdtab[0], &rdtab[0], NTAB);	/* for r_d = r_d(m) */
  spline(&vctab[NTAB], &rdtab[0], &vctab[0], NTAB);	/* for v_c = v_c(r) */
}

/*
 * GDISK: compute radial acceleration due to exponential disks.
 */

#define KMAX 1000.0
#define STEP 0.1

real gdisk(real r)
{
  real x1 = 0.5 * alpha1 * r, x2 = 0.5 * alpha2 * r, a1, a2;;

  if (epsilon1 >= 0.0)				/* compute smoothed accel.  */
    a1 = - mdisk1 * rqbe(alpha1) *
      simpson(dgdisk, r, alpha1, epsilon1, 0, KMAX*alpha1, STEP*alpha1);
  else						/* use exact expression     */
    a1 = - mdisk1 * rqbe(alpha1) *
      r * (bessi0(x1) * bessk0(x1) - bessi1(x1) * bessk1(x1)) / 2;
  if (epsilon2 >= 0.0)
    a2 = - mdisk2 * rqbe(alpha2) *
      simpson(dgdisk, r, alpha2, epsilon2, 0, KMAX*alpha2, STEP*alpha2);
  else
    a2 = - mdisk2 * rqbe(alpha2) *
      r * (bessi0(x2) * bessk0(x2) - bessi1(x2) * bessk1(x2)) / 2;
  return (a1 + a2);	       
}

/*
 * DGDISK: integrand for radial acceleration with smoothing. Note: potential
 * at z = 0 with smoothing equals potential at z = epsilon without smoothing.
 */

real dgdisk(real k, real r, real alpha, real eps)
{
  return (rexp(- k*eps) * j1(k*r) * k / rsqrt(rqbe(alpha*alpha + k*k)));
}

/*
 * SIMPSON: integrate given function using Simpson's rule.
 */

real simpson(real (*integrand)(real, real, real, real),
	     real param1, real param2, real param3,
	     real xlow, real xhigh, real step0)
{
  int nstep, i;
  real step1, x;
  double v1, v2, v4;

  nstep = 1 + 2 * (int) rceil(0.5 * (xhigh - xlow) / step0);
  step1 = (xhigh - xlow) / (nstep - 1);
  v1 = v2 = v4 = 0.0;
  for (i = 1; i <= nstep; i++) {
    x = xlow + (i - 1) * step1;
    if (i == 1 || i == nstep)
      v1 = v1 + (*integrand)(x, param1, param2, param3);
    else if (i % 2 == 0)
      v4 = v4 + (*integrand)(x, param1, param2, param3);
    else
      v2 = v2 + (*integrand)(x, param1, param2, param3);
  }
  return (step1 * (v1 + 4.0*v4 + 2.0*v2) / 3.0);
}

/*
 * LISTDISK: print table listing disk parameters.
 */

void listdisk(int nlist)
{
  int i;
  real r, phi, vcir, omega, Adisk, kappa, sigma1, sigma2,
       mu_eff, sig_r, sig_p, sig_z, vrad, vorb2, vorb;

  printf("#%5s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n",
	 "r", "vcir", "omega", "kappa", "sig_z", "sig_r", "sig_p", "vorb",
	 "Q", "rhomid", "fmax");
  for (i = 1; i <= nlist; i++) {			/* loop over radii  */
    r = (i * rcut) / ((real) nlist);
    vcir = seval(r, &rdtab[0], &vctab[0], &vctab[NTAB], NTAB);
    omega = vcir / r;
    Adisk = (omega - spldif(r, &rdtab[0], &vctab[0], &vctab[NTAB], NTAB)) / 2;
    if (omega - Adisk < 0.0)
      error("%s: kappa undefined (omega - Adisk < 0)\n"
	    "  r, omega, Adisk = %f %f %f\n", getargv0(), r, omega, Adisk);
    kappa = 2 * rsqrt(rsqr(omega) - Adisk * omega);
    sigma1 = rsqr(alpha1) * mdisk1 * rexp(- alpha1 * r) / TWO_PI;
    sigma2 = rsqr(alpha2) * mdisk2 * rexp(- alpha2 * r) / TWO_PI;
    mu_eff = (r_mu>0 ? 1 + (mu - 1) * (r / (r + r_mu)) : mu);
    sig_z = rsqrt(PI * (sigma1 + sigma2) * zdisk);
    sig_r = mu_eff * sig_z;
    sig_p = (0.5 * kappa / omega) * sig_r;
    vorb2 = rsqr(vcir) + rsqr(sig_r) * (1 - 2 * alpha1 * r) - rsqr(sig_p) +
      (r_mu>0 ? rsqr(sig_z) * r * mu_eff*(2*mu-2)*r_mu/rsqr(r+r_mu) : 0);
    vorb = rsqrt(MAX(vorb2, 0.0));
    printf("%6.4f %6.4f %6.2f %6.2f %6.4f %6.4f %6.4f %6.4f "
	   "%6.3f %6.1f %6.1f\n",
	   r, vcir, omega, kappa, sig_z, sig_r, sig_p, vorb,
	   kappa * sig_r / (3.358*(sigma1+sigma2)), sigma1/(2*zdisk),
	   sigma1 / (2*zdisk * rsqrt(rqbe(2*PI)) * sig_z*sig_r*sig_p));
  }
}

/*
 * MAKEDISK: create realization of disk.
 */

void makedisk(void)
{
  real m, r, phi, vcir, omega, Adisk, kappa, sigma1, sigma2,
       mu_eff, sig_r, sig_p, sig_z, vrad, vorb2, vorb, vphi;
  double Trr = 0.0, Tpp = 0.0, Tzz = 0.0;
  int i;
  bodyptr dp;

  for (i = 0; i < ndisk; i++) {			/* loop initializing bodies */
    m = mdtab[NTAB-1] * ((real) i + 0.5) / ndisk;
    r = seval(m, &mdtab[0], &rdtab[0], &rdtab[NTAB], NTAB);
    vcir = seval(r, &rdtab[0], &vctab[0], &vctab[NTAB], NTAB);
    omega = vcir / r;
    Adisk = (omega - spldif(r, &rdtab[0], &vctab[0], &vctab[NTAB], NTAB)) / 2;
    if (omega - Adisk < 0.0)
      error("%s: kappa undefined (omega - Adisk < 0)\n"
	    "  r, omega, Adisk = %f %f %f\n", getargv0(), r, omega, Adisk);
    kappa = 2 * rsqrt(rsqr(omega) - Adisk * omega);
    sigma1 = rsqr(alpha1) * mdisk1 * rexp(- alpha1 * r) / TWO_PI;
    sigma2 = rsqr(alpha2) * mdisk2 * rexp(- alpha2 * r) / TWO_PI;
    mu_eff = (r_mu>0 ? 1 + (mu - 1) * (r / (r + r_mu)) : mu);
    sig_z = rsqrt(PI * (sigma1 + sigma2) * zdisk);
    sig_r = mu_eff * sig_z;
    sig_p = (0.5 * kappa / omega) * sig_r;
    vorb2 = rsqr(vcir) + rsqr(sig_r) * (1 - 2 * alpha1 * r) - rsqr(sig_p) +
      (r_mu>0 ? rsqr(sig_z) * r * mu_eff*(2*mu-2)*r_mu/rsqr(r+r_mu) : 0);
    vorb = rsqrt(MAX(vorb2, 0.0));
    dp = NthBody(disk, i);			/* set up ptr to disk body  */
    Mass(dp) = mdisk1 / ndisk;
    phi = xrandom(0.0, TWO_PI);
    Pos(dp)[0] = r * rsin(phi);
    Pos(dp)[1] = r * rcos(phi);
    Pos(dp)[2] = zdisk * ratanh(xrandom(-1.0, 1.0));
    vrad = (eta > 0 ? pickdist(eta, sig_r) : grandom(0.0, sig_r));
    vphi = (eta > 0 ? pickdist(eta, sig_p) : grandom(0.0, sig_p)) + vorb;
    Vel(dp)[0] = vrad * rsin(phi) + vphi * rcos(phi);
    Vel(dp)[1] = vrad * rcos(phi) - vphi * rsin(phi);
    Vel(dp)[2] = grandom(0.0, sig_z);
    Trr += Mass(dp) * rsqr(sig_r) / 2;
    Tpp += Mass(dp) * (rsqr(vorb) + rsqr(sig_p)) / 2;
    Tzz += Mass(dp) * rsqr(sig_z) / 2;
  }
  eprintf("[%s: Trr = %f  Tpp = %f  Tzz = %f]\n", getargv0(), Trr, Tpp, Tzz);
}

/*
 * RATANH: return hyperbolic arc tangent.
 */

real ratanh(real x)
{
  return (0.5 * rlog((1.0 + x) / (1.0 - x)));
}

/*
 * PICKDIST: pick value from modified gaussian distribution.
 */

#define YMAX  2.5
#define fmap(x)  ((x) / (1 - (x)*(x)))

real pickdist(real eta, real sigma)
{
  static real eta0 = -1.0, sigcorr;
  int niter;
  real x, y, q;

  if (eta != eta0) {
    sigcorr =
      rsqrt(8 * eta /
	      (bessel_k(0.75, 1/(32*eta)) / bessel_k(0.25, 1/(32*eta)) - 1));
    eprintf("[%s: sigma correction factor = %f]\n", getargv0(), sigcorr);
    eta0 = eta;
  }
  niter = 0;
  do {
    x = xrandom(-1.0, 1.0);
    y = xrandom(0.0, YMAX);
    q = rexp(- 0.5 * rsqr(fmap(x)) - eta * rsqr(rsqr(fmap(x)))) *
      (1 + x*x) / rsqr(1 - x*x);
    if (q > YMAX)				/* should not ever happen   */
      error("%s: guess out of bounds\n  x = %f  q = %f > %f\n",
	    getargv0(), x, q, YMAX);
    niter++;
    if (niter > 1000)
      error("%s: 1000 iterations without success\n", getargv0());
  } while (y > q || x*x == 1);			/* 2nd test prevents infty  */
  return (sigcorr * sigma * fmap(x));
}

/*
 * BESSEL_K: compute modified Bessel function K of arbitrary order.
 */

real bessel_k(double nu, double x)
{
  gsl_sf_result res;
  int stat;

  stat = gsl_sf_bessel_Knu_e(nu, x, &res);
  if (stat != 0)
    error("%s.bessel_k: error status: %s\n", getargv0(), gsl_strerror(stat));
  return (res.val);
}
