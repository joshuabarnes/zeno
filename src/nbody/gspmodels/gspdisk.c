/*
 * gspdisk.c: set up an exponential disk in a gsp model.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <assert.h>

string defv[] = {		";Make exponential disk in a gsp model",
  "grav=",			";Input gsp file sets spheroid gravity.",
				";If blank, just use disk self-gravity.",
  "out=",			";Output N-body file with disk model.",
				";Contains " MassTag ", " PosTag ", " VelTag
#if defined(GAS)
				",",
				";" UinternTag ", " TypeTag ", " AuxTag
#endif
				" data.",
  "mdisk=0.1875",		";Exponential disk mass (or masses).",
				";Can take many comma-separated values;",
				";uses first to generate model,",
				";adds rest to compute rotation curve.",
  "alpha=12.0",			";Inverse radial scale length(s).",
				";Can take many comma-separated values.",
  "epsilon=-1",			";Smoothing parameter for disk(s).",
				";Can take many comma-separated values.",
				";No smoothing is used if epsilon < 0.",
#if !defined(GAS)
  "zdisk=0.01",			";Vertical scale height; sets dispersions",
  "mu=2.0",			";Ratio of radial to vertical dispersion.",
				";Constant unless r_mu > 0.",
  "r_mu=-1",			";Scale radius for mu(R) function.",
				";If r_mu > 0 then mu -> 1 as R -> 0.",
  "eta=-1",			";Velocity distribution parameter.",
				";Use Gaussian distribution if eta < 0.",
  "rcut=1.0",			";Outer disk cutoff radius",
  "nlist=256",			";Number of radii in disk listing",
#else
  "uint=0.014",                 ";Internal energy for SPH particles",
  "gamma=5.0/3.0",              ";Ratio of specific heats",
  "type=0x60",			";Body type for SPH calculation",
  "rcut=1.0",			";Outer disk cutoff radius",
#endif
  "ndisk=12288",		";Number of disk particles",
  "seed=54321",			";Seed for random number generator",
  "VERSION=2.0",		";Josh Barnes & Kelly Blumenthal  24 July 2017",
  NULL,
};

// Function prototypes.

void readgsp(void);
void writemodel(void);
void setprof(void);
double gdisk(double);
double dgdisk(double, void *);
void listdisk(int);
void makedisk(void);
double pickdist(double, double);
string *padlist(string *old, int len);

// Global parameters.

#define MAXMODEL  4
double mdisk[MAXMODEL];			// total mass of expo. disk(s)
double alpha[MAXMODEL];			// inverse radial scale length(s)
double epsilon[MAXMODEL];		// value for smoothing length(s)
int nmodel;				// number of disk models

double zdisk;				// vertical scale height of 1st disk
double mu;				// ratio of sig_r to sig_z
double r_mu;				// scale radius for mu(R) function
double eta;				// shape parameter for vel. dist.
double rcut;				// cut-off radius for disk model
double uinternal;                       // internal energy for SPH particles
double gam;                             // Ratio of specific heats
byte type;                              // body type for SPH particles

// Global arrays and data structures.

#define NTAB  (256 + 1)

double mdtab[NTAB];			// use disk mass as indp var
double rdtab[NTAB];			// radius as fcn of mass
double vctab[NTAB];			// circ. velocity as fcn of radius
gsl_interp *rm_spline;			// interpolator for r(m)
gsl_interp *vr_spline;			// interpolator for v(r)

gsprof *sphr = NULL;			// sphroid mass as fcn of radius

string bodyfields[] = {
#if !defined(GAS)
  PosTag, VelTag, MassTag, NULL,
#else
  PosTag, VelTag, MassTag, UinternTag, TypeTag, AuxTag, NULL,
#endif
};

bodyptr disk = NULL;			// array of disk particles
int ndisk;				// number of particles in the disk

int main(int argc, string argv[])
{
  string *mdisklist, *alphalist, *epsillist;

  initparam(argv, defv);
  mdisklist = burststring(getparam("mdisk"), ",");
  nmodel = xstrlen(mdisklist, sizeof(string)) - 1;
  if (nmodel < 1)
    error("%s: must specify at least one mdisk\n", getprog());
  if (nmodel > MAXMODEL)
    error("%s: can't have more than %d disks\n", getprog(), MAXMODEL);
  alphalist = padlist(burststring(getparam("alpha"), ","), nmodel);
  epsillist = padlist(burststring(getparam("epsilon"), ","), nmodel);
  if (alphalist == NULL || epsillist == NULL)
    error("%s: must specify value(s) for alpha, epsilon\n", getprog());
  for (int i = 0; i < nmodel; i++) {
    mdisk[i] = strtod(mdisklist[i], (char **) NULL);
    alpha[i] = strtod(alphalist[i], (char **) NULL);
    epsilon[i] = strtod(epsillist[i], (char **) NULL);
  }
#if !defined(GAS)
  zdisk = getdparam("zdisk");
  mu = getdparam("mu");
  r_mu = getdparam("r_mu");
  eta = getdparam("eta");
#else
  uinternal = getdparam("uint");
  gam = getdparam("gamma");
  type = (0x7f & getiparam("type"));
#endif
  rcut = getdparam("rcut");
  readgsp();
  setprof();
  layout_body(bodyfields, Precision, NDIM);
  ndisk = getiparam("ndisk");
  disk = (bodyptr) allocate(ndisk * SizeofBody);
  init_random(getiparam("seed"));
#if !defined(GAS)
  if (getiparam("nlist") > 0)
    listdisk(getiparam("nlist"));
#endif
  makedisk();
  writemodel();
  fflush(NULL);
  return 0;
}

//  readgsp: read sphr GSP from input.
//  __________________________________

void readgsp(void)
{
  stream istr;

  if (! strnull(getparam("grav"))) {
    istr = stropen(getparam("grav"), "r");
    get_history(istr);
    sphr = gsp_read(istr);
  }
}

//  writemodel: write N-body model to output.
//  _________________________________________

void writemodel(void)
{
  stream ostr;
  real tsnap = 0.0;

  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_snap(ostr, &disk, &ndisk, &tsnap, bodyfields);
  }
}

//  setprof: initialize disk tables for radius and circular velocity.
//  _________________________________________________________________

void setprof(void)
{
  double r, msphr;

  rdtab[0] = mdtab[0] = vctab[0] = 0.0;
  for (int j = 1; j < NTAB; j++) {
    rdtab[j] = r = rcut * gsl_pow_2(((double) j) / (NTAB - 1));
    mdtab[j] = 1 - exp(- alpha[0] * r) - alpha[0] * r * exp(- alpha[0] * r);
    msphr = (sphr != NULL ? gsp_mass(sphr, r) : 0.0);
    vctab[j] = sqrt(msphr / r - gdisk(r) * r);
  }
  eprintf("[%s: rcut = %8.4f/alpha  M(rcut) = %8.6f mdisk]\n",
	  getprog(), rdtab[NTAB-1] * alpha[0], mdtab[NTAB-1]);
  if ((mdtab[0] == mdtab[1]) || (mdtab[NTAB-2] == mdtab[NTAB-1]))
      error("%s: disk mass table is degenerate\n", getprog());
  rm_spline = gsl_interp_alloc(gsl_interp_akima, NTAB);
  gsl_interp_init(rm_spline, mdtab, rdtab, NTAB);
  vr_spline = gsl_interp_alloc(gsl_interp_akima, NTAB);
  gsl_interp_init(vr_spline, rdtab, vctab, NTAB);
}

//  gdisk: compute radial acceleration due to exponential disks.
//  ____________________________________________________________

#define KMAX 1000.0
#define NWKSP  1000
#define ERRTOL 1.0e-7

double gdisk(double r)
{
  static gsl_integration_workspace *wksp = NULL;
  double params[3], res, abserr, a, x, atot = 0;
  gsl_function dgdisk_func = { &dgdisk, params };

  for (int j = 0; j < nmodel; j++) {
    if (epsilon[j] >= 0.0) {			// compute smoothed accel.
      if (wksp == NULL)				// get workspace ready
	wksp = gsl_integration_workspace_alloc(NWKSP);
      params[0] = r;				// init params for dgdisk
      params[1] = alpha[j];
      params[2] = epsilon[j];
      assert(gsl_integration_qag(&dgdisk_func, 0.0, KMAX * alpha[j],
				 0.0, ERRTOL, NWKSP, GSL_INTEG_GAUSS61,
				 wksp, &res, &abserr) == GSL_SUCCESS);
    } else {					// use analytic expression
      x = 0.5 * alpha[j] * r;
      res = r * (bessel_I0(x)*bessel_K0(x) - bessel_I1(x)*bessel_K1(x)) / 2;
    }
    atot = atot - mdisk[j] * gsl_pow_3(alpha[j]) * res;
  }
  return atot;	       
}

//  dgdisk: integrand for radial acceleration with smoothing. Note: potential
//  at z = 0 with smoothing equals potential at z = epsilon without smoothing.
//  __________________________________________________________________________

double dgdisk(double k, void *params)
{
  double r = ((double *) params)[0];
  double alpha = ((double *) params)[1];
  double eps = ((double *) params)[2];
  
  return (exp(- k*eps) * bessel_J1(k*r) * k / pow(alpha*alpha + k*k, 1.5));
}

#if !defined(GAS)

//  listdisk: print table listing disk parameters.
//  ______________________________________________

void listdisk(int nlist)
{
  double r, vcir, dvdr, Omega, Aoort, kappa, Sigma1, Sigma;
  double mu_eff, sig_z, sig_r, sig_p, vorb2, vorb;

  printf("#%5s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n",
	 "r", "vcir", "Omega", "kappa", "sig_z", "sig_r", "sig_p", "vorb",
	 "Q", "rhomid", "fmax");
  for (int i = 1; i <= nlist; i++) {		// loop over table radii
    r = (i * rcut) / ((double) nlist);
    vcir = gsl_interp_eval(vr_spline, rdtab, vctab, r, NULL);
    dvdr = gsl_interp_eval_deriv(vr_spline, rdtab, vctab, r, NULL);
    Omega = vcir / r;
    Aoort = (Omega - dvdr) / 2;
    if (Omega - Aoort < 0.0)
      error("%s.listdisk: kappa undefined (Omega - Aoort < 0)\n"
	    "  r, Omega, Aoort = %f %f %f\n", getprog(), r, Omega, Aoort);
    kappa = 2 * sqrt(Omega*Omega - Aoort * Omega);
    Sigma1 = alpha[0]*alpha[0] * mdisk[0] * exp(- alpha[0] * r) / (2*M_PI);
    Sigma = Sigma1;				// sum total surface density
    for (int j = 1; j < nmodel; j++)
      Sigma += alpha[j]*alpha[j] * mdisk[j] * exp(- alpha[j] * r) / (2*M_PI);
    mu_eff = (r_mu>0 ? 1 + (mu - 1) * (r / (r + r_mu)) : mu);
    sig_z = sqrt(M_PI * Sigma * zdisk);		// model as isothermal sheet
    sig_r = mu_eff * sig_z;			// apply ratio given by user
    sig_p = (0.5 * kappa / Omega) * sig_r;	// use epicyclic relationship
    vorb2 = vcir*vcir + sig_r*sig_r * (1 - 2 * alpha[0] * r) - sig_p*sig_p +
      (r_mu>0 ? sig_z*sig_z * r * mu_eff*(2*mu-2)*r_mu/(r+r_mu)*(r+r_mu) : 0);
    vorb = sqrt(MAX(vorb2, 0.0));
    printf("%6.4f %6.4f %6.2f %6.2f %6.4f %6.4f %6.4f %6.4f "
	   "%6.3f %6.1f %6.1f\n",
	   r, vcir, Omega, kappa, sig_z, sig_r, sig_p, vorb,
	   kappa * sig_r / (3.358 * Sigma), Sigma1 / (2*zdisk),
	   Sigma1 / (2*zdisk * pow(2*M_PI, 1.5) * sig_z*sig_r*sig_p));
  }
}
#endif

//  makedisk: create realization of disk.
//  _____________________________________

void makedisk(void)
{
  bodyptr dp;
  double m, r, vcir, phi, Omega, dvdr, kappa, Sigma, K, zgas;
  double mu_eff, sig_z, sig_r, sig_p, vorb2, vorb, vrad, vphi;
  double Trr = 0.0, Tpp = 0.0, Tzz = 0.0;

  for (int i = 0; i < ndisk; i++) {		// loop initializing bodies
    dp = NthBody(disk, i);
    m = mdtab[NTAB-1] * ((double) i + 0.5) / ndisk;
    r = gsl_interp_eval(rm_spline, mdtab, rdtab, m, NULL);
    vcir = gsl_interp_eval(vr_spline, rdtab, vctab, r, NULL);
    Mass(dp) = mdisk[0] / ndisk;
    phi = xrandom(0.0, 2 * M_PI);
#if !defined(GAS)
    Omega = vcir / r;
    dvdr = gsl_interp_eval_deriv(vr_spline, rdtab, vctab, r, NULL);
    if (Omega + dvdr < 0.0)			// should NOT happen!
      error("%s.makedisk: kappa undefined because Omega + dvdr < 0\n"
	    "  r, Omega, dvdr = %f %f %f\n", getprog(), r, Omega, dvdr);
    kappa = 2 * sqrt(Omega * (Omega + dvdr) / 2);
    Sigma = 0.0;				// sum total surface density
    for (int j = 0; j < nmodel; j++)
      Sigma += alpha[j]*alpha[j] * mdisk[j] * exp(- alpha[j] * r) / (2*M_PI);
    mu_eff = (r_mu>0 ? 1 + (mu - 1) * (r / (r + r_mu)) : mu);
    sig_z = sqrt(M_PI * Sigma * zdisk);		// model as isothermal sheet
    sig_r = mu_eff * sig_z;			// apply ratio given by user
    sig_p = (0.5 * kappa / Omega) * sig_r;	// use epicyclic relationship
    vorb2 = vcir*vcir + sig_r*sig_r * (1 - 2 * alpha[0] * r) - sig_p*sig_p +
      (r_mu>0 ? sig_z*sig_z * r * mu_eff*(2*mu-2)*r_mu/(r+r_mu)*(r+r_mu) : 0);
    vorb = sqrt(MAX(vorb2, 0.0));
    vrad = (eta > 0 ? pickdist(eta, sig_r) : grandom(0.0, sig_r));
    vphi = (eta > 0 ? pickdist(eta, sig_p) : grandom(0.0, sig_p)) + vorb;
    Pos(dp)[0] = r * sin(phi);
    Pos(dp)[1] = r * cos(phi);
    Pos(dp)[2] = zdisk * atanh(xrandom(-1.0, 1.0));
    Vel(dp)[0] = vrad * sin(phi) + vphi * cos(phi);
    Vel(dp)[1] = vrad * cos(phi) - vphi * sin(phi);
    Vel(dp)[2] = grandom(0.0, sig_z);
    Trr += Mass(dp) * sig_r*sig_r / 2;
    Tpp += Mass(dp) * (vorb*vorb + sig_p*sig_p) / 2;
    Tzz += Mass(dp) * sig_z*sig_z / 2;
#else
    K = (sphr != NULL ? gsp_mass(sphr, r) / (r*r*r) : 0.0);
    if (nmodel > 1)				// stellar disk present?
      K += (mdisk[1] * alpha[1]*alpha[1] / epsilon[1]) * exp(-alpha[1]*r);
    zgas = (K > 0 ? sqrt(uinternal * (gam - 1) / K) : 0.0);
    Pos(dp)[0] = r * sin(phi);
    Pos(dp)[1] = r * cos(phi);
    Pos(dp)[2] = grandom(0, zgas);
    Vel(dp)[0] = vcir * cos(phi);
    Vel(dp)[1] = - vcir * sin(phi);
    Vel(dp)[2] = 0.0;
    Uintern(dp) = uinternal;
    Type(dp) = type;
    Aux(dp) = K;				// save force constant
    Tpp += Mass(dp) * vcir*vcir / 2;
#endif    
  }

  eprintf("[%s.makedisk: Trr = %f  Tpp = %f  Tzz = %f]\n",
	  getprog(), Trr, Tpp, Tzz);
}

//  pickdist: pick value from modified gaussian distribution.
//  _________________________________________________________

#define YMAX  2.5
#define fmap(x)  ((x) / (1 - (x)*(x)))

double pickdist(double eta, double sigma)
{
  static double eta0 = -1.0, sigcorr;
  int niter;
  double x, y, q;

  if (eta != eta0) {
    sigcorr =
      sqrt(8 * eta /
	      (bessel_K(0.75, 1/(32*eta)) / bessel_K(0.25, 1/(32*eta)) - 1));
    eprintf("[%s.pickdist: sigma correction factor = %f]\n",
	    getprog(), sigcorr);
    eta0 = eta;
  }
  niter = 0;
  do {
    x = xrandom(-1.0, 1.0);
    y = xrandom(0.0, YMAX);
    q = exp(- 0.5 * gsl_pow_2(fmap(x)) - eta * gsl_pow_4(fmap(x))) *
      (1 + x*x) / gsl_pow_2(1 - x*x);
    if (q > YMAX)				// should not ever happen
      error("%s.pickdist: guess out of bounds\n  x = %f  q = %f > %f\n",
	    getprog(), x, q, YMAX);
    niter++;
    if (niter >= 1000)
      error("%s.pickdist: 1000 iterations without success\n", getprog());
  } while (y > q || x*x == 1);			// 2nd test prevents infty
  return (sigcorr * sigma * fmap(x));
}

//  padlist: pad string list to specified length.
//  _____________________________________________

string *padlist(string *old, int newlen)
{
  int oldlen = xstrlen(old, sizeof(string)) - 1;
  string *new;
  
  if (oldlen == 0)
    return (NULL);
  if (oldlen >= newlen)
    return (old);
  new = (string *) allocate((newlen + 1) * sizeof(string));
  for (int i = 0; i <= newlen; i++)
    new[i] = (i < newlen ? old[MIN(i, oldlen - 1)] : NULL);
  return new;
}
