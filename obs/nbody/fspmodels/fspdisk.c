/*
 * FSPDISK.C: set up an exponential disk embedded in a fsp halo.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "fsp.h"

string defv[] = {		";Make exponential disk in a fsp halo",
    "fsp=",			";Input fsp for halo mass profile",
    "out=",			";Output N-body model of disk",
    "mdisk=0.1875",		";Total disk mass",
    "alpha=12.0",		";Inverse exponential scale length",
    "z0=0.01",			";Vertical scale height",
    "epsilon=-1",		";Plummer smoothing parameter.",
				";No smoothing (cf. v1.1) if < 0.",
    "mu=2.0",			";Ratio of sig_r/sig_z.",
				";Constant unless r_mu > 0.",
    "r_mu=-1",			";Scale radius for mu(R) function.",
				";If > 0 then mu -> 1 as R -> 0.",
    "eta=-1",			";Velocity distribution parameter.",
				";Use Gaussian distrib. if <= 0.",
    "rcut=1.0",			";Outer disk cutoff radius",
    "ndisk=12288",		";Number of disk particles",
    "seed=54321",		";Seed for random number generator",
    "zerocm=true",		";If TRUE, translate CM to origin",
    "VERSION=1.2",		";Josh Barnes  27 March 2002",
    NULL,
};

/* Function prototypes. */

void readfsp(void);
void writemodel(void);
void setprof(void);
real gdisk(real);
real dgdisk(real, real);
real simpson(real (*)(real, real), real, real, real, real);
void makedisk(void);
real bessel_k(real, real);
real ratanh(real);
real pickdist(real, real);

/* Global parameters. */

real mdisk;				/* total mass of exponential disk   */
real alpha;				/* inverse exponential scale length */
real z0;				/* vertical disk scale height       */
real mu;				/* ratio of sig_r to sig_z          */
real rcut;				/* cut-off radius for disk model    */
real epsilon;				/* value for smoothing parameter    */
real eta;				/* shape parameter for vel. dist.   */
real r_mu;				/* scale radius for mu(R) function  */
int ndisk;				/* number of particles in the disk  */

/* Global tables and data structures. */

#define NTAB  (256 + 1)

real mdtab[NTAB];			/* use   disk mass as indp var      */
real rdtab[4*NTAB];			/* radius as fcn of mass            */
real vctab[4*NTAB];			/* circ. velocity as fcn of radius  */

fsprof *halo = NULL;			/* def. halo mass as fcn of radius  */

string bodyfields[] = { PosTag, VelTag, MassTag, NULL };

bodyptr disk = NULL;			/* array of disk particles          */

int main(int argc, string argv[])
{
    initparam(argv, defv);
    mdisk = getdparam("mdisk");
    alpha = getdparam("alpha");
    z0 = getdparam("z0");
    mu = getdparam("mu");
    rcut = getdparam("rcut");
    epsilon = getdparam("epsilon");
    eta = getdparam("eta");
    r_mu = getdparam("r_mu");
    readfsp();
    setprof();
    layout_body(bodyfields, Precision, NDIM);
    ndisk = getiparam("ndisk");
    disk = (bodyptr) allocate(ndisk * SizeofBody);
    srandom(getiparam("seed"));
    makedisk();
    if (getbparam("zerocm"))
	snapcenter(disk, ndisk, MassField.offset);
    writemodel();
    return (0);
}

/*
 * READFSP: read halo FSP from input file.
 */

void readfsp(void)
{
    stream istr;

    if (! strnull(getparam("fsp"))) {
	istr = stropen(getparam("fsp"), "r");
	get_history(istr);
	halo = get_fsprof(istr);
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
    real r, mhalo;

    rdtab[0] = mdtab[0] = vctab[0] = 0.0;
    for (j = 1; j < NTAB; j++) {
	r = rcut * rpow(((real) j) / (NTAB - 1), 2.0);
	rdtab[j] = r;
	mdtab[j] = 1 - rexp(- alpha * r) - alpha * r * rexp(- alpha * r);
	mhalo = (halo != NULL ? mass_fsp(halo, r) : 0.0);
	vctab[j] = rsqrt(mhalo / r - gdisk(r) * r);
    }
    eprintf("[%s: rcut = %8.4f/alpha  M(rcut) = %8.6f mdisk]\n",
	    getargv0(), rdtab[NTAB-1] * alpha, mdtab[NTAB-1]);
    if ((mdtab[0] == mdtab[1]) || (mdtab[NTAB-2] == mdtab[NTAB-1]))
        error("%s: disk mass table is degenerate\n", getargv0());
    spline(&rdtab[NTAB], &mdtab[0], &rdtab[0], NTAB);	/* for r_d = r_d(m) */
    spline(&vctab[NTAB], &rdtab[0], &vctab[0], NTAB);	/* for v_c = v_c(r) */
}

/*
 * GDISK: compute radial acceleration due to exponential disk.
 */

#define KMAX 1000.0
#define STEP 0.1

real gdisk(real r)
{
    real x = 0.5 * alpha * r;

    if (epsilon >= 0.0)				/* compute smoothed accel.  */
	return (- mdisk * rqbe(alpha) *
		  simpson(dgdisk, r, 0.0, alpha * KMAX, alpha * STEP));
    else					/* use exact expression     */
	return (- mdisk * rqbe(alpha) *
		  0.5 * r * (bessi0(x) * bessk0(x) - bessi1(x) * bessk1(x)));
}

/*
 * DGDISK: integrand for radial acceleration with smoothing.
 * Note: potential at z = 0 in smoothed disk is equal to
 *       potential at z = epsilon in unsmoothed disk.
 */

real dgdisk(real k, real r)
{
    real psq = alpha*alpha + k*k;

    return (rexp(- k * epsilon) * j1(k * r) * k / (psq * rsqrt(psq)));
}

/*
 * SIMPSON: integrate given function using Simpson's rule.
 */

real simpson(real (*integrand)(real, real), real param,
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
	    v1 = v1 + (*integrand)(x, param);
	else if (i % 2 == 0)
	    v4 = v4 + (*integrand)(x, param);
	else
	    v2 = v2 + (*integrand)(x, param);
    }
    return (step1 * (v1 + 4.0*v4 + 2.0*v2) / 3.0);
}

/*
 * MAKEDISK: create realization of disk.
 */

void makedisk(void)
{
    real bfarg, cfact, m, r, phi, vcir, omega, A, kappa, sigma,
	 mu_eff, sig_r, sig_p, sig_z, vrad, vorb2, vorb, vphi;
    double Trr = 0.0, Tpp = 0.0, Tzz = 0.0;
    int i;
    bodyptr dp;

    if (eta > 0) {
	bfarg = 1 / (32 * eta);
	cfact = rsqrt(8 * eta * bessel_k(0.25, bfarg) /
		        (bessel_k(0.75, bfarg) - bessel_k(0.25, bfarg)));
	eprintf("[%s: sigma correction factor = %f]\n", getargv0(), cfact);
    }
    printf("#%5s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s\n",
	   "r", "vcir", "omega", "kappa", "sig_z", "sig_r", "sig_p",
	   "vorb", "Q", "rhomid", "fmax");
    for (i = 0; i < ndisk; i++) {		/* loop initializing bodies */
	m = mdtab[NTAB-1] * ((real) i + 0.5) / ndisk;
	r = seval(m, &mdtab[0], &rdtab[0], &rdtab[NTAB], NTAB);
	vcir = seval(r, &rdtab[0], &vctab[0], &vctab[NTAB], NTAB);
	omega = vcir / r;
	A = (omega - spldif(r, &rdtab[0], &vctab[0], &vctab[NTAB], NTAB)) / 2;
	if (omega - A < 0.0)
	    error("%s: kappa undefined (omega - A < 0)\n"
		  "  r, omega, A = %f %f %f\n", getargv0(), r, omega, A);
	kappa = 2 * rsqrt(rsqr(omega) - A * omega);
	sigma = rsqr(alpha) * mdisk * rexp(- alpha * r) / TWO_PI;
	mu_eff = (r_mu>0 ? 1 + (mu - 1) * (r / (r + r_mu)) : mu);
	sig_z = rsqrt(PI * sigma * z0);
	sig_r = mu_eff * sig_z;
	sig_p = (0.5 * kappa / omega) * sig_r;
	vorb2 = rsqr(vcir) + rsqr(sig_r) * (1 - 2 * alpha * r) - rsqr(sig_p) +
	    (r_mu>0 ? rsqr(sig_z) * r * mu_eff*(2*mu-2)*r_mu/rsqr(r+r_mu) : 0);
	vorb = rsqrt(MAX(vorb2, 0.0));
	dp = NthBody(disk, i);			/* set up ptr to disk body  */
	Mass(dp) = mdisk / ndisk;
	phi = xrandom(0.0, TWO_PI);
	Pos(dp)[0] = r * rsin(phi);
	Pos(dp)[1] = r * rcos(phi);
	Pos(dp)[2] = z0 * ratanh(xrandom(-1.0, 1.0));
	vrad = (eta > 0 ? pickdist(eta, cfact * sig_r) :
                          grandom(0.0, sig_r));
	vphi = (eta > 0 ? pickdist(eta, cfact * sig_p) + vorb :
                          grandom(0.0, sig_p) + vorb);
	Vel(dp)[0] = vrad * rsin(phi) + vphi * rcos(phi);
	Vel(dp)[1] = vrad * rcos(phi) - vphi * rsin(phi);
	Vel(dp)[2] = grandom(0.0, sig_z);
	Trr += Mass(dp) * rsqr(sig_r) / 2;
	Tpp += Mass(dp) * (rsqr(vorb) + rsqr(sig_p)) / 2;
	Tzz += Mass(dp) * rsqr(sig_z) / 2;
	if ((i + 1) % 256 == 0)
	    printf("%6.4f %6.4f %6.2f %6.2f %6.4f %6.4f %6.4f %6.4f "
		   "%6.3f %6.1f %6.1f\n",
		   r, vcir, omega, kappa, sig_z, sig_r, sig_p, vorb,
		   kappa * sig_r / (3.358*sigma), sigma / (2 * z0),
		   sigma / (2 * z0 * rsqrt(rqbe(2 * PI)) * sig_z*sig_r*sig_p));
    }
    eprintf("[%s: Trr = %f  Tpp = %f  Tzz = %f]\n", getargv0(), Trr, Tpp, Tzz);
}

/*
 * BESSEL_K: compute modified Bessel function K of arbitrary order.
 * Uses routine bessik() from "Numerical Recipies in C".
 */

real bessel_k(real ord, real arg)
{
    float bfi, bfk, bfip, bfkp;
    void bessik(float, float, float *, float *, float *, float *);

    bessik((float) arg, (float) ord, &bfi, &bfk, &bfip, &bfkp);
    return (bfk);
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
    int niter;
    real x, y, q;

    niter = 0;
    do {
	x = xrandom(-1.0, 1.0);
	y = xrandom(0.0, YMAX);
	q = rexp(- 0.5 * rsqr(fmap(x)) - eta * rsqr(rsqr(fmap(x)))) *
	      (1 + x*x) / rsqr(1 - x*x);
	if (q > YMAX)
	    error("%s: guess out of bounds\n  x = %f  q = %f > %f\n",
		  getargv0(), x, q, YMAX);
	niter++;
	if (niter > 1000)
	    error("%s: 1000 iterations without success\n", getargv0());
    } while (y > q || x*x == 1);		/* 2nd test prevents infty  */
    return (sigma * fmap(x));
}
