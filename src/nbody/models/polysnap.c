/*
 * POLYSNAP: Construct generalized polytrope.
 */
 
#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include <assert.h>
 
string defv[] = {		";Construct generalized polytrope",
    "out=???",			";Snapshot output file",
    "n=1.5",			";Binding energy index; must be > 1/2",
    "m=1.5",			";Angular momentum index; must be > -1",
    "nbody=4096",		";Number of bodies to generate",
    "seed=54321",		";Random number seed",
    "nmodel=1",			";Number of copies to generate",
    "hstep=0.01",		";Initial integration step",
    "besort=true",		";Sort particles by binding energy",
    "zerocm=true",		";Transform to center of mass coords",
    "VERSION=1.4",		";Josh Barnes	12 May 2012",
    NULL,
};

/* Procedure prototypes. */

void polymodel1(void);
void polymodel2(void);
void polysolve(real);
void initstep(real [], real);
void diffrpmw(real [], real []);
real rho(real, real);
void storestep(real []);
void fixsurface(void);
void polyscale(void);
void polymodel(void);
real rad_m(real);
real pick_psi(void);
real hfunc(real);
real pick_v(real);
real gfunc1(real);
real gfunc2(real);
void diffstep(real [], int, void (*)(real [], real []), real);
real vnpick(real (*)(real), real, real, real, string);
int berank(const void *, const void *);

extern double lgamma(double);

/* Global data. */

string bodyfields[] = { PosTag, VelTag, MassTag, PhiTag, AuxTag, NULL };

bodyptr btab;			/* body array generated below */
int nbody;			/* number of bodies in array */

real n, m;			/* polytropic indicies */
real Kprime;			/* constant in f(E, J) */
real rad1, phi1, mtot;		/* surface radius, potential, mass */

#define K	1.0		/* constant in rho(r, phi) */
#define PHI0	-1.0		/* phi(r=0) before rescaling */

#define MAXSTEP 1024		/* maximum integration steps */
real rtab[MAXSTEP];		/* radial integration points */
real ptab[MAXSTEP];		/* potential as a fn of radius */
real mtab[MAXSTEP];		/* enclosed mass as a fn of radius */
real pcoef[3*MAXSTEP];		/* coefficients for phi = phi(r) */
int nstep;			/* number of integration steps */

int main(int argc, string argv[])
{
    stream outstr;
    int nmodel;
    real tzero = 0.0;

    initparam(argv, defv);
    layout_body(bodyfields, Precision, NDIM);
    n = getdparam("n");
    m = getdparam("m");
    nbody = getiparam("nbody");
    btab = (bodyptr) allocate(nbody * SizeofBody);
    init_random(getiparam("seed"));
    nmodel = getiparam("nmodel");
    if ((! (n == 1.0 && m == -1.0)) && (! (n == 0.5 && m == -0.5)))
        polysolve(getdparam("hstep"));
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    while (--nmodel >= 0) {
	if (n == 1.0 && m == -1.0)
	    polymodel1();
	else if (n == 0.5 && m == -0.5)
	    polymodel2();
	else
	    polymodel();
	if (getbparam("besort"))
	    qsort(btab, nbody, SizeofBody, berank);
	if (getbparam("zerocm"))
            snapcenter(btab, nbody, MassField.offset);
	put_snap(outstr, &btab, &nbody, &tzero, bodyfields);
    }
    strclose(outstr);
    return (0);
}

/*
 * POLYMODEL1: generate polytrope model for n = 1, m = -1.
 */

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
    bodyfields[4] = NULL;			/* don't output Aux data    */
}

/*
 * POLYMODEL2: generate polytrope model for n = 1/2, m = -1/2.
 */

void polymodel2(void)
{
    error("%s: special case n = 1/2, m = -1/2 not implemented\n", getargv0());
    bodyfields[4] = NULL;			/* don't output Aux data    */
}

/*
 * POLYSOLVE: solve structure equation for generalized polytrope.
 */
 
void polysolve(real hstep)
{
    real h, rpmw[4], w1, w2;
 
    if (n <= 0.5 || m <= -1.0)
	error("%s: illegal value for n or m\n", getargv0());
    if (n >= 3*m + 5)
	error("%s: model would have infinite radius\n", getargv0());
    Kprime = K * (m + 1) * rpow(PI, -1.5) * rpow(2.0, - (m + 1.5)) *
	       rexp(lgamma(m + n + 1) - lgamma(m + 2) - lgamma(n - 0.5));
    h = hstep;
    do {
	rpmw[0] = 0.0;				/* init. radius             */
	rpmw[1] = PHI0;				/* init. relative potential */
	rpmw[2] = 0.0;				/* init. enclosed mass	    */
	rpmw[3] = 0.0;				/* init. binding energy     */
	nstep = 0;
	storestep(rpmw);
	while (rpmw[1] < 0.0) {			/* while rel. pot. is neg.  */
	    if (rpmw[0] == 0.0)			/* starting first step?     */
	        initstep(rpmw, h);		/* use asymp. approximation */
	    else
		diffstep(rpmw, 4, diffrpmw, h);	/* use num. integration     */
	    storestep(rpmw);
	}
	fixsurface();
	w1 = rpmw[3] - 0.5 * rsqr(mtot) / rad1;
	w2 = - (2*m + 3) / (3*m - n + 5) * rsqr(mtot) / rad1;
	eprintf("[polysolve: nstep = %3d  W = %f  error = %f]\n",
		nstep, w2, (w1 - w2)/w2);
	h = 0.5 * h;				/* refine step-size by 2    */
    } while (nstep < MAXSTEP/4);
    polyscale();
}
 
/*
 * INITSTEP: asymptotic approx. for generalized polytrope near r=0.
 */

void initstep(real rpmw[], real h)
{
    rpmw[0] = h;
    rpmw[1] = (FOUR_PI * K / ((2*m + 2) * (2*m + 3))) *
                rpow(- PHI0, n + m) * rpow(h, 2*m + 2) + PHI0;
    rpmw[2] = (FOUR_PI * K / (2*m + 3)) *
                rpow(- PHI0, n + m) * rpow(h, 2*m + 3);
    rpmw[3] = 0.5 * PHI0 * rpmw[2];
}

/*
 * DIFFRPMW: differentials of radius, potential, mass, binding energy.
 */

void diffrpmw(real drpmw[], real rpmw[])
{
    drpmw[0] = 1.0;
    drpmw[1] = rpmw[2] / rsqr(rpmw[0]);
    drpmw[2] = FOUR_PI * rsqr(rpmw[0]) * rho(rpmw[0], rpmw[1]);
    drpmw[3] = 0.5 * rpmw[1] * drpmw[2];
}

real rho(real r, real phi)
{
    return (phi < 0.0 ? K * rpow(- phi, n + m) * rpow(r, 2 * m) : 0.0);
}

/*
 * STORESTEP: tabulate results of integration.
 */

void storestep(real rpmw[])
{
    if (nstep >= MAXSTEP)
	error("%s.storestep: too many steps\n", getargv0());
    rtab[nstep] = rpmw[0];
    ptab[nstep] = rpmw[1];
    mtab[nstep] = rpmw[2];
    nstep++;
}

/*
 * FIXSURFACE: locate radius of surface by linear interpolation.
 */

void fixsurface(void)
{
    real f;

    f = -ptab[nstep-2] / (ptab[nstep-1] - ptab[nstep-2]);
    rad1 = f * rtab[nstep-1] + (1 - f) * rtab[nstep-2];
    mtot = f * mtab[nstep-1] + (1 - f) * mtab[nstep-2];
    phi1 = 0.0;					/* fixed by definition      */
}

/*
 * POLYSCALE: scale model to Henon's units.
 */

void polyscale(void)
{
    real r_henon, m_henon, p_henon, scl_r, scl_m, scl_v;
    int i;

    r_henon = (4*m + 6) / (3*m - n + 5);
    m_henon = 1.0;
    p_henon = m_henon / r_henon;
    scl_r = r_henon / rad1;
    scl_m = m_henon / mtot;
    scl_v = rsqrt(scl_m / scl_r);
    for (i = 0; i < nstep; i++) {
	rtab[i] = scl_r * rtab[i];
	mtab[i] = scl_m * mtab[i];
	ptab[i] = rsqr(scl_v) * ptab[i] - p_henon;
    }
    spline(pcoef, rtab, ptab, nstep);
    rad1 = scl_r * rad1;
    mtot = scl_m * mtot;
    phi1 = rsqr(scl_v) * phi1 - p_henon;	/* == - p_henon by def.     */
    Kprime = Kprime * scl_m / rqbe(scl_r * scl_v);
    eprintf("[polyscale: Kprime = %f]\n", Kprime);
}

/*
 * POLYMODEL: construct N-body realization of polytrope.
 */

void polymodel(void)
{
    bodyptr p;
    real x, r, pot, v, psi, vr, vp, a, E, J;
    vector rhat, vtmp, vper;

    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p)) {
	x = xrandom(0.0, mtot);
	r = rad_m(x);
	pot = seval(r, rtab, ptab, pcoef, nstep);
	v = pick_v(pot);
	psi = pick_psi();
	vr = v * rcos(psi);
	vp = v * rsin(psi);
	Mass(p) = mtot / nbody;
	pickshell(rhat, NDIM, 1.0);
	MULVS(Pos(p), rhat, r);
	pickshell(vtmp, NDIM, 1.0);
	a = dotvp(vtmp, rhat);
	MULVS(vper, rhat, -a);
	ADDV(vper, vper, vtmp);
	a = absv(vper);
	MULVS(vper, vper, vp / a);
	MULVS(Vel(p), rhat, vr);
	ADDV(Vel(p), Vel(p), vper);
	Phi(p) = pot;
	E = pot + 0.5 * rsqr(v);
	J = r * ABS(vp);
	Aux(p) = Kprime * rpow(phi1 - E, n - 1.5) * rpow(J, 2 * m);
    }
}

/*
 * RAD_M: compute radius as a function of enclosed mass.  Cannot use
 * a spline since mass increments can become very small near r_1.
 */

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
	      getargv0(), x, i, mtab[i], k, mtab[k]);
    f = (x - mtab[i]) / (mtab[k] - mtab[i]);
    return (f * rtab[k] + (1 - f) * rtab[i]);
}

/*
 * PICK_PSI: pick angle between velocity and radial vectors.
 */

real pick_psi(void)
{
    static real hmax = -1.0;
    real x, psi;

    if (hmax < 0.0)
	hmax = (m > 0 ? rpow(2.0, m) : 1);
    x = vnpick(hfunc, 0.0, 1.0, hmax, "hfunc");
    psi = racos(1 - rpow(x, 1 / (m + 1)));
    return (xrandom(-1.0, 1.0) < 0.0 ? PI - psi : psi);
}

real hfunc(real x)
{
    return (rpow(2 - rpow(x, 1 / (m + 1)), m));
}

/*
 * PICK_V: pick magnitude of velocity.
 */

real pick_v(real pot)
{
    real vmax, x;
    static real gmax = -1.0;

    vmax = rsqrt(2 * phi1 - 2 * pot);
    if (n > 1.5) {
	x = vnpick(gfunc1, 0.0, 1.0, 1.0, "gfunc1");
	return (vmax * x);
    } else {
	if (gmax < 0.0) {
	    x = (n + 4*m + 2.5) / (n + 2*m + 0.5);
	    if (0.0 <= x && x < 1.0)
		gmax = rpow(2 - x, n - 1.5) * rpow(1 - x, 2*m + 2);
	    else
		gmax = MAX(gfunc2(0.0), gfunc2(0.99));
	}
	x = vnpick(gfunc2, 0.0, 1.0, gmax, "gfunc2");
	return (vmax * (1 - rpow(x, 1 / (n - 0.5))));
    }
}

real gfunc1(real x)
{
    if (rlog10(x) * (2*m + 2) < -36.0)
	return (0.0);
    else
	return (rpow(1 - rsqr(x), n - 1.5) * rpow(x, 2*m + 2));
}

real gfunc2(real x)
{
    real y;

    y = rpow(x, 1 / (n - 0.5));
    return (rpow(2 - y, n - 1.5) * rpow(1 - y, 2*m + 2));
}

/*
 * DIFFSTEP: solve coupled difeqs using Runge-Kutta method.
 */
 
#define MAXVARS  8
 
void diffstep(real x0[], int nx, void (*diff)(real [], real []), real delta)
{
    real dxi, x1[MAXVARS], dx0[MAXVARS], dx1[MAXVARS];
    int i;
 
    (*diff)(dx1, x0);
    for (i = 0; i < nx; i++) {
        dxi = 0.5 * delta * dx1[i];
        dx0[i] = dxi;
        x1[i] = x0[i] + dxi;
    }
    (*diff)(dx1, x1);
    for (i = 0; i < nx; i++) {
        dxi = delta * dx1[i];
        dx0[i] = dx0[i] + dxi;
        x1[i] = x0[i] + 0.5 * dxi;
    }
    (*diff)(dx1, x1);
    for (i = 0; i < nx; i++) {
        dxi = delta * dx1[i];
        dx0[i] = dx0[i] + dxi;
        x1[i] = x0[i] + dxi;
    }
    (*diff)(dx1, x1);
    for (i = 0; i < nx; i++)
        x0[i] = x0[i] + (dx0[i] + 0.5 * delta * dx1[i]) / 3;
}

/*
 * VNPICK: chose from a distribution by von Neumann technique.
 * Ref: von Neumann, J. 1963. Collected Works, Vol. 5.
 */

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
	    eprintf("[vnpick: %s(%g) = %g > %s_max = %g]\n",
		    name, x, fx, name, fmax);
	    warn = TRUE;
	}
	if (fx > 1.1 * fmax)
	    error("%s.vnpick: %s(x) out of bounds\n", getargv0(), name);
	ncyc++;
	if (ncyc % NCYC == 0)
	    eprintf("[vnpick: %d cycles picking %s(x)]\n", ncyc, name);
    }
    return (x);
}

int berank(const void *p1, const void *p2)
{
    real e1, e2;

    e1 = 0.5 * rsqr(absv(Vel((bodyptr) p1))) + Phi((bodyptr) p1);
    e2 = 0.5 * rsqr(absv(Vel((bodyptr) p2))) + Phi((bodyptr) p2);
    return (e1 < e2 ? -1 : e1 > e2 ? 1 : 0);
}
