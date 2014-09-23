/*
 * GSPREALIZE.C: make N-body realization of GSP.
 *
 * Ths code has some limitations which need to be addressed: (1) Computation
 * of hmax(phi) is flawed; a fudge factor should NOT be becessary! (2) The
 * innermost radius of the density-profile gsp appears to set a lower limit
 * on the particle radii.  (3) Setting randrad=true appears to create a
 * slightly biased radial distribution.
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

string defv[] = {		";Construct N-body realization of GSP.",
				";Uses Abel integral to compute DF.",
    "gsp=???",			";Input GSP for density profile",
    "out=",			";Output SnapShot with bodies",
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
    "hfudge=1.5",		";Fudge factor for hmax(phi).",
				";This should NOT be necessary!",
    "VERSION=1.2",		";Josh Barnes  17 June 2011",
    NULL,
};

/* Prototypes for model construction.					    */

local void dfcompute(void);			/* tabulate dist. func.     */
local void intfunc(real *, real *);		/* derivatives for g(E)     */
local void hmaxinit(void);			/* init table for hmax(phi) */
local float hneg(float);			/* convenience for golden() */
local void gsprealize(void);			/* construct realization    */
local int berank(const void *, const void *);	/* compare binding energies */
local real pickspeed(real);			/* pick speed from h(v)     */
local real h_v(real, real);			/* speed distribution func  */
local real f_E(real);				/* phase-space dist. func   */

real diffstep(real *, real *, int, void (*)(real *, real *), real);
float golden(float, float, float, float (*)(float), float, float *);

/* Global data for communication between routines.			    */

local gsprof *gsp, *ggsp;		/* profiles for mass and gravity    */
local real *Etab, *Ftab, *Fcoef;	/* tables for f(E) = dF/dE          */
local int ntab;				/* number of values in tables	    */
local real *htab, *hcoef;		/* tables for hmax(phi)		    */
local bodyptr btab = NULL;		/* pointer to array of bodies	    */
local int nbody;			/* number of bodies in array	    */
local real phi0;			/* potential for hneg() function    */
local real hfudge;			/* fudge factor for hmax(phi)	    */

/* Miscellaneous parameters and constants.				    */

#define GOLDTOL 0.0001			/* relative tol. for golden()	    */
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
    dfcompute();
    if (! strnull(getparam("out"))) {
	layout_body(bodyfields, Precision, NDIM);
        nbody = getiparam("nbody");
	assert(nbody > 0);
	srandom(getiparam("seed"));
	hmaxinit();
	hfudge = getdparam("hfudge");
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
 * DFCOMPUTE: numerically integrate F(E) and tabulate result.
 */

local void dfcompute(void)
{
    real xmax, xFE[3], err, avgerr, maxerr;
    int k, nstep, i, j;

    (void) phi_gsp(ggsp, 1.0);			/* compute potential table  */
    for (k = 0; ggsp->phi[k] == ggsp->phi[k+1]; k++);
    if (k > 0)
        eprintf("[%s.dfcompute: skipping first %d values of phi]\n",
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
	xFE[2] = Etab[i];			/* pass E on to intfunc()   */
	for (j = 0; j < nstep; j++) {
	    err = diffstep(xFE, xFE, 3, intfunc, xmax / nstep);
	    avgerr += err / (ntab * nstep);
	    maxerr = MAX(maxerr, err);
	}
	Ftab[i] = xFE[1];
    }
    eprintf("[%s.dfcompute: avgerr = %f  maxerr = %f]\n",
	    getargv0(), avgerr, maxerr);
    spline(Fcoef, Etab, Ftab, ntab);
    if (getbparam("dflist")) {
	printf("#%11s  %12s  %12s\n", "E", "F(E)", "f(E)");
	for (i = 0; i < ntab; i++)
	    printf("%12.6f  %12.6g  %12.6g\n", Etab[i],
		   seval(Etab[i], Etab, Ftab, Fcoef, ntab),
		   spldif(Etab[i], Etab, Ftab, Fcoef, ntab));
    }
}    

/*
 * INTFUNC: evaluate dF/dx for dfcompute().
 */

local void intfunc(real *dxFE, real *xFE)
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
 * HMAXINIT: initalize the table used to find hmax(phi).
 */

local void hmaxinit(void)
{
    int i;
    float fmax, vtop, va, vb, vc, vmax;

    htab = (real *) allocate(ntab * sizeof(real));
    hcoef = (real *) allocate(3 * ntab * sizeof(real));
    for (i = 0; i < ntab-1; i++) {
        vtop = rsqrt(-2.0 * (Etab[i] - Etab[ntab-1]));
        va = 0.0;
	vb = (i == 0 ? 0.5 : fmax) * vtop;
	vc = vtop;
	phi0 = Etab[i];
	if (hneg(vb) < hneg(vc)) {		/* got vmax bracketed?      */
	    (void) golden(va, vb, vc, hneg, GOLDTOL, &vmax);
	    htab[i] = h_v((real) vmax, Etab[i]);
	    fmax = vmax / vc;
	} else {				/* just use top velocity    */
	    htab[i] = h_v((real) vtop, Etab[i]);
	}
    }
    htab[ntab-1] = 0.0;
    spline(hcoef, Etab, htab, ntab);
}

/*
 * HNEG: return - h(v;phi0) for minimum finder golden().
 */

local float hneg(float v)
{
    return (- h_v((real) v, phi0));
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
#ifdef BETTERSAMP
	  x = (i%2 == 0 ? xrandom(mmin, mmax) : xrandom(mmax, mmin));
	/* 
         * NOTE: swapping low and high limits on xrandom appears to reduce
	 * (but not eliminate) the bias in the radial distribution.
	 */
#else
	  x = xrandom(mmin, mmax);
#endif
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
 * BERANK: rank bodies by binding energy.
 */

local int berank(const void *a, const void *b)
{
    real Ea, Eb;

    Ea = 0.5 * dotvp(Vel((bodyptr) a), Vel((bodyptr) a)) + Phi((bodyptr) a);
    Eb = 0.5 * dotvp(Vel((bodyptr) b), Vel((bodyptr) b)) + Phi((bodyptr) b);
    return (Ea < Eb ? -1 : Ea > Eb ? 1 : 0);
}

/*
 * PICKSPEED: chose speed distributed according to h(v).
 */

local real pickspeed(real phi)
{
    real vmax, hmax, v0, h0, hv;
    int n, nwarn;

    vmax = rsqrt(-2.0 * (phi - Etab[ntab-1]));
    hmax = hfudge * seval(phi, Etab, htab, hcoef, ntab);
    nwarn = NTRIAL;
    n = 0;
    do {
	if (n > nwarn) {
	    eprintf("[%s.pickspeed: %d iterations of pickspeed]\n",
		    getargv0(), n);
	    nwarn += NTRIAL;
	}
	v0 = xrandom(0.0, vmax);
	h0 = xrandom(0.0, hmax);
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
