/*
 * GSPMODEL.C: make N-body model of GSP.
 */

#include "stdinc.h"
#include "assert.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"

string defv[] = {		";Make N-body model of GSP",
    "gsp=???",			";Input GSP for density profile",
    "out=",			";Output SnapShot with bodies",
    "grav=",		        ";Input GSP for gravity calculation",
    "beta_a=0.0",		";Anisotropy parameter: beta_a <= 1",
    "nbody=4096",		";Number of bodies to generate",
    "mcut=0.999",		";Radial cutoff in terms of total mass",
    "vcut=1.0",			";Velocity cut in terms of escape vel",
    "seed=54321",		";Usual random number seed",
    "besort=true",		";Sort particles by binding energy",
    "zerocm=true",		";Transform to center of mass coords",
    "auxvar=",			";Options are mass, rperi",
    "VERSION=2.4",		";Josh Barnes  12 May 2012",
    NULL,
};

/* Prototypes for model construction and I/O.				    */

void gspmodel(void);			/* routine to construct model	    */
real fixsigma(real, real);		/* correct sigma for vel cutoff	    */
real sigeff(real, real);		/* compute <v^2>^0.5 with cutoff    */
void picktriad(vector, vector, vector);	/* pick three orthogonal vectors    */
int berank(const void *, const void *);	/* compare binding energies         */
void setauxvar(bodyptr, int);		/* set aux to mass(r), rperi, ...   */
void readgsp(void);			/* routine to read input profiles   */
void writemodel(void);			/* routine to write output model    */

/* Global data for communication between major routines.		    */

gsprof *gsp, *ggsp;			/* profiles for mass and gravity    */
bodyptr btab = NULL;			/* pointer to array of bodies	    */
int  nbody;				/* number of bodies in array	    */

string bodyfields[] = { PosTag, VelTag, MassTag, PhiTag, AuxTag, NULL };

int main(int argc, string argv[])
{
    initparam(argv, defv);
    readgsp();
    init_random(getiparam("seed"));
    layout_body(bodyfields, Precision, NDIM);
    gspmodel();
    writemodel();
    return (0);
}

/*
 * GSPMODEL: construct realization from GSP data.
 */

void gspmodel(void)
{
    real beta_a, mcut, vcut, vfac;
    static real *sig2 = NULL;
    real r, vmax2, sigr, sig, x, y, vr, v1, v2;
    bodyptr bp;
    vector rhat, vec1, vec2, vtmp;

    beta_a = getdparam("beta_a");
    assert(beta_a <= 1.0);
    nbody = getiparam("nbody");
    assert(nbody > 0);
    mcut = getdparam("mcut");
    assert(0.0 < mcut && mcut <= 1.0);
    vcut = getdparam("vcut");
    assert(vcut > 0.0);
    if (sig2 == NULL)
        sig2 = calc_sig2_gsp(gsp, ggsp, beta_a);
    if (btab == NULL)
	btab = (bodyptr) allocate(nbody * SizeofBody);
    vfac = rsqrt(1 - beta_a);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	Mass(bp) = gsp->mtot / nbody;
	r = r_mass_gsp(gsp, xrandom(0.0, mcut * gsp->mtot));
	vmax2 = -2 * rsqr(vcut) * phi_gsp(ggsp, r);
	if (vfac > 1.0)
	    vmax2 = vmax2 / rsqr(vfac);
	sigr = rsqrt(sig2_gsp(gsp, ggsp, beta_a, sig2, r));
	sig = fixsigma(sigr, rsqrt(vmax2));
	do {
	    vr = grandom(0.0, sig);
	    v1 = grandom(0.0, sig);
	    v2 = grandom(0.0, sig);
	} while (vr*vr + v1*v1 + v2*v2 > vmax2);
	picktriad(rhat, vec1, vec2);
	MULVS(Pos(bp), rhat, r);
	MULVS(Vel(bp), rhat, vr);
	MULVS(vtmp, vec1, v1 * vfac);
	ADDV(Vel(bp), Vel(bp), vtmp);
	MULVS(vtmp, vec2, v2 * vfac);
	ADDV(Vel(bp), Vel(bp), vtmp);
	Phi(bp) = phi_gsp(ggsp, r);
	Aux(bp) = Phi(bp) + 0.5 * dotvp(Vel(bp), Vel(bp));
    }
    if (getbparam("besort"))
	qsort(btab, nbody, SizeofBody, berank);
    if (getbparam("zerocm"))
	snapcenter(btab, nbody, MassField.offset);
    if (! strnull(getparam("auxvar")))
	setauxvar(btab, nbody);
}

real fixsigma(real sig0, real vmax)
{
    real sig1, sig2, sig, sige;

    sig1 = sig0;
    sig2 = 2.0 * sig1;
    while (sig2 < 8.0 * vmax && sigeff(sig2, vmax) < sig0)
	sig2 = 2 * sig2;
    if (sigeff(sig2, vmax) < sig0)
	eprintf("[fixsigma: can\'t attain required sigma]\n");
    while ((sig2 - sig1) / (sig2 + sig1) > 0.00001) {
	sig = (sig2 + sig1) / 2.0;
	sige = sigeff(sig, vmax);
	if (sige > sig0)
	    sig2 = sig;
	else
	    sig1 = sig;
    }
    return (sig);
}

#define SQRT2   1.41421356237
#define SQRTPI  1.77245385091

real sigeff(real sig, real vmax)
{
    real x, y, z, Q;

    x = vmax / sig;
    y = (SQRTPI / SQRT2) * erf(x / SQRT2);
    z = rexp(rsqr(x) / 2);
    Q = ((-3 * x - rqbe(x)) / z + 3 * y) / (- x / z + y);
    return (sig * rsqrt(Q / 3));
}

void picktriad(vector x, vector y, vector z)
{
    real a;

    pickshell(x, NDIM, 1.0);
    pickshell(z, NDIM, 1.0);
    CROSSVP(y, x, z);
    a = absv(y);
    DIVVS(y, y, a);
    CROSSVP(z, x, y);
}

int berank(const void *a, const void *b)
{
    return (Aux((bodyptr) a) < Aux((bodyptr) b) ? -1 :
	      Aux((bodyptr) a) > Aux((bodyptr) b) ? 1 : 0);
}

#define TOL  0.0001

void setauxvar(bodyptr btab, int nbody)
{
    bodyptr bp;
    vector jvec;
    real jtot, etot, r0, r1, r;

    if (streq(getparam("auxvar"), "mass"))
	for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	    Aux(bp) = mass_gsp(ggsp, absv(Pos(bp)));
    else if (streq(getparam("auxvar"), "rperi"))
	for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	    CROSSVP(jvec, Pos(bp), Vel(bp));
	    jtot = absv(jvec);
	    etot = 0.5 * dotvp(Vel(bp), Vel(bp)) + Phi(bp);
	    r0 = 0.0;
	    r1 = absv(Pos(bp));
	    r = 0.5 * (r0 + r1);
	    while ((r1 - r0) > TOL * r) {
		if (rsqrt(2 * (etot - phi_gsp(ggsp, r))) > jtot/r)
		    r1 = r;
		else
		    r0 = r;
		r = 0.5 * (r0 + r1);
	    }
	    Aux(bp) = r;
	}
    else
	error("%s: unknown auxvar option %s\n",
	      getargv0(), getparam("auxvar"));
}

void readgsp(void)
{
    stream istr;

    istr = stropen(getparam("gsp"), "r");
    get_history(istr);
    gsp = get_gsprof(istr);
    strclose(istr);
    if (strnull(getparam("grav")))
	ggsp = gsp;
    else {
	istr = stropen(getparam("grav"), "r");
	get_history(istr);
	ggsp = get_gsprof(istr);
	strclose(istr);
    }
}

void writemodel(void)
{
    stream ostr;
    real tsnap = 0.0;
    string outfields[] = { PosTag, VelTag, MassTag, PhiTag, AuxTag, NULL };

    if (! strnull(getparam("out"))) {
	if (strnull(getparam("auxvar")))
	    outfields[3] = NULL;
	ostr = stropen(getparam("out"), "w");
	put_history(ostr);
	put_snap(ostr, &btab, &nbody, &tsnap, outfields);
	strclose(ostr);
    }
}
