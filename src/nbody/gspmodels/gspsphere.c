/*
 * GSPSPHERE.C: make SPH realization of GSP.
 */

#include "stdinc.h"
#include "assert.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"

string defv[] = {               ";Make SPH gas sphere from GSP",
    "gsp=???",                  ";Input GSP for density profile",
    "out=",                     ";Output SnapShot with bodies",
    "grav=",                    ";Input GSP for gravity calculation",
    "gamma=1.666667",           ";Ratio of specific heats",
    "nbody=4096",               ";Number of bodies to generate",
    "mcut=0.999",               ";Radial cutoff in terms of total mass",
    "seed=54321",               ";Usual random number seed",
    "zerocm=true",              ";Transform to center of mass coords",
    "VERSION=1.0",              ";Josh Barnes  12 May 2012",
    NULL,
};

/* Prototypes for model construction and I/O.                               */

void gspsphere(void);                   /* construct spherical SPH model    */
void readgsp(void);                     /* read input profile(s)	    */
void writemodel(void);                  /* write SPH model to output	    */

/* Global data for communication between major routines.                    */

gsprof *gsp, *ggsp;                     /* profiles for mass and gravity    */
bodyptr btab = NULL;                    /* pointer to array of bodies       */
int  nbody;                             /* number of bodies in array        */

string bodyfields[] = {
    PosTag, VelTag,     MassTag,
    RhoTag, EntFuncTag, UinternTag,
    NULL
};

int main(int argc, string argv[])
{
    initparam(argv, defv);
    readgsp();
    init_random(getiparam("seed"));
    layout_body(bodyfields, Precision, NDIM);
    gspsphere();
    if (getbparam("zerocm"))
        snapcenter(btab, nbody, MassField.offset);
    writemodel();
    return (0);
}

/*
 * GSPSPHERE: construct realization from GSP data.
 */

void gspsphere(void)
{
    real gamma0, mcut, r, sig2, eint = 0.0;
    static real *sig2tab = NULL;
    bodyptr bp;

    nbody = getiparam("nbody");
    assert(nbody > 0);
    gamma0 = getdparam("gamma");
    mcut = getdparam("mcut");
    assert(0.0 < mcut && mcut <= 1.0);
    if (sig2tab == NULL)
        sig2tab = calc_sig2_gsp(gsp, ggsp, 0.0);
    if (btab == NULL)
        btab = (bodyptr) allocate(nbody * SizeofBody);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
        Mass(bp) = gsp->mtot / nbody;
        r = r_mass_gsp(gsp, xrandom(0.0, mcut * gsp->mtot));
	pickshell(Pos(bp), NDIM, r);
	CLRV(Vel(bp));
	Rho(bp) = rho_gsp(gsp, r);
	sig2 = sig2_gsp(gsp, ggsp, 0.0, sig2tab, r);
        EntFunc(bp) = sig2 / rpow(Rho(bp), gamma0 - 1);
	Uintern(bp) = sig2 / (gamma0 - 1);
	eint += Mass(bp) * Uintern(bp);
    }
    eprintf("[%s: thermal energy = %f]\n", getargv0(), eint);
}

void readgsp(void)
{
    stream istr;

    istr = stropen(getparam("gsp"), "r");
    get_history(istr);
    gsp = ggsp = get_gsprof(istr);
    strclose(istr);
    if (! strnull(getparam("grav"))) {
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

    if (! strnull(getparam("out"))) {
        ostr = stropen(getparam("out"), "w");
        put_history(ostr);
        put_snap(ostr, &btab, &nbody, &tsnap, bodyfields);
        strclose(ostr);
    }
}
