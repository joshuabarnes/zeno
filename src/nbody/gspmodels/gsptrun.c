/*
 * GSPTRUN.C: truncate mass distribution of general spherical profile.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

string defv[] = {		";Truncate mass distribution of GSP",
    "in=???",			";Input file with GSP",
    "out=???",			";Output file for GSP",
    "rtrun=16.0",		";Radius at which truncation starts",
    "rescale=true",		";If true, rescale to original mass",
    "VERSION=1.0",		";Josh Barnes  12 January 1994",
    NULL,
};

void gsptrun(gsprof *, real);
void fixmass(gsprof *, real);

int main(int argc, string argv[])
{
    stream istr, ostr;
    gsprof *gsp;
    real mtot0;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    gsp = get_gsprof(istr);
    mtot0 = gsp->mtot;
    gsptrun(gsp, getdparam("rtrun"));
    if (getbparam("rescale"))
	fixmass(gsp, mtot0);
    eprintf("[%s:  beta = %f  mtot = %f]\n",
	    getargv0(), gsp->beta, gsp->mtot);
    if (gsp->mass[gsp->npoint-2] == gsp->mass[gsp->npoint-1])
	eprintf("[%s: WARNING: M(r) cannot be inverted]\n", getargv0());
    if (gsp->mass[gsp->npoint-1] == gsp->mtot)
	eprintf("[%s: mass has converged by last point]\n", getargv0());
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_gsprof(ostr, gsp);
    strclose(ostr);
    return (0);
}

/*
 * GSPTRUN: impose exponential cutoff at radius rtrun0.
 */

void gsptrun(gsprof *gsp, real rtrun0)
{
    real x, rtrun, mtrun, rstar, kstar, r;
    int i, itrun;

    x = rabs(rtrun0 - gsp->radius[0]);
    for (i = 1; i < gsp->npoint; i++)
	if (rabs(rtrun0 - gsp->radius[i]) < x) {
	    x = rabs(rtrun0 - gsp->radius[i]);
	    itrun = i;
	}
    rtrun = gsp->radius[itrun];
    mtrun = gsp->mass[itrun];
    rstar = - rtrun / (2 + gsp->beta);
    kstar = rsqr(rtrun) * rexp(- (2 + gsp->beta)) * gsp->density[itrun];
    for (i = itrun + 1; i < gsp->npoint; i++) {
	r = gsp->radius[i];
	gsp->density[i] = kstar * rexp(- r / rstar) / rsqr(r);
	gsp->mass[i] = mtrun + 4 * PI * kstar * rstar *
	    	         (rexp(2 + gsp->beta) - rexp(- r / rstar));
    }
    gsp->mtot = mtrun + 4 * PI * kstar * rstar * rexp(2 + gsp->beta);
    gsp->beta = - (2 + gsp->radius[gsp->npoint-1] / rstar);
}

/*
 * FIXMASS: rescale model to total mass mtot0.
 */

void fixmass(gsprof *gsp, real mtot0)
{
    real mscale;
    int i;

    mscale = mtot0 / gsp->mtot;
    for (i = 0; i < gsp->npoint; i++) {
	gsp->density[i] *= mscale;
	gsp->mass[i] *= mscale;
    }
    gsp->mtot = mtot0;
}
