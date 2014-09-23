/*
 * FSPTRUN.C: truncate mass distribution of finite spherical profile.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "fsp.h"

string defv[] = {		";Truncate mass distribution of FSP",
    "in=???",			";Input file with FSP",
    "out=???",			";Output file for FSP",
    "rtrun=16.0",		";Radius at which truncation starts",
    "rescale=true",		";If true, rescale to original mass",
    "VERSION=1.0",		";Josh Barnes  12 January 1994",
    NULL,
};

void fsptrun(fsprof *, real);
void fixmass(fsprof *, real);

int main(int argc, string argv[])
{
    stream istr, ostr;
    fsprof *fsp;
    real mtot0;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    fsp = get_fsprof(istr);
    mtot0 = fsp->mtot;
    fsptrun(fsp, getdparam("rtrun"));
    if (getbparam("rescale"))
	fixmass(fsp, mtot0);
    eprintf("[%s:  beta = %f  mtot = %f]\n",
	    getargv0(), fsp->beta, fsp->mtot);
    if (fsp->mass[fsp->npoint-2] == fsp->mass[fsp->npoint-1])
	eprintf("[%s: WARNING: M(r) cannot be inverted]\n", getargv0());
    if (fsp->mass[fsp->npoint-1] == fsp->mtot)
	eprintf("[%s: mass has converged by last point]\n", getargv0());
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_fsprof(ostr, fsp);
    strclose(ostr);
    return (0);
}

/*
 * FSPTRUN: impose exponential cutoff at radius rtrun0.
 */

void fsptrun(fsprof *fsp, real rtrun0)
{
    real x, rtrun, mtrun, rstar, kstar, r;
    int i, itrun;

    x = rabs(rtrun0 - fsp->radius[0]);
    for (i = 1; i < fsp->npoint; i++)
	if (rabs(rtrun0 - fsp->radius[i]) < x) {
	    x = rabs(rtrun0 - fsp->radius[i]);
	    itrun = i;
	}
    rtrun = fsp->radius[itrun];
    mtrun = fsp->mass[itrun];
    rstar = - rtrun / (2 + fsp->beta);
    kstar = rsqr(rtrun) * rexp(- (2 + fsp->beta)) * fsp->density[itrun];
    for (i = itrun + 1; i < fsp->npoint; i++) {
	r = fsp->radius[i];
	fsp->density[i] = kstar * rexp(- r / rstar) / rsqr(r);
	fsp->mass[i] = mtrun + 4 * PI * kstar * rstar *
	    	         (rexp(2 + fsp->beta) - rexp(- r / rstar));
    }
    fsp->mtot = mtrun + 4 * PI * kstar * rstar * rexp(2 + fsp->beta);
    fsp->beta = - (2 + fsp->radius[fsp->npoint-1] / rstar);
}

/*
 * FIXMASS: rescale model to total mass mtot0.
 */

void fixmass(fsprof *fsp, real mtot0)
{
    real mscale;
    int i;

    mscale = mtot0 / fsp->mtot;
    for (i = 0; i < fsp->npoint; i++) {
	fsp->density[i] *= mscale;
	fsp->mass[i] *= mscale;
    }
    fsp->mtot = mtot0;
}
