/*
 * SNAPRHOLIST: list spherically-averaged density profile.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "getparam.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";List spherical density profile",
    "in=???",			";Input file with N-body data",
    "times=all",		";Select frames averaged together",
    "mass=m",			";Expr giving particle mass",
    "nprof=256",		";Number of fine-grain radial bins",
    "rrange=0.001:10",		";Range of radii to tabulate",
    "fcrit=0.01",		";Minimum mass measurement",
    "require=",			";List of required items",
    "seed=",			";Seed for random number generator",
    "VERSION=1.1",		";Josh Barnes  14 May 2008",
    NULL,
};

string fields[] = { PosTag, MassTag, NULL };

stream execmap(string);
void setprofile(real *, real *, int, real [], bodyptr, int);
void listdensity(real *, real *, int, real [], int);

int main(int argc, string argv[])
{
    string prog, itags[MaxBodyFields];
    stream xstr;
    int nprof, nbody, nsamp;
    real *prof1, *prof2, rrange[2], tnow;
    bodyptr btab = NULL;

    initparam(argv, defv);
    layout_body(fields, Precision, NDIM);
    prog = tempnam("/tmp", "sm");
    xstr = execmap(prog);
    get_history(xstr);
    nprof = getiparam("nprof");
    prof1 = (real *) allocate((2 + nprof) * sizeof(real));
    prof2 = (real *) allocate((2 + nprof) * sizeof(real));
    setrange(rrange, getparam("rrange"));
    nsamp = 0;
    while (get_snap(xstr, &btab, &nbody, &tnow, itags, FALSE)) {
	setprofile(prof1, prof2, nprof, rrange, btab, nbody);
	nsamp++;
    }
    if (unlink(prog) != 0)
        error("%s: can't unlink %s\n", getargv0(), prog);
    if (nsamp == 0)
	error("%s: no data in input\n", getargv0());
    listdensity(prof1, prof2, nprof, rrange, nsamp);
    return (0);
}

#include <unistd.h>

void buildmap(string, string *, string *, string *, string, int);

stream execmap(string prog)
{
    int handle[2];
    char handbuf[32];
    string names[] = { MassTag, NULL }, exprs[] = { NULL, NULL };

    exprs[0] = getparam("mass");
    buildmap(prog, names, exprs, NULL, Precision, NDIM);
    pipe(handle);
    if (fork() == 0) {                           /* if this is child process */
        close(handle[0]);
        sprintf(handbuf, "-%d", handle[1]);
        execl(prog, getargv0(), getparam("in"), handbuf, getparam("times"),
              getparam("require"), MassTag "," PosTag,
	      strnull(getparam("require")) ? "true" : "false",
	      getparam("seed"), NULL);
        error("%s: execl %s failed\n", getargv0(), prog);
    }
    close(handle[1]);
    sprintf(handbuf, "-%d", handle[0]);
    return (stropen(handbuf, "r"));
}

void setprofile(real *prof1, real *prof2, int nprof, real rrange[],
		bodyptr btab, int nbody)
{
    real logrmin, logrdif, logr;
    bodyptr bp;
    int j;

    logrmin = rlog10(rrange[0]);
    logrdif = rlog10(rrange[1] / rrange[0]);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	logr = rlog10(absv(Pos(bp)));
	j = 1 + floor(nprof * (logr - logrmin) / logrdif);
	j = MAX(j, 0);
	j = MIN(j, nprof + 1);
	prof1[j] += Mass(bp);
	prof2[j] += rsqr(Mass(bp));
    }
}

void listdensity(real *prof1, real *prof2, int nprof, real rrange[], int nsamp)
{
    real logrmin, logrdif, mtot, mcrit, msum, rmid, rmsum, rad1, rad2, vol12;
    real esum;
    int i;

    logrmin = rlog10(rrange[0]);
    logrdif = rlog10(rrange[1] / rrange[0]);
    mtot = 0.0;
    for (i = 0; i <= nprof + 1; i++)
	mtot += prof1[i];
    mcrit = mtot * getdparam("fcrit");
    msum = prof1[0];
    esum = prof2[0];
    rmid = 0.707 * rrange[0];
    rmsum = rmid * prof1[0];
    rad1 = 0.0;
    printf("#%11s%12s%12s%12s%12s\n",
	   "rmin", "ravg", "rmax", "rho", "error");
    for (i = 1; i < nprof + 1; i++) {
	msum += prof1[i];
	esum += prof2[i];
	rmid = rdex(logrmin + (i - 0.5) * logrdif / nprof);
	rmsum += rmid * prof1[i];
	if (msum > mcrit) {
	    rad2 = rdex(logrmin + ((real) i) * logrdif / nprof);
	    vol12 = FRTHRD_PI * (rqbe(rad2) - rqbe(rad1)) * nsamp;
	    printf("%#12.5g%#12.5g%#12.5g%#12.5g%#12.5g\n",
		   rad1, rmsum / msum, rad2, msum / vol12,
		   rsqrt(esum) / vol12);
	    msum = 0.0;
	    esum = 0.0;
	    rmsum = 0.0;
	    rad1 = rad2;
	}
    }
}
