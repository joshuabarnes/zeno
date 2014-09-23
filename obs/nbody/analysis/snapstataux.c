/*
 * SNAPSTATAUX.C: calculate statistics for Aux variable of input snapshots.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "getparam.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {			";Statistics of Aux data",
    "in=???",				";Input file name",
    "time=all",				";Time to analyze",
    "VERSION=2.0",			";Josh Barnes  14 June 1999",
    NULL,
};

string bodytags[] = { AuxTag, NULL };

void snapstat(bodyptr, int);
int cmpaux(const void *, const void *);

int main(int argc, string argv[])
{
    stream istr;
    bodyptr btab = NULL;
    int nbody;
    real tnow;
    string times, intags[MaxBodyFields];

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    layout_body(bodytags, Precision, NDIM);
    times = getparam("time");
    if (! get_snap_t(istr, &btab, &nbody, &tnow, intags, FALSE, times) ||
	! set_member(intags, AuxTag))
	    error("%s: %s data missing\n", getargv0(), AuxTag);
    snapstat(btab, nbody);
    return (0);
}

void snapstat(bodyptr btab, int nbody)
{
    real avg1, avg2, avg3, avg4;
    bodyptr bp;

    avg1 = avg2 = avg3 = avg4 = 0.0;
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	avg1 += Aux(bp) / nbody;
	avg2 += rsqr(Aux(bp)) / nbody;
	avg3 += rqbe(Aux(bp)) / nbody;
	avg4 += rsqr(rsqr(Aux(bp))) / nbody;
    }
    printf("%12s  %12s  %12s  %12s  %12s\n",
	   "values", "average", "r.m.s.", "avg^3", "avg^4");
    printf("%12d  %12.5f  %12.5f  %12.5f  %12.5f\n",
	   nbody, avg1, rsqrt(avg2), rcbrt(avg3), rsqrt(rsqrt(avg4)));
    if (nbody > 1) {
	qsort(btab, nbody, SizeofBody, cmpaux);
	printf("%12s  %12s  %12s  %12s  %12s\n",
	       "minimum", "1st quart", "median", "3rd quart", "maximum");
	printf("%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n",
	       Aux(NthBody(btab, 0)), 
	       Aux(NthBody(btab, nbody/4)),
	       Aux(NthBody(btab, nbody/2 - 1)),
	       Aux(NthBody(btab, 3*nbody/4 - 1)),
	       Aux(NthBody(btab, nbody - 1)));
    }
}

int cmpaux(const void *a, const void *b)
{
    return (Aux((bodyptr) a) < Aux((bodyptr) b) ? -1 :
	      Aux((bodyptr) a) > Aux((bodyptr) b) ? 1 : 0);
}
