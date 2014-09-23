/*
 * SNAPINTERPOLATE.C: interpolate frames in snapshot stream.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {
    "in=???",
    "out=???",
    "VERSION=1",
    NULL,
};

string bodytags[] = { PosTag, NULL, };

int main(int argc, string argv[])
{
    stream istr, ostr;
    bodyptr btab1 = NULL, btab2 = NULL, btab, btmp;
    int nbody1, nbody2, i;
    real tsnap1, tsnap2, tsnap;
    string itags1[MaxBodyFields], itags2[MaxBodyFields];

    initparam(argv, defv);
    layout_body(bodytags, Precision, NDIM);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    if (! get_snap(istr, &btab1, &nbody1, &tsnap1, itags1, FALSE))
	error("%s: first input failed\n", getargv0());
    if (! set_member(itags1, PosTag))
	error("%s: %s data missing\n", getargv0(), PosTag);
    put_snap(ostr, &btab1, &nbody1, &tsnap1, itags1);
    btab = (bodyptr) allocate(nbody1 * SizeofBody);
    while (get_snap(istr, &btab2, &nbody2, &tsnap2, itags2, FALSE)) {
	if (! set_member(itags2, PosTag))
	    error("%s: %s data missing\n", getargv0(), PosTag);
	if (nbody2 != nbody1)
	    error("%s: nbody changed between frames\n", getargv0());
	tsnap = (tsnap1 + tsnap2) / 2;
	for (i = 0; i < nbody1; i++) {
	    CLRV(Pos(NthBody(btab, i)));
	    ADDMULVS(Pos(NthBody(btab, i)), Pos(NthBody(btab1, i)), 0.5);
	    ADDMULVS(Pos(NthBody(btab, i)), Pos(NthBody(btab2, i)), 0.5);
	}
	put_snap(ostr, &btab, &nbody1, &tsnap, itags2);
	btmp = btab1;
	btab1 = btab2;
	btab2 = btmp;
	tsnap1 = tsnap2;
	put_snap(ostr, &btab1, &nbody1, &tsnap1, itags2);
    }
    return (0);
}
