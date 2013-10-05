/*
 * SNAPFORCE.C: compute forces by direct-sum algorithm.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Compute forces by direct-sum",
    "in=???",			";N-body snapshot input file",
    "out=???",			";N-body snapshot output file",
    "nforce=",			";Number of bodies with forces",
    "eps=0.01",			";Plummer smoothing parameter",
    "VERSION=1.3",		";Josh Barnes  3 September 2002",
    NULL,
};

string bodytags[MaxBodyFields] = {
    PosTag, VelTag, MassTag, PhiTag, AccTag, NULL
};

local void sum1force(bodyptr, int, int, real);

int main(int argc, string argv[])
{
    stream istr, ostr;
    real eps2, tnow;
    int nforce, nbody, i;
    bodyptr btab = NULL;
    string intags[MaxBodyFields];

    initparam(argv, defv);
    layout_body(bodytags, Precision, NDIM);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    eps2 = rsqr(getdparam("eps"));
    nforce = -1;				/* use nforce as flag...    */
    while (get_snap(istr, &btab, &nbody, &tnow, intags, FALSE)) {
        if (nforce == -1 && ! set_member(intags, MassTag))
	    error("%s: Mass data missing from 1st snapshot\n", getargv0());
	if (! set_member(intags, PosTag))
	    error("%s: %s data missing\n", getargv0(), PosTag);
	if (nforce == -1)
	   nforce = (strnull(getparam("nforce")) ?
		       nbody : MIN(getiparam("nforce"), nbody));
	for (i = 0; i < nforce; i++)
	    sum1force(btab, nbody, i, eps2);
	put_snap(ostr, &btab, &nforce, &tnow, bodytags);
    }
    return (0);
}

local void sum1force(bodyptr btab, int nbody, int i0, real eps2)
{
    bodyptr bp0, bp;
    double phi0, acc0[NDIM];
    vector pos0, dr;
    real dr2, drab, mri, mr3i;
    int j;

    bp0 = NthBody(btab, i0);
    phi0 = 0.0;
    CLRV(acc0);
    SETV(pos0, Pos(bp0));
    for (j = 0; j < nbody; j++)
	if (j != i0) {
	    bp = NthBody(btab, j);
	    DOTPSUBV(dr2, dr, Pos(bp), pos0);
	    dr2 += eps2;
	    drab = rsqrt(dr2);
	    mri = Mass(bp) / drab;
	    phi0 -= mri;
	    mr3i = mri / dr2;
	    ADDMULVS(acc0, dr, mr3i);
	}
    Phi(bp0) = phi0;
    SETV(Acc(bp0), acc0);
}
