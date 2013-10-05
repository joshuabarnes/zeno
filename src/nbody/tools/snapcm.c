/*
 * SNAPCM.C: read a snapshot file and calculate center-of-mass.
 */

#include "stdinc.h"
#include "getparam.h"
#include "filestruct.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"

string defv[] = {	       	";Calculate center-of-mass coords",
    "in=???",			";Input file of N-body snapshots",
    "VERSION=2.3",		";Josh Barnes  29 March 2000",
    NULL,
};

string bodytags[] = { PosTag, VelTag, MassTag, NULL };

int main(int argc, string argv[])
{
    stream istr;
    bodyptr btab = NULL, bp;
    int nbody;
    real tnow, mtot = 0.0;
    vector cmr, cmv;
    string intags[MaxBodyFields];
    bool firstloop = TRUE;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    layout_body(bodytags, Precision, NDIM);
    while (get_snap(istr, &btab, &nbody, &tnow, intags, FALSE)) {
	if (set_member(intags, MassTag)) {
	    mtot = 0.0;
	    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
		mtot += Mass(bp);
	}
	if (firstloop)
	    printf("#%9s %9s %9s %9s %9s %9s %9s\n",
		   "time", "xcm", "ycm", "zcm", "vxcm", "vycm", "vzcm");
	if (set_member(intags, PosTag) && set_member(intags, VelTag)) {
	    if (mtot == 0.0)
		error("%s: mass data missing\n", getargv0());
	    snapcmpos(cmr, btab, nbody, MassField.offset);
	    snapcmvel(cmv, btab, nbody, MassField.offset);
	    printf(" %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
		   tnow, cmr[0], cmr[1], cmr[2], cmv[0], cmv[1], cmv[2]);
	}
	firstloop = FALSE;
    }
    return (0);
}
