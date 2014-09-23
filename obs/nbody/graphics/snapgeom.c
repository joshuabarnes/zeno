/*
 * SNAPGEOM.C: extract disk geometry and tracks.
 */

#include "stdinc.h"
#include "getparam.h"
#include "filestruct.h"
#include "vectmath.h"
#include "phatbody.h"

string defv[] = {	       	";Extract disk geometry and tracks",
    "key=???",			";Keyed initial data input file",
    "out=???",			";Encounter geometry output file",
    "VERSION=1.0",		";Josh Barnes  19 May 1998",
    NULL,
};

string bodytags[] = { PosTag, VelTag, PhiTag, KeyTag, NULL };

void geometry(bodyptr, int);
void spinvect(vector, bodyptr, int);
void halocent(bodyptr, int);
void writegeom(stream);

int main(int argc, string argv[])
{
    stream istr, ostr;
    bodyptr btab = NULL;
    int nbody;
    real tnow;
    string intags[MaxBodyFields];

    initparam(argv, defv);
    istr = stropen(getparam("key"), "r");
    get_history(istr);
    layout_body(bodytags, Precision, NDIM);
    if (! get_snap(istr, &btab, &nbody, &tnow, intags, FALSE))
	error("%s: get_snap failed\n", getargv0());
    if (! (set_member(intags, PosTag) && set_member(intags, VelTag)))
	error("%s: coordinate data missing\n", getargv0());
    if (! set_member(intags, PhiTag))
	error("%s: potential data missing\n", getargv0());
    if (! set_member(intags, KeyTag))
	error("%s: key data missing\n", getargv0());
    geometry(btab, nbody);
    halocent(btab, nbody);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    writegeom(ostr);
    strclose(ostr);
    return(0);
}

#define MAXGAL  12				/* max. galaxies tracked    */
#define MAXSEG  33				/* max. segments per galaxy */
#define TOTSEG  (MAXGAL * MAXSEG)		/* max. segments in total   */

int ngal;					/* actual no. of galaxies   */
int nseg[MAXGAL+1];				/* segments for each galaxy */
int iseg[TOTSEG+1];				/* indicies of all segments */
vector jseg[TOTSEG];				/* spin vects of all segs.  */
vector rgal[MAXGAL];				/* position of each galaxy  */

void geometry(bodyptr btab, int nbody)
{
    int ns, i1, i2, kval;

    ngal = 0;					/* count number of galaxies */
    nseg[ngal] = ns = 0;			/* and numbers of segments  */
    for (i1 = 0; i1 < nbody; i1 = i2) {		/* loop over all bodies     */
	kval = Key(NthBody(btab, i1));		/* get key for this segment */
	i2 = i1;				/* scan ahead from here...  */
	while (i2 < nbody && kval == Key(NthBody(btab, i2)))
	    i2++;				/* to find end of segment   */
	if (kval == 0) {			/* start of halo segmant?   */
	    CLRV(jseg[ns]);			/* zero halo spin vector    */
	    ngal++;				/* count another galaxy     */
	    nseg[ngal] = ns;			/* init galaxy seg. count   */
	} else					/* start of disk segment?   */
	    spinvect(jseg[ns], NthBody(btab, i1), i2 - i1);
						/* record disk spin vector  */
	nseg[ngal]++;				/* count another segment    */
	iseg[ns++] = i1;			/* record start of segment  */
    }
    iseg[ns] = nbody;				/* record end of last seg.  */
}

void spinvect(vector jvec, bodyptr btab, int nb)
{
    vector cmr, cmv, r, v, jtmp;
    bodyptr bp;
    real jmag;

    CLRV(cmr);
    CLRV(cmv);
    for (bp = btab; bp < NthBody(btab, nb); bp = NextBody(bp)) {
	ADDMULVS(cmr, Pos(bp), 1.0 / nb);
	ADDMULVS(cmv, Vel(bp), 1.0 / nb);
    }
    CLRV(jvec);
    for (bp = btab; bp < NthBody(btab, nb); bp = NextBody(bp)) {
	SUBV(r, Pos(bp), cmr);
	SUBV(v, Vel(bp), cmv);
	CROSSVP(jtmp, v, r);
	jmag = absv(jtmp);
	ADDMULVS(jvec, jtmp, 1.0 / jmag);
    }
    jmag = absv(jvec);
    DIVVS(jvec, jvec, jmag);
}

void halocent(bodyptr btab, int nbody)
{
    int n, i1, i2, i;
    real phimin;

    for (n = 0; n < ngal; n++) {
	i1 = iseg[ nseg[n] ];
	i2 = iseg[nseg[n]+1];
	phimin = Phi(NthBody(btab, i1));
	SETV(rgal[n], Pos(NthBody(btab, i1)));
	for (i = i1 + 1; i < i2; i++)
	    if (Phi(NthBody(btab, i)) < phimin) {
		phimin = Phi(NthBody(btab, i));
		SETV(rgal[n], Pos(NthBody(btab, i)));
	    }
    }
}

void writegeom(stream ostr)
{
    int ns = nseg[ngal];

    put_set(ostr, "EncounterGeometry");
    put_data(ostr, "NGalaxy", IntType, &ngal, 0);
    put_data(ostr, "NSegment", IntType, nseg, ngal + 1, 0);
    put_data(ostr, "ISegment", IntType, iseg, ns + 1, 0);
    put_data(ostr, "JSegment", RealType, jseg, ns, NDIM, 0);
    put_data(ostr, "RGalaxy", RealType, rgal, ngal, NDIM, 0);
    put_tes(ostr, "EncounterGeometry");
}
