/*
 * snaptestdisk.c: set up a uniform-density test disk in a spherical galaxy.
 */

#include <string.h>
#include "stdinc.h"
#include "strset.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "snapcenter.h"

string defv[] = {       ";Make test disk in spherical system",
    "in=???",           ";Spheroid N-body input file",
    "out=",             ";Resulting N-body output file",
    "rmin=0.1",         ";Inner radius of test disk",
    "rmax=0.4",         ";Outer radius of test disk",
    "eps=0.0",          ";Softening for rotation curve",
    "sigma=0.0",	";Random velocity dispersion",
    "ndisk=2048",       ";Number of disk particles",
    "seed=54321",       ";Seed for random number generator",
    "nosphr=false",	";If true, only output disk",
    "zerocm=false",     ";If true, translate to centroid",
    "VERSION=2.1",      ";Josh Barnes  12 May 2012",
    NULL,
};

string bodyfields[] = { PosTag, VelTag, MassTag, NULL };

bodyptr spheroid = NULL, galaxy = NULL;

int nspheroid, ngalaxy;

void readsphr(string);
void sphrprof(void);
void testdisk(void);
void writegal(string);

int main(int argc, string argv[])
{
  initparam(argv, defv);
  layout_body(bodyfields, Precision, NDIM);
  readsphr(getparam("in"));
  sphrprof();
  testdisk();
  if (! strnull(getparam("out")))
    writegal(getparam("out"));
  return (0);
}

/*
 * READSPHR: read spheroid model from input.
 */

void readsphr(string name)
{
    stream instr;
    real tsnap;
    string infields[MaxBodyFields];

    instr = stropen(name, "r");
    get_history(instr);
    if (! get_snap(instr, &spheroid, &nspheroid, &tsnap, infields, FALSE))
        error("readsphr in %s: no data in input file\n", getargv0());
    strclose(instr);
    if (! (set_member(infields, MassTag) && set_member(infields, PosTag)))
        error("readsphr in %s: required data missing\n", getargv0());
}

/*
 * WRITEGAL: write galaxy model to output.
 */

void writegal(string name)
{
    stream outstr;
    real tsnap = 0.0;

    outstr = stropen(name, "w");
    put_history(outstr);
    put_snap(outstr, &galaxy, &ngalaxy, &tsnap, bodyfields);
    strclose(outstr);
}

/*
 * SPHRPROF: process spheroid to generate mass profile table [rsph, msph].
 */

#define NTAB  257                               /* length of profile table */

real rsph[NTAB];                                /* radii of selected bodies */
real msph[NTAB];                                /* mass in rsph, inclusive */
real mcof[3*NTAB];                              /* spline coefs for mass */

int rankrad(const void *, const void *);

void sphrprof(void)
{
    bodyptr *rsort;
    int skip, j, i;
    real r;

    if (nspheroid > 1) {
	rsort = (bodyptr *) allocate(nspheroid * sizeof(bodyptr));
	for (i = 0; i < nspheroid; i++)
	    rsort[i] = NthBody(spheroid, i);
	qsort(rsort, nspheroid, sizeof(bodyptr), rankrad);
        skip = (int) rceil(MAX(((real) nspheroid) / (NTAB - 1), 1.0));
	msph[0] = rsph[0] = 0.0;
	for (j = 1; j < NTAB; j++) {
	    i = skip * j - 1;
	    rsph[j] = (i < nspheroid ? absv(Pos(rsort[i])) : 1.0 + rsph[j-1]);
	    msph[j] = 0.0;
	}
	for (i = 0; i < nspheroid; i++) {
	    j = i / skip + 1;
	    if (j > NTAB-1)
		error("%s.sphrprof: table overflow: i = %d  skip = %d\n",
			  getargv0(), i, skip);
	    msph[j] = msph[j] + Mass(rsort[i]);
	}
	for (j = 1; j < NTAB; j++)
	    msph[j] += msph[j - 1];
    } else {
        for (j = 0; j < NTAB; j++) {
	    rsph[j] = 1.0 * j;
	    msph[j] = Mass(NthBody(spheroid, 0));
	}
    }
    eprintf("[%s.sphrprof: rsph = %f %f %f ... %f %f]\n", getargv0(),
	    rsph[0], rsph[1], rsph[2], rsph[NTAB-2], rsph[NTAB-1]);
    spline(mcof, rsph, msph, NTAB);
}

int rankrad(const void *x, const void *y)
{
    real rxsq = dotvp(Pos(*((bodyptr *) x)), Pos(*((bodyptr *) x)));
    real rysq = dotvp(Pos(*((bodyptr *) y)), Pos(*((bodyptr *) y)));

    return (rxsq < rysq ? -1 : rxsq > rysq ? 1 : 0);
}

/*
 * TESTDISK: use tabulated spheroid mass profile to make a uniform
 * density test disk.  Note that the softening correction is only
 * exact for a central point mass.
 */

void testdisk(void)
{
    int ndisk, i;
    real rmin2, rmax2, eps2, sigma, r_i, theta_i, m_i, v_i;
    bodyptr gp, sp;

    ndisk = getiparam("ndisk");
    ngalaxy = ndisk + (getbparam("nosphr") ? 0 : nspheroid);
    galaxy = (bodyptr) allocate(ngalaxy * SizeofBody);
    rmin2 = rsqr(getdparam("rmin"));
    rmax2 = rsqr(getdparam("rmax"));
    eps2 = rsqr(getdparam("eps"));
    sigma = getdparam("sigma");
    init_random(getiparam("seed"));
    for (i = 0; i < ndisk; i++) {			/* build disk       */
        gp = NthBody(galaxy, i);
        Mass(gp) = 0.0;                                 /* set mass to zero */
        r_i = rsqrt(rmin2 + i * (rmax2 - rmin2) / (ndisk - 1.0));
        theta_i = xrandom(0.0, TWO_PI);
        Pos(gp)[0] = r_i * rsin(theta_i);               /* set positions    */
        Pos(gp)[1] = r_i * rcos(theta_i);
        Pos(gp)[2] = 0.0;
        if (r_i < rsph[NTAB-1])
            m_i = seval(r_i, &rsph[0], &msph[0], &msph[NTAB], NTAB);
        else
            m_i = msph[NTAB-1];
	v_i = rsqrt(MAX(m_i, 0.0) * r_i*r_i / rpow(r_i*r_i + eps2, 1.5));
							/* set velocities   */
        Vel(gp)[0] = grandom(  v_i * rcos(theta_i), sigma);
        Vel(gp)[1] = grandom(- v_i * rsin(theta_i), sigma);
        Vel(gp)[2] = grandom(                  0.0, sigma);
    }
    if (! getbparam("nosphr"))
	for (i = 0; i < nspheroid; i++) {		/* append spheroid  */
	    sp = NthBody(spheroid, i);
	    gp = NthBody(galaxy, ndisk + i);
	    memcpy(gp, sp, SizeofBody);
	}
    if (getbparam("zerocm"))
        snapcenter(galaxy, ngalaxy, MassField.offset);
}

