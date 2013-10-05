/*
 * MAKESNAP: make a configuration of bodies (formerly MKCONFIG).
 */

#include "stdinc.h"
#include <string.h>
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "snapcenter.h"

string defv[] = {		";Make a configuration of bodies",
    "out=???",			";Output snapshot file name",
    "shape=shell",		";shell,ball,cusp,gauss,grid,line,ring",
    "radius=1.0",		";Characteristic or maximum radius",
    "vshape=shell",		";Like shape, but in velocity space",
    "speed=1.0",		";Characteristic or maximum velocity",
    "mass=1.0",			";Total mass of configuration",
    "nbody=4096",		";Number of bodies in configuration",
    "copies=1",			";Number of copies to generate",
    "seed=123",			";Seed for random numbers",
    "zerocm=false",		";Translate center of mass to origin",
    "VERSION=3.1",		";Josh Barnes  11 May 2012",
    NULL,
};

string bodytags[] = { PosTag, VelTag, MassTag, NULL };

int nbody;

bodyptr btab = NULL;

void initbody(real);
void mkconfig(string, real, int);
void makeshell(real, int);
void makeball(real, int);
void makecusp(real, int);
void makegauss(real, int);
void makegrid(real, int);
void makeline(real, int);
void makering(real, int);

int main(int argc, string argv[])
{
    stream outstr;
    int n;
    real tzero = 0.0;

    initparam(argv, defv);
    nbody = getiparam("nbody");
    if (nbody < 1)
	error("%s: nbody = %d is absurd!\n", getargv0());
    layout_body(bodytags, Precision, NDIM);
    btab = (bodyptr) allocate(nbody * SizeofBody);
    init_random(getiparam("seed"));
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    for (n = 1; n <= getiparam("copies"); n++) {
	initbody(getdparam("mass"));
	mkconfig(getparam("shape"), getdparam("radius"), PosField.offset);
	mkconfig(getparam("vshape"), getdparam("speed"), VelField.offset);
	if (getbparam("zerocm"))
	    snapcenter(btab, nbody, MassField.offset);
	put_snap(outstr, &btab, &nbody, &tzero, bodytags);
    }
    strclose(outstr);
    return (0);
}

void initbody(real mass)
{
    bodyptr bp;

    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	Mass(bp) = mass / nbody;
	CLRV(Pos(bp));
	CLRV(Vel(bp));
    }
}

void mkconfig(string shape, real radius, int voff)
{
    if (streq(shape, "shell"))
	makeshell(radius, voff);
    else if (streq(shape, "ball"))
	makeball(radius, voff);
    else if (streq(shape, "cusp"))
	makecusp(radius, voff);
    else if (streq(shape, "gauss"))
	makegauss(radius, voff);
    else if (streq(shape, "grid"))
	makegrid(radius, voff);
    else if (streq(shape, "line"))
	makeline(radius, voff);
    else if (streq(shape, "ring"))
	makering(radius, voff);
    else
	error("%s: shape %s unknown\n", getargv0(), shape);
}

void makeshell(real radius, int voff)
{
    bodyptr bp;

    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	pickshell(SelectVect(bp, voff), NDIM, radius);
}

void makeball(real radius, int voff)
{
    bodyptr bp;

    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	pickball(SelectVect(bp, voff), NDIM, radius);
}

void makecusp(real radius, int voff)
{
    bodyptr bp;

    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	pickshell(SelectVect(bp, voff), NDIM, xrandom(0.0, radius));
}

void makegauss(real radius, int voff)
{
    bodyptr bp;

    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	SelectVect(bp, voff)[0] = grandom(0.0, radius);
	SelectVect(bp, voff)[1] = grandom(0.0, radius);
	SelectVect(bp, voff)[2] = grandom(0.0, radius);
    }
}

void makegrid(real radius, int voff)
{
    int nside, j, k;
    bodyptr bp;
    
    nside = (int) rsqrt((real) nbody);
    if (nside*nside != nbody)
	error("%s: nbody not square\n", getargv0());
    bp = btab;
    for (j = 0; j < nside; j++)
	for (k = 0; k < nside; k++) {
	    SelectVect(bp, voff)[0] = radius * (2.0 * j / (nside - 1.0) - 1.0);
	    SelectVect(bp, voff)[1] = radius * (2.0 * k / (nside - 1.0) - 1.0);
	    SelectVect(bp, voff)[2] = 0.0;
	    bp = NextBody(bp);
	}
}

void makeline(real radius, int voff)
{
    bodyptr bp;
    real x = 0.0;

    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	SelectVect(bp, voff)[0] = x;
	SelectVect(bp, voff)[1] = 0.0;
	SelectVect(bp, voff)[2] = 0.0;
	x += radius / (nbody - 1);
    }
}

void makering(real radius, int voff)
{
    bodyptr bp;
    real theta = 0.0;

    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	SelectVect(bp, voff)[0] = radius * rsin(theta);
	SelectVect(bp, voff)[1] = radius * rcos(theta);
	SelectVect(bp, voff)[2] = 0.0;
	theta += TWO_PI / nbody;
    }
}
