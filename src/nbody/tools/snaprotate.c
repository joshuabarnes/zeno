/*
 * SNAPROTATE.C: read a snapshot file and rotate particle coordinates.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"

#include <string.h>
#include <ctype.h>

string defv[] = {		";Rotate N-body configuration",
    "in=???",			";Input snapshot file",
    "out=???",			";Output snapshot file",
    "times=all",		";Range of times to rotate",
    "order=xyz",		";Order to do rotations",
    "thetax=0.0",		";Angle of rotation about x (in deg)",
    "thetay=0.0",		";Angle of rotation about y",
    "thetaz=0.0",		";Angle of rotation about z",
    "invert=false",		";If true, invert transformation",
    "produce=*",		";List of items to produce",
    "VERSION=3.0",		";Josh Barnes  6 January 1999",
    NULL,
};

void snaprotate(bodyptr, int, string *, string, bool, real, real, real);
void rotate(bodyptr, int, string *, char, real, real, real);
void xmatrix(matrix, real), ymatrix(matrix, real), zmatrix(matrix, real);
void rotatevec(vector, matrix);

int main(int argc, string argv[])
{
    stream istr, ostr;
    string times, *produce, iotags[MaxBodyFields];
    bodyptr btab = NULL;
    int nbody;
    real tsnap;
    bool expand;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    times = getparam("times");
    expand = streq(getparam("produce"), "*");
    if (! expand) {
	produce = burststring(getparam("produce"), ", ");
	layout_body(produce, Precision, NDIM);
    }
    while (get_snap_t(istr, &btab, &nbody, &tsnap, iotags, expand, times)) {
	snaprotate(btab, nbody, iotags, getparam("order"),
		   getbparam("invert"), getdparam("thetax"),
		   getdparam("thetay"), getdparam("thetaz"));
	put_snap(ostr, &btab, &nbody, &tsnap, iotags);
    }
    strclose(ostr);
    return (0);
}

/*
 * SNAPROTATE: perform rotation about specified axes.
 */

void snaprotate(bodyptr btab, int nbody, string *tags, string order,
		bool invert, real thetax, real thetay, real thetaz)
{
    if (strlen(order) != 3)
	error("%s: order must be 3 chars long\n", getargv0());
    if (! invert) {
	rotate(btab, nbody, tags, order[0], thetax, thetay, thetaz);
	rotate(btab, nbody, tags, order[1], thetax, thetay, thetaz);
	rotate(btab, nbody, tags, order[2], thetax, thetay, thetaz);
    } else {
	rotate(btab, nbody, tags, order[2], -thetax, -thetay, -thetaz);
	rotate(btab, nbody, tags, order[1], -thetax, -thetay, -thetaz);
	rotate(btab, nbody, tags, order[0], -thetax, -thetay, -thetaz);
    }
}

void rotate(bodyptr btab, int nbody, string *tags, char axis,
	    real thetax, real thetay, real thetaz)
{
    matrix rmat;
    bodyptr bp;

    switch (tolower(axis)) {
      case 'x':
	xmatrix(rmat, thetax);
	break;
      case 'y':
	ymatrix(rmat, thetay);
	break;
      case 'z':
	zmatrix(rmat, thetaz);
	break;
      default:
	error("%s: unknown axis %c\n", getargv0(), axis);
    }
    if (set_member(tags, PosTag))
	for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	    rotatevec(Pos(bp), rmat);
    if (set_member(tags, VelTag))
	for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	    rotatevec(Vel(bp), rmat);
    if (set_member(tags, AccTag))
	for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	    rotatevec(Acc(bp), rmat);
    if (set_member(tags, AuxVecTag))
	for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	    rotatevec(AuxVec(bp), rmat);
}

#define DEG2RAD 0.017453

void xmatrix(matrix rmat, real theta)
{
    real s = rsin(DEG2RAD * theta), c = rcos(DEG2RAD * theta);

    rmat[0][0] = 1.0;    rmat[0][1] = 0.0;    rmat[0][2] = 0.0;
    rmat[1][0] = 0.0;    rmat[1][1] =  c ;    rmat[1][2] =  s ;
    rmat[2][0] = 0.0;    rmat[2][1] = -s ;    rmat[2][2] =  c ;
}

void ymatrix(matrix rmat, real theta)
{
    real s = rsin(DEG2RAD * theta), c = rcos(DEG2RAD * theta);

    rmat[0][0] =  c ;    rmat[0][1] = 0.0;    rmat[0][2] = -s ;
    rmat[1][0] = 0.0;    rmat[1][1] = 1.0;    rmat[1][2] = 0.0;
    rmat[2][0] =  s ;    rmat[2][1] = 0.0;    rmat[2][2] =  c ;
}

void zmatrix(matrix rmat, real theta)
{
    real s = rsin(DEG2RAD * theta), c = rcos(DEG2RAD * theta);

    rmat[0][0] =  c ;    rmat[0][1] =  s ;    rmat[0][2] = 0.0;
    rmat[1][0] = -s ;    rmat[1][1] =  c ;    rmat[1][2] = 0.0;
    rmat[2][0] = 0.0;    rmat[2][1] = 0.0;    rmat[2][2] = 1.0;
}

void rotatevec(vector vec, matrix mat)
{
    vector tmp;

    MULMV(tmp, mat, vec);
    SETV(vec, tmp);
}
