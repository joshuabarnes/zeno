/*
 * SNAPTURN.C: Read one snapshot, produce stream of rotated snapshots.
 */

#include "stdinc.h"
#include "string.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Produce stream of rotated snapshots",
    "in=???",			";Input snapshot file",
    "out=???",			";Stream of output snapshots",
    "time=all",			";Time value for snapshot",
    "nsnap=361",		";Number of snapshots output",
    "angle=1.0",		";Rotation angle between snaps",
    "axis=z",			";Axis of rotation: x, y, z",
    "produce=*",		";List of items to produce",
    "VERSION=1.2",		";Josh Barnes  12 April 2000",
    NULL,
};

void snaprotate(bodyptr, bodyptr, int, real, char);

int main(int argc, string argv[])
{
    stream istr, ostr;
    string time, *produce, iotags[MaxBodyFields];
    bodyptr btab1 = NULL, btab2 = NULL;
    int nbody, i;
    real tnow, angle;
    bool expand;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    time = getparam("time");
    expand = streq(getparam("produce"), "*");
    if (! expand) {
	produce = burststring(getparam("produce"), ", ");
	layout_body(produce, Precision, NDIM);
    }
    if (! get_snap(istr, &btab1, &nbody, &tnow, iotags, expand, time))
        error("%s: no data in input file\n", getargv0());
    strclose(istr);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    btab2 = (bodyptr) allocate(nbody * SizeofBody);
    angle = (PI / 180.0) * getdparam("angle");
    for (i = 0; i < getiparam("nsnap"); i++) {
	snaprotate(btab2, btab1, nbody, i * angle, getparam("axis")[0]);
	put_snap(ostr, &btab2, &nbody, &tnow, iotags);
    }
    strclose(ostr);
    return (0);
}

void snaprotate(bodyptr btab2, bodyptr btab1, int nbody,
		real theta, char axis)
{
    matrix rmat;
    int i;
    bodyptr bp1, bp2;

    SETMI(rmat);
    switch (axis) {
      case 'x':
      case 'X':
	rmat[1][1] =   rcos(theta);
	rmat[2][1] = - rsin(theta);
	rmat[1][2] =   rsin(theta);
	rmat[2][2] =   rcos(theta);
	break;
      case 'y':
      case 'Y':
	rmat[2][2] =   rcos(theta);
	rmat[0][2] = - rsin(theta);
	rmat[2][0] =   rsin(theta);
	rmat[0][0] =   rcos(theta);
	break;
      case 'z':
      case 'Z':
	rmat[0][0] =   rcos(theta);
	rmat[1][0] = - rsin(theta);
	rmat[0][1] =   rsin(theta);
	rmat[1][1] =   rcos(theta);
	break;
      default:
	error("%s: unknown axis %c\n", getargv0(), axis);
    }
    for (i = 0; i < nbody; i++) {
	bp1 = NthBody(btab1, i);
	bp2 = NthBody(btab2, i);
	memcpy(bp2, bp1, SizeofBody);
	MULMV(Pos(bp2), rmat, Pos(bp1));
	MULMV(Vel(bp2), rmat, Vel(bp1));
    }
}
