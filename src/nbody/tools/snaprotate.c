/*
 * snaprotate.c: read a snapshot file and rotate particle coordinates.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include <string.h>
#include <ctype.h>

string defv[] = {		";Rotate N-body configuration",
  "in=???",			";Input snapshot file",
  "out=???",			";Output snapshot file",
  "times=all",			";Range of times to rotate",
  "order=xyz",			";Order to do rotations",
  "thetax=0.0",			";Angle of rotation about x (in deg)",
  "thetay=0.0",			";Angle of rotation about y",
  "thetaz=0.0",			";Angle of rotation about z",
  "invert=false",		";If true, invert transformation",
  "vectors=" PosTag "," VelTag "," AccTag "," AuxVecTag,
				";List of vectors to rotate",
  "produce=*",			";List of items to produce",
  "VERSION=3.1",		";Josh Barnes  2 February 2015",
  NULL,
};

//  snapshot: prototype of high-level snapshot structure, and some macros.
//  ______________________________________________________________________

typedef struct {
  bodyptr bodies;
  int nbody;
  real time;
} snapshot;

#define get_snapshot_timed(str, snap, tags, exp, times)			\
  get_snap_t(str, &(snap).bodies, &(snap).nbody, &(snap).time,		\
             tags, exp, times)

#define put_snapshot(str, snap, tags)					\
  put_snap(str, &(snap).bodies, &(snap).nbody, &snap.time, tags)

#define for_all_bodies(bp, snap)					\
  for (bp = (snap).bodies;						\
       bp < NthBody((snap).bodies, (snap).nbody);			\
       bp = NextBody(bp))

void snaprotate(snapshot *snap, string *tags, string *vecs, string order,
		bool invert, real thetax, real thetay, real thetaz);

void rotate(snapshot *snap, string *tags, string *vecs, char axis,
	    real thetax, real thetay, real thetaz);

void xmatrix(matrix rmat, real theta);
void ymatrix(matrix rmat, real theta);
void zmatrix(matrix rmat, real theta);
void rotatevec(vector vec, matrix mat);

int main(int argc, string argv[])
{
  stream istr, ostr;
  string times, *vecs, *produce, iotags[MaxBodyFields];
  bool expand;
  snapshot snap = { NULL, 0, 0.0 };

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  times = getparam("times");
  vecs = burststring(getparam("vectors"), ", ");
  expand = streq(getparam("produce"), "*");
  if (! expand) {
    produce = burststring(getparam("produce"), ", ");
    layout_body(produce, Precision, NDIM);
  }
  while (get_snapshot_timed(istr, snap, iotags, expand, times)) {
    eprintf("[%s: rotating time %f]\n", getprog(), snap.time);
    snaprotate(&snap, iotags, vecs, getparam("order"), getbparam("invert"),
	       getdparam("thetax"), getdparam("thetay"), getdparam("thetaz"));
    put_snapshot(ostr, snap, iotags);
    skip_history(istr);
  }
  strclose(ostr);
  return (0);
}

//  snaprotate: perform rotation about specified axes.
//  __________________________________________________

void snaprotate(snapshot *snap, string *tags, string *vecs, string order,
		bool invert, real thetax, real thetay, real thetaz)
{
  if (strlen(order) != 3)
    error("%s: order must be 3 chars long\n", getargv0());
  if (! invert) {
    rotate(snap, tags, vecs, order[0], thetax, thetay, thetaz);
    rotate(snap, tags, vecs, order[1], thetax, thetay, thetaz);
    rotate(snap, tags, vecs, order[2], thetax, thetay, thetaz);
  } else {
    rotate(snap, tags, vecs, order[2], -thetax, -thetay, -thetaz);
    rotate(snap, tags, vecs, order[1], -thetax, -thetay, -thetaz);
    rotate(snap, tags, vecs, order[0], -thetax, -thetay, -thetaz);
  }
}

//  rotate: perform rotation about one axis.
//  ________________________________________

void rotate(snapshot *snap, string *tags, string *vecs, char axis,
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
  if (set_member(tags, PosTag) && set_member(vecs, PosTag))
    for_all_bodies(bp, *snap)
      rotatevec(Pos(bp), rmat);
  if (set_member(tags, VelTag) && set_member(vecs, VelTag))
    for_all_bodies(bp, *snap)
      rotatevec(Vel(bp), rmat);
  if (set_member(tags, AccTag) && set_member(vecs, AccTag))
    for_all_bodies(bp, *snap)
      rotatevec(Acc(bp), rmat);
  if (set_member(tags, AuxVecTag) && set_member(vecs, AuxVecTag))
    for_all_bodies(bp, *snap)
      rotatevec(AuxVec(bp), rmat);
}

#define DEG2RAD  (PI / 180.0)

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
