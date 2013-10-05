/*
 * snapregister.c: register snapshot centroids.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Register snapshot centroids",
    "in=???",			";Input snapshot, with " TypeTag " data",
    "out=???",			";Output snapshot, scaled and rotated",
    "sphrtype=???",		";" TypeTag " value tagging spheroids",
    "ncent=1024",		";No. of bodies used to find centroids",
    "dist=2.0",			";Scaled distance between centroids",
    "angle=0.0",		";Centroid position angle wrt X-axis",
    "VERSION=1.0",		";Josh Barnes  19 December 2012",
    NULL,
};

void findcenters(vector cent1, vector cent2,
		 bodyptr btab, int nbody, int sphrtype, int ncent);

void snapregister(bodyptr btab, int nbody, vector cent1, vector cent2,
		  real dist, real angle);

int main(int argc, string argv[])
{
  stream istr, ostr = NULL;
  bodyptr btab = NULL;
  int nbody;
  real tnow;
  string iotags[MaxBodyFields];
  vector cent1, cent2;

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  while (get_snap(istr, &btab, &nbody, &tnow, iotags, TRUE)) {
    if (ostr == NULL && ! set_member(iotags, TypeTag))
      error("%s: %s input data missing\n", getargv0(), TypeTag);
    if (! set_member(iotags, PosTag))
      error("%s: %s input data missing\n", getargv0(), PosTag);
    if (! set_member(iotags, VelTag))
      error("%s: %s input data missing\n", getargv0(), VelTag);
    findcenters(cent1, cent2, btab, nbody,
		getiparam("sphrtype"), getiparam("ncent"));
    snapregister(btab, nbody, cent1, cent2,
		 getdparam("dist"), (PI/180) * getdparam("angle"));
    if (ostr == NULL) {
      ostr = stropen(getparam("out"), "w");
      put_history(ostr);
    }
    put_snap(ostr, &btab, &nbody, &tnow, iotags);
  }
  return (0);
}

void findcenters(vector cent1, vector cent2,
		 bodyptr btab, int nbody, int sphrtype, int ncent)
{
  int block = 0, ncent1 = 0, ncent2 = 0;
  bodyptr bp;

  CLRV(cent1);
  CLRV(cent2);
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    if (Type(bp) == sphrtype) {
      if (block % 2 == 0)
	block++;
      if ((block == 1 ? ncent1 : ncent2) < ncent) {
	(* (block == 1 ? &ncent1 : &ncent2))++;
	if (block == 1) {
	  ADDV(cent1, cent1, Pos(bp));
	} else {
	  ADDV(cent2, cent2, Pos(bp));
	}
      }
    } else {
      if (block % 2 == 1)
	block++;
    }
  }
  eprintf("[%s: ncent1 = %d  ncent2 = %d]\n", getprog(), ncent1, ncent2);
  if (ncent1 == 0 || ncent2 == 0)
    error("%s: ncent1 = %d  ncent2 = %d\n", getprog(), ncent1, ncent2);
  DIVVS(cent1, cent1, ncent1);
  DIVVS(cent2, cent2, ncent2);
}

void snapregister(bodyptr btab, int nbody, vector cent1, vector cent2, 
		  real dist1, real angle1)
{
  vector xyvec, xyoff, tmp;
  int block = 0, ncent1 = 0, ncent2 = 0;
  bodyptr bp;
  real dist0, angle0, s, c;

  SUBV(xyvec, cent1, cent2);
  xyvec[2] = 0.0;
  ADDV(xyoff, cent1, cent2);
  xyoff[2] = 0.0;
  DIVVS(xyoff, xyoff, 2.0);
  dist0 = absv(xyvec);
  angle0 = ratan2(xyvec[1], xyvec[0]);
  eprintf("[%s: dist0 = %f  angle0 = %f]\n", getprog(),
	  dist0, (180/PI) * angle0);
  if (dist0 == 0.0)
    error("%s: projected seperation vector vanishes\n", getprog());
  s = rsin(angle1 - angle0);
  c = rcos(angle1 - angle0);
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    SUBV(tmp, Pos(bp), xyoff);
    MULVS(tmp, tmp, dist1 / dist0);
    Pos(bp)[0] = c * tmp[0] - s * tmp[1];
    Pos(bp)[1] = s * tmp[0] + c * tmp[1];
    Pos(bp)[2] = tmp[2];
    SETV(tmp, Vel(bp));
    Vel(bp)[0] = c * tmp[0] - s * tmp[1];
    Vel(bp)[1] = s * tmp[0] + c * tmp[1];
    Vel(bp)[2] = tmp[2];
  }
}
