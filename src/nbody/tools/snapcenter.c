/*
 * snapcenter.c: transform snapshot file to coordinates centered
 * on weighted subset of particles.
 */

#include "stdinc.h"
#include "getparam.h"
#include "strset.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "buildmap.h"
#include "snapcenter.h"
#include <unistd.h>

string defv[] = {		";Center snapshot on weighted bodies",
  "in=???",                     ";Input snapshot file name",
  "out=???",                    ";Output snapshot file name",
  "times=all",                  ";Range of times to process",
  "weight=1.0",			";C language expression weighting bodies.",
				";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "coords=" PosTag "," VelTag,	";Coordinate data to center",
  "require=",			";Input items required",
  "produce=",			";Output items produced",
  "passall=true",		";If true, pass on input data",
  "seed=",			";Seed for random number generator",
  "VERSION=2.3",		";Josh Barnes  9 Sep 2014",
  NULL,
};

stream execmap(string);				// start snapmap process
void del_tag(string *, string *, string);	// remove tag from list

string names[2] = { "Weight",  NULL };
string exprs[2] = { NULL,      NULL };
string types[2] = { RealType,  NULL };

#define WeightField  phatbody[NewBodyFields+0]
#define Weight(b)  SelectReal(b, WeightField.offset)

int main(int argc, string argv[])
{
  string prog, coords, itags[MaxBodyFields], otags[MaxBodyFields];
  stream xstr, ostr;
  bodyptr btab = NULL, bp;
  int nbody;
  real tnow;
  vector cmpos, cmvel, cmacc;

  initparam(argv, defv);
  exprs[0] = getparam("weight");
  prog = tempnam("/tmp", "sm");
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog);
  if (get_tag_ok(xstr, "History"))
    skip_item(xstr);
  get_history(xstr);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  coords = getparam("coords");
  new_field(&WeightField, RealType, "Weight");
  new_field(&WeightField + 1, NULL, NULL);
  while (get_snap(xstr, &btab, &nbody, &tnow, itags, TRUE)) {
    if (scanopt(coords, PosTag) && set_member(itags, PosTag)) {
      snapcmpos(cmpos, btab, nbody, WeightField.offset);
      for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	SUBV(Pos(bp), Pos(bp), cmpos);
      }
      eprintf("[%s: centroid position: %f,%f,%f]\n", getprog(),
	      cmpos[0], cmpos[1], cmpos[2]);
    }
    if (scanopt(coords, VelTag) && set_member(itags, VelTag)) {
      snapcmvel(cmvel, btab, nbody, WeightField.offset);
      for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	SUBV(Vel(bp), Vel(bp), cmvel);
      }
      eprintf("[%s: centroid velocity: %f,%f,%f]\n", getprog(),
	      cmvel[0], cmvel[1], cmvel[2]);
    }
    if (scanopt(coords, AccTag) && set_member(itags, AccTag)) {
      snapcmacc(cmacc, btab, nbody, WeightField.offset);
      for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	SUBV(Acc(bp), Acc(bp), cmacc);
      }
      eprintf("[%s: cen. acceleration: %f,%f,%f]\n", getprog(),
	      cmacc[0], cmacc[1], cmacc[2]);
    }
    del_tag(otags, itags, "Weight");
    put_snap(ostr, &btab, &nbody, &tnow, otags);
  }
  strclose(ostr);
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getprog(), prog);
  return (0);
}

//  execmap: start snapmap subprocess, and return snapmap output stream.
//  ____________________________________________________________________

stream execmap(string prog)
{
  int handle[2];
  char handbuf[32], produce[512];

  pipe(handle);
  if (fork() == 0) {				// if this is child process
    close(handle[0]);
    sprintf(handbuf, "-%d", handle[1]);
    sprintf(produce, "%s,Weight", getparam("produce"));
    execl(prog, getprog(), getparam("in"), handbuf, getparam("times"),
	  getparam("require"), produce, getparam("passall"),
	  getparam("seed"), NULL);
    error("%s: execl %s failed\n", getprog(), prog);
  }
  close(handle[1]);
  sprintf(handbuf, "-%d", handle[0]);
  return (stropen(handbuf, "r"));
}

void del_tag(string *olist, string *ilist, string tag)
{
  string *op, *ip;

  for (op = olist, ip = ilist; *ip != NULL; ip++)
    if (! streq(*ip, tag))
      *op++ = *ip;
  *op = NULL;
}
