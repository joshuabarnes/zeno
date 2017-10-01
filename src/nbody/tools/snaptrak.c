/*
 * snaptrak.c: track centroids of user-specified particle groups.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "buildmap.h"
#include <unistd.h>

string defv[] = {		";Track centroids of specified groups",
  "in=???",			";Input snapshot file name",
  "out=???",			";Output snapshot file name",
  "group=???",			";Expression (C code) for group membership",
				";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "times=all",			";Range of times to process",
  "seed=",			";Generator seed for random values",
  "VERSION=2.2",                ";Josh Barnes  26 July 2017",
  NULL,
};

void snaptrak(bool usemass);			// track groups of bodies
local stream execmap(string);			// start snapmap process

#define GroupField  phatbody[NewBodyFields+0]
#define Group(b)    SelectInt(b, GroupField.offset)
#define GroupTag    "Group"

string names[2] = { GroupTag,  NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { IntType,  NULL };

string btags[] = { PosTag, VelTag, KeyTag, MassTag, GroupTag, NULL};
string otags[] = { PosTag, VelTag, KeyTag, MassTag, NULL};

bodyptr bodytab = NULL, traktab = NULL;
int nbody = 0, ntrak = 0;
real tbody;

int main(int argc, string argv[])
{
  string prog, itags[MaxBodyFields];
  stream xstr, ostr;
  int nold = -1;
  bool usemass;

  initparam(argv, defv);
  exprs[0] = getparam("group");
  prog = mktemp((string) copxstr("/tmp/sm_XXXXXX", sizeof(char)));
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog);
  if (get_tag_ok(xstr, "History"))
    skip_item(xstr);
  get_history(xstr);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  new_field(&GroupField, IntType, GroupTag);
  new_field(&GroupField + 1, NULL, NULL);
  layout_body(btags, Precision, NDIM);
  while (get_snap(xstr, &bodytab, &nbody, &tbody, itags, FALSE, NULL)) {
    if (nold == -1) {				// reading first frame?
      usemass = set_member(itags, MassTag);	// remember if mass present
      if (! usemass) {				// if no masses available
	otags[3] = NULL;			// don't output mass data
	eprintf("[%s: warning: assuming equal masses]\n", getprog());
      }
    }
    snaptrak(usemass);
    put_snap(ostr, &traktab, &ntrak, &tbody, otags);
    if (ntrak != nold)
      eprintf("[%s: wrote %d groups at t = %f]\n",
	      getprog(), ntrak, tbody);
    nold = ntrak;
  }
  strclose(ostr);
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getprog(), prog);
  return 0;
}

void snaptrak(bool usemass)
{
  int nzero = 0;
  
  if (traktab == NULL) {
    ntrak = 0;
    for (bodyptr bp = bodytab; bp < NthBody(bodytab,nbody); bp = NextBody(bp))
      ntrak = MAX(ntrak, Group(bp));
    eprintf("[%s: allocating %d groups]\n", getprog(), ntrak);
    traktab = (bodyptr) allocate(ntrak * SizeofBody);
  }
  for (bodyptr gp = traktab; gp < NthBody(traktab,ntrak); gp = NextBody(gp)) {
    Mass(gp) = 0.0;
    CLRV(Pos(gp));
    CLRV(Vel(gp));
    Key(gp) = 0;
  }
  for (bodyptr bp = bodytab; bp < NthBody(bodytab,nbody); bp = NextBody(bp)) {
    if (Group(bp) > 0) {
      if (Group(bp) > ntrak)
	error("%s: cant expand group array\n", getprog());
      gp = NthBody(traktab, Group(bp) - 1);
      ADDMULVS(Pos(gp), Pos(bp), (usemass ? Mass(bp) : 1.0));
      ADDMULVS(Vel(gp), Vel(bp), (usemass ? Mass(bp) : 1.0));
      Key(gp)++;
      Mass(gp) += (usemass ? Mass(bp) : 1.0);
    }
  }
  for (bodyptr gp = traktab; gp < NthBody(traktab,ntrak); gp = NextBody(gp))
    if (Mass(gp) != 0.0) {
      DIVVS(Pos(gp), Pos(gp), Mass(gp));
      DIVVS(Vel(gp), Vel(gp), Mass(gp));
    } else
      nzero++;
  if (nzero > 0)
    eprintf("[%s: %d groups have zero mass]\n", getprog(), nzero);
}

//  execmap: start snapmap subprocess, and return snapmap output stream.
//  ____________________________________________________________________

stream execmap(string prog)
{
  int handle[2];
  char handbuf[32];

  pipe(handle);
  if (fork() == 0) {                            // if this is child process
    close(handle[0]);
    sprintf(handbuf, "-%d", handle[1]);
    execl(prog, getprog(), getparam("in"), handbuf, getparam("times"),
	  PosTag "," VelTag,
	  PosTag "," VelTag "," GroupTag,
	  "true", getparam("seed"), NULL);
    error("%s: execl %s failed\n", getprog(), prog);
  }
  close(handle[1]);
  sprintf(handbuf, "-%d", handle[0]);
  return (stropen(handbuf, "r"));
}
