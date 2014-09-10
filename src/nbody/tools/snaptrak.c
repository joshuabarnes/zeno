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
  "group=???",			";C language expression for group membership",
				";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "times=all",			";Range of times to process",
  "seed=",			";Seed for random number generator",
  "VERSION=2.1",                ";Josh Barnes  9 Sep 2014",
  NULL,
};

local void snaptrak(void);			// track groups of bodies
local stream execmap(string);			// start snapmap process

string names[2] = { "Group",  NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { IntType,  NULL };

#define GroupField  phatbody[NewBodyFields+0]
#define Group(b)  SelectInt(b, GroupField.offset)

string btags[] = { MassTag, PosTag, VelTag, KeyTag, "Group", NULL};
string otags[] = { MassTag, PosTag, VelTag, KeyTag, NULL};

bodyptr bodytab = NULL, traktab = NULL;
int nbody, ntrak;
real tbody;

int main(int argc, string argv[])
{
  string prog, itags[MaxBodyFields];
  stream xstr, ostr;
  int nold = -1;

  initparam(argv, defv);
  exprs[0] = getparam("group");
  prog = tempnam("/tmp", "sm");
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog);
  if (get_tag_ok(xstr, "History"))
    skip_item(xstr);
  get_history(xstr);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  new_field(&GroupField, IntType, "Group");
  new_field(&GroupField + 1, NULL, NULL);
  layout_body(btags, Precision, NDIM);
  while (get_snap(xstr, &bodytab, &nbody, &tbody, itags, FALSE)) {
    snaptrak();
    put_snap(ostr, &traktab, &ntrak, &tbody, otags);
    if (ntrak != nold)
      eprintf("[%s: wrote %d groups at t = %f]\n",
	      getprog(), ntrak, tbody);
    nold = ntrak;
  }
  strclose(ostr);
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getprog(), prog);
  return (0);
}

void snaptrak(void)
{
  bodyptr bp, gp;
  int nzero;
  
  if (traktab == NULL) {
    ntrak = 0;
    for (bp = bodytab; bp < NthBody(bodytab, nbody); bp = NextBody(bp))
      ntrak = MAX(ntrak, Group(bp));
    eprintf("[%s: allocating %d groups]\n", getprog(), ntrak);
    traktab = (bodyptr) allocate(ntrak * SizeofBody);
  }
  for (gp = traktab; gp < NthBody(traktab, ntrak); gp = NextBody(gp)) {
    Mass(gp) = 0.0;
    CLRV(Pos(gp));
    CLRV(Vel(gp));
    Key(gp) = 0;
  }
  for (bp = bodytab; bp < NthBody(bodytab, nbody); bp = NextBody(bp)) {
    if (Group(bp) > ntrak)
      error("snaptrak: cant expand group array\n");
    if (Group(bp) > 0) {
      gp = NthBody(traktab, Group(bp) - 1);
      Mass(gp) += Mass(bp);
      ADDMULVS(Pos(gp), Pos(bp), Mass(bp));
      ADDMULVS(Vel(gp), Vel(bp), Mass(bp));
      Key(gp)++;
    }
  }
  nzero = 0;
  for (gp = traktab; gp < NthBody(traktab, ntrak); gp = NextBody(gp))
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
	  MassTag "," PosTag "," VelTag,
	  MassTag "," PosTag "," VelTag ",Group",
	  "true", getparam("seed"), NULL);
    error("%s: execl %s failed\n", getprog(), prog);
  }
  close(handle[1]);
  sprintf(handbuf, "-%d", handle[0]);
  return (stropen(handbuf, "r"));
}
