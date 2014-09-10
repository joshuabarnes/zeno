/*
 * snapcons.c: construct compound n-body system.
 */

#include "stdinc.h"
#include "strset.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include <string.h>

string defv[] = {		";Construct compound N-body system",
  "in=???",			";Input file with series of models",
  "out=???",			";Output file for composite model",
  "frame=",			";Input file for CM coordinates.",
				";If blank, use zero offsets.",
  "produce=" PosTag "," VelTag "," MassTag,
				";List of output data items",
  "VERSION=2.2",		";Josh Barnes  5 January 2013",
  NULL,
};

typedef struct _snapshot {
  bodyptr bodies;
  int nbody;
  real time;
  struct _snapshot *link;
} snapshot;

#define get_snapshot(str, snap, tags, ext)				\
  get_snap(str, &(snap).bodies, &(snap).nbody, &(snap).time, tags, ext)

int main(int argc, string argv[])
{
  string *produce, snap_tags[MaxBodyFields], frame_tags[MaxBodyFields];
  snapshot *snaplist = NULL, *ssp, frame = { NULL, 0, 0.0, NULL };
  stream instr, outstr;
  int n, nbody, i;
  bodyptr btab, bp, sp;
  bool roff, voff;

  initparam(argv, defv);
  produce = burststring(getparam("produce"), ", ");
  layout_body(produce, Precision, NDIM);
  instr = stropen(getparam("in"), "r");
  get_history(instr);
  snaplist = ssp = (snapshot *) allocate(sizeof(snapshot));
  n = nbody = 0;
  while (get_snapshot(instr, *ssp, snap_tags, FALSE)) {
    if (! set_equal(snap_tags, produce))
      error("%s: input file lacks required items\n", getprog());
    n++;
    nbody += ssp->nbody;
    get_history(instr);
    ssp->link = (snapshot *) allocate(sizeof(snapshot));
    ssp = ssp->link;
  }
  strclose(instr);
  eprintf("[%s: read %d frames, %d bodies]\n", getprog(), n, nbody);
  if (! strnull(getparam("frame"))) {
    instr = stropen(getparam("frame"), "r");
    get_history(instr);
    if (! get_snapshot(instr, frame, frame_tags, FALSE))
      error("%s: no snapshot in frame file\n", getprog());
    strclose(instr);
    if (frame.nbody < n)
      error("%s: need at least %d points in frame\n", getprog(), n);
  }
  roff = (frame.nbody > 0) && set_member(frame_tags, PosTag);
  voff = (frame.nbody > 0) && set_member(frame_tags, VelTag);
  bp = btab = (bodyptr) allocate(nbody * SizeofBody);
  ssp = snaplist;
  for (i = 0; i < n; i++) {
    sp = ssp->bodies;
    while (sp < NthBody(ssp->bodies, ssp->nbody)) {
      memmove((void *) bp, (void *) sp, (size_t) SizeofBody);
      if (roff) {
	ADDV(Pos(bp), Pos(bp), Pos(NthBody(frame.bodies, i)));
      }
      if (voff) {
	ADDV(Vel(bp), Vel(bp), Vel(NthBody(frame.bodies, i)));
      }
      bp = NextBody(bp);
      sp = NextBody(sp);
    }
    ssp = ssp->link;
  }
  outstr = stropen(getparam("out"), "w");
  put_history(outstr);
  put_snap(outstr, &btab, &nbody, &frame.time, produce);
  strclose(outstr);
  return (0);
}
