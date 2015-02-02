/*
 * SNAPSTACK: combine two snapshot frames.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

#include <string.h>

string defv[] = {		";Stack two N-body systems",
  "in1=???",			";1st N-body snapshot input file",
  "in2=???",			";2nd N-body snapshot input file",
  "out=???",			";N-body snapshot output file",
  "deltar=0.0,0.0,0.0",		";Position offset vector",
  "deltav=0.0,0.0,0.0",		";Velocity offset vector",
  "produce=*",			";List of output items",
  "VERSION=2.0",		";Josh Barnes  26 November 1997",
  NULL,
};

void checktags(string *, string *);
void snapstack(bodyptr, bodyptr, int, bodyptr, int, string *);
void setvect(vector, string *);

int main(int argc, string argv[])
{
  stream istr, ostr;
  string itags1[MaxBodyFields], itags2[MaxBodyFields];
  string *produce, *otags;
  bodyptr btab1 = NULL, btab2 = NULL, btab;
  int nbody1, nbody2, nbody;
  real tsnap;

  initparam(argv, defv);
  if (! streq(getparam("produce"), "*")) {
    produce = burststring(getparam("produce"), ", ");
    layout_body(produce, Precision, NDIM);
  } else
    produce = NULL;
  istr = stropen(getparam("in1"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab1, &nbody1, &tsnap, itags1, produce == NULL))
    error("%s: no data in 1st input file\n");
  istr = stropen(getparam("in2"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab2, &nbody2, &tsnap, itags2, FALSE))
    error("%s: no data in 2nd input file\n");
  otags = set_inter(itags1, itags2);
  if (produce != NULL)
    checktags(otags, produce);
  eprintf("[%s: nbody = %d + %d]\n", getargv0(), nbody1, nbody2);
  nbody = nbody1 + nbody2;
  btab = (bodyptr) allocate(SizeofBody * nbody);
  snapstack(btab, btab1, nbody1, btab2, nbody2, otags);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  put_snap(ostr, &btab, &nbody, &tsnap, otags);
  strclose(ostr);
  return (0);
}

void checktags(string *otags, string *produce)
{
  string *missing, *mp;
  char buf[512];

  missing = set_diff(produce, otags);
  if (set_length(missing) > 0) {
    buf[0] = (char) NULL;
    for (mp = missing; *mp != NULL; mp++) {
      if (mp != missing)
	strcat(buf, ",");
      strcat(buf, *mp);
    }
    error("%s: missing %s data\n", getargv0(), buf);
  }
  free(missing);
}

void snapstack(bodyptr btab, bodyptr bt1, int nb1, bodyptr bt2, int nb2,
	       string *tags)
{
  vector deltar, deltav;
  bodyptr bp;
  int i;
  
  setvect(deltar, burststring(getparam("deltar"), ", "));
  setvect(deltav, burststring(getparam("deltav"), ", "));
  for (i = 0; i < nb1; i++) {
    bp = NthBody(btab, i);
    memcpy(bp, NthBody(bt1, i), SizeofBody);
    if (set_member(tags, PosTag)) {
      ADDMULVS(Pos(bp), deltar, 0.5);
    }
    if (set_member(tags, VelTag)) {
      ADDMULVS(Vel(bp), deltav, 0.5);
    }
  }
  for (i = 0; i < nb2; i++) {
    bp = NthBody(btab, i + nb1);
    memcpy(bp, NthBody(bt2, i), SizeofBody);
    if (set_member(tags, PosTag)) {
      ADDMULVS(Pos(bp), deltar, -0.5);
    }
    if (set_member(tags, VelTag)) {
      ADDMULVS(Vel(bp), deltav, -0.5);
    }
  }
}

void setvect(vector vect, string *values)
{
  int i;

  for (i = 0; i < NDIM; i++) {
    if (values[i] == NULL)
      error("%s: not enough values\n", getargv0());
    vect[i] = atof(values[i]);
  }
}
