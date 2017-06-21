/*
 * snapavg.c: calculate averages or standard deviations over time.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "getparam.h"
#include "filestruct.h"
#include "phatbody.h"
#include "strset.h"
#include "buildmap.h"
#include <unistd.h>

string defv[] = {		";Find averages or std. devs. over time.",
				";Results are output using Aux field.",
  "in=???",			";Input file name",
  "out=???",			";Output file name",
  "times=all",			";Times to analyze",
  "value=???",			";C language expression for value to average.",
			        ";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "stddev=false",		";Find standard deviation instead of average",
  "require=",			";Input items required",
  "produce=",                   ";Output items produced",
  "passall=true",               ";If true, pass on input data",
  "seed=",			";Seed for random number generator",
  "VERSION=1.1",		";Josh Barnes  24 May 2015",
  NULL,
};

stream execmap(string);				// start snapmap process 

string names[2] = { "Value",  NULL };
string exprs[2] = { NULL,      NULL };
string types[2] = { RealType,  NULL };

#define ValueField  phatbody[NewBodyFields+0]
#define Value(b)  SelectReal(b, ValueField.offset)

string btags[] = { AuxTag, "Value", NULL };

int main(int argc, string argv[])
{
  string prog, itags[MaxBodyFields], *otags = NULL;
  stream xstr, ostr;
  bodyptr btab = NULL;
  int nbody, nsnap, i;
  real tnow;
  double *sum = NULL, *sum2 = NULL, std2;

  initparam(argv, defv);
  exprs[0] = getparam("value");
  prog = mktemp((string) copxstr("/tmp/sm_XXXXXX", sizeof(char)));
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog);
  if (get_tag_ok(xstr, "History"))
    skip_item(xstr);
  get_history(xstr);
  new_field(&ValueField, RealType, "Value");
  new_field(&ValueField + 1, NULL, NULL);
  layout_body(btags, Precision, NDIM);
  nsnap = 0;
  while (get_snap(xstr, &btab, &nbody, &tnow, itags, TRUE, NULL)) {
    if (otags == NULL)
      otags = set_diff(set_union(itags, set_cons(AuxTag, NULL)), names);
    if (sum == NULL)
      sum = (double *) allocate(nbody * sizeof(double));
    for (i = 0; i < nbody; i++)
      sum[i] += Value(NthBody(btab, i));
    if (getbparam("stddev")) {
      if (sum2 == NULL)
	sum2 = (double *) allocate(nbody * sizeof(double));
      for (i = 0; i < nbody; i++)
	sum2[i] += rsqr(Value(NthBody(btab, i)));
    }
    nsnap++;
  }
  if (nsnap == 0)
    error("%s: no frames in input\n", getargv0());
  eprintf("[%s: averaging %d frames]\n", getargv0(), nsnap);
  if (! getbparam("stddev")) {
    for (i = 0; i < nbody; i++)
      Aux(NthBody(btab, i)) = sum[i] / nsnap;
  } else {
    for (i = 0; i < nbody; i++) {
      std2 = sum2[i] / nsnap - (sum[i] / nsnap) * (sum[i] / nsnap);
      Aux(NthBody(btab, i)) = sqrt(MAX(std2, 0.0));
    }
  }
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getargv0(), prog);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  put_snap(ostr, &btab, &nbody, &tnow, otags);
  return (0);
}

//  execmap: start snapmap subprocess, and return snapmap output stream.
//  ____________________________________________________________________

stream execmap(string prog)
{
  int handle[2];
  char handbuf[32], produce[512];

  pipe(handle);
  if (fork() == 0) {                            // if this is child process
    close(handle[0]);
    sprintf(handbuf, "-%d", handle[1]);
    sprintf(produce, "%s,%s", getparam("produce"), names[0]);
    execl(prog, getargv0(), getparam("in"), handbuf, getparam("times"),
	  getparam("require"), produce, getparam("passall"), 
	  getparam("seed"), NULL);
    error("%s: execl %s failed\n", getargv0(), prog);
  }
  close(handle[1]);
  sprintf(handbuf, "-%d", handle[0]);
  return (stropen(handbuf, "r"));
}
