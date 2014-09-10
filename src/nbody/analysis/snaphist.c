/*
 * snaphist.c: make histogram of snapshot data.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "getparam.h"
#include "filestruct.h"
#include "phatbody.h"
#include "buildmap.h"
#include <unistd.h>

string defv[] = {		";Construct histogram of particle data",
  "in=???",			";Input snapshot file name",
  "times=all",			";Times to combine in analysis",
  "value=???",			";C language expression for value to bin.",
			        ";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "weight=",			";C language expression for weight factor.",
				";If blank, do unweighted enumeration.",
  "require=",			";Input items required for computation",
  "nbins=32",			";Number of bins in histogram",
  "range=0.0:1.0",		";Range of values binned in histogram",
  "plotstyle=false",		";If true, format output for plotting",
  "seed=",			";Seed for random number generator",
  "VERSION=2.3",		";Josh Barnes  9 Sep 2014",
  NULL,
};

void loadihist(int *, int, real *, bodyptr, int);
void listihist(int *, int, real *, bool);
void loadrhist(real *, int, real *, bodyptr, int);
void listrhist(real *, int, real *, bool);
stream execmap(string, bool);

string names[3] = { "Value",  "Weight",  NULL };
string exprs[3] = { NULL,     NULL,      NULL };
string types[3] = { RealType, RealType,  NULL };

#define ValueField  phatbody[NewBodyFields+0]
#define WeightField phatbody[NewBodyFields+1]
#define Value(b)    SelectReal(b, ValueField.offset)
#define Weight(b)   SelectReal(b, WeightField.offset)

int main(int argc, string argv[])
{
  string prog, itags[MaxBodyFields];
  stream xstr;
  int nbins, *ibins, nbody;
  real range[2], *rbins, tnow;
  bodyptr btab = NULL;

  initparam(argv, defv);
  exprs[0] = getparam("value");
  if (! strnull(getparam("weight")))
    exprs[1] = getparam("weight");
  else
    names[1] = exprs[1] = types[1] = NULL;
  prog = tempnam("/tmp", "sm");
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog, exprs[1] != NULL);
  get_history(xstr);
  if (exprs[1] == NULL) {
    new_field(&ValueField, RealType, "Value");
    new_field(&ValueField + 1, NULL, NULL);
  } else {
    new_field(&ValueField, RealType, "Value");
    new_field(&WeightField, RealType, "Weight");
    new_field(&WeightField + 1, NULL, NULL);
  }
  layout_body(names, Precision, NDIM);
  nbins = getiparam("nbins");
  setrange(range, getparam("range"));
  if (exprs[1] == NULL) {
    ibins = (int *) allocate((nbins + 2) * sizeof(int));
    while (get_snap(xstr, &btab, &nbody, &tnow, itags, FALSE))
      loadihist(ibins, nbins, range, btab, nbody);
    listihist(ibins, nbins, range, getbparam("plotstyle"));
  } else {
    rbins = (real *) allocate((nbins + 2) * sizeof(real));
    while (get_snap(xstr, &btab, &nbody, &tnow, itags, FALSE))
      loadrhist(rbins, nbins, range, btab, nbody);
    listrhist(rbins, nbins, range, getbparam("plotstyle"));
  }
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getargv0(), prog);
  return (0);
}

void loadihist(int *bins, int nbins, real *range, bodyptr btab, int nbody)
{
  bodyptr bp;
  int i;

  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    i = 1 + nbins * (Value(bp) - range[0]) / (range[1] - range[0]);
    i = MAX(i, 0);
    i = MIN(i, nbins + 1);
    bins[i]++;
  }
}

void listihist(int *bins, int nbins, real *range, bool plotstyle)
{
  real binwidth;
  int i;

  binwidth = (range[1] - range[0]) / nbins;
  printf("#  count_low: %d  count_high: %d\n", bins[0], bins[nbins+1]);
  if (! plotstyle) {
    printf("#%11s %11s %11s\n", "min", "max", "count");
    for (i = 1; i <= nbins; i++)
      printf(" %#11.5g %#11.5g %11d\n",
	     range[0] + binwidth * (i-1), range[0] + binwidth * i,
	     bins[i]);
  } else {
    printf("#%11s %11s\n", "value", "count");
    for (i = 1; i <= nbins; i++) {
      printf(" %#11.5g %11d\n", range[0] + binwidth * (i-1), bins[i]);
      printf(" %#11.5g %11d\n", range[0] + binwidth * i, bins[i]);
    }
  }
}

void loadrhist(real *bins, int nbins, real *range, bodyptr btab, int nbody)
{
  bodyptr bp;
  int i;
  
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    i = 1 + nbins * (Value(bp) - range[0]) / (range[1] - range[0]);
    i = MAX(i, 0);
    i = MIN(i, nbins + 1);
    bins[i] += Weight(bp);
  }
}

void listrhist(real *bins, int nbins, real *range, bool plotstyle)
{
  real binwidth;
  int i;

  binwidth = (range[1] - range[0]) / nbins;
  printf("#  sum_low: %g  sum_high: %g\n", bins[0], bins[nbins+1]);
  if (! plotstyle) {
    printf("#%11s %11s %11s\n", "min", "max", "sum");
    for (i = 1; i <= nbins; i++)
      printf(" %#11.5g %#11.5g %11.5g\n",
	     range[0] + binwidth * (i-1), range[0] + binwidth * i,
	     bins[i]);
  } else {
    printf("#%11s %11s\n", "value", "count");
    for (i = 1; i <= nbins; i++) {
      printf(" %#11.5g %11.5g\n", range[0] + binwidth * (i-1), bins[i]);
      printf(" %#11.5g %11.5g\n", range[0] + binwidth * i,     bins[i]);
    }
  }
}

//  execmap: start snapmap subprocess, and return snapmap output stream.
//  ____________________________________________________________________

stream execmap(string prog, bool weightflag)
{
  int handle[2];
  char handbuf[32];

  pipe(handle);
  if (fork() == 0) {                            // if this is child process
    close(handle[0]);
    sprintf(handbuf, "-%d", handle[1]);
    execl(prog, getargv0(), getparam("in"), handbuf, getparam("times"),
	  getparam("require"), weightflag ? "Value,Weight" : "Value",
	  strnull(getparam("require")) ? "true" : "false",
	  getparam("seed"), NULL);
    error("%s: execl %s failed\n", getargv0(), prog);
  }
  close(handle[1]);
  sprintf(handbuf, "-%d", handle[0]);
  return (stropen(handbuf, "r"));
}
