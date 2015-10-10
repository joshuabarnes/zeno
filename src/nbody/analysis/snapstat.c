/*
 * snapstat.c: compute statistics of snapshot data.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "getparam.h"
#include "filestruct.h"
#include "phatbody.h"
#include "buildmap.h"
#include <unistd.h>

string defv[] = {		";Compute statistics of snapshot data",
  "in=???",			";Input file: one or more snapshots",
  "times=all",			";Range of times to analyze",
  "value=???",			";C language expression for value analyzed.",
			        ";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "require=",			";Input data items required",
  "options=avg,med",		";Statistical measures to tabulate.",
				";Possible choices are one or more of:",
				";  avg: list average, rms, etc.",
				";  sum: list sum, sum of squares, etc.",
				";  med: list median, quartiles, & limits.",
				";  oct: list 1st through 7th octiles.",
				";  OCT: list octiles & limits (wide!).",
				";  time: prepend current time to line.",
  "seed=",			";Random number seed for value computation",
  "formats=%11s,%#11.5g,%11d,%9s,%#9.3g",
				";Format specifications for output table",
  "output=",			";Output file; uses stdout if blank",
  "stream=",			";Output stream copies input snapshots",
  "VERSION=2.5",		";Josh Barnes  26 June 2015",
  NULL,
};

void snapstat(stream, bodyptr, int, real, string, string*);
void snapavg(stream, bodyptr, int, real, bool, bool, string*);
void snapsum(stream, bodyptr, int, real, bool, bool, string*);
void snapmed(stream, bodyptr, int, real, bool, bool, string*);
void snapoct(stream, bodyptr, int, real, bool, bool, string*);
void snapOCT(stream, bodyptr, int, real, bool, bool, string*);

int cmpvalue(const void *, const void *);

stream execmap(string);				// start snapmap process

#define ValueTag    "Value"
#define ValueField  phatbody[NewBodyFields+0]
#define Value(b)    SelectReal(b, ValueField.offset)

string names[2] = { ValueTag, NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { RealType, NULL };

int main(int argc, string argv[])
{
  string  *formats, prog, itags[MaxBodyFields], *otags;
  stream xstr, outstr, strstr;
  bodyptr btab = NULL;
  int nbody;
  real tnow;

  initparam(argv, defv);
  formats = burststring(getparam("formats"), ",");
  if (xstrlen(formats, sizeof(string)) != 6)
    error("%s: formats arg needs 5 specifications\n", getprog());
  exprs[0] = getparam("value");
  prog = mktemp((string) copxstr("/tmp/sm_XXXXXX", sizeof(char)));
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog);
  get_history(xstr);
  new_field(&ValueField, RealType, ValueTag);
  new_field(&ValueField + 1, NULL, NULL);
  layout_body(names, Precision, NDIM);
  outstr = (streq(getparam("output"), "") ?
	    stdout : stropen(getparam("output"), "w"));
  strstr = (streq(getparam("stream"), "") ?
	    NULL : stropen(getparam("stream"), "w"));
  while (get_snap(xstr, &btab, &nbody, &tnow, itags, strstr != NULL)) {
    snapstat(outstr, btab, nbody, tnow, getparam("options"), formats);
    if (strstr != NULL) {
      otags = set_diff(itags, names);
      put_snap(strstr, &btab, &nbody, &tnow, otags);
      free(otags);
    }
    fflush(NULL);
  }
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getargv0(), prog);
  return (0);
}

void snapstat(stream outstr, bodyptr btab, int nbody, real tnow,
	      string opts, string *fmts)
{
  static bool showhead = TRUE;
  bool showtime = scanopt(opts, "time");
  int nopt;

  nopt = 0;
  if (scanopt(opts, "avg")) {
    snapavg(outstr, btab, nbody, tnow, showtime, showhead, fmts);
    nopt++;
  }
  if (scanopt(opts, "sum")) {
    snapsum(outstr, btab, nbody, tnow, showtime, showhead, fmts);
    nopt++;
  }
  if (scanopt(opts, "med")) {
    snapmed(outstr, btab, nbody, tnow, showtime, showhead, fmts);
    nopt++;
  }
  if (scanopt(opts, "oct")) {
    snapoct(outstr, btab, nbody, tnow, showtime, showhead, fmts);
    nopt++;
  }
  if (scanopt(opts, "OCT")) {
    snapOCT(outstr, btab, nbody, tnow, showtime, showhead, fmts);
    nopt++;
  }
  if (nopt == 0)
    error("%s: options %s not recognized\n", getargv0(), opts);
  showhead = (nopt > 1);
}

void snapavg(stream outstr, bodyptr btab, int nbody, real tnow,
	     bool showtime, bool showhead, string *fmts)
{
  double avg1, avg2, avg3, avg4;
  bodyptr bp;
  char fmtbuf[512];

  avg1 = avg2 = avg3 = avg4 = 0.0;
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    avg1 += Value(bp) / nbody;
    avg2 += rsqr(Value(bp)) / nbody;
    avg3 += rqbe(Value(bp)) / nbody;
    avg4 += rsqr(rsqr(Value(bp))) / nbody;
  }
  if (showtime) {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s %s\n",
	       fmts[0], fmts[0], fmts[0], fmts[0], fmts[0], fmts[0]);
      fprintf(outstr, fmtbuf, "time",
	      "values", "average", "r.m.s.", "avg^3", "avg^4");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s %s\n",
	     fmts[1], fmts[2], fmts[1], fmts[1], fmts[1], fmts[1]);
    fprintf(outstr, fmtbuf, tnow,
	    nbody, avg1, rsqrt(avg2), rcbrt(avg3), rsqrt(rsqrt(avg4)));
  } else {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s\n",
	       fmts[0], fmts[0], fmts[0], fmts[0], fmts[0]);
      fprintf(outstr, fmtbuf,
	      "values", "average", "r.m.s.", "avg^3", "avg^4");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s\n",
	     fmts[2], fmts[1], fmts[1], fmts[1], fmts[1]);
    fprintf(outstr, fmtbuf,
	    nbody, avg1, rsqrt(avg2), rcbrt(avg3), rsqrt(rsqrt(avg4)));
  }
}

void snapsum(stream outstr, bodyptr btab, int nbody, real tnow,
	     bool showtime, bool showhead, string *fmts)
{
  double sum1, sum2, sum3, sum4;
  bodyptr bp;
  char fmtbuf[512];

  sum1 = sum2 = sum3 = sum4 = 0.0;
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    sum1 += Value(bp);
    sum2 += rsqr(Value(bp));
    sum3 += rqbe(Value(bp));
    sum4 += rsqr(rsqr(Value(bp)));
  }
  if (showtime) {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s %s\n",
	       fmts[0], fmts[0], fmts[0], fmts[0], fmts[0], fmts[0]);
      fprintf(outstr, fmtbuf, "time",
	      "values", "sum", "sum^2", "sum^3", "sum^4");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s %s\n",
	     fmts[1], fmts[2], fmts[1], fmts[1], fmts[1], fmts[1]);
    fprintf(outstr, fmtbuf, tnow, nbody, sum1, sum2, sum3, sum4);
  } else {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s\n",
	       fmts[0], fmts[0], fmts[0], fmts[0], fmts[0]);
      fprintf(outstr, fmtbuf, "values", "sum", "sum^2", "sum^3", "sum^4");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s\n",
	     fmts[2], fmts[1], fmts[1], fmts[1], fmts[1]);
    fprintf(outstr, fmtbuf, nbody, sum1, sum2, sum3, sum4);
  }
}

void snapmed(stream outstr, bodyptr btab, int nbody, real tnow,
	     bool showtime, bool showhead, string *fmts)
{
  real *valarr = (real*) allocate(nbody * sizeof(real));
  char fmtbuf[512];

  for (int i = 0; i < nbody; i++)
    valarr[i] = Value(NthBody(btab, i));
  qsort(valarr, nbody, sizeof(real), cmpvalue);
  if (showtime) {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s %s\n",
	       fmts[0], fmts[0], fmts[0], fmts[0], fmts[0], fmts[0]);
      fprintf(outstr, fmtbuf, "time",
	      "minimum", "1st quart", "median", "3rd quart", "maximum");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s %s\n",
	     fmts[1], fmts[1], fmts[1], fmts[1], fmts[1], fmts[1]);
    fprintf(outstr, fmtbuf, tnow,
	    valarr[0], 
	    valarr[nbody/4],
	    valarr[nbody/2],
	    valarr[3*nbody/4],
	    valarr[nbody - 1]);
  } else {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s\n",
	       fmts[0], fmts[0], fmts[0], fmts[0], fmts[0]);
      fprintf(outstr, fmtbuf,
	      "minimum", "1st quart", "median", "3rd quart", "maximum");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s\n",
	     fmts[1], fmts[1], fmts[1], fmts[1], fmts[1]);
    fprintf(outstr, fmtbuf,
	    valarr[0], 
	    valarr[nbody/4],
	    valarr[nbody/2],
	    valarr[3*nbody/4],
	    valarr[nbody - 1]);
  }
  free(valarr);
}

void snapoct(stream outstr, bodyptr btab, int nbody, real tnow,
	     bool showtime, bool showhead, string *fmts)
{
  real *valarr = (real*) allocate(nbody * sizeof(real));
  char fmtbuf[512];

  for (int i = 0; i < nbody; i++)
    valarr[i] = Value(NthBody(btab, i));
  qsort(valarr, nbody, sizeof(real), cmpvalue);
  if (showtime) {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s %s %s %s\n", fmts[3],
	       fmts[3], fmts[3], fmts[3], fmts[3], fmts[3], fmts[3], fmts[3]);
      fprintf(outstr, fmtbuf, "time",
	      "oct1", "oct2", "oct3", "oct4", "oct5", "oct6", "oct7");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s %s %s %s\n", fmts[4],
	     fmts[4], fmts[4], fmts[4], fmts[4], fmts[4], fmts[4], fmts[4]);
    fprintf(outstr, fmtbuf, tnow,
	    valarr[(1*nbody)/8],
	    valarr[(2*nbody)/8],
	    valarr[(3*nbody)/8],
	    valarr[(4*nbody)/8],
	    valarr[(5*nbody)/8],
	    valarr[(6*nbody)/8],
	    valarr[(7*nbody)/8]);
  } else {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s %s %s\n",
	       fmts[3], fmts[3], fmts[3], fmts[3], fmts[3], fmts[3], fmts[3]);
      fprintf(outstr, fmtbuf,
	      "oct1", "oct2", "oct3", "oct4", "oct5", "oct6", "oct7");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s %s %s\n",
	     fmts[4], fmts[4], fmts[4], fmts[4], fmts[4], fmts[4], fmts[4]);
    fprintf(outstr, fmtbuf,
	    valarr[(1*nbody)/8],
	    valarr[(2*nbody)/8],
	    valarr[(3*nbody)/8],
	    valarr[(4*nbody)/8],
	    valarr[(5*nbody)/8],
	    valarr[(6*nbody)/8],
	    valarr[(7*nbody)/8]);
  }
  free(valarr);
}

void snapOCT(stream outstr, bodyptr btab, int nbody, real tnow,
	     bool showtime, bool showhead, string *fmts)
{
  real *valarr = (real*) allocate(nbody * sizeof(real));
  char fmtbuf[512];

  for (int i = 0; i < nbody; i++)
    valarr[i] = Value(NthBody(btab, i));
  qsort(valarr, nbody, sizeof(real), cmpvalue);
  if (showtime) {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s %s %s %s %s %s\n",
	       fmts[3], fmts[3], fmts[3], fmts[3], fmts[3], fmts[3], fmts[3],
	       fmts[3], fmts[3], fmts[3]);
      fprintf(outstr, fmtbuf, "time", "min", "oct1", "oct2", "oct3", "oct4",
	      "oct5", "oct6", "oct7", "max");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s %s %s %s %s %s\n",
	     fmts[4], fmts[4], fmts[4], fmts[4], fmts[4], fmts[4], fmts[4],
	     fmts[4], fmts[4], fmts[4]);
    fprintf(outstr, fmtbuf, tnow,
	    valarr[0],
	    valarr[(1*nbody)/8],
	    valarr[(2*nbody)/8],
	    valarr[(3*nbody)/8],
	    valarr[(4*nbody)/8],
	    valarr[(5*nbody)/8],
	    valarr[(6*nbody)/8],
	    valarr[(7*nbody)/8],
	    valarr[nbody - 1]);
  } else {
    if (showhead) {
      snprintf(fmtbuf, sizeof(fmtbuf), "#%s %s %s %s %s %s %s %s %s\n",
	       fmts[3], fmts[3], fmts[3], fmts[3], fmts[3], fmts[3],
	       fmts[3], fmts[3], fmts[3]);
      fprintf(outstr, fmtbuf, "min", "oct1", "oct2", "oct3", "oct4", "oct5",
	      "oct6", "oct7", "max");
    }
    snprintf(fmtbuf, sizeof(fmtbuf), " %s %s %s %s %s %s %s %s %s\n",
	     fmts[4], fmts[4], fmts[4], fmts[4], fmts[4], fmts[4],
	     fmts[4], fmts[4], fmts[4]);
    fprintf(outstr, fmtbuf,
	    valarr[0],
	    valarr[(1*nbody)/8],
	    valarr[(2*nbody)/8],
	    valarr[(3*nbody)/8],
	    valarr[(4*nbody)/8],
	    valarr[(5*nbody)/8],
	    valarr[(6*nbody)/8],
	    valarr[(7*nbody)/8],
	    valarr[nbody - 1]);
  }
  free(valarr);
}

int cmpvalue(const void *a, const void *b)
{
  real aval = *((real*) a), bval = *((real*) b);
  return (aval < bval ? -1 : aval > bval ? 1 : 0);
}

stream execmap(string prog)
{
  int handle[2];
  char handbuf[32];

  pipe(handle);
  if (fork() == 0) {                           // if this is child process
    close(handle[0]);
    sprintf(handbuf, "-%d", handle[1]);
    execl(prog, getargv0(), getparam("in"), handbuf, getparam("times"),
	  getparam("require"), ValueTag,
	  strnull(getparam("require")) ? "true" : "false",
	  getparam("seed"), NULL);
    error("%s: execl %s failed\n", getargv0(), prog);
  }
  close(handle[1]);
  sprintf(handbuf, "-%d", handle[0]);
  return (stropen(handbuf, "r"));
}
