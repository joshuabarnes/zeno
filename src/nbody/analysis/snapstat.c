/*
 * snapstat.c: calculate statistics for input snapshots.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "getparam.h"
#include "filestruct.h"
#include "phatbody.h"
#include "buildmap.h"
#include <unistd.h>

string defv[] = {		";Analyze statistics of SnapShot data.",
				";Uses doubles for sums and averages.",
  "in=???",			";Input file name",
  "times=all",			";Times to analyze",
  "value=???",			";C language expression for analyzed value.",
			        ";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "require=",			";Input items required",
  "options=avg,med",		";Statistical measures to tabulate.",
				";Possible choices are one or more of:",
				";avg: list average, rms, etc.",
				";sum: list sum, sum of squares, etc.",
				";med: list median, quartiles, & limits.",
				";oct: list 1st through 7th octiles.",
				";OCT: list octiles & limits (wide!).",
				";time: prepend current time to line.",
  "seed=",			";Seed for random number generator",
  "VERSION=2.3",		";Josh Barnes  14 May 2008",
  NULL,
};

void snapstat(bodyptr, int, real, string);
void snapavg(bodyptr, int, real, bool, bool);
void snapsum(bodyptr, int, real, bool, bool);
void snapmed(bodyptr, int, real, bool, bool);
void snapoct(bodyptr, int, real, bool, bool);
void snapOCT(bodyptr, int, real, bool, bool);
int cmpvalue(const void *, const void *);
stream execmap(string);				// start snapmap process

string names[2] = { "Value",  NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { RealType, NULL };

#define ValueField  phatbody[NewBodyFields+0]
#define Value(b)  SelectReal(b, ValueField.offset)

int main(int argc, string argv[])
{
  string prog, itags[MaxBodyFields];
  stream xstr;
  bodyptr btab = NULL;
  int nbody;
  real tnow;

  initparam(argv, defv);
  exprs[0] = getparam("value");
  prog = tempnam("/tmp", "sm");
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog);
  get_history(xstr);
  new_field(&ValueField, RealType, "Value");
  new_field(&ValueField + 1, NULL, NULL);
  layout_body(names, Precision, NDIM);
  while (get_snap(xstr, &btab, &nbody, &tnow, itags, FALSE))
    snapstat(btab, nbody, tnow, getparam("options"));
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getargv0(), prog);
  return (0);
}

void snapstat(bodyptr btab, int nbody, real tnow, string options)
{
    static bool showhead = TRUE;
    bool showtime = scanopt(options, "time");
    int nopt;

    nopt = 0;
    if (scanopt(options, "avg")) {
	snapavg(btab, nbody, tnow, showtime, showhead);
	nopt++;
    }
    if (scanopt(options, "sum")) {
	snapsum(btab, nbody, tnow, showtime, showhead);
	nopt++;
    }
    if (scanopt(options, "med")) {
	snapmed(btab, nbody, tnow, showtime, showhead);
	nopt++;
    }
    if (scanopt(options, "oct")) {
	snapoct(btab, nbody, tnow, showtime, showhead);
	nopt++;
    }
    if (scanopt(options, "OCT")) {
	snapOCT(btab, nbody, tnow, showtime, showhead);
	nopt++;
    }
    if (nopt == 0)
      error("%s: options %s not recognized\n", getargv0(), options);
    showhead = (nopt > 1);
}

void snapavg(bodyptr btab, int nbody, real tnow, bool showtime, bool showhead)
{
    double avg1, avg2, avg3, avg4;
    bodyptr bp;

    avg1 = avg2 = avg3 = avg4 = 0.0;
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	avg1 += Value(bp) / nbody;
	avg2 += rsqr(Value(bp)) / nbody;
	avg3 += rqbe(Value(bp)) / nbody;
	avg4 += rsqr(rsqr(Value(bp))) / nbody;
    }
    if (showhead) {
	printf(showtime ? "#%11s " : "#", "time");
	printf("%11s %11s %11s %11s %11s\n",
	       "values", "average", "r.m.s.", "avg^3", "avg^4");
    }
    if (showtime)
	printf("%#12.5g", tnow);
    printf(" %11d %#11.5g %#11.5g %#11.5g %#11.5g\n",
	   nbody, avg1, rsqrt(avg2), rcbrt(avg3), rsqrt(rsqrt(avg4)));
}

void snapsum(bodyptr btab, int nbody, real tnow, bool showtime, bool showhead)
{
    double sum1, sum2, sum3, sum4;
    bodyptr bp;

    sum1 = sum2 = sum3 = sum4 = 0.0;
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	sum1 += Value(bp);
	sum2 += rsqr(Value(bp));
	sum3 += rqbe(Value(bp));
	sum4 += rsqr(rsqr(Value(bp)));
    }
    if (showhead) {
	printf(showtime ? "#%11s " : "#", "time");
	printf("%11s %11s %11s %11s %11s\n",
	       "values", "sum", "sum^2", "sum^3", "sum^4");
    }
    if (showtime)
	printf("%#12.5g", tnow);
    printf(" %11d %#11.5g %#11.5g %#11.5g %#11.5g\n",
	   nbody, sum1, sum2, sum3, sum4);
}

void snapmed(bodyptr btab, int nbody, real tnow, bool showtime, bool showhead)
{
    if (showhead) {
	printf(showtime ? "#%11s " : "#", "time");
	printf("%11s %11s %11s %11s %11s\n",
	       "minimum", "1st quart", "median", "3rd quart", "maximum");
    }
    qsort(btab, nbody, SizeofBody, cmpvalue);
    if (showtime)
	printf("%#12.5g", tnow);
    printf(" %#11.5g %#11.5g %#11.5g %#11.5g %#11.5g\n",
	   Value(NthBody(btab, 0)), 
	   Value(NthBody(btab, nbody/4)),
	   Value(NthBody(btab, nbody/2)),
	   Value(NthBody(btab, 3*nbody/4)),
	   Value(NthBody(btab, nbody - 1)));
}

void snapoct(bodyptr btab, int nbody, real tnow, bool showtime, bool showhead)
{
    if (showhead) {
	printf(showtime ? "#%7s " : "#", "time");
	printf("%9s %9s %9s %9s %9s %9s %9s\n",
	       "oct1", "oct2", "oct3", "oct4", "oct5", "oct6", "oct7");
    }
    qsort(btab, nbody, SizeofBody, cmpvalue);
    if (showtime)
	printf("%#8.5g", tnow);
    printf(" %#9.3g %#9.3g %#9.3g %#9.3g %#9.3g %#9.3g %#9.3g\n",
	   Value(NthBody(btab, (1*nbody)/8)),
	   Value(NthBody(btab, (2*nbody)/8)),
	   Value(NthBody(btab, (3*nbody)/8)),
	   Value(NthBody(btab, (4*nbody)/8)),
	   Value(NthBody(btab, (5*nbody)/8)),
	   Value(NthBody(btab, (6*nbody)/8)),
	   Value(NthBody(btab, (7*nbody)/8)));
}

void snapOCT(bodyptr btab, int nbody, real tnow, bool showtime,	bool showhead)
{
    if (showhead) {
	printf(showtime ? "#%11s " : "#", "time");
	printf("%11s %11s %11s %11s %11s %11s %11s %11s %11s\n", 
	       "min", "oct1", "oct2", "oct3", "oct4",
	       "oct5", "oct6", "oct7", "max");
    }
    qsort(btab, nbody, SizeofBody, cmpvalue);
    if (showtime)
	printf("%#12.5g", tnow);
    printf(" %#11.5g %#11.5g %#11.5g %#11.5g %#11.5g"
	   " %#11.5g %#11.5g %#11.5g %#11.5g\n",
	   Value(NthBody(btab, 0)),
	   Value(NthBody(btab, (1*nbody)/8)),
	   Value(NthBody(btab, (2*nbody)/8)),
	   Value(NthBody(btab, (3*nbody)/8)),
	   Value(NthBody(btab, (4*nbody)/8)),
	   Value(NthBody(btab, (5*nbody)/8)),
	   Value(NthBody(btab, (6*nbody)/8)),
	   Value(NthBody(btab, (7*nbody)/8)),
	   Value(NthBody(btab, nbody - 1)));
}

int cmpvalue(const void *a, const void *b)
{
    return (Value((bodyptr) a) < Value((bodyptr) b) ? -1 :
	      Value((bodyptr) a) > Value((bodyptr) b) ? 1 : 0);
}

stream execmap(string prog)
{
    int handle[2];
    char handbuf[32];

    pipe(handle);
    if (fork() == 0) {                           /* if this is child process */
        close(handle[0]);
        sprintf(handbuf, "-%d", handle[1]);
        execl(prog, getargv0(), getparam("in"), handbuf, getparam("times"),
              getparam("require"), "Value",
	      strnull(getparam("require")) ? "true" : "false",
	      getparam("seed"), NULL);
        error("%s: execl %s failed\n", getargv0(), prog);
    }
    close(handle[1]);
    sprintf(handbuf, "-%d", handle[0]);
    return (stropen(handbuf, "r"));
}
