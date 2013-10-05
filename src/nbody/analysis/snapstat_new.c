/*
 * snapstat.c: calculate statistics for input snapshots.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "getparam.h"
#include "filestruct.h"
#define SafeSelect
#include "phatbody.h"
#include <string.h>
#include <dlfcn.h>

string defv[] = {		";Analyze statistics of SnapShot data.",
				";Use doubles for sums and averages.",
    "in=???",			";Input file name",
    "times=all",		";Times to analyze",
    "value=x",			";Computed value to analyze",
    "require=*",		";Input items required",
    "options=avg,med",		";Statistical measures to tabulate.",
				";Possible choices are one or more of:",
				";avg: list average, rms, etc.",
				";sum: list sum, sum of squares, etc.",
				";med: list median, quartiles, & limits.",
				";oct: list 1st through 7th octiles.",
				";OCT: list octiles & limits (wide!).",
				";time: prepend current time to line.",
    "seed=",			";Seed for random number generator",
    "VERSION=3.0",		";Josh Barnes  31 July 2012",
    NULL,
};

#define ValueField  phatbody[NewBodyFields+0]
#define Value(b)  SelectReal(b, ValueField.offset)

string names[2] = { "Value",  NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { RealType, NULL };

void snapstat(bodyptr, int, real, string);
void snapavg(bodyptr, int, real, bool, bool);
void snapsum(bodyptr, int, real, bool, bool);
void snapmed(bodyptr, int, real, bool, bool);
void snapoct(bodyptr, int, real, bool, bool);
void snapOCT(bodyptr, int, real, bool, bool);
int cmpvalue(const void *, const void *);

local void checktags(string *tags, string *req)
{
  string *rp;

  for (rp = req; *rp != NULL; rp++)
    if (! set_member(tags, *rp))
      error("%s: input data %s missing\n", getargv0(), *rp);
}

local void snapmap(bodyptr btab, int nbody, real tnow,
		   void (*mapfunc)(bodyptr, bodyptr, real, int, int))
{
  bodyptr cp, bp;
  int i;

  cp = (bodyptr) allocate(SizeofBody);		// alloc space for copy
  for (i = 0; i < nbody; i++) {			// loop over all bodies
    bp = NthBody(btab, i);			// get address of each
    memcpy(cp, bp, SizeofBody);			// make copy of contents
    (*mapfunc)(bp, cp, tnow, i, nbody);		// perform transformation
  }
  free(cp);
}

int main(int argc, string argv[])
{
  stream istr;
  string times, *require, prog, err, itags[MaxBodyFields];
  bool expand, firstloop = TRUE;
  void *handle;
  char buf[256];
  void (*mapfunc)(bodyptr, bodyptr, real, int, int);
  bodyptr btab = NULL;
  int nbody;
  real tnow;

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  times = getparam("times");
  exprs[0] = getparam("value");

  new_field(&ValueField, RealType, "Value");
  new_field(&ValueField + 1, NULL, NULL);
  layout_body(names, Precision, NDIM);
  expand = streq(getparam("require"), "*");
  if (! expand) {
    require = burststring(getparam("require"), ", ");
    layout_body(require, Precision, NDIM);
  }

  sprintf(buf, "sm_%s.so", getparam("value"));
  handle = dlopen(buf, RTLD_LAZY);
  if (! handle)
    error("%s: error opening %s\n%s\n", getprog(), buf, dlerror());
  mapfunc = dlsym(handle, "computemap");
  if ((err = dlerror()) != NULL)
    error("%s: error finding %s\n%s\n", getprog(), "computemap", err);

  // prog = tempnam("/tmp", "sm");

  while (get_snap_t(istr, &btab, &nbody, &tnow, itags, expand, times)) {
    if (firstloop && ! expand)
      checktags(itags, require);
    snapmap(btab, nbody, tnow, mapfunc);
    snapstat(btab, nbody, tnow, getparam("options"));
    firstloop = FALSE;
  }

  // if (unlink(prog) != 0)
  //   error("%s: can't unlink %s\n", getargv0(), prog);

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
