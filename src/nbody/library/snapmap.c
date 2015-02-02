/*
 * snapmap.c: program to implement arbitrary body transformation.
 */

#include "stdinc.h"
#include "getparam.h"
#include "mathfns.h"
#include "vectdefs.h"
#include "phatbody.h"

string defv[] = {
  "in=???",			";N-body snapshot input file",
  "out=???",			";N-body snapshot output file",
  "times=all",			";Range of times to process",
  "require=",			";List of input data required",
  "produce=",			";List of output data produced ",
  "passall=false",		";If true, output all data",
  "seed=",			";Seed for random number generator",
  "HISTORY=",			";History from client program",
  "VERSION=1.3",		";Josh Barnes  8 Sep 2014",
  NULL,
};

local void definebody(string *, string *);	// layout all body fields
local void snapmap(bodyptr, int, real);		// perform client's mapping
local void checktags(string *, string *);	// catch missing input data
local void jointags(string *, string *);	// include absent tags

//  extendbody, computemap, computetime: these routines are generated
//  by the client.  They are compiled seperately and linked with this
//  code to produce an instance of a snapshot transformation program.
//  _________________________________________________________________

void extendbody(void);				// add nonstandard fields
void computemap(bodyptr, bodyptr, real, int, int);	// map one body
real computetime(real, int);			// compute new time

int main(int argc, string argv[])
{
  stream istr, ostr;
  string *require, *produce, iotags[MaxBodyFields];
  bool firstsnap = TRUE;
  bodyptr btab = NULL;
  int nbody;
  real tnow, tnew;

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  require = burststring(getparam("require"), ", ");
  produce = burststring(getparam("produce"), ", ");
  extendbody();					// add client's extensions
  definebody(require, produce);			// set up all body fields
  if (! strnull(getparam("seed")))
    init_random(getiparam("seed"));
  while (get_snap_t(istr, &btab, &nbody, &tnow, iotags,
		    firstsnap && getbparam("passall"), getparam("times"))) {
    if (firstsnap)
      checktags(iotags, require);		// check for required data
    snapmap(btab, nbody, tnow);			// map particle array
    tnew = computetime(tnow, nbody);		// and snapshot time value
    if (getbparam("passall")) {
      jointags(iotags, produce);		// form union w/ input set
      put_snap(ostr, &btab, &nbody, &tnew, iotags);
    } else
      put_snap(ostr, &btab, &nbody, &tnew, produce);
    fflush(ostr);
    firstsnap = FALSE;
    get_history(istr);				// gobble any history items
  }
  return (0);
}	

//  definebody: construct list of body fields and layout body structure.
//  ____________________________________________________________________

local void definebody(string *require, string *produce)
{
  string tags[MaxBodyFields] = { PosTag, VelTag, };

  if ((set_member(require, PosTag) || set_member(produce, PosTag)) &&
      (set_member(require, VelTag) || set_member(produce, VelTag)))
    tags[2] = NULL;				// keep Pos, Vel at front
  else
    tags[0] = NULL;				// use arbitrary tag order
  jointags(tags, require);
  jointags(tags, produce);
  layout_body(tags, Precision, NDIM);
}

//  snapmap: perform mapping on particle array, one body at a time.
//  _______________________________________________________________

local void snapmap(bodyptr btab, int nbody, real tnow)
{
  bodyptr cp, bp;
  int i;

  cp = (bodyptr) allocate(SizeofBody);		// create temporary space
  for (i = 0; i < nbody; i++) {			// loop over all bodies
    bp = NthBody(btab, i);			// get address of each
    memcpy(cp, bp, SizeofBody);			// make copy of body
    computemap(bp, cp, tnow, i, nbody);		// perform transformation
  }
  free(cp);
}

//  checktags: make sure all required input fields are present.
//  ___________________________________________________________

local void checktags(string *iotags, string *require)
{
  string *rp;

  for (rp = require; *rp != NULL; rp++)
    if (! set_member(iotags, *rp))
      error("%s: input data %s missing\n", getargv0(), *rp);
}

//  jointags: merge tags into existing list.
//  ________________________________________

local void jointags(string *list, string *more)
{
  string *lp, *mp;

  for (lp = list; *lp != NULL; lp++)		// find end of current list
    continue;
  for (mp = more; *mp != NULL; mp++)		// loop over new tag list
    if (! set_member(list, *mp)) {
      if (lp == list + MaxBodyFields - 1)
	error("%s: field list overflow\n", getargv0());
      *lp++ = *mp;
      *lp = NULL;
    }
}
