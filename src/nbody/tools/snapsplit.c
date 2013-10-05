/*
 * SNAPSPLIT.C: Split snapshot into samples.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Split snapshot into samples",
    "in=???",			";Snapshot input file",
    "out=???",			";Snapshot output stream",
    "nbody=1024",		";Nubmer per output sample",
    "produce=*",		";List of output items",
    "VERSION=1.0",		";Josh Barnes  26 October 2004",
    NULL,
};

int main(int argc, string argv[])
{
  stream istr, ostr;
  string *produce, iotags[MaxBodyFields];
  bodyptr btab = NULL, bptr;
  int nbody, ninp, nout, nleft;
  real tsnap;
  bool expand;

  initparam(argv, defv);
  expand = streq(getparam("produce"), "*");
  if (! expand) {
    produce = burststring(getparam("produce"), ", ");
    layout_body(produce, Precision, NDIM);
  }
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab, &ninp, &tsnap, iotags, expand))
    error("%s: snapshot input failed\n", getargv0());
  if (! (expand || set_equal(iotags, produce)))
    eprintf("[%s: outputting only %d fields]\n",
	    getargv0(), set_length(iotags));
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  nbody = getiparam("nbody");
  nleft = ninp;
  bptr = btab;
  while (nleft > 0) {
    nout = MIN(nleft, nbody);
    if (nout < nbody)
      eprintf("[%s outputting only %d bodies]\n", getargv0(), nout);
    put_snap(ostr, &bptr, &nout, &tsnap, iotags);
    nleft -= nout;
    bptr = NthBody(bptr, nout);
  }
  strclose(ostr);
  return (0);
}
