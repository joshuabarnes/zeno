/*
 * setrange.c: parse range a expression and assign values.  The range
 * expression has the form
 *
 *	range = <value>:<value> | <value> | <value>:
 *	value = <number> | <number>/<number>
 *
 * In a range expression with only one <value>, the other is taken
 * to be zero; <value> is interpreted as 0:<value>, while <value>:
 * is interpreted as <value>:0.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include <stdlib.h>

void setrange(real *rval, string rexp)
{
  string rptr;
    
  rval[0] = strtod(rexp, &rptr);
  if (*rptr == '/')
    rval[0] = rval[0] / strtod(rptr + 1, &rptr);
  if (*rptr == ':') {
    rval[1] = strtod(rptr + 1, &rptr);
    if (*rptr == '/')
      rval[1] = rval[1] / strtod(rptr + 1, &rptr);
  } else {
    rval[1] = rval[0];
    rval[0] = 0.0;
  }
  if (*rptr != (char) NULL)
    eprintf("[%s.setrange: WARNING: ignoring trailing \"%s\"]\n",
	    getprog(), rptr);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
    "range=1/64:3.2e1",
    NULL,
};

main(int argc, string *argv)
{
  real range[2];

  initparam(argv, defv);
  setrange(range, getparam("range"));
  printf("%s: %s -> %f,%f\n", getprog(), getparam("range"),
	 range[0], range[1]);
}

#endif
