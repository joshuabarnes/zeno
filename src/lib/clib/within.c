/*
 * within.c: determine if a floating point number is within specified range.
 */

#include "stdinc.h"
#include "mathfns.h"
#include <string.h>

#define COMMA  ','
#define COLON  ':'

//  within: determine if a floating point number is within specified
//  range list, passed in the form "<range>_1,<range>_2,..." where
//  each <range> is either one floating point number, two numbers
//  "<low>:<high>", or three numbers "<low>:<high>:<step>".  To allow
//  for small uncertainties in the values tested, floating-point
//  comparison is done with a specified fuzzyness parameter.  Now
//  accepts <range> of form
//  _________________________________________________________________


bool within(double val, string ranges, double fuzz)
{
  char *endptr, *rngptr, *comptr, *colptr;
  double rngval, rngmax, rngstep, x;

  endptr = ranges + strlen(ranges);		// point to terminating NULL
  for (rngptr = ranges; rngptr != endptr; ) {	// loop over each range
    comptr = strchr(rngptr, COMMA);		// point to end of range
    comptr = (comptr==NULL ? endptr : comptr);	// use final NULL for last
    colptr = strchr(rngptr, COLON);		// scan ahead for colon
    colptr = (colptr > comptr ? NULL : colptr);	// keep if in this range
    if (colptr == NULL) {			// single value specified?
      rngval = atof(rngptr);			// get value as given
      if (ABS(val - rngval) < fuzz)		// make fuzzy comparison
	return (TRUE);				// done if value matches
    } else {					// handle true range
      rngval = atof(rngptr);			// set low end of range
      rngmax = atof(colptr+1);			// set high end of range
      if (rngval - fuzz <= val &&		// within fuzzy range?
	  val <= rngmax + fuzz) {
	colptr = strchr(colptr+1, COLON);	// look for another colon
	if (colptr == NULL || colptr > comptr)	// if none in this range
	  return (TRUE);			// know value is in range
	else {					// handle stepped range
	  rngstep = atof(colptr+1);		// get step value
	  x = (val - rngval) / rngstep;		// compute step number
	  if (ABS(x-round(x)) <= fuzz/rngstep)	// close enough to integer?
	    return (TRUE);			// count it as a match
	}
      }
    }
    rngptr = (*comptr == COMMA ? comptr+1 : comptr);
						// advance to next range
  }
  return (FALSE);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
    "val=1.0",
    "range=0.5:0.7,1.0,1.2:1.3",
    "fuzz=0.001",
    NULL,
};

int main(int argc, string argv[])
{
  initparam(argv, defv);
  if (within(getdparam("val"), getparam("range"), getdparam("fuzz")))
    printf("within returns TRUE\n");
  else
    printf("within returns FALSE\n");
  return 0;
}

#endif

#if 0

bool within(double val, string range, double fuzz)
{
  char *endptr, *subptr, *sepptr, *colptr;
  double sublow, subhi;

  endptr = range + strlen(range);		// point to term. NULL
  for (subptr = range; subptr != endptr; ) {	// for each subrange
    sepptr = strchr(subptr, COMMA);		// pnt to subrange end
    if (sepptr == NULL)				// last subrange listed?
      sepptr = endptr;				// fix up subend ptr
    colptr = strchr(subptr, COLON);		// scan subrange for
    if (colptr > sepptr)			// in another subrange?
      colptr = NULL;				// then dont use it
    sublow = atof(subptr) - fuzz/2.0;		// set low end of range
    if (colptr != NULL)				// high end specified?
      subhi = atof(colptr+1) + fuzz/2.0;	// set high end
    else
      subhi = sublow + fuzz;			// just use low end
    if (sublow <= val && val <= subhi)		// within subrange?
      return (TRUE);
    subptr = sepptr;				// advance subrange ptr
    if (*subptr == COMMA)			// more ranges to do?
      subptr++;					// move on to next
  }
  return (FALSE);
}

#endif

