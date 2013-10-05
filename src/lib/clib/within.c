/*
 * within.c: determine if a floating point number is within specified range,
 * represented as a string of the form "<subrange1>,<subrange2>,..." where
 * each <range> is either a single floating point number, or a pair of
 * numbers seperated by a ":".  To allow for small uncertainties in the
 * values tested, floating-point comparison is done with a specified
 * fuzzyness parameter.
 */

#include "stdinc.h"
#include <string.h>

bool within(double val, string range, double fuzz)
{
    char *endptr, *subptr, *sepptr, *colptr;
    double sublow, subhi;

    endptr = range + strlen(range);		// point to term. NULL
    for (subptr = range; subptr != endptr; ) {	// for each subrange
	sepptr = strchr(subptr, ',');		// pnt to subrange end
	if (sepptr == NULL)			// last subrange listed?
	    sepptr = endptr;			// fix up subend ptr
	colptr = strchr(subptr, ':');		// scan subrange for
	if (colptr > sepptr)			// in another subrange?
	    colptr = NULL;			// then dont use it
	sublow = atof(subptr) - fuzz/2.0;	// set low end of range
	if (colptr != NULL)			// high end specified?
	    subhi = atof(colptr+1) + fuzz/2.0;	// set high end
	else
	    subhi = sublow + fuzz;		// just use low end
	if (sublow <= val && val <= subhi)	// within subrange?
	    return (TRUE);
	subptr = sepptr;			// advance subrange ptr
	if (*subptr == ',')			// more ranges to do?
	    subptr++;				// move on to next
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

void main(int argc, string argv[])
{
    initparam(argv, defv);
    if (within(getdparam("val"), getparam("range"), getdparam("fuzz")))
	printf("within returns TRUE\n");
    else
	printf("within returns FALSE\n");
}

#endif
