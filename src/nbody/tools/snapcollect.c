/*
 * SNAPCOLLECT: read several frames from a snapshot input file,
 * and output single frame with the latest data items.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Collect info from several frames",
    "in=???",			";N-body snapshot input file",
    "out=???",			";N-body snapshot output file",
    "times=all",		";Range of times to collect",
    "produce=*",		";List of output items",
    "VERSION=2.0",		";Josh Barnes  4 November 1997",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr, ostr;
    string times, *produce, itags[MaxBodyFields], atags[MaxBodyFields];
    bodyptr btab = NULL;
    int nbody, ntag, i;
    real tsnap;
    bool expand;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    times = getparam("times");
    expand = streq(getparam("produce"), "*");
    if (! expand) {
	produce = burststring(getparam("produce"), ", ");
	layout_body(produce, Precision, NDIM);
    }
    atags[ntag = 0] = NULL;
    while (get_snap_t(istr, &btab, &nbody, &tsnap, itags, expand, times)) {
	for (i = 0; itags[i] != NULL; i++)
	    if (! set_member(atags, itags[i])) {
		atags[ntag++] = itags[i];
		atags[ntag] = NULL;
	    }
	get_history(istr);
    }
    if (! expand) {
	for (i = 0; produce[i] != 0; i++)
	    if (! set_member(atags, produce[i]))
		error("%s: field %s missing\n", getargv0(), produce[i]);
    } else
	produce = atags;
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_snap(ostr, &btab, &nbody, &tsnap, produce);
    strclose(ostr);
    return (0);
}
