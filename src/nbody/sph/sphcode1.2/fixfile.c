#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {
    "in=???",
    "out=???",
    "VERSION=1",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr, ostr;
    bodyptr btab = NULL;
    int nbody;
    real tlast, tnow;
    string intags[MaxBodyFields];

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    tlast = -1.0;
    while (get_snap(istr, &btab, &nbody, &tnow, intags, TRUE))
	if (tnow > tlast) {
	    put_snap(ostr, &btab, &nbody, &tnow, intags);
	    tlast = tnow;
	}
    strclose(ostr);
    exit(0);
}
