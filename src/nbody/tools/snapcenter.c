/*
 * SNAPCENTER: transform snapshot file to coordinates centered
 * on some weighted subset of particles.
 */

#include "stdinc.h"
#include "getparam.h"
#include "strset.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "snapcenter.h"

string defv[] = {		";Center snapshot on weighted bodies",
    "in=???",                   ";Input snapshot file name",
    "out=???",                  ";Output snapshot file name",
    "times=all",                ";Range of times to process",
    "weight=1.0",		";Expression used to weight bodies",
    "coords=" PosTag "," VelTag,";Coordinate data to center",
    "require=",			";Input items required",
    "produce=",			";Output items produced",
    "passall=true",		";If true, pass on input data",
    "seed=",			";Seed for random number generator",
    "VERSION=2.3",		";Josh Barnes  26 October 2011",
    NULL,
};

void buildmap(string, string *, string *, string *, string, int);

string names[2] = { "Weight",  NULL };
string exprs[2] = { NULL,      NULL };
string types[2] = { RealType,  NULL };

stream execmap(string);				/* start snapmap process    */
void del_tag(string *, string *, string);	/* remove tag from list     */

#define WeightField  phatbody[NewBodyFields+0]
#define Weight(b)  SelectReal(b, WeightField.offset)

int main(int argc, string argv[])
{
    string prog, coords, itags[MaxBodyFields], otags[MaxBodyFields];
    stream xstr, ostr;
    bodyptr btab = NULL, bp;
    int nbody;
    real tnow;
    vector cmpos, cmvel, cmacc;

    initparam(argv, defv);
    prog = tempnam("/tmp", "sm");
    exprs[0] = getparam("weight");
    buildmap(prog, names, exprs, types, Precision, NDIM);
    xstr = execmap(prog);
    if (get_tag_ok(xstr, "History"))
	skip_item(xstr);
    get_history(xstr);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    new_field(&WeightField, RealType, "Weight");
    new_field(&WeightField + 1, NULL, NULL);
    coords = getparam("coords");
    while (get_snap(xstr, &btab, &nbody, &tnow, itags, TRUE)) {
	if (scanopt(coords, PosTag) && set_member(itags, PosTag)) {
	    snapcmpos(cmpos, btab, nbody, WeightField.offset);
	    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
		SUBV(Pos(bp), Pos(bp), cmpos);
	    }
	}
	if (scanopt(coords, VelTag) && set_member(itags, VelTag)) {
	    snapcmvel(cmvel, btab, nbody, WeightField.offset);
	    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
		SUBV(Vel(bp), Vel(bp), cmvel);
	    }
	}
	if (scanopt(coords, AccTag) && set_member(itags, AccTag)) {
	    snapcmacc(cmacc, btab, nbody, WeightField.offset);
	    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
		SUBV(Acc(bp), Acc(bp), cmacc);
	    }
	}
	del_tag(otags, itags, "Weight");
	put_snap(ostr, &btab, &nbody, &tnow, otags);
    }
    strclose(ostr);
    if (unlink(prog) != 0)
        error("%s: can't unlink %s\n", getargv0(), prog);
    return (0);
}

#include <unistd.h>

stream execmap(string prog)
{
    int handle[2];
    char handbuf[32], produce[512];

    pipe(handle);
    if (fork() == 0) {                           /* if this is child process */
        close(handle[0]);
        sprintf(handbuf, "-%d", handle[1]);
	sprintf(produce, "%s,Weight", getparam("produce"));
        execl(prog, getargv0(), getparam("in"), handbuf, getparam("times"),
              getparam("require"), produce, getparam("passall"),
	      getparam("seed"), NULL);
        error("%s: execl %s failed\n", getargv0(), prog);
    }
    close(handle[1]);
    sprintf(handbuf, "-%d", handle[0]);
    return (stropen(handbuf, "r"));
}

void del_tag(string *olist, string *ilist, string tag)
{
    string *op, *ip;

    for (op = olist, ip = ilist; *ip != NULL; ip++)
	if (! streq(*ip, tag))
	    *op++ = *ip;
    *op = NULL;
}
