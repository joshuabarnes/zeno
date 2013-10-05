/*
 * SNAPTRAK.C: track centroids of user-specified particle groups.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Track centroids of specified groups",
    "in=???",			";Input snapshot file name",
    "out=???",			";Output snapshot file name",
    "group=???",		";Expression for group membership",
    "times=all",		";Range of times to process",
    "seed=",			";Seed for random number generator",
    "VERSION=2.1",              ";Josh Barnes  12 April 2011",
    NULL,
};

string btags[] = { MassTag, PosTag, VelTag, KeyTag, "Group", NULL};
string otags[] = { MassTag, PosTag, VelTag, KeyTag, NULL};

void buildmap(string, string *, string *, string *, string, int);

string names[2] = { "Group",  NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { IntType,  NULL };

local stream execmap(string);			/* start snapmap process    */
local void snaptrak(void);			/* track groups of bodies   */

#define GroupField  phatbody[NewBodyFields+0]
#define Group(b)  SelectInt(b, GroupField.offset)

bodyptr bodytab = NULL, traktab = NULL;
int nbody, ntrak;
real tbody;

int main(int argc, string argv[])
{
    string prog, itags[MaxBodyFields];
    stream xstr, ostr;
    int nold = -1;

    initparam(argv, defv);
    prog = tempnam("/tmp", "sm");
    exprs[0] = getparam("group");
    buildmap(prog, names, exprs, types, Precision, NDIM);
    xstr = execmap(prog);
    if (get_tag_ok(xstr, "History"))
	skip_item(xstr);
    get_history(xstr);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    new_field(&GroupField, IntType, "Group");
    new_field(&GroupField + 1, NULL, NULL);
    layout_body(btags, Precision, NDIM);
    while (get_snap(xstr, &bodytab, &nbody, &tbody, itags, FALSE)) {
	snaptrak();
	put_snap(ostr, &traktab, &ntrak, &tbody, otags);
	if (ntrak != nold)
	  eprintf("[%s: wrote %d groups at t = %f]\n",
		  getargv0(), ntrak, tbody);
	nold = ntrak;
    }
    strclose(ostr);
    if (unlink(prog) != 0)
        error("%s: can't unlink %s\n", getargv0(), prog);
    return (0);
}

#include <unistd.h>

#define StdTags  MassTag "," PosTag "," VelTag

stream execmap(string prog)
{
    int handle[2];
    char handbuf[32];

    pipe(handle);
    if (fork() == 0) {                           /* if this is child process */
        close(handle[0]);
        sprintf(handbuf, "-%d", handle[1]);
        execl(prog, getargv0(), getparam("in"), handbuf, getparam("times"),
	      StdTags, StdTags ",Group", "true", getparam("seed"), NULL);
        error("%s: execl %s failed\n", getargv0(), prog);
    }
    close(handle[1]);
    sprintf(handbuf, "-%d", handle[0]);
    return (stropen(handbuf, "r"));
}

void snaptrak(void)
{
    bodyptr bp, gp;
    int nzero;

    if (traktab == NULL) {
	ntrak = 0;
	for (bp = bodytab; bp < NthBody(bodytab, nbody); bp = NextBody(bp))
	    ntrak = MAX(ntrak, Group(bp));
	eprintf("[%s: allocating %d groups]\n", getargv0(), ntrak);
	traktab = (bodyptr) allocate(ntrak * SizeofBody);
    }
    for (gp = traktab; gp < NthBody(traktab, ntrak); gp = NextBody(gp)) {
	Mass(gp) = 0.0;
	CLRV(Pos(gp));
	CLRV(Vel(gp));
	Key(gp) = 0;
    }
    for (bp = bodytab; bp < NthBody(bodytab, nbody); bp = NextBody(bp)) {
	if (Group(bp) > ntrak)
	    error("snaptrak: cant expand group array\n");
	if (Group(bp) > 0) {
	    gp = NthBody(traktab, Group(bp) - 1);
	    Mass(gp) += Mass(bp);
	    ADDMULVS(Pos(gp), Pos(bp), Mass(bp));
	    ADDMULVS(Vel(gp), Vel(bp), Mass(bp));
	    Key(gp)++;
	}
    }
    nzero = 0;
    for (gp = traktab; gp < NthBody(traktab, ntrak); gp = NextBody(gp))
	if (Mass(gp) != 0.0) {
	    DIVVS(Pos(gp), Pos(gp), Mass(gp));
	    DIVVS(Vel(gp), Vel(gp), Mass(gp));
	} else
	    nzero++;
    if (nzero > 0)
	eprintf("[%s: %d groups have zero mass]\n", getargv0(), nzero);
}
