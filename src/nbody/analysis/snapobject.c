/*
 * SNAPOBJECT.C: identify objects using percolation algorithm.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "strset.h"
#include "phatbody.h"

string defv[] = {		";Identify objects using percolation",
    "in=???",			";Input file with N-body system",
    "out=???",			";Output file with identifications",
    "bcrit=0.1",		";Percolation length parameter",
    "nmin=12",			";Minimum members per listed group",
    "subkey=false",		";If TRUE, link bodies with key != 0",
    "VERSION=1.1",		";Josh Barnes  17 June 1999",
    NULL,
};

#define LinkTag    "Link"
#define LinkField  phatbody[NewBodyFields+0]
#define Link(b)    (*((void **)((byte *)(b) + (LinkField.offset))))

string bodytags[] = {  PosTag, KeyTag, LinkTag, NULL };

bodyptr btab = NULL;		/* array of bodies to process */
int nbody;			/* number of bodies in btab */
real tnow = 0.0;		/* time body data is defined */

void findobj(real, int, bool);
void linkpart(bodyptr, bodyptr);
bodyptr toplink(bodyptr);

int main(int argc, string argv[])
{
    stream istr, ostr;
    string itags[MaxBodyFields];

    initparam(argv, defv);
    new_field(&LinkField, IntType, LinkTag);	/* use int's worth of space */
    new_field(&LinkField + 1, NULL, NULL);
    layout_body(bodytags, Precision, NDIM);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    if (! get_snap(istr, &btab, &nbody, &tnow, itags, TRUE))
	error("%s: snapshot input failed\n", getargv0());
    if (! set_member(itags, PosTag))
	error("%s: %s data missing\n", getargv0(), PosTag);
    if (getbparam("subkey") && ! set_member(itags, KeyTag))
	error("%s: %s data missing\n", getargv0(), KeyTag);
    findobj(getdparam("bcrit"),	getiparam("nmin"), getbparam("subkey"));
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_snap(ostr, &btab, &nbody, &tnow,
	     set_union(itags, set_cons(KeyTag, NULL)));
    strclose(ostr);
    return (0);
}

/*
 * FINDOBJ: identify objects by percolation.
 */

void findobj(real bcrit, int nmin, bool subkey)
{
    bodyptr p, q;
    int nobj;

    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p)) {
	Link(p) = p;
	Key(p) = (subkey && Key(p) == 0 ? -1 : 0);
    }
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
	if (Key(p) > -1)
	    for (q = btab; q < p; q = NextBody(q))
		if (Key(q) > -1 && distv(Pos(p), Pos(q)) < bcrit)
		    linkpart(p, q);
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p)) {
	Link(p) = toplink(p);
	Key(Link(p))++;
    }
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
	if (Key(Link(p)) < nmin)
	    Link(p) = NULL;
    nobj = 0;
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
	if (Link(p) == p)
	    Key(p) = ++nobj;
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
	Key(p) = (Link(p) != NULL ? Key(Link(p)) : 0);
    eprintf("[%s: %d objects found at t = %.3f]\n", getargv0(), nobj, tnow);
}

/*
 * LINKPART: put particles in same equiv. class.
 */

void linkpart(bodyptr p, bodyptr q)
{
    p = toplink(p);
    q = toplink(q);
    if (p < q)
	Link(q) = p;
    else
	Link(p) = q;
}

/*
 * TOPLINK: follow link chain to self-referent particle.
 */

bodyptr toplink(bodyptr p)
{
    while (Link(p) != p)
	p = Link(p);
    return (p);
}
