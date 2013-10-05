/*
 * SNAPSORT.C: sort bodies by value.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {               ";Sort bodies by value",
    "in=???",                   ";Input snapshot file name",
    "out=???",                  ";Output snapshot file name",
    "times=all",                ";Range of times to process",
    "value=???",		";Computed value to sort by",
    "require=",			";Input items required",
    "produce=",			";Output items produced",
    "passall=true",		";If true, pass on input data",
    "seed=",			";Seed for random number generator",
    "VERSION=2.1",              ";Josh Barnes  14 May 2008",
    NULL,
};

void buildmap(string, string *, string *, string *, string, int);

string names[2] = { "Value",  NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { RealType, NULL };

stream execmap(string);				/* start snapmap process    */
int cmpvalue(const void *, const void *);	/* compare values of bodies */
void del_tag(string *, string *, string);	/* remove tag from list     */

#define ValueField  phatbody[NewBodyFields+0]
#define Value(b)  SelectReal(b, ValueField.offset)

int main(int argc, string argv[])
{
    string prog, itags[MaxBodyFields], otags[MaxBodyFields];
    stream xstr, ostr;
    bodyptr btab = NULL;
    int nbody;
    real tnow;

    initparam(argv, defv);
    prog = tempnam("/tmp", "sm");
    exprs[0] = getparam("value");
    buildmap(prog, names, exprs, types, Precision, NDIM);
    xstr = execmap(prog);
    if (get_tag_ok(xstr, "History"))
	skip_item(xstr);
    get_history(xstr);
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    new_field(&ValueField, RealType, "Value");
    new_field(&ValueField + 1, NULL, NULL);
    while (get_snap(xstr, &btab, &nbody, &tnow, itags, TRUE)) {
	qsort(btab, nbody, SizeofBody, cmpvalue);
	del_tag(otags, itags, "Value");
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
	sprintf(produce, "%s,Value", getparam("produce"));
        execl(prog, getargv0(), getparam("in"), handbuf, getparam("times"),
              getparam("require"), produce, getparam("passall"),
	      getparam("seed"), NULL);
        error("%s: execl %s failed\n", getargv0(), prog);
    }
    close(handle[1]);
    sprintf(handbuf, "-%d", handle[0]);
    return (stropen(handbuf, "r"));
}

int cmpvalue(const void *a, const void *b)
{
    return (Value((bodyptr) a) < Value((bodyptr) b) ? -1 :
	      Value((bodyptr) a) > Value((bodyptr) b) ? 1 : 0);
}

void del_tag(string *olist, string *ilist, string tag)
{
    string *op, *ip;

    for (op = olist, ip = ilist; *ip != NULL; ip++)
	if (! streq(*ip, tag))
	    *op++ = *ip;
    *op = NULL;
}
