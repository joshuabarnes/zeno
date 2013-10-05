/*
 * SNAPKEY.C: set key variable to given expression.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "bodytags.h"

string defv[] = {               ";Set key variable in snapshot",
    "in=???",                   ";Input snapshot file name",
    "out=???",                  ";Output snapshot file name",
    "times=all",                ";Range of times to process",
    "key=???",                  ";Expression for key value",
    "require=",			";List of input items necessary",
    "produce=" KeyTag,          ";List of output items produced",
    "passall=true",		";If true, pass on input data",
    "seed=",			";Seed for random number generator",
    "VERSION=2.1",              ";Josh Barnes  14 May 2008",
    NULL,
};

void buildmap(string, string *, string *, string *, string, int);
void execmap(string);

int main(int argc, string argv[])
{
    string prog;
    string names[2] = { KeyTag, NULL };
    string exprs[2] = { NULL,   NULL };

    initparam(argv, defv);
    prog = tempnam("/tmp", "sm");
    exprs[0] = getparam("key");
    buildmap(prog, names, exprs, NULL, Precision, NDIM);
    execmap(prog);
    if (unlink(prog) != 0)
        error("%s: can't unlink %s\n", getargv0(), prog);
    return (0);
}

#include <unistd.h>
#include <sys/wait.h>

void execmap(string prog)
{
    int mappid, mapstat;
    char histbuf[512];

    mappid = fork();
    if (mappid == 0) {                          /* if this is child process */
	sprintf(histbuf, "HISTORY=%s", *ask_history());
        execl(prog, getargv0(), getparam("in"), getparam("out"),
	      getparam("times"), getparam("require"), getparam("produce"),
	      getparam("passall"), getparam("seed"), histbuf, NULL);
        error("%s: execl %s failed\n", getargv0(), prog);
    }
    while (wait(&mapstat) != mappid)
        eprintf("[%s: waiting on subprocess %d]\n", getargv0(), mappid);
}
