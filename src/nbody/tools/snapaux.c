/*
 * snapaux.c: set aux variable to given expression.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "bodytags.h"
#include "buildmap.h"
#include <unistd.h>
#include <sys/wait.h>

string defv[] = {               ";Set aux value for each body in snapshot",
  "in=???",                     ";Input snapshot file name",
  "out=???",                    ";Output snapshot file name",
  "times=all",                  ";Range of times to process",
  "aux=???",                    ";C language expression for aux value.",
				";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "require=",			";List of input items necessary",
  "produce=" AuxTag,            ";List of output items produced",
  "passall=true",		";If true, pass on input data",
  "seed=",			";Seed for random number generator",
  "VERSION=2.1",                ";Josh Barnes  9 Sep 2014",
  NULL,
};

void execmap(string);

int main(int argc, string argv[])
{
  string prog, names[2] = { AuxTag, NULL }, exprs[2] = { NULL, NULL };

  initparam(argv, defv);
  prog = tempnam("/tmp", "sm");
  exprs[0] = getparam("aux");
  buildmap(prog, names, exprs, NULL, NULL, Precision, NDIM, TRUE);
  execmap(prog);
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getargv0(), prog);
  return (0);
}

//  execmap: use compiled program to map input snap to output snap.
//  _______________________________________________________________

void execmap(string prog)
{
  int mappid, mapstat;
  char histbuf[512];

  mappid = fork();
  if (mappid == 0) {                            // if this is child process
    sprintf(histbuf, "HISTORY=%s", *ask_history());
    execl(prog, getargv0(), getparam("in"), getparam("out"),
	  getparam("times"), getparam("require"), getparam("produce"),
	  getparam("passall"), getparam("seed"), histbuf, NULL);
    error("%s: execl %s failed\n", getargv0(), prog);
  }
  while (wait(&mapstat) != mappid)
    eprintf("[%s: waiting on subprocess %d]\n", getargv0(), mappid);
}
