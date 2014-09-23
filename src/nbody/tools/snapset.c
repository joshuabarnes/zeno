/*
 * snapset.c: set body variables to given expressions.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "bodytags.h"
#include "buildmap.h"
#include <unistd.h>
#include <sys/wait.h>
#include <assert.h>

string defv[] = {               ";Set values for each body in snapshot.",
				";Expressions are written in C language.",
				";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "in=???",                     ";Input snapshot file name",
  "out=???",                    ";Output snapshot file name",
  "times=all",                  ";Range of times to process",
  "x=",				";Expression for X position",
  "y=",				";Expression for Y position",
  "z=",				";Expression for Z position",
  "vx=",			";Expression for X velocity",
  "vy=",			";Expression for Y velocity",
  "vz=",			";Expression for Z velocity",
  "m=",				";Expression for mass",
  "phi=",                       ";Expression for potential",
  "ax=",			";Expression for X acceleration",
  "ay=",			";Expression for Y acceleration",
  "az=",			";Expression for Z acceleration",
  "smooth=",			";Expression for smoothing length",
  "rho=",			";Expression for density",
  "entf=",			";Expression for entropy function",
  "uint=",			";Expression for internal energy",
  "udot=",			";Expression for derivative of uint",
  "udotrad=",			";Expression for radiative losses",
  "udotvis=",			";Expression for viscous heating",
  "tau=",			";Expression for optical depth",
  "type=",			";Expression for body type",
  "birth=",			";Expression for birth date",
  "death=",			";Expression for death date",
  "key=",			";Expression for key value",
  "aux=",                       ";Expression for aux value",
  "auxvx=",                     ";Expression for X aux vector",
  "auxvy=",                     ";Expression for Y aux vector",
  "auxvz=",                     ";Expression for Z aux vector",
  "tmap=",			";Expression used to map time.",
				";Bound variables are: "
				  SNAPMAP_TIME_VARS ".",
  "require=",                   ";List of input items necessary",
  "produce=",                   ";List of output items produced",
  "passall=true",		";If true, pass on all input data",
  "seed=",			";Seed for random number generator",
  "VERSION=2.8",                ";Josh Barnes  9 Sep 2014",
  NULL,
};

void execmap(string prog);

int main(int argc, string argv[])
{
  string *mdtab, *names, *exprs, prog;
  int i, j;

  initparam(argv, defv);
  mdtab = getmapdefs();
  for (i = 0; mdtab[i] != NULL; i++)
    continue;
  assert((i & 1) == 0);
  names = (string *) allocate(sizeof(string *) * i / 2);
  exprs = (string *) allocate(sizeof(string *) * i / 2);
  for (i = j = 0; mdtab[i] != NULL; i += 2)
    if (getparamstat(mdtab[i]) & ARGPARAM) {	// if ident has been assigned
      exprs[j] = getparam(mdtab[i]);		// list given value as expr
      names[j] = mdtab[i+1];			// and access macro as name
      j++;
    }
  exprs[j] = names[j] = NULL;
  eprintf("[%s: defined %d variable%s]\n", getargv0(), j, j!=1 ? "s" : "");
  prog = tempnam("/tmp", "sm");
  buildmap(prog, names, exprs, NULL,
	   strnull(getparam("tmap")) ? NULL : getparam("tmap"),
	   Precision, NDIM, TRUE);
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
  if (mappid == 0) {				// if this is child process
    sprintf(histbuf, "HISTORY=%s", *ask_history());
    execl(prog, getargv0(), getparam("in"), getparam("out"),
	  getparam("times"), getparam("require"), getparam("produce"),
	  getparam("passall"), getparam("seed"), histbuf, NULL);
    error("%s: execl %s failed\n", getargv0(), prog);
  }
  while (wait(&mapstat) != mappid)
    eprintf("[%s: waiting on subprocess %d]\n", getargv0(), mappid);
}
