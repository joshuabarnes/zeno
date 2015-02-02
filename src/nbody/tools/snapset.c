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
				";Expressions (C code) may use these values:",
				  SNAPMAP_BODY_VARS ".",
  "in=???",                     ";Input snapshot file name",
  "out=???",                    ";Output snapshot file name",
  "times=all",                  ";Range of times to process",
  "x=",				";Expression for x position",
  "y=",				";Expression for y position",
  "z=",				";Expression for z position",
  "vx=",			";Expression for x velocity",
  "vy=",			";Expression for y velocity",
  "vz=",			";Expression for z velocity",
  "ax=",			";Expression for x acceleration",
  "ay=",			";Expression for y acceleration",
  "az=",			";Expression for z acceleration",
  "m=",				";Expression for mass",
  "phi=",                       ";Expression for potential",
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
  "keyarr=",			";Expressions for key array",
  "aux=",                       ";Expression for aux value",
  "auxvx=",                     ";Expression for X aux vector",
  "auxvy=",                     ";Expression for Y aux vector",
  "auxvz=",                     ";Expression for Z aux vector",
  "auxarr=",			";Expressions for aux array",
  "t=",				";Expression used to map time.",
				";May use these values: "
				  SNAPMAP_TIME_VARS ".",
  "require=",                   ";List of input items necessary",
  "produce=",                   ";List of output items produced",
  "passall=true",		";If true, pass on all input data",
  "seed=",			";Generator seed for random values",
  "VERSION=2.9",                ";Josh Barnes  2 February 2015",
  NULL,
};

void execmap(string prog);

int main(int argc, string argv[])
{
  string *mdtab, *names, *exprs, prog;
  int nexp = 0, i, j;

  initparam(argv, defv);
  mdtab = getmapdefs();				// get list of mapping vars
  for (i = 0; mdtab[i] != NULL; i += 2)
    if (getparamstat(mdtab[i]) & ARGPARAM)	// if var has assigned value
      nexp++;
  eprintf("[%s: %scounted %d variable assignments]\n", getprog(),
	  nexp > 0 ? "" : "warning: ", nexp);
  names = (string *) allocate(sizeof(string *) * (nexp + 1));
  exprs = (string *) allocate(sizeof(string *) * (nexp + 1));
  for (i = j = 0; mdtab[i] != NULL; i += 2)
    if (getparamstat(mdtab[i]) & ARGPARAM) {	// if var has assigned value
      exprs[j] = getparam(mdtab[i]);		// list value given as expr
      names[j] = mdtab[i+1];			// and name of access macro
      j++;
    }
  exprs[j] = names[j] = NULL;
  prog = mktemp((string) copxstr("/tmp/sm_XXXXXX", sizeof(char)));
  buildmap(prog, names, exprs, NULL,
	   strnull(getparam("t")) ? NULL : getparam("t"),
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
