/*
 * snapsift.c: select bodies obeying predicate.
 */

#include "stdinc.h"
#include "string.h"
#include "getparam.h"
#include "vectdefs.h"
#include "filestruct.h"
#include "phatbody.h"
#include "buildmap.h"
#include <unistd.h>

string defv[] = {               ";Select bodies obeying predicate",
  "in=???",                     ";Input snapshot file name",
  "out=???",                    ";Output snapshot file name",
  "times=all",                  ";Range of times to process",
  "sieve=???",                  ";C language predicate to select bodies.",
				";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "require=",			";Input items required",
  "produce=",			";Output items produced",
  "passall=true",		";If true, pass on input data",
  "seed=",			";Seed for random number generator",
  "VERSION=2.1",                ";Josh Barnes  9 Sep 2014",
  NULL,
};

void snapsift(bodyptr, int *, bodyptr, int);	// sift array of bodies
stream execmap(string);				// start snapmap process
void del_tag(string *, string *, string);	// remove tag from list

string names[2] = { "Sieve",  NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { BoolType, NULL };

#define SieveField  phatbody[NewBodyFields+0]
#define Sieve(b)  SelectBool(b, SieveField.offset)

int main(int argc, string argv[])
{
  string prog, itags[MaxBodyFields], otags[MaxBodyFields];
  stream xstr, ostr;
  bodyptr btab = NULL;
  int nbody, nout, nold = -1;
  real tnow;

  initparam(argv, defv);
  exprs[0] = getparam("sieve");
  prog = tempnam("/tmp", "sm");
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog);			// start map process running
  if (get_tag_ok(xstr, "History"))
    skip_item(xstr);
  get_history(xstr);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  new_field(&SieveField, BoolType, "Sieve");
  new_field(&SieveField + 1, NULL, NULL);
  while (get_snap(xstr, &btab, &nbody, &tnow, itags, TRUE)) {
    snapsift(btab, &nout, btab, nbody);	// can sift in place b/c map process
					// writes ALL fields of each snapshot
    del_tag(otags, itags, "Sieve");
    if (nout != nold)
      eprintf("[%s: writing %d bodies at time %f]\n", getprog(), nout, tnow);
    put_snap(ostr, &btab, &nout, &tnow, otags);
    nold = nout;
  }
  strclose(ostr);
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getprog(), prog);
  return (0);
}

void snapsift(bodyptr bout, int *nout, bodyptr btab, int nbody)
{
  bodyptr src = btab, dst = bout;

  *nout = 0;
  while (src < NthBody(btab, nbody)) {
    if (Sieve(src)) {
      memcpy(dst, src, SizeofBody);
      dst = NextBody(dst);
      (*nout)++;
    }
    src = NextBody(src);
  }
}

//  execmap: start snapmap subprocess, and return snapmap output stream.
//  ____________________________________________________________________

stream execmap(string prog)
{
  int handle[2];
  char handbuf[32], produce[512];

  pipe(handle);
  if (fork() == 0) {				// if this is child process
    close(handle[0]);
    sprintf(handbuf, "-%d", handle[1]);
    sprintf(produce, "%s,Sieve", getparam("produce"));
    execl(prog, getprog(), getparam("in"), handbuf, getparam("times"),
	  getparam("require"), produce, getparam("passall"),
	  getparam("seed"), NULL);
    error("%s: execl %s failed\n", getprog(), prog);
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
