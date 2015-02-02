/*
 * snapsort.c: sort bodies by value.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "filestruct.h"
#include "phatbody.h"
#include "buildmap.h"
#include <unistd.h>

string defv[] = {               ";Sort bodies by computed value",
  "in=???",                     ";Input snapshot file name",
  "out=???",                    ";Output snapshot file name",
  "times=all",                  ";Range of times to process",
  "value=???",			";Expression (C code) for sort value.",
				";May use these values (if given in input):",
				  SNAPMAP_BODY_VARS ".",
  "require=",			";Input items required",
  "produce=",			";Output items produced",
  "passall=true",		";If true, pass on input data",
  "seed=",			";Generator seed for random values",
  "VERSION=2.1",                ";Josh Barnes  2 February 2015",
  NULL,
};

int cmpvalue(const void *, const void *);	// compare values of bodies
stream execmap(string);				// start snapmap process
void del_tag(string *, string *, string);	// remove tag from list

string names[2] = { "Value",  NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { RealType, NULL };

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
  exprs[0] = getparam("value");
  prog = mktemp((string) copxstr("/tmp/sm_XXXXXX", sizeof(char)));
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
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
    error("%s: can't unlink %s\n", getprog(), prog);
  return (0);
}

int cmpvalue(const void *a, const void *b)
{
  return (Value((bodyptr) a) < Value((bodyptr) b) ? -1 :
	  Value((bodyptr) a) > Value((bodyptr) b) ? 1 : 0);
}

//  execmap: start snapmap subprocess, and return snapmap output stream.
//  ____________________________________________________________________

stream execmap(string prog)
{
  int handle[2];
  char handbuf[32], produce[512];

  pipe(handle);
  if (fork() == 0) {                            // if this is child process
    close(handle[0]);
    sprintf(handbuf, "-%d", handle[1]);
    sprintf(produce, "%s,Value", getparam("produce"));
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
