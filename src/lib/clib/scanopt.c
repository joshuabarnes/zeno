/*
 * scanopt.c: scan string of the form "word1,word2,..." for match. Warning:
 * words must be separated by exactly one comma -- no spaces allowed!
 */

#include "stdinc.h"

bool scanopt(string opt, string key)
{
  char *op, *kp;

  op = (char *) opt;				// start scan of options
  while (*op != (char) NULL) {			// loop over words
    kp = key;
    while ((*op != ',' ? *op : (char) NULL) == *kp) {
						// compare with key
      if (*kp++ == (char) NULL)
	return (TRUE);				// found keyword
      op++;
    }
    while (*op != (char) NULL && *op++ != ',')	// scan for next word
      continue;
  }
  return (FALSE);				// keyword not found
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
    "opt=foo,bar,fum",
    "key=foo",
    NULL,
};

void main(int argc, string argv[])
{
  string opt, key;

  initparam(argv, defv);
  opt = getparam("opt");
  key = getparam("key");
  printf("scanopt(\"%s\", \"%s\") returns %s\n",
	 opt, key, scanopt(opt, key) ? "true" : "false");
}

#endif
