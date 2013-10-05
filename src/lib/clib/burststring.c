/*
 * burststring.c: break a string of the form "word1, word2, ..." into
 * seperate strings "word1", "word2", ... and return them in an
 * extended-string (ie, NULL-terminated sequence of pointers).
 */

#include "stdinc.h"
#include "getparam.h"
#include <string.h>

#define MWRD  64	// max words in list
#define MSTR 256	// max chars per word

string *burststring(string lst, string sep)
{
  string wrdbuf[MWRD], *wp;
  char strbuf[MSTR], *sp, *lp;

  wp = wrdbuf;
  sp = strbuf;
  lp = lst;
  do {							// scan over string
    if (*lp == (char) NULL ||
	  strchr(sep, *lp) != NULL) {			// is this a sep?
      if (sp > strbuf) {				// and got a word?
	*sp = (char) NULL;
	*wp++ = (string) copxstr(strbuf, sizeof(char));
	if (wp == &wrdbuf[MWRD])			// no room in buf?
	  error("%s.burststring: too many words\n", getprog());
	sp = strbuf;					// ready for next
      }
    } else {						// part of word
      *sp++ = *lp;					// so copy it over
      if (sp == &strbuf[MSTR])				// no room left?
	error("%s.burststring: word too long\n", getprog());
    }
  } while (*lp++ != (char) NULL);			// until list ends
  *wp = NULL;						// end word list
  return ((string *) copxstr(wrdbuf, sizeof(string)));
}

#ifdef TESTBED

string defv[] = {
    "lst=foo, bar,waldo ,",
    "sep= ,",
    NULL,
};

main(argc, argv)
int argc;
string argv[];
{
    string getparam(), lst, sep, *wrds;

    initparam(argv, defv);
    lst = getparam("lst");
    sep = getparam("sep");
    wrds = burststring(lst, sep);
    while (*wrds != NULL)
	printf("\"%s\"  ", *wrds++);
    printf("\n");
}

#endif
