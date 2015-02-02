/*
 * burststring.c: break a string of the form "word1, word2, ..." into
 * separate strings "word1", "word2", ... and return them in an
 * extended-string (ie, NULL-terminated sequence of pointers).
 */

#include "stdinc.h"
#include "getparam.h"
#include <string.h>

string *burststring(string lst, string sep)
{
  char *lp, *fp;
  int nwrd = 0;
  string *wlst, *wp;
  
  for (lp = lst, fp = NULL; *lp != NULL; lp++)	// count words in list
    if (strchr(sep, *lp) == NULL) {		// part of a word?
      if (fp == NULL) {				// and first char?
	fp = lp;				// save ptr to start
	nwrd++;					// count another word
      }
    } else					// not part of word
      fp = NULL;				// change scan state
  wp = wlst = (string *) allocate((nwrd + 1) * sizeof(string));
						// get space for list
  for (lp = lst, fp = NULL; *lp != NULL; lp++)	// store words in list
    if (strchr(sep, *lp) == NULL) {		// part of a word?
      if (fp == NULL)				// and first char?
	fp = lp;				// save ptr to start
    } else {					// not part of word
      if (fp != NULL) {				// but prev char was
	*wp = (char *) allocate((1 + lp - fp) * sizeof(char));
						// get space for word
	(void) strncpy(*wp, fp, lp - fp);	// copy word to space
	(*wp++)[lp - fp] = NULL;		// make sure it ended
      }
      fp = NULL;				// change scan state
    }
  if (fp != NULL) {				// handle final word
    *wp = (char *) allocate((1 + lp - fp) * sizeof(char));
    (void) strncpy(*wp, fp, lp - fp);
    (*wp++)[lp - fp] = NULL;
  }
  *wp = NULL;					// terminate word list
  return (wlst);
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
