/*
 * burststring.c: break a string of the form "word1<sep>word2<sep>..."
 * into separate strings "word1", "word2", ... and return them in an
 * xstr of strings (ie, NULL-terminated sequence of string pointers).
 * Words may have length zero if two separator characters are adjacent
 * (or if the string begins or ends with a separator).  However, if a
 * blank space counts as a separator, then zero-length words are not
 * returned; thus, if str = "spam, eggs" and sep = ", ", the result is
 * ["spam", "eggs", NULL] instead of ["spam", "", "eggs", NULL].  This
 * seems more likely to yield the right behavior in most applications.
 */

#include "stdinc.h"
#include "getparam.h"
#include <string.h>

string *burststring(string str, string sep)
{
  char *sp, *fp;
  int nsep;
  string *wrdlst, *wp;
  bool keepnull = (strchr(sep, ' ') == NULL);

  nsep = 0;					// init separator count 
  for (sp = str; *sp != (char) NULL; sp++)	// loop over string
    if (strchr(sep, *sp) != NULL)		// identify separator chars
      nsep++;					// and count them up
  wp = wrdlst = (string *) allocate((nsep + 2) * sizeof(string));
						// reserve adequate storage
						// (some may not be used...)
  for (fp = sp = str; *fp != (char) NULL; sp++)	// scan along string
    if (strchr(sep, *sp) != NULL) {		// found next separator
						// (or end, if *sp == NULL)
      if (sp > fp || keepnull) {		// check length of word
	*wp = (string) allocate(MAX(sp - fp, 1) * sizeof(char));
						// get room for word
	(void) strncpy(*wp++, fp, sp - fp);	// make copy of text
      }
      fp = (*sp == (char) NULL ? sp : sp + 1);	// advance to next word
    }
  if ((strnull(str) || strchr(sep, *(fp - 1)) != NULL) && keepnull)
						// if last char was a sep
    *wp = (string) allocate(1 * sizeof(char));	// tack on an empty word
  return (wrdlst);
}

#ifdef TESTBED

string defv[] = {
    "lst=foo,bar,spam;eggs",
    "sep=,;",
    NULL,
};

int main(int argc, string argv[])
{
  string getparam(), lst, sep, *wrds;

  initparam(argv, defv);
  lst = getparam("lst");
  sep = getparam("sep");
  wrds = burststring(lst, sep);
  while (*wrds != NULL)
    printf("\"%s\"  ", *wrds++);
  printf("\n");
  return (0);
}

#endif
