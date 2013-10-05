/*
 * stropen.c: stream-open enhancement of fopen function.
 */

#include "stdinc.h"
#include "getparam.h"
#include <sys/types.h>
#include <sys/stat.h>

//  ___________________________________________________________________
//  stropen: open a stdio stream like fopen(), with these extensions:
//  1.  error checking: stropen() can only return NULL if mode == "r?";
//  2.  data protection: files can only be clobbered if mode == "w!";
//  3.  names of the form "-" map to stdin/stdout, depending on mode;
//  4.  names of the form "-num" open a stream to read/write f.d. num.

stream stropen(string name, string xmode)
{
    bool readflag, testflag;
    char mode[2];
    int fds;
    stream res;
    struct stat buf;

    if (! (streq(xmode, "r") || streq(xmode, "r?") ||
	     streq(xmode, "w") || streq(xmode, "w!") ||
	       streq(xmode, "a")))
	error("%s.stropen: illegal xmode \"%s\"\n",
	      getprog(), xmode);
    readflag = (xmode[0] == 'r');
    testflag = (xmode[1] == '?');
    mode[0] = xmode[0];
    mode[1] = (char) NULL;
    if (name[0] == '-') {			
	if (streq(name, "-")) {
	    fds = dup(fileno(readflag ? stdin : stdout));
	    if (fds == -1)
		error("%s.stropen: cannot dup %s\n",
		      getprog(), readflag ? "stdin" : "stdout");
	} else if (sscanf(name, "-%d", &fds) != 1)
	    error("%s.stropen: bad f.d. number \"%s\"\n",
		  getprog(), name);
	res = fdopen(fds, mode);
	if (res == NULL && !testflag)
	    error("%s.stropen: cannot open f.d. %d for %s\n",
		  getprog(), fds, readflag ? "input" : "output");
    } else {
	if (streq(xmode, "w") && stat(name, &buf) == 0)
	    error("%s.stropen: file \"%s\" already exists\n",
		  getprog(), name);
	res = fopen(name, mode);
	if (res == NULL && !testflag)
	    error("%s.stropen: cannot open file \"%s\" for %s\n",
		  getprog(), name, readflag ? "input" : "output");
    }
    return (res);
}

#ifdef TESTBED

string defv[] = {
    "name=foo.bar",
    "mode=w",
    "text=boo hoo foo",
    NULL,
};

void main(int argc, string argv[])
{
    string name, mode, text;
    stream str;
    char buf[128];

    initparam(argv, defv);
    name = getparam("name");
    mode = getparam("mode");
    text = getparam("text");
    str = stropen(name, mode);
    if (str == NULL)
	error("%s: stropen() returned NULL\n", getprog());
    if (streq(mode, "r") || streq(mode, "r?")) {
	while (fgets(buf, 127, str) != NULL)
	    printf("%s", buf);
    } else {
	sprintf(buf, "%s\n", text);
	fputs(buf, str);
    }
    fclose(str);
}

#endif
