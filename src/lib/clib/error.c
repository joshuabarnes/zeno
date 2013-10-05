/* 
 * error.c: routines to report errors, etc.
 */

#include "stdinc.h"
#include <stdarg.h>
#include <string.h>

local stream errstr = NULL;			// use to report errors

//  _______________________________________________________
//  set_error_stream: direct error message to given stream.

void set_error_stream(stream str)
{
  errstr = str;
}

//  ______________________________
//  error: scream and die quickly.

void error(string fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);			// printf msg to std error
  va_end(ap);
  fflush(stderr);
  if (errstr != NULL) {				// error logging enabled
    va_start(ap, fmt);
    vfprintf(errstr, fmt, ap);
    va_end(ap);
    fflush(errstr);
  }
  exit(1);					// quit with error status
}

//  ____________________________________
//  fatal: scream and die an ugly death.

void fatal(string fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);			// printf msg to std error
  va_end(ap);
  fflush(stderr);
  if (errstr == NULL)				// default to error log
    errstr = fopen("fatal_error.log", "a");
  if (errstr != NULL) {				
    va_start(ap, fmt);
    vfprintf(errstr, fmt, ap);			// printf msg to error str
    va_end(ap);
    fflush(errstr);
  }
  abort();					// quit, leave core image
}

//  _______________________________________________________________________
//  eprintf: print messages and warnings.  Uses "ZENO_MSG_OPTION" env. var.
//  to control printing; "warn" or "none" suppress some or all output.


void eprintf(string fmt, ...)
{
  static string msgopt = NULL;
  va_list ap;

  if (msgopt == NULL)				// if NULL, get env. value
    msgopt = getenv("ZENO_MSG_OPTION");
  if ((msgopt == NULL) ||
      (strne(msgopt, "none") &&
       (strne(msgopt, "warn") || 
	(strstr(fmt, "warn") != NULL || strstr(fmt, "WARN") != NULL)))) {
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);			// printf msg to std error
    va_end(ap);
    fflush(stderr);				// drain std error buffer
  }
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
  "logfile=",
  "fatal=false",
  NULL,
};

int main(int argc, string argv[])
{
  initparam(argv, defv);
  if (!strnull(getparam("logfile")))
    set_error_stream(fopen(getparam("logfile"), "w"));
  eprintf("[%s: foo=%f  bar=%d  fum=\"%s\"]\n",
	  getprog(), 2.7183, 16384, "frodo");
  eprintf("[%s: WARNING: foo=%f  bar=%d  fum=\"%s\"]\n",
	  getprog(), 3.1415, 32768, "bilbo");
  eprintf("[%s: warning: foo=%f  bar=%d  fum=\"%s\"]\n",
	  getprog(), 3.1415, 32768, "bilbo");
  if (! getbparam("fatal"))
    error("error: foo=%f  bar=%d  fum=\"%s\"\n", 0.5772, 65536, "waldo");
  else
    fatal("error: foo=%f  bar=%d  fum=\"%s\"\n", 0.5772, 65536, "waldo");
  return (0);
}

#endif
