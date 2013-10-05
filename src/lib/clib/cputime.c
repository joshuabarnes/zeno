/*
 * cputime.c: compute total process CPU time in minutes.
 */

#include "stdinc.h"
#include "getparam.h"

#if !defined(OBSOLETE_CODE)

#include <sys/time.h>
#include <sys/resource.h>

double cputime(void)
{
  struct rusage buf;

  if (getrusage(RUSAGE_SELF, &buf) == -1)
    error("%s.cputime: getrusage() call failed\n", getprog());
  return ((buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1000000.0 +
	   buf.ru_stime.tv_sec + buf.ru_stime.tv_usec / 1000000.0) / 60.0);
}

#else

#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>

double cputime(void)
{
    struct tms buffer;

    if (times(&buffer) == -1)
	error("%s.cputime: times() call failed\n", getprog());
    return ((buffer.tms_utime + buffer.tms_stime) / (60.0 * HZ));
}

#endif
