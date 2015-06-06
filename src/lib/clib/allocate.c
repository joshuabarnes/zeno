/*
 * allocate.c: memory allocation and zeroing, with error checking.
 */

#include "stdinc.h"
#include "getparam.h"

void *allocate(size_t nbyte)
{
  void *mem;

  mem = calloc(nbyte, 1);			// use calloc to zero mem
  if (mem == NULL)
    error("%s.allocate: calloc failed (nbyte = %lu)\n", getprog(), nbyte);
  return (mem);
}
