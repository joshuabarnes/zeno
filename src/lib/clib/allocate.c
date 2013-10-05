/*
 * allocate.c: memory allocation and zeroing, with error checking.
 */

#include "stdinc.h"
#include "getparam.h"

void *allocate(int nb)
{
  void *mem;

  if (nb <= 0)
    error("%s.allocate: absurd request (nb = %d)\n", getprog(), nb);
  mem = calloc(nb, 1);			// use calloc to zero mem
  if (mem == NULL)
    error("%s.allocate: not enough memory (nb = %d)\n", getprog(), nb);
  return (mem);
}
