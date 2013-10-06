/*
 * extstring.c: an extension of the standard C concept of character
 * strings to strings of n-byte (nonzero) values, terminated by a
 * marker of n zero bytes.
 */

#include "stdinc.h"
#include "getparam.h"

#define MAXLEN	1024				// input buffer length

//  getxstr: read extended string from input stream, return pointer.
//  ________________________________________________________________

void *getxstr(stream inpt, int nbyt)
{
  byte buf[MAXLEN], *bp;
  bool lpflg;
  int i, ch;

  bp = &buf[0];					// set point into buffer
  do {						// loop reading data in
    lpflg = FALSE;				// init loop flag
    for (i = 0; i < nbyt; i++) {		// loop over block of bytes
      ch = getc(inpt);				// input next byte
      if (bp > &buf[MAXLEN-1])			// detect overflow error
	error("%s.getxstr: buffer overflow\n", getprog());
      *bp = ch != EOF ? ch : (int) NULL;	// map EOF to NULL
      if (*bp++ != (byte) NULL)			// got a byte of real data?
	lpflg = TRUE;				// repeat for next block
    }
  } while (lpflg);				// until a block of NULLs
  return (copxstr(&buf[0], nbyt));		// make copy and return it
}

//  putxstr: write extended string to output stream, return TRUE on success.
//  ________________________________________________________________________

bool putxstr(stream outp, void *xspt, int nbyt)
{
  byte *bp, c;
  int n;

  bp = (byte *) xspt;				// init pointer to string
  n = nbyt * xstrlen(xspt, nbyt);		// get length in bytes
  while (--n >= 0) {				// loop writing out bytes
    c = *bp++;					// get byte to output
    putc(0377 & c, outp);			// and write eight bits out
    if (ferror(outp))				// did output call fail?
      return (FALSE);				// then so did we
  }
  return (TRUE);				// return sign of success
}

//  copxstr: make copy of extended string in allocated memory, return pointer.
//  __________________________________________________________________________

void *copxstr(void *xspt, int nbyt)
{
  byte *sp, *dp, *dest;
  int n;

  sp = (byte *) xspt;				// init pointer to string
  n = nbyt * xstrlen(xspt, nbyt);		// get length in bytes
  dp = dest = (byte *) allocate(n);		// allocate new storage
  while (--n >= 0)				// loop over bytes
    *dp++ = *sp++;				// copy each in turn
  return (dest);				// return copy string
}

//  xstrlen: count number of values (including null) in extended string.
//  ____________________________________________________________________

int xstrlen(void *xspt, int nbyt)
{
  byte *bp;
  int nval, i;
  bool lpflg;

  bp = (byte *) xspt;				// init pointer to string
  nval = 0;					// init count of values
  do {						// loop over values
    nval++;					// count one more
    lpflg = FALSE;				// init loop flag
    for (i = 0; i < nbyt; i++)			// loop over bytes
      if (*bp++ != (byte) NULL)			// a byte of data?
	lpflg = TRUE;				// set loop flag
  } while (lpflg);				// until a NULL value
  return (nval);				// return total count
}

//  xstreq: determine if extended strings are equal, return TRUE if so.
//  ___________________________________________________________________

bool xstreq(void *xp1, void *xp2, int nbyt)
{
  byte *bp1, *bp2;
  int n;

  bp1 = (byte *) xp1;				// init pointers to strings
  bp2 = (byte *) xp2;
  n = nbyt * xstrlen(xp1, nbyt);		// get length in bytes
  while (--n >= 0)				// loop over bytes
    if (*bp1++ != *bp2++)			// bytes not equal?
      return (FALSE);				// then strs unequal
  return (TRUE);				// indicate equality
}

#ifdef TESTBED

main(int argc, string argv[])
{
  int i;
  long lstr[32], *lcop;
  stream opt, ipt;

  for (i = 0; i < 32; i++)
    lstr[i] = (i < 21 ? 12345 + 512 * i : (int) NULL);
  printf("xstrlen(lstr, %d) == %d\n",
	 sizeof(long), xstrlen(lstr, sizeof(long)));
  lcop = (long *) copxstr(lstr, sizeof(long));
  printf("xstrlen(lcop, %d) == %d\n",
	 sizeof(long), xstrlen(lcop, sizeof(long)));
  printf("xstreq(lstr, lcop, %d) == %d\n",
	 sizeof(long), xstreq(lstr, lcop, sizeof(long)));
  printf("changing *lcop\n");
  *lcop = -1;
  printf("xstreq(lstr, lcop, %d) == %d\n",
	 sizeof(long), xstreq(lstr, lcop, sizeof(long)));
  opt = stropen("foobar.dat", "w!");
  printf("putxstr(opt, lstr, %d) == %d\n",
	 sizeof(long), putxstr(opt, lstr, sizeof(long)));
  printf("putxstr(opt, lcop, %d) == %d\n",
	 sizeof(long), putxstr(opt, lcop, sizeof(long)));
  fclose(opt);
  ipt = stropen("foobar.dat", "r");
  lcop = (long *) getxstr(ipt, sizeof(long));
  if (lcop == NULL)
    printf("getxstr(ipt, %d) failed\n", sizeof(long));
  printf("xstreq(lstr, lcop, %d) == %d\n",
	 sizeof(long), xstreq(lstr, lcop, sizeof(long)));
}

#endif
