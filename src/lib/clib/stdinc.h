/* 
 * stdinc.h: standard include file for Zeno C programs.
 * Copyright (c) 2017  Josh Barnes  Honolulu, HI.
 */

#ifndef _stdinc_h
#define _stdinc_h

//  Always include stdio.h, stdlib.h, and string.h.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//  NULL: value for null pointers, normally defined by stdio.h.
//  ___________________________________________________________

#if !defined(NULL)
#  define NULL 0L
#endif

//  local: synonym for static declares an object as local to a source file.
//  _______________________________________________________________________

#define local static

//  bool, TRUE, FALSE: standard names for logical values.
//  _____________________________________________________

typedef short int bool;

#if !defined(TRUE)
#  define TRUE  ((bool) 1)
#  define FALSE ((bool) 0)
#endif

//  byte: name a handy chunk of bits.
//  _________________________________

typedef unsigned char byte;

//  string: null-terminated char array.
//  ___________________________________

typedef char *string;

//  stream: more elegant synonym for FILE *.
//  ________________________________________

typedef FILE *stream;			// note: stdio.h is included above

//  real, realptr: Compile-time precision specification.  Options are:
//      DOUBLEPREC:     everything (variables & functions) is double.
//      MIXEDPREC:      user values are float, -lm functions are double.
//      SINGLEPREC:     everything (variables & functions) is float.
//  See <mathfns.h> for a list of real-valued functions.  If single
//  precision library functions are not availible then use MIXEDPREC
//  instead of SINGLEPREC.
//  ____________________________________________________________________

//  Default precision is SINGLEPREC on LINUX and SGI, and MIXEDPREC on Sun.

#if !defined(MIXEDPREC) && !defined(SINGLEPREC) && !defined(DOUBLEPREC)
#  if !defined(SUN)
#    define SINGLEPREC
#  else
#    define MIXEDPREC
#  endif
#endif

#if defined(DOUBLEPREC)
#  undef SINGLEPREC
#  undef MIXEDPREC
   typedef double real, *realptr;
#  define Precision "DOUBLEPREC"
#  define REALFMT "%lf"
#endif

#if defined(MIXEDPREC)
#  undef DOUBLEPREC
#  undef SINGLEPREC
   typedef float *realptr, real;
#  define Precision "MIXEDPREC"
#  define REALFMT "%f"
#endif

#if defined(SINGLEPREC)
#  undef DOUBLEPREC
#  undef MIXEDPREC
   typedef float real, *realptr;
#  define Precision "SINGLEPREC"
#  define REALFMT "%f"
#endif

//  PI, etc.: mathematical constants.
//  _________________________________

#define PI         3.14159265358979323846
#define TWO_PI     6.28318530717958647693
#define FOUR_PI   12.56637061435917295385
#define HALF_PI    1.57079632679489661923
#define FRTHRD_PI  4.18879020478639098462

//  streq, strne: string-equality macros. strnull: test empty string.
//  Note that string.h should be included before these are used.
//  _________________________________________________________________

#define streq(x,y) (strcmp((x), (y)) == 0)
#define strne(x,y) (strcmp((x), (y)) != 0)
#define strnull(x) (strcmp((x), "") == 0)

//  ABS: returns the absolute value of its argument.
//  MAX: returns the argument with the highest value.
//  MIN: returns the argument with the lowest value.
//  _________________________________________________

#define ABS(x)   (((x)<0)?-(x):(x))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

//  Prototypes for misc. functions in Clib.

void *allocate(size_t nbyte);		// alloc, zero, & check for errors

string *burststring(string, string);	// burst string into token list

double cputime(void);			// returns CPU time in minutes

void set_error_stream(stream);		// send error msgs to given stream
void error(string, ...);		// complain about error and exit
void fatal(string, ...);		// complain about error and abort
void eprintf(string, ...);		// print message to stderr	

void *getxstr(stream, int);		// read extended string
bool putxstr(stream, void *, int);	// write extended string	
void *copxstr(void *, int);		// make copy of extended string
int xstrlen(void *, int);		// find length of extended string
bool xstreq(void *, void *, int);	// compare extended strings

bool scanopt(string, string);		// scan options for keyword

stream stropen(string, string);		// arguments are much like fopen

void get_history(stream instr);		// read history data from stream
void put_history(stream outstr);	// write history data to stream
void add_history(string str);		// append item to history data
void skip_history(stream instr);
string *ask_history(void);		// return vector of history entries

bool within(double, string, double);	// test if value lies in range(s)

#endif  // ! _stdinc_h
