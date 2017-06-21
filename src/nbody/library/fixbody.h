/*
 * fixbody.h: Macro definitions and function prototypes for
 * fixed-offset interface to phat body structures.
 */

#ifndef _fixbody_h
#define _fixbody_h

#include "strset.h"
#include "phatstruct.h"
#include "bodytags.h"

//  nextbody, prevbody, nthbody: macros for stepping through body arrays.
//  _____________________________________________________________________

#define NextBody(bp)   ((bp) + 1)
#define PrevBody(bp)   ((bp) - 1)
#define NthBody(bp,n)  ((bp) + (n))

//  bodyoffset: offset of named field.
//  __________________________________

#define BodyOffset(f)  ((int) &f((bodyptr) 0))

//  define_body, define_body_offset: structure definition routines.
//  _______________________________________________________________

void define_body(int, string, int);
void define_body_offset(string, int);

//  Snapshot I/O functions; actual name depends on precision.
//  _________________________________________________________

#if defined(THREEDIM)
#  if defined(SINGLEPREC) || defined(MIXEDPREC)
#    define put_snap    f3put_snap
#    define get_snap    f3get_snap
#  else
#    define put_snap    d3put_snap
#    define get_snap    d3get_snap
#  endif
#endif

void put_snap(stream, bodyptr *, int *, real *, string *);
bool get_snap(stream, bodyptr *, int *, real *, string *, bool, string);

//  phatbody: array of phat structure fields for body structures.
//  _____________________________________________________________

#define MaxBodyFields  32			// total number of fields
#define NewBodyFields  22			// index of 1st free field

extern ps_field phatbody[MaxBodyFields];	// describe phat bodies

#endif /* ! _fixbody_h */
