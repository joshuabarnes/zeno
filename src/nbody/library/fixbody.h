/*
 * FIXBODY.H: Macro definitions and function prototypes for
 * fixed-offset interface to phat body structures.
 */

#ifndef _fixbody_h
#define _fixbody_h

#include "strset.h"
#include "phatstruct.h"
#include "bodytags.h"

/*
 * NEXTBODY, PREVBODY, NTHBODY: macros for stepping through body arrays.
 */

#define NextBody(bp)   ((bp) + 1)
#define PrevBody(bp)   ((bp) - 1)
#define NthBody(bp,n)  ((bp) + (n))

/*
 * BODYOFFSET: offset of named field.
 */

#define BodyOffset(f)  ((int) &f((bodyptr) 0))

/*
 * DEFINE_BODY, DEFINE_BODY_OFFSET: structure definition routines.
 */

void define_body(int, string, int);
void define_body_offset(string, int);

/*
 * Snapshot I/O functions.
 */

#if defined(THREEDIM)
#  if defined(SINGLEPREC) || defined(MIXEDPREC)
#    define put_snap    f3put_snap
#    define get_snap    f3get_snap
#    define get_snap_t  f3get_snap_t
#  else
#    define put_snap    d3put_snap
#    define get_snap    d3get_snap
#    define get_snap_t  d3get_snap_t
#  endif
#endif

void put_snap(stream, bodyptr *, int *, real *, string *);
bool get_snap(stream, bodyptr *, int *, real *, string *, bool);
bool get_snap_t(stream, bodyptr *, int *, real *, string *, bool, string);

/*
 * PHATBODY: array of phat structure fields for body structures.
 */

#define MaxBodyFields  25			/* total number of fields   */
#define NewBodyFields  16			/* index of 1st free field  */

extern ps_field phatbody[MaxBodyFields];	/* describe phat bodies     */

#endif /* ! _fixbody_h */
