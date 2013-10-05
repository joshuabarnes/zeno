/*
 * PHATBODY.H: Macro definitions and prototypes for phat body structures.
 */

#ifndef _phatbody_h
#define _phatbody_h

#include "strset.h"
#include "phatstruct.h"
#include "bodytags.h"

/*
 * BODYPTR: name for pointer to a phat body.
 */

typedef void *bodyptr;

/*
 * LAYOUT_BODY: structure definition routine.
 */

void layout_body(string *, string, int);

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
 * PHATBODY: array of fields for phat body structures.
 */

/* Define names for standard fields in phat body structure. */

#define BodyField      phatbody[ 0]		/* 0th field defines bodies */
#define PosField       phatbody[ 1]
#define VelField       phatbody[ 2]
#define MassField      phatbody[ 3]
#define PhiField       phatbody[ 4]
#define AccField       phatbody[ 5]
#define SmoothField    phatbody[ 6]
#define RhoField       phatbody[ 7]
#define EntFuncField   phatbody[ 8]
#define UinternField   phatbody[ 9]
#define UdotIntField   phatbody[10]
#define UdotRadField   phatbody[11]
#define UdotVisField   phatbody[12]
#define TauField       phatbody[13]
#define TypeField      phatbody[14]
#define BirthField     phatbody[15]
#define DeathField     phatbody[16]
#define KeyField       phatbody[17]
#define AuxField       phatbody[18]
#define AuxVecField    phatbody[19]

#define NewBodyFields 20			/* index of 1st free field  */
#define MaxBodyFields 32			/* total number of fields   */

extern ps_field phatbody[MaxBodyFields];	/* describe phat bodies     */

/*
 * SIZEOFBODY: number of bytes per body.
 */

#define SizeofBody (BodyField.length)

/*
 * NEXTBODY, PREVBODY, NTHBODY: macros for stepping through body arrays.
 */

#define NextBody(bp)  ((bodyptr) ((byte *)(bp) + SizeofBody))
#define PrevBody(bp)  ((bodyptr) ((byte *)(bp) - SizeofBody))
#define NthBody(bp,n) ((bodyptr) ((byte *)(bp) + SizeofBody * (n)))

/*
 * GoodOffset: Check offset before using, if SafeSelect is defined.
 * If it is, referencing an nonexistent body component will produce
 * an intelligible error message; if not, the result is undefined.
 */

#if defined(SafeSelect)
#  include "getparam.h"
#  define GoodOffset(field) \
	    (field.offset != BadOffset ? field.offset : \
	       (error("%s: %s undefined\n", getargv0(), field.name), 0))
#else
#  define GoodOffset(field) (field.offset)
#endif

/*
 * Accessor macros for body components.
 */

#define Pos(b)       SelectVect(b, GoodOffset(PosField))
#define PosX(b)      (Pos(b)[0])
#define PosY(b)      (Pos(b)[1])
#define PosZ(b)      (Pos(b)[2])
#define Vel(b)       SelectVect(b, GoodOffset(VelField))
#define VelX(b)      (Vel(b)[0])
#define VelY(b)      (Vel(b)[1])
#define VelZ(b)      (Vel(b)[2])
#define Mass(b)      SelectReal(b, GoodOffset(MassField))
#define Phi(b)       SelectReal(b, GoodOffset(PhiField))
#define Acc(b)       SelectVect(b, GoodOffset(AccField))
#define AccX(b)      (Acc(b)[0])
#define AccY(b)      (Acc(b)[1])
#define AccZ(b)      (Acc(b)[2])
#define Smooth(b)    SelectReal(b, GoodOffset(SmoothField))
#define Rho(b)       SelectReal(b, GoodOffset(RhoField))
#define EntFunc(b)   SelectReal(b, GoodOffset(EntFuncField))
#define Uintern(b)   SelectReal(b, GoodOffset(UinternField))
#define UdotInt(b)   SelectReal(b, GoodOffset(UdotIntField))
#define UdotRad(b)   SelectReal(b, GoodOffset(UdotRadField))
#define UdotVis(b)   SelectReal(b, GoodOffset(UdotVisField))
#define Tau(b)       SelectReal(b, GoodOffset(TauField))
#define Type(b)      SelectByte(b, GoodOffset(TypeField))
#define Birth(b)     SelectReal(b, GoodOffset(BirthField))
#define Death(b)     SelectReal(b, GoodOffset(DeathField))
#define Key(b)       SelectInt(b,  GoodOffset(KeyField))
#define Aux(b)       SelectReal(b, GoodOffset(AuxField))
#define AuxVec(b)    SelectVect(b, GoodOffset(AuxVecField))
#define AuxVecX(b)   (AuxVec(b)[0])
#define AuxVecY(b)   (AuxVec(b)[1])
#define AuxVecZ(b)   (AuxVec(b)[2])

#endif /* ! _phatbody_h */
