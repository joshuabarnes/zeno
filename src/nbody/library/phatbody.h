/*
 * phatbody.h: Macro definitions and prototypes for phat body structures.
 */

#ifndef _phatbody_h
#define _phatbody_h

#include "strset.h"
#include "phatstruct.h"
#include "bodytags.h"

//  bodyptr: name for pointer to a phat body.
//  _________________________________________

typedef void *bodyptr;

//  layout_body: structure definition routine.
//  __________________________________________

void layout_body(string *, string, int);

//  Snapshot i/o functions.
//  _______________________

#if defined(threedim)
#  if defined(singleprec) || defined(mixedprec)
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

//  phatbody: array of fields for phat body structures.
//  ___________________________________________________

#define MaxBodyFields 32			// total number of fields

extern ps_field phatbody[MaxBodyFields];	// describe phat bodies

//  Define names for standard fields in phat body structure.

#define BodyField      phatbody[ 0]		// 0th field defines bodies
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
#define NewBodyFields		20		// index of 1st free field

//  SizeofBody: number of bytes per body.
//  _____________________________________

#define SizeofBody (BodyField.length)

//  NextBody, PrevBody, NthBody: macros for stepping through body arrays.
//  _____________________________________________________________________

#define NextBody(bp)  ((bodyptr) ((byte *)(bp) + SizeofBody))
#define PrevBody(bp)  ((bodyptr) ((byte *)(bp) - SizeofBody))
#define NthBody(bp,n) ((bodyptr) ((byte *)(bp) + SizeofBody * (n)))

//  SafeOffset: If SafeSelect is defined, referencing an nonexistent
//  body component will produce an intelligible error message, else
//  if not, the result is undefined.
//  ________________________________________________________________

#if defined(SafeSelect)
#  include "getparam.h"
#  define SafeOffset(field) \
	    (field.offset != BadOffset ? field.offset :	\
	      (error("%s: %s undefined\n", getargv0(), field.name), 0))
#else
#  define SafeOffset(field) (field.offset)
#endif

//  Accessor macros for body components.
//  ____________________________________

#define Pos(b)       SelectVect(b, SafeOffset(PosField))
#define PosX(b)      (Pos(b)[0])
#define PosY(b)      (Pos(b)[1])
#define PosZ(b)      (Pos(b)[2])
#define Vel(b)       SelectVect(b, SafeOffset(VelField))
#define VelX(b)      (Vel(b)[0])
#define VelY(b)      (Vel(b)[1])
#define VelZ(b)      (Vel(b)[2])
#define Mass(b)      SelectReal(b, SafeOffset(MassField))
#define Phi(b)       SelectReal(b, SafeOffset(PhiField))
#define Acc(b)       SelectVect(b, SafeOffset(AccField))
#define AccX(b)      (Acc(b)[0])
#define AccY(b)      (Acc(b)[1])
#define AccZ(b)      (Acc(b)[2])
#define Smooth(b)    SelectReal(b, SafeOffset(SmoothField))
#define Rho(b)       SelectReal(b, SafeOffset(RhoField))
#define EntFunc(b)   SelectReal(b, SafeOffset(EntFuncField))
#define Uintern(b)   SelectReal(b, SafeOffset(UinternField))
#define UdotInt(b)   SelectReal(b, SafeOffset(UdotIntField))
#define UdotRad(b)   SelectReal(b, SafeOffset(UdotRadField))
#define UdotVis(b)   SelectReal(b, SafeOffset(UdotVisField))
#define Tau(b)       SelectReal(b, SafeOffset(TauField))
#define Type(b)      SelectByte(b, SafeOffset(TypeField))
#define Birth(b)     SelectReal(b, SafeOffset(BirthField))
#define Death(b)     SelectReal(b, SafeOffset(DeathField))
#define Key(b)       SelectInt(b,  SafeOffset(KeyField))
#define Aux(b)       SelectReal(b, SafeOffset(AuxField))
#define AuxVec(b)    SelectVect(b, SafeOffset(AuxVecField))
#define AuxVecX(b)   (AuxVec(b)[0])
#define AuxVecY(b)   (AuxVec(b)[1])
#define AuxVecZ(b)   (AuxVec(b)[2])

#endif /* ! _phatbody_h */
