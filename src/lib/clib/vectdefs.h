/*
 * vectdefs.h: definitions from vectmath.h, separated for programs which
 * need to define vectors without loading the whole mess of definitions.
 */

#ifndef _vectdefs_h
#define _vectdefs_h

#if !defined(NDIM) && !defined(TWODIM) && !defined(THREEDIM)
#define THREEDIM			      // specify default dimensions
#endif

#if defined(THREEDIM) || (NDIM==3)
#undef  TWODIM
#define NDIM 3
#define THREEDIM
#endif

#if defined(TWODIM) || (NDIM==2)
#undef  THREEDIM
#define NDIM 2
#define TWODIM
#endif

typedef real vector[NDIM];
typedef real matrix[NDIM][NDIM];

//  ____________________________________________________________________
//  Scalar-valued vector functions.  These pass the number of dimensions
//  as the last argument so they will work with any value of NDIM.

#define dotvp(v,u) (_dotvp((real*)(v),(real*)(u),NDIM))
#define absv(v)    (_absv((real*)(v),NDIM))
#define distv(v,u) (_distv((real*)(v),(real*)(u),NDIM))
#define tracem(p)  (_tracem((real*)(p),NDIM))

#if defined(MIXEDPREC)
#define _dotvp  _mdotvp
#define _absv   _mabsv
#define _distv  _mdistv
#define _tracem _mtracem
#endif

#if defined(SINGLEPREC)
#define _dotvp  _fdotvp
#define _absv   _fabsv
#define _distv  _fdistv
#define _tracem _ftracem
#endif

#if defined(DOUBLEDPREC)
#define _dotvp  _ddotvp
#define _absv   _dabsv
#define _distv  _ddistv
#define _tracem _dtracem
#endif

real _dotvp(real *, real *, int);
real _absv(real *, int);
real _distv(real *, real *, int);
real _tracem(real *, int);

//  ______________________________________________________
//  VectType: Datatypes-style type for 3-D or 2-D vectors.

#ifdef THREEDIM
#define VectType  (RealType RealType RealType)
#else
#define VectType  (RealType RealType)
#endif

#endif  // ! _vectdefs_h
