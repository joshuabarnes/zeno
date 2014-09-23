/*
 * snapcenter.h: routines to compute center-of-mass coordinates.
 */

#if defined(THREEDIM)
#if defined(SINGLEPREC) || defined(MIXEDPREC)
#define snapcmpos   f3snapcmpos
#define snapcmvel   f3snapcmvel
#define snapcmacc   f3snapcmacc
#define snapcenter  f3snapcenter
#else
#define snapcmpos   d3snapcmpos
#define snapcmvel   d3snapcmvel
#define snapcmacc   d3snapcmacc
#define snapcenter  d3snapcenter
#endif
#endif

void snapcmpos(vector, bodyptr, int, int);

void snapcmvel(vector, bodyptr, int, int);

void snapcenter(bodyptr, int, int);
