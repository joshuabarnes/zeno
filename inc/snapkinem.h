/*
 * snapkinem.h: routines to compute kinematic parameters.
 */

#if defined(THREEDIM)
#if defined(SINGLEPREC) || defined(MIXEDPREC)
#define snapke     f3snapke
#define snappe     f3snappe
#define snapamvec  f3snapamvec
#define snapketen  f3snapketen
#define snappeten  f3snappeten
#define snapmiten  f3snapmiten
#else
#define snapke     d3snapke
#define snappe     d3snappe
#define snapamvec  d3snapamvec
#define snapketen  d3snapketen
#define snappeten  d3snappeten
#define snapmiten  d3snapmiten
#endif
#endif

real snapke(bodyptr, int, int);

real snappe(bodyptr, int, int);

void snapamvec(vector, bodyptr, int, int);

void snapketen(matrix, bodyptr, int, int);

void snappeten(matrix, bodyptr, int, int);

void snapmiten(matrix, bodyptr, int, int);
