/****************************************************************************/
/* SMOOTH.H: include file for SPH smoothing routines in smooth.c.           */
/* Copyright (c) 2000 by Joshua E. Barnes, Santa Barbara, CA.               */
/* Adapted from SMOOTH V2.01; see http://www-hpcc.astro.washington.edu/     */
/****************************************************************************/

#ifndef _smooth_h
#define _smooth_h

/* 
 * Misc. constants.
 */

#define MAXLEV   20			/* max levels in timestep hierarchy */
#define NCOEFS    5			/* coefficients in smoothing func   */
#define NRPROF   11			/* bins in smoothing ball profile   */
#define EXTLIST  10                     /* extra room on neighbor lists     */

/*
 * PQNODE: item in priority que.
 */

typedef struct _pqnode {
    real pqkey;
    int pind;
    struct _pqnode *pqLoser;
    struct _pqnode *pqFromInt;
    struct _pqnode *pqFromExt;
    struct _pqnode *pqWinner;
} pqnode;

/*
 * SMCONTEXT, SMXPTR: context for smoothing process.
 */

typedef struct {
    kdxptr kd;				/* KD context from init_kdtree	    */
    int  nsmooth;			/* bodies within smoothing sphere   */
    real coefs[NCOEFS];			/* coefficients for SPH kernel	    */
    int  rprof[NRPROF];			/* counts of bodies within sphere   */
    real freqmax;			/* maximum integration frequency    */
    real freqsum;			/* sum of integration frequencies   */
    pqnode *pque;			/* priority que for neighbor search */
    pqnode *pqhead;			/* most distant neighbor in que     */
    real *r2ball;			/* array of smoothing radii	    */
    real *r2list;			/* list of neighbor distances^2     */
    int *inlist;			/* list of neighbor indicies        */
    real cpustart;			/* CPU time at start of smoothing   */
} smcontext, *smxptr;

/*
 * WSMOOTH, DWSMOOTH: macros to compute kernel and gradient.
 */

#define WSmooth(wsm, wsc, x2, cf)					\
{									\
    real _x = rsqrt(x2);                                                \
    if (x2 < 1)                                                         \
        wsm = wsc * (cf[0] + cf[1]*_x + cf[2]*x2 + cf[3]*_x*x2);        \
    else                                                                \
        wsm = wsc * cf[4] * (2-_x)*(2-_x)*(2-_x);                       \
}

#define dWSmooth(dsm, dsc, x2, cf)					\
{					                                \
    real _x = rsqrt(x2);                                                \
    if (x2 < 1)                                                         \
        dsm = dsc * (cf[1] + 2*cf[2]*_x + 3*cf[3]*x2);                  \
    else                                                                \
        dsm = dsc * cf[4] * -3 * (2-_x)*(2-_x);                         \
}

/*
 * Smoothing procedure prototypes.
 */

smxptr init_smooth(kdxptr, int, real);
void smooth(smxptr, void (*)(smxptr, int, int));
void resmooth(smxptr, void (*)(smxptr, int, int));
void report_smooth(smxptr, stream, string);
void finish_smooth(smxptr);

#endif  /* ! _smooth_h */
