/****************************************************************************/
/* SMOOTH.C: routines to perform smoothing for SPH.  Public routines:       */
/* init_smooth, finish_smooth, density, derivatives, report_smooth.         */
/* Copyright (c) 2000 by Joshua E. Barnes, Santa Barbara, CA.               */
/* Adapted from SMOOTH V2.01; see http://www-hpcc.astro.washington.edu/     */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "sphdefs.h"
#include "kdtree.h"
#include "smooth.h"

/*
 * PQINIT, PQBUILD, PQREPLACE: macros to manipulate priority queues.
 */

#define PQInit(pq, n)							\
{									\
    int _j;								\
    for (_j = 0; _j < (n); _j++) {					\
	(pq)[_j].pqFromInt = (_j < 2 ? NULL : &(pq)[_j >> 1]);		\
	(pq)[_j].pqFromExt = &(pq)[(_j+(n)) >> 1];			\
    }									\
}

#define PQBuild(pq, n, q)						\
{									\
    int _j, _i;								\
    pqnode *_n1, *_n2;							\
    for (_j = (n) - 1; _j>0; _j--) {					\
        _i = (_j << 1);							\
        _n1 = (_i < (n) ? (pq)[_i].pqWinner : &(pq)[_i-(n)]);		\
        _i++;								\
        _n2 = (_i < (n) ? (pq)[_i].pqWinner : &(pq)[_i-(n)]);		\
        (pq)[_j].pqLoser = (_n1->pqkey < _n2->pqkey ? _n1 : _n2);	\
        (pq)[_j].pqWinner = (_n1->pqkey < _n2->pqkey ? _n2 : _n1);	\
    }									\
    (q) = (pq)[1].pqWinner;						\
}

#define PQReplace(q)							\
{									\
    pqnode *_n1, *_n2;							\
    _n1 = (q)->pqFromExt;						\
    while (_n1) {							\
        if (_n1->pqLoser->pqkey > (q)->pqkey) {				\
            _n2 = _n1->pqLoser;						\
            _n1->pqLoser = (q);						\
            (q) = _n2;							\
        }								\
        _n1 = _n1->pqFromInt;						\
    }									\
}

/*
 * INTERSECT: macro to determine if node intersects search ball.
 */

#define Intersect(node, r2ball, pos, done)				\
{									\
    real _dxl, _dxr, _dyl, _dyr, _dzl, _dzr, _dr2;			\
    _dxl = node.bnd.minb[0] - pos[0];					\
    _dxr = pos[0] - node.bnd.maxb[0];					\
    if (_dxl > 0.0) {							\
        _dr2 = _dxl*_dxl;						\
        if (_dr2 > r2ball) goto done;					\
    } else if (_dxr > 0.0) {						\
        _dr2 = _dxr*_dxr;						\
        if (_dr2 > r2ball) goto done;					\
    } else								\
        _dr2 = 0.0;							\
    _dyl = node.bnd.minb[1] - pos[1];					\
    _dyr = pos[1] - node.bnd.maxb[1];					\
    if (_dyl > 0.0) {							\
        _dr2 += _dyl*_dyl;						\
        if (_dr2 > r2ball) goto done;					\
    } else if (_dyr > 0.0) {						\
        _dr2 += _dyr*_dyr;						\
        if (_dr2 > r2ball) goto done;					\
    }									\
    _dzl = node.bnd.minb[2] - pos[2];					\
    _dzr = pos[2] - node.bnd.maxb[2];					\
    if (_dzl > 0.0) {							\
        _dr2 += _dzl*_dzl;						\
        if (_dr2 > r2ball) goto done;					\
    } else if (_dzr > 0.0) {						\
        _dr2 += _dzr*_dzr;						\
        if (_dr2 > r2ball) goto done;					\
    }									\
}

/* 
 * Prototypes for local routines.
 */

local void ball_search(smxptr, real, real *);
local int  ball_gather(smxptr, real, real *);

/*
 * INIT_SMOOTH: create smoothing context and initialize arrays.
 */

smxptr init_smooth(kdxptr kd, int nsmooth, real slope)
{
    smxptr sm;
    int i, listsize = nsmooth + EXTLIST;

    sm = (smxptr) allocate(sizeof(smcontext));
    sm->kd = kd;
    sm->nsmooth = nsmooth;
    sm->coefs[0] =  1.0 - slope * (14.0/45.0);
    sm->coefs[1] = slope;
    sm->coefs[2] = -1.5 - slope * (31.0/30.0);
    sm->coefs[3] = 0.75 + slope * (7.0/20.0);
    sm->coefs[4] = 0.25 + slope * (1.0/180.0);
    for (i = 0; i < NRPROF; i++)
        sm->rprof[i] = 0;
    sm->pque = (pqnode *) allocate(nsmooth * sizeof(pqnode));
    sm->r2ball = (real *) allocate(kd->ngas * sizeof(real));
    sm->r2list = (real *) allocate(listsize * sizeof(real));
    sm->inlist = (int *) allocate(listsize * sizeof(int));
    PQInit(sm->pque, nsmooth);
    sm->cpustart = cputime();
    return (sm);
}

/*
 * FINISH_SMOOTH: clean up after smoothing operations.
 */

void finish_smooth(smxptr sm)
{
    free(sm->r2ball);
    free(sm->pque);
    free(sm->r2list);
    free(sm->inlist);
    free(sm);
}

/*
 * SMOOTH: construct neighbor lists and invoke smoothing function.
 */

void smooth(smxptr sm, void (*smproc)(smxptr, int, int))
{
    kdnode *ntab = sm->kd->ntab;
    bodyptr *bptr = sm->kd->bptr;
    int nsmooth = sm->nsmooth, pi, pj, pnext, pscan, nmisc, nball, cell, j, k;
    pqnode *pque = sm->pque, *pqend = sm->pque + nsmooth, *pq;
    real h2, xb2;

    for (pi = 0; pi < sm->kd->ngas; ++pi) {     /* loop over gas bodies     */
        ClrFlag(bptr[pi], INQUE);               /* reset in-que flag        */
        ClrFlag(bptr[pi], DONE);		/* and done flag            */
    }
    pj = 0;                                     /* start at zeroth body     */
    for (pq = pque; pq < pqend; ++pq) {         /* and preload priority que */
        pq->pind = pj;                          /* place each body in que   */
        SetFlag(bptr[pj++], INQUE);             /* and flag it as such      */
    }
    pnext = 0;                                  /* set next body to process */
    pscan = 0;                                  /* set high-water mark      */
    nmisc = 0;					/* count miscounts of que   */
    for ( ; ; ) {                               /* loop processing bodies   */
        if (! Done(bptr[pnext]))                /* is pnext not yet done?   */
            pi = pnext;                         /* do it using current que  */
        else {                                  /* else restart at next one */
            while (pscan < sm->kd->ngas && Done(bptr[pscan]))
                ++pscan;                        /* scan for one to do       */
            if (pscan == sm->kd->ngas)          /* exit loop if all done    */
                break;
            pi = pscan;                         /* else do this one         */
            for (pq = pque; pq < pqend; ++pq)   /* loop over priority que   */
		ClrFlag(bptr[pq->pind], INQUE);	/* remove all old bodies    */
            cell = KDROOT;                      /* start at root of tree    */
            while (cell < sm->kd->nsplit) {     /* find bucket holding pi   */
                if (Pos(bptr[pi])[ntab[cell].dim] < ntab[cell].split)
                    cell = Lower(cell);         /* follow lower branch      */
                else
                    cell = Upper(cell);         /* follow upper branch      */
            }
            pj = MIN(ntab[cell].first, sm->kd->ngas - nsmooth);
						/* start with local bucket  */
            for (pq = pque; pq < pqend; ++pq) { /* and reload priority que  */
                pq->pind = pj;                  /* place local body in que  */
		SetFlag(bptr[pj++], INQUE);     /* and flag it as such      */
            }
        }
        for (pq = pque; pq < pqend; ++pq) {     /* loop over current que    */
	    DISTSQV(pq->pqkey, Pos(bptr[pi]), Pos(bptr[pq->pind]));
						/* store dist to each body  */
        }
        PQBuild(pque, nsmooth, sm->pqhead);     /* build que of candidates  */
        ball_search(sm, sm->pqhead->pqkey, Pos(bptr[pi]));
						/* and search for neighbors */
        sm->r2ball[pi] = sm->pqhead->pqkey;     /* record dist to outer one */
	SetFlag(bptr[pi], DONE);		/* flag this body as done   */
        nball = 0;                              /* init count of neighbors  */
        pnext = pi;                             /* look for next body to do */
        h2 = sm->pqhead->pqkey;                 /* within radius of ball    */
        for (pq = pque; pq < pqend; ++pq)       /* loop over bodies in que  */
            if (pq != sm->pqhead) {             /* if not most distant one  */
		pj = pq->pind;
                sm->r2list[nball] = pq->pqkey;  /* store dist^2 in r2list   */
                sm->inlist[nball++] = pj;	/* store index in inlist    */
                if (! Done(bptr[pj]) && pq->pqkey < h2) {
						/* not yet done, & closer?  */
                    pnext = pj;                 /* then make it the next    */
                    h2 = pq->pqkey;             /* and keep track of dist   */
                }
            }
        if (nball != nsmooth - 1)
	    nmisc++;
	Smooth(bptr[pi]) = 0.5 * rsqrt(sm->r2ball[pi]);
        (*smproc)(sm, pi, nball);
        for (j = 0; j < nball; j++)
            if (sm->inlist[j] != pi) {
                xb2 = sm->r2list[j] / sm->r2ball[pi];
                k = 0;
                while (k < NRPROF && xb2 <= 1.0) {
                    sm->rprof[k]++;
                    k++;
                    xb2 = 2 * xb2;
                }
            }
    }
#if defined(NOTEMISCOUNT)
    if (nmisc != 0)
        eprintf("[%s.smooth: %d miscount%s]\n", getargv0(), nmisc,
		nmisc > 1 ? "s" : "");
#endif
}

/*
 * RESMOOTH: reconstruct neighbor lists and invoke smoothing function.
 */

void resmooth(smxptr sm, void (*smproc)(smxptr, int, int))
{
    bodyptr *bptr = sm->kd->bptr;
    int nmisc, pi, nball, i;

    nmisc = 0;
    for (pi = 0; pi < sm->kd->ngas; ++pi) {
	nball = ball_gather(sm, sm->r2ball[pi], Pos(bptr[pi]));
	if (nball != sm->nsmooth - 1)
	    nmisc++;
	(*smproc)(sm, pi, nball);
    }
#if defined(NOTEMISCOUNT)
    if (nmisc != 0)
        eprintf("[%s.resmooth: %d miscount%s]\n", getargv0(), nmisc,
		nmisc > 1 ? "s" : "");
#endif
}

/*
 * REPORT_SMOOTH: output report on smoothing operations.
 */

void report_smooth(smxptr sm, stream ostr, string options)
{
    bodyptr *bptr = sm->kd->bptr;
    real rbmin, rbsum, xi[NRPROF-1];
    int *rprof = sm->rprof, nlev[MAXLEV+1], nupd, i;

    rbmin = Smooth(bptr[0]);
    rbsum = 0.0;
    nupd = 0;
    for (i = 0; i <= MAXLEV; i++)
	nlev[i] = 0;
    for (i = 0; i < sm->kd->ngas; i++) {
        if (Update(bptr[i])) {
	    nlev[NewLevel(bptr[i])]++;
            rbmin = MIN(rbmin, Smooth(bptr[i]));
            rbsum += Smooth(bptr[i]);
            nupd++;
        }
    }
    for (i = 0; i <= NRPROF - 2; i++)
	xi[i] = -1 + rpow(8.0, i/2.0) * (rprof[i] - rprof[i+1]) /
                       ((1 - rsqrt(0.125)) * rprof[0]);
    fprintf(ostr, "\n%10s %5s %5s %5s %11s %9s %9s %9s %7s\n",
            "xi8", "xi6", "xi4", "xi2", "hmin", "havg",
	    "freqmax", "freqavg", "CPUsph"); 
    fprintf(ostr, "%10.2f %5.2f %5.2f %5.2f %11.4f %9.4f %9.2f %9.2f %7.3f\n",
	    xi[8], xi[6], xi[4], xi[2], rbmin, rbsum / nupd,
	    sm->freqmax, sm->freqsum / nupd, cputime() - sm->cpustart);
    if (scanopt(options, "corrfunc")) {
	fprintf(ostr, "\n     ");
	for (i = NRPROF - 2; i >= 0; i--)
	    fprintf(ostr, "   xi%d", i);
	fprintf(ostr, "\n     ");
	for (i = NRPROF - 2; i >= 0; i--)
	    fprintf(ostr, "%6.2f", xi[i]);
	fprintf(ostr, "\n");
    }	
    if (scanopt(options, "levelhist")) {
	fprintf(ostr, "\n     ");
	for (i = 0; i <= MAXLEV; i++)
	    if (nlev[i] != 0)
		fprintf(ostr, i<10 ? "  lev%d" : " lev%d", i);
	fprintf(ostr, "\n     ");
	for (i = 0; i <= MAXLEV; i++)
	    if (nlev[i] != 0)
		fprintf(ostr, "%6d", nlev[i]);
	fprintf(ostr, "\n");
    }	
    fflush(ostr);
}

/*
 * BALL_SEARCH: given approximate priority que, construct final version.
 */

local void ball_search(smxptr sm, real r2ball, real *ri)
{
    kdnode *ntab = sm->kd->ntab;
    bodyptr *bptr = sm->kd->bptr;
    pqnode *pq = sm->pqhead;
    int cell, cp, ct, pj;
    real dist2;

    cell = KDROOT;                              /* start at root of tree    */
    while (cell < sm->kd->nsplit) {             /* descend to local bucket  */
        if (ri[ntab[cell].dim] < ntab[cell].split)
            cell = Lower(cell);
        else
            cell = Upper(cell);
    }
    for (pj = ntab[cell].first; pj <= ntab[cell].last; ++pj)
	if (! InQue(bptr[pj])) {		/* in bucket, but not que?  */
	    DISTSQV(dist2, ri, Pos(bptr[pj]));  /* compute dist^2 to center */
	    if (dist2 < r2ball) {		/* within current ball?     */
		ClrFlag(bptr[pq->pind], INQUE);	/* drop furthest from que   */
		SetFlag(bptr[pj], INQUE);	/* and add this one to que  */
		pq->pqkey = dist2;              /* store its distance       */
		pq->pind = pj;                  /* and its index            */
		PQReplace(pq);                  /* move to rightful place   */
		r2ball = pq->pqkey;             /* adopt new search radius  */
	    }
	}
    while (cell != KDROOT) {			/* scan back toward root    */
        cp = Sibling(cell);
        ct = cp;
        SetNext(ct);
	do {
            Intersect(ntab[cp], r2ball, ri, GetNextCell);
						/* got intersection to test */
            if (cp < sm->kd->nsplit) {          /* not yet down to bucket?  */
                cp = Lower(cp);
                continue;
            } else                              /* scan bucket for winners  */
                for (pj = ntab[cp].first; pj <= ntab[cp].last; ++pj)
		    if (! InQue(bptr[pj])) {	/* not already in the que?  */
			DISTSQV(dist2, ri, Pos(bptr[pj]));
			if (dist2 < r2ball) {	/* but within current ball? */
			    ClrFlag(bptr[pq->pind], INQUE);
			    SetFlag(bptr[pj], INQUE);
			    pq->pqkey = dist2;
			    pq->pind = pj;
			    PQReplace(pq);
			    r2ball = pq->pqkey;
			}
		    }
          GetNextCell:
            SetNext(cp);
        } while (cp != ct);
        cell = Parent(cell);			/* climb down towards root  */
    }
    sm->pqhead = pq;
}

/*
 * BALL_GATHER: search tree for all bodies strictly within search radius.
 */

local int ball_gather(smxptr sm, real r2ball, real *ri)
{
    kdnode *ntab = sm->kd->ntab;
    bodyptr *bptr = sm->kd->bptr;
    int cp, nball, pj;
    real dist2;

    nball = 0;
    cp = KDROOT;
    do {
        Intersect(ntab[cp], r2ball, ri, GetNextCell);
						/* got intersection to test */
        if (cp < sm->kd->nsplit) {
            cp = Lower(cp);
            continue;
        } else {
            for (pj = ntab[cp].first; pj <= ntab[cp].last; ++pj) {
		DISTSQV(dist2, ri, Pos(bptr[pj]));
                if (dist2 < r2ball) {
                    sm->r2list[nball] = dist2;
                    sm->inlist[nball++] = pj;
                }
            }
        }
      GetNextCell:
        SetNext(cp);
    } while (cp != KDROOT);
    if (nball > sm->nsmooth + EXTLIST)
	error("%s: gathered list overflow\n", getargv0());
    return (nball);
}
