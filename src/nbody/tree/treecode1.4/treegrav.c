/****************************************************************************/
/* TREEGRAV.C: routines to compute gravity. Public routines: gravcalc().    */
/* Copyright (c) 2000 by Joshua E. Barnes, Honolulu, HI.                    */
/****************************************************************************/
 
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "treedefs.h"
 
/* Local routines to perform force calculations. */
 
local void walktree(nodeptr *, nodeptr *, cellptr, cellptr,
		    nodeptr, real, vector);
local bool accept(nodeptr, real, vector);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector);
local void gravsum(bodyptr, cellptr, cellptr);
local void sumnode(cellptr, cellptr, vector, real *, vector);
local void sumcell(cellptr, cellptr, vector, real *, vector);
 
/* Lists of active nodes and interactions. */
 
#ifndef FACTIVE
#  define FACTIVE  1.0				/* active list fudge factor */
#endif
 
local int actmax;                               /* length as allocated      */
local int acttot;                               /* actual active length     */

local nodeptr *active = NULL;                   /* list of nodes tested     */

local cellptr interact = NULL;                  /* list of interactions     */

/* <A NAME="gravcalc"></A>
 * GRAVCALC: perform force calculation on all particles.
 */
 
void gravcalc(void)
{
    double cpustart;
    vector rmid;
 
    if (active == NULL) {                       /* if this is the 1st call  */
        actmax = FACTIVE * 216 * (tdepth + 1);  /* estimate list length     */
#ifndef QUICKSCAN
	if (theta > 0.1)
	    actmax = actmax * rpow(theta,-2.5);	/* allow for opening angle  */
	else
	    actmax = 5.0 * ncell;		/* guess total node count   */
#endif
        active = (nodeptr *) allocate(actmax * sizeof(nodeptr));
        interact = (cellptr) allocate(actmax * sizeof(cell));
    }
    cpustart = cputime();			/* record time, less alloc  */
    acttot = nfcalc = nbbcalc = nbccalc = 0;    /* zero cumulative counters */
    active[0] = (nodeptr) root;                 /* initialize active list   */
    CLRV(rmid);                                 /* set center of root cell  */
    walktree(active, active + 1, interact, interact + actmax,
             (nodeptr) root, rsize, rmid);      /* scan tree, update forces */
    cpuforce = cputime() - cpustart;		/* store CPU time w/o alloc */
}

/* <A NAME="walktree"></A>
 * WALKTREE: do a complete walk of the tree, building the interaction
 * list level-by-level and computing the resulting force on each body.
 */
 
local void walktree(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                    nodeptr p, real psize, vector pmid)
{
  nodeptr *np, *ap, q;
  int actsafe;
 
  if (Update(p)) {				/* are new forces needed?   */
    np = nptr;					/* start new active list    */
    actsafe = actmax - NSUB;			/* leave room for NSUB more */
    for (ap = aptr; ap < nptr; ap++) {		/* loop over active nodes   */
      if (Type(*ap) == CELL) {			/* is this node a cell?     */
	if (accept(*ap, psize, pmid)) {		/* does it pass the test?   */
	  if (Mass(*ap) > 0.0) {		/* and contribute to field? */
	    Mass(cptr) = Mass(*ap);		/* copy to interaction list */
	    SETV(Pos(cptr), Pos(*ap));
	    SETM(Quad(cptr), Quad(*ap));
	    cptr++;				/* and bump cell array ptr  */
	  }
	} else {				/* else it fails the test   */
	  if (np - active >= actsafe)		/* check list has room      */
	    fatal("walktree: active list overflow\n");
	  for (q = More(*ap); q != Next(*ap); q = Next(q))
						/* loop over all subcells   */
	    *np++= q;				/* put on new active list   */
	}
      } else					/* else this node is a body */
	if (*ap != p && Mass(*ap) > 0.0) {	/* if not self-interaction  */
	  --bptr;				/* bump body array ptr      */
	  Mass(bptr) = Mass(*ap);		/* and copy data to array   */
	  SETV(Pos(bptr), Pos(*ap));
	}
    }
    acttot = MAX(acttot, np - active);		/* keep track of max active */
    if (np != nptr) {				/* if new actives listed    */
      walksub(nptr, np, cptr, bptr, p, psize, pmid);
						/* then visit next level    */
    } else {					/* else no actives left, so */
      if (Type(p) != BODY)			/* must have found a body   */
	fatal("walktree: recursion terminated with cell\n"
	      "  p = 0x%x  psize   = %.8f  Mass(p) = %g\n"
	      "  pmid =   (%.8f,%.8f,%.8f)\n"
	      "  Pos(p) = (%.8f,%.8f,%.8f)\n",
	      (int) p, psize, Mass(p),
	      pmid[0], pmid[1], pmid[2],
	      Pos(p)[0], Pos(p)[1], Pos(p)[2]);
      gravsum((bodyptr) p, cptr, bptr);		/* sum force on the body    */
    }
  }
}

#ifdef QUICKSCAN
 
/* <A NAME="accept"></A>
 * ACCEPT: quick criterion accepts any cell not touching cell p.
 */
 
local bool accept(nodeptr c, real psize, vector pmid)
{
    real p15, dk;
 
    p15 = ((real) 1.5) * psize;                 /* premultiply cell size    */
    dk = Pos(c)[0] - pmid[0];                   /* find distance to midpnt  */
    if (ABS(dk) > p15)                          /* if c does not touch p    */
        return (TRUE);                          /* then accept interaction  */
    dk = Pos(c)[1] - pmid[1];                   /* find distance to midpnt  */
    if (ABS(dk) > p15)                          /* if c does not touch p    */
        return (TRUE);                          /* then accept interaction  */
    dk = Pos(c)[2] - pmid[2];                   /* find distance to midpnt  */
    if (ABS(dk) > p15)                          /* if c does not touch p    */
        return (TRUE);                          /* then accept interaction  */
    return (FALSE);                             /* else do not accept it    */
}
 
#else
 
/*
 * ACCEPT: standard criterion accepts cell if its critical radius
 * does not intersect cell p, and also imposes above condition.
 */
 
local bool accept(nodeptr c, real psize, vector pmid)
{
    real dmax, dsq, dk;
    int k;
 
    dmax = psize;                               /* init maximum distance    */
    dsq = 0.0;                                  /* and squared min distance */
    for (k = 0; k < NDIM; k++) {                /* loop over space dims     */
        dk = Pos(c)[k] - pmid[k];               /* form distance to midpnt  */
        if (dk < 0)                             /* and get absolute value   */
            dk = - dk;
        if (dk > dmax)                          /* keep track of max value  */
            dmax = dk;
        dk -= ((real) 0.5) * psize;             /* allow for size of cell   */
        if (dk > 0)
            dsq += dk * dk;                     /* sum min dist to cell ^2  */
    }
    return (dsq > Rcrit2(c) &&                  /* test angular criterion   */
              dmax > ((real) 1.5) * psize);     /* and adjacency criterion  */
}
 
#endif

/* <A NAME="walksub"></A>
 * WALKSUB: test next level's active list against subnodes of p.
 */
 
/* #define PTRACE  0x361c70 */
/* #define PTRACE  0x3eff910 */
/* #define PTRACE  0x3eff960 */

local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid)
{
  nodeptr q;
  int k;
  vector qmid;
 
#if PTRACE
  if (PTRACE == (int) p)
    eprintf("\nwalksub(0x%x, 0x%x, 0x%x, 0x%x, p = 0x%x,\n"
	    "        psize = %.8f, pmid = [%.8f,%.8f,%.8f])\n",
	    (int) nptr, (int) np, (int) cptr, (int) bptr,
	    (int) p, psize, pmid[0], pmid[1], pmid[2]);
#endif
  if (Type(p) == CELL) {                        /* fanout over descendents  */
#if PTRACE
    if (PTRACE == (int) p)
      eprintf("  p is a cell; Pos: [%.8f,%.8f,%.8f]  Mass: %.8f\n",
	      Pos(p)[0], Pos(p)[1], Pos(p)[2], Mass(p));
#endif
    for (q = More(p); q != Next(p); q = Next(q)) {
                                                /* loop over all subcells   */
#if PTRACE
      if (PTRACE == (int) q)
	eprintf("parent of 0x%x is 0x%x\n", (int) q, (int) p);
      if (PTRACE == (int) p)
	eprintf("    %s 0x%x: Pos: [%.8f,%.8f,%.8f]  Mass: %.8f\n"
		"      offset/psize = [%.8f,%.8f,%.8f]\n",
		Type(q) == CELL ? "cell" : "body", (int) q,
		Pos(q)[0], Pos(q)[1], Pos(q)[2], Mass(q), 
		(Pos(q)[0] - pmid[0])/psize,
		(Pos(q)[1] - pmid[1])/psize,
		(Pos(q)[2] - pmid[2])/psize);
#endif      
      for (k = 0; k < NDIM; k++)                /* set subcell's midpoint   */
	qmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - psize : psize) / 4;
#if PTRACE
      if (PTRACE == (int) p)
	eprintf("    qmid = [%.8f,%.8f,%.8f]\n", qmid[0], qmid[1], qmid[2]);
#endif
      walktree(nptr, np, cptr, bptr, q, psize / 2, qmid);
                                                /* recurse on subcell       */
    }
  } else {                                      /* extend "virtual" tree    */
    for (k = 0; k < NDIM; k++)                  /* set virtual midpoint     */
      qmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - psize : psize) / 4;
    walktree(nptr, np, cptr, bptr, p, psize / 2, qmid);
                                                /* and search next level    */
  }
}

/* <A NAME="gravsum"></A>
 * GRAVSUM: compute gravitational field at body p0.
 */
 
local void gravsum(bodyptr p0, cellptr cptr, cellptr bptr)
{
    vector pos0, acc0;
    real phi0;
 
    SETV(pos0, Pos(p0));                        /* copy position of body    */
    phi0 = 0.0;                                 /* init total potential     */
    CLRV(acc0);                                 /* and total acceleration   */
    if (usequad)                                /* if using quad moments    */
        sumcell(interact, cptr, pos0, &phi0, acc0);
                                                /* sum cell forces w quads  */
    else                                        /* not using quad moments   */
        sumnode(interact, cptr, pos0, &phi0, acc0);
                                                /* sum cell forces wo quads */
    sumnode(bptr, interact + actmax, pos0, &phi0, acc0);
                                                /* sum forces from bodies   */
    Phi(p0) = phi0;                             /* store total potential    */
    SETV(Acc(p0), acc0);                        /* and total acceleration   */
    nfcalc++;                                   /* update counters          */
    nbbcalc += interact + actmax - bptr;
    nbccalc += cptr - interact;
}

/* <A NAME="sumnode"></A>
 * SUMNODE: add up body-node interactions.
 */
 
local void sumnode(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
#ifndef IBMAIX
    real eps2, dr2, drab, phi_p, mr3i;
#else
    double eps2, dr2, drab, phi_p, mr3i;
#endif
    vector dr;
 
    eps2 = eps * eps;                           /* avoid extra multiplys    */
    for (p = start; p < finish; p++) {          /* loop over node list      */
        DOTPSUBV(dr2, dr, Pos(p), pos0);        /* compute separation       */
                                                /* and distance squared     */
        dr2 += eps2;                            /* add standard softening   */
        drab = rsqrt(dr2);                      /* form scalar "distance"   */
        phi_p = Mass(p) / drab;                 /* get partial potential    */
        *phi0 -= phi_p;                         /* decrement tot potential  */
        mr3i = phi_p / dr2;                     /* form scale factor for dr */
        ADDMULVS(acc0, dr, mr3i);               /* sum partial acceleration */}
    
}

/* <A NAME="sumcell"></A>
 * SUMCELL: add up body-cell interactions.
 */
 
local void sumcell(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i, drqdr, dr5i, phi_q;
    vector dr, qdr;
 
    eps2 = eps * eps;
    for (p = start; p < finish; p++) {          /* loop over node list      */
        DOTPSUBV(dr2, dr, Pos(p), pos0);        /* do mono part of force    */
        dr2 += eps2;
        drab = rsqrt(dr2);
        phi_p = Mass(p) / drab;
        mr3i = phi_p / dr2;
        DOTPMULMV(drqdr, qdr, Quad(p), dr);     /* do quad part of force    */
        dr5i = ((real) 1.0) / (dr2 * dr2 * drab);
        phi_q = ((real) 0.5) * dr5i * drqdr;
        *phi0 -= phi_p + phi_q;                 /* add mono and quad pot    */
        mr3i += ((real) 5.0) * phi_q / dr2;
        ADDMULVS2(acc0, dr, mr3i, qdr, -dr5i);  /* add mono and quad acc    */
    }
}
