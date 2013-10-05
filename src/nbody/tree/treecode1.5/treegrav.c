/****************************************************************************/
/* TREEGRAV.C: routines to compute gravity. Public routines: gravcalc().    */
/* Copyright (c) 2012 by Joshua E. Barnes, Honolulu, HI.                    */
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
 
  if (active == NULL) {				/* if this is the 1st call  */
    actmax = FACTIVE * 216 * (tdepth + 1);	/* estimate list length     */
#ifndef QUICKSCAN
    if (theta > 0.1)
      actmax = actmax * rpow(theta,-2.5);	/* allow for opening angle  */
    else
      actmax = 5.0 * ncell;			/* guess total node count   */
#endif
    active = (nodeptr *) allocate(actmax * sizeof(nodeptr));
    interact = (cellptr) allocate(actmax * sizeof(cell));
  }
  cpustart = cputime();				/* record time, less alloc  */
  acttot = nfcalc = nbbcalc = nbccalc = 0;	/* zero cumulative counters */
  active[0] = (nodeptr) root;			/* initialize active list   */
  CLRV(rmid);					/* set center of root cell  */
  walktree(active, active + 1, interact, interact + actmax,
	   (nodeptr) root, rsize, rmid);	/* scan tree, update forces */
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
#ifdef TESTFORCE
	    More(cptr) = More(*ap);
	    Next(cptr) = Next(*ap);
#endif
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
 
  p15 = ((real) 1.5) * psize;			/* premultiply cell size    */
  dk = Pos(c)[0] - pmid[0];			/* find distance to midpnt  */
  if (ABS(dk) > p15)				/* if c does not touch p    */
    return (TRUE);				/* then accept interaction  */
  dk = Pos(c)[1] - pmid[1];			/* find distance to midpnt  */
  if (ABS(dk) > p15)				/* if c does not touch p    */
    return (TRUE);				/* then accept interaction  */
  dk = Pos(c)[2] - pmid[2];			/* find distance to midpnt  */
  if (ABS(dk) > p15)				/* if c does not touch p    */
    return (TRUE);				/* then accept interaction  */
  return (FALSE);				/* else do not accept it    */
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
 
local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid)
{
  nodeptr q;
  int k;
  vector qmid;
 
  if (Type(p) == CELL) {                        /* fanout over descendents  */
    for (q = More(p); q != Next(p); q = Next(q)) {
                                                /* loop over all subcells   */
      for (k = 0; k < NDIM; k++)                /* set subcell's midpoint   */
	qmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - psize : psize) / 4;
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
 
  SETV(pos0, Pos(p0));                          /* copy position of body    */
  phi0 = 0.0;                                   /* init total potential     */
  CLRV(acc0);                                   /* and total acceleration   */
  if (usequad)                                  /* if using quad moments    */
    sumcell(interact, cptr, pos0, &phi0, acc0); /* sum cell forces w quads  */
  else                                          /* not using quad moments   */
    sumnode(interact, cptr, pos0, &phi0, acc0); /* sum cell forces wo quads */
  sumnode(bptr, interact + actmax, pos0, &phi0, acc0);
                                                /* sum forces from bodies   */
  Phi(p0) = phi0;                               /* store total potential    */
  SETV(Acc(p0), acc0);                          /* and total acceleration   */
  nfcalc++;                                     /* update counters          */
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
  real eps2 = eps*eps, dr2, drab, phi_p, mr3i;
  vector dr;
 
  for (p = start; p < finish; p++) {            /* loop over node list      */
    DOTPSUBV(dr2, dr, Pos(p), pos0);            /* compute separation       */
                                                /* and distance squared     */
    dr2 += eps2;                                /* add standard softening   */
    drab = rsqrt(dr2);                          /* form scalar "distance"   */
    phi_p = Mass(p) / drab;                     /* get partial potential    */
    *phi0 -= phi_p;                             /* decrement tot potential  */
    mr3i = phi_p / dr2;                         /* form scale factor for dr */
    ADDMULVS(acc0, dr, mr3i);                   /* sum partial acceleration */}
}

#ifndef TESTFORCE

/* <A NAME="sumcell"></A>
 * SUMCELL: add up body-cell interactions.
 */
 
local void sumcell(cellptr start, cellptr finish, vector pos0,
		   real *phi0, vector acc0)
{
  cellptr p;
  real eps2 = eps*eps, dr2, drab, phi_p, mr3i;
  real drqdr, dr5i, phi_q;
  vector dr, qdr;
 
  for (p = start; p < finish; p++) {            /* loop over node list      */
    DOTPSUBV(dr2, dr, Pos(p), pos0);            /* do mono part of force    */
    dr2 += eps2;
    drab = rsqrt(dr2);
    phi_p = Mass(p) / drab;
    mr3i = phi_p / dr2;
    DOTPMULMV(drqdr, qdr, Quad(p), dr);         /* do quad part of force    */
    dr5i = ((real) 1.0) / (dr2 * dr2 * drab);
    phi_q = ((real) 0.5) * dr5i * drqdr;
    *phi0 -= phi_p + phi_q;                     /* add mono and quad pot    */
    mr3i += ((real) 5.0) * phi_q / dr2;
    ADDMULVS2(acc0, dr, mr3i, qdr, -dr5i);      /* add mono and quad acc    */
  }
}

#endif

#ifdef TESTFORCE

local void exactgrav(cellptr p, vector pos0, real *phip, vector accp);
local real smoothmass(cellptr p, real rad);

/* 
 * SUMCELL: add up body-cell interactions.
 */
 
local void sumcell(cellptr start, cellptr finish, vector pos0,
		   real *phi0, vector acc0)
{
  cellptr p;
  real eps2 = eps*eps, dr2, drab, phi_p, mr3i;
  real drqdr, dr5i, phi_q, phi_ap, phi_ex;
  vector dr, qdr, acc_ap, acc_ex;
 
  for (p = start; p < finish; p++) {            /* loop over node list      */
    DOTPSUBV(dr2, dr, Pos(p), pos0);            /* do mono part of force    */
    dr2 += eps2;
    drab = rsqrt(dr2);
    phi_p = Mass(p) / drab;
    mr3i = phi_p / dr2;
    DOTPMULMV(drqdr, qdr, Quad(p), dr);         /* do quad part of force    */
    dr5i = ((real) 1.0) / (dr2 * dr2 * drab);
    phi_q = ((real) 0.5) * dr5i * drqdr;
    *phi0 -= phi_p + phi_q;                     /* add mono and quad pot    */
    mr3i += ((real) 5.0) * phi_q / dr2;
    ADDMULVS2(acc0, dr, mr3i, qdr, -dr5i);      /* add mono and quad acc    */

    phi_ap = -(phi_p + phi_q);
    CLRV(acc_ap);
    ADDMULVS2(acc_ap, dr, mr3i, qdr, -dr5i);
    exactgrav(p, pos0, &phi_ex, acc_ex);
    printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	   pos0[0], Mass(p), distv(pos0, Pos(p)),
	   Pos(p)[0], Pos(p)[1], Pos(p)[2],
	   (phi_ap - phi_ex) / phi_ex, distv(acc_ap, acc_ex) / absv(acc_ex));

  }
}

local void exactgrav(cellptr p0, vector pos0, real *phi0, vector acc0)
{
  nodeptr p;
  real eps2 = eps * eps, dr2, drab, phi_p, mr3i;
  vector dr;

  *phi0 = 0.0;
  CLRV(acc0);
  p = More(p0);
  while (p != Next(p0))
    if (Type(p) == CELL)
      p = More(p);
    else {
      DOTPSUBV(dr2, dr, Pos(p), pos0);
      dr2 += eps2;
      drab = rsqrt(dr2);
      phi_p = Mass(p) / drab;
      *phi0 -= phi_p;
      mr3i = phi_p / dr2;
      ADDMULVS(acc0, dr, mr3i);
      p = Next(p);
    }
}

local real smoothmass(cellptr p, real rad)
{
  nodeptr q;
  real mtot = 0, eps2 = eps * eps;
  double ri;

  q = More(p);
  while (q != Next(p))
    if (Type(q) == CELL)
      q = More(q);
    else {
      ri = distv(Pos(p), Pos(q));
      if (ri > 0.0) 
	mtot += Mass(q) *
	  ((rad + ri + eps2/ri) / sqrt((rad + ri)*(rad + ri) + eps2) +
	   (rad - ri - eps2/ri) / sqrt((rad - ri)*(rad - ri) + eps2)) / 2;
      else
	mtot += Mass(q) / rpow(1 + eps2/rad*rad, 1.5);
      q = Next(q);
    }
  return (mtot);
}

#endif
