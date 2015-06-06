/*
 * treegrav.c: routines to compute gravity. Public routines: gravcalc().
 * Copyright (c) 2015 by Joshua E. Barnes, Kyoto, Japan.
 */
 
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "treedefs.h"
 
//  Local routines to perform force calculations.
 
local void walktree(nodeptr *, nodeptr *, cellptr, cellptr,
		    nodeptr, real, vector);
local bool accept(nodeptr, real, vector);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector);
local void gravsum(bodyptr, cellptr, cellptr);
local void sumnode(cellptr, cellptr, vector, real *, vector);
local void sumcell(cellptr, cellptr, vector, real *, vector);
 
//  Lists of active nodes and interactions.
 
#if !defined(FACTIVE)
#  define FACTIVE  1.0				// active list fudge factor
#endif
 
local int actmax;                               // length as allocated
local int acttot;                               // actual active length

local nodeptr *active = NULL;                   // list of nodes tested
local cellptr interact = NULL;                  // list of interactions

//  gravcalc: perform force calculation on all particles.
//  _____________________________________________________
 
void gravcalc(void)
{
  double cpustart;
  vector rmid;
 
  if (active == NULL) {				// if this is the 1st call
    actmax = FACTIVE * 216 * (tdepth + 1);	// estimate list length
#if !defined(QUICKSCAN)
    if (theta > 0.1)
      actmax = actmax * rpow(theta,-2.5);	// allow for opening angle
    else
      actmax = 5.0 * ncell;			// guess total node count
#endif
    active = (nodeptr *) allocate(actmax * sizeof(nodeptr));
    interact = (cellptr) allocate(actmax * sizeof(cell));
  }
  cpustart = cputime();				// record time, less alloc
  acttot = nfcalc = nbbcalc = nbccalc = 0;	// zero cumulative counters
  active[0] = (nodeptr) root;			// initialize active list
  CLRV(rmid);					// set center of root cell
  walktree(active, active + 1, interact, interact + actmax,
	   (nodeptr) root, rsize, rmid);	// scan tree, update forces
  cpuforce = cputime() - cpustart;		// store CPU time w/o alloc
}

//  walktree: do a complete walk of the tree, building the interaction
//  list level-by-level and computing the resulting force on each body.
//  ___________________________________________________________________
 
local void walktree(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                    nodeptr p, real psize, vector pmid)
{
  nodeptr *np, *ap, q;
  int actsafe;
  matrix trQM;
 
  if (Update(p)) {				// new forces needed in node?
    np = nptr;					// start new active list
    actsafe = actmax - NSUB;			// leave room for NSUB more
    for (ap = aptr; ap < nptr; ap++) {		// loop over active nodes
      if (Type(*ap) == CELL) {			// is this node a cell?
	if (accept(*ap, psize, pmid)) {		// does it pass the test?
	  if (Mass(*ap) > 0.0) {		// and contribute to field?
	    Mass(cptr) = Mass(*ap);		// copy to interaction list
	    SETV(Pos(cptr), Pos(*ap));
#if defined(SOFTCORR)
	    TRACEM(Trace(cptr), Quad(*ap));	// save trace in copy
	    SETMI(trQM);
	    MULMS(trQM, trQM, Trace(cptr)/3);
	    SUBM(Quad(cptr), Quad(*ap), trQM);	// store traceless moment
#else
	    SETM(Quad(cptr), Quad(*ap));	// copy traceless moment
#endif
	    cptr++;				// and bump cell array ptr
	  }
	} else {				// this cell fails the test
	  if (np - active >= actsafe)		// make sure list has room
	    fatal("%s.walktree: active list overflow\n", getprog());
	  for (q = More(*ap); q != Next(*ap); q = Next(q))
						// loop over all subcells
	    *np++= q;				// put them on active list
	}
      } else					// else this node is a body
	if (*ap != p && Mass(*ap) > 0.0) {	// not self-interaction?
	  --bptr;				// bump body array ptr
	  Mass(bptr) = Mass(*ap);		// and copy data to array
	  SETV(Pos(bptr), Pos(*ap));
	}
    }
    acttot = MAX(acttot, np - active);		// keep track of max active
    if (np != nptr) {				// if new actives were added
      walksub(nptr, np, cptr, bptr, p, psize, pmid);
						// then visit next level
    } else {					// else no actives left
      if (Type(p) != BODY)			// make sure we got a body
	fatal("%s.walktree: recursion terminated with cell\n"
	      "  p = 0x%x  psize   = %.8f  Mass(p) = %g\n"
	      "  pmid =   (%.8f,%.8f,%.8f)\n  Pos(p) = (%.8f,%.8f,%.8f)\n",
	      getprog(), (int) p, psize, Mass(p),
	      pmid[0], pmid[1], pmid[2], Pos(p)[0], Pos(p)[1], Pos(p)[2]);
      gravsum((bodyptr) p, cptr, bptr);		// sum force on this body
    }
  }
}

#if defined(QUICKSCAN)
 
//  accept: quick criterion accepts any cell not touching cell p.
//  _____________________________________________________________
 
local bool accept(nodeptr c, real psize, vector pmid)
{
  real p15, dk;
 
  p15 = ((real) 1.5) * psize;			// premultiply cell size
  dk = Pos(c)[0] - pmid[0];			// find distance to midpnt
  if (ABS(dk) > p15)				// if c does not touch p
    return (TRUE);				// then accept interaction
  dk = Pos(c)[1] - pmid[1];			// find distance to midpnt
  if (ABS(dk) > p15)				// if c does not touch p
    return (TRUE);				// then accept interaction
  dk = Pos(c)[2] - pmid[2];			// find distance to midpnt
  if (ABS(dk) > p15)				// if c does not touch p
    return (TRUE);				// then accept interaction
  return (FALSE);				// else do not accept it
}
 
#else
 
//  accept: standard criterion accepts cell if its critical radius
//  does not intersect cell p, and also imposes above condition.
//  ______________________________________________________________
 
local bool accept(nodeptr c, real psize, vector pmid)
{
  real dmax, dsq, dk;
  int k;
 
  dmax = psize;                                 // init maximum distance
  dsq = 0.0;                                    // and squared min distance
  for (k = 0; k < NDIM; k++) {                  // loop over space dims
    dk = Pos(c)[k] - pmid[k];                   // form distance to midpnt
    if (dk < 0)                                 // and get absolute value
      dk = - dk;
    if (dk > dmax)                              // keep track of max value
      dmax = dk;
    dk -= ((real) 0.5) * psize;                 // allow for size of cell
    if (dk > 0)
      dsq += dk * dk;                           // sum (min dist to cell)^2
  }
  return (dsq > Rcrit2(c) &&                    // test angular criterion
	  dmax > ((real) 1.5) * psize);         // and adjacency criterion
}
 
#endif

//  walksub: test next level's active list against subnodes of p.
//  _____________________________________________________________
 
local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid)
{
  nodeptr q;
  int k;
  vector qmid;
 
  if (Type(p) == CELL) {                        // fanout over descendents
    for (q = More(p); q != Next(p); q = Next(q)) {
                                                // loop over all subcells
      for (k = 0; k < NDIM; k++)                // set subcell's midpoint
	qmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - psize : psize) / 4;
      walktree(nptr, np, cptr, bptr, q, psize / 2, qmid);
                                                // recurse on subcell
    }
  } else {                                      // extend "virtual" tree
    for (k = 0; k < NDIM; k++)                  // set virtual midpoint
      qmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - psize : psize) / 4;
    walktree(nptr, np, cptr, bptr, p, psize / 2, qmid);
                                                // and search next level
  }
}

//  gravsum: compute gravitational field at body p0.
//  ________________________________________________
 
local void gravsum(bodyptr p0, cellptr cptr, cellptr bptr)
{
  vector pos0, acc0;
  real phi0;
 
  SETV(pos0, Pos(p0));                          // copy position of body
  phi0 = 0.0;                                   // init total potential
  CLRV(acc0);                                   // and total acceleration
  if (usequad)                                  // if using quad moments
    sumcell(interact, cptr, pos0, &phi0, acc0); // sum cell forces w quads
  else                                          // not using quad moments
    sumnode(interact, cptr, pos0, &phi0, acc0); // sum cell forces wo quads
  sumnode(bptr, interact + actmax, pos0, &phi0, acc0);
                                                // sum forces from bodies
  Phi(p0) = phi0;                               // store total potential
  SETV(Acc(p0), acc0);                          // and total acceleration
  nfcalc++;                                     // update counters
  nbbcalc += interact + actmax - bptr;
  nbccalc += cptr - interact;
}

//  sumnode: add up body-node interactions.
//  _______________________________________
 
local void sumnode(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
  real eps2, dr2, dr2i, dr1i, mdr1i, mdr3i;
  vector dr;
 
  eps2 = eps * eps;				// premultiply softening
  for (cellptr p = start; p < finish; p++) {	// loop over node list
    DOTPSUBV(dr2, dr, Pos(p), pos0);		// compute sep. vector
						// and square of distance
    dr2i = ((real) 1.0) / (dr2 + eps2);		// perform only division
    dr1i = rsqrt(dr2i);				// set inverse soft distance
    mdr1i = Mass(p) * dr1i;			// form partial potential
    mdr3i = mdr1i * dr2i;			// form scale factor for dr
    *phi0 -= mdr1i;				// sum potential
    ADDMULVS(acc0, dr, mdr3i);			// sum acceleration
  }
}

//  sumcell: add up body-cell interactions, including quad moments.
//  _______________________________________________________________

local void sumcell(cellptr start, cellptr finish, vector pos0,
		   real *phi0, vector acc0)
{
  real eps2, eps2thrd, dr2, dr2i, dr1i, mdr1i, mdr3i, qdr2, dr5i, phi_q;
  vector dr, qdr;
 
  eps2 = eps * eps;				// premultiply softening
  eps2thrd = eps2 / 3.0;			// predivide for soft corr
  for (cellptr p = start; p < finish; p++) {	// loop over node list
    DOTPSUBV(dr2, dr, Pos(p), pos0);		// do mono part of force
    dr2i = ((real) 1.0) / (dr2 + eps2);		// perform only division
    dr1i = rsqrt(dr2i);				// set inverse soft distance
    mdr1i = Mass(p) * dr1i;			// form mono potential
    mdr3i = mdr1i * dr2i;			// get scale factor for dr
    DOTPMULMV(qdr2, qdr, Quad(p), dr);		// do quad part of force
#if defined(SOFTCORR)
    qdr2 -= eps2thrd * Trace(p);		// apply Keigo's correction
#endif
    dr5i = ((real) 3.0) * dr2i * dr2i * dr1i;	// factor 3 saves a multiply
    phi_q = ((real) 0.5) * dr5i * qdr2;		// form quad potential
    mdr3i += ((real) 5.0) * phi_q * dr2i;	// adjust radial term
    *phi0 -= mdr1i + phi_q;			// sum mono and quad pot
    ADDMULVS2(acc0, dr, mdr3i, qdr, - dr5i);	// sum mono and quad acc
  }
}
