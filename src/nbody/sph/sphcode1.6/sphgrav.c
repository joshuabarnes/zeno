/*
 * sphgrav.c: routines to compute gravity. Public routines: gravforce(),
 * report_force(). Copyright (c) 2016 by Joshua E. Barnes, Honolulu, HI.
 */
 
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "sphdefs.h"
 
//  Local routines to perform force calculations.
//  _____________________________________________

local void walkgrav(nodeptr *, nodeptr *, cellptr, cellptr,
		    nodeptr, real, vector);
local bool accept(nodeptr, real, vector);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector);
local void gravsum(bodyptr, cellptr, cellptr);
local void sumnode(cellptr, cellptr, vector, real *, vector);
local void sumcell(cellptr, cellptr, vector, real *, vector);
 
//  Lists of active nodes and interactions.
//  _______________________________________
 
#if !defined(FACTIVE)
#  define FACTIVE  1.0				// active list fudge factor
#endif
 
local int actlen;                               // length as allocated

local nodeptr *active;                          // list of nodes tested
local cellptr interact;                         // list of interactions

//  gravforce: perform force calculation on all particles.
//  ______________________________________________________
 
void gravforce(void)
{
  double cpustart;
  vector rmid;

  actlen = FACTIVE * 216 * (tdepth + 1);	// estimate list length
  actlen = actlen * rpow(theta, -2.5);		// allow for opening angle
  active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
  interact = (cellptr) allocate(actlen * sizeof(cell));
  cpustart = cputime();				// record time, less alloc
  actmax = nfcalc = nbbcalc = nbccalc = 0;	// zero cumulative counters
  active[0] = (nodeptr) root;			// initialize active list
  CLRV(rmid);					// set center of root cell
  walkgrav(active, active + 1, interact, interact + actlen,
	   (nodeptr) root, rsize, rmid);	// scan tree, update forces
  cpuforce = cputime() - cpustart;		// store CPU time w/o alloc
  free(active);					// free up temp storage
  free(interact);
}

//  report_force: print statistics on tree construction and force calculation.
//  __________________________________________________________________________

void gravreport(stream ostr, int nbody)
{
  fprintf(ostr, "\n%9s %9s %9s %9s %14s %14s %9s\n",
	  "rsize", "tdepth", "ftree",
	  "actmax", "nbbint", "nbcint", "CPUgfc");
  fprintf(ostr, "%9.1f %9d %9.3f %9d %14ld %14ld %9.3f\n",
	  rsize, tdepth, (nbody + ncell - 1) / ((real) ncell),
	  actmax, nbbcalc, nbccalc, cpuforce);
  fflush(ostr);
}

//  walkgrav: do a complete walk of the tree, building the interaction
//  list level-by-level and computing the resulting force on each body.
//  ___________________________________________________________________
 
local void walkgrav(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                    nodeptr p, real psize, vector pmid)
{
  nodeptr *np, *ap, q;
  int actsafe;
  matrix trQM;
 
  if (Update(p)) {				// are new forces needed?
    np = nptr;					// start new active list
    actsafe = actlen - NSUB;			// leave room for NSUB more 
    for (ap = aptr; ap < nptr; ap++)		// loop over active nodes
      if (Cell(*ap)) {				// is this node a cell?
	if (accept(*ap, psize, pmid)) {		// does it pass the test?
	  Mass(cptr) = Mass(*ap);		// copy to interaction list
	  SETV(Pos(cptr), Pos(*ap));
#if defined(NOSOFTCORR)
	    SETM(Quad(cptr), Quad(*ap));	// store traceless moment
#else
	    TRACEM(Trace(cptr), Quad(*ap));	// save value of trace
	    SETMI(trQM);
	    MULMS(trQM, trQM, Trace(cptr)/3.0);	// scale unit matrix by trace
	    SUBM(Quad(cptr), Quad(*ap), trQM);	// compute traceless moment
#endif
	  cptr++;				// and bump cell array ptr
	} else {				// else it fails the test
	  if (np - active >= actsafe)		// check list has room
	    fatal("%s.walkgrav: active list overflow\n", getprog());
	  for (q = More(*ap); q != Next(*ap); q = Next(q))
						// loop over all subcells
	    *np++= q;				// put on new active list
	}
      } else					// else this node is a body
	if (*ap != p) {				// if not self-interaction
	  --bptr;				// bump body array ptr
	  Mass(bptr) = Mass(*ap);		// and copy data to array
	  SETV(Pos(bptr), Pos(*ap));
	}
    actmax = MAX(actmax, np - active);		// keep track of max active
    if (np != nptr)				// if new actives listed
      walksub(nptr, np, cptr, bptr, p, psize, pmid);
						// then visit next level
    else {					// else no actives left, so
      if (! Body(p))				// must have found a body
	error("%s.walkgrav: recursion terminated with cell\n", getprog());
      gravsum((bodyptr) p, cptr, bptr);		// sum force on the body
    }
  }
}

//  accept: standard criterion accepts cell if its critical radius
//  does not intersect cell p, and also imposes above condition.
//  _______________________________________________________________
 
local bool accept(nodeptr c, real psize, vector pmid)
{
  real dmax, dsq, dk;
  int k;
 
  dmax = psize;					// init maximum distance
  dsq = 0.0;					// and squared min distance
  for (k = 0; k < NDIM; k++) {			// loop over space dims
    dk = Pos(c)[k] - pmid[k];			// form distance to midpnt
    if (dk < 0)					// and get absolute value
      dk = - dk;
    if (dk > dmax)				// keep track of max value
      dmax = dk;
    dk -= ((real) 0.5) * psize;			// allow for size of cell
    if (dk > 0)
      dsq += dk * dk;				// sum min dist to cell ^2
  }
  return (dsq > Rcrit2(c) &&			// test angular criterion
	  dmax > ((real) 1.5) * psize);		// and adjacency criterion
}

//  walksub: test next level's active list against subnodes of p.
//  _____________________________________________________________
 
local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid)
{
  real poff;
  nodeptr q;
  int k;
  vector nmid;
 
  poff = psize / 4;				// precompute mid. offset
  if (Cell(p)) {				// fanout over descendents
    for (q = More(p); q != Next(p); q = Next(q)) {
						// loop over all subcells
      for (k = 0; k < NDIM; k++)		// locate each's midpoint
	nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
      walkgrav(nptr, np, cptr, bptr, q, psize / 2, nmid);
						// recurse on subcell
    }
  } else {					// extend virtual tree
    for (k = 0; k < NDIM; k++)			// locate next midpoint
      nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
    walkgrav(nptr, np, cptr, bptr, p, psize / 2, nmid);
						// and search next level
  }
}

//  gravsum: compute gravitational field at body p0.
//  ________________________________________________
 
local void gravsum(bodyptr p0, cellptr cptr, cellptr bptr)
{
  vector pos0, acc0;
  real phi0;
 
  SETV(pos0, Pos(p0));				// copy position of body
  phi0 = 0.0;					// init total potential
  CLRV(acc0);					// and total acceleration
  if (usequad)					// if using quad moments?
    sumcell(interact, cptr, pos0, &phi0, acc0);	// sum cell forces w quads
  else						// not using quad moments?
    sumnode(interact, cptr, pos0, &phi0, acc0);	// sum using monopole only
  sumnode(bptr, interact + actlen, pos0, &phi0, acc0);
						// sum forces from bodies
  Phi(p0) = phi0;				// store total potential
  SETV(Acc(p0), acc0);				// and total acceleration
  nfcalc++;					// update counters
  nbbcalc += interact + actlen - bptr;
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
 
//  sumcell: add up body-cell interactions.
//  _______________________________________
 
local void sumcell(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
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
#if !defined(NOSOFTCORR)
    qdr2 -= eps2thrd * Trace(p);		// apply Keigo's correction
#endif
    dr5i = ((real) 3.0) * dr2i * dr2i * dr1i;	// factor 3 saves a multiply
    phi_q = ((real) 0.5) * dr5i * qdr2;		// form quad potential
    mdr3i += ((real) 5.0) * phi_q * dr2i;	// adjust radial term
    *phi0 -= mdr1i + phi_q;			// sum mono and quad pot
    ADDMULVS2(acc0, dr, mdr3i, qdr, - dr5i);	// sum mono and quad acc
  }
}
