/*
 * sphload.c: routines to create tree.  Public routines: maketree().
 * Copyright (c) 2016 by Joshua E. Barnes, Honolulu, HI.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "sphdefs.h"

//  Local routines and variables to perform tree construction.
//  __________________________________________________________
 
local void newtree(void);			// flush existing tree
local cellptr makecell(void);			// create an empty cell
local void expandbox(bodyptr, int);		// set size of root cell
local void loadbody(bodyptr);			// load body into tree
local int subindex(bodyptr, cellptr);		// compute subcell index
local void hackcofm(cellptr, real, int);	// find centers of mass
local void setrcrit(cellptr, vector, real);	// set cell's crit. radius
local void threadtree(nodeptr, nodeptr);	// set next and more links
local void hackquad(cellptr);			// compute quad moments
 
local bool bh86, sw94, theff;			// use alternate criteria
local nodeptr freecell = NULL;			// list of free cells

#define MAXLEVEL  48				// max depth of tree 

local int cellhist[MAXLEVEL];			// count cells by level
local int subnhist[MAXLEVEL];			// count subnodes by level

//  maketree: initialize tree structure for hierarchical force calculation.
//  _______________________________________________________________________
 
void maketree(bodyptr btab, int nbody)
{
  double cpustart;
  bodyptr p;
  int i;
 
  cpustart = cputime();				// record time at start
  newtree();					// flush existing tree, etc
  root = makecell();				// allocate the root cell
  CLRV(Pos(root));				// initialize the midpoint
  expandbox(btab, nbody);			// and expand cell to fit
  for (p = btab; p < btab+nbody; p++)		// loop over all bodies
    loadbody(p);				// insert each into tree
  bh86 = scanopt(options, "bh86");		// set flags for alternate
  sw94 = scanopt(options, "sw94");		// ...cell opening criteria
  theff = scanopt(options, "theta-eff");
  if ((bh86 && sw94) || (sw94 && theff) || (theff && bh86))
						// allow just one at a time
    error("%s.maketree: pick one of bh86, sw94, theta-eff\n", getprog());
  tdepth = 0;					// init count of levels
  for (i = 0; i < MAXLEVEL; i++)		// and init tree histograms
    cellhist[i] = subnhist[i] = 0;
  hackcofm(root, rsize, 0);			// find c-of-m coords, etc
  threadtree((nodeptr) root, NULL);		// add next and more links
  if (usequad)					// if including quad terms
    hackquad(root);				// find quadrupole moments
  cputree = cputime() - cpustart;		// store elapsed CPU time
}

//  newtree: reclaim cells in tree, prepare to build new one.
//  _________________________________________________________
 
local void newtree(void)
{
  static bool firstcall = TRUE;
  nodeptr p;
 
  if (! firstcall) {				// if cells to reclaim
    p = (nodeptr) root;				// start with the root
    while (p != NULL)				// loop scanning tree
      if (Cell(p)) {				// if we found a cell
	Next(p) = freecell;			// then save existing list
	freecell = p;				// and add it to the front
	p = More(p);				// then scan down tree
      } else					// else, skip over bodies
	p = Next(p);				// by going on to the next
  } else					// else nothing to reclaim
    firstcall = FALSE;				// so just note it
  root = NULL;					// flush existing tree
  ncell = 0;					// reset cell count
}
 
//  makecell: return pointer to free cell.
//  ______________________________________
 
local cellptr makecell(void)
{
  cellptr c;
  int i;
 
  if (freecell == NULL)				// if no free cells left
    c = (cellptr) allocate(sizeof(cell));	// then allocate a new one
  else {					// else use existing cell
    c = (cellptr) freecell;			// take the one in front
    freecell = Next(c);				// and go on to next one
  }
  Type(c) = CELL;				// initialize node type
  Flags(c) = 0;					// and force update flag
  for (i = 0; i < NSUB; i++)			// loop over subcells
    Subp(c)[i] = NULL;				// and empty each one
  ncell++;					// count one more cell
  return (c);					// return pointer to cell
}

//  expandbox: find range of coordinate values (with respect to root)
//  and expand root cell to fit.  The size is doubled at each step to
//  take advantage of exact representation of powers of two.
//  _________________________________________________________________
 
local void expandbox(bodyptr btab, int nbody)
{
  real dmax, d;
  bodyptr p;
  int k;
 
  dmax = 0.0;					// keep track of max value
  for (p = btab; p < btab+nbody; p++)		// loop over all bodies
    for (k = 0; k < NDIM; k++) {		// and over all dimensions
      d = rabs(Pos(p)[k] - Pos(root)[k]);	// find distance to midpnt
      if (d > dmax)				// if bigger than old one
	dmax = d;				// store new max value
    }
  while (rsize < 2 * dmax)			// loop until value fits
    rsize = 2 * rsize;				// doubling box each time
}

//  loadbody: descend tree and insert body p in appropriate place.
//  ______________________________________________________________
 
local void loadbody(bodyptr p)
{
  cellptr q, c;
  int qind, k;
  real qsize;
 
  q = root;					// start with tree root
  qind = subindex(p, q);			// get index of subcell
  qsize = rsize;				// keep track of cell size
  while (Subp(q)[qind] != NULL) {		// loop descending tree
    if (Body(Subp(q)[qind])) {			// if reached a "leaf"
      if (Pos(p)[0] == Pos(Subp(q)[qind])[0] &&
	  Pos(p)[1] == Pos(Subp(q)[qind])[1] &&
	  Pos(p)[2] == Pos(Subp(q)[qind])[2])
	fatal("%s.loadbody: two bodies have identical position (%a,%a,%a)\n",
	      getprog(), Pos(p)[0], Pos(p)[1], Pos(p)[2]);
      c = makecell();				// then allocate new cell
      for (k = 0; k < NDIM; k++)		// and initialize midpoint
	Pos(c)[k] = Pos(q)[k] +			// offset from parent
	  (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;
      Subp(c)[subindex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
						// put body in cell
      Subp(q)[qind] = (nodeptr) c;		// link cell in tree
    }
    q = (cellptr) Subp(q)[qind];		// advance to next level
    qind = subindex(p, q);			// get index to examine
    qsize = qsize / 2;				// shrink current cell
  }
  Subp(q)[qind] = (nodeptr) p;			// found place, store p
}

//  subindex: compute subcell index for body p in cell q.
//  _____________________________________________________
 
local int subindex(bodyptr p, cellptr q)
{
  int ind, k;
 
  ind = 0;					// accumulate subcell index
  for (k = 0; k < NDIM; k++)			// loop over dimensions
    if (Pos(q)[k] <= Pos(p)[k])			// if beyond midpoint
      ind += NSUB >> (k + 1);			// then skip over subcells
  return (ind);
}

//  hackcofm: descend tree finding center-of-mass coordinates;
//  also sets critical cell radii, if appropriate.
//  __________________________________________________________
 
local void hackcofm(cellptr p, real psize, int lev)
{
  vector cmpos, tmpv;
  nodeptr q;
 
  if (lev >= MAXLEVEL)				// limit depth of tree
    fatal("%s.hackcofm: tree depth exceeds maximum!\n", getprog());
  tdepth = MAX(tdepth, lev);			// remember maximum level
  cellhist[lev]++;				// count cells by level
  Mass(p) = 0.0;				// init cell's total mass
  CLRV(cmpos);					// and center of mass pos
  for (int i = 0; i < NSUB; i++)		// loop over the subnodes
    if ((q = Subp(p)[i]) != NULL) {		// skipping the NULLs
      subnhist[lev]++;				// count existing subnodes
      if (Cell(q))				// and if node is a cell
	hackcofm((cellptr) q, psize/2, lev+1);	// then do the same for it
      Flags(p) |= Flags(q);			// propagate flags upward
      Mass(p) += Mass(q);			// accumulate total mass
      MULVS(tmpv, Pos(q), Mass(q));		// weight position by mass
      ADDV(cmpos, cmpos, tmpv);			// and sum c-of-m position
    }
  if (Mass(p) > 0.0) {				// usually, cell has mass
    DIVVS(cmpos, cmpos, Mass(p));		// so find c-of-m position
  } else {					// but if no mass inside
    SETV(cmpos, Pos(p));			// use geo. center for now
  }
  for (int k = 0; k < NDIM; k++)		// check center-of-mass
    if (! ((Pos(p)[k] - psize/2 <= cmpos[k]) &&
	   (cmpos[k] <= Pos(p)[k] + psize/2))) {
      eprintf("[%s.hackcofm: warning: CM out of bounds: psize = %a  lev = %d\n"
	    "  k = %d:  (%.7a <= %.7a <= %.7a) fails]\n", getprog(), psize,
	    lev, k, Pos(p)[k] - psize/2, cmpos[k], Pos(p)[k] + psize/2);
      for (int i = 0; i < NSUB; i++)
	if ((q = Subp(p)[i]) != NULL)
	  eprintf("[%s.hackcofm: warning: subcell %d: %s  mass = %a"
		  "  pos = (%.6a,%.6a,%.6a)]\n", getprog(), i,
		  Cell(q) ? "cell" : "body", Mass(q),
		  Pos(q)[0], Pos(q)[1], Pos(q)[2]);
    }
  setrcrit(p, cmpos, psize);			// set critical radius
  SETV(Pos(p), cmpos);				// and center-of-mass pos
}

//  setrcrit: assign critical radius for cell p, using center-of-mass
//  position cmpos and cell size psize.
//  _________________________________________________________________

local void setrcrit(cellptr p, vector cmpos, real psize)
{
  real bmax2, d;
  int k;

if (theta == 0.0)				// if exact calculation
    Rcrit2(p) = rsqr(2 * rsize);		// then always open cells
  else if (sw94) {				// if using S&W's criterion
    bmax2 = 0.0;				// compute max distance^2
    for (k = 0; k < NDIM; k++) {		// loop over dimensions
      d = cmpos[k] - Pos(p)[k] + psize/2;	// get dist from corner
      bmax2 += rsqr(MAX(d, psize - d));		// and sum max distance^2
    }
    Rcrit2(p) = bmax2 / rsqr(theta);		// use max dist from cm
  } else if (bh86) {				// if using old criterion
    Rcrit2(p) = rsqr(psize / theta);		// then use size of cell
  } else if (theff) {				// if using effective theta
    DISTV(d, cmpos, Pos(p));			// find offset from center
    Rcrit2(p) = rsqr(psize / (theta * (eps + psize) / (2*eps + psize)) + d);
  } else {					// else use new criterion
    DISTV(d, cmpos, Pos(p));			// find offset from center
    Rcrit2(p) = rsqr(psize / theta + d);	// use size plus offset
  }
}

//  threadtree: do a recursive treewalk starting from node p,
//  with next stop n, installing Next and More links.
//  _________________________________________________________
 
local void threadtree(nodeptr p, nodeptr n)
{
  int ndesc, i;
  nodeptr desc[NSUB+1];
 
  Next(p) = n;					// set link to next node
  if (Cell(p)) {				// if descendents to thread
    ndesc = 0;					// start counting them
    for (i = 0; i < NSUB; i++)			// loop over all subcells
      if (Subp(p)[i] != NULL)			// if this one is occupied
	desc[ndesc++] = Subp(p)[i];		// then store it in table
    More(p) = desc[0];				// set link to 1st one
    desc[ndesc] = n;				// thread last one to next
    for (i = 0; i < ndesc; i++)			// loop over descendents
      threadtree(desc[i], desc[i+1]);		// and thread them together
  }
}

//  hackquad: descend tree, evaluating quadrupole moments.  Note that this
//  routine is coded so that the Subp() and Quad() components of a cell can
//  share the same memory locations.
//  _______________________________________________________________________
 
local void hackquad(cellptr p)
{
  int ndesc, i;
  nodeptr desc[NSUB], q;
  vector dr;
  matrix Qtmp, trQM;
  real trQ;

  ndesc = 0;					// count occupied subnodes
  for (i = 0; i < NSUB; i++)			// loop over all subnodes
    if (Subp(p)[i] != NULL)			// if this one's occupied
      desc[ndesc++] = Subp(p)[i];		// copy it to safety
  CLRM(Quad(p));				// init quadrupole moment
  for (i = 0; i < ndesc; i++) {			// loop over real subnodes
    q = desc[i];				// access ech one in turn
    if (Cell(q))				// if it's also a cell
      hackquad((cellptr) q);			// then process it first
    SUBV(dr, Pos(q), Pos(p));			// find displacement vect.
    OUTVP(Qtmp, dr, dr);			// form outer prod. of dr
    MULMS(Qtmp, Qtmp, Mass(q));			// scale by mass of subnode
#if defined(NOSOFTCORR)
    TRACEM(trQ, Qtmp);				// form trace (= dot prod.)
    SETMI(trQM);				// init unit matrix
    MULMS(trQM, trQM, trQ/3.0);			// and scale by trace
    SUBM(Qtmp, Qtmp, trQM);			// get traceless quad moment
#endif
    if (Cell(q))				// if subnode is cell
      ADDM(Qtmp, Qtmp, Quad(q));		// then include its moment
    ADDM(Quad(p), Quad(p), Qtmp);		// sum cell's quad moment
  }
}

