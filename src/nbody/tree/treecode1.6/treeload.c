/*
 * treeload.c: routines to create tree.	 Public routines: maketree().
 * Copyright (c) 2015 by Joshua E. Barnes, Kyoto, Japan.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "treedefs.h"

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

#define MAXLEVEL  48				// max height of tree

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
  newtree();					// flush existing tree
  root = makecell();				// allocate the root cell
  CLRV(Pos(root));				// initialize the midpoint
  expandbox(btab, nbody);			// expand to fit bodies
  for (p = btab; p < btab+nbody; p++)		// loop over all bodies
    loadbody(p);				// insert each into tree
  bh86 = scanopt(options, "bh86");		// set flags for alternate
  sw94 = scanopt(options, "sw94");		// cell opening criteria
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
      if (Type(p) == CELL) {			// if we found a cell
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
  Update(c) = FALSE;				// and force update flag
  for (i = 0; i < NSUB; i++)			// loop over subcells
    Subp(c)[i] = NULL;				// and empty each one
  ncell++;					// count cells in use
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
    if (Type(Subp(q)[qind]) == BODY) {		// found a body in subcell?
      if (Pos(p)[0] == Pos(Subp(q)[qind])[0] &&
	  Pos(p)[1] == Pos(Subp(q)[qind])[1] &&
	  Pos(p)[2] == Pos(Subp(q)[qind])[2])
	fatal("%s.loadbody: two bodies have identical position (%a,%a,%a)\n",
	      getprog(), Pos(p)[0], Pos(p)[1], Pos(p)[2]);
      c = makecell();				// allocate cell for both
      for (k = 0; k < NDIM; k++)		// and initialize midpoint
	Pos(c)[k] = Pos(q)[k] +	(Pos(p)[k]<Pos(q)[k] ? - qsize : qsize) / 4;
						// set offset from parent
      Subp(c)[subindex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
						// put old body in cell
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
  int i, k;
  nodeptr q;
  real phalf = psize / 2.0;
 
  if (lev >= MAXLEVEL)				// limit depth of tree
    fatal("%s.hackcofm: tree depth exceeds maximum!\n", getprog());
  tdepth = MAX(tdepth, lev);			// remember maximum level
  cellhist[lev]++;				// count cells by level
  Mass(p) = 0.0;				// init cell's total mass
  Ndesc(p) = 0;					// zero descendents count
  CLRV(cmpos);					// and center of mass pos
  for (i = 0; i < NSUB; i++)			// loop over the subnodes
    if ((q = Subp(p)[i]) != NULL) {		// skipping empty ones
      subnhist[lev]++;				// count existing subnodes
      if (Type(q) == CELL)			// and if node is a cell
	hackcofm((cellptr) q, phalf, lev+1);	// then do the same for it
      Update(p) |= Update(q);			// propagate update request
      Mass(p) += Mass(q);			// accumulate total mass
      Ndesc(p) += (Type(q) == CELL ? Ndesc(q) : 1);
      ADDMULVS(cmpos, Pos(q), Mass(q));		// and center of mass posn
    }
  if (Mass(p) > 0.0) {				// usually, cell has mass
    DIVVS(cmpos, cmpos, Mass(p));		// so find c-of-m position
  } else {					// but if no mass inside
    SETV(cmpos, Pos(p));			// use centroid position
  }
  SUBV(tmpv, cmpos, Pos(p));
  if (ABS(tmpv[0])>phalf || ABS(tmpv[1])>phalf || ABS(tmpv[2])>phalf)
    eprintf("[%s.hackcofm: WARNING: center of mass outside cell:\n"
	    "  psize = %14.8g  lev = %d\n"
	    "  midp = (%14.8g,%14.8g,%14.8g)\n"
	    "  cofm = (%14.8g,%14.8g,%14.8g)]\n", getprog(), psize, lev,
	    Pos(p)[0], Pos(p)[1], Pos(p)[2], cmpos[0], cmpos[1], cmpos[2]);
#if !defined(QUICKSCAN)
  setrcrit(p, cmpos, psize);			// set critical radius
#endif
  SETV(Pos(p), cmpos);				// and center-of-mass pos
}

#if !defined(QUICKSCAN)

//  setrcrit: assign critical radius for cell p, using center-of-mass
//  position cmpos and cell size psize.
//  _________________________________________________________________

local void setrcrit(cellptr p, vector cmpos, real psize)
{
  real bmax2, d;
  int k;

  if (theta == 0.0) {				// if exact calculation
    Rcrit2(p) = rsqr(2 * rsize);		// then always open cells
  } else if (sw94) {				// if using S&W's criterion
    bmax2 = 0.0;				// compute max distance^2
    for (k = 0; k < NDIM; k++) {		// loop over dimensions
      d = cmpos[k] - Pos(p)[k] + psize/2; 	// get dist from corner
      bmax2 += rsqr(MAX(d, psize - d));		// and sum max distance^2
    }
    Rcrit2(p) = bmax2 / rsqr(theta);		// use max dist from cm
  } else if (bh86) {				// if using old criterion
    Rcrit2(p) = rsqr(psize / theta);		// just use size of cell
  } else if (theff) {				// if using effective theta
    DISTV(d, cmpos, Pos(p));			// find offset from center
    Rcrit2(p) = rsqr(psize / (theta * (eps + psize) / (2*eps + psize)) + d);
  } else {					// else use new criterion
    DISTV(d, cmpos, Pos(p));			// find offset from center
    Rcrit2(p) = rsqr(psize / theta + d);	// use size plus offset
  }
}

#endif

//  threadtree: do a recursive treewalk starting from node p,
//  with next stop n, installing Next and More links.
//  _________________________________________________________
 
local void threadtree(nodeptr p, nodeptr n)
{
  int ndesc, i;
  nodeptr desc[NSUB+1];
 
  Next(p) = n;					// set link to next node
  if (Type(p) == CELL) {			// if descendents to thread
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
 
  ndesc = 0;					// list descendent subnodes
  for (i = 0; i < NSUB; i++)			// loop over all subnodes
    if (Subp(p)[i] != NULL)			// if this one's occupied
      desc[ndesc++] = Subp(p)[i];		// copy it to safety
  CLRM(Quad(p));				// init quadrupole moment
  for (i = 0; i < ndesc; i++) {			// loop over real subnodes
    q = desc[i];				// access ech one in turn
    if (Type(q) == CELL)			// if it's also a cell
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
    if (Type(q) == CELL)			// if subnode is cell
      ADDM(Qtmp, Qtmp, Quad(q));		// then include its moment
    ADDM(Quad(p), Quad(p), Qtmp);		// sum cell's quad moment
  }
}
