/*
 * treedefs.h: include file for hierarchical force calculation routines.
 * These definitions are needed for treeload.c and treegrav.c, but this
 * file does not provide definitions for other parts of the N-body code.
 * Copyright (c) 2015 by Joshua E. Barnes, Kyoto, Japan.
 */
 
#ifndef _treedefs_h
#define _treedefs_h

//  Body and cell data structures are used to represent the tree.  During
//  tree construction, descendent pointers are stored in the subp arrays:
// 
//           +-------------------------------------------------------------+
//  root --> | CELL: mass, pos, next, rcrit2, more, subp:[/,o,/,/,/,/,o,/] |
//           +----------------------------------------------|---------|----+
//                                                          |         |
//      +---------------------------------------------------+         |
//      |                                                             |
//      |    +--------------------------------------+                 |
//      +--> | BODY: mass, pos, next, vel, acc, phi |                 |
//           +--------------------------------------+                 |
//                                                                    |
//      +-------------------------------------------------------------+
//      |
//      |    +-------------------------------------------------------------+
//      +--> | CELL: mass, pos, next, rcrit2, more, subp:[o,/,/,o,/,/,o,/] |
//           +--------------------------------------------|-----|-----|----+
//                                                       etc   etc   etc
 
//  node: data common to body and cell structures.
//  ______________________________________________

typedef struct _node {
  short type;                   // code for node type
  bool update;		        // status in force calc
  struct _node *next;		// link to next force calc
  real mass;                    // total mass of node
  vector pos;                   // position of node
} node, *nodeptr;
 
#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Next(x)   (((nodeptr) (x))->next)
#define Mass(x)   (((nodeptr) (x))->mass)
#define Pos(x)    (((nodeptr) (x))->pos)

#define BODY 01                 // type code for bodies
#define CELL 02                 // type code for cells

//  body: data structure used to represent particles.
//  _________________________________________________
 
typedef struct {
  node bodynode;                // data common to all nodes
  vector vel;                   // velocity of body
  vector acc;                   // acceleration of body
  real phi;                     // potential at body
} body, *bodyptr;

#define Vel(x)    (((bodyptr) (x))->vel)
#define Acc(x)    (((bodyptr) (x))->acc)
#define Phi(x)    (((bodyptr) (x))->phi)

//  cell: structure used to represent internal nodes of tree.
//  _________________________________________________________
  
#define NSUB (1 << NDIM)        // subcells per cell
 
typedef struct {
  node cellnode;                // data common to all nodes
#if !defined(QUICKSCAN)
  real rcrit2;                  // critical c-of-m radius^2
#endif
  union {
    nodeptr more;		// link to first descendent
    real trace;			// saved trace for softening correction
  } mort;
  union {
    nodeptr subp[NSUB];         // descendents of cell
    matrix quad;                // quad. moment of cell
  } sorq;
} cell, *cellptr;
 
#if !defined(QUICKSCAN)
#define Rcrit2(x) (((cellptr) (x))->rcrit2)
#endif

#define More(x)   (((cellptr) (x))->mort.more)
#define Trace(x)  (((cellptr) (x))->mort.trace)
#define Subp(x)   (((cellptr) (x))->sorq.subp)
#define Quad(x)   (((cellptr) (x))->sorq.quad)

//  global: pseudo-keyword for storage class.
//  _________________________________________

#if !defined(global)
#  define global extern
#endif
 
//  Parameters for tree construction and force calculation.
//  _______________________________________________________
 
#if !defined(QUICKSCAN)
global real theta;              // force accuracy parameter
#endif

global string options;          // various option keywords
 
global bool usequad;		// use quadrupole corrections

global real eps;                // density smoothing parameter

//  Tree construction: interface and global results.
//  ________________________________________________
 
void maketree(bodyptr, int);	// construct tree structure

global cellptr root;            // pointer to root cell
global real rsize;              // side-length of root cell
global int ncell;		// count of cells in tree
global int tdepth;		// count of levels in tree
global real cputree;		// CPU time to build tree	
 
//  Force calculation: interface and diagnostics.
//  _____________________________________________

void gravcalc(void);		// update force on bodies

global int nfcalc;		// total force calculations
global long nbbcalc;		// total body-body interactions
global long nbccalc;		// total body-cell interactions
global real cpuforce;		// CPU time for force calc

#endif // ! _treedefs_h */
