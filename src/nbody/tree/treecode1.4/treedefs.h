/****************************************************************************/
/* TREEDEFS.H: include file for hierarchical force calculation routines.    */
/* These definitions are needed for treeload.c and treegrav.c, but this     */
/* file does not provide definitions for other parts of the N-body code.    */
/* Copyright (c) 2000 by Joshua E. Barnes, Honolulu, HI.                    */
/****************************************************************************/
 
#ifndef _treedefs_h
#define _treedefs_h

/*
 * Body and cell data structures are used to represent the tree.  During
 * tree construction, descendent pointers are stored in the subp arrays:
 *
 *          +-------------------------------------------------------------+
 * root --> | CELL: mass, pos, next, rcrit2, more, subp:[/,o,/,/,/,/,o,/] |
 *          +----------------------------------------------|---------|----+
 *                                                         |         |
 *     +---------------------------------------------------+         |
 *     |                                                             |
 *     |    +--------------------------------------+                 |
 *     +--> | BODY: mass, pos, next, vel, acc, phi |                 |
 *          +--------------------------------------+                 |
 *                                                                   |
 *     +-------------------------------------------------------------+
 *     |
 *     |    +-------------------------------------------------------------+
 *     +--> | CELL: mass, pos, next, rcrit2, more, subp:[o,/,/,o,/,/,o,/] |
 *          +--------------------------------------------|-----|-----|----+
 *                                                      etc   etc   etc
 *
 */
 
/* <A NAME="node"></A>
 * NODE: data common to BODY and CELL structures.
 */

typedef struct _node {
    short type;                 /* code for node type */
    bool update;		/* status in force calc */
    struct _node *next;		/* link to next force calc */   
    real mass;                  /* total mass of node */
    vector pos;                 /* position of node */
} node, *nodeptr;
 
#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Next(x)   (((nodeptr) (x))->next)
#define Mass(x)   (((nodeptr) (x))->mass)
#define Pos(x)    (((nodeptr) (x))->pos)

#define BODY 01                 /* type code for bodies */
#define CELL 02                 /* type code for cells */

/* <A NAME="body"></A>
 * BODY: data structure used to represent particles.
 */
 
typedef struct {
    node bodynode;              /* data common to all nodes */
    vector vel;                 /* velocity of body */
    vector acc;                 /* acceleration of body */
    real phi;                   /* potential at body */
} body, *bodyptr;

#define Vel(x)    (((bodyptr) (x))->vel)
#define Acc(x)    (((bodyptr) (x))->acc)
#define Phi(x)    (((bodyptr) (x))->phi)

/* <A NAME="cell"></A>
 * CELL: structure used to represent internal nodes of tree.
 */
  
#define NSUB (1 << NDIM)        /* subcells per cell */
 
typedef struct {
    node cellnode;              /* data common to all nodes */
#ifndef QUICKSCAN
    real rcrit2;                /* critical c-of-m radius^2 */
#endif
    union {
	nodeptr more;		/* link to first descendent */   
	real dlast;		/* dist from last force sum */
    } mord;
    union {
	nodeptr subp[NSUB];     /* descendents of cell */
	matrix quad;            /* quad. moment of cell */
    } sorq;
} cell, *cellptr;
 
#ifndef QUICKSCAN
#define Rcrit2(x) (((cellptr) (x))->rcrit2)
#endif

#define More(x)   (((cellptr) (x))->mord.more)
#define Dlast(x)  (((cellptr) (x))->mord.dlast)
#define Subp(x)   (((cellptr) (x))->sorq.subp)
#define Quad(x)   (((cellptr) (x))->sorq.quad)

/* <A NAME="global"></A>
 * GLOBAL: pseudo-keyword for storage class.
 */

#ifndef global
#  define global extern
#endif
 
/*
 * Parameters for tree construction and force calculation.
 */
 
#ifndef QUICKSCAN
global real theta;                      /* force accuracy parameter         */
#endif

global string options;                  /* various option keywords          */
 
global bool usequad;		        /* use quadrupole corrections       */

global real eps;                        /* density smoothing parameter      */

/*
 * Tree construction.
 */
 
void maketree(bodyptr, int);		/* construct tree structure         */

global cellptr root;                    /* pointer to root cell             */
global real rsize;                      /* side-length of root cell         */
global int ncell;			/* count of cells in tree           */
global int tdepth;			/* count of levels in tree          */
global real cputree;			/* CPU time to build tree	    */
 
/*
 * Force calculation.
 */

void gravcalc(void);			/* update force on bodies           */

global int nfcalc;			/* total force calculations         */
global int nbbcalc;			/* total body-body interactions     */
global int nbccalc;			/* total body-cell interactions     */
global real cpuforce;			/* CPU time for force calc	    */

#endif /* ! _treedefs_h */
