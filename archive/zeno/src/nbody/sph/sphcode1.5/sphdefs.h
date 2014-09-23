/****************************************************************************/
/* SPHDEFS.H: include file for hierarchical SPH force calculation routines. */
/* Copyright (c) 2009 by Joshua E. Barnes, Honolulu, Hawai'i.               */
/****************************************************************************/
 
#ifndef _sphdefs_h
#define _sphdefs_h
 
/*
 * NODE: data common to BODY and CELL structures.
 */

typedef struct _node {
    byte type;                  /* code for type of node */
    byte flags;			/* status in force calculation */
    byte curlevel;		/* level in timestep hierarchy */
    byte newlevel;		/* level set by latest force calc. */
    struct _node *next;		/* link to next in threaded tree */   
    real mass;                  /* total mass of node */
    vector pos;                 /* position of node */
} node, *nodeptr;
 
/* Accessor macros for nodes. */

#define Type(x)      (((nodeptr) (x))->type)
#define Flags(x)     (((nodeptr) (x))->flags)
#define CurLevel(x)  (((nodeptr) (x))->curlevel)
#define NewLevel(x)  (((nodeptr) (x))->newlevel)
#define Next(x)      (((nodeptr) (x))->next)
#define Mass(x)      (((nodeptr) (x))->mass)
#define Pos(x)       (((nodeptr) (x))->pos)

/* Bits used to construct type of a node. Note that low 4 bits are ignored. */

#define CELL  0x80		/* node is cell of oct-tree */
#define BODY  0x40		/* node is any kind of body */
#define GAS   0x20		/* node is smooth hydro body */
#define STAR  0x10		/* node has become a star */

#define Cell(x)  ((Type(x) & CELL) != 0)
#define Body(x)  ((Type(x) & BODY) != 0)
#define Gas(x)   ((Type(x) & GAS) != 0)
#define Star(x)  ((Type(x) & STAR) != 0)

/* Bits governing gravity and hydro updates. */

#define INCLUDE  0x01		/* include body in tree construction */
#define UPDATE   0x02		/* update body, or bodies within cell */
#define INQUE    0x04		/* body is listed in current priority que */
#define DONE     0x08		/* smoothing operation is complete */
#define SURFACE  0x10		/* body has optically thin neighbors */

#define Include(x)  ((Flags(x) & INCLUDE) != 0)
#define Update(x)   ((Flags(x) & UPDATE) != 0)
#define InQue(x)    ((Flags(x) & INQUE) != 0)
#define Done(x)     ((Flags(x) & DONE) != 0)
#define Surface(x)  ((Flags(x) & SURFACE) != 0)

#define SetFlag(x,f)  (Flags(x) |= (f))
#define ClrFlag(x,f)  (Flags(x) &= ~(f))

/*
 * BODY: structure representing both gas and collisionless bodies.
 */
 
typedef struct {
    node bodynode;              /* data common to all nodes		    */
    vector vel;                 /* velocity of body			    */
    vector vmid;		/* velocity at midpoint			    */
    vector acc;                 /* acceleration of body			    */
    real smooth;		/* smoothing length			    */
    real phi;                   /* gravitational potential at body	    */
    real rho;			/* interpolated gas density		    */
#if defined(ENTROPY)
    real entf;			/* thermodynamic variable: a(s)		    */
#else
    real uint;			/* thermodynamic variable: u		    */
#endif
    real udotint;		/* internal change in specific energy	    */
#if defined(RADIATING)
    real udotrad;		/* radiative change in specific energy	    */
#endif
#if defined(COMPVISC)
    real udotvis;		/* viscous change in specific energy	    */
#endif
    real press;			/* gas pressure				    */
    real freq;			/* inverse of time-step			    */
#if defined(DIFFUSING) || defined(OPAQUE)
    real tau;			/* optical depth estimate		    */
#endif
#if defined(STARFORM) || defined(MASSLOSS)
    real birth;			/* time of star formation event		    */
#endif
#if defined(MASSLOSS)
    real death;			/* time particle converts to gas	    */
#endif
} body, *bodyptr;

/* Accessor macros for bodies.  Note that these are defined for all fields, */
/* even though only some of these fields may exist in any given version.    */

#define Vel(x)       (((bodyptr) (x))->vel)
#define Vmid(x)      (((bodyptr) (x))->vmid)
#define Acc(x)       (((bodyptr) (x))->acc)
#define Phi(x)       (((bodyptr) (x))->phi)
#define Smooth(x)    (((bodyptr) (x))->smooth)
#define Rho(x)       (((bodyptr) (x))->rho)
#define EntFunc(x)   (((bodyptr) (x))->entf)
#define Uintern(x)   (((bodyptr) (x))->uint)
#define UdotInt(x)   (((bodyptr) (x))->udotint)
#define UdotRad(x)   (((bodyptr) (x))->udotrad)
#define UdotVis(x)   (((bodyptr) (x))->udotvis)
#define Press(x)     (((bodyptr) (x))->press)
#define Freq(x)      (((bodyptr) (x))->freq)
#define Tau(x)       (((bodyptr) (x))->tau)
#define Birth(x)     (((bodyptr) (x))->birth)
#define Death(x)     (((bodyptr) (x))->death)

#define NthBody(bp,n)  ((bp) + (n))

/*
 * CELL: structure used to represent internal nodes of oct-tree.
 */
  
#define NSUB (1 << NDIM)        /* subcells per cell */
 
typedef struct {
    node cellnode;              /* data common to all nodes */
    real rcrit2;                /* critical c-of-m radius^2 */
    nodeptr more;		/* link to first descendent */   
    union {
	nodeptr subp[NSUB];     /* descendents of cell */
	matrix quad;            /* quad. moment of cell */
    } sorq;
} cell, *cellptr;
 
/* Accessor macros for cells. */

#define Rcrit2(x)    (((cellptr) (x))->rcrit2)
#define More(x)      (((cellptr) (x))->more)
#define Subp(x)      (((cellptr) (x))->sorq.subp)
#define Quad(x)      (((cellptr) (x))->sorq.quad)

/*
 * GLOBAL: pseudo-keyword for storage class.
 */

#if !defined(global)
#define global extern
#endif
 
/*
 * Parameters for tree construction and force calculation.
 */
 
global real theta;                      /* force accuracy parameter         */
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

void gravforce(void);			/* update force on bodies           */
void gravreport(stream, int);		/* report on force calculation      */

global int actmax;			/* maximum length of active list    */
global int nfcalc;			/* total force calculations         */
global long nbbcalc;			/* total body-body interactions     */
global long nbccalc;			/* total body-cell interactions     */
global real cpuforce;			/* CPU time for force calc	    */

#endif /* ! _sphdefs_h */
