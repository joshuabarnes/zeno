/****************************************************************************/
/* SPHGRAV.C: routines to compute gravity. Public routines: gravforce(),    */
/* report_force(). Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.    */
/****************************************************************************/
 
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "sphdefs.h"
 
/* Local routines to perform force calculations. */
 
local void walkgrav(nodeptr *, nodeptr *, cellptr, cellptr,
		    nodeptr, real, vector);
local bool accept(nodeptr, real, vector);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector);
local void gravsum(bodyptr, cellptr, cellptr);
local void sumnode(cellptr, cellptr, vector, real *, vector);
local void sumcell(cellptr, cellptr, vector, real *, vector);
 
/* Lists of active nodes and interactions. */
 
#if !defined(FACTIVE)
#define FACTIVE  1				/* active list fudge factor */
#endif
 
local int actlen;                               /* length as allocated      */

local nodeptr *active;                          /* list of nodes tested     */
local cellptr interact;                         /* list of interactions     */

/*
 * GRAVFORCE: perform force calculation on all particles.
 */
 
void gravforce(void)
{
    double cpustart;
    vector rmid;
 
    actlen = FACTIVE * 216 * tdepth;		/* estimate list length     */
    actlen = actlen * rpow(theta, -2.5);	/* allow for opening angle  */
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
    interact = (cellptr) allocate(actlen * sizeof(cell));
    cpustart = cputime();			/* record time, less alloc  */
    actmax = nfcalc = nbbcalc = nbccalc = 0;    /* zero cumulative counters */
    active[0] = (nodeptr) root;                 /* initialize active list   */
    CLRV(rmid);                                 /* set center of root cell  */
    walkgrav(active, active + 1, interact, interact + actlen,
             (nodeptr) root, rsize, rmid);      /* scan tree, update forces */
    cpuforce = cputime() - cpustart;		/* store CPU time w/o alloc */
    free(active);				/* free up temp storage	    */
    free(interact);
}

/*
 * REPORT_FORCE: print staristics on tree construction and force calculation.
 */

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

/*
 * WALKGRAV: do a complete walk of the tree, building the interaction
 * list level-by-level and computing the resulting force on each body.
 */
 
local void walkgrav(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                    nodeptr p, real psize, vector pmid)
{
    nodeptr *np, *ap, q;
    int actsafe;
 
    if (Update(p)) {				/* are new forces needed?   */
	np = nptr;				/* start new active list    */
	actsafe = actlen - NSUB;                /* leave room for NSUB more */
	for (ap = aptr; ap < nptr; ap++)        /* loop over active nodes   */
	    if (Cell(*ap)) {                    /* is this node a cell?     */
		if (accept(*ap, psize, pmid)) {	/* does it pass the test?   */
		    Mass(cptr) = Mass(*ap);	/* copy to interaction list */
		    SETV(Pos(cptr), Pos(*ap));
		    SETM(Quad(cptr), Quad(*ap));
		    cptr++;			/* and bump cell array ptr  */
		} else {			/* else it fails the test   */
		    if (np - active >= actsafe) /* check list has room      */
			error("walkgrav: active list overflow\n");
		    for (q = More(*ap); q != Next(*ap); q = Next(q))
						/* loop over all subcells   */
			*np++= q;               /* put on new active list   */
		}
	    } else                              /* else this node is a body */
		if (*ap != p) {                 /* if not self-interaction  */
		    --bptr;			/* bump body array ptr      */
		    Mass(bptr) = Mass(*ap);	/* and copy data to array   */
		    SETV(Pos(bptr), Pos(*ap));
		}
	actmax = MAX(actmax, np - active);	/* keep track of max active */
	if (np != nptr)                         /* if new actives listed    */
	    walksub(nptr, np, cptr, bptr, p, psize, pmid);
						/* then visit next level    */
	else {                                  /* else no actives left, so */
	    if (! Body(p))                      /* must have found a body   */
		error("walkgrav: recursion terminated with cell\n");
	    gravsum((bodyptr) p, cptr, bptr);   /* sum force on the body    */
	}
    }
}

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

/*
 * WALKSUB: test next level's active list against subnodes of p.
 */
 
local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid)
{
    real poff;
    nodeptr q;
    int k;
    vector nmid;
 
    poff = psize / 4;                           /* precompute mid. offset   */
    if (Cell(p)) {                              /* fanout over descendents  */ 
        for (q = More(p); q != Next(p); q = Next(q)) {
                                                /* loop over all subcells   */
            for (k = 0; k < NDIM; k++)          /* locate each's midpoint   */
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walkgrav(nptr, np, cptr, bptr, q, psize / 2, nmid);
                                                /* recurse on subcell       */
        }
    } else {                                    /* extend virtual tree      */
        for (k = 0; k < NDIM; k++)              /* locate next midpoint     */
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walkgrav(nptr, np, cptr, bptr, p, psize / 2, nmid);
                                                /* and search next level    */
    }
}

/*
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
    sumnode(bptr, interact + actlen, pos0, &phi0, acc0);
                                                /* sum forces from bodies   */
    Phi(p0) = phi0;                             /* store total potential    */
    SETV(Acc(p0), acc0);                        /* and total acceleration   */
    nfcalc++;                                   /* update counters          */
    nbbcalc += interact + actlen - bptr;
    nbccalc += cptr - interact;
}

/*
 * SUMNODE: add up body-node interactions.
 */
 
local void sumnode(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i;
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
        ADDMULVS(acc0, dr, mr3i);               /* sum partial acceleration */
    }
}
 
/*
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
