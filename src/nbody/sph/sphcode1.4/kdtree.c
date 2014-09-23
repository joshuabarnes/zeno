/****************************************************************************/
/* KDTREE.C: routines to build kd-tree for SPH.				    */
/* Public routines: init_kdtree(), finish_kdtree(), build_kdtree().	    */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/* Adapted from SMOOTH V2.01; see http://www-hpcc.astro.washington.edu/     */
/****************************************************************************/

#include "stdinc.h"
#include "getparam.h"
#include "mathfns.h"
#include "vectdefs.h"
#include <assert.h>

#include "sphdefs.h"
#include "kdtree.h"

local void set_bounds(bound *, bodyptr *, int, int);
local int  median_index(bodyptr *, int, int, int);
local void combine_nodes(kdnode *, kdnode *, kdnode *);
local void upward_pass(kdxptr, int);

kdxptr init_kdtree(bodyptr btab, int nbody, int ngas)
{
    kdxptr kd;
    int i, j;

    kd = (kdxptr) allocate(sizeof(kdcontext));
    kd->ngas = ngas;
    kd->bptr = (bodyptr *) allocate(ngas * sizeof(bodyptr));
    for (i = j = 0; i < nbody; i++)
	if (Gas(NthBody(btab, i)))
	    kd->bptr[j++] = NthBody(btab, i);
    assert(j == ngas);
    set_bounds(&kd->bnd, kd->bptr, 0, kd->ngas - 1);
    return (kd);
}

void finish_kdtree(kdxptr kd)
{
    free(kd->bptr);
    free(kd->ntab);
    free(kd);
}

void build_kdtree(kdxptr kd, int nbucket)
{
    int k, n, i, d, m, j, ct;
    kdnode *ntab;

    n = kd->ngas;
    k = 1;
    while (n > nbucket) {
	n = n>>1;
	k = k<<1;
    }
    kd->nnode = k<<1;
    kd->nsplit = k;
    ntab = kd->ntab = (kdnode *) allocate(kd->nnode * sizeof(kdnode));
    ntab[KDROOT].first = 0;			/* initialize root node	    */
    ntab[KDROOT].last = kd->ngas-1;
    ntab[KDROOT].bnd = kd->bnd;
    i = KDROOT;
    ct = KDROOT;
    SetNext(ct);
    for ( ; ; ) {				/* loop splitting nodes	    */
	if (i < kd->nsplit) {
	    d = 0;				/* find longest dimension   */
	    for (j=1; j<3; ++j) {
		if (ntab[i].bnd.maxb[j]-ntab[i].bnd.minb[j] > 
		      ntab[i].bnd.maxb[d]-ntab[i].bnd.minb[d])
		    d = j;
	    }
	    m = median_index(kd->bptr, d, ntab[i].first, ntab[i].last);
	    ntab[i].dim = d;
	    ntab[i].split = Pos(kd->bptr[m])[d];
	    ntab[Lower(i)].bnd = ntab[i].bnd;
	    ntab[Lower(i)].bnd.maxb[d] = ntab[i].split;
	    ntab[Lower(i)].first = ntab[i].first;
	    ntab[Lower(i)].last = m-1;
	    ntab[Upper(i)].bnd = ntab[i].bnd;
	    ntab[Upper(i)].bnd.minb[d] = ntab[i].split;
	    ntab[Upper(i)].first = m;
	    ntab[Upper(i)].last = ntab[i].last;
	    i = Lower(i);
	} else {
	    ntab[i].dim = -1;
	    SetNext(i);
	    if (i == ct) break;
	}
    }
    upward_pass(kd, KDROOT);
}

/*
 * SET_BOUNDS: compute bounds from body pointers in specified range.
 */

local void set_bounds(bound *bndptr, bodyptr *bptr, int lo, int hi)
{
    int k, i;
    bound bnd;

    for (k = 0; k < 3; ++k)			/* initialize bounds	    */
	bnd.maxb[k] =  bnd.minb[k] = Pos(bptr[lo])[k];
    for (i = lo + 1; i <= hi; ++i) {		/* find actual bounds	    */
	for (k = 0; k < 3; ++k) {
	    if (bnd.minb[k] > Pos(bptr[i])[k]) 
		bnd.minb[k] = Pos(bptr[i])[k];
	    else if (bnd.maxb[k] < Pos(bptr[i])[k])
		bnd.maxb[k] = Pos(bptr[i])[k];
	}
    }
    *bndptr = bnd;				/* store actual bounds	    */
}

/*
 * MEDIAN_INDEX: partly sort body pointers in specified range,
 * and return index of median, using JST's median algorithm.
 */

#define SwapBody(b1,b2)  { bodyptr _tmp; _tmp = b1; b1 = b2; b2 = _tmp; }

local int median_index(bodyptr *p, int d, int lo, int hi)
{
    int i, j, m;
    real f;

    m = j = (lo + hi) / 2;
    while (lo < hi) {
	m = (lo + hi) / 2;
	f = Pos(p[m])[d];
	SwapBody(p[m], p[hi]);
	i = hi - 1;
	m = lo;
	while (Pos(p[m])[d] < f)
	     ++m;
	while (m < i) {
	    while (Pos(p[i])[d] >= f)
		if (--i == m)
		    break;
	    SwapBody(p[m], p[i]);
	    --i;
	    while (Pos(p[m])[d] < f)
		++m;
	}
	SwapBody(p[m], p[hi]);
        if (j <= m)
	    hi = m - 1;
        if (j >= m)
	    lo = m + 1;
    }
    return (m);
}

/*
 * UPWARD_PASS: Adjust bounds of each node to fit bodies exactly.
 */

local void upward_pass(kdxptr kd, int cell)
{
    kdnode *ntab = kd->ntab;

    if (ntab[cell].dim != -1) {			/* not a terminal node?	    */
	upward_pass(kd, Lower(cell));
	upward_pass(kd, Upper(cell));
	combine_nodes(&ntab[cell], &ntab[Lower(cell)], &ntab[Upper(cell)]);
    } else					/* scan bodies in node	    */
	set_bounds(&ntab[cell].bnd, kd->bptr,
		   ntab[cell].first, ntab[cell].last);
}

/*
 * COMBINE_NODES: combine bounds of two nodes.
 */

local void combine_nodes(kdnode *pout, kdnode *p1, kdnode *p2)
{
    int k;

    for (k = 0; k < 3; ++k) {
	pout->bnd.minb[k] = MIN(p2->bnd.minb[k], p1->bnd.minb[k]);
	pout->bnd.maxb[k] = MAX(p2->bnd.maxb[k], p1->bnd.maxb[k]);
    }
}
