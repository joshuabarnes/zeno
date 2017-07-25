#ifndef _kdtree_h
#define _kdtree_h

/*
 * BOUND: structure representing the bounds of a kd node.
 */

typedef struct {
    real minb[3];
    real maxb[3];
} bound;

/*
 * KDNODE: structure representing a node of a kd-tree.
 */

typedef struct {
    bound bnd;				/* min, max bounds of kd-cell	    */
    int dim;				/* splitting dimension		    */
    real split;				/* median coordinate value	    */
    int first;				/* index of first body in node	    */
    int last;				/* index of last body in node	    */
} kdnode;

/*
 * KDROOT, ..., Sibling: macros for indicies into kd tree.
 */

#define KDROOT		1
#define Lower(i)	((i)<<1)
#define Upper(i)	(((i)<<1)+1)
#define Parent(i)	((i)>>1)
#define Sibling(i) 	(((i)&1)?(i)-1:(i)+1)

#define SetNext(i) {				\
    while (i&1)					\
        i=i>>1;					\
    ++i;					\
}

/*
 * KDCONTEXT, KDXPTR: structure to represent context of KD tree.
 */

typedef struct {
    int ngas;
    bodyptr *bptr;
    bound bnd;
    int nnode;
    int nsplit;
    kdnode *ntab;
} kdcontext, *kdxptr;

kdxptr init_kdtree(bodyptr, int, int);
void build_kdtree(kdxptr, int);
void finish_kdtree(kdxptr);

#endif  /* ! _kdtree_h */
