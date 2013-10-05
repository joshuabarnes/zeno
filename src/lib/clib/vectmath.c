/*
 * vectmath.c: source code for vector/matrix operations.
 * Joshua Barnes  2 Aug 1996   Honolulu, HI.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"

real _dotvp(real *v, real *u, int n)
{
    real s;

    s = 0.0;
    while (--n >= 0)
	s += (*v++) * (*u++);
    return (s);
}

real _absv(real *v, int n)
{
    real s;

    s = 0.0;
    while (--n >= 0) {
	s += (*v) * (*v);
	v++;
    }
    return (rsqrt(s));
}

real _distv(real *v, real *u, int n)
{
    real d, s;

    s = 0.0;
    while (--n >= 0) {
	d = (*v++) - (*u++);
	s += d * d;
    }
    return (rsqrt(s));
}

real _tracem(real *p, int n)
{
    real s;
    int i;

    s = 0.0;
    for (i = n; --i >= 0; ) {
	s += (*p);
	p += (n + 1);			// next diag. element
    }
    return (s);
}
