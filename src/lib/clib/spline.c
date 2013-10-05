/*
 * spline.c: functions to evaluate cubic spline approximations.
 * Reference: Forsythe, Malcolm & Moler, "Computer Methods for
 *	      Mathematical Computations", pp. 76-79.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"

local void splsub(real *, real *, real *, real *, real *, int);

//  ___________________________________
//  spline: compute cubic spline coefs.

void spline(real *coef, real *x, real *y, int n)
{
    splsub(&coef[0], &coef[n], &coef[2*n], x, y, n);
}

local void splsub(real *b, real *c, real *d, real *x, real *y, int n)
{
    int i;
    real t1, tn, t;

    if (n < 3)
	error("%s.spline: n lt 3\n", getprog());
    d[0] = x[1] - x[0];
    c[1] = (y[1] - y[0]) / d[0];
    for (i = 1; i <= n-2; i++) {
	d[i] = x[i+1] - x[i];
	b[i] = 2.0 * (d[i-1] + d[i]);
	c[i+1] = (y[i+1] - y[i]) / d[i];
	c[i] = c[i+1] - c[i];
    }
    b[ 0 ] = - d[ 0 ];
    b[n-1] = - d[n-2];
    if (n == 3)
	c[0] = c[n-1] = 0.0;
    else {
	t1 = c[ 2 ] / (x[ 3 ] - x[ 1 ]) - c[ 1 ] / (x[ 2 ] - x[ 0 ]);
	tn = c[n-2] / (x[n-1] - x[n-3]) - c[n-3] / (x[n-2] - x[n-4]);
	c[ 0 ] =   t1 * d[ 0 ]*d[ 0 ] / (x[ 3 ] - x[ 0 ]);
	c[n-1] = - tn * d[n-2]*d[n-2] / (x[n-1] - x[n-4]);
    }
    for (i = 1; i <= n-1; i++) {
	t = d[i-1] / b[i-1];
	b[i] = b[i] - t * d[i-1];
	c[i] = c[i] - t * c[i-1];
    }
    c[n-1] = c[n-1] / b[n-1];
    for (i = n-2; i >= 0; i--)
	c[i] = (c[i] - d[i] * c[i+1]) / b[i];
    b[n-1] = (y[n-1] - y[n-2]) / d[n-2] + d[n-2] * (c[n-2] + 2 * c[n-1]);
    for (i = 0; i <= n-2; i++) {
	b[i] = (y[i+1] - y[i]) / d[i] - d[i] * (c[i+1] + 2 * c[i]);
	d[i] = (c[i+1] - c[i]) / d[i];
	c[i] = 3 * c[i];
    }
    c[n-1] = 3 * c[n-1];
    d[n-1] = d[n-2];
}

//  ___________________________________________
//  seval: evaluate cubic spline interpolation.

real seval(real x0, real *x, real *y, real *coef, int n)
{
    int i, j, k;
    real u;

    i = 0;
    k = n;
    while (i+1 < k) {
	j = (i + k) / 2;
	if (x[j] <= x0)
	    i = j;
	else
	    k = j;
    }
    u = x0 - x[i];
    return (y[i] + u * (coef[i] + u * (coef[n+i] + u * coef[2*n+i])));
}

//  ____________________________________________
//  spldif: evaluate derivative of cubic spline.

real spldif(real x0, real *x, real *y, real *coef, int n)
{
    int i, j, k;
    real u;

    i = 0;
    k = n;
    while (i+1 < k) {
	j = (i + k) / 2;
	if (x[j] <= x0)
	    i = j;
	else
	    k = j;
    }
    u = x0 - x[i];
    return (coef[i] + u * (2.0*coef[i+n] + u * 3.0*coef[i+2*n]));
}

#ifdef TESTBED

#define N       11

#define fun(x0) (rsqrt(x0) * rcos(PI*x0))
#define fpr(x0) (0.5*rcos(PI*x0) / rsqrt(x0) - rsqrt(x0) * PI*rsin(PI*x0))

void main(int argc, string argv[])
{
    int i;
    real x[N], y[N], coef[3*N];
    double x0;

    for (i = 0; i < N; i++) {
	x[i] = i / (N - 1.0);
	y[i] = fun(x[i]);
    }
    spline(coef, x, y, N);
    printf("\n%12s%12s%12s%12s%12s\n", "x", "y", "b", "c", "d");
    for (i = 0; i < N; i++) {
	printf("%12.6f%12.6f%12.6f%12.6f%12.6f\n",
	       x[i], y[i], coef[i], coef[N+i], coef[2*N+i]);
    }
    printf("\n");
    for ( ; ; ) {
	printf("x0: ");
	scanf("%lf", &x0);
	printf("        %12s\t%12s\n", "f(x)", "f'(x)");
	printf("exact   %12.6f\t%12.6f\n", fun(x0), fpr(x0));
	printf("spline  %12.6f\t%12.6f\n",
	       seval(x0, x, y, coef, N),
	       spldif(x0, x, y, coef, N));
    }
}

#endif
