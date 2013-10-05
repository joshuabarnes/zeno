/*
 * DIFFSTEP.C: integrate a set of 1st order differential equations.
 */

#include "stdinc.h"
#include "mathfns.h"

/*
 * DIFFSTEP: take a double step using Runge-Kutta integrator,
 * while checking errors by comparison with a single step.
 */

void rkstep(real *, real *, int, void (*)(real *, real *), real);

real diffstep(real *x1, real *x0, int n, void (*f)(real *, real *), real h)
{
    real *y1, *y2, toterr;
    int i;

    y1 = (real *) allocate(2 * n * sizeof(real));
    y2 = y1 + n;
    rkstep(y1, x0, n, f, h);
    rkstep(y2, x0, n, f, h / 2);
    rkstep(y2, y2, n, f, h / 2);
    toterr = 0.0;
    for (i = 0; i < n; i++)
	if (y2[i] - x0[i] != 0.0)
	    toterr += rabs((y2[i] - y1[i]) / (y2[i] - x0[i]));
    for (i = 0; i < n; i++)
	x1[i] = y2[i];
    free((void *) y1);
    return (toterr);
}

/*
 * RKSTEP: take one step using Runge-Kutta 4th order integrator.
 */

void rkstep(real *y, real *x, int n, void (*f)(real *, real *), real h)
{
    real *x0, *k1, *k2, *k3, *k4;
    int i;

    x0 = (real *) allocate(5 * n * sizeof(real));
    (*f)(k1 = x0 + n, x);
    for (i = 0; i < n; i++)
	x0[i] = x[i] + h * k1[i] / 2;
    (*f)(k2 = k1 + n, x0);
    for (i = 0; i < n; i++)
	x0[i] = x[i] + h * k2[i] / 2;
    (*f)(k3 = k2 + n, x0);
    for (i = 0; i < n; i++)
	x0[i] = x[i] + h * k3[i];
    (*f)(k4 = k3 + n, x0);
    for (i = 0; i < n; i++)
	y[i] = x[i] + h * ((k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6);
    free((void *) x0);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
    "nstep=16",
    NULL,
};

void func(real *, real *);

main(int argc, string argv[])
{
    int nstep, i;
    real x[3], toterr;

    initparam(argv, defv);
    nstep = getiparam("nstep");
    x[0] = 0.0; x[1] = 0.0; x[2] = 1.0;
    printf("%8s  %10s  %10s  %10s\n", "t", "x", "y", "error");
    printf("%8.6f  %10.7f  %10.7f\n", x[0], x[1], x[2]);
    for (i = 0; i < nstep; i++) {
	toterr = diffstep(x, x, 3, func, PI / nstep);
	printf("%8.6f  %10.7f  %10.7f  %10.7f\n", x[0], x[1], x[2], toterr);
    }
    exit(0);
}

void func(real *dx, real *x)
{
    dx[0] = 1.0;
    dx[1] = x[2];
    dx[2] = - x[1];
}

#endif
