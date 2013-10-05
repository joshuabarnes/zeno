/*
 * EXPTEST.C: test integration of fields of an exponential disk.
 */

#include <stdinc.h>
#include <mathfns.h>
#include <getparam.h>
#include <filestruct.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "r=1.0",			/* radius of field point */
    "z=1.0",			/* height of field point */
    "mdisk=1.0",		/* total disk mass */
    "rdisk=1.0",		/* disk scale length */
    "step=0.01",		/* integration step */
    "kmax=10.0",		/* upper limit of k integeral */
    "VERSION=1",
    NULL,
};

real mass;
real alpha;
real phi(real s0, real z0, real step, real kmax);
real a_s(real s0, real z0, real step, real kmax);
real a_z(real s0, real z0, real step, real kmax);

int main(int argc, string *argv)
{
    real s, z, step, kmax;

    initparam(argv, defv);
    s = getdparam("r");
    z = getdparam("z");
    mass = getdparam("mdisk");
    alpha = 1.0 / getdparam("rdisk");
    step = getdparam("step");
    kmax = getdparam("kmax");
    printf("a.rad: %12.6f\n", a_s(s, z, step, kmax));
    printf("a.ver: %12.6f\n", a_z(s, z, step, kmax));
    printf("phi:   %12.6f\n", phi(s, z, step, kmax));
    return (0);
}

real s, z;		/* point at which to evaluate field */

real simpson(real (*func)(real x), real xlow, real xhigh, real step);
real d_phi(real k);
real d_a_s(real k);
real d_a_z(real k);

/*
 * PHI, A_S, A_Z: compute potential, acceleration.
 */

real phi(real s0, real z0, real step, real kmax)
{
    s = s0;
    z = z0;
    return (- mass * rqbe(alpha) * simpson(d_phi, 0.0, kmax, step));
}

real a_s(real s0, real z0, real step, real kmax)
{
    s = s0;
    z = z0;
    return (- mass * rqbe(alpha) * simpson(d_a_s, 0.0, kmax, step));
}

real a_z(real s0, real z0, real step, real kmax)
{
    s = s0;
    z = z0;
    return (- mass * rqbe(alpha) * simpson(d_a_z, 0.0, kmax, step));
}

/*
 * D_PHI, D_A_S, D_A_Z: integrands for potential, accelerations.
 */

real d_phi(real k)
{
    real psq;

    psq = alpha*alpha + k*k;
    return (rexp(- k*z) * j0(k*s) / (psq * rsqrt(psq)));
}

real d_a_s(real k)
{
    real psq;

    psq = alpha*alpha + k*k;
    return (rexp(- k*z) * j1(k*s) * k / (psq * rsqrt(psq)));
}

real d_a_z(real k)
{
    real psq;

    psq = alpha*alpha + k*k;
    return (rexp(- k*z) * j0(k*s) * k / (psq * rsqrt(psq)));
}

/*
 * SIMPSON: integrate given function using Simpson's rule.
 */

real simpson(real (*func)(real x), real xlow, real xhigh, real step)
{
    int nstep, i;
    real step1, x, v;
    double v1, v2, v4;

    nstep = 1 + 2 * (int) rceil(0.5 * (xhigh - xlow) / step);
    step1 = (xhigh - xlow) / (nstep - 1);
    v1 = v2 = v4 = 0.0;
    for (i = 1; i <= nstep; i++) {
	x = xlow + (i - 1) * step1;
	v = (*func)(x);
	if (i == 1 || i == nstep)
	    v1 = v1 + v;
	else if (i % 2 == 0)
	    v4 = v4 + v;
	else
	    v2 = v2 + v;
    }
    return (step1 * (v1 + 4.0*v4 + 2.0*v2) / 3.0);
}
