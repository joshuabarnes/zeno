/*
 * HALOFSP.C: generate profile tables for Navarro et al model.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "assert.h"
#include "fsp.h"

/*
 * HALOFSP_E: initialize Navarro et al model with exponential taper.
 */

fsprof *halofsp_e(real m_a, real a, real b, int np, real rmin, real rmax)
{
    fsprof *fsp;
    real *rtab, *dtab, *mtab, lgrs;
    real mu, rho_b, gam, m_b, r_i;
    int i;

    assert(m_a > 0.0 && a > 0.0 && a < b);
    fsp = (fsprof *) allocate(sizeof(fsprof));
    rtab = (real *) allocate(np * sizeof(real));
    dtab = (real *) allocate(np * sizeof(real));
    mtab = (real *) allocate(np * sizeof(real));
    lgrs = rlog2(rmax / rmin) / (np - 1);
    mu = m_a / (rlog(2.0) - 0.5);
    rho_b = (mu / FOUR_PI) / (b * rsqr(a + b));
    gam = b / (b + a) - 0.5;
    m_b = mu * (rlog((b + a) / a) - b / (a + b));
    eprintf("[mu = %f  rho_b = %f]\n", mu, rho_b);
    eprintf("[gamma = %f  m_b = %f]\n", gam, m_b);
    for (i = 0; i < np; i++) {
	r_i = rtab[i] = rmin * rexp2(lgrs * i);
	if (r_i <= b) {
	    dtab[i] = (mu / FOUR_PI) / (r_i * rsqr(a + r_i));
	    mtab[i] = mu * (log1p(r_i / a) - r_i / (a + r_i));
	} else {
	    dtab[i] = rho_b * rsqr(b/r_i) *
	                rexp(-2 * gam * (r_i/b - 1));
	    mtab[i] = m_b + (TWO_PI / gam) * rqbe(b) * rho_b *
	                (1 - rexp(-2 * gam * (r_i/b - 1)));
	}
    }
    fsp->npoint = np;
    fsp->radius = rtab;
    fsp->density = dtab;
    fsp->mass = mtab;
    fsp->alpha = -1.0;
    fsp->beta = -2 * gam * rmax / b - 2;
    fsp->mtot = m_b + (TWO_PI / gam) * rqbe(b) * rho_b;
    eprintf("[beta = %f  mtot = %f]\n", fsp->beta, fsp->mtot);
    return (fsp);
}

extern double erf(double);

/*
 * HALOFSP_G: initialize Navarro et al model with gaussian taper.
 */

fsprof *halofsp_g(real m_a, real a, real b, int np, real rmin, real rmax)
{
    fsprof *fsp;
    real *rtab, *dtab, *mtab, lgrs;
    real mu, rho_b, gam, m_b, r_i, ghalf, pi3half;
    int i;

    assert(m_a > 0.0 && a > 0.0 && a < b);
    fsp = (fsprof *) allocate(sizeof(fsprof));
    rtab = (real *) allocate(np * sizeof(real));
    dtab = (real *) allocate(np * sizeof(real));
    mtab = (real *) allocate(np * sizeof(real));
    lgrs = rlog2(rmax / rmin) / (np - 1);
    mu = m_a / (rlog(2.0) - 0.5);
    rho_b = (mu / FOUR_PI) / (b * rsqr(a + b));
    gam = b / (b + a) - 0.5;
    m_b = mu * (rlog((b + a) / a) - b / (a + b));
    eprintf("[mu = %f  rho_b = %f]\n", mu, rho_b);
    eprintf("[gamma = %f  m_b = %f]\n", gam, m_b);
    ghalf = rsqrt(gam);
    pi3half = rsqrt(rqbe(PI));
    for (i = 0; i < np; i++) {
	r_i = rtab[i] = rmin * rexp2(lgrs * i);
	if (r_i <= b) {
	    dtab[i] = (mu / FOUR_PI) / (r_i * rsqr(a + r_i));
	    mtab[i] = mu * (log1p(r_i / a) - r_i / (a + r_i));
	} else {
	    dtab[i] = rho_b * rsqr(b/r_i) *
	                rexp(- gam * (rsqr(r_i/b) - 1));
	    mtab[i] = m_b + (2 * pi3half / ghalf) * rqbe(b) * rho_b *
	                rexp(gam) * (erf(ghalf * r_i/b) - erf(ghalf));
	}
    }
    fsp->npoint = np;
    fsp->radius = rtab;
    fsp->density = dtab;
    fsp->mass = mtab;
    fsp->alpha = -1.0;
    fsp->beta = -2 * gam * rsqr(rmax / b) - 2;
    fsp->mtot = m_b + (2 * pi3half / ghalf) * rqbe(b) * rho_b *
	          rexp(gam) * (1.0 - erf(ghalf));
    eprintf("[beta = %f  mtot = %f]\n", fsp->beta, fsp->mtot);
    return (fsp);
}

extern float gammp(), gammln();

/*
 * HALOFSP_SW: initialize Navarro et al model with fast
 * exponential taper (Springel & White 1998, astro-ph/9807320).
 */

fsprof *halofsp_sw(real m_a, real a, real b, int np, real rmin, real rmax)
{
    fsprof *fsp;
    real *rtab, *dtab, *mtab, lgrs;
    real mu, rho_b, gam, m_b, gscale, r_i;
    int i;

    assert(m_a > 0.0 && a > 0.0 && a < b);
    fsp = (fsprof *) allocate(sizeof(fsprof));
    rtab = (real *) allocate(np * sizeof(real));
    dtab = (real *) allocate(np * sizeof(real));
    mtab = (real *) allocate(np * sizeof(real));
    lgrs = rlog2(rmax / rmin) / (np - 1);
    mu = m_a / (rlog(2.0) - 0.5);
    rho_b = (mu / FOUR_PI) / (b * rsqr(a + b));
    gam = (b / a) - (a + 3 * b) / (a + b);
    m_b = mu * (rlog((b + a) / a) - b / (a + b));
    gscale = rexp(gam * rlog(a/b) + b/a + gammln(3 + gam));
    eprintf("[mu = %f  rho_b = %f]\n", mu, rho_b);
    eprintf("[gamma = %f  m_b = %f]\n", gam, m_b, gscale);
    for (i = 0; i < np; i++) {
	rtab[i] = r_i = rmin * rexp2(lgrs * i);
	if (r_i <= b) {
	    dtab[i] = (mu / FOUR_PI) / (r_i * rsqr(a + r_i));
	    mtab[i] = mu * (log1p(r_i / a) - r_i / (a + r_i));
	} else {
	    dtab[i] = rho_b * rpow(r_i / b, gam) * rexp(- (r_i - b) / a);
	    mtab[i] = m_b + FOUR_PI * rho_b * rqbe(a) * gscale *
			  (gammp(3 + gam, r_i/a) - gammp(3 + gam, b/a));
	}
     }
    fsp->npoint = np;
    fsp->radius = rtab;
    fsp->density = dtab;
    fsp->mass = mtab;
    fsp->alpha = -1.0;
    fsp->beta = gam - rmax / a;
    fsp->mtot = m_b + FOUR_PI * rho_b * rqbe(a) * gscale *
		    (1.0 - gammp(3 + gam, b/a));
    eprintf("[beta = %f  mtot = %f]\n", fsp->beta, fsp->mtot);
    return (fsp);
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Generate profile for halo model",
    "out=???",			";Output file for profile tables",
    "m_a=1.0",			";Mass within radius a",
    "a=1.0",			";Radial scale of model",
    "b=4.0",			";Radius to begin taper",
    "taper=exp",		";Tapering: exp, gauss, or sw",
    "npoint=257",		";Number of points in tables",
    "rrange=1/4096:16",		";Range of radii tabulated",
    "VERSION=1.3",		";Josh Barnes  10 April 2002",
    NULL,
};

int main(int argc, string argv[])
{
    real rrange[2];
    fsprof *fsp;
    stream ostr;

    initparam(argv, defv);
    setrange(rrange, getparam("rrange"));
    if (streq(getparam("taper"), "exp"))
	fsp = halofsp_e(getdparam("m_a"), getdparam("a"), getdparam("b"),
			getiparam("npoint"), rrange[0], rrange[1]);
    else if (streq(getparam("taper"), "gauss"))
	fsp = halofsp_g(getdparam("m_a"), getdparam("a"), getdparam("b"),
			getiparam("npoint"), rrange[0], rrange[1]);
    else if (streq(getparam("taper"), "sw"))
	fsp = halofsp_sw(getdparam("m_a"), getdparam("a"), getdparam("b"),
			 getiparam("npoint"), rrange[0], rrange[1]);
    else
	error("%s: unknown taper option %s\n", getargv0(), getparam("taper"));
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_fsprof(ostr, fsp);
    strclose(ostr);
    return (0);
}

#endif
