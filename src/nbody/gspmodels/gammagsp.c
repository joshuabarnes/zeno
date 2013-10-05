/*
 * GAMMAGSP.C: generate profile tables for Dehnen models.
 * See Dehnen, W. 1993, MNRAS, 265, 250.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "gsp.h"

/*
 * GAMMAGSP: initialize tables for Gamma model.
 */

gsprof *gammagsp(real gam, real mtot, real a, int np, real rmin, real rmax)
{
  gsprof *gsp;
  real *rtab, *dtab, *mtab, lgrs;
  int i;

  gsp = (gsprof *) allocate(sizeof(gsprof));
  rtab = (real *) allocate(np * sizeof(real));
  dtab = (real *) allocate(np * sizeof(real));
  mtab = (real *) allocate(np * sizeof(real));
  lgrs = rlog2(rmax / rmin) / (np - 1);
  for (i = 0; i < np; i++) {
    rtab[i] = rmin * rexp2(lgrs * i);
    dtab[i] = ((3 - gam) / FOUR_PI) * mtot * a /
                (rpow(rtab[i], gam) * rpow(rtab[i] + a, 4 - gam));
    mtab[i] = mtot * rpow(rtab[i] / (rtab[i] + a), 3 - gam);
  }
  gsp->npoint = np;
  gsp->radius = rtab;
  gsp->density = dtab;
  gsp->mass = mtab;
  /*
  gsp->alpha = - (gam + (4 - gam) * rmin / (rmin + a));
  gsp->beta = - (gam + (4 - gam) * rmax / (rmax + a));
  */
  gsp->alpha = - gam;
  gsp->beta = -4.0;
  gsp->mtot = mtot;
  return (gsp);
}

#ifdef UTILITY

#include "getparam.h"

string defv[] = {		";Generate profile for gamma model",
    "out=???",			";Output file for profile tables",
    "gamma=1.0",		";Structure parameter: 0 < gamma < 3",
    "mtot=1.0",			";Total mass of model",
    "a=1.0",			";Length scale of model",
    "npoint=257",		";Number of points in tables",
    "rrange=1/4096:16",		";Range of radii tabulated",
    "VERSION=1.2",		";Josh Barnes  4 February 1995",
    NULL,
};

int main(int argc, string argv[])
{
  real rrange[2];
  gsprof *gsp;
  stream ostr;

  initparam(argv, defv);
  setrange(rrange, getparam("rrange"));
  gsp = gammagsp(getdparam("gamma"), getdparam("mtot"), getdparam("a"),
		 getiparam("npoint"), rrange[0], rrange[1]);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  put_gsprof(ostr, gsp);
  strclose(ostr);
  return (0);
}

#endif

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Test representation of gamma model",
    "gamma=1.0",		";Structure parameter: 0 < gamma < 3",
    "mtot=1.0",			";Total mass of model",
    "a=1.0",			";Length scale of model",
    "npoint=257",		";Number of points in tables",
    "rrange=1/256:256",		";Range of radii tabulated",
    "VERSION=1.0",		";Josh Barnes  15 July 2011",
    NULL,
};

int main(int argc, string argv[])
{
  real gam, mtot, a, rrange[2], lgrs, r, d, m;
  gsprof *gsp;
  int i;

  initparam(argv, defv);
  gam = getdparam("gamma");
  mtot = getdparam("mtot");
  a = getdparam("a");
  setrange(rrange, getparam("rrange"));
  gsp = gammagsp(gam, mtot, a, getiparam("npoint"), rrange[0], rrange[1]);
  lgrs = rlog2(1000.0 / 0.001) / (6144 - 1);
  for (i = 0; i < 6144; i++) {
    r = 0.001 * rexp2(lgrs * i);
    d = ((3-gam) / FOUR_PI) * mtot * a / (rpow(r, gam) * rpow(r + a, 4-gam));
    m = mtot * rpow(r / (r + a), 3 - gam);
    printf("%12.6f  %12.5e  %12.5e  %12.5e  %12.5e\n", r, d, m,
	   (rho_gsp(gsp, r) - d) / d, (mass_gsp(gsp, r) - m) / m);

  }
  return (0);
}

#endif
