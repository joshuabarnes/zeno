/*
 * gammagsp.c: generate profile tables for gamma models.
 * See Dehnen, W. 1993, MNRAS, 265, 250.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

//  Functions defining gamma models.  First argument holds parameters.

local double density(void *p, double r);
local double mass(void *p, double r);
local double potential(void *p, double r);
local double radius(void *p, double m);

//  pars: structure used to pass parameters to model functions.
//  ___________________________________________________________

typedef struct {
  double gam;
  double mtot;
  double ascl;
} pars;

//  gsp_gamma: initialize tables for gamma model.
//  _____________________________________________

gsprof *gsp_gamma(double gam, double mtot, double ascl,
		  int np, double rmin, double rmax)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  double lgrs;
  pars p = { gam, mtot, ascl };

  gsp->npoint = np;
  gsp->radius  = (double *) allocate(np * sizeof(double));
  gsp->density = (double *) allocate(np * sizeof(double));
  gsp->mass    = (double *) allocate(np * sizeof(double));
  // gsp->phi     = (double *) allocate(np * sizeof(double));
  lgrs = log2(rmax / rmin) / (np - 1);
  eprintf("[%s.gsp_gamma: lgrs = %f (%a)]\n", getprog(), lgrs, lgrs);
  for (int i = 0; i < np; i++) {
    gsp->radius[i]  = rmin * exp2(lgrs * i);
    gsp->density[i] = density(&p, gsp->radius[i]);
    gsp->mass[i]    = mass(&p, gsp->radius[i]);
    if (i > 0 && gsp->mass[i] == gsp->mass[i-1])
      error("%s.gsp_gamma: mass degenerate (i = %d)\n", getprog(), i);
    // gsp->phi[i]     = potential(&p, gsp->radius[i]);
  }
  // gsp->alpha = - (gam + (4 - gam) * rmin / (rmin + ascl));
  // gsp->beta = - (gam + (4 - gam) * rmax / (rmax + ascl));
  gsp->alpha = - gam;				// true asymptotic slopes
  gsp->beta = -4.0;				// yield higher accuracy
  gsp->mtot = mtot;
  gsp_test_rad(gsp, gsp_rho, density, &p, "density");
  gsp_test_rad(gsp, gsp_mass, mass, &p, "mass");
  gsp_test_mass(gsp, gsp_mass_rad, radius, &p, "radius");
  gsp_test_rad(gsp, gsp_phi, potential, &p, "potential");
  return gsp;
}

//  Access macros for parameters.

#define _gam   (((pars *) p)->gam)
#define _mtot  (((pars *) p)->mtot)
#define _ascl  (((pars *) p)->ascl)

//  density: compute space density as function of radius.
//  _____________________________________________________

local double density(void *p, double r)
{
  return (((3 - _gam) / (4 * M_PI)) * _mtot * _ascl /
	  (rpow(r, _gam) * pow(r + _ascl, 4 - _gam)));
}

//  mass: compute enclosed mass as function of radius.
//  __________________________________________________

local double mass(void *p, double r)
{
  return (_mtot * pow(r / (r + _ascl), 3 - _gam));
}

//  potential: compute potential as function of radius.
//  ___________________________________________________

local double potential(void *p, double r)
{
  if (_gam != 2.0)
    return ((_mtot/_ascl) * (1 - pow(r / (r + _ascl), 2-_gam)) / (_gam-2));
  else
    return ((_mtot/_ascl) * log(r / (r + _ascl)));
}

//  radius: compute radius as function of enclosed mass.
//  ____________________________________________________

local double radius(void *p, double m)
{
  double xi = pow(m / _mtot, 1.0 / (3 - _gam));
  return (_ascl * xi / (1 - xi));
}

#ifdef UTILITY

string defv[] = {		";Generate profile for gamma model",
  "out=",			";Output file for profile tables",
  "gamma=1.0",			";Structure parameter: 0 < gamma < 3",
  "mtot=1.0",			";Total mass of model",
  "a=1.0",			";Length scale of model",
  "npoint=1281",		";Number of points in tables",
  "rrange=1/1024:1024",		";Range of radii tabulated",
  "VERSION=2.0",		";Josh Barnes  28 May 2017",
  NULL,
};

int main(int argc, string argv[])
{
  real rrange[2];
  gsprof *gsp;
  stream ostr;

  initparam(argv, defv);
  setrange(rrange, getparam("rrange"));
  gsp = gsp_gamma(getdparam("gamma"), getdparam("mtot"), getdparam("a"),
		  getiparam("npoint"), rrange[0], rrange[1]);
  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    gsp_write(ostr, gsp);
  }
  fflush(NULL);
  return 0;
}

#endif
