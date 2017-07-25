/*
 * plumgsp.c: generate profile tables for Plummer model.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "gsp.h"

//  Functions defining Plummer models.  First argument holds parameters.

local double density(void *p, double r);
local double gradient(void *p, double r);
local double mass(void *p, double r);
local double potential(void *p, double r);

//  pars: structure used to pass parameters.
//  ________________________________________

typedef struct {
  double mtot;
  double ascl;
} pars;

//  gsp_plum: initialize tables for Plummer model.
//  ______________________________________________

gsprof *gsp_plum(double mtot, double ascl,
		  int np, double rmin, double rmax)
{
  gsprof *gsp = (gsprof *) allocate(sizeof(gsprof));
  double lgrs;
  pars p = { mtot, ascl };

  gsp->npoint = np;
  gsp->radius  = (double *) allocate(np * sizeof(double));
  gsp->density = (double *) allocate(np * sizeof(double));
  gsp->mass    = (double *) allocate(np * sizeof(double));
  // gsp->phi     = (double *) allocate(np * sizeof(double));
  lgrs = log2(rmax / rmin) / (np - 1);
  eprintf("[%s.gsp_plum: lgrs = %f (%a)]\n", getprog(), lgrs, lgrs);
  for (int i = 0; i < np; i++) {
    gsp->radius[i]  = rmin * exp2(lgrs * i);
    gsp->density[i] = density(&p, gsp->radius[i]);
    gsp->mass[i]    = mass(&p, gsp->radius[i]);
    if (i > 0 && gsp->mass[i] == gsp->mass[i-1])
      error("%s.gsp_plum: mass degenerate (i = %d)\n", getprog(), i);
    // gsp->phi[i]     = potential(&p, gsp->radius[i]);
  }
  gsp->alpha = 0.0;
  gsp->beta = -5.0;
  gsp->mtot = mtot;
  gsp_test_rad(gsp, gsp_rho, density, &p, "density");
  gsp_test_rad(gsp, gsp_grad, gradient, &p, "gradient");
  gsp_test_rad(gsp, gsp_mass, mass, &p, "mass");
  gsp_test_rad(gsp, gsp_phi, potential, &p, "potential");
  return gsp;
}

//  Access macros for parameters.

#define _mtot  (((pars *) p)->mtot)
#define _ascl  (((pars *) p)->ascl)

//  density: compute space density as function of radius.
//  _____________________________________________________

local double density(void *p, double r)
{
  return ((3 / (4 * M_PI)) * _mtot * _ascl*_ascl /
	  pow(r*r + _ascl*_ascl, 2.5));

}

//  gradient: compute density gradient as function of radius.
//  _________________________________________________________

local double gradient(void *p, double r)
{
  return (- (15 / (4 * M_PI)) * _mtot * _ascl*_ascl * r /
	  pow(r*r + _ascl*_ascl, 3.5));
}

//  mass: compute enclosed mass as function of radius.
//  __________________________________________________

local double mass(void *p, double r)
{
  return (_mtot * r*r*r / pow(r*r + _ascl*_ascl, 1.5));
}

//  potential: compute potential as function of radius.
//  ___________________________________________________

local double potential(void *p, double r)
{
  return (- _mtot / sqrt(r*r + _ascl*_ascl));
}

#ifdef UTILITY

string defv[] = {		";Generate profile for Plummer model",
  "out=",			";Output file for profile tables",
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
  gsp = gsp_plum(getdparam("mtot"), getdparam("a"),
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
