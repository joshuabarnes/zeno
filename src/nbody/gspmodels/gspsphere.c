/*
 * gspsphere.c: make SPH realization of GSP.
 */

#include "stdinc.h"
#include "assert.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"

string defv[] = {               ";Make SPH gas sphere from GSP",
  "gsp=???",                    ";Input GSP for density profile",
  "out=",                       ";Output SnapShot with bodies",
  "grav=",                      ";Input GSP for gravity calculation",
  "gamma=5/3",                  ";Ratio of specific heats",
  "nbody=16384",                ";Number of bodies to generate",
  "mcut=0.999",                 ";Radial cutoff in terms of total mass",
  "seed=54321",                 ";Usual random number seed",
  "zerocm=false",               ";Transform to center of mass coords",
  "type=0x60",			";Specify body type for SPH calculation",
  "VERSION=2.0",                ";Josh Barnes  24 June 2017",
  NULL,
};

//  Prototypes for model construction and I/O.

void gspsphere(void);                   // construct spherical SPH model
void readgsp(void);                     // read input profile(s)
void writemodel(void);                  // write SPH model to output

//  Global data for communication between major routines.

gsprof *gsp, *ggsp;                     // profiles for mass and gravity
bodyptr btab = NULL;                    // pointer to array of bodies
int nbody;                              // number of bodies in array

string bodyfields[] = {
  PosTag, VelTag, MassTag, RhoTag, EntFuncTag, UinternTag, TypeTag, NULL
};

int main(int argc, string argv[])
{
  initparam(argv, defv);
  readgsp();
  init_random(getiparam("seed"));
  layout_body(bodyfields, Precision, NDIM);
  gspsphere();
  if (getbparam("zerocm"))
      snapcenter(btab, nbody, MassField.offset);
  writemodel();
  fflush(NULL);
  return 0;
}

//  gspsphere: construct realization from GSP data.
//  _______________________________________________

void gspsphere(void)
{
  double gam, mcut, r, sig2, eint = 0.0;
  byte type = (0x7f & getiparam("type"));

  nbody = getiparam("nbody");
  assert(nbody > 0);
  gam = getdparam("gamma");
  mcut = getdparam("mcut");
  assert(0.0 < mcut && mcut <= 1.0);
  if (btab == NULL) {
    btab = (bodyptr) allocate(nbody * SizeofBody);
    gsp_calc_sig2(gsp, ggsp, 0.0);
  }
  for (bodyptr bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    Mass(bp) = gsp->mtot / nbody;
    r = gsp_mass_rad(gsp, xrandom(0.0, mcut * gsp->mtot));
    pickshell(Pos(bp), NDIM, r);
    CLRV(Vel(bp));
    Rho(bp) = gsp_rho(gsp, r);
    sig2 = gsp_sig2(gsp, r);
    EntFunc(bp) = sig2 / pow(Rho(bp), gam - 1);
    Uintern(bp) = sig2 / (gam - 1);
    Type(bp) = type;
    eint += Mass(bp) * Uintern(bp);
  }
  eprintf("[%s: thermal energy = %f]\n", getprog(), eint);
}

void readgsp(void)
{
  stream istr;

  istr = stropen(getparam("gsp"), "r");
  get_history(istr);
  gsp = ggsp = gsp_read(istr);
  if (! strnull(getparam("grav"))) {
    istr = stropen(getparam("grav"), "r");
    get_history(istr);
    ggsp = gsp_read(istr);
  }
}

void writemodel(void)
{
  stream ostr;
  real tsnap = 0.0;
  
  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_snap(ostr, &btab, &nbody, &tsnap, bodyfields);
  }
}
