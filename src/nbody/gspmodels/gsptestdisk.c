/*
 * gsptestdisk.c: set up a test-particle disk embedded in a gsp spheroid.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"
#include "gsp.h"

string defv[] = {		";Make test disk in a gsp spheroid",
  "gsp=???",			";Input gsp for mass profile",
  "out=???",			";Output N-body model of disk",
  "alpha=12.0",			";Inverse exponential scale length",
  "rcut=1.0",			";Outer disk cutoff radius",
  "model=0",			";Select model for disk distribution.",
				";model=-3:  quad. surface density,",
				";model=-2:  linear surface density,",
				";model=-1:  const. surface density,",
				";model=0:   normal exponential disk,",
				";model=1,2: biased exponential disk.",
  "ndisk=12288",		";Number of disk particles",
  "randspin=false",		";If true, generate test-particle ball",
  "seed=54321",			";Seed for random number generator",
  "VERSION=1.1",		";Josh Barnes  28 August 2014",
  NULL,
};

// Function prototypes.

void readgsp(void);
void writemodel(void);
void setprof(void);
void makedisk(bool randspin);
void xmatrix(matrix rmat, real theta);
void zmatrix(matrix rmat, real theta);

// Global parameters.

real alpha;				// inverse exponential scale length
real rcut;				// cut-off radius for disk model
int model;				// selects model to generate

// Global tables and data structures.

#define NTAB  (256 + 1)

real mdtab[NTAB];			// use disk mass as indp var
real rdtab[4*NTAB];			// radius as fcn of mass
real vctab[4*NTAB];			// circ. velocity as fcn of radius

gsprof *spheroid;			// spheroid mass as fcn of radius

string bodyfields[] = { PosTag, VelTag, AuxVecTag, NULL };

bodyptr disk = NULL;			// array of disk particles
int ndisk;				// number of particles in the disk

int main(int argc, string argv[])
{
  initparam(argv, defv);
  alpha = getdparam("alpha");
  rcut = getdparam("rcut");
  model = (alpha <= 0.0 ? -1 : getiparam("model"));
  ndisk = getiparam("ndisk");
  init_random(getiparam("seed"));
  readgsp();
  setprof();
  layout_body(bodyfields, Precision, NDIM);
  disk = (bodyptr) allocate(ndisk * SizeofBody);
  makedisk(getbparam("randspin"));
  if (! getbparam("randspin"))			// if spins not random
    bodyfields[2] = NULL;			// don't write AuxVec field
  writemodel();
  return (0);
}

//  readgsp: read spheroid gsp from input file.
//  ___________________________________________

void readgsp(void)
{
  stream istr;

  istr = stropen(getparam("gsp"), "r");
  get_history(istr);
  spheroid = get_gsprof(istr);
  strclose(istr);
}

//  writemodel: write N-body model to output file.
//  ______________________________________________

void writemodel(void)
{
  stream ostr;
  real tsnap = 0.0;

  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  put_snap(ostr, &disk, &ndisk, &tsnap, bodyfields);
  strclose(ostr);
}

//  setprof: initialize disk tables for radius and circular velocity.
//  _________________________________________________________________

void setprof(void)
{
  int j;
  real r, x;

  rdtab[0] = mdtab[0] = vctab[0] = 0.0;
  for (j = 1; j < NTAB; j++) {
    r = rcut * rpow(((real) j) / (NTAB - 1), 2.0);
    rdtab[j] = r;
    x = alpha * r;
    switch (model) {
      case -3:
        mdtab[j] = rsqr(rsqr(r)) / rsqr(rsqr(rcut));
	break;
      case -2:
        mdtab[j] = rqbe(r) / rqbe(rcut);
	break;
      case -1:
        mdtab[j] = rsqr(r) / rsqr(rcut);
	break;
      case 0:
	mdtab[j] = 1 - rexp(-x) - x * rexp(-x);
	break;
      case 1:
	mdtab[j] = (2 - 2 * rexp(-x) - (2*x + x*x) * rexp(-x)) / 2;
	break;
      case 2:
	mdtab[j] = (6 - 6 * rexp(-x) - (6*x + 3*x*x + x*x*x) * rexp(-x)) / 6;
	break;
      default:
	error("%s: bad choice for model\n", getargv0());
    }
    vctab[j] = rsqrt(mass_gsp(spheroid, r) / r);
  }
  if (model != -1)
    eprintf("[%s: rcut = %8.4f/alpha  M(rcut) = %8.6f mdisk]\n",
	    getargv0(), rdtab[NTAB-1] * alpha, mdtab[NTAB-1]);
  if ((mdtab[0] == mdtab[1]) || (mdtab[NTAB-2] == mdtab[NTAB-1]))
    error("%s: disk mass table is degenerate\n", getargv0());
  spline(&rdtab[NTAB], &mdtab[0], &rdtab[0], NTAB);	// for r_d = r_d(m)
  spline(&vctab[NTAB], &rdtab[0], &vctab[0], NTAB);	// for v_c = v_c(r)
}

//  makedisk: create realization of disk.
//  _____________________________________

void makedisk(bool randspin)
{
  int i;
  bodyptr bp;
  real m, r, v, phi;
  matrix xmat, zmat;
  vector tmpv;

  for (i = 0; i < ndisk; i++) {			// loop initializing bodies
    bp = NthBody(disk, i);			// set ptr to body number i
    m = mdtab[NTAB-1] * ((real) i + 0.5) / ndisk;
    r = seval(m, &mdtab[0], &rdtab[0], &rdtab[NTAB], NTAB);
    v = seval(r, &rdtab[0], &vctab[0], &vctab[NTAB], NTAB);
    phi = xrandom(0.0, TWO_PI);
    Pos(bp)[0] = r * rsin(phi);
    Pos(bp)[1] = r * rcos(phi);
    Pos(bp)[2] = 0.0;
    Vel(bp)[0] = v * rcos(phi);
    Vel(bp)[1] = - v * rsin(phi);
    Vel(bp)[2] = 0.0;
    pickshell(AuxVec(bp), NDIM, 1.0);
    if (randspin) {
      xmatrix(xmat, racos(AuxVec(bp)[2]));
      zmatrix(zmat, ratan2(AuxVec(bp)[0], AuxVec(bp)[1]));
      MULMV(tmpv, xmat, Pos(bp));
      MULMV(Pos(bp), zmat, tmpv);
      MULMV(tmpv, xmat, Vel(bp));
      MULMV(Vel(bp), zmat, tmpv);
    }
  }
}

void xmatrix(matrix rmat, real theta)
{
  real s = rsin(theta), c = rcos(theta);

  rmat[0][0] = 1.0;    rmat[0][1] = 0.0;    rmat[0][2] = 0.0;
  rmat[1][0] = 0.0;    rmat[1][1] =  c ;    rmat[1][2] =  s ;
  rmat[2][0] = 0.0;    rmat[2][1] = -s ;    rmat[2][2] =  c ;
}

void zmatrix(matrix rmat, real theta)
{
  real s = rsin(theta), c = rcos(theta);

  rmat[0][0] =  c ;    rmat[0][1] =  s ;    rmat[0][2] = 0.0;
  rmat[1][0] = -s ;    rmat[1][1] =  c ;    rmat[1][2] = 0.0;
  rmat[2][0] = 0.0;    rmat[2][1] = 0.0;    rmat[2][2] = 1.0;
}
