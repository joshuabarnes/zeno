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
#include <gsl/gsl_math.h>

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
  "VERSION=2.0",		";Josh Barnes  22 June 2017",
  NULL,
};

// Function prototypes.

void readgsp(void);
void writemodel(void);
void setprof(int model, double alpha, double rcut);
void makedisk(bool randspin);
void xmatrix(matrix rmat, double theta);
void zmatrix(matrix rmat, double theta);

// Global arrays and data structures.

#define NTAB  (256 + 1)

double mdtab[NTAB];			// use disk mass as indp var
double rdtab[NTAB];			// radius as fcn of mass
double vctab[NTAB];			// circ. velocity as fcn of radius
gsl_interp *rm_spline;			// interpolator for r(m)
gsl_interp *vr_spline;			// interpolator for v(r)

gsprof *spheroid;			// spheroid mass as fcn of radius

string bodyfields[] = { PosTag, VelTag, AuxVecTag, NULL };

bodyptr disk = NULL;			// array of disk particles
int ndisk;				// number of particles in the disk

int main(int argc, string argv[])
{
  double alpha, rcut;
  int model;

  initparam(argv, defv);
  alpha = getdparam("alpha");
  rcut = getdparam("rcut");
  model = (alpha <= 0.0 ? -1 : getiparam("model"));
  ndisk = getiparam("ndisk");
  init_random(getiparam("seed"));
  readgsp();
  setprof(model, alpha, rcut);
  layout_body(bodyfields, Precision, NDIM);
  disk = (bodyptr) allocate(ndisk * SizeofBody);
  makedisk(getbparam("randspin"));
  if (! getbparam("randspin"))			// if spins not random
    bodyfields[2] = NULL;			// don't write AuxVec field
  writemodel();
  fflush(NULL);
  return 0;
}

//  readgsp: read spheroid gsp from input file.
//  ___________________________________________

void readgsp(void)
{
  stream istr;

  istr = stropen(getparam("gsp"), "r");
  get_history(istr);
  spheroid = gsp_read(istr);
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
}

//  setprof: initialize disk tables for radius and circular velocity.
//  _________________________________________________________________

void setprof(int model, double alpha, double rcut)
{
  int j;
  double r, x;

  rdtab[0] = mdtab[0] = vctab[0] = 0.0;
  for (j = 1; j < NTAB; j++) {
    r = rcut * pow(((double) j) / (NTAB - 1), 2.0);
    rdtab[j] = r;
    x = alpha * r;
    switch (model) {
      case -3:
        mdtab[j] = gsl_pow_4(r / rcut);
	break;
      case -2:
        mdtab[j] = gsl_pow_3(r / rcut);
	break;
      case -1:
        mdtab[j] = gsl_pow_2(r / rcut);
	break;
      case 0:
	mdtab[j] = 1 - exp(-x) - x * exp(-x);
	break;
      case 1:
	mdtab[j] = (2 - 2 * exp(-x) - (2*x + x*x) * exp(-x)) / 2;
	break;
      case 2:
	mdtab[j] = (6 - 6 * exp(-x) - (6*x + 3*x*x + x*x*x) * exp(-x)) / 6;
	break;
      default:
	error("%s: bad choice for model\n", getprog());
    }
    vctab[j] = sqrt(gsp_mass(spheroid, r) / r);
  }
  if (model > -1)
    eprintf("[%s: rcut = %8.4f/alpha  M(rcut) = %8.6f*mdisk]\n",
	    getprog(), rdtab[NTAB-1] * alpha, mdtab[NTAB-1]);
  if ((mdtab[0] == mdtab[1]) || (mdtab[NTAB-2] == mdtab[NTAB-1]))
    error("%s: disk mass table is degenerate\n", getprog());
  rm_spline = gsl_interp_alloc(gsl_interp_akima, NTAB);
  gsl_interp_init(rm_spline, mdtab, rdtab, NTAB);
  vr_spline = gsl_interp_alloc(gsl_interp_akima, NTAB);
  gsl_interp_init(vr_spline, rdtab, vctab, NTAB);
}

//  makedisk: create realization of disk.
//  _____________________________________

void makedisk(bool randspin)
{
  int i;
  bodyptr bp;
  double m, r, v, phi;
  matrix xmat, zmat;
  vector tmpv;

  for (i = 0; i < ndisk; i++) {			// loop initializing bodies
    bp = NthBody(disk, i);			// set ptr to body number i
    m = mdtab[NTAB-1] * ((double) i + 0.5) / ndisk;
    r = gsl_interp_eval(rm_spline, mdtab, rdtab, m, NULL);
    v = gsl_interp_eval(vr_spline, rdtab, vctab, r, NULL);
    phi = xrandom(0.0, 2 * M_PI);
    Pos(bp)[0] = r * sin(phi);
    Pos(bp)[1] = r * cos(phi);
    Pos(bp)[2] = 0.0;
    Vel(bp)[0] = v * cos(phi);
    Vel(bp)[1] = - v * sin(phi);
    Vel(bp)[2] = 0.0;
    pickshell(AuxVec(bp), NDIM, 1.0);
    if (randspin) {
      xmatrix(xmat, acos(AuxVec(bp)[2]));
      zmatrix(zmat, atan2(AuxVec(bp)[0], AuxVec(bp)[1]));
      MULMV(tmpv, xmat, Pos(bp));
      MULMV(Pos(bp), zmat, tmpv);
      MULMV(tmpv, xmat, Vel(bp));
      MULMV(Vel(bp), zmat, tmpv);
    }
  }
}

void xmatrix(matrix rmat, double theta)
{
  real s = sin(theta), c = cos(theta);

  rmat[0][0] = 1.0;    rmat[0][1] = 0.0;    rmat[0][2] = 0.0;
  rmat[1][0] = 0.0;    rmat[1][1] =  c ;    rmat[1][2] =  s ;
  rmat[2][0] = 0.0;    rmat[2][1] = -s ;    rmat[2][2] =  c ;
}

void zmatrix(matrix rmat, double theta)
{
  real s = sin(theta), c = cos(theta);

  rmat[0][0] =  c ;    rmat[0][1] =  s ;    rmat[0][2] = 0.0;
  rmat[1][0] = -s ;    rmat[1][1] =  c ;    rmat[1][2] = 0.0;
  rmat[2][0] = 0.0;    rmat[2][1] = 0.0;    rmat[2][2] = 1.0;
}
