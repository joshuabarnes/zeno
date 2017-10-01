/*
 * testcode.c: follow test body orbits in given potential.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "gsp.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#if defined(MPICODE)
#  include <mpi.h>
#endif

#if !(defined(NBDGRAV) || defined(GSPGRAV) || defined(HMPGRAV))
#define NBDGRAV
#endif

string defv[] = {
#if defined(NBDGRAV)
#if !defined(ZOMGRAV)
				";Follow orbits in N-body potential",
#else
				";Follow orbits in zombie potential",
#endif
#elif defined(GSPGRAV)
				";Follow orbits in GSP potential",
#elif defined(HMPGRAV)
				";Follow orbits in Hernquist model potential",
#endif
  "in=???",			";Initial conditions for test bodies",
  "out=???",			";Output stream of test body snapshots",
#if defined(NBDGRAV)
  "grav=???",			";N-body configuration generating potential",
  "xyzscale=1.0,1.0,1.0",	";N-body scale factors along X, Y, Z",
  "eps=0.025",			";Softening parameter for force calculation",
#endif
#if defined(GSPGRAV)
  "gsp=???",			";General spherical profile for potential",
#endif
#if defined(ZOMGRAV)
  "gsp=???",			";General spherical profile for zombies",
#endif
#if defined(HMPGRAV)
  "M=1.0",			";Total mass of Hernquist model",
  "a=1.0",			";Major axis of Hernquist model",
  "b=1.0",			";Minor axis of Hernquist model",
  "tol=1.0e-6",			";Error tolerance for integrals",
#endif
  "tstop=2.0",			";Duration of integration",
  "dtime=1/32",			";Leap-frog integration time-step",
  "decrit=1.0e-5",		";Initial energy change tolerance",
  "outputs=" PosTag "," VelTag,	";Per-body data arrays written to output",
  "dtout=1/4",			";Time interval between output frames",
  "confmt=",			";Save configuration at unit times",
  "VERSION=2.2",		";Josh Barnes  22 July 2017",
  NULL,
};

//  Force calculation and other routines.
//  _____________________________________

local void mainloop(real epot0);

local void forcecalc(void);

local void sumforces(bodyptr btab, int nbody, bodyptr gtab, int ngrav,
		     real xyzscale[3], real eps2);

local void gspforces(bodyptr btab, int nbody, gsprof *ggsp);

local void hmpforces(bodyptr btab, int nbody, real M, real a, real b,
		     real tol);

local void initcode(void);

local void saveconfig(string confmt);

//  Global data and parameters.  Each version uses a different subset.
//  __________________________________________________________________

int nproc = 1, iproc = 0;		// number of processes, process index

bodyptr btab = NULL;			// array of test bodies to integrate

int nbody = 0;				// number of test bodies to integrate

real tnow, tstop, dtime, decrit0, tout;	// integration parameters

stream ostr;				// output stream, always open

string *otags;				// list of body tags to output

bodyptr gtab = NULL;			// array of bodies for N-body pot.

int ngrav = 0;				// number of bodies for N-body pot.

real xyzscale[NDIM], eps;		// N-body potential parameters

gsprof *ggsp = NULL;			// GSP for general spherical pot.

real Mhmp, ahmp, bhmp, tol;		// spheroidal Hernquist model params.

//  Define field for inital binding energy.
//  _______________________________________

#define EinitPBF  phatbody[NewBodyFields+0]	// add initial E to bodies
#define Einit(b)  SelectReal(b,EinitPBF.offset)	// accessor macro for above
#define EinitTag  "Einit"			// field name for above

string bodytags[MaxBodyFields] = {
  PosTag, VelTag, MassTag, PhiTag, AccTag, EinitTag, NULL
};

int main(int argc, string argv[])
{
  real epot0;

#if defined(MPICODE)
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
#endif
  initparam(argv, defv);
  initcode();
  forcecalc();
  epot0 = 0.0;					// use as energy scale
  for (bodyptr bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    epot0 += Phi(bp) / nbody;			// compute avg. potential
    Einit(bp) = Phi(bp) + dotvp(Vel(bp), Vel(bp)) / 2;
  }
  eprintf("[%s[%d]: initial average potential = %g]\n", getprog(), iproc, epot0);
  if (iproc == 0) {
    put_snap(ostr, &btab, &nbody, &tnow, otags);
    tout = tnow + getdparam("dtout");		// set next output time
    fflush(NULL);
  }
  mainloop(epot0);
  if (ggsp != NULL)
    gsp_free(ggsp);
#if defined(MPICODE)
  MPI_Finalize();
#endif
  return 0;
}

//  mainloop: perform simulation with test bodies and optional zombies.
//  ___________________________________________________________________

local void mainloop(real epot0)
{
  real decrit, demin, demax, derms, de2avg, denow;

  decrit = decrit0;
  demin = demax = derms = 0.0;			// track maximum errors
  while (tnow < tstop) {			// enter main loop
    for (bodyptr bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      ADDMULVS(Vel(bp), Acc(bp), 0.5 * dtime);	// step velocities by dt/2
      ADDMULVS(Pos(bp), Vel(bp), dtime);	// step positions by dt
    }
#if defined(ZOMGRAV)
    for (bodyptr gp = gtab; gp < NthBody(gtab, ngrav); gp = NextBody(gp)) {
      ADDMULVS(Vel(gp), Acc(gp), 0.5 * dtime);
      ADDMULVS(Pos(gp), Vel(gp), dtime);
    }
#endif
    forcecalc();
    for (bodyptr bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
      ADDMULVS(Vel(bp), Acc(bp), 0.5 * dtime);	// step velocities by dt/2
#if defined(ZOMGRAV)
    for (bodyptr gp = gtab; gp < NthBody(gtab, ngrav); gp = NextBody(gp))
      ADDMULVS(Vel(gp), Acc(gp), 0.5 * dtime);
#endif
    tnow = tnow + dtime;			// complete time-step
    de2avg = 0.0;
    for (bodyptr bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      denow = (dotvp(Vel(bp), Vel(bp))/2 + Phi(bp) - Einit(bp)) / ABS(epot0);
      demin = MIN(demin, denow);
      demax = MAX(demax, denow);
      de2avg += rsqr(denow) / nbody;
    }
    derms = MAX(derms, rsqrt(de2avg));
    if (demin < -decrit || demax > decrit) {
      if (iproc == 0)
	eprintf("[%s: warning: energy error exceeds %.4e at time = %-12.8f\n"
		" min,max,rms = %e,%e,%e  threshold now %.4e]\n", getprog(),
		decrit, tnow, demin, demax, rsqrt(de2avg), decrit * sqrt(2.0));
      decrit = decrit * sqrt(2.0);
    }
    if (iproc == 0 && tout <= tnow) {
      put_snap(ostr, &btab, &nbody, &tnow, otags);
      tout = tout + getdparam("dtout");
    }
    if (tnow == floor(tnow) && ! strnull(getparam("confmt")))
      saveconfig(getparam("confmt"));
    fflush(NULL);
  }
  if (iproc == 0)
    eprintf("[%s: %senergy error: min,max,rms = %.6g,%.6g,%.6g]\n", getprog(),
	    (decrit == decrit0 ? "" : "warning: "), demin, demax, derms);
}

//  forcecalc: compute forces on test bodies (and zombie bodies);
//  if running in parallel, distribute forces to all processes.
//  _____________________________________________________________

local void forcecalc(void)
{
  real *apbuf[nproc], *apptr;

#if defined(NBDGRAV)
  sumforces(btab, nbody, gtab, ngrav, xyzscale, eps*eps);
#if defined(MPICODE)
  for (int n = 0; n < nproc; n++)		// alloc buffers for all proc.
    apbuf[n] = (real *) allocate(4 * sizeof(real) * nbody);
  apptr = apbuf[iproc];				// copy data to local buffer
  for (bodyptr bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    *apptr++ = Acc(bp)[0];
    *apptr++ = Acc(bp)[1];
    *apptr++ = Acc(bp)[2];
    *apptr++ = Phi(bp);
  }
  for (int n = 0; n < nproc; n++) {		// share data among processes
#if defined(DOUBLEPREC)
    MPI_Bcast((void *) apbuf[n], 4 * nbody, MPI_DOUBLE, n, MPI_COMM_WORLD);
#else
    MPI_Bcast((void *) apbuf[n], 4 * nbody, MPI_FLOAT, n, MPI_COMM_WORLD);
#endif
  }
  for (bodyptr bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    CLRV(Acc(bp));
    Phi(bp) = 0.0;
  }
  for (int n = 0; n < nproc; n++) {		// sum data in strict sequence
    apptr = apbuf[n];
    for (bodyptr bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      Acc(bp)[0] += *apptr++ / nproc;		// average accelerations
      Acc(bp)[1] += *apptr++ / nproc;
      Acc(bp)[2] += *apptr++ / nproc;
      Phi(bp)    += *apptr++ / nproc;		// average potential
    }
  }
  for (int n = 0; n < nproc; n++)
    free(apbuf[n]);
#endif
#if defined(ZOMGRAV)
  gspforces(gtab, ngrav, ggsp);			// compute forces on zombies
#endif
#elif defined(GSPGRAV)
  gspforces(btab, nbody, ggsp);
#elif defined(HMPGRAV)
  hmpforces(btab, nbody, Mhmp, ahmp, bhmp, tol);
#endif
}

//  sumforces: compute forces on test bodies due to massive bodies;
//  massive body distribution can be rescaled along X,Y,Z axies.
//  _______________________________________________________________

local void sumforces(bodyptr btab, int nbody, bodyptr gtab, int ngrav,
		     real xyzscale[3], real eps2)
{
  int i;
  bodyptr bp, gp;
  real phiB[nbody], mG, dr2, dr2i, dr1i, mdr1i, mdr3i;
  vector accB[nbody], posB[nbody], posG, dr;

  for (i = 0, bp = btab; i < nbody; i++, bp = NextBody(bp)) {
    phiB[i] = 0.0;
    CLRV(accB[i]);
    SETV(posB[i], Pos(bp));
  }
  for (gp = gtab; gp < NthBody(gtab, ngrav); gp = NextBody(gp)) {
    for (int k = 0; k < 3; k++)
      posG[k] = xyzscale[k] * Pos(gp)[k];
    mG = Mass(gp);
    for (i = 0; i < nbody; i++) {
      DOTPSUBV(dr2, dr, posG, posB[i]);
      dr2i = ((real) 1.0) / (dr2 + eps2);
      dr1i = rsqrt(dr2i);
      mdr1i = mG * dr1i;
      mdr3i = mdr1i * dr2i;
      phiB[i] -= mdr1i;
      ADDMULVS(accB[i], dr, mdr3i);
    }
  }
  for (i = 0, bp = btab; i < nbody; i++, bp = NextBody(bp)) {
    Phi(bp) = phiB[i];
    SETV(Acc(bp), accB[i]);
  }
}

//  gspforces: compute forces on test or massive bodies due to GSP.
//  _______________________________________________________________

local void gspforces(bodyptr btab, int nbody, gsprof *ggsp)
{
  bodyptr bp;
  real r, mr3i;

  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    r = absv(Pos(bp));
#if !defined(ZOMGRAV)
    Phi(bp) = gsp_phi(ggsp, r);
#endif
    mr3i = gsp_mass(ggsp, r) / rqbe(r);
    MULVS(Acc(bp), Pos(bp), - mr3i);
  }
}

//  a2, b2, R, z: access macros for parameter block passed to integrands.
//  _____________________________________________________________________

#define a2(par)  (((double *) par)[0])
#define b2(par)  (((double *) par)[1])
#define R(par)   (((double *) par)[2])
#define z(par)   (((double *) par)[3])

//  int_Phi, int_aR, int_az: Integrands for potential, radial, and vertical
//  acceleration.  Coded entirely in double precision to improve behavior
//  of integration routines, which expect accurate integrands.
//  _______________________________________________________________________

local double intPhi(double u, void *par)
{
  double Delta, mu;

  Delta = (a2(par) + u) * sqrt(b2(par) + u);
  mu = sqrt(R(par)*R(par) / (a2(par) + u) + z(par)*z(par) / (b2(par) + u));
  return (0.5 / (Delta * (1+mu)*(1+mu)));
}

local double int_aR(double u, void *par)
{
  double Delta, mu;

  Delta = (a2(par) + u) * sqrt(b2(par) + u);
  mu = sqrt(R(par)*R(par) / (a2(par) + u) + z(par)*z(par) / (b2(par) + u));
  return (R(par) / ((a2(par) + u) * Delta * mu * (1+mu)*(1+mu)*(1+mu)));
}

local double int_az(double u, void *par)
{
  double Delta, mu;

  Delta = (a2(par) + u) * sqrt(b2(par) + u);
  mu = sqrt(R(par)*R(par) / (a2(par) + u) + z(par)*z(par) / (b2(par) + u));
  return (z(par) / ((b2(par) + u) * Delta * mu * (1+mu)*(1+mu)*(1+mu)));
}

//  hmpforces: compute test body forces due to spheroidal Hernquist model.
//  ______________________________________________________________________

local void hmpforces(bodyptr btab, int nbody, real M, real a, real b,
		     real tol)
{
  bodyptr bp;
  double r, mr3i, params[4], phi0, aR0, az0, abserr[3];
  static gsl_integration_workspace *wksp = NULL;
  gsl_function FPhi, F_aR, F_az;
  static double maxerr = 0.0;
  int stat[3];

  if (a == b) {					// spherical case is easy!
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      r = absv(Pos(bp));
      Phi(bp) = - M / (a + r);
      mr3i = M * rsqr(r / (a + r)) / rqbe(r);
      MULVS(Acc(bp), Pos(bp), - mr3i);
    }
  } else {					// flattened case is harder
    if (wksp == NULL) {				// on first call, initialze
      wksp = gsl_integration_workspace_alloc(1000);
      gsl_set_error_handler_off();		// handle errors below
    }
    FPhi.function = &intPhi;
    F_aR.function = &int_aR;
    F_az.function = &int_az;
    FPhi.params = F_aR.params = F_az.params = params;
    a2(params) = rsqr(a);
    b2(params) = rsqr(b);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      R(params) = rsqrt(rsqr(Pos(bp)[0]) + rsqr(Pos(bp)[1]));
      z(params) = Pos(bp)[2];
      stat[0] = gsl_integration_qagiu(&FPhi, 0.0, tol, 0.0,
				      1000, wksp, &phi0, &abserr[0]);
      stat[1] = gsl_integration_qagiu(&F_aR, 0.0, tol, 0.0,
				      1000, wksp, &aR0, &abserr[1]);
      stat[2] = gsl_integration_qagiu(&F_az, 0.0, tol, 0.0,
				      1000, wksp, &az0, &abserr[2]);
      if (stat[0] || stat[1] || stat[2])	// any errors reported?
	for (int i = 0; i < 3; i++)
	  if (stat[i] != 0 && abserr[i] > maxerr) {
	    eprintf("[%s.hmpforces: warning: %s  abserr[%d] = %g]\n",
		    getprog(), gsl_strerror(stat[i]), i+1, abserr[i]);
	    maxerr = abserr[i];			// adjust reporting threshold
	  }
      Acc(bp)[0] = - M * (Pos(bp)[0] / R(params)) * aR0;
      Acc(bp)[1] = - M * (Pos(bp)[1] / R(params)) * aR0;
      Acc(bp)[2] = - M * az0;
      Phi(bp)    = - M * phi0;
    }
  }
}

//  initcode: read initial conditions and initialize parameters.
//  ____________________________________________________________

local void initcode(void)
{
  stream istr, gstr;
  string btags[MaxBodyFields], gtags[MaxBodyFields];
  real tgrav;

  new_field(&EinitPBF, RealType, EinitTag);	// define initial energy field
  layout_body(bodytags, Precision, NDIM);	// layout required fields
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! (get_snap(istr, &btab, &nbody, &tnow, btags, TRUE, NULL) &&
	 set_member(btags, PosTag) && set_member(btags, VelTag)))
    error("%s: can't read input snapshot data\n", getprog());
#if defined(NBDGRAV)
  gstr = stropen(getparam("grav"), "r");
  get_history(gstr);
  for (int i = 0; i < iproc; i++)		// if running in parallel
    if (! skip_item(gstr))			// skip to correct grav frame
      error("%s[%d]: unexpected EOF reading grav data\n", getprog(), iproc);
  if (! (get_snap(gstr, &gtab, &ngrav, &tgrav, gtags, FALSE, NULL) &&
	 set_member(gtags, MassTag) && set_member(gtags, PosTag)))
    error("%s[%d]: can't read gravity snapshot data\n", getprog(), iproc);
  if (sscanf(getparam("xyzscale"), REALFMT "," REALFMT "," REALFMT,
	     &xyzscale[0], &xyzscale[1], &xyzscale[2]) != 3)
    error("%s: can't parse xyzscale\n", getprog());
  eps = getdparam("eps");
#endif
#if defined(GSPGRAV) || defined(ZOMGRAV)
  gstr = stropen(getparam("gsp"), "r");
  get_history(gstr);
  ggsp = gsp_read(gstr);			// read GSP for grav. field
#endif
#if defined(HMPGRAV)
  Mhmp = getdparam("M");
  ahmp = getdparam("a");
  bhmp = getdparam("b");
  tol = getdparam("tol");
#endif
  tstop = getdparam("tstop");
  dtime = getdparam("dtime");
  decrit0 = getdparam("decrit");
  otags = burststring(getparam("outputs"), ",");
  if (iproc == 0) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
  }
}

//  saveconfig: checkpoint current test (and zombie) body configurations;
//  in parallel, confmt has two %d specs. to keep processes separate.
//  _____________________________________________________________________

local void saveconfig(string confmt)
{
  string cfnm;
  stream cstr;

#if !defined(MPICODE)
  asprintf(&cfnm, confmt, ((int) floor(tnow)) % 2);
#else
  asprintf(&cfnm, confmt, iproc, ((int) floor(tnow)) % 2);
#endif
  eprintf("[%s[%d].saveconfig: writing configuration to \"%s\"]\n",
	  getprog(), iproc, cfnm);
  cstr = stropen(cfnm, "w!");
  put_history(cstr);
  put_snap(cstr, &btab, &nbody, &tnow, bodytags);
#if defined(ZOMGRAV)
  put_snap(cstr, &gtab, &ngrav, &tnow, bodytags);
#endif
  strclose(cstr);
  free(cfnm);
}
