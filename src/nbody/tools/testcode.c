/*
 * testcode.c: integrate test bodies in given potential.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "gsp.h"

#if !(defined(GSPGRAV) || defined(HQMGRAV))
#define NBDGRAV
#endif

string defv[] = {
#if defined(NBDGRAV)
				";Integrate test bodies in N-body potential",
#elif defined(GSPGRAV)
				";Integrate test bodies in GSP potential",
#elif defined(HQMGRAV)
				";Integrate test bodies in Hernquist model",
#endif
  "in=???",			";Initial conditions for test bodies",
  "out=???",			";Output stream of test body snapshots",
#if defined(NBDGRAV)
  "grav=???",			";N-body configuration generating potential",
  "eps=0.025",			";Softening parameter for force calculation",
  "frozen=true",		";If FALSE, read new grav snap for each step",
#elif defined(GSPGRAV)
  "grav=???",			";General spherical profile for potential",
#elif defined(HQMGRAV)
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
  "VERSION=1.1",		";Josh Barnes  9 June 2015",
  NULL,
};

#if defined(NBDGRAV)
local void sumforces(bodyptr btab, int nbody, bodyptr gtab, int ngrav,
		     real eps2);
#elif defined(GSPGRAV)
local void gspforces(bodyptr btab, int nbody, gsprof *gravgsp);
#elif defined(HQMGRAV)
local void hqmforces(bodyptr btab, int nbody, real M, real a, real b,
		     real tol);
#endif

#define EinitPBF  phatbody[NewBodyFields+0]	// add initial E to bodies
#define Einit(b)  SelectReal(b,EinitPBF.offset)	// accessor macro for above
#define EinitTag  "Einit"			// field name for above

string bodytags[MaxBodyFields] = {
  PosTag, VelTag, MassTag, PhiTag, AccTag, EinitTag, NULL
};

int main(int argc, string argv[])
{
  stream istr, ostr, gstr;
  real tnow, tgrav, eps2, tstop, dtime, tout, Mhqm, ahqm, bhqm, tol;
  real decrit, epot0, demin, demax, derms, de2avg, enow, denow;
  int nbody, ngrav;
  bodyptr btab = NULL, gtab = NULL, bp;
  string bdtags[MaxBodyFields], grtags[MaxBodyFields], *optags;
  gsprof *gravgsp = NULL;
  bool decrit_inc = FALSE;

  initparam(argv, defv);
  new_field(&EinitPBF, RealType, EinitTag);	// define initial energy field
  layout_body(bodytags, Precision, NDIM);	// layout necessary fields
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab, &nbody, &tnow, bdtags, TRUE))
    error("%s: can't read input snapshot\n", getprog());
  if (! (set_member(bdtags, PosTag) && set_member(bdtags, VelTag)))
    error("%s: required data missing from input snapshot\n", getprog());
#if defined(NBDGRAV)
  gstr = stropen(getparam("grav"), "r");
  get_history(gstr);
  if (! get_snap(gstr, &gtab, &ngrav, &tgrav, grtags, FALSE))
    error("%s: can't read gravity snapshot\n", getprog());
  if (! (set_member(grtags, MassTag) && set_member(grtags, PosTag)))
    error("%s: required data missing from gravity snapshot\n", getprog());
  eps2 = rsqr(getdparam("eps"));
#elif defined(GSPGRAV)
  gstr = stropen(getparam("grav"), "r");
  get_history(gstr);
  gravgsp = get_gsprof(gstr);			// read GSP for grav. field
#elif defined(HQMGRAV)
  Mhqm = getdparam("M");
  ahqm = getdparam("a");
  bhqm = getdparam("b");
  tol = getdparam("tol");
#endif
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  tstop = getdparam("tstop");
  dtime = getdparam("dtime");
  decrit = getdparam("decrit");
  optags = burststring(getparam("outputs"), ",");
#if defined(NBDGRAV)
  sumforces(btab, nbody, gtab, ngrav, eps2);	// prime the pump...
#elif defined(GSPGRAV)
  gspforces(btab, nbody, gravgsp);
#elif defined(HQMGRAV)
  hqmforces(btab, nbody, Mhqm, ahqm, bhqm, tol);
#endif

  epot0 = 0.0;					// use as energy scale
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    epot0 += Phi(bp) / nbody;			// compute avg. potential
    Einit(bp) = Phi(bp) + dotvp(Vel(bp), Vel(bp)) / 2;
  }
  eprintf("[%s: initial average potential = %g]\n", getprog(), epot0);
  put_snap(ostr, &btab, &nbody, &tnow, optags);
  fflush(NULL);
  tout = tnow + getdparam("dtout");
  demin = demax = derms = 0.0;			// track maximum errors
  while (tnow < tstop) {			// enter main loop
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      ADDMULVS(Vel(bp), Acc(bp), 0.5 * dtime);	// step velocities by dt/2
      ADDMULVS(Pos(bp), Vel(bp), dtime);	// step positions by dt
    }
    tnow = tnow + dtime;			// step time to new value
#if defined(NBDGRAV)
    if (! getbparam("frozen"))
      if (! get_snap(gstr, &gtab, &ngrav, &tgrav, grtags, TRUE))
	error("%s: can't read gravity snapshot\n", getprog());
    sumforces(btab, nbody, gtab, ngrav, eps2);	// get new accelerations
#elif defined(GSPGRAV)
    gspforces(btab, nbody, gravgsp);
#elif defined(HQMGRAV)
    hqmforces(btab, nbody, Mhqm, ahqm, bhqm, tol);
#endif
    de2avg = 0.0;
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      ADDMULVS(Vel(bp), Acc(bp), 0.5 * dtime);	// step velocities by dt/2
      enow = 0.5 * dotvp(Vel(bp), Vel(bp)) + Phi(bp);
      denow = (enow - Einit(bp)) / ABS(epot0);	// compute rel. energy change
      demin = MIN(demin, denow);
      demax = MAX(demax, denow);
      de2avg += rsqr(denow) / nbody;
    }
    derms = MAX(derms, rsqrt(de2avg));
    if (demin < -decrit || demax > decrit) {
      eprintf("[%s: warning: energy error exceeds %.4e at time = %-12.8f\n"
	      " min,max,rms = %.6g,%.6g,%.6g  threshold now %.4e]\n",
	      getprog(), decrit, tnow,
	      demin, demax, rsqrt(de2avg), decrit * rsqrt(2.0));
      decrit = decrit * rsqrt(2.0);
      decrit_inc = TRUE;
    }
    if (tout <= tnow) {
      put_snap(ostr, &btab, &nbody, &tnow, optags);
      tout = tout + getdparam("dtout");
    }
    fflush(NULL);
  }
  eprintf(decrit_inc ?
	  "[%s: WARNING: energy error: min,max,rms = %.6g,%.6g,%.6g]\n" :
	  "[%s: energy error: min,max,rms = %.6g,%.6g,%.6g]\n",
	  getprog(), demin, demax, derms);
  return (0);
}

#if defined(NBDGRAV)

local void sumforces(bodyptr btab, int nbody, bodyptr gtab, int ngrav,
		     real eps2)
{
  bodyptr bp, gp;
  double phi0, acc0[NDIM];
  vector pos0, dr;
  real dr2, dr2i, dr1i, mdr1i, mdr3i;

  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    phi0 = 0.0;
    CLRV(acc0);
    SETV(pos0, Pos(bp));
    for (gp = gtab; gp < NthBody(gtab, ngrav); gp = NextBody(gp)) {
      DOTPSUBV(dr2, dr, Pos(gp), pos0);
      dr2i = ((real) 1.0) / (dr2 + eps2);
      dr1i = rsqrt(dr2i);
      mdr1i = Mass(gp) * dr1i;
      mdr3i = mdr1i * dr2i;
      phi0 -= mdr1i;
      ADDMULVS(acc0, dr, mdr3i);
    }
    Phi(bp) = phi0;
    SETV(Acc(bp), acc0);
  }
}

#endif

#if defined(GSPGRAV)

local void gspforces(bodyptr btab, int nbody, gsprof *gravgsp)
{
  bodyptr bp;
  real r, mr3i;

  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    r = absv(Pos(bp));
    Phi(bp) = phi_gsp(gravgsp, r);
    mr3i = mass_gsp(gravgsp, r) / rqbe(r);
    MULVS(Acc(bp), Pos(bp), -mr3i);
  }
}

#endif

#if defined(HQMGRAV)

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define a2(par)  (((double *) par)[0])
#define b2(par)  (((double *) par)[1])
#define R(par)   (((double *) par)[2])
#define z(par)   (((double *) par)[3])

local double intPhi(double u, void *par)
{
  real Delta, mu;

  Delta = (a2(par) + u) * rsqrt(b2(par) + u);
  mu = rsqrt(rsqr(R(par)) / (a2(par) + u) + rsqr(z(par)) / (b2(par) + u));
  return (0.5 / (Delta * rsqr(1 + mu)));
}

local double int_aR(double u, void *par)
{
  real Delta, mu;

  Delta = (a2(par) + u) * rsqrt(b2(par) + u);
  mu = rsqrt(rsqr(R(par)) / (a2(par) + u) + rsqr(z(par)) / (b2(par) + u));
  return (R(par) / ((a2(par) + u) * Delta * mu * rqbe(1 + mu)));
}

local double int_az(double u, void *par)
{
  real Delta, mu;

  Delta = (a2(par) + u) * rsqrt(b2(par) + u);
  mu = rsqrt(rsqr(R(par)) / (a2(par) + u) + rsqr(z(par)) / (b2(par) + u));
  return (z(par) / ((b2(par) + u) * Delta * mu * rqbe(1 + mu)));
}

local void hqmforces(bodyptr btab, int nbody, real M, real a, real b,
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
	    eprintf("[%s.hqmforces: warning: %s  abserr[%d] = %g]\n",
		    getprog(), gsl_strerror(stat[i]), i+1, abserr[i]);
	    maxerr = abserr[i];			// adjust reporting threshold
	  }
      Phi(bp) = - M * phi0;
      Acc(bp)[0] = - M * (Pos(bp)[0] / R(params)) * aR0;
      Acc(bp)[1] = - M * (Pos(bp)[1] / R(params)) * aR0;
      Acc(bp)[2] = - M * az0;
    }
  }
}

#endif
