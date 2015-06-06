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

string defv[] = {
#if defined(GSPGRAV)
				";Integrate test bodies in GSP potential",
#elif defined(HQMGRAV)
				";Integrate test bodies in Hernquist model",
#else
				";Integrate test bodies in N-body potential",
#endif
  "in=???",			";Initial conditions for test bodies",
  "out=???",			";Output stream of test body snapshots",
#if defined(GSPGRAV)
  "grav=???",			";General spherical profile for potential",
#elif defined(HQMGRAV)
  "M=1.0",			";Total mass of Hernquist model",
  "a=1.0",			";Major axis of Hernquist model",
  "b=1.0",			";Minor axis of Hernquist model",
#else
  "grav=???",			";N-body configuration generating potential",
  "eps=0.025",			";Softening parameter for force calculation",
  "frozen=true",		";If FALSE, read new grav snap for each step",
#endif
  "tstop=2.0",			";Duration of integration",
  "dtime=1/32",			";Leap-frog integration time-step",
  "decrit=1.0e-5",		";Initial energy change tolerance",
  "outputs=" PosTag "," VelTag,	";Per-body data arrays written to output",
  "dtout=1/4",			";Time interval between output frames",
  "VERSION=1.0",		";Josh Barnes  23 May 2015",
  NULL,
};

#define EinitPBF  phatbody[NewBodyFields+0]	// add initial E to bodies
#define Einit(b)  SelectReal(b,EinitPBF.offset)	// accessor macro for above
#define EinitTag  "Einit"			// field name for above

string bodytags[MaxBodyFields] = {
  PosTag, VelTag, MassTag, PhiTag, AccTag, EinitTag, NULL
};

local void sumforces(bodyptr btab, int nbody, bodyptr gtab, int ngrav,
		     real eps2);

local void gspforces(bodyptr btab, int nbody, gsprof *gravgsp);

#if defined(HQMGRAV)
local void hqmforces(bodyptr btab, int nbody, real M, real a, real b);
#endif

int main(int argc, string argv[])
{
  stream istr, ostr, gstr;
  real tnow, tgrav, eps2, tstop, dtime, tout;
  real decrit, epot0, demin, demax, derms, de2avg, enow, denow;
  int nbody, ngrav;
  bodyptr btab = NULL, gtab = NULL, bp;
  string bdtags[MaxBodyFields], grtags[MaxBodyFields], *optags;
  gsprof *gravgsp = NULL;

  initparam(argv, defv);
  new_field(&EinitPBF, RealType, EinitTag);	// define initial energy field
  layout_body(bodytags, Precision, NDIM);	// layout necessary fields
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab, &nbody, &tnow, bdtags, TRUE))
    error("%s: can't read input snapshot\n", getprog());
  if (! (set_member(bdtags, PosTag) && set_member(bdtags, VelTag)))
    error("%s: required data missing from input snapshot\n", getprog());
#if defined(GSPGRAV)
  gstr = stropen(getparam("grav"), "r");
  get_history(gstr);
  gravgsp = get_gsprof(gstr);			// read GSP for grav. field
#elif !defined(HQMGRAV)
  gstr = stropen(getparam("grav"), "r");
  get_history(gstr);
  if (! get_snap(gstr, &gtab, &ngrav, &tgrav, grtags, FALSE))
    error("%s: can't read gravity snapshot\n", getprog());
  if (! (set_member(grtags, MassTag) && set_member(grtags, PosTag)))
    error("%s: required data missing from gravity snapshot\n", getprog());
  eps2 = rsqr(getdparam("eps"));
#endif
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  tstop = getdparam("tstop");
  dtime = getdparam("dtime");
  decrit = getdparam("decrit");
  optags = burststring(getparam("outputs"), ",");
#if defined(GSPGRAV)
  gspforces(btab, nbody, gravgsp);		// prime the pump...
#elif defined(HQMGRAV)
  hqmforces(btab, nbody, getdparam("M"), getdparam("a"), getdparam("b"));
#else
  sumforces(btab, nbody, gtab, ngrav, eps2);
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
#if defined(GSPGRAV)
    gspforces(btab, nbody, gravgsp);		// get new accelerations
#elif defined(HQMGRAV)
    hqmforces(btab, nbody, getdparam("M"), getdparam("a"), getdparam("b"));
#else
    if (! getbparam("frozen"))
      if (! get_snap(gstr, &gtab, &ngrav, &tgrav, grtags, TRUE))
	error("%s: can't read gravity snapshot\n", getprog());
    sumforces(btab, nbody, gtab, ngrav, eps2);
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
    }
    if (tout <= tnow) {
      put_snap(ostr, &btab, &nbody, &tnow, optags);
      tout = tout + getdparam("dtout");
    }
    fflush(NULL);
  }
  eprintf("[%s: energy error: min,max,rms = %.6g,%.6g,%.6g]\n",
	  getprog(), demin, demax, derms);
  return (0);
}

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

local double phiint(double u, void *par)
{
  real Delta, mu;

  Delta = (a2(par) + u) * rsqrt(b2(par) + u);
  mu = rsqrt(rsqr(R(par)) / (a2(par) + u) + rsqr(z(par)) / (b2(par) + u));
  return (0.5 / (Delta * rsqr(1 + mu)));
}

local double aRint(double u, void *par)
{
  real Delta, mu;

  Delta = (a2(par) + u) * rsqrt(b2(par) + u);
  mu = rsqrt(rsqr(R(par)) / (a2(par) + u) + rsqr(z(par)) / (b2(par) + u));
  return (R(par) / ((a2(par) + u) * Delta * mu * rqbe(1 + mu)));
}

local double azint(double u, void *par)
{
  real Delta, mu;

  Delta = (a2(par) + u) * rsqrt(b2(par) + u);
  mu = rsqrt(rsqr(R(par)) / (a2(par) + u) + rsqr(z(par)) / (b2(par) + u));
  return (z(par) / ((b2(par) + u) * Delta * mu * rqbe(1 + mu)));
}

#define AEPS  1.0e-6
#define REPS  1.0e-6

local void hqmforces(bodyptr btab, int nbody, real M, real a, real b)
{
  bodyptr bp;
  double r, mr3i, params[4], phiI, phiE, aRI, aRE, azI, azE;
  static gsl_integration_workspace *wksp = NULL;
  gsl_function phiF, aRF, azF;

  if (a == b) {					// spherical case is easy!
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      r = absv(Pos(bp));
      Phi(bp) = - M / (a + r);
      mr3i = M * rsqr(r / (a + r)) / rqbe(r);
      MULVS(Acc(bp), Pos(bp), - mr3i);
    }
  } else {					// flattened case is harder
    if (wksp == NULL) {				// first call; initialze
      wksp = gsl_integration_workspace_alloc(1000);
      gsl_set_error_handler_off();		// live dangerously
    }
    phiF.function = &phiint;
    aRF.function = &aRint;
    azF.function = &azint;
    phiF.params = aRF.params = azF.params = params;
    a2(params) = rsqr(a);
    b2(params) = rsqr(b);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
      R(params) = rsqrt(rsqr(Pos(bp)[0]) + rsqr(Pos(bp)[1]));
      z(params) = Pos(bp)[2];
      gsl_integration_qagiu(&phiF, 0.0, AEPS, REPS, 1000, wksp, &phiI, &phiE);
      gsl_integration_qagiu(&aRF, 0.0, AEPS, REPS, 1000, wksp, &aRI, &aRE);
      gsl_integration_qagiu(&azF, 0.0, AEPS, REPS, 1000, wksp, &azI, &azE);
      Phi(bp) = - M * phiI;
      Acc(bp)[0] = - M * (Pos(bp)[0] / R(params)) * aRI;
      Acc(bp)[1] = - M * (Pos(bp)[1] / R(params)) * aRI;
      Acc(bp)[2] = - M * azI;
    }
  }
}

#endif
