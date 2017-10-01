/*
 * sphcode.c: hierarchical SPH/N-body simulation code family.
 * Copyright (c) 2017 by Joshua E. Barnes, Honolulu, Hawai'i.
 */

#define global					// don't default to extern
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "datatypes.h"
#include "gsp.h"
#include "sphcode.h"
#include "kdtree.h"
#include "smooth.h"
#include "fixbody.h"

//  Describe code variant on basis of compile-time flags.
//  _____________________________________________________

#if defined(ENTROPY)
#  define THERMOPT		";Entropy formulation,",
#  if defined(ADIABATIC)
#    define GASOPT		";Adiabatic gas (constant entf),",
#  else
#    define GASOPT		";Non-adiabatic gas (shock heating),",
#  endif
#else
#  define THERMOPT		";Energy formulation,",
#  if defined(ISOTHERMAL)
#    define GASOPT		";Isothermal gas (constant uint),",
#  elif defined(RADIATING)
#    if defined(DIFFUSING)
#      define GASOPT		";Gas radiates (diffusion model),",
#    elif defined(OPAQUE)
#      define GASOPT		";Gas radiates (opaque model),",
#    else
#      define GASOPT		";Gas radiates (optically thin),",
#    endif
#  elif defined(CONDUCTING)
#    define GASOPT		";Gas conducts heat,",
#  else
#    define GASOPT		";Ideal gas (shock heating),",
#  endif
#endif

#if defined(STARFORM) && defined(MASSLOSS)
#  define STAROPT               ";Star formation and mass loss,",
#elif defined(STARFORM) && defined(COMPVISC)
#  define STAROPT		";Star formation (using udotvis),",
#elif defined(STARFORM)
#  define STAROPT		";Star formation (using udot),",
#else
#  define STAROPT
#endif

#if defined(GRAVITY) && !defined(NOSOFTCORR)
#  define GRAVOPT		";Self-consistent gravity (soft corr).",
#elif defined(GRAVITY) && defined(NOSOFTCORR)
#  define GRAVOPT		";Self-consistent gravity (no soft corr).",
#elif defined(EXTGRAV)
#  define GRAVOPT		";External gravitational field.",
#elif defined(NOACCEL)
#  define GRAVOPT		";No accelerations computed.",
#endif

//  Command-line parameters and defaults.
//  _____________________________________

string defv[] = {		";SPH/N-body simulation code.",
				THERMOPT GASOPT STAROPT GRAVOPT
  "in=",			";Input file with initial conditions",
  "out=",			";Output file patern for SPH frames",
  "save=",			";Write state file as code runs",
  "restore=",			";Continue run from state file",
  "gamma=5/3",			";Ratio of specific heats",
#if defined(RADIATING)
  "uradpk=1.0",			";Uinternal at peak of cooling curve",
  "lambpk=1.0",			";Peak cooling rate at unit density",
#endif
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
  "uintmax=0.0",		";Enforce maximum uint if gt zero",
#endif
#if defined(DIFFUSING)
  "sigmastar=1.0",		";Stefan-Boltzmann parameter.",
				";Units are flux / specific_energy^4.",
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
  "opacity=1.0",		";Radiation opacity parameter",
#endif
#if defined(CONDUCTING)
  "conduct=1.0",		";Thermal conduction parameter",
#endif
#if defined(STARFORM)
  "cstar=0.0,0.0,0.0",		";Star formation law constants",
  "nstar=1.0,1.0,1.0",		";Star formation density indicies",
  "mstar=0.0,0.0,0.0",		";Star formation shock indicies",
  "starlog=starlog.txt",	";Star formation and evolution log",
  "starseed=12345",		";Random number seed for star form.",
#if defined(MASSLOSS)
  "tau_ml=0.01",		";Mass-loss timescale (2.5 Myr)",
  "beta_ml=12.0",		";Mass-loss power-law index",
#endif
#endif
  "alpha=1.0",			";Bulk viscosity parameter",
  "beta=2.0",			";vN-R viscosity parameter",
  "nsmooth=40",			";Bodies in smoothing volume",
  "nbucket=16",			";Bodies in leaves of KD tree",
  "slope=0.0",			";Kernel slope at origin",
  "courant=0.25",		";Courant condition parameter",
  "dtime=1/256",		";Basic integration timestep",
  "fdrag=0.0",			";Velocity damping factor (1/t)",
#if defined(GRAVITY)
  "eps=0.01",			";Gravitational smoothing length",
  "usequad=true",		";If false, don't use quad moments",
  "theta=1.0",			";Force accuracy parameter",
#elif defined(EXTGRAV)
  "gravgsp=",			";Input GSP for external gravity",
#endif
  "options=",			";Assorted flags for simulation control.",
				";Options: corrfunc, levelhist, new-tout,",
				";nolimit, fixstep, lockstep, reset-time,",
				";bh86, sw94, theta-eff.",
  "outputs="
#if defined(ENTROPY) && !defined(ADIABATIC)
             EntFuncTag ","
#elif !defined(ENTROPY) && !defined(ISOTHERMAL)
             UinternTag ","
#endif
             PosTag "," VelTag, ";Particle data written to output file.",
				";First frame also has mass and gas state.",

  "tstop=2.0",			";Time to stop integration",
  "dtout=1/16",			";Data output timestep",
  "testbody=65536",		";Number of gas particles for test data",
  "testseed=123",		";Random number seed for test data",
  "testuint=0.05",		";Internal energy for test data",
  "stream=",			";Output file pattern for snapshot stream",
  "log=",			";Output file name for simulation log",
				";Defaults to stdout unless stream=-",
  "VERSION=1.6",		";Joshua Barnes  3 May 2017",
  NULL,
};

//  Local procedure prototypes.
//  ___________________________

local void startrun(void);			// initialize system state
local void newrun(void);			// get params & input data
local void testdata(void);			// generate test data
local void oldrun(void);			// restore state & params
local void initstep(void);			// init multistep algorithm
local void macrostep(void);			// advance by basic step
local void microstep(real *);			// advance by minimum step
local void update(int, real);			// update forces and levels
local void evaluate(bool);			// do force calculation
local void gravcalc(bool);			// gravitational part
local void sphcalc(bool);			// pressure, etc part
local void sph_density(smxptr, int, int);	// compute SPH densities
local void sph_depth(smxptr, int, int);		// estimate optical depth
local void sph_derivs(smxptr, int, int);	// compute SPH derivatives
local real coolfunc(bodyptr);			// compute cooling function
local void setlevels(smxptr);			// adjust timestep levels
local void starform(real);			// handle star formation
local void massloss(real);			// handle star destruction

//  main: toplevel routine for hierarchical SPH/N-body code.
//  ________________________________________________________

int main(int argc, string argv[])
{
  initparam(argv, defv);			// initialize param access
  headline = defv[0] + 1;			// set ident. message
  startrun();					// get params & input data
  startoutput(defv + 1);			// open streams for output
  if (nstep == 0)				// starting new calculation?
    initstep();					// then prime the integrator
  outputdata();					// report initial diagsnostics
  while (tstop - tnow > 0.01 * dtime) {		// while not past tstop
    macrostep();				// advance step by step
    outputdata();				// do output after each
  }
  return (0);
}

//  startrun: startup hierarchical N-body code.
//  ___________________________________________

local void startrun(void)
{
  define_body(sizeof(body), Precision, NDIM);	// setup interface to dynbody
  define_body_offset(TypeTag, BodyOffset(Type));
  define_body_offset(PosTag, BodyOffset(Pos));
  define_body_offset(VelTag, BodyOffset(Vel));
  define_body_offset(MassTag, BodyOffset(Mass));
  define_body_offset(SmoothTag, BodyOffset(Smooth));
  define_body_offset(PhiTag, BodyOffset(Phi));
  define_body_offset(AccTag, BodyOffset(Acc));
  define_body_offset(RhoTag, BodyOffset(Rho));
#if defined(ENTROPY)
  define_body_offset(EntFuncTag, BodyOffset(EntFunc));
#else
  define_body_offset(UinternTag, BodyOffset(Uintern));
#endif
  define_body_offset(UdotIntTag, BodyOffset(UdotInt));
#if defined(RADIATING)
  define_body_offset(UdotRadTag, BodyOffset(UdotRad));
#endif
#if defined(COMPVISC)
  define_body_offset(UdotVisTag, BodyOffset(UdotVis));
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
  define_body_offset(TauTag, BodyOffset(Tau));
#endif
#if defined(STARFORM) || defined(MASSLOSS)
  define_body_offset(BirthTag, BodyOffset(Birth));
#endif
#if defined(MASSLOSS)
  define_body_offset(DeathTag, BodyOffset(Death));
#endif
  infile = getparam("in");			// set up I/O file names
  outpatn = getparam("out");
  savepatn = getparam("save");
  restfile = getparam("restore");
  strmpatn = getparam("stream");
  logfile = getparam("log");
  if (! strnull(infile) && ! strnull(restfile))
      error("%s: can't read both input and restore files\n", getprog());
  if (strnull(restfile)) {			// no restore file; start run
    if (strnull(infile))
      testdata();				// make up initial data
    else
      inputdata();				// read inital data
    newrun();
  } else {
    restorestate();				// read in state file
    oldrun();
  }
}

//  testdata: generate data for test case: collapse of 1/r density profile
//  from isothermal initial state (Evrard 1988, Hernquist & Katz 1989).
//  ______________________________________________________________________

local void testdata(void)
{
  bodyptr p;
  real uint, gam0, r;

  nbody = getiparam("testbody");
  init_random(getiparam("testseed"));
  uint = getdparam("testuint");
  gam0 = getdparam("gamma");			// use user-specified gamma
  btab = (bodyptr) allocate(nbody * sizeof(body));
  for (p = btab; p < btab+nbody; p++) {
    Type(p) = BODY | GAS;
    Mass(p) = 1.0 / nbody;
    r = rsqrt(xrandom(0.0, 1.0));		// make 1/r density profile
    pickshell(Pos(p), NDIM, r);
    CLRV(Vel(p));
#if defined(ENTROPY)
    EntFunc(p) = (gam0 - 1) * uint * rpow(2 * PI * r, gam0 - 1);
#else
    Uintern(p) = uint;
#endif
  }
  ngas = nbody;
  tnow = 0.0;
}

//  RFMT3: scanf format for three real parameters.
//  ______________________________________________

#if defined(DOUBLEPREC)
#  define RFMT3  "%lf,%lf,%lf"
#else
#  define RFMT3  "%f,%f,%f"
#endif

//  newrun: set command-line parameters and read in input data.
//  ___________________________________________________________

local void newrun(void)
{
  gamma0 = getdparam("gamma");			// get input parameters
#if defined(RADIATING)
  uradpk = getdparam("uradpk");
  lambpk = getdparam("lambpk");
#endif
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
  uintmax = getdparam("uintmax");
#endif
#if defined(DIFFUSING)
  sigmastar = getdparam("sigmastar");
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
  opacity = getdparam("opacity");
#endif
#if defined(CONDUCTING)
  conduct = getdparam("conduct");
#endif
#if defined(STARFORM)
  init_random(getiparam("starseed"));
  if (sscanf(getparam("cstar"), RFMT3, cstar, cstar+1, cstar+2) != 3 ||
      sscanf(getparam("nstar"), RFMT3, nstar, nstar+1, nstar+2) != 3 ||
      sscanf(getparam("mstar"), RFMT3, mstar, mstar+1, mstar+2) != 3)
    error("%s: can't parse star formation parameters\n", getprog());
#if defined(MASSLOSS)
  tau_ml = getdparam("tau_ml");
  beta_ml = getdparam("beta_ml");
#endif
#endif
  alpha = getdparam("alpha");
  beta = getdparam("beta");
  nsmooth = getiparam("nsmooth");
  nbucket = getiparam("nbucket");
  slope0 = getdparam("slope");
  courant = getdparam("courant");
  dtime = getdparam("dtime");
  fdrag = getdparam("fdrag");
#if defined(GRAVITY)
  eps = getdparam("eps");
  usequad = getbparam("usequad");
  theta = getdparam("theta");
#elif defined(EXTGRAV)
  if (! strnull(getparam("gravgsp")))
    inputgrav(getparam("gravgsp"));		// get external grav field
  else
    gravgsp = NULL;
#endif
  options = getparam("options");
  outputs = getparam("outputs");
  tstop = getdparam("tstop");
  dtout = getdparam("dtout");
  if (scanopt(options, "reset-time"))		// reset starting time?
    tnow = 0.0;
  tout = tnow;					// schedule first output
  rsize = 1.0;					// start root w/ unit cube
  nstep = 0;					// begin counting steps
  eradiate = 0.0;				// zero radiated energy
}

//  oldrun: restore state from file and set new parameters.
//  _______________________________________________________

local void oldrun(void)
{
  if (getparamstat("gamma") & ARGPARAM)		// reset explicit params
    gamma0 = getdparam("gamma");
#if defined(RADIATING)
  if (getparamstat("uradpk") & ARGPARAM)
    uradpk = getdparam("uradpk");
  if (getparamstat("lambpk") & ARGPARAM)
    lambpk = getdparam("lambpk");
#endif
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
  if (getparamstat("uintmax") & ARGPARAM)
    uintmax = getdparam("uintmax");
#endif
#if defined(DIFFUSING)
  if (getparamstat("sigmastar") & ARGPARAM)
    sigmastar = getdparam("sigmastar");
#endif    
#if defined(DIFFUSING) || defined(OPAQUE)
  if (getparamstat("opacity") & ARGPARAM)
    opacity = getdparam("opacity");
#endif
#if defined(CONDUCTING)
  if (getparamstat("conduct") & ARGPARAM)
    conduct = getdparam("conduct");
#endif
#if defined(STARFORM)
  if (getparamstat("starseed") & ARGPARAM)
    init_random(getiparam("starseed"));
  if (getparamstat("cstar") & ARGPARAM)
    if (sscanf(getparam("cstar"), RFMT3, cstar, cstar+1, cstar+2) != 3)
      error("%s: can't parse star formation parameters\n", getprog());
  if (getparamstat("nstar") & ARGPARAM)
    if (sscanf(getparam("nstar"), RFMT3, nstar, nstar+1, nstar+2) != 3)
      error("%s: can't parse star formation parameters\n", getprog());
  if (getparamstat("mstar") & ARGPARAM)
    if (sscanf(getparam("mstar"), RFMT3, mstar, mstar+1, mstar+2) != 3)
      error("%s: can't parse star formation parameters\n", getprog());
#if defined(MASSLOSS)
  if (getparamstat("tau_ml") & ARGPARAM)
    tau_ml = getdparam("tau_ml");
  if (getparamstat("beta_ml") & ARGPARAM)
    beta_ml = getdparam("beta_ml");
#endif
#endif
  if (getparamstat("alpha") & ARGPARAM)
    alpha = getdparam("alpha");
  if (getparamstat("beta") & ARGPARAM)
    beta = getdparam("beta");
  if (getparamstat("nsmooth") & ARGPARAM)
    nsmooth = getiparam("nsmooth");
  if (getparamstat("nbucket") & ARGPARAM)
    nbucket = getiparam("nbucket");
  if (getparamstat("slope") & ARGPARAM)
    slope0 = getdparam("slope");
  if (getparamstat("dtime") & ARGPARAM)
    eprintf("[%s: warning: cannot change timestep]\n", getprog());
  if (getparamstat("courant") & ARGPARAM)
    courant = getdparam("courant");
  if (getparamstat("fdrag") & ARGPARAM)
    fdrag = getdparam("fdrag");

#if defined(GRAVITY)
  if (getparamstat("eps") & ARGPARAM)
    eps = getdparam("eps");
  if (getparamstat("usequad") & ARGPARAM)
    usequad = getbparam("usequad");
  if (getparamstat("theta") & ARGPARAM)
    theta = getdparam("theta");
#endif
  if (getparamstat("options") & ARGPARAM)
    options = getparam("options");
  if (getparamstat("outputs") & ARGPARAM)
    outputs = getparam("outputs");
  if (getparamstat("tstop") & ARGPARAM)
    tstop = getdparam("tstop");
  if (getparamstat("dtout") & ARGPARAM)
    dtout = getdparam("dtout");
  if (scanopt(options, "new-tout"))		// if output time reset
    tout = tnow + dtout;			// then offset from now
}

//  initstep: prepare for multi-step integration.
//  _____________________________________________

local void initstep(void)
{
  real ftmp, dt;
  bodyptr p;

  for (p = btab; p < btab+nbody; p++)	{	// loop over all bodies
    Freq(p) = 1.0 / dtime;			// set default freqency
    Flags(p) = INCLUDE | UPDATE;		// do full force evaluation
  }
  evaluate(TRUE);				// initial force calculation
  levmax = 0;
  for (p = btab; p < btab+nbody; p++) {		// loop over all bodies
    CurLevel(p) = NewLevel(p);			// init multi-step level
    levmax = MAX(levmax, CurLevel(p));		// determine maximum level
    dt = dtime / (2 << CurLevel(p));		// compute half timestep
    SETV(Vmid(p), Vel(p));			// init midpoint velocity
    ADDMULVS(Vmid(p), Acc(p), dt);		// forward by half a step
  }
}

//  macrostep: advance system by dtime.
//  ___________________________________

local void macrostep(void)
{
  real fstep;
  int lev, n, m;

#if defined(TRACESTEP)
  if (logstream != NULL)
    fprintf(logstream, "\n    %12s%8s%8s\n", "fstep", "level", "update");
#endif
  fstep = 0.0;					// count fractional steps
  microstep(&fstep);				// advance system variables
  while (fstep < 1.0) {				// loop until step complete
    lev = levmax;				// find level of next step
    n = (1 << levmax) * fstep;			// form microstep count
    for (m = n ^ (n-1); m > 1; m = m >> 1)
      lev--;
    update(lev, fstep);				// update bodies at level
    microstep(&fstep);				// advance variables again
  }
  if (fstep > 1.0)
    error("%s: fstep = %f\n", getprog(), fstep);
  nstep++;					// advance macrostep count
  update(0, fstep);				// and update all bodies
}

//  microstep: advance system by short time-step.
//  _____________________________________________

local void microstep(real *fstep)
{
  real dt, edotrad = 0.0, udotrad;
  vector vtmp;
  
  dt = dtime / (1 << levmax);			// find current time-step
  for (bodyptr p = btab; p < btab+nbody; p++) {
    SETV(vtmp, Vel(p));
    ADDMULVS(vtmp, Acc(p), 0.5 * dt);
    ADDMULVS(Pos(p), vtmp, dt);
    ADDMULVS(Vel(p), Acc(p), dt);
    if (Gas(p)) {
#if defined(ENTROPY)
#  if defined(ADIABATIC)
      udotrad = - UdotInt(p);
#  elif defined(RADIATING)
      udotrad = UdotRad(p);
#  else
      udotrad = 0.0;
#  endif
      EntFunc(p) += dt * (UdotInt(p) + udotrad) *
	(gamma0-1) / rpow(Rho(p), gamma0-1);
      if (EntFunc(p) < 0.0)
	error("%s: negative entropy function; timestep = %g\n"
	      "  body %d:  entf[now] = %g  udot = %g,%g  freq = %g (%d)\n",
	      getprog(), dt, (int) (p - btab), EntFunc(p) - dt *
	      (UdotInt(p) + udotrad) * (gamma0-1) / rpow(Rho(p), gamma0-1),
	      UdotInt(p), udotrad, Freq(p), CurLevel(p));
#else
#  if defined(ISOTHERMAL)
      udotrad = - UdotInt(p);
#  elif defined(RADIATING)
      udotrad = UdotRad(p);
#  else
      udotrad = 0.0;
#  endif
      Uintern(p) += dt * (UdotInt(p) + udotrad);
      if (Uintern(p) < 0.0)
	error("%s: negative internal energy; timestep = %g\n"
	      "  body %d:  uint = %g  udot = %g,%g  freq = %g (%d)\n",
	      getprog(), dt, (int) (p - btab),
	      Uintern(p) - dt * (UdotInt(p) + udotrad),
	      UdotInt(p), udotrad, Freq(p), CurLevel(p));
#  if !defined(ISOTHERMAL)
      if (uintmax > 0 && Uintern(p) > uintmax) {
	udotrad += (uintmax - Uintern(p)) / dt;
	Uintern(p) = uintmax;
      }
#  endif 
#endif
      edotrad -= Mass(p) * udotrad;
    }
  }
  eradiate += dt * edotrad;
#if defined(STARFORM)
  starform(dt);					// form dt worth of stars
#endif
#if defined(MASSLOSS)
  massloss(dt);					// kill dt worth of stars
#endif
  tnow += dt;
  *fstep += 1.0 / (1 << levmax);
}

//  update: recompute derivatives for all bodies at current level or higher.
//  ________________________________________________________________________

local void update(int lev, real fstep)
{
  int nupdate, oldlev, newlev;
  bodyptr p;
  real dt1, dt2, x;
  vector vdif, aold;

  nupdate = 0;					// init count of updates
  for (p = btab; p < btab+nbody; p++) {		// loop over body table
    Flags(p) = INCLUDE | (CurLevel(p) < lev ? 0 : UPDATE);
						// flag bodies to update
    if (Update(p))				// is this one flagged?
      nupdate++;				// count updating bodies
  }
#if defined(TRACESTEP)
  fprintf(logstream, "    %12.5f%8d%8d\n", fstep, lev, nupdate);
  fflush(logstream);
#endif
  evaluate(fstep == 1.0);			// update flagged bodies
  levmax = lev;					// keep floor on levmax
  for (p = btab; p < btab+nbody; p++)		// loop over all bodies...
    if (Update(p)) {				// which need updating
      oldlev = CurLevel(p);			// remember previous level
      newlev = MAX(NewLevel(p), MAX(oldlev - 1, lev));
						// enforce limits on level
      CurLevel(p) = newlev;			// update level of body
      levmax = MAX(levmax, CurLevel(p));	// and max level of system
      dt1 = dtime / (2 << oldlev);		// use half previous step
      dt2 = dtime / (2 << newlev);		// plus half current step
      SUBV(vdif, Vel(p), Vmid(p));
      DIVVS(aold, vdif, dt1);			// recover previous accel.
      SETV(Vel(p), Vmid(p));
      ADDMULVS(Vel(p), Acc(p), 0.75 * dt1);
      ADDMULVS(Vel(p), aold, 0.25 * dt1);
      x = (3 + dt2 / dt1) / 4;			// set interpolation factor
      ADDMULVS(Vmid(p), Acc(p), x * (dt1 + dt2));
      ADDMULVS(Vmid(p), aold, (1 - x) * (dt1 + dt2));
    }
}

//  evaluate: calculate forces and time derivatives.
//  ________________________________________________

local void evaluate(bool report)
{
  bodyptr p;

  if (report)					// log output requested?
    outputhead();				// announce force calc.
  gravcalc(report);				// do gravity calculation
  sphcalc(report);				// do SPH calculation
  if (fdrag != 0.0)				// are drag forces enabled?
    for (p = btab; p < btab+nbody; p++)		// loop over all bodies
      if (Gas(p) && Update(p)) {		// pick updated gas bodies
	ADDMULVS(Acc(p), Vel(p), - fdrag);
						// subtract drag force
      }
}

//  gravcalc: calculate gravitational forces.
//  _________________________________________

local void gravcalc(bool report)
{
  bodyptr p;
  real r, mr3i;

#if defined(GRAVITY)
  maketree(btab, nbody);			// construct tree structure
  gravforce();					// compute self-gravity
  if (logstream != NULL && report)		// log output requested?
    gravreport(logstream, nbody);		// print gravity stats.
#elif defined(EXTGRAV)
  for (p = btab; p < btab+nbody; p++)		// loop over all bodies
    if (Update(p) && gravgsp != NULL) {		// update in extern field?
      r = absv(Pos(p));				// get distance from origin
      mr3i = - gsp_mass(gravgsp, r) / (r*r*r);
      MULVS(Acc(p), Pos(p), mr3i);		// set extern acc and phi
      Phi(p) = gsp_phi(gravgsp, r);
    } else if (Update(p)) {			// update in zero field?
      CLRV(Acc(p));				// zero grav. acceleration
      Phi(p) = 0.0;				// zero grav. potential
    }
#elif defined(NOACCEL)
  for (p = btab; p < btab+nbody; p++)
    if (Update(p)) {
      CLRV(Acc(p));
      Phi(p) = 0.0;
    }
#endif
}

//  sphcalc: calculate pressure forces and other SPH parameters.
//  ____________________________________________________________

local void sphcalc(bool report)
{
  kdxptr kd;
  smxptr sm;
  bodyptr p;

  kd = init_kdtree(btab, nbody, ngas);		// prepare for kd tree
  build_kdtree(kd, nbucket);			// do tree construction
  sm = init_smooth(kd, nsmooth, slope0);	// prepare for smoothing
  for (p = btab; p < btab+nbody; p++)		// loop over gas bodies
    if (Gas(p))
      Rho(p) = 0.0;				// prepare to sum density
  smooth(sm, sph_density);			// compute gas density
#if defined(DIFFUSING) || defined(OPAQUE)
  resmooth(sm, sph_depth);			// estimate optical depth
#endif
  for (p = btab; p < btab+nbody; p++)		// loop over all gas bodies
    if (Gas(p)) {
#if defined(ENTROPY)
      Press(p) = EntFunc(p) * rpow(Rho(p), gamma0);
						// set presure using a(s)
#else
      Press(p) = (gamma0 - 1) * Uintern(p) * Rho(p);
						// set presure using uint
#endif
      UdotInt(p) = 0.0;				// prepare to sum udot
#if defined(COMPVISC)
      UdotVis(p) = 0.0;
#endif
    }
  resmooth(sm, sph_derivs);			// compute time derivatives
  for (p = btab; p < btab+nbody; p++)		// loop over all gas bodies
    if (Gas(p)) {
#if defined(RADIATING)
#if defined(DIFFUSING) || defined(OPAQUE)
      UdotRad(p) = - (Tau(p) < 1.0 ? coolfunc(p) : 0.0);
#if defined(DIFFUSING)
      if (Surface(p))
	UdotRad(p) -= sigmastar * rsqr(rsqr(Uintern(p))) /
                (2 * Smooth(p) * Rho(p));
#endif
#else
      UdotRad(p) = - coolfunc(p);
#endif
      if (Update(p))
	Freq(p) = MAX(Freq(p), (ABS(UdotRad(p)) / Uintern(p)) / courant);
#endif	// RADIATING
#if defined(CONDUCTING) || (defined(RADIATING) && defined(DIFFUSING))
      if (Update(p))
	Freq(p) = MAX(Freq(p), (ABS(UdotInt(p)) / Uintern(p)) / courant);
#endif
    }
  setlevels(sm);				// set levels using freqs
  if (logstream != NULL && report)		// log output requested?
    sphreport(logstream, sm, options);		// print smoothing stats.
  finish_smooth(sm);				// deallocate smooth data
  finish_kdtree(kd);				// deallocate kdtree data
}

//  sph_density: sum smoothed density for body and its neighbors.
//  _____________________________________________________________

local void sph_density(smxptr sm, int pi, int nball)
{
  bodyptr bi = sm->kd->bptr[pi], bj;
  real hinv2, wsc, rhinv2, wsm;
  int j;

  hinv2 = 4 / sm->r2ball[pi];
  wsc = 0.5 * rsqrt(hinv2) * hinv2 / PI;
  for (j = 0; j < nball; ++j) {
    bj = sm->kd->bptr[sm->inlist[j]];
    rhinv2 = sm->r2list[j] * hinv2;
    WSmooth(wsm, wsc, rhinv2, sm->coefs);
    Rho(bi) += wsm * Mass(bj);
    Rho(bj) += wsm * Mass(bi);
  }
}

#if defined(DIFFUSING) || defined(OPAQUE)

//  sph_depth: estimate optical depth using gradient of density.
//  ____________________________________________________________

local void sph_depth(smxptr sm, int pi, int nball)
{
  bodyptr bi = sm->kd->bptr[pi], bj;
  real hinv2, dwsc, rhinv2, dwsm, dwsmrinv, fij;
  vector gradrho, rij, gradW;
  int j;

  hinv2 = 4 / sm->r2ball[pi];
  dwsc = hinv2 * hinv2 / PI;
  CLRV(gradrho);
  for (j = 0; j < nball; ++j) {
    bj = sm->kd->bptr[sm->inlist[j]];
    if (bi != bj) {
      rhinv2 = sm->r2list[j] * hinv2;
      dWSmooth(dwsm, dwsc, rhinv2, sm->coefs);
      dwsmrinv = dwsm / sqrt(sm->r2list[j]);
      SUBV(rij, Pos(bi), Pos(bj));
      fij = Rho(bj) / Rho(bi) - 1.0;
      ADDMULVS(gradrho, rij, Mass(bj) * fij * dwsmrinv);
    }
  }
  Tau(bi) = opacity * rsqr(Rho(bi)) / absv(gradrho);
}

#endif

//  sph_derivs: sum SPH accelerations and internal energy derivatives.
//  __________________________________________________________________

#define ETA2  0.01				// prevent divergence

#if defined(TRUE)
#define TauSWITCH(t)  \
  (t<0.8 ? 0 : t<1.0 ? 12.5*rsqr(t-0.8) : t<1.2 ? 1 - 12.5*rsqr(1.2-t) : 1)
#else
#define TauSWITCH(t)  \
  (t<1.0 ? 0 : 1)
#endif

local void sph_derivs(smxptr sm, int pi, int nball)
{
  bodyptr bi = sm->kd->bptr[pi], bj;
  real hinv2, dwsc, epsi, hsi, csi, divv, mumax, rhinv2, dwsm, dwsmrinv;
  real epsj, rv, hsavg, csavg, rhoavg, muvis, pivis, ascale, freqcou;
  real chii, chij, chiavg, *r2list = sm->r2list;
  int *inlist = sm->inlist, j, lev;
  vector rij, vij;

  hinv2 = 4 / sm->r2ball[pi];
  dwsc = 0.5 * hinv2 * hinv2 / PI;
  epsi = Press(bi) / (Rho(bi) * Rho(bi));
  hsi = rsqrt(sm->r2ball[pi]) / 2;
  csi = rsqrt(gamma0 * Press(bi) / Rho(bi));
  divv = mumax = 0.0;
  for (j = 0; j < nball; j++) {
    bj = sm->kd->bptr[inlist[j]];
    if (bi != bj) {
      rhinv2 = r2list[j] * hinv2;
      dWSmooth(dwsm, dwsc, rhinv2, sm->coefs);
      dwsmrinv = dwsm / rsqrt(r2list[j]);
      epsj = Press(bj) / (Rho(bj) * Rho(bj));
      SUBV(rij, Pos(bi), Pos(bj));
      SUBV(vij, Vel(bi), Vel(bj));
      DOTVP(rv, rij, vij);
      hsavg = (hsi + rsqrt(sm->r2ball[inlist[j]]) / 2) / 2;
      muvis = rv / (r2list[j] / hsavg + ETA2 * hsavg);
      mumax = MAX(mumax, ABS(muvis));
      if (rv < 0) {                       	// compute artificial visc.
	csavg = (csi + rsqrt(gamma0 * Press(bj) / Rho(bj))) / 2;
	rhoavg = (Rho(bi) + Rho(bj)) / 2;
	pivis = muvis * (beta * muvis - alpha * csavg) / rhoavg;
	UdotInt(bi) += Mass(bj) * 0.5 * pivis * rv * dwsmrinv;
	UdotInt(bj) += Mass(bi) * 0.5 * pivis * rv * dwsmrinv;
						// add shock heating term
#if defined(COMPVISC)
	UdotVis(bi) += Mass(bj) * 0.5 * pivis * rv * dwsmrinv;
	UdotVis(bj) += Mass(bi) * 0.5 * pivis * rv * dwsmrinv;
#endif
	ascale = (epsi + epsj + pivis) * dwsmrinv;
      } else
	ascale = (epsi + epsj) * dwsmrinv;
#if !defined(NOACCEL)
      if (Update(bi)) {
	ADDMULVS(Acc(bi), rij, - Mass(bj) * ascale);
      }
      if (Update(bj)) {
	ADDMULVS(Acc(bj), rij, Mass(bi) * ascale);
      }
#endif

#if !defined(ENTROPY)
      UdotInt(bi) += Mass(bj) * epsi * rv * dwsmrinv;
      UdotInt(bj) += Mass(bi) * epsj * rv * dwsmrinv;
						// include PdV work
#if defined(RADIATING) && defined(DIFFUSING)
      chii = TauSWITCH(Tau(bi)) * (16.0/3.0) *
	       (sigmastar / opacity) * rqbe(Uintern(bi)) / Rho(bi);
      chij = TauSWITCH(Tau(bj)) * (16.0/3.0) *
	       (sigmastar / opacity) * rqbe(Uintern(bj)) / Rho(bj);
      chiavg = 0.5 * (chii + chij);
      UdotInt(bi) += Mass(bj) * (chiavg / (Rho(bi) * Rho(bj))) *
	(Uintern(bi) - Uintern(bj)) * dwsm *
	rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
      UdotInt(bj) += Mass(bi) * (chiavg / (Rho(bi) * Rho(bj))) *
	(Uintern(bj) - Uintern(bi)) * dwsm *
	rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
      if ((Tau(bi) > 1.0) != (Tau(bj) > 1.0))	// neighbors straddle surface?
	SetFlag(bi, SURFACE);
#endif  // RADIATING && DIFFUSING
#if defined(CONDUCTING)
      UdotInt(bi) += Mass(bj) * (conduct / (Rho(bi) * Rho(bj))) *
	(Uintern(bi) - Uintern(bj)) * dwsm *
	rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
      UdotInt(bj) += Mass(bi) * (conduct / (Rho(bi) * Rho(bj))) *
	(Uintern(bj) - Uintern(bi)) * dwsm *
	rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
#endif  // CONDUCTING
#endif  // ENTROPY
      divv -= 2 * Mass(bj) * rv * dwsmrinv / Rho(bi);
    }
  }
  if (Update(bi))				// set hydro freq value
    Freq(bi) = ((csi + 1.2 * (alpha*csi + beta*mumax)) / hsi + ABS(divv)) /
                 courant;
}

//  coolfunc: evaluate simplified cooling function.
//  _______________________________________________

local real coolfunc(bodyptr p)
{
#if defined(RADIATING)
  real logx, f1, f2, f3, udotrad;

  logx = log10(Uintern(p) / uradpk);
  f1 = -0.17 * rsqr(logx + 1.23) - 3.14;
  f2 = -1.88 * rsqr(rsqr(logx));
  f3 = (logx > 0.97 ? 0.0 : -2.0 * rsqr(rsqr(logx - 0.97)))  - 1.60;
  udotrad = lambpk * Rho(p) * (rdex(f1) + rdex(f2) + rdex(f3));
  return (udotrad);
#else
  return (0.0);
#endif
}    

//  setlevels: set NewLevel values from computed Frequencies.
//  _________________________________________________________

local void setlevels(smxptr sm)
{
  bodyptr p;
  int lev, nmax;
  
  sm->freqmax = sm->freqsum = 0.0;
  for (p = btab; p < btab+nbody; p++) {
    if (Gas(p) && Update(p)) {
      sm->freqmax = MAX(sm->freqmax, Freq(p));
      sm->freqsum += Freq(p);
      lev = rceil(rlog2(MAX(dtime * Freq(p), 1.0)));
      if (lev > MAXLEV) {
	eprintf("[%s: courant condition violated; level = %d]\n",
		getprog(), lev);
	lev = MAXLEV;
      }
      NewLevel(p) = lev;
    }
    if (! Gas(p))
      NewLevel(p) = 0;
  }
  if (! scanopt(options, "nolimit")) {		// nolimit option not set?
    lev = 0;
    for (p = btab; p < btab+nbody; p++)
      if (NewLevel(p) > lev) {			// found higher level?
	lev = NewLevel(p);
	nmax = 1;
      } else if (NewLevel(p) == lev)
	nmax++;
    if (nmax < nsmooth/4)			// less than 25% have max?
      for (p = btab; p < btab+nbody; p++)
	if (NewLevel(p) == lev)			// demote them one level
	  NewLevel(p)--;
  }
  if (scanopt(options, "fixstep")) {		// fixstep option set?
    if (scanopt(options, "lockstep"))
      error("%s: fixstep and lockstep incompatable\n", getprog());
    for (p = btab; p < btab+nbody; p++)
      NewLevel(p) = 0;
  }
  if (scanopt(options, "lockstep")) {		// lockstep option set?
    lev = 0;					// place floor on level
    for (p = btab; p < btab+nbody; p++)
      lev = MAX(NewLevel(p), lev);		// determine maximum level
    for (p = btab; p < btab+nbody; p++)
      NewLevel(p) = lev;			// assign max level to all
  }
}

#if defined(STARFORM)

stream starlog = NULL;

//  starform: handle star formation as a random process.
//  ____________________________________________________

local void starform(real dt)
{
  real rho, udot, rate;

  for (bodyptr p = btab; p < btab+nbody; p++)
    if (Gas(p)) {				// select all gas bodies
      rho = Rho(p);
#if defined(COMPVISC)
      udot = UdotVis(p);			// viscosity udot is >= 0
#else	
      udot = MAX(UdotInt(p), 0.0);		// total udot may be < 0
#endif
      rate = 0.0;				// rate has units of 1/time
      for (int k = 0; k < 3; k++)
	if (cstar[k] > 0)			// sum rate for all rules
	  rate += cstar[k] * rpow(rho, nstar[k] - 1.0) * rpow(udot, mstar[k]);
      if (xrandom(0.0, 1.0) < rate * dt) {	// a star is born!
	if (starlog == NULL)
	  starlog = stropen(getparam("starlog"), "w!");
	fprintf(starlog, "%8d  born  %10.8f  %10.8f\n",
		(int) (p - btab), tnow, rate);
	Type(p) = (Type(p) & ~GAS) | STAR;	// switch gas to star
	ngas--;					// count one less gas body
	eradiate += Mass(p) * Uintern(p);	// dump internal energy
	Birth(p) = tnow;			// record birth date
#if defined(MASSLOSS)
	Death(p) = tnow + tau_ml * rpow(xrandom(0.0, 1.0), - beta_ml);
						// cut thread of lifetime
#endif
      }
    }
  if (starlog != NULL)
    fflush(starlog);
}

#endif

#if defined(STARFORM) && defined(MASSLOSS)

//  massloss: handle stellar mass loss as a random process.
//  _______________________________________________________

local void massloss(real dt)
{
  bodyptr p;

  for (p = btab; p < btab+nbody; p++)
    if (Star(p)) {				// select all star bodies
      if (Death(p) < tnow + dt) {		// star expires during step?
	if (starlog == NULL)
	  starlog = stropen(getparam("starlog"), "w!");
	fprintf(starlog, "%8d  expd  %10.8f\n", (int) (p - btab), Death(p));
	Type(p) = (Type(p) & ~STAR) | GAS;	// switch star to gas
	ngas++;					// count one more gas body
	eradiate -= Mass(p) * Uintern(p);	// restore internal energy
	Birth(p) = Death(p) = 0;		// zero out star dates
      }
    }
  if (starlog != NULL)
    fflush(starlog);
}

#endif
