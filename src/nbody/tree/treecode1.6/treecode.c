/*
 * treecode.c: hierarchical N-body simulation code.
 * Copyright (c) 2017 by Joshua E. Barnes, Honolulu, Hawaii.
 */

#define global				// don't default global to extern

#include "treecode.h"

//  Default values for input parameters.
//  ____________________________________

string defv[] = {
#if defined(QUICKSCAN)
				";Hierarchical N-body code (quick scan)",
#elif defined(EXTGRAV)
				";Hierarchical N-body code (+ ext grav)",
#else
				";Hierarchical N-body code"
#  if defined(NOSOFTCORR)
                                " (no soft corr)",
#  else
                                " (soft corr)",
#  endif
#endif
    "in=",			";Input file with initial conditions",
    "out=",			";Output file pattern for N-body frames.",
				";Directive (eg, %d) formats step number.",
    "restore=",			";Continue run from existing state file",
    "save=",			";State file pattern, written each step.",
				";Directive (eg, %d) alternates 0 and 1.",
    "dtime=1/32",		";Timestep for leapfrog integration",
    "eps=0.025",		";Smoothing length for force calculation",
#if !defined(QUICKSCAN)
    "theta=1.0",		";Hierarchical force accuracy parameter",
#endif
    "usequad=true",		";If false, don't use quadrapole moments",
    "options=",			";Comma-separated list of program options.",
				";new-tout: reschedule output times.",
				";reset-time: set time to zero.",
				";bh86,sw94,theta-eff: opening criteria.",
    "outputs=" PosTag "," VelTag,
				";Body data arrays written to output.",
				";Other choices: "
				  MassTag "," PhiTag "," AccTag ".",
    "tstop=2.0",		";Time to stop integration",
    "dtout=1/4",		";Time interval between data outputs",
    "nstatic=0",		";Number of static bodies in array.",
				";If pos (neg), count from start (end).",
    "nbody=65536",		";Number of bodies generated for test run.",
				";If no input file given, make Plummer model.",
    "seed=123",			";Random number seed for test run",
    "stream=",			";Output file pattern for frame stream",
    "log=",			";Output file name for calculation log.",
				";Defaults to stdout unless stream=-.",
#if defined(EXTGRAV)
    "gravgsp=",			";Input GSP for external gravity field",
#endif
    "VERSION=1.6",		";Joshua Barnes  3 May 2017",
    NULL,
};

//  Prototypes for local procedures.
//  ________________________________

local void startrun(void);			// initialize system state
local void newrun(void);			// start new simulation
local void oldrun(void);			// resume old simulation
local void testdata(void);			// generate test data
local void stepsystem(void);			// advance by one time-step
local void treeforce(void);			// do force calculation

#if defined(EXTGRAV)
#include "gsp.h"
gsprof *gravgsp = NULL;				// GSP for external field
#endif

//  main: toplevel routine for hierarchical N-body code.
//  ____________________________________________________

int main(int argc, string argv[])
{
  initparam(argv, defv);			// initialize param access
  headline = defv[0] + 1;			// use default headline
  startrun();					// get params & input data
  startoutput();				// activate output code
  if (nstep == 0) {				// if data just initialized
    treeforce();				// calculate initial forces
    output();					// generate initial output
  }
  if (dtime != 0.0)				// if time steps requested
    while (tstop - tnow > 0.01 * dtime) {	// while not past tstop
      stepsystem();				// advance step by step
      output();					// output results each time
    }
  return (0);					// end with proper status
}

//  startrun: startup hierarchical N-body code.
//  ___________________________________________

local void startrun(void)
{
  bodyptr p1, p2, p;
  stream gravstr;

  define_body(sizeof(body), Precision, NDIM);	// setup phat body struct
  define_body_offset(PosTag,  BodyOffset(Pos));
  define_body_offset(VelTag,  BodyOffset(Vel));
  define_body_offset(MassTag, BodyOffset(Mass));
  define_body_offset(PhiTag,  BodyOffset(Phi));
  define_body_offset(AccTag,  BodyOffset(Acc));
#ifdef TESTBED
  define_body_offset(KeyTag,  BodyOffset(Key));
  define_body_offset(AuxTag,  BodyOffset(Aux));
#endif
  infile = getparam("in");			// set I/O file names
  outfile = getparam("out");
  savefile = getparam("save");
  if (strnull(getparam("restore")))		// starting a new run?
    newrun();
  else						// else resume old run
    oldrun();
  if (ABS(nstatic) > nbody)			// check nstatic is OK
    error("%s: absurd value for nstatic\n", getargv0());
  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);
  testcalc = TRUE;				// determine type of calc:
  for (p = p1; p < p2; p++)
    testcalc = testcalc && (Mass(p) == 0);	// look for dynamic masses
  strfile = getparam("stream");
  logfile = getparam("log");
#if defined(EXTGRAV)
  if (! strnull(getparam("gravgsp"))) {		// was GSP file given?
    gravstr = stropen(getparam("gravgsp"), "r");
    get_history(gravstr);
    gravgsp = get_gsprof(gravstr);		// read external field GSP
    strclose(gravstr);
  }
#endif
}

//  newrun, oldrun: initialize parameters and bodies.
//  _________________________________________________

local void newrun(void)
{
  eps = getdparam("eps");			// get input parameters
  dtime = getdparam("dtime");
  nstatic = getiparam("nstatic");
#if !defined(QUICKSCAN)
  theta = getdparam("theta");
#endif
  usequad = getbparam("usequad");
  tstop = getdparam("tstop");
  dtout = getdparam("dtout");
  options = getparam("options");
  outputs = getparam("outputs");
  if (! strnull(infile))			// if data file was given
    inputdata();				// then read inital data
  else {					// else make initial data
    nbody = getiparam("nbody");			// get number of bodies
    init_random(getiparam("seed"));		// set random number gen.
    testdata();					// and make plummer model
  }
  rsize = 1.0;					// start root w/ unit cube
  nstep = 0;					// begin counting steps
  tout = tnow;					// schedule first output
}

local void oldrun(void)
{
  restorestate(getparam("restore"));		// read in old state file
  if (getparamstat("eps") & ARGPARAM)		// was eps given new value?
    eps = getdparam("eps");			// use command line value
  if (getparamstat("nstatic") & ARGPARAM)	// likewise for others...
    nstatic = getiparam("nstatic");
#if !defined(QUICKSCAN)
  if (getparamstat("theta") & ARGPARAM)
    theta = getdparam("theta");
#endif
  if (getparamstat("usequad") & ARGPARAM)
    usequad = getbparam("usequad");
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

//  testdata: generate Plummer model initial conditions for test runs,
//  scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
//  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
//  __________________________________________________________________

#define MFRAC  0.999				// cut off 1-MFRAC of mass

local void testdata(void)
{
  real rsc, vsc, r, v, x, y;
  bodyptr p;

  if (nbody < 1)				// check for silly values
    error("%s: absurd value for nbody\n", getargv0());
  bodytab = (bodyptr) allocate(nbody * sizeof(body));
						// alloc space for bodies
  rsc = (3 * PI) / 16;				// set length scale factor
  vsc = rsqrt(1.0 / rsc);			// find speed scale factor
  for (p = bodytab; p < bodytab+nbody; p++) {	// loop over particles
    Type(p) = BODY;				// tag as a body
    Mass(p) = 1.0 / nbody;			// set masses equal
    x = xrandom(0.0, MFRAC);			// pick enclosed mass
    r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1);	// find enclosing radius
    pickshell(Pos(p), NDIM, rsc * r);		// pick position vector
    do {					// select from fn g(x)
      x = xrandom(0.0, 1.0);			// for x in range 0:1
      y = xrandom(0.0, 0.1);			// max of g(x) is 0.092
    } while (y > x*x * rpow(1 - x*x, 3.5));	// using von Neumann tech
    v = x * rsqrt(2.0 / rsqrt(1 + r*r));	// find resulting speed
    pickshell(Vel(p), NDIM, vsc * v);		// pick velocity vector
  }
  tnow = 0.0;					// set elapsed model time
}

//  stepsystem: advance N-body system using simple leap-frog.
//  _________________________________________________________

local void stepsystem(void)
{
  bodyptr p1, p2, p;

  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);
  for (p = p1; p < p2; p++) {			// loop over body range	
    ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);	// advance v by 1/2 step
    ADDMULVS(Pos(p), Vel(p), dtime);		// advance r by 1 step
  }
  treeforce();
  for (p = p1; p < p2; p++) {			// loop over body range	
    ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);	// advance v by 1/2 step
  }
  nstep++;					// count another time step
  tnow = tnow + dtime;				// finally, advance time
}

//  treeforce: supervise force calculation.
//  _______________________________________

local void treeforce(void)
{
  bodyptr p1, p2, p;
  real r, mr3i;

  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);
  for (p = bodytab; p < bodytab+nbody; p++)	// loop over all bodies
    Update(p) = (testcalc ? p1 <= p && p < p2 : TRUE);
						// flag bodies to update
  maketree(bodytab, nbody);			// construct tree structure
  gravcalc();					// compute current forces
  forcereport();				// print force statistics
#if defined(EXTGRAV)
  for (p = bodytab; p < bodytab+nbody; p++)	// loop over all bodies
    if (Update(p) && gravgsp != NULL) {		// update in extern field?
      r = absv(Pos(p));				// get distance from origin
      mr3i = - mass_gsp(gravgsp, r) / (r*r*r);
      ADDMULVS(Acc(p), Pos(p), mr3i);		// add extern acc and phi
      Phi(p) += phi_gsp(gravgsp, r);
    }
#endif
}
