/*
 * treeio.c: I/O routines for hierarchical N-body code.
 * Public routines: inputdata(), startoutput(), forcereport(),
 * output(), savestate(), restorestate().
 * Copyright (c) 2015 by Joshua E. Barnes, Kyoto, Japan.
 */

#include "treecode.h"
#include "filestruct.h"
#include <sys/types.h>
#include <sys/stat.h>

//  Prototypes for local routines and variables.
//  ____________________________________________

local void diagnostics(void);			// eval N-body diagnostics

local stream logstr;				// log stream, or stdout
local bool forcehead;				// force headers printed?
local real mtot;		                // total mass of system
local real etot[3];				// Etot, KE, PE of system
local matrix keten;				// kinetic energy tensor
local matrix peten;				// potential energy tensor
local vector cmpos;				// center of mass position
local vector cmvel;				// center of mass velocity
local vector amvec;				// angular momentum vector

//  inputdata: read initial conditions from input file.
//  ___________________________________________________

void inputdata(void)
{
  stream instr;
  string intags[MaxBodyFields];
  bodyptr p;

  bodytab = NULL;				// request new input data
  instr = stropen(infile, "r");			// open input stream
  get_history(instr);				// read file history data
  if (! get_snap(instr, &bodytab, &nbody, &tnow, intags, FALSE, NULL))
    error("%s.inputdata: no data in input file\n", getprog());
  strclose(instr);				// close input stream
  if (! set_member(intags, MassTag) || ! set_member(intags, PosTag) ||
      ! set_member(intags, VelTag))
    error("%s.inputdata: essential data missing\n", getprog());
  if (scanopt(options, "reset-time"))		// reset starting time?
    tnow = 0.0;
  for (p = bodytab; p < bodytab+nbody; p++)	// loop over new bodies
    Type(p) = BODY;				// initializing body type
}

//  startoutput: begin output to log file.
//  ______________________________________

void startoutput(void)
{
  if (! strnull(logfile))
    logstr = stropen(logfile, "w!");		// always open file if given
  else
    logstr = (streq(strfile, "-") ? NULL : stdout);	// don't hog stdout
  if (logstr != NULL) {
    fprintf(logstr, "\n%s\n", headline);	// print headline, params
    fprintf(logstr, "\n%9s%9s%9s%9s", "nbody", "dtime", "nstatic", "eps");
#if !defined(QUICKSCAN)
    fprintf(logstr, "%9s", "theta");
#endif
    fprintf(logstr, "%9s%9s%9s\n", "usequad", "dtout", "tstop");
    fprintf(logstr, "%9d%9.5f%9d%9.4f", nbody, dtime, nstatic, eps);
#if !defined(QUICKSCAN)
    fprintf(logstr, "%9.2f", theta);
#endif
    fprintf(logstr, "%9s%9.5f%9.4f\n",
	    usequad ? "true" : "false", dtout, tstop);
    if (! strnull(options))			// print options, if any
      fprintf(logstr, "\n\toptions: %s\n", options);
  }
  if (! strnull(savefile))			// was state file given?
    savestate(savefile);			// save initial data
  forcehead = FALSE;				// require header printing
  fflush(NULL);					// empty all output buffers
}

//  forcereport: print staristics on tree construction and force calculation.
//  _________________________________________________________________________

void forcereport(void)
{
  if (logstr != NULL) {
    if (! forcehead)				// no force header printed?
      fprintf(logstr, "\n    %8s%8s%8s%8s%12s%12s%8s\n",
	      "rsize", "tdepth", "ftree",
	      "nfcalc", "nbbtot", "nbctot", "CPUfc");
    fprintf(logstr, "    %8.1f%8d%8.3f%8d%12ld%12ld%8.3f\n",
	    rsize, tdepth, (nbody + ncell - 1) / ((real) ncell),
	    nfcalc, nbbcalc, nbccalc, cpuforce);
    forcehead = TRUE;
    fflush(NULL);				// empty all output buffers
  }
}

//  output: compute diagnostics and output binary data.
//  ___________________________________________________

void output(void)
{
  real teff;
  int n;
  string outtags[MaxBodyFields];
  char namebuf[256];
  struct stat buf;
  stream outstr;
  string streamtags[] = { MassTag, PosTag, VelTag, PhiTag, AccTag, NULL };

  diagnostics();				// compute std diagnostics
  if (logstr != NULL) {
    fprintf(logstr, "\n    %8s%8s%8s%8s%8s%8s%8s%8s\n",
	    "time", "|T+U|", "T", "-U", "-T/U", "|Vcom|", "|Jtot|", "CPUtot");
    fprintf(logstr, "    %8.3f%8.4f%8.4f%8.4f%8.5f%8.5f%8.4f%8.2f\n",
	    tnow, ABS(etot[0]), etot[1], -etot[2], -etot[1]/etot[2],
	    absv(cmvel), absv(amvec), cputime());
  }
  teff = tnow + dtime/8;			// anticipate slightly...
  if (! strnull(outfile) && teff >= tout) {	// time for data output?
    n = 0;
    if (scanopt(outputs, PosTag))		// if listed in outputs
      outtags[n++] = PosTag;			// include tag in list
    if (scanopt(outputs, VelTag))
      outtags[n++] = VelTag;
    if (scanopt(outputs, MassTag) || (nstep == 0))
      outtags[n++] = MassTag;
    if (scanopt(outputs, PhiTag))
      outtags[n++] = PhiTag;			// select potential data
    if (scanopt(outputs, AccTag))
      outtags[n++] = AccTag;			// select acceleration data
#ifdef TESTBED
    outtags[n++] = KeyTag;
    outtags[n++] = AuxTag;
#endif
    outtags[n] = NULL;
    sprintf(namebuf, outfile, nstep);		// make up output file name
    if (stat(namebuf, &buf) != 0) {		// no output file exists?
      outstr = stropen(namebuf, "w");		// create & open for output
      put_history(outstr);			// write out hiatory data
    } else					// else file already exists
      outstr = stropen(namebuf, "a");		// reopen and append output
    put_snap(outstr, &bodytab, &nbody, &tnow, outtags);
    strclose(outstr);				// close up output file
    if (logstr != NULL)
      fprintf(logstr, "\n\tdata output to file %s at time %f\n",
	      namebuf, tnow);
    tout += dtout;				// schedule next output
  }
  if (! strnull(strfile)) {
    if (! streq(strfile, "-")) {
      sprintf(namebuf, strfile, nstep);
      if (stat(namebuf, &buf) != 0) {
	outstr = stropen(namebuf, "w");
	put_history(outstr);
      } else
	outstr = stropen(namebuf, "a");
    } else
      outstr = stdout;
    put_snap(outstr, &bodytab, &nbody, &tnow, streamtags);
    fflush(outstr);
    if (! streq(strfile, "-"))
      strclose(outstr);
  }
  if (! strnull(savefile))			// was state file given?
    savestate(savefile);			// save data for restart
  forcehead = FALSE;				// insure headers print
  fflush(NULL);					// empty all output buffers
}

//  diagnostics: compute set of dynamical diagnostics.
//  __________________________________________________

local void diagnostics(void)
{
  bodyptr p1, p2, p;
  real mp, velsq;
  vector tmpv;
  matrix tmpt;

  mtot = 0.0;					// zero total mass
  etot[1] = etot[2] = 0.0;			// zero total KE and PE
  CLRM(keten);					// zero ke tensor
  CLRM(peten);					// zero pe tensor
  CLRV(amvec);					// zero am vector
  CLRV(cmpos);					// zero c. of m. position
  CLRV(cmvel);					// zero c. of m. velocity
  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);
  for (p = p1; p < p2; p++) {			// loop over body range	
    mp = (testcalc ? 1.0 / (nbody - ABS(nstatic)) : Mass(p));
						// use eq. mass in testcalc
    mtot += mp;					// sum particle masses
    DOTVP(velsq, Vel(p), Vel(p));		// square vel vector
    etot[1] += 0.5 * mp * velsq;		// sum current KE
    etot[2] += (testcalc ? 1.0 : 0.5) * mp * Phi(p);
						// and PE, weighted right
    MULVS(tmpv, Vel(p), 0.5 * mp);		// sum 0.5 m v_i v_j
    OUTVP(tmpt, tmpv, Vel(p));
    ADDM(keten, keten, tmpt);
    MULVS(tmpv, Pos(p), mp);			// sum m r_i a_j
    OUTVP(tmpt, tmpv, Acc(p));
    ADDM(peten, peten, tmpt);
    CROSSVP(tmpv, Vel(p), Pos(p));		// sum angular momentum
    MULVS(tmpv, tmpv, mp);
    ADDV(amvec, amvec, tmpv);
    MULVS(tmpv, Pos(p), mp);			// sum cm position
    ADDV(cmpos, cmpos, tmpv);
    MULVS(tmpv, Vel(p), mp);			// sum cm momentum
    ADDV(cmvel, cmvel, tmpv);
  }
  etot[0] = etot[1] + etot[2];			// sum KE and PE
  DIVVS(cmpos, cmpos, mtot);        		// normalize cm coords
  DIVVS(cmvel, cmvel, mtot);
}

//  savestate: write current state to external file.
//  ________________________________________________

void savestate(string pattern)
{
  char namebuf[256];
  stream str;

  sprintf(namebuf, pattern, nstep & 1);		// construct alternate name
  str = stropen(namebuf, "w!");			// open state output file
  put_string(str, "program", getprog());
  put_string(str, "version", getversion());
  put_string(str, "headline", headline);	// save control parameters
  put_data(str, "dtime", RealType, &dtime, 0);
  put_data(str, "nstatic", IntType, &nstatic, 0);
#if !defined(QUICKSCAN)
  put_data(str, "theta", RealType, &theta, 0);
#endif
  put_data(str, "usequad", BoolType, &usequad, 0);
  put_data(str, "eps", RealType, &eps, 0);
  put_string(str, "options", options);
  put_string(str, "outputs", outputs);
  put_data(str, "tstop", RealType, &tstop, 0);
  put_data(str, "dtout", RealType, &dtout, 0);
  put_data(str, "tnow", RealType, &tnow, 0);	// save state variables
  put_data(str, "tout", RealType, &tout, 0);
  put_data(str, "nstep", IntType, &nstep, 0);
  put_data(str, "rsize", RealType, &rsize, 0);
  put_data(str, "nbody", IntType, &nbody, 0);
  put_data(str, "bodytab", AnyType, bodytab, nbody, sizeof(body), 0);
  strclose(str);
}

//  restorestate: restore state from external file.
//  _______________________________________________

void restorestate(string file)
{
  stream str;
  string program, version;

  str = stropen(file, "r");			// open state input file
  program = get_string(str, "program");
  version = get_string(str, "version");
  if (! streq(program, getprog()) ||		// check program, version
      ! streq(version, getversion()))
    eprintf("[%s: warning: state file may be outdated]\n", getprog());
  headline = get_string(str, "headline");	// read control parameters
  get_data(str, "dtime", RealType, &dtime, 0);
  get_data(str, "nstatic", IntType, &nstatic, 0);
#if !defined(QUICKSCAN)
  get_data(str, "theta", RealType, &theta, 0);
#endif
  get_data(str, "usequad", BoolType, &usequad, 0);
  get_data(str, "eps", RealType, &eps, 0);
  options = get_string(str, "options");
  outputs = get_string(str, "outputs");
  get_data(str, "tstop", RealType, &tstop, 0);
  get_data(str, "dtout", RealType, &dtout, 0);
  get_data(str, "tnow", RealType, &tnow, 0);	// read state variables
  get_data(str, "tout", RealType, &tout, 0);
  get_data(str, "nstep", IntType, &nstep, 0);
  get_data(str, "rsize", RealType, &rsize, 0);
  get_data(str, "nbody", IntType, &nbody, 0);
  bodytab = (bodyptr) allocate(nbody * sizeof(body));
  get_data(str, "bodytab", AnyType, bodytab, nbody, sizeof(body), 0);
  strclose(str);
}
