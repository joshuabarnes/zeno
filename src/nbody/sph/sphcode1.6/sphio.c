/*
 * sphio.c: I/O routines for hierarchical SPH/N-body code.
 * Copyright (c) 2016 by Joshua E. Barnes, Honolulu, Hawai'i.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "strset.h"
#include "filestruct.h"
#include "gsp.h"
#include "sphcode.h"
#include "fixbody.h"

#include <sys/types.h>
#include <sys/stat.h>

//  Prototypes for local routines.
//  ______________________________

local void diagnostics(void);			// eval N-body diagnostics
local void writestream(void);			// dump body data to stream

//  Diagnostic parameters.
//  ______________________

local real mtot;		                // total mass of system
local double etot[4];				// system E, Ek, Ei, Ep
local matrix keten;				// kinetic energy tensor
local matrix peten;				// potential energy tensor
local vector cmpos;				// center of mass position
local vector cmvel;				// center of mass velocity
local vector amvec;				// angular momentum vector

//  inputdata: read initial conditions from input file.
//  ___________________________________________________

#define CHECKTAGS(tags, x) \
  if (! set_member(tags, x)) error("%s: %s data missing\n", getprog(), x)

void inputdata(void)
{
  stream instr;
  string intags[MaxBodyFields];

  instr = stropen(infile, "r");			// open input stream
  get_history(instr);				// read file history data
  btab = NULL;					// be sure to alloc space
  if (! get_snap(instr, &btab, &nbody, &tnow, intags, FALSE, NULL))
    error("%s.inputdata: no data in input file\n", getprog());
  strclose(instr);				// close input stream
  CHECKTAGS(intags, TypeTag);
  CHECKTAGS(intags, MassTag);
  CHECKTAGS(intags, PosTag);
  CHECKTAGS(intags, VelTag);
#if defined(ENTROPY)
  CHECKTAGS(intags, EntFuncTag);
#else
  CHECKTAGS(intags, UinternTag);
#endif
#if defined(MASSLOSS)
  CHECKTAGS(intags, BirthTag);
  CHECKTAGS(intags, DeathTag);
#endif
  ngas = 0;					// count gas bodies
  for (bodyptr p = btab; p < btab+nbody; p++) {	// check each body's type
    if (Cell(p) || ! Body(p))
      error("%s.inputdata: only body types allowed in input\n", getprog());
    if (Gas(p) && Star(p))
      error("%s.inpudata: gas & star types mutually exclusive\n", getprog());
    if (Gas(p)) {
      ngas++;
#if defined(ENTROPY)
      if (EntFunc(p) <= 0)
	error("%s.inputdata: entropy must be positive definite\n", getprog());
#else
      if (Uintern(p) <= 0)
	error("%s.inputdata: energy must be positive definite\n", getprog());
#endif
    }
    if (Star(p)) {
#if defined(MASSLOSS)
      if (Birth(p) > tnow)
	error("%s.inputdata: stars must already be born\n", getprog());
      if (Death(p) < tnow)
	error("%s.inputdata: stars must still be alive\n", getprog());
      if (Birth(p) >= Death(p))
	error("%s.inputdata: stars must be born before dying\n", getprog());
#endif
    }
  }
}

//  inputgrav: read GSP for external grav field, if defined.
//  ________________________________________________________

void inputgrav(string gspfile)
{
  stream gspstr;

#if defined(EXTGRAV)
  gspstr = stropen(gspfile, "r");
  get_history(gspstr);
  gravgsp = get_gsprof(gspstr);			// read external field GSP
  strclose(gspstr);
#endif
}

//  startoutput: begin output to log file and stream file.
//  __________________________________________________

void startoutput(string *opts)
{
  int ncol = 10;

  if (! strnull(logfile))			// explicit log file given?
    logstream = (strne(logfile, "/dev/null") ? stropen(logfile, "w!") : NULL);
  else						// use stdout unless busy
    logstream = (strne(strmpatn, "-") ? stdout : NULL);
  set_error_stream(logstream);			// set up error log

  if (logstream != NULL) {
    fprintf(logstream, "\n%s  [v%s]\n", headline, getversion());
    for (int i = 0; opts[i][0] == ';'; i++)
      fprintf(logstream, "    %s\n", opts[i]+1);
    fprintf(logstream, "\n %9s %9s %9s %9s %9s %9s %9s\n",
	    "nbody", "ngas", "nsmooth", "courant", "dtime", "dtout", "tstop");
    fprintf(logstream, " %9d %9d %9d %9.4f %9.6f %9.6f %9.4f\n\n",
	    nbody, ngas, nsmooth, courant, dtime, dtout, tstop);
    fprintf(logstream, " %9s", "gamma");
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
    fprintf(logstream, " %9s", "uintmax");
#endif
#if defined(RADIATING)
    fprintf(logstream, " %9s %9s", "uradpk", "lambpk");
#endif
#if defined(DIFFUSING)
    fprintf(logstream, " %9s", "sigmastar");
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    fprintf(logstream, " %9s", "opacity");
#endif
#if defined(CONDUCTING)
    fprintf(logstream, " %9s", "conduct");
#endif
#if defined(STARFORM)
    fprintf(logstream, " %9s %9s %9s", "cstar", "nstar", "mstar");
#if defined(MASSLOSS)
    fprintf(logstream, " %9s %9s", "tau_ml", "beta_ml");
#endif
#endif
    fprintf(logstream, "\n");
    fprintf(logstream, " %9.4f", gamma0);
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
    fprintf(logstream, " %9.4g", uintmax);
    ncol += 10;
#endif

#if defined(RADIATING)
    fprintf(logstream, " %9.4g %9.4g", uradpk, lambpk);
    ncol += 20;
#endif
#if defined(DIFFUSING)
    fprintf(logstream, " %9.4g", sigmastar);
    ncol += 10;
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    fprintf(logstream, " %9.4g", opacity);
    ncol += 10;
#endif
#if defined(CONDUCTING)
    fprintf(logstream, " %9.4g", conduct);
    ncol += 10;
#endif
#if defined(STARFORM)
    fprintf(logstream, " %9.4g %9.2f %9.2f", cstar[0], nstar[0], mstar[0]);
#if defined(MASSLOSS)
    fprintf(logstream, " %9.4g %9.2f", tau_ml, beta_ml);
#endif
#endif
    fprintf(logstream, "\n");
    for (int k = 1; k < 3; k++)
      if (cstar[k] > 0) {
	for (int i = 0; i < ncol; i++)
	  fputc(' ', logstream);
	fprintf(logstream, " %9.4g %9.2f %9.2f\n",
		cstar[k], nstar[k], mstar[k]);
      }
    fprintf(logstream, "\n");
#if defined(GRAVITY)
    fprintf(logstream, " %9s %9s %9s", "eps", "usequad", "theta");
#endif
    fprintf(logstream, " %9s %9s %9s %9s\n",
	    "alpha", "beta", "slope", "fdrag");
#if defined(GRAVITY)
    fprintf(logstream, " %9.4f %9s %9.2f",
	    eps, usequad ? "true" : "false", theta);
#endif
    fprintf(logstream, " %9.2f %9.2f %9.4f %9.4f\n",
	    alpha, beta, slope0, fdrag);
    if (! strnull(options))			// print options, if any
      fprintf(logstream, "\n\toptions: %s\n", options);
    fflush(logstream);
  }
  if (! strnull(savepatn))			// was state file given?
    savestate();				// save initial data
}

//  outputhead: announce beginning of new time-step.
//  ________________________________________________

void outputhead(void)
{
  if (logstream != NULL) {
    fprintf(logstream,
	    "\n____________________________________________________\n");
    fprintf(logstream, "time: %.8f    nstep: %d    cputime: %.2f\n",
	    tnow, nstep, cputime());
    fflush(logstream);
  }
}

//  outputdata: output diagnostics and binary data.
//  _______________________________________________

void outputdata(void)
{
  static string *outtags = NULL;
  int n;
  string reqtags[MaxBodyFields], *alltags;
  char outbuf[256];
  struct stat outstat;
  stream outstr;

  diagnostics();				// compute all diagnostics
  if (logstream != NULL) {
    fprintf(logstream, "\n%9s %9s %9s %9s %9s %9s %9s %9s\n",
	    "Time", "Etot", "Eint", "Ekin", "Epot", "Erad",
	    "|Jtot|", "|Vcom|");
    fprintf(logstream, "%9.4f %9.5f %9.5f %9.5f %9.5f %9.5f %9.6f %9.6f\n",
	    tnow, etot[0], etot[1], etot[2], etot[3], eradiate,
	    absv(amvec), absv(cmvel));
  }
  if (outtags == NULL)				// list requested data
    outtags = burststring(outputs, ", \t\n");
  if (! strnull(outpatn) && tnow + (dtime > 0 ? dtime/8 : 0) >= tout) {
    n = 0;					// build output list
    if (nstep == 0) {				// first output of run?
      reqtags[n++] = TypeTag;			// output type codes first
      reqtags[n++] = MassTag;			// output mass data next
#if defined(ENTROPY)
      reqtags[n++] = EntFuncTag;		// always output entropies
#else
      reqtags[n++] = UinternTag;		// always output energies
#endif
      reqtags[n++] = PosTag;			// always output positions
      reqtags[n++] = VelTag;			// always output velocities
    }
    reqtags[n] = NULL;				// terminate required list
    alltags = set_union(reqtags, outtags);	// list all output data
    sprintf(outbuf, outpatn, nstep);		// make up output file name
    if (stat(outbuf, &outstat) != 0) {		// no output file exists?
      outstr = stropen(outbuf, "w");		// create & open for output
      put_history(outstr);			// write out history data
    } else					// else file already exists
      outstr = stropen(outbuf, "a");		// reopen and append output
    put_snap(outstr, &btab, &nbody, &tnow, alltags);
    strclose(outstr);				// close up output file
    free(alltags);
    if (logstream != NULL)
      fprintf(logstream, "\n\tdata output to file %s at time %f\n",
	      outbuf, tnow);
    tout += dtout;				// schedule next output
  }
  if (logstream != NULL)
    fflush(logstream);
  if (! strnull(savepatn))			// was state file given?
    savestate();				// save data for restart
  writestream();
}

//  diagnostics: compute dynamic and thermodynamic diagnostics.
//  ___________________________________________________________

local void diagnostics(void)
{
  bodyptr p;
  real uint, velsq;
  vector tmpv;
  matrix tmpt;

  mtot = 0;					// zero total mass
  etot[1] = etot[2] = etot[3] = 0;		// zero energy totals
  CLRM(keten);					// zero ke tensor
  CLRM(peten);					// zero pe tensor
  CLRV(amvec);					// zero am vector
  CLRV(cmpos);					// zero c. of m. position
  CLRV(cmvel);					// zero c. of m. velocity
  for (p = btab; p < btab+nbody; p++) {		// loop over all particles
    mtot += Mass(p);				// sum particle masses
    if (Gas(p)) {
#if defined(ENTROPY)
      uint = EntFunc(p) * rpow(Rho(p), gamma0 - 1) / (gamma0 - 1);
#else
      uint = Uintern(p);
#endif
      etot[1] += Mass(p) * uint;		// sum internal energies
    }
    DOTVP(velsq, Vel(p), Vel(p));		// square vel vector
    etot[2] += 0.5 * Mass(p) * velsq;		// sum current KE
#if defined(GRAVITY)
    etot[3] += 0.5 * Mass(p) * Phi(p);		// and current PE
#elif defined(EXTGRAV)
    etot[3] += Mass(p) * Phi(p);		// and external PE
#endif
    MULVS(tmpv, Vel(p), 0.5 * Mass(p));		// sum 0.5 m v_i v_j
    OUTVP(tmpt, tmpv, Vel(p));
    ADDM(keten, keten, tmpt);
    MULVS(tmpv, Pos(p), Mass(p));		// sum m r_i a_j
    OUTVP(tmpt, tmpv, Acc(p));
    ADDM(peten, peten, tmpt);
    CROSSVP(tmpv, Vel(p), Pos(p));		// sum angular momentum
    MULVS(tmpv, tmpv, Mass(p));
    ADDV(amvec, amvec, tmpv);
    MULVS(tmpv, Pos(p), Mass(p));		// sum cm position
    ADDV(cmpos, cmpos, tmpv);
    MULVS(tmpv, Vel(p), Mass(p));		// sum cm momentum
    ADDV(cmvel, cmvel, tmpv);
  }
  etot[0] = etot[1] + etot[2] + etot[3] + eradiate;
						// sum total energy
  DIVVS(cmpos, cmpos, mtot);        		// normalize cm coords
  DIVVS(cmvel, cmvel, mtot);
}

//  writestream: output full particle data to stream file.
//  ______________________________________________________

local void writestream(void)
{
  char strmbuf[256];
  struct stat strmstat;
  stream strmstr;
  static string strmtags[] = {
    TypeTag, PosTag, VelTag, MassTag, SmoothTag, PhiTag, AccTag, RhoTag,
#if defined(ENTROPY)
    EntFuncTag,
#else
    UinternTag,
#endif
    UdotIntTag,
#if defined(RADIATING)
    UdotRadTag,
#endif
#if defined(COMPVISC)
    UdotVisTag,
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    TauTag,
#endif
#if defined(STARFORM) || defined(MASSLOSS)
    BirthTag,
#endif
#if defined(MASSLOSS)
    DeathTag,
#endif
    NULL,
  };
  if (! strnull(strmpatn)) {
    if (strne(strmpatn, "-")) {
      sprintf(strmbuf, strmpatn, nstep);
      if (stat(strmbuf, &strmstat) != 0) {
	strmstr = stropen(strmbuf, "w");
	put_history(strmstr);
      } else
	strmstr = stropen(strmbuf, "a");
    } else
      strmstr = stdout;
    put_snap(strmstr, &btab, &nbody, &tnow, strmtags);
    fflush(strmstr);
    if (strne(strmpatn, "-"))
      strclose(strmstr);
  }
}

//  savestate: write current state to disk file.
//  ____________________________________________

void savestate(void)
{
  char savebuf[256];
  stream str;
  static int randsize;
  static void *randstate = NULL;

  sprintf(savebuf, savepatn, nstep & 1);	// generate state file name
  str = stropen(savebuf, "w!");			// open state output file
  put_string(str, "program", getprog());
  put_string(str, "version", getversion());
  put_string(str, "headline", headline);
  put_data(str, "gamma", RealType, &gamma0, 0);	// save input parameters
#if defined(RADIATING)
  put_data(str, "uradpk", RealType, &uradpk, 0);
  put_data(str, "lambpk", RealType, &lambpk, 0);
#endif
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
  put_data(str, "uintmax", RealType, &uintmax, 0);
#endif
#if defined(CONDUCTING)
  put_data(str, "conduct", RealType, &conduct, 0);
#endif
#if defined(DIFFUSING)
  put_data(str, "sigmastar", RealType, &sigmastar, 0);
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
  put_data(str, "opacity", RealType, &opacity, 0);
#endif
#if defined(STARFORM)
  put_data(str, "cstar", RealType, cstar, 3, 0);
  put_data(str, "nstar", RealType, nstar, 3, 0);
  put_data(str, "mstar", RealType, mstar, 3, 0);
#if defined(MASSLOSS)
  put_data(str, "tau_ml", RealType, &tau_ml, 0);
  put_data(str, "beta_ml", RealType, &beta_ml, 0);
#endif
#endif
  put_data(str, "alpha", RealType, &alpha, 0);
  put_data(str, "beta", RealType, &beta, 0);
  put_data(str, "nsmooth", IntType, &nsmooth, 0);
  put_data(str, "nbucket", IntType, &nbucket, 0);
  put_data(str, "slope", RealType, &slope0, 0);
  put_data(str, "courant", RealType, &courant, 0);
  put_data(str, "dtime", RealType, &dtime, 0);
  put_data(str, "fdrag", RealType, &fdrag, 0);
#if defined(GRAVITY)
  put_data(str, "eps", RealType, &eps, 0);
  put_data(str, "usequad", BoolType, &usequad, 0);
  put_data(str, "theta", RealType, &theta, 0);
#elif defined(EXTGRAV)
  if (gravgsp != NULL)
    put_gsprof(str, gravgsp);
#endif
  put_string(str, "options", options);
  put_string(str, "outputs", outputs);
  put_data(str, "tstop", RealType, &tstop, 0);
  put_data(str, "dtout", RealType, &dtout, 0);
  put_data(str, "nstep", IntType, &nstep, 0);	// save current state vars
  put_data(str, "levmax", IntType, &levmax, 0);
  put_data(str, "rsize", RealType, &rsize, 0);
  put_data(str, "eradiate", RealType, &eradiate, 0);
  put_data(str, "tout", RealType, &tout, 0);
  put_data(str, "tnow", RealType, &tnow, 0);
#if defined(STARFORM)
  get_random_state(&randsize, &randstate);
  put_data(str, "randstate", ByteType, randstate, randsize, 0);
#endif
  put_data(str, "nbody", IntType, &nbody, 0);
  put_data(str, "ngas", IntType, &ngas, 0);
  put_data(str, "btab", AnyType, btab, nbody, sizeof(body), 0);
  strclose(str);
}

//  restorestate: restore state from disk file.
//  ___________________________________________

void restorestate(void)
{
  stream str;
  string program, version;
  int randsize;
  void *randstate;

  str = stropen(restfile, "r");			// open state input file
  program = get_string(str, "program");
  version = get_string(str, "version");
  if (strne(program, getprog()) || strne(version, getversion()))
    eprintf("[%s: warning: input state file may be outdated]\n",
	    getprog());
  headline = get_string(str, "headline");
  get_data(str, "gamma", RealType, &gamma0, 0);	// read program parameters
#if defined(RADIATING)
  get_data(str, "uradpk", RealType, &uradpk, 0);
  get_data(str, "lambpk", RealType, &lambpk, 0);
#endif
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
  if (get_tag_ok(str, "uintmax"))
    get_data(str, "uintmax", RealType, &uintmax, 0);
  else
    uintmax = 0;
#endif
#if defined(CONDUCTING)
  get_data(str, "conduct", RealType, &conduct, 0);
#endif
#if defined(DIFFUSING)
  get_data(str, "sigmastar", RealType, &sigmastar, 0);
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
  get_data(str, "opacity", RealType, &opacity, 0);
#endif
#if defined(STARFORM)
  get_data(str, "cstar", RealType, cstar, 3, 0);
  get_data(str, "nstar", RealType, nstar, 3, 0);
  get_data(str, "mstar", RealType, mstar, 3, 0);
#if defined(MASSLOSS)
  get_data(str, "tau_ml", RealType, &tau_ml, 0);
  get_data(str, "beta_ml", RealType, &beta_ml, 0);
#endif
#endif
  get_data(str, "alpha", RealType, &alpha, 0);
  get_data(str, "beta", RealType, &beta, 0);
  get_data(str, "nsmooth", IntType, &nsmooth, 0);
  get_data(str, "nbucket", IntType, &nbucket, 0);
  get_data(str, "slope", RealType, &slope0, 0);
  get_data(str, "courant", RealType, &courant, 0);
  get_data(str, "dtime", RealType, &dtime, 0);
  get_data(str, "fdrag", RealType, &fdrag, 0);
#if defined(GRAVITY)
  get_data(str, "eps", RealType, &eps, 0);
  get_data(str, "usequad", BoolType, &usequad, 0);
  get_data(str, "theta", RealType, &theta, 0);
#elif defined(EXTGRAV)
  if (get_tag_ok(str, "GeneralSphericalProfile"))
    gravgsp = get_gsprof(str);			// set external field GSP
  else
    gravgsp = NULL;				// else no external field
#endif
  options = get_string(str, "options");
  outputs = get_string(str, "outputs");
  get_data(str, "tstop", RealType, &tstop, 0);
  get_data(str, "dtout", RealType, &dtout, 0);
  get_data(str, "nstep", IntType, &nstep, 0);	// read current state vars
  get_data(str, "levmax", IntType, &levmax, 0);
  get_data(str, "rsize", RealType, &rsize, 0);
  get_data(str, "eradiate", RealType, &eradiate, 0);
  get_data(str, "tout", RealType, &tout, 0);
  get_data(str, "tnow", RealType, &tnow, 0);
#if defined(STARFORM)
  randsize = get_length(str, "randstate");
  randstate = allocate(randsize);
  get_data(str, "randstate", ByteType, randstate, randsize, 0);
  set_random_state(&randsize, &randstate);
  free(randstate);
#endif
  get_data(str, "nbody", IntType, &nbody, 0);
  get_data(str, "ngas", IntType, &ngas, 0);
  btab = (bodyptr) allocate(nbody * sizeof(body));
  get_data(str, "btab", AnyType, btab, nbody, sizeof(body), 0);
  strclose(str);
}
