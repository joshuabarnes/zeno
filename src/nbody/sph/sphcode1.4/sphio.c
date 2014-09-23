/*
 * sphio.c: I/O routines for hierarchical SPH/N-body code.
 * Copyright (c) 2012 by Joshua E. Barnes, Honolulu, Hawai'i.
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

// Prototypes for local routines.

local void diagnostics(void);			// eval N-body diagnostics
local void writetrace(void);			// dump body data to trace

// Output state variables.

local real mtot;		                // total mass of system
local double etot[4];				// system E, Ek, Ei, Ep
local matrix keten;				// kinetic energy tensor
local matrix peten;				// potential energy tensor
local vector cmpos;				// center of mass position
local vector cmvel;				// center of mass velocity
local vector amvec;				// angular momentum vector

// INPUTDATA: read initial conditions from input file.

void inputdata(void)
{
  stream istr;
  string intags[MaxBodyFields];
  bodyptr p;

  istr = stropen(infile, "r");			// open input stream
  get_history(istr);				// read file history data
  btab = NULL;					// be sure to alloc space
  if (! get_snap(istr, &btab, &nbody, &tnow, intags, FALSE))
    error("%s: no data in input file\n", getargv0());
  strclose(istr);				// close input stream
  if (! set_member(intags, TypeTag))
    error("%s: type data missing\n", getargv0());
  if (! set_member(intags, MassTag))
    error("%s: mass data missing\n", getargv0());
  if (! set_member(intags, PosTag))
    error("%s: position data missing\n", getargv0());
  if (! set_member(intags, VelTag))
    error("%s: velocity data missing\n", getargv0());
#if defined(ENTROPY)
  if (! set_member(intags, EntFuncTag))
    error("%s: entropy function data missing\n", getargv0());
#else
  if (! set_member(intags, UinternTag))
    error("%s: internal energy data missing\n", getargv0());
#endif
#if defined(MASSLOSS)
  if (! set_member(intags, BirthTag))
    error("%s: birth data missing\n", getargv0());
  if (! set_member(intags, DeathTag))
    error("%s: death data missing\n", getargv0());
#endif
  ngas = 0;					// count gas bodies
  for (p = btab; p < btab+nbody; p++) {		// check each body's type
    if (Cell(p) || ! Body(p))
      error("%s: input must contain bodies, not cells\n", getargv0());
    if (Gas(p) && Star(p))
      error("%s: gas and star types mutually exclusive\n", getargv0());
    if (Gas(p)) {
      ngas++;
#if defined(ENTROPY)
      if (EntFunc(p) <= 0)
	error("%s: entropy function must be positive definite\n", getargv0());
#else
      if (Uintern(p) <= 0)
	error("%s: internal energy must be positive definite\n", getargv0());
#endif
    }
    if (Star(p)) {
#if defined(MASSLOSS)
      if (Birth(p) > tnow)
	error("%s: stars must be born before run starts\n", getargv0());
      if (Death(p) < tnow)
	error("%s: stars must be alive when run starts\n", getargv0());
      if (Birth(p) >= Death(p))
	error("%s: stars must be born before they die\n", getargv0());
#endif
    }
  }
}

// INPUTGRAV: read GSP for external grav field, if defined.

void inputgrav(string gspfile)
{
  stream gstr;

#if defined(EXTGRAV)
  gstr = stropen(gspfile, "r");
  get_history(gstr);
  gravgsp = get_gsprof(gstr);			// read external field GSP
  strclose(gstr);
#endif
}

// STARTOUT: begin output to log file and trace file.

void startout(string defv[])
{
  int i;

  headline = defv[0] + 1;			// set ident. message
  fprintf(logstr, "\n%s  [v%s]\n", headline, getversion());
  for (i = 1; defv[i][0] == ';'; i++)
    fprintf(logstr, "    %s\n", defv[i]+1);
  fprintf(logstr, "\n %9s %9s %9s %9s %9s %9s %9s\n",
	  "nbody", "ngas", "nsmooth", "courant", "dtime", "dtout", "tstop");
  fprintf(logstr, " %9d %9d %9d %9.4f %9.6f %9.6f %9.4f\n\n",
	  nbody, ngas, nsmooth, courant, dtime, dtout, tstop);
  fprintf(logstr, " %9s", "gamma");
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
  fprintf(logstr, " %9s", "uintmax");
#endif
#if defined(RADIATING)
  fprintf(logstr, " %9s %9s", "uradpk", "lambpk");
#endif
#if defined(DIFFUSING)
  fprintf(logstr, " %9s", "sigmastar");
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
  fprintf(logstr, " %9s", "opacity");
#endif
#if defined(CONDUCTING)
  fprintf(logstr, " %9s", "conduct");
#endif
#if defined(STARFORM)
  fprintf(logstr, " %9s %9s %9s", "starprob", "rhoindx", "udotindx");
#if defined(MASSLOSS)
  fprintf(logstr, " %9s %9s", "tau_ml", "beta_ml");
#endif
#endif
  fprintf(logstr, "\n");
  fprintf(logstr, " %9.4f", gamma0);
#if !defined(ENTROPY) && !defined(ISOTHERMAL)
  fprintf(logstr, " %9.4g", uintmax);
#endif
#if defined(RADIATING)
  fprintf(logstr, " %9.4g %9.4g", uradpk, lambpk);
#endif
#if defined(DIFFUSING)
  fprintf(logstr, " %9.4g", sigmastar);
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
  fprintf(logstr, " %9.4g", opacity);
#endif
#if defined(CONDUCTING)
  fprintf(logstr, " %9.4g", conduct);
#endif
#if defined(STARFORM)
  fprintf(logstr, " %9.4g %9.2f %9.2f", starprob, rhoindx, udotindx);
#if defined(MASSLOSS)
  fprintf(logstr, " %9.4g %9.2f", tau_ml, beta_ml);
#endif
#endif
  fprintf(logstr, "\n\n");
#if defined(GRAVITY)
  fprintf(logstr, " %9s %9s %9s", "eps", "usequad", "theta");
#endif
  fprintf(logstr, " %9s %9s %9s %9s\n", "alpha", "beta", "slope", "fdrag");
#if defined(GRAVITY)
  fprintf(logstr, " %9.4f %9s %9.2f", eps, usequad ? "true" : "false", theta);
#endif
  fprintf(logstr, " %9.2f %9.2f %9.4f %9.4f\n", alpha, beta, slope0, fdrag);
  if (! strnull(options))			// print options, if any
    fprintf(logstr, "\n\toptions: %s\n", options);
  fflush(logstr);
  if (! strnull(savefile))			// was state file given?
    savestate();				// save initial data
  if (! strnull(tracefile)) {
    tracestr = stropen(tracefile, "a");
    put_history(tracestr);
  } else
    tracestr = NULL;
}

// OUTPUTHEAD: announce beginning of new time-step.

void outputhead(void)
{
  fprintf(logstr, "\n____________________________________________________\n");
  fprintf(logstr, "time: %.8f    nstep: %d    cputime: %.2f\n",
	  tnow, nstep, cputime());
  fflush(logstr);
}

// OUTPUTDATA: output diagnostics and binary data.

void outputdata(void)
{
  static string *outtags = NULL;
  int n;
  string reqtags[MaxBodyFields], *alltags;
  char namebuf[256];
  struct stat buf;
  stream outstr;

  diagnostics();				// compute all diagnostics
  fprintf(logstr, "\n%9s %9s %9s %9s %9s %9s %9s %9s\n",
	  "Time", "Etot", "Eint", "Ekin", "Epot", "Erad",
	  "|Jtot|", "|Vcom|");
  fprintf(logstr, "%9.4f %9.5f %9.5f %9.5f %9.5f %9.5f %9.6f %9.6f\n",
	  tnow, etot[0], etot[1], etot[2], etot[3], eradiate,
	  absv(amvec), absv(cmvel));
  if (outtags == NULL)				// list requested data
    outtags = burststring(outputs, ", \t\n");
  if (! strnull(outfile) && tnow + (dtime > 0 ? dtime/8 : 0) >= tout) {
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
    sprintf(namebuf, outfile, nstep);		// make up output file name
    if (stat(namebuf, &buf) != 0) {		// no output file exists?
      outstr = stropen(namebuf, "w");		// create & open for output
      put_history(outstr);			// write out history data
    } else					// else file already exists
      outstr = stropen(namebuf, "a");		// reopen and append output
    put_snap(outstr, &btab, &nbody, &tnow, alltags);
    strclose(outstr);				// close up output file
    free(alltags);
    fprintf(logstr, "\n\tdata output to file %s at time %f\n",
	    namebuf, tnow);
    tout += dtout;				// schedule next output
  }
  fflush(logstr);
  if (! strnull(savefile))			// was state file given?
    savestate();				// save data for restart
  if (tracestr != NULL)
    writetrace();
}

// DIAGNOSTICS: compute dynamic and thermodynamic diagnostics.

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

//  ______________________________________________________
//  writetrace: output full particle data to trace stream.

local void writetrace(void)
{
  static string tracetags[] = {
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

  put_snap(tracestr, &btab, &nbody, &tnow, tracetags);
  fflush(tracestr);    
}

// SAVESTATE: write current state to disk file.

void savestate(void)
{
  char namebuf[256];
  stream str;
  static int randsize;
  static void *randstate = NULL;

  sprintf(namebuf, savefile, nstep & 1);	// generate state file name
  str = stropen(namebuf, "w!");			// open state output file
  put_string(str, "program", getargv0());
  put_string(str, "version", getversion());
  put_string(str, "headline", headline);
  put_data(str, "gamma", RealType, &gamma0, 0);
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
  put_data(str, "starprob", RealType, &starprob, 0);
  put_data(str, "rhoindx", RealType, &rhoindx, 0);
  put_data(str, "udotindx", RealType, &udotindx, 0);
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

// RESTORESTATE: restore state from disk file.

void restorestate(void)
{
  stream str;
  string program, version;
  int randsize;
  void *randstate;

  str = stropen(restfile, "r");			// open state input file
  program = get_string(str, "program");
  version = get_string(str, "version");
  if (! streq(program, getargv0()) ||		// check program, version
      ! streq(version, getversion()))
    eprintf("[%s: input state file may be outdated]\n", getargv0());
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
  get_data(str, "starprob", RealType, &starprob, 0);
  get_data(str, "rhoindx", RealType, &rhoindx, 0);
  get_data(str, "udotindx", RealType, &udotindx, 0);
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
