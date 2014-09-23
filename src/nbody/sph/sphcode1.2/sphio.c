/****************************************************************************/
/* SPHIO.C: I/O routines for hierarchical SPH code. Public routines:        */
/* inputdata(), startoutput(), outputdata(), savestate(), restorestate().   */
/* Copyright (c) 2009 by Joshua E. Barnes, Honolulu, Hawai'i.               */
/****************************************************************************/

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

/*
 * Prototypes for local routines.
 */

local void diagnostics(void);			/* eval N-body diagnostics  */

/*
 * Output state variables.
 */

local real mtot;		                /* total mass of system     */
local real etot[4];				/* system E, Ek, Ei, Ep     */
local matrix keten;				/* kinetic energy tensor    */
local matrix peten;				/* potential energy tensor  */
local vector cmpos;				/* center of mass position  */
local vector cmvel;				/* center of mass velocity  */
local vector amvec;				/* angular momentum vector  */

/*
 * INPUTDATA: read initial conditions from input file.
 */

void inputdata(void)
{
    stream istr, gstr;
    string intags[MaxBodyFields];
    bodyptr p;

    btab = NULL;				/* request new input data   */
    istr = stropen(infile, "r");		/* open input stream        */
    get_history(istr);				/* read file history data   */
    if (! get_snap(istr, &btab, &nbody, &tnow, intags, FALSE))
	error("inputdata: no data in input file\n");
    strclose(istr);				/* close input stream       */
    if (! set_member(intags, MassTag))
	error("inputdata: mass data missing\n");
    if (! set_member(intags, PosTag))
	error("inputdata: position data missing\n");
    if (! set_member(intags, VelTag))
	error("inputdata: velocity data missing\n");
#if defined(ENTROPY)
    if (! set_member(intags, EntFuncTag))
	error("inputdata: entropy function data missing\n");
#else
    if (! set_member(intags, UinternTag))
	error("inputdata: internal energy data missing\n");
#endif
#if defined(STARDEATH)
    if (! set_member(intags, BirthTag))
      error("inputdata: birth data missing\n");
#endif
    ngas = 0;					/* start gas body count     */
    for (p = btab; p < btab+nbody; p++) {	/* initialize node type     */
#if defined(ENTROPY)
	Type(p) = BODY | (EntFunc(p) >= 0 ? GAS : 0);
						/* gas has nonnegative a(s) */
#else
	Type(p) = BODY | (Uintern(p) >= 0 ? GAS : 0);
						/* gas has nonnegative uint */
#endif
	if (Gas(p))				/* was gas bit set?	    */
	    ngas++;				/* then count gas body	    */
    }
#if defined(EXTGRAV)
    if (! strnull(gspfile)) {			/* was GSP file given?      */
	gstr = stropen(gspfile, "r");
	get_history(gstr);
	gravgsp = get_gsprof(gstr);		/* set external field GSP   */
	strclose(gstr);
    } else
	gravgsp = NULL;				/* else no external field   */
#endif
    if (scanopt(options, "reset-time"))		/* reset starting time?     */
	tnow = 0.0;
}

/*
 * STARTOUTPUT: begin output to log file.
 */

void startoutput(stream ostr, string defv[])
{
    int i;

    fprintf(ostr, "\n%s  [v%s]\n", headline, getversion());
    for (i = 1; defv[i][0] == ';'; i++)
        fprintf(ostr, "    %s\n", defv[i]+1);
    fprintf(ostr, "\n %9s %9s %9s %9s %9s %9s %9s\n", "nbody", "ngas",
	    "nsmooth", "courant", "dtime", "dtout", "tstop");
    fprintf(ostr, " %9d %9d %9d %9.4f %9.6f %9.6f %9.4f\n\n", nbody, ngas,
	    nsmooth, courant, dtime, dtout, tstop);
    fprintf(ostr, " %9s", "gamma");
#if defined(RADIATING)
    fprintf(ostr, " %9s %9s", "uradmax", "lambmax");
#endif
#if defined(DIFFUSING)
    fprintf(ostr, " %9s", "sigmastar");
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    fprintf(ostr, " %9s", "opacity");
#endif
#if defined(CONDUCTING)
    fprintf(ostr, " %9s", "conduct");
#endif
#if defined(STARFORM)
    fprintf(ostr, " %9s %9s %9s", "starprob", "rhoindx", "udotindx");
#endif
    fprintf(ostr, "\n");
    fprintf(ostr, " %9.4f", gamma0);
#if defined(RADIATING)
    fprintf(ostr, " %9.4g %9.4g", uradmax, lambmax);
#endif
#if defined(DIFFUSING)
    fprintf(ostr, " %9.4g", sigmastar);
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    fprintf(ostr, " %9.4g", opacity);
#endif
#if defined(CONDUCTING)
    fprintf(ostr, " %9.4g", conduct);
#endif
#if defined(STARFORM)
    fprintf(ostr, " %9.4g %9.2f %9.2f", starprob, rhoindx, udotindx);
#endif
    fprintf(ostr, "\n\n");
#if defined(GRAVITY)
    fprintf(ostr, " %9s %9s", "eps", "usequad");
#if !defined(QUICKSCAN)
    fprintf(ostr, " %9s", "theta");
#endif
#endif
    fprintf(ostr, " %9s %9s %9s %9s\n", "alpha", "beta", "slope", "fdrag");
#if defined(GRAVITY)
    fprintf(ostr, " %9.4f %9s", eps, usequad ? "true" : "false");
#if !defined(QUICKSCAN)
    fprintf(ostr, " %9.2f", theta);
#endif
#endif
    fprintf(ostr, " %9.2f %9.2f %9.4f %9.4f\n", alpha, beta, slope0, fdrag);
    if (! strnull(options))			/* print options, if any    */
        fprintf(ostr, "\n\toptions: %s\n", options);
    fflush(ostr);
    if (! strnull(savefile))			/* was state file given?    */
	savestate(savefile);			/* save initial data        */
}

/*
 * OUTPUTHEAD: announce beginning of new time-step.
 */

void outputhead(stream ostr)
{
    fprintf(ostr, "\n");
    fprintf(ostr, "########## %9s %9s    ################"
	    "############################\n", "Time", "Nstep");
    fprintf(ostr, "########## %9.4f %9d    ################"
	    "############################\n", tnow, nstep);
    fflush(ostr);
}

/*
 * OUTPUTDATA: compute diagnostics and output binary data.
 */

void outputdata(stream ostr)
{
    real teff;
    int n;
    string reqtags[MaxBodyFields], *outtags, *alltags;
    char namebuf[256];
    struct stat buf;
    stream outstr;

    diagnostics();				/* compute all diagnostics  */
    fprintf(ostr, "\n %9s %9s %9s %9s %9s %9s %9s %7s\n",
	   "Etot", "Eint", "Ekin", "Epot", "Erad",
	    "|Jtot|", "|Vcom|", "CPUtot");
    fprintf(ostr, " %9.5f %9.5f %9.5f %9.5f %9.5f %9.6f %9.6f %7.1f\n",
	   etot[0], etot[1], etot[2], etot[3], eradiate,
	    absv(amvec), absv(cmvel), cputime());
    teff = tnow + (dtime > 0 ? dtime/8 : 0);	/* anticipate slightly...   */
    if (! strnull(outfile) && teff >= tout) {	/* time for data output?    */
	n = 0;					/* list required data	    */
	if (scanopt(outputs, MassTag) || (nstep == 0))
	    reqtags[n++] = MassTag;		/* output mass data         */
	reqtags[n++] = PosTag;			/* output position data	    */
	reqtags[n++] = VelTag;			/* output velocity data	    */
#if defined(ENTROPY)
#  if defined(ADIABATIC)
	if (scanopt(outputs, UinternTag) || (nstep == 0))
	    reqtags[n++] = EntFuncTag;
#  else
	reqtags[n++] = EntFuncTag;		/* output entropy data	    */
#  endif
#else
#  if defined(ISOTHERMAL)
	if (scanopt(outputs, UinternTag) || (nstep == 0))
	    reqtags[n++] = UinternTag;		/* output energy data	    */
#  else
	reqtags[n++] = UinternTag;		/* output energy data	    */
#  endif
#endif
	reqtags[n] = NULL;			/* terminate required list  */
	outtags = burststring(outputs, ", ");	/* list requested data	    */
	alltags = set_union(reqtags, outtags);	/* list all output data	    */
	sprintf(namebuf, outfile, nstep);	/* make up output file name */
	if (stat(namebuf, &buf) != 0) {		/* no output file exists?   */
	    outstr = stropen(namebuf, "w");     /* create & open for output */
	    put_history(outstr);		/* write out hiatory data   */
	} else					/* else file already exists */
	    outstr = stropen(namebuf, "a");	/* reopen and append output */
	put_snap(outstr, &btab, &nbody, &tnow, alltags);
	strclose(outstr);			/* close up output file     */
	free(outtags);				/* collect tag lists	    */
	free(alltags);
	fprintf(ostr, "\n\tdata output to file %s at time %f\n",
		namebuf, tnow);
	tout += dtout;				/* schedule next output     */
    }
    fflush(ostr);
    if (! strnull(savefile))			/* was state file given?    */
	savestate(savefile);			/* save data for restart    */
}

/*
 * DIAGNOSTICS: compute dynamical and thermodynamical diagnostics.
 */

local void diagnostics(void)
{
    register bodyptr p;
    real uint, velsq;
    vector tmpv;
    matrix tmpt;

    mtot = 0;					/* zero total mass          */
    etot[1] = etot[2] = etot[3] = 0;		/* zero energy totals       */
    CLRM(keten);				/* zero ke tensor           */
    CLRM(peten);				/* zero pe tensor           */
    CLRV(amvec);				/* zero am vector           */
    CLRV(cmpos);				/* zero c. of m. position   */
    CLRV(cmvel);				/* zero c. of m. velocity   */
    for (p = btab; p < btab+nbody; p++) {	/* loop over all particles  */
	mtot += Mass(p);                        /* sum particle masses      */
	if (Gas(p)) {
#if defined(ENTROPY)
	    uint = EntFunc(p) * rpow(Rho(p), gamma0 - 1) / (gamma0 - 1);
#else
	    uint = Uintern(p);
#endif
	    etot[1] += Mass(p) * uint;		/* sum internal energies    */
	}
	DOTVP(velsq, Vel(p), Vel(p));		/* square vel vector        */
	etot[2] += 0.5 * Mass(p) * velsq;	/* sum current KE           */
#if defined(GRAVITY)
	etot[3] += 0.5 * Mass(p) * Phi(p);	/* and current PE           */
#elif defined(EXTGRAV)
	etot[3] += Mass(p) * Phi(p);		/* and external PE          */
#endif
	MULVS(tmpv, Vel(p), 0.5 * Mass(p));	/* sum 0.5 m v_i v_j        */
	OUTVP(tmpt, tmpv, Vel(p));
	ADDM(keten, keten, tmpt);
	MULVS(tmpv, Pos(p), Mass(p));		/* sum m r_i a_j            */
	OUTVP(tmpt, tmpv, Acc(p));
	ADDM(peten, peten, tmpt);
	CROSSVP(tmpv, Vel(p), Pos(p));		/* sum angular momentum     */
	MULVS(tmpv, tmpv, Mass(p));
	ADDV(amvec, amvec, tmpv);
	MULVS(tmpv, Pos(p), Mass(p));		/* sum cm position          */
	ADDV(cmpos, cmpos, tmpv);
	MULVS(tmpv, Vel(p), Mass(p));		/* sum cm momentum          */
	ADDV(cmvel, cmvel, tmpv);
    }
    etot[0] = etot[1] + etot[2] + etot[3] + eradiate;
						/* sum total energy         */
    DIVVS(cmpos, cmpos, mtot);        		/* normalize cm coords      */
    DIVVS(cmvel, cmvel, mtot);
}

/*
 * SAVESTATE: write current state to disk file.
 */

void savestate(string pattern)
{
    char namebuf[256];
    stream str;

    sprintf(namebuf, pattern, nstep & 1);	/* generate state file name */
    str = stropen(namebuf, "w!");		/* open state output file   */
    put_string(str, "program", getargv0());
    put_string(str, "version", getversion());
    put_string(str, "headline", headline);
    put_data(str, "gamma", RealType, &gamma0, 0);
#if defined(RADIATING)
    put_data(str, "uradmax", RealType, &uradmax, 0);
    put_data(str, "lambmax", RealType, &lambmax, 0);
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
#if !defined(QUICKSCAN)
    put_data(str, "theta", RealType, &theta, 0);
#endif
#elif defined(EXTGRAV)
    put_string(str, "gspfile", gspfile);
    if (! strnull(gspfile))
	put_gsprof(str, gravgsp);
#endif
    put_string(str, "options", options);
    put_string(str, "outputs", outputs);
    put_data(str, "tstop", RealType, &tstop, 0);
    put_data(str, "dtout", RealType, &dtout, 0);
    put_data(str, "nstep", IntType, &nstep, 0);	/* save current state vars  */
    put_data(str, "levmax", IntType, &levmax, 0);
    put_data(str, "rsize", RealType, &rsize, 0);
    put_data(str, "eradiate", RealType, &eradiate, 0);
    put_data(str, "tout", RealType, &tout, 0);
    put_data(str, "tnow", RealType, &tnow, 0);
#if defined(STARFORM)
    put_data(str, "randstate", ByteType, randstate, RANDSIZE, 0);
#endif
    put_data(str, "nbody", IntType, &nbody, 0);
    put_data(str, "ngas", IntType, &ngas, 0);
    put_data(str, "btab", AnyType, btab, nbody, sizeof(body), 0);
    strclose(str);
}

/*
 * RESTORESTATE: restore state from disk file.
 */

void restorestate(string file)
{
    stream str, gstr;
    string program, version;

    str = stropen(file, "r");			/* open state input file    */
    program = get_string(str, "program");
    version = get_string(str, "version");
    if (! streq(program, getargv0()) ||		/* check program, version   */
	  ! streq(version, getversion()))
	eprintf("[%s: input state file may be outdated]\n", getargv0());
    headline = get_string(str, "headline");
    get_data(str, "gamma", RealType, &gamma0, 0);
#if defined(RADIATING)
    put_data(str, "uradmax", RealType, &uradmax, 0);
    put_data(str, "lambmax", RealType, &lambmax, 0);
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
#if !defined(QUICKSCAN)
    get_data(str, "theta", RealType, &theta, 0);
#endif
#elif defined(EXTGRAV)
    gspfile = get_string(str, "gspfile");
    if (! strnull(gspfile))			/* was GSP file given?      */
	gravgsp = get_gsprof(str);		/* set external field GSP   */
    else
	gravgsp = NULL;				/* else no external field   */
#endif
    options = get_string(str, "options");
    outputs = get_string(str, "outputs");
    get_data(str, "tstop", RealType, &tstop, 0);
    get_data(str, "dtout", RealType, &dtout, 0);
    get_data(str, "nstep", IntType, &nstep, 0);	/* read current state vars  */
    get_data(str, "levmax", IntType, &levmax, 0);
    get_data(str, "rsize", RealType, &rsize, 0);
    get_data(str, "eradiate", RealType, &eradiate, 0);
    get_data(str, "tout", RealType, &tout, 0);
    get_data(str, "tnow", RealType, &tnow, 0);
#if defined(STARFORM)
    get_data(str, "randstate", ByteType, randstate, RANDSIZE, 0);
    (void) setstate(randstate);
#endif
    get_data(str, "nbody", IntType, &nbody, 0);
    get_data(str, "ngas", IntType, &ngas, 0);
    btab = (bodyptr) allocate(nbody * sizeof(body));
    get_data(str, "btab", AnyType, btab, nbody, sizeof(body), 0);
    strclose(str);
}
