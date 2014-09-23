/****************************************************************************/
/* SPHCODE.C: hierarchical SPH/N-body code.                                 */
/* Copyright (c) 2009 by Joshua E. Barnes, Honolulu, Hawai'i.               */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "datatypes.h"
#include "gsp.h"

#define global					/* don't default to extern  */
#include "sphcode.h"
#include "kdtree.h"
#include "smooth.h"
#include "fixbody.h"

/* Describe code variant on basis of compile-time flags. */

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

#if defined(STARFORM) && defined(COMPVISC)
#  define STAROPT		";Star formation (using udotvis),",
#elif defined(STARFORM)
#  define STAROPT		";Star formation (using udot),",
#elif defined(STARDEATH)
#  define STAROPT               ";Star death (recycling),",
#else
#  define STAROPT
#endif

#if defined(GRAVITY)
#  define DYNOPT		";Self-consistent gravity.",
#elif defined(EXTGRAV)
#  define DYNOPT		";External gravitational field.",
#elif defined(NOACCEL)
#  define DYNOPT		";No accelerations computed.",
#endif

/* Command-line parameters and defaults. */

string defv[] = {		";SPH/N-body simulation code.",
				THERMOPT
				GASOPT
				STAROPT
				DYNOPT
    "in=",			";Input file for initial conditions",
    "out=",			";Output file for N-body frames",
    "gamma=1.666667",		";Ratio of specific heats",
#if defined(RADIATING)
    "uradmax=1.0",		";Uinternal at max. of cooling curve",
    "lambmax=1.0",		";Max. cooling rate at unit density",
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
    "starprob=0.01",		";Star form. prob. per unit time",
    "rhoindx=0.0",		";Power-law index for density",
    "udotindx=0.0",		";Power-law index for dissipation",
    "starseed=12345",		";Random seed for star formation",
#endif
    "alpha=1.0",		";Bulk viscosity parameter",
    "beta=2.0",			";vN-R viscosity parameter",
    "nsmooth=40",		";Bodies in smoothing volume",
    "nbucket=16",		";Bodies in leaves of KD tree",
    "slope=0.0",		";Kernel slope at origin",
    "courant=0.25",		";Courant condition parameter",
    "dtime=1/64",		";Basic integration timestep",
    "fdrag=0.0",		";Velocity damping factor (1/t)",
#if defined(GRAVITY)
    "eps=0.025",		";Gravitational smoothing length",
    "usequad=false",		";If true, use quad moments",
#if !defined(QUICKSCAN)
    "theta=1.0",		";Force accuracy parameter",
#endif
#elif defined(EXTGRAV)
    "gravgsp=",			";Input GSP for external gravity",
#endif
    "options=",			";Various control options.",
				";Choices include: corrfunc, levelhist,",
				";new-tout, nolimit, fixstep, lockstep,",
				";reset-time, bh86, and sw94.",
    "outputs=",			";Optional particle data written to output.",
				";Key dynamical variables always included.",
    "tstop=2.0",		";Time to stop integration",
    "dtout=1/16",		";Data output timestep",
    "save=",			";Write state file as code runs",
    "restore=",			";Continue run from state file",
    "VERSION=1.2x",		";Joshua Barnes  02 December 2009",
    NULL,
};

/*
 * Prototypes for local procedures.
 */

local void startrun(void);			/* initialize system state  */
local void setupbody(void);			/* set offsets for dynbody  */
local void startnewrun(void);			/* get params & input data  */
local void startoldrun(void);			/* restore state & params   */
local void initstep(void);			/* init multistep algorithm */
local void macrostep(void);			/* advance by basic step    */
local void microstep(real *);			/* advance by minimum step  */
local void update(int, real);			/* update forces and levels */
local real coolfunc(bodyptr);			/* compute cooling function */
local void evaluate(bool);			/* do force calculation     */
local void gravcalc(bool);			/* gravitational part       */
local void sphcalc(bool);			/* pressure, etc part       */
local void sum_density(smxptr, int, int);	/* compute SPH densities    */
local void est_tau(smxptr, int, int);		/* estimate optical depth   */
local void sum_derivs(smxptr, int, int);	/* compute SPH derivatives  */
local void setlevels(smxptr);			/* adjust timestep levels   */
local void starform(real);			/* handle star formation    */
local void stardeath(real);			/* handle star destruction  */

/*
 * MAIN: toplevel routine for hierarchical SPH/N-body code.
 */

int main(int argc, string argv[])
{
    initparam(argv, defv);			/* initialize param access  */
    headline = defv[0] + 1;			/* set ident. message	    */
    startrun();					/* get params & input data  */
    startoutput(stdout, defv);			/* activate output code     */
    if (strnull(restfile))			/* start of new run?	    */
	initstep();				/* prime the integrator     */
    outputdata(stdout);				/* report initial diags.    */
    while (tstop - tnow > 0.01 * dtime) {	/* while not past tstop     */
	macrostep();				/* advance step by step     */
	outputdata(stdout);			/* do output after each	    */
    }
    return (0);
}

/*
 * STARTRUN: startup hierarchical N-body code.
 */

local void startrun(void)
{
    infile = getparam("in");			/* set up I/O file names    */
    outfile = getparam("out");
    savefile = getparam("save");
    restfile = getparam("restore");
    setupbody();				/* interface to dynbody     */
    if (! strnull(infile))
	startnewrun();
    else if (! strnull(restfile))
	startoldrun();
    else
	error("%s: must supply file to input or restore\n", getargv0());
    if (scanopt(options, "new-tout"))		/* if output time reset     */
	tout = tnow + dtout;			/* then offset from now     */
}

/*
 * SETUPBODY: inform dynamic body routines of relevant body fields.
 */

local void setupbody(void)
{
    define_body(sizeof(body), Precision, NDIM);
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
#if defined(STARFORM) || defined(STARDEATH)
    define_body_offset(BirthTag, BodyOffset(Birth));
#endif
}

/*
 * STARTNEWRUN: set command-line parameters and read in input data.
 */

local void startnewrun(void)
{
    real dt1, dt2;

    gamma0 = getdparam("gamma");		/* get input parameters     */
#if defined(RADIATING)
    uradmax = getdparam("uradmax");
    lambmax = getdparam("lambmax");
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
    starprob = getdparam("starprob");
    rhoindx = getdparam("rhoindx");
    udotindx = getdparam("udotindx");
    (void) initstate(getiparam("starseed"), randstate, RANDSIZE);
#endif
    alpha = getdparam("alpha");
    beta = getdparam("beta");
    nsmooth = getiparam("nsmooth");
    nbucket = getiparam("nbucket");
    slope0 = getdparam("slope");
    courant = getdparam("courant");
    dtime = (sscanf(getparam("dtime"), "%f/%f", &dt1, &dt2) == 2 ?
	     dt1 / dt2 : getdparam("dtime"));
    fdrag = getdparam("fdrag");
#if defined(GRAVITY)
    eps = getdparam("eps");
    usequad = getbparam("usequad");
#if !defined(QUICKSCAN)
    theta = getdparam("theta");
#endif
#elif defined(EXTGRAV)
    gspfile = getparam("gravgsp");
#endif
    options = getparam("options");
    outputs = getparam("outputs");
    tstop = getdparam("tstop");
    dtout = (sscanf(getparam("dtout"), "%f/%f", &dt1, &dt2) == 2 ?
	     dt1 / dt2 : getdparam("dtout"));
    inputdata();				/* now read inital data     */
    tout = tnow;				/* schedule first output    */
    rsize = 1.0;				/* start root w/ unit cube  */
    nstep = 0;					/* begin counting steps     */
    eradiate = 0.0;				/* zero radiated energy     */
}

/*
 * STARTOLDRUN: restore state from file and set new parameters.
 */

local void startoldrun(void)
{
    real dt1, dt2;

    restorestate(restfile);			/* read in state file       */
    if (getparamstat("gamma") & ARGPARAM)	/* reset explicit params    */
	gamma0 = getdparam("gamma");
#if defined(RADIATING)
    if (getparamstat("uradmax") & ARGPARAM)
	uradmax = getdparam("uradmax");
    if (getparamstat("lambmax") & ARGPARAM)
	lambmax = getdparam("lambmax");
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
    if (getparamstat("starprob") & ARGPARAM)
	starprob = getdparam("starprob");
    if (getparamstat("rhoindx") & ARGPARAM)
	rhoindx = getdparam("rhoindx");
    if (getparamstat("udotindx") & ARGPARAM)
	udotindx = getdparam("udotindx");
    if (getparamstat("starseed") & ARGPARAM)
        (void) initstate(getiparam("starseed"), randstate, RANDSIZE);
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
    if (getparamstat("courant") & ARGPARAM)
	courant = getdparam("courant");
    if (getparamstat("fdrag") & ARGPARAM)
	fdrag = getdparam("fdrag");
#if defined(GRAVITY)
    if (getparamstat("eps") & ARGPARAM)
	eps = getdparam("eps");
    if (getparamstat("usequad") & ARGPARAM)
	usequad = getbparam("usequad");
#if !defined(QUICKSCAN)
    if (getparamstat("theta") & ARGPARAM)
	theta = getdparam("theta");
#endif
#endif
    if (getparamstat("options") & ARGPARAM)
	options = getparam("options");
    if (getparamstat("outputs") & ARGPARAM)
	outputs = getparam("outputs");
    if (getparamstat("tstop") & ARGPARAM)
	tstop = getdparam("tstop");
    if (getparamstat("dtout") & ARGPARAM)
	dtout = (sscanf(getparam("dtout"), "%f/%f", &dt1, &dt2) == 2 ?
		 dt1 / dt2 : getdparam("dtout"));
}

/*
 * INITSTEP: prepare for multi-step integration.
 */

local void initstep(void)
{
    real ftmp, dt;
    bodyptr p;

    for (p = btab; p < btab+nbody; p++)	{	/* loop over all bodies     */
        Frequency(p) = 1.0 / dtime;		/* set default freqency     */
	Flags(p) = INCLUDE | UPDATE;		/* request force evaluation */
    }
    evaluate(TRUE);				/* do initial force calc.   */
    levmax = 0;
    for (p = btab; p < btab+nbody; p++) {	/* loop over all bodies     */
	CurLevel(p) = NewLevel(p);		/* init multi-step level    */
	levmax = MAX(levmax, CurLevel(p));	/* determine maximum level  */
	dt = dtime / (2 << CurLevel(p));	/* compute half timestep    */
	SETV(Vmid(p), Vel(p));			/* init midpoint velocity   */
	ADDMULVS(Vmid(p), Acc(p), dt);		/* forward by half a step   */
    }
}

/*
 * MACROSTEP: advance system by dtime.
 */

local void macrostep(void)
{
    real fstep;
    int lev, n, m;

#if defined(TRACESTEP)
    fprintf(stdout, "\n    %12s%8s%8s\n", "fstep", "level", "update");
#endif
    fstep = 0.0;				/* count fractional steps   */
    microstep(&fstep);				/* advance system variables */
    while (fstep < 1.0) {			/* loop until step complete */
	lev = levmax;				/* find level of next step  */
	n = (1 << levmax) * fstep;		/* form microstep count     */
	for (m = n ^ (n-1); m > 1; m = m >> 1)
	    lev--;
	update(lev, fstep);			/* update bodies at level   */
	microstep(&fstep);			/* advance variables again  */
    }
    if (fstep > 1.0)
	error("%s: fstep = %f\n", getargv0(), fstep);
    nstep++;					/* advance macrostep count  */
    update(0, fstep);				/* and update all bodies    */
}

/*
 * MICROSTEP: advance system by short time-step.
 */

local void microstep(real *fstep)
{
    real dt, edotrad, udotrad;
    bodyptr p;
    vector vtmp;

    dt = dtime / (1 << levmax);			/* find current time-step   */
    edotrad = 0.0;
    for (p = btab; p < btab+nbody; p++) {
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
			    (gamma0 - 1) / rpow(Rho(p), gamma0 - 1);
	    if (EntFunc(p) < 0.0)
		error("%s: negative entropy function\n", getargv0());
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
		error("%s: negative internal energy\n", getargv0());
#endif
	    edotrad -= Mass(p) * udotrad;
	}
    }
    eradiate += dt * edotrad;
#if defined(STARFORM)
    starform(dt);				/* form dt worth of stars   */
#endif
#if defined(STARDEATH)
    stardeath(dt);				/* kill dt worth of stars   */
#endif
    tnow += dt;
    *fstep += 1.0 / (1 << levmax);
}

/*
 * UPDATE: recompute derivatives for all bodies at current level or higher.
 */

local void update(int lev, real fstep)
{
    int nupdate, oldlev, newlev;
    bodyptr p;
    real dt1, dt2, x;
    vector vdif, aold;

    nupdate = 0;				/* init count of updates    */
    for (p = btab; p < btab+nbody; p++) {	/* loop over body table	    */
	Flags(p) = INCLUDE | (CurLevel(p) < lev ? 0 : UPDATE);
						/* flag bodies to update    */
	if (Update(p))				/* is this one flagged?     */
	    nupdate++;				/* count updating bodies    */
    }
#if defined(TRACESTEP)
    fprintf(stdout, "    %12.5f%8d%8d\n", fstep, lev, nupdate);
#endif
    fflush(stdout);
    evaluate(fstep == 1.0);			/* update flagged bodies    */
    levmax = lev;				/* keep floor on levmax     */
    for (p = btab; p < btab+nbody; p++)		/* loop over all bodies...  */
	if (Update(p)) {			/* which need updating      */
	    oldlev = CurLevel(p);		/* remember previous level  */
	    newlev = MAX(NewLevel(p), MAX(oldlev - 1, lev));
						/* enforce limits on level  */
	    CurLevel(p) = newlev;		/* update level of body     */
	    levmax = MAX(levmax, CurLevel(p));	/* and max level of system  */
	    dt1 = dtime / (2 << oldlev);	/* use half previous step   */
	    dt2 = dtime / (2 << newlev);	/* plus half current step   */
	    SUBV(vdif, Vel(p), Vmid(p));
	    DIVVS(aold, vdif, dt1);		/* recover previous accel.  */
	    SETV(Vel(p), Vmid(p));
	    ADDMULVS(Vel(p), Acc(p), 0.75 * dt1);
	    ADDMULVS(Vel(p), aold, 0.25 * dt1);
	    x = (3 + dt2 / dt1) / 4;		/* set interpolation factor */
	    ADDMULVS(Vmid(p), Acc(p), x * (dt1 + dt2));
	    ADDMULVS(Vmid(p), aold, (1 - x) * (dt1 + dt2));
	}
}

/*
 * EVALUATE: calculate forces and time derivatives.
 */

local void evaluate(bool report)
{
    bodyptr p;

    if (report)					/* log output requested?    */
	outputhead(stdout);			/* announce force calc.     */
    gravcalc(report);				/* do gravity calculation   */
    sphcalc(report);				/* do SPH calculation       */
    if (fdrag != 0.0)				/* are drag forces enabled? */
        for (p = btab; p < btab+nbody; p++)	/* loop over all bodies     */
            if (Gas(p) && Update(p)) {		/* pick updated gas bodies  */
		ADDMULVS(Acc(p), Vel(p), - fdrag);
						/* subtract drag force      */
	}
}

/*
 * GRAVCALC: calculate gravitational forces.
 */

local void gravcalc(bool report)
{
    bodyptr p;
    real r, mrinv3;

#if defined(GRAVITY)
    maketree(btab, nbody);			/* construct tree structure */
    gravforce();				/* compute self-gravity     */
    if (report)					/* log output requested?    */
	report_force(stdout, nbody);		/* print gravity stats.     */
#elif defined(EXTGRAV)
    for (p = btab; p < btab+nbody; p++)		/* loop over all bodies     */
	if (Update(p) && gravgsp != NULL) {	/* update in extern field?  */
	    r = absv(Pos(p));			/* get distance from origin */
	    mrinv3 = - mass_gsp(gravgsp, r) / (r*r*r);
	    MULVS(Acc(p), Pos(p), mrinv3);	/* set extern acc and phi   */
	    Phi(p) = phi_gsp(gravgsp, r);
	} else if (Update(p)) {			/* update in zero field?    */
	    CLRV(Acc(p));			/* zero grav. acceleration  */
	    Phi(p) = 0.0;			/* zero grav. potential	    */
	}
#elif defined(NOACCEL)
    for (p = btab; p < btab+nbody; p++)
	if (Update(p)) {
	    CLRV(Acc(p));
	    Phi(p) = 0.0;
	}
#endif
}

/*
 * SPHCALC: calculate pressure forces and other SPH parameters.
 */

local void sphcalc(bool report)
{
    kdxptr kd;
    smxptr sm;
    bodyptr p;
    real freqrad;

    kd = init_kdtree(btab, nbody, ngas);	/* prepare for kd tree	    */
    build_kdtree(kd, nbucket);			/* do tree construction	    */
    sm = init_smooth(kd, nsmooth, slope0);	/* prepare for smoothing    */
    for (p = btab; p < btab+nbody; p++)		/* loop over all bodies     */
        if (Gas(p)) {				/* which represent gas      */
	    Rho(p) = UdotInt(p) = 0.0;		/* reset cumulative values  */
#if defined(RADIATING)
	    UdotRad(p) = 0.0;
#endif
#if defined(COMPVISC)
	    UdotVis(p) = 0.0;
#endif
	}
    smooth(sm, sum_density);			/* compute smoothed density */
#if defined(DIFFUSING) || defined(OPAQUE)
    resmooth(sm, est_tau);
#endif
    for (p = btab; p < btab+nbody; p++)		/* loop over all bodies	    */
	if (Gas(p)) {				/* which represent gas	    */
#if defined(ENTROPY)
	    Press(p) = EntFunc(p) * rpow(Rho(p), gamma0);
						/* set presure using a(s)   */
#else
	    Press(p) = (gamma0 - 1) * Uintern(p) * Rho(p);
						/* set presure using uint   */
#endif
	}
    resmooth(sm, sum_derivs);			/* compute SPH derivatives  */
#if defined(RADIATING)
    for (p = btab; p < btab+nbody; p++)		/* loop over all bodies	    */
	if (Gas(p)) {				/* which represent gas	    */
#if defined(DIFFUSING) || defined(OPAQUE)
	    UdotRad(p) -= (Tau(p) < 1.0 ? coolfunc(p) : 0.0);
#if defined(DIFFUSING)
	    if (Surface(p))
		UdotRad(p) -= sigmastar * rsqr(rsqr(Uintern(p))) /
			        (2 * Smooth(p) * Rho(p));
#endif
#else
	    UdotRad(p) -= coolfunc(p);
#endif
	    freqrad = (ABS(UdotRad(p)) / Uintern(p)) / courant;
	    Frequency(p) = MAX(Frequency(p), freqrad);
	}
#endif
    setlevels(sm);				/* set timestep levels      */
    if (report)					/* log output requested?    */
	report_smooth(sm, stdout, options);	/* print smoothing stats.   */
    finish_smooth(sm);				/* deallocate smooth data   */
    finish_kdtree(kd);				/* deallocate kdtree data   */
}

/*
 * SUM_DENSITY: sum smoothed density for body and its neighbors.
 */

local void sum_density(smxptr sm, int pi, int nball)
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

/*
 * EST_TAU: estimate optical depth using gradient of density.
 */

local void est_tau(smxptr sm, int pi, int nball)
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

/*
 * SUM_DERIVS: sum SPH accelerations and internal energy derivatives.
 */

local void sum_derivs(smxptr sm, int pi, int nball)
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
            if (rv < 0) {                       /* compute artificial visc. */
                csavg = (csi + rsqrt(gamma0 * Press(bj) / Rho(bj))) / 2;
                rhoavg = (Rho(bi) + Rho(bj)) / 2;
                pivis = muvis * (beta * muvis - alpha * csavg) / rhoavg;
#if defined(ENTROPY)
		UdotInt(bi) += Mass(bj) * 0.5 * pivis * rv * dwsmrinv;
		UdotInt(bj) += Mass(bi) * 0.5 * pivis * rv * dwsmrinv;
#endif
            } else
                pivis = 0.0;
#if !defined(NOACCEL)
            ascale = (epsi + epsj + pivis) * dwsmrinv;
	    if (Update(bi)) {
		ADDMULVS(Acc(bi), rij, - Mass(bj) * ascale);
	    }
	    if (Update(bj)) {
		ADDMULVS(Acc(bj), rij, Mass(bi) * ascale);
	    }
#endif

#if !defined(ENTROPY)
	    UdotInt(bi) += Mass(bj) * (epsi + 0.5 * pivis) * rv * dwsmrinv;
	    UdotInt(bj) += Mass(bi) * (epsj + 0.5 * pivis) * rv * dwsmrinv;
#  if defined(COMPVISC)
	    UdotVis(bi) += Mass(bj) * 0.5 * pivis * rv * dwsmrinv;
	    UdotVis(bj) += Mass(bi) * 0.5 * pivis * rv * dwsmrinv;
#  endif
#  if defined(RADIATING) && defined(DIFFUSING)
	    chii = (16.0/3.0) * (sigmastar / opacity) *
		     (Tau(bi) > 1.0 ? rqbe(Uintern(bi)) / Rho(bi) : 0.0);
	    chij = (16.0/3.0) * (sigmastar / opacity) *
		     (Tau(bj) > 1.0 ? rqbe(Uintern(bj)) / Rho(bj) : 0.0);
	    chiavg = 0.5 * (chii + chij);
	    UdotInt(bi) += Mass(bj) * (chiavg / (Rho(bi) * Rho(bj))) *
			     (Uintern(bi) - Uintern(bj)) * dwsm *
			       rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
	    UdotInt(bj) += Mass(bi) * (chiavg / (Rho(bi) * Rho(bj))) *
			     (Uintern(bj) - Uintern(bi)) * dwsm *
			       rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
	    if ((Tau(bi) > 1.0) != (Tau(bj) > 1.0))
		SetFlag(bi, SURFACE);
#  endif  /* RADIATING && DIFFUSING */
#  if defined(CONDUCTING)
	    UdotInt(bi) += Mass(bj) * (conduct / (Rho(bi) * Rho(bj))) *
			     (Uintern(bi) - Uintern(bj)) * dwsm *
			       rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
	    UdotInt(bj) += Mass(bi) * (conduct / (Rho(bi) * Rho(bj))) *
			     (Uintern(bj) - Uintern(bi)) * dwsm *
			       rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
#  endif  /* CONDUCTING */
#endif  /* ! ENTROPY */
            divv -= 2 * Mass(bj) * rv * dwsmrinv / Rho(bi);
        }
    }
    if (Update(bi))
	Frequency(bi) = ((csi + 1.2 * (alpha*csi + beta*mumax)) / hsi +
			   ABS(divv)) / courant;
}

#if defined(RADIATING)

/*
 * COOLFUNC: evaluate standard cooling function.
 */

local real coolfunc(bodyptr p)
{
    real logx, f1, f2, f3, udotrad;

    logx = log10(Uintern(p) / uradmax);
    f1 = -0.17 * rsqr(logx + 1.23) - 3.14;
    f2 = -1.88 * rsqr(rsqr(logx));
    f3 = (logx > 0.97 ? 0.0 : -2.0 * rsqr(rsqr(logx - 0.97)))  - 1.60;
    udotrad = lambmax * Rho(p) * (rdex(f1) + rdex(f2) + rdex(f3));
    return (udotrad);
}    

#endif

/*
 * SETLEVELS: set NewLevel values from computed Frequencies.
 */

local void setlevels(smxptr sm)
{
    bodyptr p;
    real freqmin;
    int lev, nmax;

    sm->freqmax = sm->freqsum = 0.0;
    for (p = btab; p < btab+nbody; p++) {
        if (Gas(p) && Update(p)) {
#if defined(CONDUCTING)
	    freqmin = (ABS(UdotInt(p)) / Uintern(p)) / courant;
	    if (Frequency(p) < freqmin)
	        Frequency(p) = freqmin;
#endif
	    sm->freqmax = MAX(sm->freqmax, Frequency(p));
	    sm->freqsum += Frequency(p);
	    lev = rceil(rlog2(MAX(dtime * Frequency(p), 1.0)));
	    if (lev > MAXLEV) {
	        eprintf("[%s: courant condition violated; level = %d]\n",
		        getargv0(), lev);
		lev = MAXLEV;
	    }
	    NewLevel(p) = lev;
	}
	if (! Gas(p))
	    NewLevel(p) = 0;
    }
    if (! scanopt(options, "nolimit")) {	/* nolimit option not set?  */
	lev = 0;
	for (p = btab; p < btab+nbody; p++)
	    if (NewLevel(p) > lev) {		/* found higher level?      */
		lev = NewLevel(p);
		nmax = 1;
	    } else if (NewLevel(p) == lev)
		nmax++;
	if (nmax < nsmooth/4)			/* less than 25% have max?  */
	    for (p = btab; p < btab+nbody; p++)
		if (NewLevel(p) == lev)		/* demote them one level    */
		    NewLevel(p)--;
    }
    if (scanopt(options, "fixstep")) {		/* fixstep option set?      */
	if (scanopt(options, "lockstep"))
	    error("%s: fixstep and lockstep incompatable\n", getargv0());
	for (p = btab; p < btab+nbody; p++)
	    NewLevel(p) = 0;
    }
    if (scanopt(options, "lockstep")) {		/* lockstep option set?     */
	lev = 0;				/* place floor on level     */
	for (p = btab; p < btab+nbody; p++)
	    lev = MAX(NewLevel(p), lev);	/* determine maximum level  */
	for (p = btab; p < btab+nbody; p++)
	    NewLevel(p) = lev;			/* assign max level to all  */
    }
}

stream starlog = NULL;

#if defined(STARFORM)

/*
 * STARFORM: handle star formation as a random process.
 */

local void starform(real dt)
{
  bodyptr p;
  real udot, prob;

  if (starlog == NULL)
    starlog = stropen("starlog.txt", "w!");
  for (p = btab; p < btab+nbody; p++)
    if (Gas(p)) {
#if defined(COMPVISC)
      udot = UdotVis(p);
#else	
      udot = MAX(UdotInt(p), 0.0);
#endif
      prob = starprob * rpow(Rho(p), rhoindx) * rpow(udot, udotindx) * dt;
      if (xrandom(0.0, 1.0) < prob) {		/* a star is born!          */
	Type(p) = BODY | STAR;
	eradiate += Mass(p) * Uintern(p);
	Birth(p) = tnow;
	ngas--;
	fprintf(starlog, "%8d  %10.8f  form  %10.8f\n", p - btab, tnow, prob);
      }
    }
  fflush(starlog);
}

#endif

#if defined(STARDEATH)

#define AlphaSD  (1.0/12.0)
#define TauSD    0.01

local real deathrate(real tau)
{
  if (tau < 0)
    error("%s.deathrate: tau = %f < 0\n", getargv0(), tau);
  if (tau < TauSD)
    return (0.0);
  return ((AlphaSD / TauSD) * rpow(tau / TauSD, - (AlphaSD + 1.0)));
}

local void stardeath(real dt)
{
  bodyptr p;

  if (starlog == NULL)
    starlog = stropen("starlog.txt", "w!");
  for (p = btab; p < btab+nbody; p++)
    if (Birth(p) > 0 && xrandom(0.0, 1.0) < dt * deathrate(tnow - Birth(p))) {
      Type(p) = BODY | GAS;
      eradiate -= Mass(p) * Uintern(p);
      ngas++;
      fprintf(starlog, "%8d  %10.8f  dead  %10.8f\n",
	      p - btab, tnow, dt * deathrate(tnow - Birth(p)));
      Birth(p) = 0;
    } else if (Birth(p) < 0 && -Birth(p) < tnow + dt) {
      if (Uintern(p) > 0)
	error("%s.stardeath: doomed star %d has positive uint\n",
	      getargv0(), p - btab);
      Type(p) = BODY | GAS;
      Uintern(p) = -Uintern(p);
      eradiate -= Mass(p) * Uintern(p);
      ngas++;
      fprintf(starlog, "%8d  %10.8f  doom\n", p - btab, tnow);
      Birth(p) = 0;
    }
  fflush(starlog);
}

#endif
