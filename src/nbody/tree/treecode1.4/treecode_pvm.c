/****************************************************************************/
/* TREECODE_PVM.C: hierarchical N-body code for PVM.                        */
/* Copyright (c) 2000 by Joshua E. Barnes, Honolulu, HI.                    */
/****************************************************************************/

#define global					/* don't default to extern  */
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "datatypes.h"
#include "treecode.h"
#include "fixbody.h"
#include "snapcenter.h"
#include "pvm3.h"
#include <sys/time.h>
#include <sys/resource.h>

/* Default values for input parameters. */

string defv[] = {		";Hierarchical PVM N-body code "
#ifdef QUICKSCAN
				    "(quick scan)",
#else
				    "(theta scan)",
#endif
    "in=",			";Input file with initial conditions",
    "out=",			";Output file pattern for N-body frames.",
				";Replace %d with step number of frame.",
    "restore=",			";Continue run from state file",
    "save=",			";Write state file as code runs",
    "dtime=1/32",		";Leapfrog integration timestep",
    "nstatic=0",		";Number of static bodies in array.",
				";If pos (neg), count from start (end).",
    "eps=0.025",		";Density field smoothing length",
#ifndef QUICKSCAN
    "theta=1.0",		";Force accuracy parameter",
#endif
    "usequad=false",		";If true, use quad moments",
    "options=",			";Various control options.",
				";Choices: new-tout,reset-time,bh86,sw94.",
    "outputs=" PosTag "," VelTag,
				";Data fields to output.",
				";Others: " MassTag "," PhiTag "," AccTag ".",
    "tstop=2.0",		";Time to stop integration",
    "dtout=1/4",		";Data output timestep",
    "nbody=16384",		";Number of particles (for test run).",
				";If no input file, make Plummer model.",
    "seed=123",			";Random number seed (for test run)",
    "ntrace=0",			";Number of trace bodies output each step",
    "trace=",			";Output file pattern for trace frames",
    "ntask=???",		";Number of tasks to run",
    "slack=0.95",		";Slack for leading task",
    "VERSION=1.4",		";Joshua Barnes  14 April 2001",
    NULL,
};

/* Prototypes for local procedures. */

local void startrun(void);			/* initialize treecode      */
local void loadstate(void);			/* setup params and bodies  */
local void testdata(void);			/* generate test data       */
local void sharestate(void);			/* send state to helpers    */
local void stepsystem(void);			/* advance by one time-step */
local void treeforce(void);			/* do force calculation     */
local void sharestats(void);			/* exchange gravity stats.  */
local void sharegravity(void);			/* exchange gravity data    */
local void stoprun(void);			/* shut down PVM tasks      */

/* Definitions and parameters for parallel calculation. */

#define NiceValue  10			/* parameter for nice option	    */

#define ProgState  0			/* initial program state message    */
#define GravStats  1			/* force calc. statistics message   */
#define AccData    2			/* body acceleration data message   */
#define PhiPrompt  3			/* request potential data message   */
#define PhiData    4			/* body potential data message      */

local int ntask;			/* number of cooperating tasks	    */
local int jtask;			/* index of this task [0:ntask-1]   */
local int *taskid;			/* task identification numbers      */
local int *nload;			/* task workloads (no. of bodies)   */
local real *twall;			/* task clock time per step         */
local real slack;			/* extra time for leading process   */

/* Synonyms for PVM routines. */

#ifndef DOUBLEPREC
#  define pvm_pkreal  pvm_pkfloat
#  define pvm_upkreal pvm_upkfloat
#else
#  define pvm_pkreal  pvm_pkdouble
#  define pvm_upkreal pvm_upkdouble
#endif

#define pvm_pkbool  pvm_pkshort
#define pvm_upkbool pvm_upkshort

/*
 * MAIN: toplevel routine for hierarchical N-body code.
 */

int main(int argc, string argv[])
{
    initparam(argv, defv);			/* initialize param access  */
    headline = defv[0] + 1;			/* use default headline     */
    startrun();					/* get params & input data  */
    if (jtask == 0)				/* if process is leader     */
        startoutput();				/* activate output code     */
    if (nstep == 0) {				/* if data just initialized */
	treeforce();				/* calculate initial forces */
	if (jtask == 0)				/* if process is leader     */
	    output();				/* generate initial output  */
    }
    if (dtime != 0.0)				/* if time steps requested  */
	while (tstop - tnow > 0.01 * dtime) {	/* while not past tstop     */
	    stepsystem();			/* advance step by step     */
	    if (jtask == 0)			/* if process is leader     */
	        output();			/* output results of step   */
	}
    stoprun();					/* finish parallel section  */
    return (0);					/* end with proper status   */
}

/*
 * STARTRUN: startup parallel hierarchical N-body code.
 */

local void startrun(void)
{
    int mytid, spawnflag, j, ntot;
    string spawnargs[3];
    char arg1[65];
    bodyptr p1, p2, p;

    define_body(sizeof(body), Precision, NDIM);	/* setup phat body struct   */
    define_body_offset(PosTag,  BodyOffset(Pos));
    define_body_offset(VelTag,  BodyOffset(Vel));
    define_body_offset(MassTag, BodyOffset(Mass));
    define_body_offset(PhiTag,  BodyOffset(Phi));
    define_body_offset(AccTag,  BodyOffset(Acc));
    ntask = getiparam("ntask");			/* get number of tasks      */
    taskid = (int *) allocate(sizeof(int) * ntask);
    nload = (int *) allocate(sizeof(int) * ntask);
    twall = (real *) allocate(sizeof(real) * ntask);
    mytid = pvm_mytid();			/* get ident of this task   */
    if (pvm_parent() < 0) {			/* if process is leader     */
        infile = getparam("in");		/* init I/O file names      */
	outfile = getparam("out");
	savefile = getparam("save");
	tracefile = getparam("trace");
	ntrace = getiparam("ntrace");		/* init misc parameters     */
	slack = getdparam("slack");
        loadstate();				/* initialize state data    */
	spawnargs[0] = getargv0();		/* set up args for spawn    */
	sprintf(spawnargs[1] = arg1, "ntask=%d", ntask);
	spawnargs[2] = NULL;
	spawnflag = (scanopt(options, "use-host") ?
		       PvmTaskDefault : PvmTaskHost + PvmHostCompl);
	if (pvm_spawn(getargv0(), spawnargs, spawnflag, ".",
		      ntask - 1, &taskid[1]) != ntask - 1)
	    error("%s: can't start %d helpers\n", getargv0(), ntask - 1);
        taskid[0] = mytid;			/* drop our tid into place  */
    }
    pvm_setopt(PvmRoute, PvmRouteDirect);	/* route messages directly  */
    sharestate();				/* send state to helpers    */
    for (jtask = 0; taskid[jtask] != mytid; jtask++)
						/* scan taskid[] for self   */
        if (jtask == ntask - 1)			/* complain if not found    */
	    fatal("%s: can't find taskid = %d\n", getargv0(), mytid);
    ntot = 0;					/* init sum of workloads    */
    for (j = 0; j < ntask; j++)	{		/* loop over all tasks      */
        nload[j] = nbody / ntask;		/* assign workload to each  */
        ntot += nload[j];			/* sum overall workload     */
    }
    nload[ntask-1] += nbody - ntot;		/* last task gets remainder */
    testcalc = TRUE;				/* determine type of calc.  */
    p1 = bodytab + MAX(nstatic, 0);		/* set dynamic body range   */
    p2 = bodytab + nbody + MIN(nstatic, 0);
    for (p = p1; p < p2; p++)			/* loop over dynamic bodies */
        testcalc = testcalc && (Mass(p) == 0);	/* look for massive ones    */
    if (scanopt(options, "nice"))
        if (setpriority(PRIO_PROCESS, 0, NiceValue) == -1)
	    fatal("%s: task %d can\'t set priority\n", getargv0(), jtask);
}

/*
 * LOADSTATE: initialize parameters and bodies.
 */

local void loadstate(void)
{
    real dt1, dt2;

    if (strnull(getparam("restore"))) {
	eps = getdparam("eps");			/* get input parameters     */
	dtime = (sscanf(getparam("dtime"), "%f/%f", &dt1, &dt2) == 2 ?
		 dt1 / dt2 : getdparam("dtime"));
	nstatic = getiparam("nstatic");
#ifndef QUICKSCAN
	theta = getdparam("theta");
#endif
	usequad = getbparam("usequad");
	tstop = getdparam("tstop");
	dtout = (sscanf(getparam("dtout"), "%f/%f", &dt1, &dt2) == 2 ?
		 dt1 / dt2 : getdparam("dtout"));
	options = getparam("options");
	outputs = getparam("outputs");
	if (! strnull(infile))			/* if data file was given   */
	    inputdata();			/* then read inital data    */
	else {					/* else make initial data   */
	    nbody = getiparam("nbody");		/* get number of bodies     */
	    srandom(getiparam("seed"));		/* set random number gen.   */
	    testdata();				/* and make plummer model   */
	}
	rsize = 1.0;				/* start root w/ unit cube  */
	nstep = 0;				/* begin counting steps     */
	tout = tnow;				/* schedule first output    */
    } else {
	restorestate(getparam("restore"));	/* read in old state file   */
	if (getparamstat("eps") & ARGPARAM)	/* was eps given new value? */
	    eps = getdparam("eps");		/* use command line value   */
	if (getparamstat("dtime") & ARGPARAM)
	    error("%s: can't change dtime\n", getargv0());
	if (getparamstat("nstatic") & ARGPARAM)
	    nstatic = getiparam("nstatic");
    #ifndef QUICKSCAN
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
	    dtout = (sscanf(getparam("dtout"), "%f/%f", &dt1, &dt2) == 2 ?
		     dt1 / dt2 : getdparam("dtout"));
	if (scanopt(options, "new-tout"))	/* if output time reset     */
	    tout = tnow;			/* then offset from now     */
    }
    if (ABS(nstatic) > nbody)			/* check nstatic is OK      */
      error("%s: absurd value for nstatic\n", getargv0());
}

/*
 * TESTDATA: generate Plummer model initial conditions for test runs,
 * scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 */

#define MFRAC  0.999				/* cut off 1-MFRAC of mass  */

local void testdata(void)
{
    real rsc, vsc, r, v, x, y;
    bodyptr p;

    if (nbody < 1)				/* check for silly values   */
	error("%s: absurd value for nbody\n", getargv0());
    bodytab = (bodyptr) allocate(nbody * sizeof(body));
						/* alloc space for bodies   */
    rsc = (3 * PI) / 16;			/* and length scale factor  */
    vsc = rsqrt(1.0 / rsc);			/* find speed scale factor  */
    for (p = bodytab; p < bodytab+nbody; p++) {	/* loop over particles      */
	Type(p) = BODY;				/* tag as a body            */
	Mass(p) = 1.0 / nbody;			/* set masses equal         */
	x = xrandom(0.0, MFRAC);		/* pick enclosed mass       */
	r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1);	/* find enclosing radius    */
	pickshell(Pos(p), NDIM, rsc * r);	/* pick position vector     */
	do {					/* select from fn g(x)      */
	    x = xrandom(0.0, 1.0);		/* for x in range 0:1       */
	    y = xrandom(0.0, 0.1);		/* max of g(x) is 0.092     */
	} while (y > x*x * rpow(1 - x*x, 3.5));	/* using von Neumann tech   */
	v = x * rsqrt(2.0 / rsqrt(1 + r*r));	/* find resulting speed     */
	pickshell(Vel(p), NDIM, vsc * v);	/* pick velocity vector     */
    }
    snapcenter(bodytab, nbody, BodyOffset(Mass));
						/* subtract cm coordinates  */
    tnow = 0.0;					/* set elapsed model time   */
}

/*
 * SHARESTATE: share state information with helpers.
 */

local void sharestate(void)
{
    int nchr;

    if (pvm_parent() < 0) {			/* if task is the leader    */
	pvm_initsend(PvmDataInPlace);
	pvm_pkint(taskid, ntask, 1);
	pvm_pkreal(&dtime, 1, 1);
	pvm_pkint(&nstatic, 1, 1);
#ifndef QUICKSCAN
	pvm_pkreal(&theta, 1, 1);
#endif
	pvm_pkbool(&usequad, 1, 1);
	pvm_pkreal(&eps, 1, 1);
	nchr = xstrlen(options, sizeof(char));
	pvm_pkint(&nchr, 1, 1);
	pvm_pkbyte(options, nchr, 1);
	pvm_pkreal(&tstop, 1, 1);
	pvm_pkreal(&dtout, 1, 1);
	pvm_pkreal(&tnow, 1, 1);
	pvm_pkreal(&tout, 1, 1);
	pvm_pkint(&nstep, 1, 1);
	pvm_pkreal(&rsize, 1, 1);
	pvm_pkreal(&slack, 1, 1);
	pvm_pkint(&nbody, 1, 1);
	pvm_pkbyte((char *) bodytab, nbody * sizeof(body), 1);
	pvm_mcast(taskid, ntask, ProgState);
    } else {					/* else task is a helper    */
	pvm_recv(pvm_parent(), ProgState);
	pvm_upkint(taskid, ntask, 1);
	pvm_upkreal(&dtime, 1, 1);
	pvm_upkint(&nstatic, 1, 1);
#ifndef QUICKSCAN
	pvm_upkreal(&theta, 1, 1);
#endif
	pvm_upkbool(&usequad, 1, 1);
	pvm_upkreal(&eps, 1, 1);
	pvm_upkint(&nchr, 1, 1);
	options = (char *) allocate(nchr * sizeof(char));
	pvm_upkbyte(options, nchr, 1);
	pvm_upkreal(&tstop, 1, 1);
	pvm_upkreal(&dtout, 1, 1);
	pvm_upkreal(&tnow, 1, 1);
	pvm_upkreal(&tout, 1, 1);
	pvm_upkint(&nstep, 1, 1);
	pvm_upkreal(&rsize, 1, 1);
	pvm_upkreal(&slack, 1, 1);
	pvm_upkint(&nbody, 1, 1);
	bodytab = (bodyptr) allocate(nbody * sizeof(body));
	pvm_upkbyte((char *) bodytab, nbody * sizeof(body), 1);
    }
}

/*
 * STEPSYSTEM: advance N-body system using simple leap-frog.
 */

local void stepsystem(void)
{
    bodyptr p1, p2, p;

    p1 = bodytab + MAX(nstatic, 0);		/* set dynamic body range   */
    p2 = bodytab + nbody + MIN(nstatic, 0);
    for (p = p1; p < p2; p++) {			/* loop over body range	    */
	ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);	/* advance v by 1/2 step    */
	ADDMULVS(Pos(p), Vel(p), dtime);	/* advance r by 1 step      */
    }
    treeforce();
    for (p = p1; p < p2; p++) {			/* loop over body range	    */
	ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);  /* advance v by 1/2 step    */
    }
    nstep++;					/* count another time step  */
    tnow = tnow + dtime;			/* finally, advance time    */
}

/*
 * TREEFORCE: parallel force calculation.
 */

local void treeforce(void)
{
    struct timeval t1, t2;
    bodyptr p;
    int j, i, ntot;
    real speedtot, speedj;

    (void) gettimeofday(&t1, NULL);		/* start wall-clock timer   */
    p = bodytab;				/* scan through body array  */
    for (j = 0; j < ntask; j++)			/* loop over all tasks      */
        for (i = 0; i < nload[j]; i++) {	/* loop over each workload  */
	    Update(p) = (j == jtask ? TRUE : FALSE);
						/* flag bodies to update    */
	    p = NextBody(p);			/* and go on to next body   */
	}
    maketree(bodytab, nbody);			/* construct tree structure */
    gravcalc();					/* compute grav. forces     */
    (void) gettimeofday(&t2, NULL);		/* stop wall-clock timer    */
    twall[jtask] = (t2.tv_sec -  t1.tv_sec) + (t2.tv_usec - t1.tv_usec)/1.0e6;
						/* and store elapsed time   */
    sharestats();				/* exchange timing data     */
    sharegravity();				/* exchange computed forces */
    if (jtask == 0)				/* if process is leader     */
        forcereport();				/* print force statistics   */
    speedtot = slack * nload[0] / twall[0];	/* init sum to lead's speed */
    for (j = 1; j < ntask; j++)			/* loop over all helpers    */
        speedtot += nload[j] / twall[j];	/* add up speed of each     */
    ntot = 0;					/* init sum of workloads    */
    for (j = 0; j < ntask; j++) {		/* loop over all tasks      */
	speedj = (j == 0 ? slack : 1) *  nload[j] / twall[j];
						/* estimate speed of each   */
        nload[j] = (nload[j] + (int) rfloor(nbody * speedj / speedtot)) / 2;
						/* scale load by rel. speed */
	ntot += nload[j];			/* sum overall workload     */
    }
    nload[ntask-1] += nbody - ntot;		/* last task gets remainder */
    if (nload[ntask-1] < 0)			/* check for unlikely error */
	fatal("%s: load balance error\n", getargv0());
}

/*
 * SHARESTATS: broadcast and receive force calculation statistics.
 */

local void sharestats(void)
{
    int j, nfcj, nbbj, nbcj;
    real cpufj;

    if (jtask == 0) {				/* if this task handles I/O */
        printf("\n\t%8s%8s%12s%12s%8s%8s\n",
	       "task", "nfcalc", "nbbcalc", "nbccalc", "twall", "tcpu");
	printf("\t%8d%8d%12d%12d%8.3f%8.3f\n",
	       jtask, nfcalc, nbbcalc, nbccalc, twall[jtask] / 60, cpuforce);
    }
    pvm_initsend(PvmDataRaw);			/* start statistics message */
    pvm_pkint(&nfcalc, 1, 1);			/* count forces computed    */
    pvm_pkint(&nbbcalc, 1, 1);			/* count body-body inter.   */
    pvm_pkint(&nbccalc, 1, 1);			/* count body-cell inter.   */
    pvm_pkreal(&twall[jtask], 1, 1);		/* elapsed wall-clock time  */
    pvm_pkreal(&cpuforce, 1, 1);		/* CPU time consumed        */
    pvm_mcast(taskid, ntask, GravStats);	/* send statistics message  */
    for (j = 0; j < ntask; j++)			/* loop over all tasks      */
        if (j != jtask) {			/* which are not this task  */
	    pvm_recv(taskid[j], GravStats);	/* get statistics message   */
	    pvm_upkint(&nfcj, 1, 1);		/* get forces computed.     */
	    pvm_upkint(&nbbj, 1, 1);		/* get body-body inter.     */
	    pvm_upkint(&nbcj, 1, 1);		/* get body-cell inter.     */
	    pvm_upkreal(&twall[j], 1, 1);	/* get wall-clock time      */
	    pvm_upkreal(&cpufj, 1, 1);		/* get CPU time consumed    */
	    nfcalc += nfcj;			/* sum force calculations   */
	    nbbcalc += nbbj;			/* sum body-body interact.  */
	    nbccalc += nbcj;			/* sum body-cell interact.  */
	    cpuforce += cpufj;			/* sum CPU time consumed    */
	    if (jtask == 0)			/* if this task handles I/O */
	        printf("\t%8d%8d%12d%12d%8.3f%8.3f\n",
		       j, nfcj, nbbj, nbcj, twall[j] / 60, cpufj);
	}
}

/*
 * SHAREGRAVITY: broadcast and receive gravity data.
 */

local void sharegravity(void)
{
    bodyptr p;
    int i, j, npot;

    p = bodytab;				/* scan throgh body array   */
    for (j = 0; j < ntask; j++)			/* loop over all tasks      */
        if (j == jtask) {			/* if this task is current  */
	    pvm_initsend(PvmDataRaw);		/* prepare to send data     */
	    for (i = 0; i < nload[j]; i++) {	/* loop over workload       */
		pvm_pkreal(Acc(p), NDIM, 1);	/* pack acceleration data   */
		p = NextBody(p);		/* and go on to next body   */
	    }
	    pvm_mcast(taskid, ntask, AccData);	/* now send gravity data    */
	} else {				/* another task goes next   */
	    pvm_recv(taskid[j], AccData);	/* receive its gravity data */
	    for (i = 0; i < nload[j]; i++) {	/* loop over its workload   */
		pvm_upkreal(Acc(p), NDIM, 1);	/* unpack acceleration data */
		p = NextBody(p);		/* and go on to next body   */
	    }
	}
    p = NthBody(bodytab, nload[0]);		/* skip leader's bodies     */
    if (jtask == 0) {				/* if this task is leader   */
	for (j = 1; j < ntask; j++) {		/* loop over all helpers    */
	    pvm_initsend(PvmDataRaw);		/* prepare to prompt helper */
	    pvm_send(taskid[j], PhiPrompt);	/* send prompts in order    */
	    pvm_recv(taskid[j], PhiData);	/* get phi data from helper */
	    for (i = 0; i < nload[j]; i++) {	/* loop over the workload   */
	        pvm_upkreal(&Phi(p), 1, 1);	/* upack potential data     */
		p = NextBody(p);		/* step on to next body     */
	    }
	}
    } else {					/* if this task is a helper */
	pvm_recv(taskid[0], PhiPrompt);		/* await prompt from leader */
        pvm_initsend(PvmDataRaw);		/* prepare to send data     */
	while (p < bodytab+nbody) {		/* loop over all bodies     */
	    if (Update(p))			/* if this body was updated */
		pvm_pkreal(&Phi(p), 1, 1);	/* pack its potential data  */
	    p = NextBody(p);			/* and go on to next body   */
	}
	pvm_send(taskid[0], PhiData);		/* send message to leader   */
    }
}

/*
 * STOPRUN: shutdown parallel N-body code.
 */

local void stoprun(void)
{
    pvm_exit();
}
