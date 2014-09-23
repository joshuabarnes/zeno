/****************************************************************************/
/* TREEPEEK.C: peek at contents of treecode state file.                     */
/* Copyright (c) 2000 by Joshua E. Barnes, Santa Barbara, CA.               */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "filestruct.h"
#define global					/* don't default to extern  */
#include "treecode.h"
#include "fixbody.h"

#include <sys/types.h>
#include <sys/stat.h>

string defv[] = {		";Peek at contents of state file",
    "in=???",			";Treecode simulation state file",
    "out=???",			";Output file with N-body frame",
    "outputs=" PosTag "," VelTag,
				";Data fields to output",
    "VERSION=1.1",		";Joshua Barnes  24 November 2000",
    NULL,
};

/* Prototypes for local procedures. */

local void setupbody(void);			/* set offsets for dynbody  */

/*
 * MAIN: toplevel routine for hierarchical N-body code.
 */

int main(int argc, string argv[])
{
    initparam(argv, defv);			/* initialize param access  */
    setupbody();				/* interface to dynbody     */
    restorestate(getparam("in"));
    output();
    return (0);					/* end with proper status   */
}

/*
 * SETUPBODY: inform dynamic body routines of relevant body fields.
 */

local void setupbody(void)
{
    define_body(sizeof(body), Precision, NDIM);
    define_body_offset(PosTag,  BodyOffset(Pos));
    define_body_offset(VelTag,  BodyOffset(Vel));
    define_body_offset(MassTag, BodyOffset(Mass));
    define_body_offset(PhiTag,  BodyOffset(Phi));
    define_body_offset(AccTag,  BodyOffset(Acc));
}

/*
 * OUTPUT: write out snapshot file.
 */

void output(void)
{
    string *outputs;
    char namebuf[256];
    struct stat buf;
    stream outstr;

    outputs = burststring(getparam("outputs"), ", ");
    sprintf(namebuf, getparam("out"), nstep);	/* make up output file name */
    if (stat(namebuf, &buf) != 0) {		/* no output file exists?   */
        outstr = stropen(namebuf, "w");         /* create & open for output */
	put_history(outstr);			/* write out hiatory data   */
    } else					/* else file already exists */
        outstr = stropen(namebuf, "a");		/* reopen and append output */
    put_snap(outstr, &bodytab, &nbody, &tnow, outputs);
    strclose(outstr);				/* close up output file     */
}

/*
 * RESTORESTATE: restore state from disk file.
 */

void restorestate(string file)
{
    stream str;
    string program, version;
    real freq;

    str = stropen(file, "r");			/* open state input file    */
    program = get_string(str, "program");
    version = get_string(str, "version");
    eprintf("[%s: reading state file from %s (v%s)]\n",
	    getargv0(), program, version);
    headline = get_string(str, "headline");
    if (streq(version, "1.3")) {
        get_data(str, "freq", RealType, &freq, 0);
        dtime = 1.0 / freq;
    } else
        get_data(str, "dtime", RealType, &dtime, 0);
    get_data(str, "nstatic", IntType, &nstatic, 0);
    if (get_tag_ok(str, "theta"))
        get_data(str, "theta", RealType, &theta, 0);
    get_data(str, "usequad", BoolType, &usequad, 0);
    get_data(str, "eps", RealType, &eps, 0);
    options = get_string(str, "options");
    if (streq(version, "1.3"))
        outputs = PosTag "," VelTag;
    else
        outputs = get_string(str, "outputs");
    get_data(str, "tstop", RealType, &tstop, 0);
    if (streq(version, "1.3")) {
        get_data(str, "freqout", RealType, &freq, 0);
        dtout = 1.0 / freq;
    } else
        get_data(str, "dtout", RealType, &dtout, 0);
    get_data(str, "tnow", RealType, &tnow, 0);
    get_data(str, "tout", RealType, &tout, 0);
    get_data(str, "nstep", IntType, &nstep, 0);
    get_data(str, "rsize", RealType, &rsize, 0);
    get_data(str, "nbody", IntType, &nbody, 0);
    bodytab = (bodyptr) allocate(nbody * sizeof(body));
    get_data(str, "bodytab", AnyType, bodytab, nbody, sizeof(body), 0);
    strclose(str);
    eprintf("[%s: read %d bodies at time %f, step %d]\n",
	    getargv0(), nbody, tnow, nstep);
}
