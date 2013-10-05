/*
 * SNAPASCII.C: convert binary N-body file to ASCII format.
 *         "I did it my way" --  F. Sinatra.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "bodytags.h"

#ifndef TIPSY
#  define PURPOSE		";Convert binary N-body file to ascii"
#  define RFORMAT		"%21.14E"
#else
#  define PURPOSE		";Convert binary N-body file to tipsy"
#  define RFORMAT		"%14.7E"
#endif

string defv[] = {		PURPOSE,
    "in=???",                   ";Input file name",
    "out=???",                  ";Output file name",
    "options=" MassTag "," PosTag "," VelTag,
				";Others are " PhiTag "," AuxTag "," KeyTag,
    "iformat=  %d",		";Output format for integers",
    "rformat= " RFORMAT,	";Output format for floating-point",
    "VERSION=2.3",		";Josh Barnes  11 June 1998",
    NULL,
};

local string options, iformat, rformat;

local stream instr, outstr;

local void snapascii(void);
local void out_header(int, real);
local void out_idata(int *, int);
local void out_rdata(real *, int);
local void out_vdata(real *, int);
local void out_phase(real *, int);
local void out_int(int);
local void out_3int(int, int, int);
local void out_real(real);
local void out_vect(vector);

int main(int argc, string argv[])
{
    initparam(argv, defv);
    instr = stropen(getparam("in"), "r");
    get_history(instr);
    outstr = stropen(getparam("out"), "w");
    options = getparam("options");
    iformat = getparam("iformat");
    rformat = getparam("rformat");
    snapascii();
    fflush(outstr);
    return (0);
}

local void snapascii(void)
{
    int nbody, *ibuf;
    real tsnap, *rbuf, *vbuf;

    while (get_tag_ok(instr, SnapShotTag)) {	/* loop reading data frames */
	get_set(instr, SnapShotTag);
	get_set(instr, ParametersTag);
	if (get_tag_ok(instr, TimeTag))
	    get_data(instr, TimeTag, RealType, &tsnap, 0);
	else
	    tsnap = 0.0;
	get_data(instr, NBodyTag, IntType, &nbody, 0);
	get_tes(instr, ParametersTag);
	out_header(nbody, tsnap);
	get_set(instr, ParticlesTag);
        if (scanopt(options, MassTag)) {	/* mass data to convert?    */
	    rbuf = (real *) allocate(nbody * sizeof(real));
	    get_data(instr, MassTag, RealType, rbuf, nbody, 0);
	    out_rdata(rbuf, nbody);
	    free(rbuf);
        }
        if (scanopt(options, PosTag)) {		/* position data?           */
	    vbuf = (real *) allocate(nbody * NDIM * sizeof(real));
	    get_data(instr, PosTag, RealType, vbuf, nbody, NDIM, 0);
	    out_vdata(vbuf, nbody);
	    free(vbuf);
	}
        if (scanopt(options, VelTag)) {		/* velocity data?           */
	    vbuf = (real *) allocate(nbody * NDIM * sizeof(real));
	    get_data(instr, VelTag, RealType, vbuf, nbody, NDIM, 0);
	    out_vdata(vbuf, nbody);
	    free(vbuf);
	}
        if (scanopt(options, PhaseTag)) {	/* phase-space data?        */
	    vbuf = (real *) allocate(nbody * 2 * NDIM * sizeof(real));
	    get_data(instr, PhaseTag, RealType, vbuf, nbody, 2, NDIM, 0);
	    out_vdata(vbuf, nbody);
	    free(vbuf);
	}
        if (scanopt(options, PhiTag)) {		/* potential data?          */
	    rbuf = (real *) allocate(nbody * sizeof(real));
	    get_data(instr, PhiTag, RealType, rbuf, nbody, 0);
	    out_rdata(rbuf, nbody);
	    free(rbuf);
        }
        if (scanopt(options, AuxTag)) {		/* auxiliary data?          */
	    rbuf = (real *) allocate(nbody * sizeof(real));
	    get_data(instr, AuxTag, RealType, rbuf, nbody, 0);
	    out_rdata(rbuf, nbody);
	    free(rbuf);
        }
        if (scanopt(options, KeyTag)) {		/* key data?                */
	    ibuf = (int *) allocate(nbody * sizeof(int));
	    get_data(instr, KeyTag, IntType, ibuf, nbody, 0);
	    out_idata(ibuf, nbody);
	    free(ibuf);
        }
	get_tes(instr, ParticlesTag);
	get_tes(instr, SnapShotTag);
    }
}

local void out_header(int nbody, real tsnap)
{
    eprintf("[%s: writing %d bodies at time %f]\n", getargv0(), nbody, tsnap);
#ifndef TIPSY
    out_int(nbody);
#else
    out_3int(nbody, 0, 0);
#endif
    out_int(NDIM);
    out_real(tsnap);
}

local void out_idata(int *ip, int nbody)
{
    int i;

    for (i = 0; i < nbody; i++) {
	out_int(*ip);
	ip++;
    }
}

local void out_rdata(real *rp, int nbody)
{
    int i;

    for (i = 0; i < nbody; i++) {
	out_real(*rp);
	rp++;
    }
}

local void out_vdata(real *cptr, int nbody)
{
    real *cp;
    int i, k;

#ifndef TIPSY
    cp = cptr;
    for (i = 0; i < nbody; i++) {
	out_vect(cp);
	cp += NDIM;
    }
#else
    for (k = 0; k < NDIM; k++) {
	cp = cptr + k;
	for (i = 0; i < nbody; i++) {
	    out_real(*cp);
	    cp += NDIM;
	}
    }
#endif
}

local void out_phase(real *cptr, int nbody)
{
    real *cp;
    int k, i;

#ifndef TIPSY
    for (k = 0; k <= 1; k++) {
	cp = cptr + NDIM * k;
	for (i = 0; i < nbody; i++) {
	    out_vect(cp);
	    cp += 2 * NDIM;
	}
    }
#else
    for (k = 0; k < 6; k++) {
	cp = cptr + k;
	for (i = 0; i < nbody; i++) {
	    out_real(*cp);
	    cp += 2 * NDIM;
	}
    }
#endif
}

#define MAXBUF  256

local void out_int(int ival)
{
    char fmtbuf[MAXBUF];

    sprintf(fmtbuf, "%s\n", iformat);
    fprintf(outstr, fmtbuf, ival);
}

local void out_3int(int i1, int i2, int i3)
{
    char fmtbuf[MAXBUF];

    sprintf(fmtbuf, "%s%s%s\n", iformat, iformat, iformat);
    fprintf(outstr, fmtbuf, i1, i2, i3);
}

local void out_real(real rval)
{
    char fmtbuf[MAXBUF];

    sprintf(fmtbuf, "%s\n", rformat);
    fprintf(outstr, fmtbuf, rval);
}

local void out_vect(vector vec)
{
    char fmtbuf[MAXBUF];

    sprintf(fmtbuf, "%s%s%s\n", rformat, rformat, rformat);
    fprintf(outstr, fmtbuf, vec[0], vec[1], vec[2]);
}
