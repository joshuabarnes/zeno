/****************************************************************************/
/* TREECODE.H: define various things for treecode.c and treeio.c.           */
/* Copyright (c) 2000 by Joshua E. Barnes, Honolulu, HI.                    */
/****************************************************************************/

#ifndef _treecode_h
#define _treecode_h

#include "treedefs.h"

/*
 * Parameters, state variables, and diagnostics for N-body integration.
 */

global string infile;			/* file name for snapshot input     */
global string outfile;			/* file pattern for snapshot output */
global string savefile;			/* file pattern for state output    */
global real dtime;			/* basic integration timestep       */
global real dtout;			/* data output timestep             */
global real tstop;			/* time to stop calculation         */
global string outputs;			/* list of data to output           */
global string headline;			/* message describing calculation   */
global real tnow;			/* current value of time            */
global real tout;			/* time of next output              */
global int nstep;			/* number of time-steps             */
global int nbody;			/* number of bodies in system       */
global int nstatic;			/* number of static bodies          */
global bodyptr bodytab;			/* points to array of bodies        */
global int ntrace;			/* total number of trace points     */
global string tracefile;		/* file pattern for trace output    */
global bool testcalc;			/* indicate test-particle calc.     */

/*
 * Prototypes for I/O routines.
 */

void inputdata(void);			/* read initial data file           */
void startoutput(void);			/* open files for output            */
void forcereport(void);			/* report on force calculation      */
void output(void);			/* perform output operation         */
void savestate(string);			/* save system state		    */
void restorestate(string);		/* restore system state             */

#endif /* ! _treecode_h */
