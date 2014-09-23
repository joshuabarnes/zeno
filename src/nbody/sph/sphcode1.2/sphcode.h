/****************************************************************************/
/* SPHCODE.H: define global variables for sphcode.c and sphio.c.            */
/* Copyright (c) 2003 by Joshua E. Barnes, Honolulu, Hawaii.                */
/****************************************************************************/

#ifndef _sphcode_h
#define _sphcode_h

#include "sphdefs.h"

/* Parameters, state variables, and diagnostics for SPH integration. */

#define RANDSIZE 32			/* size of random number state      */
#define ETA2     0.01			/* prevents divergence in viscosity */

global string infile;			/* file name for snapshot input     */
global string outfile;			/* file name for snapshot output    */
global string savefile;			/* file name for state output       */
global string restfile;			/* file name for state input	    */
global string gspfile;			/* file name for extern grav input  */
global real gamma0;			/* ratio of specific heats	    */
global real uradmax;			/* uint at max of cooling curve     */
global real lambmax;			/* max. cooling rate at rho = 1     */
global real sigmastar;			/* Stefan-Boltzmann parameter       */
global real opacity;			/* radiation opacity parameter      */
global real conduct;			/* heat conduction parameter        */
global real starprob;			/* star formation probability	    */
global real rhoindx;			/* power-law index for density	    */
global real udotindx;			/* power-law index for dissipation  */
global real alpha;			/* bulk viscosity parameter	    */
global real beta;			/* vN-R viscosity parameter	    */
global int nsmooth;			/* bodies in smoothing volume	    */
global int nbucket;			/* bodies in leaves of kd tree	    */
global real slope0;			/* slope of kernel at origin	    */
global real dtime;			/* basic integration timestep       */
global real courant;			/* Courant condition parameter      */
global real fdrag;			/* velocity damping parameter	    */
global gsprof *gravgsp;			/* GSP for external grav. field	    */
global string outputs;			/* additional fields to output      */
global real tstop;			/* time to stop calculation         */
global real dtout;			/* data output timestep             */
global string headline;			/* message identifying program	    */
global int nstep;			/* number of time-steps             */
global int levmax;			/* max. level in multistep scheme   */
global real eradiate;			/* energy lost/gained via radiation */
global real tnow;			/* current value of time            */
global real tout;			/* time of next output              */
global char randstate[RANDSIZE];	/* state of random number generator */
global int nbody;			/* total number of bodies in system */
global int ngas;			/* number of SPH bodies in system   */
global bodyptr btab;			/* pointer to array of bodies       */

/* Prototypes for I/O routines. */

void inputdata(void);			/* read initial data file           */
void startoutput(stream, string []);	/* print parameters to log stream   */
void outputhead(stream);		/* announce start of time-step      */
void outputdata(stream);		/* perform output operation         */
void savestate(string);			/* save system state		    */
void restorestate(string);		/* restore system state             */

#endif /* ! _sphcode_h */
