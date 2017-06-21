/*
 * sphcode.h: define global variables for sphcode.c and sphio.c.
 * Copyright (c) 2016 by Joshua E. Barnes, Honolulu, Hawaii.
 */

#ifndef _sphcode_h
#define _sphcode_h

#include "sphdefs.h"

//  Parameters, state variables, and diagnostics for SPH integration.

global string infile;			// file name for snapshot input
global string outpatn;			// file pattern for snapshot output
global string savepatn;			// file pattern for state output
global string restfile;			// file name for state input
global string strmpatn;			// file pattern for simulation stream
global string logfile;			// name for log output file
global stream logstream;		// open stream for log output

global real gamma0;			// ratio of specific heats
global real uradpk;			// uint at peak of cooling curve
global real lambpk;			// peak cooling rate at rho = 1
global real uintmax;			// maximum internal energy (if > 0)
global real sigmastar;			// Stefan-Boltzmann parameter
global real opacity;			// radiation opacity parameter
global real conduct;			// heat conduction parameter
global real cstar[3];			// star formation law constants
global real nstar[3];			// star formation density indicies
global real mstar[3];			// star formation shock indicies
global real tau_ml;			// timescale for mass-loss to start
global real beta_ml;			// index for mass-loss power-law
global real alpha;			// bulk viscosity parameter
global real beta;			// vN-R viscosity parameter
global int nsmooth;			// bodies in smoothing volume
global int nbucket;			// bodies in leaves of kd tree
global real slope0;			// slope of kernel at origin
global real dtime;			// basic integration timestep
global real courant;			// Courant condition parameter
global real fdrag;			// velocity damping parameter
global gsprof *gravgsp;			// GSP for external grav. field
global string outputs;			// additional fields to output
global real tstop;			// time to stop calculation
global real dtout;			// data output timestep
global string headline;			// message identifying program
global int nstep;			// number of time-steps
global int levmax;			// max. level in multistep scheme
global real eradiate;			// energy lost/gained via radiation
global real tnow;			// current value of time
global real tout;			// time of next output
global int nbody;			// total number of bodies in system
global int ngas;			// number of SPH bodies in system
global bodyptr btab;			// pointer to array of bodies

//  Prototypes for I/O routines.

void inputdata(void);			// read initial data file
void inputgrav(string);			// read external GSP file
void startoutput(string *opts);		// open streams, print header
void outputhead(void);			// announce start of time-step
void outputdata(void);			// perform output operation
void savestate(void);			// save system state	
void restorestate(void);		// restore system state

#endif // ! _sphcode_h
