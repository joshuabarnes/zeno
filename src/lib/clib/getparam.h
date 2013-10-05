/*
 * getparam.h: include file for command-line processing.
 */

#ifndef _getparam_h
#define  _getparam_h

void initparam(string *argv, string *defv);	// initialize param package

string getparam(string name);			// return value as string

int getiparam(string name);			// return value as int

bool getbparam(string name);			// return value as boolean

double getdparam(string name);			// return value as double

//  Macros to obtain name and version of program.
//  _____________________________________________

#define getprog()     (getparam("argv0"))	// return name of program
#define getargv0()    (getparam("argv0"))	// return name of program

#define getversion()  (getparam("VERSION"))	// return version as string

//  Function to inquire about parameter status.
//  ___________________________________________

int getparamstat(string name);			// return status flags

//  Bits for parameter status flags.
//  ________________________________

#define DEFPARAM	001			// param has default value

#define REQPARAM	002			// must be given a value

#define ARGPARAM	004			// value reset by argument

#define INPARAM		010			// param used for input

#define OUTPARAM	020			// param used for output

#endif  // ! _getparam_h
