/*
 * buildmap.h: interface to buildmap routine.
 */

#ifndef _buildmap_h
#define _buildmap_h

//  buildmap: generate extendbody, computemap, and computetime routines;
//  compile and link with snapmap main program.
//  ____________________________________________________________________

void buildmap(string prog,		// name of executable
	      string *names,		// access macros for body variables
	      string *exprs,		// C expressions for body variables
	      string *types,		// types for new varaibles, if any
	      string tmap,		// C expression for snapshot time
	      string prec,		// numerical precision 
	      int ndim,			// number of space dimensions
	      bool cleanup);		// if TRUE, delete prog.c file

//  getmapdefs: return pointer to mapdefs table.
//  ____________________________________________

string *getmapdefs(void);

//  SNAPMAP_BODY_VARS: list of variables and parameters which can appear in
//  expressions for body components; formatted for use in defv[] string.
//  _______________________________________________________________________

#define SNAPMAP_BODY_VARS  ";x,y,z,vx,vy,vz,ax,ay,az,m,phi,smooth,rho,", \
			   ";entf,uint,udot,udotrad,udotvis,tau,type,", \
                           ";birth,death,key,keyarr,auxvx,auxvy,auxvz,", \
			   ";auxarr,r,R,v,vr,vt,etot,jx,jy,yz,jtot;t,i,n"

//  SNAPMAP_TIME_VARS: list of parameters which can appear in expression
//  for snapshot time; formatted for use in defv[] string.

#define SNAPMAP_TIME_VARS  "t,n"

#endif /* ! _buildmap_h */
