/*
 * vectmath.h: include file for vector/matrix operations.
 */

#ifndef _vectmath_h
#define _vectmath_h

#include "vectdefs.h"

//  __________________
//  Vector operations.

#define CLRV(v)			/* CLeaR Vector */			\
{									\
    int _i;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = 0.0;							\
}

#define UNITV(v,j)		/* UNIT Vector */			\
{									\
    int _i;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = (_i == (j) ? 1.0 : 0.0);				\
}

#define SETV(v,u)		/* SET Vector */			\
{ 									\
    int _i; 								\
    for (_i = 0; _i < NDIM; _i++) 					\
	(v)[_i] = (u)[_i]; 						\
}

#ifdef THREEDIM

#define DOTVP(s,v,u)		/* DOT Vector Product */		\
{									\
    (s) = (v)[0]*(u)[0] + (v)[1]*(u)[1] + (v)[2]*(u)[2];		\
}

#else

#define DOTVP(s,v,u)		/* DOT Vector Product */		\
{									\
    int _i;								\
    (s) = 0.0;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(s) += (v)[_i] * (u)[_i];					\
}

#endif

#define ABSV(s,v)		/* ABSolute value of a Vector */	\
{									\
    real _tmp;		                                                \
    int _i;								\
    _tmp = 0.0;								\
    for (_i = 0; _i < NDIM; _i++)					\
	_tmp += (v)[_i] * (v)[_i];					\
    (s) = rsqrt(_tmp);                                                  \
}

#ifdef THREEDIM

#define ADDV(v,u,w)		/* ADD Vector */			\
{									\
    (v)[0] = (u)[0] + (w)[0];						\
    (v)[1] = (u)[1] + (w)[1];						\
    (v)[2] = (u)[2] + (w)[2];						\
}

#define SUBV(v,u,w)		/* SUBtract Vector */			\
{									\
    (v)[0] = (u)[0] - (w)[0];						\
    (v)[1] = (u)[1] - (w)[1];						\
    (v)[2] = (u)[2] - (w)[2];						\
}

#define MULVS(v,u,s)		/* MULtiply Vector by Scalar */		\
{									\
    (v)[0] = (u)[0] * (s);						\
    (v)[1] = (u)[1] * (s);						\
    (v)[2] = (u)[2] * (s);						\
}

#define DIVVS(v,u,s)		/* DIVide Vector by Scalar */		\
{									\
    (v)[0] = (u)[0] / (s);						\
    (v)[1] = (u)[1] / (s);						\
    (v)[2] = (u)[2] / (s);						\
}

#else

#define ADDV(v,u,w)		/* ADD Vector */			\
{									\
    int _i;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = (u)[_i] + (w)[_i];					\
}

#define SUBV(v,u,w)		/* SUBtract Vector */			\
{									\
    int _i;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = (u)[_i] - (w)[_i];					\
}

#define MULVS(v,u,s)		/* MULtiply Vector by Scalar */		\
{									\
    int _i;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = (u)[_i] * (s);					\
}

#define DIVVS(v,u,s)		/* DIVide Vector by Scalar */		\
{									\
    int _i;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = (u)[_i] / (s);					\
}

#endif

#define DISTV(s,u,v)		/* DISTance between Vectors */		\
{									\
    real _tmp;                                                		\
    int _i;								\
    _tmp = 0.0;								\
    for (_i = 0; _i < NDIM; _i++)					\
	_tmp += ((u)[_i]-(v)[_i]) * ((u)[_i]-(v)[_i]);		        \
    (s) = rsqrt(_tmp);                                                  \
}

#ifdef THREEDIM

#define DISTSQV(s,u,v)		/* DISTance SQuared between Vectors */	\
{									\
    real _dx, _dy, _dz;							\
    _dx = (u)[0] - (v)[0];						\
    _dy = (u)[1] - (v)[1];						\
    _dz = (u)[2] - (v)[2];						\
    (s) = _dx*_dx + _dy*_dy + _dz*_dz;					\
}

#else

#define DISTSQV(s,u,v)		/* DISTance SQuared between Vectors */	\
{									\
    int _i;								\
    (s) = 0.0;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(s) += ((u)[_i]-(v)[_i]) * ((u)[_i]-(v)[_i]);		        \
}

#endif

#ifdef TWODIM

#define CROSSVP(s,v,u)		/* CROSS Vector Product */		\
{									\
    (s) = (v)[0]*(u)[1] - (v)[1]*(u)[0];				\
}

#endif

#ifdef THREEDIM

#define CROSSVP(v,u,w)		/* CROSS Vector Product */		\
{									\
    (v)[0] = (u)[1]*(w)[2] - (u)[2]*(w)[1];				\
    (v)[1] = (u)[2]*(w)[0] - (u)[0]*(w)[2];				\
    (v)[2] = (u)[0]*(w)[1] - (u)[1]*(w)[0];				\
}

#endif

//  __________________
//  Matrix operations.

#define CLRM(p)			/* CLeaR Matrix */			\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = 0.0;						\
}

#define SETMI(p)		/* SET Matrix to Identity */		\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (_i == _j ? 1.0 : 0.0);			\
}

#define SETM(p,q)		/* SET Matrix */			\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j];					\
}

#define TRANM(p,q)		/* TRANspose Matrix */			\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_j][_i];					\
}

#define ADDM(p,q,r)		/* ADD Matrix */			\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j] + (r)[_i][_j];			\
}

#define SUBM(p,q,r)		/* SUBtract Matrix */			\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j] - (r)[_i][_j];			\
}

#define MULM(p,q,r)		/* Multiply Matrix */			\
{									\
    int _i, _j, _k;							\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++) {					\
	    (p)[_i][_j] = 0.0;						\
	    for (_k = 0; _k < NDIM; _k++)				\
		(p)[_i][_j] += (q)[_i][_k] * (r)[_k][_j];		\
	}								\
}

#define MULMS(p,q,s)		/* MULtiply Matrix by Scalar */		\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)				        \
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j] * (s);				\
}

#define DIVMS(p,q,s)		/* DIVide Matrix by Scalar */		\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (q)[_i][_j] / (s);				\
}

#define MULMV(v,p,u)		/* MULtiply Matrix by Vector */		\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++) {					\
	(v)[_i] = 0.0;							\
	for (_j = 0; _j < NDIM; _j++)					\
	    (v)[_i] += (p)[_i][_j] * (u)[_j];				\
    }									\
}

#define MULMV1(v,p,u)		/* MULtiply Matrix by Vector */		\
  {				/* args v, u may be identical */	\
    int _i, _j;								\
    vector _tmp;							\
    for (_i = 0; _i < NDIM; _i++) {					\
	_tmp[_i] = 0.0;							\
	for (_j = 0; _j < NDIM; _j++)					\
	    _tmp[_i] += (p)[_i][_j] * (u)[_j];				\
    }									\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = _tmp[_i];						\
}

#define OUTVP(p,v,u)		/* OUTer Vector Product */		\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (v)[_i] * (u)[_j];				\
}

#define TRACEM(s,p)		/* TRACE of Matrix */			\
{									\
    int _i;								\
    (s) = 0.0;								\
    for (_i = 0.0; _i < NDIM; _i++)					\
	(s) += (p)[_i][_i];						\
}

//  ____________________________
//  Enhancements for tree codes.

#ifdef THREEDIM

#define DOTPSUBV(s,v,u,w)	/* SUB Vectors, form DOT Prod */	\
{									\
    (v)[0] = (u)[0] - (w)[0];    (s)  = (v)[0] * (v)[0];		\
    (v)[1] = (u)[1] - (w)[1];    (s) += (v)[1] * (v)[1];		\
    (v)[2] = (u)[2] - (w)[2];    (s) += (v)[2] * (v)[2];		\
}

#define DOTPMULMV(s,v,p,u)	/* MUL Mat by Vect, form DOT Prod */	\
{									\
    DOTVP(v[0], p[0], u);    (s)  = (v)[0] * (u)[0];			\
    DOTVP(v[1], p[1], u);    (s) += (v)[1] * (u)[1];			\
    DOTVP(v[2], p[2], u);    (s) += (v)[2] * (u)[2];			\
}

#define ADDMULVS(v,u,s)		/* MUL Vect by Scalar, ADD to vect */	\
{									\
    (v)[0] += (u)[0] * (s);						\
    (v)[1] += (u)[1] * (s);						\
    (v)[2] += (u)[2] * (s);						\
}

#define ADDMULVS2(v,u,s,w,r)	/* 2 times MUL V by S, ADD to vect */	\
{									\
    (v)[0] += (u)[0] * (s) + (w)[0] * (r);				\
    (v)[1] += (u)[1] * (s) + (w)[1] * (r);				\
    (v)[2] += (u)[2] * (s) + (w)[2] * (r);				\
}

#endif

//  ________________________
//  Misc. impure operations.

#define SETVS(v,s)		/* SET Vector to Scalar */		\
{									\
    int _i;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = (s);							\
}

#define ADDVS(v,u,s)		/* ADD Vector and Scalar */		\
{									\
    int _i;								\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = (u)[_i] + (s);					\
}

#define SETMS(p,s)		/* SET Matrix to Scalar */		\
{									\
    int _i, _j;								\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (p)[_i][_j] = (s);						\
}

#endif  // ! _vectmath_h
