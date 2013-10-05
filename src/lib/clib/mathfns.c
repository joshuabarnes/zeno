/*
 * mathfns.c: utility routines for various sorts of math operations.
 * Most these functions work with real values, meaning that they can
 * handle either floats or doubles, depending on compiler switches.
 */

#include "stdinc.h"
#include "mathfns.h"

#ifndef TESTBED

//  __________________________________
//  rsqr, rqbe: compute x*x and x*x*x.

real rsqr(real x)
{
    return (x * x);
}

real rqbe(real x)
{
    return (x * x * x);
}

#if defined(SINGLEPREC)

//  __________________________
//  fcbrt: floating cube root.

float fcbrt(float x)
{
    return ((float) cbrt((double) x));
}

#endif

//  ___________________________
//  rdex: inverse log base ten.

real rdex(real x)
{
    return (rexp(M_LN10 * x));
}

#if defined(SINGLEPREC) || !defined(MACOSX)

//  ___________________________________________
//  RLOG2, REXP2: log, inverse log to base two.

real rlog2(real x)
{
    return (rlog(x) / M_LN2);
}

real rexp2(real x)
{
    return (rexp(M_LN2 * x));
}

#endif

#endif // ! TESTBED

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
    "x=2.0",
    "y=1.0",
    "VERSION=0",
    NULL,
};

main(int argc, string argv[])
{
    real x, y;

    initparam(argv, defv);
    x = getdparam("x");
    y = getdparam("y");
#if defined(MIXEDPREC)
    printf("Mixed Precision Version\n");
#endif
#if defined(SINGLEPREC)
    printf("Single Precision Version\n");
#endif
#if defined(DOUBLEPREC)
    printf("Double Precision Version\n");
#endif
    printf("rsqrt(%4.2f)        = %12.8f\n", x, rsqrt(x));
    printf("rcbrt(%4.2f)        = %12.8f\n", x, rcbrt(x));
    printf("rsqr(%4.2f)         = %12.8f\n", x, rsqr(x));
    printf("rqbe(%4.2f)         = %12.8f\n", x, rqbe(x));
    printf("rsin(%4.2f)         = %12.8f\n", x, rsin(x));
    printf("rcos(%4.2f)         = %12.8f\n", x, rcos(x));
    printf("rtan(%4.2f)         = %12.8f\n", x, rtan(x));
    printf("rasin(%4.2f)        = %12.8f\n", x, rasin(x));
    printf("racos(%4.2f)        = %12.8f\n", x, racos(x));
    printf("ratan(%4.2f)        = %12.8f\n", x, ratan(x));
    printf("ratan2(%4.2f,%5.2f) = %12.8f\n", x, y, ratan2(x, y));
    printf("rlog(%4.2f)         = %12.8f\n", x, rlog(x));
    printf("rexp(%4.2f)         = %12.8f\n", x, rexp(x));
    printf("rlog2(%4.2f)        = %12.8f\n", x, rlog2(x));
    printf("rexp2(%4.2f)        = %12.8f\n", x, rexp2(x));
    printf("rlog10(%4.2f)       = %12.8f\n", x, rlog10(x));
    printf("rdex(%4.2f)         = %12.8f\n", x, rdex(x));
    printf("rsinh(%4.2f)        = %12.8f\n", x, rsinh(x));
    printf("rcosh(%4.2f)        = %12.8f\n", x, rcosh(x));
    printf("rtanh(%4.2f)        = %12.8f\n", x, rtanh(x));
    printf("rpow(%4.2f,%5.2f)   = %12.8f\n", x, y, rpow(x, y));
    printf("rabs(%4.2f)         = %12.8f\n", x, rabs(x));
    printf("rfloor(%4.2f)       = %12.8f\n", x, rfloor(x));
    printf("rceil(%4.2f)        = %12.8f\n", x, rceil(x));
}

#endif
