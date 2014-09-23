/*
 * mathfns.h: header file for system and zeno math functions; assumes
 * role of <math.h>.  Defines real-valued synonyms for system
 * functions (eg, rsqrt for square root) and zeno functions (eg,
 * seval), depending on precision switch (MIXEDPREC, SINGLEPREC, or
 * DOUBLEPREC).
 */

#ifndef _mathfns_h
#define _mathfns_h

#include <math.h>

//  System math functions.  Use double-precision versions in both
//  double and mixed precision.

#if defined(DOUBLEPREC) || defined(MIXEDPREC)
#define rsqrt    sqrt
#define rcbrt    cbrt
#define rsin     sin
#define rcos     cos
#define rtan     tan
#define rasin    asin
#define racos    acos
#define ratan    atan
#define ratan2   atan2
#define rlog     log
#define rexp     exp
#define rlog10   log10
#define rsinh    sinh
#define rcosh    cosh
#define rtanh    tanh
#define rpow     pow
#define rabs     fabs
#define rfloor   floor
#define rceil    ceil
#endif

//  System math functions.  Single precision is supplied in two forms,
//  depending on the naming convention.

#if defined(SINGLEPREC)

#if defined(LINUX) || defined(MACOSX)

#define rsqrt    sqrtf
#define rsin     sinf
#define rcos     cosf
#define rtan     tanf
#define rasin    asinf
#define racos    acosf
#define ratan    atanf
#define ratan2   atan2f
#define rlog     logf
#define rexp     expf
#define rlog10   log10f
#define rsinh    sinhf
#define rcosh    coshf
#define rtanh    tanhf
#define rpow     powf
#define rabs     fabsf
#define rfloor   floorf
#define rceil    ceilf

#else

#define rsqrt    fsqrt
#define rsin     fsin
#define rcos     fcos
#define rtan     ftan
#define rasin    fasin
#define racos    facos
#define ratan    fatan
#define ratan2   fatan2
#define rlog     flog
#define rexp     fexp
#define rlog10   log10f
#define rsinh    fsinh
#define rcosh    fcosh
#define rtanh    ftanh
#define rpow     powf
#define rabs     fabsf
#define rfloor   ffloor
#define rceil    fceil

#endif

#endif

//  Functions in mathfns.c; invoked just like those above.

#if defined(DOUBLEPREC) || defined(MIXEDPREC)

#define rsqr     sqr
#define rqbe     qbe
#define rdex     dex
#define rlog2    log2
#define rexp2    exp2

double sqr(double);
double qbe(double);
double dex(double);

#if !defined(MACOSX)
double log2(double);
double exp2(double);
#endif

#else

#define rsqr     fsqr
#define rqbe     fqbe
#define rcbrt    fcbrt
#define rdex     fdex
#define rlog2    flog2
#define rexp2    fexp2

float fsqr(float);
float fqbe(float);
float fcbrt(float);
float fdex(float);
float flog2(float);
float fexp2(float);

#endif

//  Bessel functions available only in double precision.

double bessel_I0(double x);
double bessel_I1(double x);
double bessel_K0(double x);
double bessel_K1(double x);
double bessel_I(double nu, double x);
double bessel_K(double nu, double x);

//  Random number functions return double-precision values.

void init_random(unsigned long seed);
double xrandom(double xl, double xh);
double grandom(double mean, double sdev);
void get_random_state(int *nb, void **st);
void set_random_state(int *nb, void **st);

//  Functions which traffic in pointers to real values must be
//  provided in all three variants.  See also "vectdefs.h".

#if defined(DOUBLEPREC)
#define pickshell dpickshell
#define pickball  dpickball
#define pickbox   dpickbox
#define setrange  dsetrange
#define spline    dspline
#define seval     dseval
#define spldif    dspldif
#endif

#if defined(MIXEDPREC)
#define pickshell mpickshell
#define pickball  mpickball
#define pickbox   mpickbox
#define setrange  msetrange
#define spline    mspline
#define seval     mseval
#define spldif    mspldif
#endif

#if defined(SINGLEPREC)
#define pickshell fpickshell
#define pickball  fpickball
#define pickbox   fpickbox
#define setrange  fsetrange
#define spline    fspline
#define seval     fseval
#define spldif    fspldif
#endif

//  Functions in pickpnt.c

void pickshell(real *, int, real);
void pickball(real *, int, real);
void pickbox(real *, int, real);

//  Functions in spline.c

void spline(real *, real *, real *, int);
real seval(real, real *, real *, real *, int);
real spldif(real, real *, real *, real *, int);

//  Functions in setrange.c

void setrange(real *, string);

#endif  // ! _mathfns_h
