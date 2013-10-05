/*
 * EXTFORCE.C: interpolate axisymmetric field.
 * Defines: initextfield(), extforce().
 */

#include <stdinc.h>
#include <mathfns.h>
#include <getparam.h>
#include <filestruct.h>

/*
 * ARAD, AVER, PHI: tables of field values, evaluated at points
 *
 *	r = scale * ir*ir,	z = scale * iz*iz,
 *
 * for index values 0 <= ir,iz < NTAB.
 */

#define	NTAB	65		/* size of field tabulation array */

local real scale;		/* spacing of finest grid-point */

local real arad[NTAB][NTAB];	/* horizontal acceleration */

local real aver[NTAB][NTAB];	/* vertical acceleration */

local real phi[NTAB][NTAB];	/* gravitational potential */

/*
 * INITEXTFIELD: read force-field tables.
 */

void initextfield(string name)
{
    stream str;
    int  ntab;

    str = stropen(name, "r");
    get_history(str);
    get_set(str, "DiskTable");
    get_data(str, "ntab", IntType, &ntab, 0);
    if (ntab != NTAB)
	error("initextfield in %s: expected ntab = %d, got %d\n",
	      getargv0(), NTAB, ntab);
    get_data(str, "scale", RealType, &scale, 0);
    get_data(str, "arad", RealType, arad, ntab, ntab, 0);
    get_data(str, "aver", RealType, aver, ntab, ntab, 0);
    get_data(str, "phi", RealType, phi, ntab, ntab, 0);
    get_tes(str, "DiskTable");
    strclose(str);
}

/*
 * INTERP: linear interpolation macro.
 */

#define INTERP(x,i,j,f,g)						\
		((1 - g) * ((1 - f) * x[ i ][ j ] + f * x[i+1][ j ]) +	\
		    g    * ((1 - f) * x[ i ][j+1] + f * x[i+1][j+1]))

/*
 * EXTFORCE: find force and potential of axysymmetric field.
 */

void extforce(real *ar, real *az, real *ph, real r, real z, 
	      real mext, real rext, real eps)
{
    real r1, z1, fr, fz, ar1, az1, ph1;
    int ir, iz;

    r1 = r / rext;				/* get dimensonless radius */
    z1 = rsqrt(z*z + eps*eps) / rext;		/* and softened height */
    ir = (int) rfloor(rsqrt(r1 / scale));	/* find table indices */
    iz = (int) rfloor(rsqrt(z1 / scale));
    if (ir < NTAB-1 && iz < NTAB-1) {		/* within tabulated range? */
	fr = (r1/scale - ir*ir) / (2*ir + 1);	/*   compute remainders */
	fz = (z1/scale - iz*iz) / (2*iz + 1);
	ar1 = INTERP(arad, ir, iz, fr, fz);	/*   interp. field values */
	az1 = INTERP(aver, ir, iz, fr, fz);
	ph1 = INTERP(phi,  ir, iz, fr, fz);
    } else {					/* use 1/r approximation */
	ph1 = 1.0 / rsqrt(r1*r1 + z1*z1);	/*   compute potential */
	ar1 = r1 * ph1 * ph1 * ph1;		/*   radial acceleration */
	az1 = z1 * ph1 * ph1 * ph1;		/*   vertical acceleration */
    }
    if (z < 0.0)				/* check sign of actual z */
	az1 = - az1;				/*   and invert if needed */
    *ph = - ph1 * mext / rext;			/* scale to physical units */
    *ar = - ar1 * mext / (rext * rext);
    *az = - az1 * mext / (rext * rext);
}

#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "r=1.0",
    "z=1.0",
    "mext=1.0",
    "rext=1.0",
    "eps=0.0",
    "ext=exptab.dat",
    "VERSION=1",
    NULL,
};

int main(int argc, string *argv)
{
    real ar, az, ph;

    initparam(argv, defv);
    initextfield(getparam("ext"));
    extforce(&ar, &az, &ph, getdparam("r"), getdparam("z"),
	     getdparam("mext"), getdparam("rext"), getdparam("eps"));
    printf("a.rad: %12.6f\n", ar);
    printf("a.ver: %12.6f\n", az);
    printf("phi  : %12.6f\n", ph);
    return (0);
}

#endif
