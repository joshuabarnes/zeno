/*
 * FSP.C: routines for operating on finite spherical profiles.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "filestruct.h"
#include "fsp.h"

/*
 * RHO_FSP: evaluate density at given radius.
 */

real rho_fsp(fsprof *fsp, real r)
{
    int n = fsp->npoint - 1;

    if (r < 0)
        error("%s.rho_fsp: undefined for r = %g\n", getargv0(), r);
    if (r < fsp->radius[0])
	return (fsp->density[0] * rpow(r / fsp->radius[0], fsp->alpha));
    else if (r > fsp->radius[n])
	return (fsp->density[n] * rpow(r / fsp->radius[n], fsp->beta));
    if (fsp->dr_coef == NULL) {
	fsp->dr_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
	spline(fsp->dr_coef, fsp->radius, fsp->density, n + 1);
    }
    return (seval(r, fsp->radius, fsp->density, fsp->dr_coef, n + 1));
}

/*
 * DRHO_FSP: evaluate radial derivative of density at given radius.
 */

real drho_fsp(fsprof *fsp, real r)
{
    int n = fsp->npoint - 1;

    if (r < 0)
        error("%s.drho_fsp: undefined for r = %g\n", getargv0(), r);
    if (r < fsp->radius[0])
	return (fsp->density[0] * rpow(r / fsp->radius[0], fsp->alpha) *
		  fsp->alpha / r);
    else if (r > fsp->radius[n])		/* less glitches at large r */
	return (fsp->density[n] * rpow(r / fsp->radius[n], fsp->beta) *
		  fsp->beta / r);
    if (fsp->dr_coef == NULL) {
	fsp->dr_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
	spline(fsp->dr_coef, fsp->radius, fsp->density, n + 1);
    }
    return (spldif(r, fsp->radius, fsp->density, fsp->dr_coef, n + 1));
}

/*
 * MASS_FSP: evaluate enclosed mass at given radius.
 */

real mass_fsp(fsprof *fsp, real r)
{
    int n = fsp->npoint - 1;

    if (r < 0)
        error("%s.mass_fsp: undefined for r = %g\n", getargv0(), r);
    if (r < fsp->radius[0])
	return (fsp->mass[0] * rpow(r / fsp->radius[0], 3 + fsp->alpha));
    else if (r > fsp->radius[n])
	return (fsp->mtot - (fsp->mtot - fsp->mass[n]) *
		  rpow(r / fsp->radius[n], 3 + fsp->beta));
    if (fsp->mr_coef == NULL) {
	fsp->mr_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
	spline(fsp->mr_coef, fsp->radius, fsp->mass, n + 1);
    }
    return (seval(r, fsp->radius, fsp->mass, fsp->mr_coef, n + 1));
}

/*
 * R_MASS_FSP: evaluate radius enclosing given mass.
 */

real r_mass_fsp(fsprof *fsp, real m)
{
    int n = fsp->npoint - 1, k;

    if (m < 0 || m > fsp->mtot)
        error("%s.r_mass_fsp: undefined for m = %g\n", getargv0(), m);
    if (m < fsp->mass[0])
	return (fsp->radius[0] * rpow(m / fsp->mass[0], 1/(3+fsp->alpha)));
    else if (m > fsp->mass[n])
	return (fsp->radius[n] *
		  rpow((fsp->mtot - m) / (fsp->mtot - fsp->mass[n]),
		       1 / (3 + fsp->beta)));
    for (k = 0; fsp->mass[k] == fsp->mass[k+1]; k++);
    while (fsp->mass[n] == fsp->mass[n-1] || fsp->mass[n-1] == fsp->mass[n-2])
        n--;
    if (fsp->rm_coef == NULL) {
	fsp->rm_coef = (real *) allocate(3 * fsp->npoint * sizeof(real));
	if (k > 0 || n < fsp->npoint - 1)
	    eprintf("[%s.r_mass_fsp: using values %d to %d]\n",
		    getargv0(), k, n);
	spline(fsp->rm_coef, fsp->mass + k, fsp->radius + k, n + 1 - k);
    }
    return (seval(m, fsp->mass + k, fsp->radius + k, fsp->rm_coef, n + 1 - k));
}

/*
 * GET_FSPROF: read profile tables from input stream.
 */

fsprof *get_fsprof(stream istr)
{
    fsprof *fsp;

    fsp = (fsprof *) allocate(sizeof(fsprof));
    get_set(istr, "FiniteSphericalProfile");
    get_data(istr, "Npoint", IntType, &fsp->npoint, 0);
    fsp->radius = (real *) allocate(fsp->npoint * sizeof(real));
    fsp->density = (real *) allocate(fsp->npoint * sizeof(real));
    fsp->mass = (real *) allocate(fsp->npoint * sizeof(real));
    get_data(istr, "Radius", RealType, fsp->radius, fsp->npoint, 0);
    get_data(istr, "Density", RealType, fsp->density, fsp->npoint, 0);
    get_data(istr, "Mass", RealType, fsp->mass, fsp->npoint, 0);
    get_data(istr, "Alpha", RealType, &fsp->alpha, 0);
    get_data(istr, "Beta", RealType, &fsp->beta, 0);
    get_data(istr, "Mtot", RealType, &fsp->mtot, 0);
    get_tes(istr, "FiniteSphericalProfile");
    return (fsp);
}

/*
 * PUT_FSPROF: write profile tables to output stream.
 */

void put_fsprof(stream ostr, fsprof *fsp)
{
    put_set(ostr, "FiniteSphericalProfile");
    put_data(ostr, "Npoint", IntType, &fsp->npoint, 0);
    put_data(ostr, "Radius", RealType, fsp->radius, fsp->npoint, 0);
    put_data(ostr, "Density", RealType, fsp->density, fsp->npoint, 0);
    put_data(ostr, "Mass", RealType, fsp->mass, fsp->npoint, 0);
    put_data(ostr, "Alpha", RealType, &fsp->alpha, 0);
    put_data(ostr, "Beta", RealType, &fsp->beta, 0);
    put_data(ostr, "Mtot", RealType, &fsp->mtot, 0);
    put_tes(ostr, "FiniteSphericalProfile");
}

#ifdef TESTBED

/*
 * MAIN: test profile evaluation and input/output routines.
 */

#include "getparam.h"

string defv[] = {		";Test finite spherical profile",
    "in=???",			";Input file with FSP",
    "out=",			";Output file for FSP",
    "npoint=17",		";Number of points to list",
    "r0inv=256.0",		";1/radius of first point",
    "lgrstep=1.0",		";Log2 increment in radius",
    "VERSION=1.4",		";Josh Barnes  20 November 2000",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr, ostr;
    fsprof *fsp;
    int np, i;
    real r0, lgrs, r;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    fsp = get_fsprof(istr);
    if (! strnull(getparam("out"))) {
	ostr = stropen(getparam("out"), "w");
	put_history(ostr);
	put_fsprof(ostr, fsp);
	strclose(ostr);
    }
    np = getiparam("npoint");
    r0 = 1.0 / getdparam("r0inv");
    lgrs = getdparam("lgrstep");
    printf("%12s%12s%12s%12s%12s%12s\n",
	   "radius", "log rho", "drho/dr", "mass", "mtot-mass", "radius(m)");
    for (i = 0; i < np; i++) {
	r = r0 * rpow(2.0, lgrs * i);
	printf("%12.5f%12.7f%12.3e%12.8f%12.8f%12.5f\n", r,
	       rlog10(rho_fsp(fsp, r)), drho_fsp(fsp, r),
	       mass_fsp(fsp, r), fsp->mtot - mass_fsp(fsp, r),
	       r_mass_fsp(fsp, mass_fsp(fsp, r)));
    }
    return (0);
}

#endif
