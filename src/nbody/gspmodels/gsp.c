/*
 * GSP.C: routines for operating on general spherical profiles.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "filestruct.h"
#include "gsp.h"

#if (!defined(LINUX) && !defined(MACOSX))
#include <ieeefp.h>
#endif

/*
 * RHO_GSP: evaluate density at given radius.
 */

real rho_gsp(gsprof *gsp, real r)
{
    int n = gsp->npoint - 1;

    if (r < 0)
        error("%s.rho_gsp: undefined for r = %g\n", getargv0(), r);
    if (r < gsp->radius[0])
	return (gsp->density[0] * rpow(r / gsp->radius[0], gsp->alpha));
    else if (r > gsp->radius[n])
	return (gsp->density[n] * rpow(r / gsp->radius[n], gsp->beta));
    if (gsp->dr_coef == NULL) {
	gsp->dr_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
	spline(gsp->dr_coef, gsp->radius, gsp->density, n + 1);
    }
    return (seval(r, gsp->radius, gsp->density, gsp->dr_coef, n + 1));
}

/*
 * DRHO_GSP: evaluate radial derivative of density at given radius.
 */

real drho_gsp(gsprof *gsp, real r)
{
    int n = gsp->npoint - 1;

    if (r < 0)
        error("%s.drho_gsp: undefined for r = %g\n", getargv0(), r);
    if (r < gsp->radius[0])
	return (gsp->density[0] * rpow(r / gsp->radius[0], gsp->alpha) *
		  gsp->alpha / r);
    else if (r > gsp->radius[n])		/* less glitches at large r */
	return (gsp->density[n] * rpow(r / gsp->radius[n], gsp->beta) *
		  gsp->beta / r);
    if (gsp->dr_coef == NULL) {
	gsp->dr_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
	spline(gsp->dr_coef, gsp->radius, gsp->density, n + 1);
    }
    return (spldif(r, gsp->radius, gsp->density, gsp->dr_coef, n + 1));
}

/*
 * MASS_GSP: evaluate enclosed mass at given radius.
 */

real mass_gsp(gsprof *gsp, real r)
{
    int n = gsp->npoint - 1;

    if (r < 0)
        error("%s.mass_gsp: undefined for r = %g\n", getargv0(), r);
    if (r < gsp->radius[0])
	return (gsp->mass[0] * rpow(r / gsp->radius[0], 3 + gsp->alpha));
    else if (r > gsp->radius[n])
	return (gsp->mtot - (gsp->mtot - gsp->mass[n]) *
		  rpow(r / gsp->radius[n], 3 + gsp->beta));
    if (gsp->mr_coef == NULL) {
	gsp->mr_coef = (real *) allocate(3 * (n + 1) * sizeof(real));
	spline(gsp->mr_coef, gsp->radius, gsp->mass, n + 1);
    }
    return (seval(r, gsp->radius, gsp->mass, gsp->mr_coef, n + 1));
}

/*
 * R_MASS_GSP: evaluate radius enclosing given mass.
 */

real r_mass_gsp(gsprof *gsp, real m)
{
    int n = gsp->npoint - 1, k = 0;

    if (m < 0 || m > gsp->mtot)
        error("%s.r_mass_gsp: undefined for m = %g\n", getargv0(), m);
    if (m < gsp->mass[0])
	return (gsp->radius[0] * rpow(m / gsp->mass[0], 1/(3+gsp->alpha)));
    else if (m > gsp->mass[n])
	return (gsp->radius[n] *
		  rpow((gsp->mtot - m) / (gsp->mtot - gsp->mass[n]),
		       1 / (3 + gsp->beta)));
    while (gsp->mass[k] == gsp->mass[k+1])
        k++;
    while (gsp->mass[n] == gsp->mass[n-1] ||
	   gsp->mass[n-1] == gsp->mass[n-2])
        n--;
    if (gsp->rm_coef == NULL) {
	gsp->rm_coef = (real *) allocate(3 * gsp->npoint * sizeof(real));
	if (k > 0 || n < gsp->npoint - 1)
	    eprintf("[%s.r_mass_gsp: spline range %d to %d]\n",
		    getargv0(), k, n);
	spline(gsp->rm_coef, gsp->mass + k, gsp->radius + k, n + 1 - k);
	if (isnan((double) gsp->rm_coef[0]))
	    error("%s.r_mass_gsp: spline fit undefined\n", getargv0());
    }
    return (seval(m, gsp->mass + k, gsp->radius + k, gsp->rm_coef, n + 1 - k));
}

/*
 * GET_GSPROF: read profile tables from input stream.
 */

gsprof *get_gsprof(stream istr)
{
    string gsptag = "GeneralSphericalProfile";
    gsprof *gsp;

    if (get_tag_ok(istr, "FiniteSphericalProfile"))	/* old-style file?  */
	gsptag = "FiniteSphericalProfile";		/* use old tag      */
    gsp = (gsprof *) allocate(sizeof(gsprof));
    get_set(istr, gsptag);
    get_data(istr, "Npoint", IntType, &gsp->npoint, 0);
    gsp->radius = (real *) allocate(gsp->npoint * sizeof(real));
    gsp->density = (real *) allocate(gsp->npoint * sizeof(real));
    gsp->mass = (real *) allocate(gsp->npoint * sizeof(real));
    get_data(istr, "Radius", RealType, gsp->radius, gsp->npoint, 0);
    get_data(istr, "Density", RealType, gsp->density, gsp->npoint, 0);
    get_data(istr, "Mass", RealType, gsp->mass, gsp->npoint, 0);
    get_data(istr, "Alpha", RealType, &gsp->alpha, 0);
    get_data(istr, "Beta", RealType, &gsp->beta, 0);
    get_data(istr, "Mtot", RealType, &gsp->mtot, 0);
    get_tes(istr, gsptag);
    return (gsp);
}

/*
 * PUT_GSPROF: write profile tables to output stream.
 */

void put_gsprof(stream ostr, gsprof *gsp)
{
    put_set(ostr, "GeneralSphericalProfile");
    put_data(ostr, "Npoint", IntType, &gsp->npoint, 0);
    put_data(ostr, "Radius", RealType, gsp->radius, gsp->npoint, 0);
    put_data(ostr, "Density", RealType, gsp->density, gsp->npoint, 0);
    put_data(ostr, "Mass", RealType, gsp->mass, gsp->npoint, 0);
    put_data(ostr, "Alpha", RealType, &gsp->alpha, 0);
    put_data(ostr, "Beta", RealType, &gsp->beta, 0);
    put_data(ostr, "Mtot", RealType, &gsp->mtot, 0);
    put_tes(ostr, "GeneralSphericalProfile");
}

/*
 * FREE_GSPROF: deallocate a gsp and associated tables.
 */

void free_gsprof(gsprof *gsp)
{
  if (gsp->radius != NULL)
    free((void *) gsp->radius);
  if (gsp->density != NULL)
    free((void *) gsp->density);
  if (gsp->mass != NULL)
    free((void *) gsp->mass);
  if (gsp->phi != NULL)
    free((void *) gsp->phi);
  if (gsp->dr_coef != NULL)
    free((void *) gsp->dr_coef);
  if (gsp->mr_coef != NULL)
    free((void *) gsp->mr_coef);
  if (gsp->rm_coef != NULL)
    free((void *) gsp->rm_coef);
  if (gsp->pr_coef != NULL)
    free((void *) gsp->pr_coef);
  if (gsp->rp_coef != NULL)
    free((void *) gsp->rp_coef);
  free((void *) gsp);
}

#ifdef TESTBED

/*
 * MAIN: test profile evaluation and input/output routines.
 */

#include "getparam.h"

string defv[] = {		";Test general spherical profile",
    "in=???",			";Input file with GSP",
    "out=",			";Output file for GSP",
    "npoint=17",		";Number of points to list",
    "r0inv=256.0",		";1/radius of first point",
    "lgrstep=1.0",		";Log2 increment in radius",
    "VERSION=1.4",		";Josh Barnes  11 July 2010",
    NULL,
};

int main(int argc, string argv[])
{
    stream istr, ostr;
    gsprof *gsp;
    int np, i;
    real r0, lgrs, r;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    gsp = get_gsprof(istr);
    if (! strnull(getparam("out"))) {
	ostr = stropen(getparam("out"), "w");
	put_history(ostr);
	put_gsprof(ostr, gsp);
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
	       rlog10(rho_gsp(gsp, r)), drho_gsp(gsp, r),
	       mass_gsp(gsp, r), gsp->mtot - mass_gsp(gsp, r),
	       r_mass_gsp(gsp, mass_gsp(gsp, r)));
    }
    free_gsprof(gsp);
    return (0);
}

#endif
