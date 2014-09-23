/*
 * GSP.H: structure for profile of general spherical model.
 */

typedef struct {
    int npoint;			/* number of points in radial arrays        */
    real *radius;		/* radii defining model profile		    */
    real *density;		/* density tabulated at each radius	    */
    real *mass;			/* enclosed mass within each radius	    */
    real *phi;			/* grav. potential at each radius           */
    real alpha;			/* density power-law index for small r	    */
    real beta;			/* density power-law index for large r	    */
    real mtot;			/* total mass of model (within infinity)    */
    real *dr_coef;		/* coefs for spline-fit to density(radius)  */
    real *mr_coef;		/* coefs for spline-fit to mass(radius)	    */
    real *rm_coef;		/* coefs for spline-fit to radius(mass)	    */
    real *pr_coef;		/* coefs for spline-fit to phi(radius)      */
    real *rp_coef;		/* coefs for spline-fit to radius(phii)      */
} gsprof;

/* Low-level GSP functions; see gsp.c for code. */

real rho_gsp(gsprof *, real);

real drho_gsp(gsprof *, real);

real mass_gsp(gsprof *, real);

real r_mass_gsp(gsprof *, real);

gsprof *get_gsprof(stream);

void put_gsprof(stream, gsprof *);

void free_gsprof(gsprof *gsp);

/* High-level GSP functions. */

void calc_phi_gsp(gsprof *, string);

real phi_gsp(gsprof *, real);

real r_phi_gsp(gsprof *, real);

real *calc_sig2_gsp(gsprof *, gsprof *, real);

real sig2_gsp(gsprof *, gsprof *, real, real *, real);

/* Constructor functions for GSP models. */

gsprof *expdgsp(real, real, real, int, real, real);

gsprof *gammagsp(real, real, real, int, real, real);

gsprof *halogsp_e(real, real, real, int, real, real);

gsprof *halogsp_g(real, real, real, int, real, real);

gsprof *halogsp_sw(real, real, real, int, real, real);

gsprof *isothgsp(real, real, real, int, real, real);

gsprof *plumgsp(real, real, int, real, real);

gsprof *polygsp(real, real, real, int);

/* Transformation functions for GSP models. */

gsprof *gspsmooth(gsprof *gsp, real eps, real kappa, string trace);
