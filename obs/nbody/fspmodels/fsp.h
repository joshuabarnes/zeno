/*
 * FSP.H: structure for profile of finite-mass, spherical model.
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
} fsprof;

/* Low-level FSP functions; see fsp.c for code. */

real rho_fsp(fsprof *, real);

real drho_fsp(fsprof *, real);

real mass_fsp(fsprof *, real);

real r_mass_fsp(fsprof *, real);

fsprof *get_fsprof(stream);

void put_fsprof(stream, fsprof *);

/* High-level FSP functions. */

real phi_fsp(fsprof *, real);

real r_phi_fsp(fsprof *, real);

real *calc_sig2_fsp(fsprof *, fsprof *, real);

real sig2_fsp(fsprof *, fsprof *, real, real *, real);

/* Constructor functions for FSP models. */

fsprof *expdfsp(real, real, int, real, real);

fsprof *gammafsp(real, real, real, int, real, real);

fsprof *halofsp_e(real, real, real, int, real, real);

fsprof *halofsp_g(real, real, real, int, real, real);

fsprof *halofsp_sw(real, real, real, int, real, real);

fsprof *isothfsp(real, real, real, int, real, real);

fsprof *plumfsp(real, real, int, real, real);

fsprof *polyfsp(real, real, real, int);

