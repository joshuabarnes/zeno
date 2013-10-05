/*
 * KEPSNAP.C: program to set up initial data for a 2-body orbit.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Compute initial data for 2-body orbit",
    "out=",			";Output snapshot with initial data",
    "mass1=1.25",		";Mass of 1st body",
    "mass2=1.25",		";Mass of 2nd body",
    "r_peri=0.5",		";Separation along X at pericenter.",
				";Use negative value to reverse Jz,",
				";zero value for head-on encounter.",
    "eccent=0.6",		";Eccentricity of orbit (if r_peri != 0).",
				";0 => circle; 1 => parabola.",
    "etot=1.0",			";Total energy of orbit (if r_peri == 0)",
    "t_peri=1.0",		";Time until pericenter",
    "nsteps=2048",		";LF integration steps",
    "thetaz=0.0",		";Angle to rotate results about Z",
    "VERSION=1.0",		";Josh Barnes  31 October 2002",
    NULL,
};

void start_point(vector *, vector *, real *, real, real, real, int);
void linear_orbit(vector *, vector *, real *, real, real);
void newton_acc(vector *, vector *, real *);
real total_energy(vector *, vector *, real *);
void rotate_vectors(vector *, int, real);

string bodytags[] = { PosTag, VelTag, MassTag, NULL };

#define DEG2RAD 0.017453

#define SETVC(v,x,y,z)		/* SET Vector to Components */		\
{									\
    (v)[0] = (x);							\
    (v)[1] = (y);							\
    (v)[2] = (z);							\
}    

#define ADDMULVS(v,u,s)         /* MUL Vect by Scalar, ADD to vect */   \
{                                                                       \
    (v)[0] += (u)[0] * (s);                                             \
    (v)[1] += (u)[1] * (s);                                             \
    (v)[2] += (u)[2] * (s);                                             \
}

int main(int argc, string argv[])
{
    real mass[2], eccent, r_peri, t_peri, tzero = 0.0;
    vector pos[2], vel[2], dr, dv;
    bodyptr btab;
    stream outstr;
    int nbody = 2;

    initparam(argv, defv);
    mass[0] = getdparam("mass1");
    mass[1] = getdparam("mass2");
    if (mass[0] < 0 || mass[1] < 0)
	error("%s: negative mass\n", getargv0());
    eccent = getdparam("eccent");
    if (eccent < 0)
	error("%s: negative eccentricity\n", getargv0());
    r_peri = getdparam("r_peri");
    t_peri = getdparam("t_peri");
    if (ABS(r_peri) > 0)
        start_point(pos, vel, mass, eccent, r_peri, t_peri,
		    getiparam("nsteps"));
    else
        linear_orbit(pos, vel, mass, getdparam("etot"), t_peri);
    rotate_vectors(pos, 2, getdparam("thetaz"));
    rotate_vectors(vel, 2, getdparam("thetaz"));
    SUBV(dr, pos[1], pos[0]);
    SUBV(dv, vel[1], vel[0]);
    eprintf("[%s: delta_pos: %f,%f,%f (%f)]\n",
	    getargv0(), dr[0], dr[1], dr[2], absv(dr));
    eprintf("[%s: delta_vel: %f,%f,%f (%f)]\n",
	    getargv0(), dv[0], dv[1], dv[2], absv(dv));
    if (! strnull(getparam("out"))) {
	layout_body(bodytags, Precision, NDIM);
	btab = (bodyptr) allocate(2 * SizeofBody);
	Mass(NthBody(btab, 0)) = mass[0];
	Mass(NthBody(btab, 1)) = mass[1];
	SETV(Pos(NthBody(btab, 0)), pos[0]);
	SETV(Pos(NthBody(btab, 1)), pos[1]);
	SETV(Vel(NthBody(btab, 0)), vel[0]);
	SETV(Vel(NthBody(btab, 1)), vel[1]);
	outstr = stropen(getparam("out"), "w");
	put_history(outstr);
	put_snap(outstr, &btab, &nbody, &tzero, bodytags);
	strclose(outstr);
    }
    return (0);
}

void start_point(vector *r, vector *v, real *m,
		 real eccent, real r_peri, real t_peri, int nsteps)
{
    real mtot, v_peri, dt;
    int n;
    vector a[2];

    mtot = m[0] + m[1];
    v_peri = rsqrt((1 + eccent) * mtot / ABS(r_peri));
    if (eccent < 1.0)
        eprintf("[%s: v_peri = %f  r_apoc = %f]\n", getargv0(),
		v_peri, ABS(r_peri) * (1 + eccent) / (1 - eccent));
    else
        eprintf("[%s: v_peri = %f]\n", getargv0(), v_peri);
    SETVC(r[0],  r_peri * m[1]/mtot, 0.0, 0.0);	/* set position at pericent */
    SETVC(r[1], -r_peri * m[0]/mtot, 0.0, 0.0);
    SETVC(v[0], 0.0, -v_peri * m[1]/mtot, 0.0);	/* set velocity at pericent */
    SETVC(v[1], 0.0,  v_peri * m[0]/mtot, 0.0);
    eprintf("[%s: total energy: %f", getargv0(), total_energy(r, v, m));
    dt = - t_peri / nsteps;
    for (n = 0; n < nsteps; n++) {		/* step back from pericent  */
	newton_acc(a, r, m);
	ADDMULVS(v[0], a[0], 0.5 * dt);
	ADDMULVS(v[1], a[1], 0.5 * dt);
	ADDMULVS(r[0], v[0], dt);
	ADDMULVS(r[1], v[1], dt);
	newton_acc(a, r, m);
	ADDMULVS(v[0], a[0], 0.5 * dt);
	ADDMULVS(v[1], a[1], 0.5 * dt);
    }
    eprintf(" -> %f]\n", total_energy(r, v, m));
}

void linear_orbit(vector *r, vector *v, real *m, real etot, real t_peri)
{
    real A, B, eta, t_eta, r12, v12;
    int ncyc = 0;

    if (etot > 0) {
	A = m[0] * m[1] / (2 * etot);
	B = rsqrt(rqbe(m[0] * m[1]) / (8 * (m[0] + m[1]) * rqbe(etot)));
	eta = MAX(asinh((double) (t_peri / B)), 0.1);
	t_eta = B * (rsinh(eta) - eta);
	while (ABS(t_eta - t_peri) > 1.0e-4 * t_peri) {
	    eta = eta - (t_eta - t_peri) / (B * (rcosh(eta) - 1));
	    t_eta = B * (rsinh(eta) - eta);
	    if (++ncyc > 50)
	        error("%s: convergence failed\n", getargv0());
	}
	r12 = A * (rcosh(eta) - 1);
	v12 = (A / B) * rsinh(eta) / (rcosh(eta) - 1);
    } else
        error("%s: etot <= 0 not implemented for r_peri == 0\n", getargv0());
    SETVC(r[0],  r12 * m[1]/ (m[0] + m[1]), 0.0, 0.0);
    SETVC(r[1], -r12 * m[0]/ (m[0] + m[1]), 0.0, 0.0);
    SETVC(v[0], -v12 * m[1]/ (m[0] + m[1]), 0.0, 0.0);
    SETVC(v[1],  v12 * m[0]/ (m[0] + m[1]), 0.0, 0.0);
    eprintf("[%s: eta = %f  t_peri = %f  etot = %f]\n",
	    getargv0(), eta, t_eta, total_energy(r, v, m));
}

void newton_acc(vector *a, vector *r, real *m)
{
    vector dr;
    real dr2, dr3;

    SUBV(dr, r[0], r[1]);
    DOTVP(dr2, dr, dr);
    dr3 = dr2 * rsqrt(dr2);
    MULVS(a[0], dr, -m[1] / dr3);
    MULVS(a[1], dr,  m[0] / dr3);
}

real total_energy(vector *r, vector *v, real *m)
{
    real dr, v2, e_tot;

    DISTV(dr, r[0], r[1]);
    e_tot = - m[0] * m[1] / dr;
    DOTVP(v2, v[0], v[0]);
    e_tot += 0.5 * m[0] * v2;
    DOTVP(v2, v[1], v[1]);
    e_tot += 0.5 * m[1] * v2;
    return (e_tot);
}    

void rotate_vectors(vector *vec, int nvec, real theta)
{
    real s, c;
    int i;
    vector vtmp;

    s = rsin(DEG2RAD * theta);
    c = rcos(DEG2RAD * theta);
    for (i = 0; i < nvec; i++) {
	SETV(vtmp, vec[i]);
	vec[i][0] = c * vtmp[0] + s * vtmp[1];
	vec[i][1] = c * vtmp[1] - s * vtmp[0];
    }
}
