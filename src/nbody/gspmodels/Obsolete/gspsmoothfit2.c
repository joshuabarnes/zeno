/*
 * GSPSMOOTHFIT2.C: Fit smoothed gsp to gravitational field.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "assert.h"
#include "gsp.h"

string defv[] = {		";Fit smoothed gsp to gravitational field",
    "gsp=???",			";GSP file with unsmoothed mass profile",
    "in=???",			";Snapshot file with pos, phi, acc data",
    "eps=0.005:0.015",		";Plummer smoothing length range",
    "neps=21",			";Number of eps values to try",
    "eta=1.0:3.0",		";Interpolation parameter range",
    "ncycle=20",		";Number of binary search cycles",
    "zerophi=true",		";If false, find zero of acceleration",
    "VERSION=1.0",		";Josh Barnes  14 July 2010",
    NULL,
};

string bodyfields[] = { PosTag, PhiTag, AccTag, AuxTag, NULL };

void eval_fit(real res[], gsprof *gsp0, real eps, real eta,
	      bodyptr btab, int nbody);

gsprof *gspsmooth2(gsprof *gsp, real rho_0, real eta, string trace);

int main(int argc, string argv[])
{
  stream gstr, istr;
  gsprof *gsp0;
  bodyptr btab = NULL, p;
  int nbody, kzero, i, j;
  real tnow, eps[2], eta[2], eps1, res0[4], res1[4], res2[4], eta1;
  string intags[MaxBodyFields];

  initparam(argv, defv);
  layout_body(bodyfields, Precision, NDIM);
  gstr = stropen(getparam("gsp"), "r");
  get_history(gstr);
  gsp0 = get_gsprof(gstr);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab, &nbody, &tnow, intags, TRUE))
    error("%s: snapshot input failed\n", getargv0());
  if (! (set_member(intags, PosTag) &&
	   set_member(intags, PhiTag) && set_member(intags, AccTag)))
    error("%s: necessary data missing\n", getargv0());
  if (! set_member(intags, AuxTag))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = 1.0;
  setrange(eps, getparam("eps"));
  kzero = (getbparam("zerophi") ? 0 : 2);
  printf("#%11s %11s %11s %11s %11s %11s\n",
	 "eps", "eta", "dphiavg", "dphirms", "daccavg", "daccrms");
  for (i = 0; i < getiparam("neps"); i++) {
    eps1 = eps[0] + (eps[1] - eps[0]) * ((real) i) / (getiparam("neps") - 1.0);
    setrange(eta, getparam("eta"));
    eval_fit(res0, gsp0, eps1, eta[0], btab, nbody);
    eval_fit(res2, gsp0, eps1, eta[1], btab, nbody);
    if (res0[kzero] < 0.0 && res2[kzero] > 0.0) {
      for (j = 0; j < getiparam("ncycle"); j++) {
	eta1 = (eta[0] + eta[1]) / 2.0;
	eval_fit(res1, gsp0, eps1, eta1, btab, nbody);
	eta[res1[kzero] < 0.0 ? 0 : 1] = eta1;
      }
      printf(" %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n",
	     eps1, eta1, res1[0], res1[1], res1[2], res1[3]);
    } else
      eprintf("[%s: eps1 = %f  res0[%d] = %f  res2[%d] = %f]\n",
	      getargv0(), eps1, kzero, res0[kzero], kzero, res2[kzero]);
  }
  return (0);
}

void eval_fit(real res[], gsprof *gsp0, real eps, real eta,
	      bodyptr btab, int nbody)
{
  gsprof *gsp1;
  bodyptr p;
  double wsum, dphisum, dphi2sum, daccsum, dacc2sum;
  real rho_0, r, phi1, acc1;

  rho_0 = - gsp0->alpha * rho_gsp(gsp0, eps);
					/* use semi-bogus expr for rho(0) */
  gsp1 = gspsmooth2(gsp0, rho_0, eta, NULL);
  calc_phi_gsp(gsp1, NULL);
  wsum = dphisum = dphi2sum = daccsum = dacc2sum = 0.0;
  for (p = btab; p < NthBody(btab, nbody); p = NextBody(p)) {
    r = absv(Pos(p));
    phi1 = phi_gsp(gsp1, r);
    acc1 = - mass_gsp(gsp1, r) / (r * r);
    wsum += Aux(p);
    dphisum += Aux(p) * (Phi(p) - phi1);
    dphi2sum += Aux(p) * rsqr(Phi(p) - phi1);
    daccsum += Aux(p) * (dotvp(Acc(p), Pos(p)) / r - acc1);
    dacc2sum += Aux(p) * rsqr(dotvp(Acc(p), Pos(p)) / r - acc1);
  }
  res[0] = dphisum / wsum;
  res[1] = rsqrt(dphi2sum / wsum);
  res[2] = daccsum / wsum;
  res[3] = rsqrt(dacc2sum / wsum);
  free_gsprof(gsp1);
}

/*
 * GSPSMOOTH2: new version based on mass-interpolation formula,
 * from Barnes (2011).
 */

gsprof *gspsmooth2(gsprof *gsp, real rho_0, real eta, string trace)
{
  gsprof *sgsp;
  int i;
  real r_i, mass_0, mass_r, rho_r;

  assert(rho_0 > 0.0 && eta > 0.0 && gsp->alpha <= 0.0);
  if (trace != NULL)
    eprintf("[%s:  rho_0 = %f  eta = %f]\n", trace, rho_0, eta);
  sgsp = (gsprof *) allocate(sizeof(gsprof));
  sgsp->npoint = gsp->npoint;
  sgsp->radius = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->density = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->mass = (real *) allocate(sgsp->npoint * sizeof(real));
  sgsp->alpha = 0.0;
  sgsp->beta = gsp->beta;
  sgsp->mtot = gsp->mtot;
  for (i = 0; i < gsp->npoint; i++) {
    r_i = gsp->radius[i];
    mass_0 = (4 * PI / 3.0) * rqbe(r_i) * rho_0;
    mass_r = rpow(rpow(gsp->mass[i], -eta) + rpow(mass_0, -eta), -1/eta);
    rho_r = rpow(mass_r / gsp->mass[i], eta+1) * gsp->density[i] +
            rpow(mass_r / mass_0, eta+1) * rho_0;
    sgsp->radius[i] = r_i;
    sgsp->density[i] = rho_r;
    sgsp->mass[i] = mass_r;
  }
  return (sgsp);
}
