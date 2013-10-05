/*
 * GSPSMOOTHFIT.C: Fit smoothed gsp to gravitational field.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "gsp.h"

string defv[] = {		";Fit smoothed gsp to gravitational field",
    "gsp=???",			";GSP file with unsmoothed mass profile",
    "in=???",			";Snapshot file with pos, phi, acc data",
    "eps=0.005:0.015",		";Plummer smoothing length range",
    "neps=21",			";Number of eps values to try",
    "kappa=1.0:3.0",		";Interpolation parameter range",
    "ncycle=20",		";Number of binary search cycles",
    "zerophi=true",		";If false, find zero of acceleration",
    "VERSION=1.0",		";Josh Barnes  14 July 2010",
    NULL,
};

string bodyfields[] = { PosTag, PhiTag, AccTag, AuxTag, NULL };

void eval_fit(real res[], gsprof *gsp0, real eps, real kappa,
	      bodyptr btab, int nbody);

int main(int argc, string argv[])
{
  stream gstr, istr;
  gsprof *gsp0;
  bodyptr btab = NULL, p;
  int nbody, kzero, i, j;
  real tnow, eps[2], kappa[2], eps1, res0[4], res1[4], res2[4], kappa1;
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
	 "eps", "kappa", "dphiavg", "dphirms", "daccavg", "daccrms");
  for (i = 0; i < getiparam("neps"); i++) {
    eps1 = eps[0] + (eps[1] - eps[0]) * ((real) i) / (getiparam("neps") - 1.0);
    setrange(kappa, getparam("kappa"));
    eval_fit(res0, gsp0, eps1, kappa[0], btab, nbody);
    eval_fit(res2, gsp0, eps1, kappa[1], btab, nbody);
    if (res0[kzero] < 0.0 && res2[kzero] > 0.0) {
      for (j = 0; j < getiparam("ncycle"); j++) {
	kappa1 = (kappa[0] + kappa[1]) / 2.0;
	eval_fit(res1, gsp0, eps1, kappa1, btab, nbody);
	kappa[res1[kzero] < 0.0 ? 0 : 1] = kappa1;
      }
      printf(" %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n",
	     eps1, kappa1, res1[0], res1[1], res1[2], res1[3]);
    } else
      eprintf("[%s: eps1 = %f  res0[%d] = %f  res2[%d] = %f]\n",
	      getargv0(), eps1, kzero, res0[kzero], kzero, res2[kzero]);
  }
  return (0);
}

void eval_fit(real res[], gsprof *gsp0, real eps, real kappa,
	      bodyptr btab, int nbody)
{
  gsprof *gsp1;
  bodyptr p;
  double wsum, dphisum, dphi2sum, daccsum, dacc2sum;
  real r, phi1, acc1;

  gsp1 = gspsmooth(gsp0, eps, kappa, NULL);
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
