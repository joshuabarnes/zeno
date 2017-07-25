/*
 * gspeval.c: evaluate GSP at particle positions.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "gsp.h"

string defv[] = {		";Evaluate GSP at particle positions",
  "gsp=???",			";Input GSP file",
  "in=???",			";Input snapshot file",
  "out=???",			";Output snapshot file",
  "option=rho",			";Other choices: grad, mass, phi",
  "VERSION=2.0",		";Josh Barnes  21 June 2017",
  NULL,
};

int main(int argc, string argv[])
{
  string bodyfields[] = { AuxTag, NULL }, intags[MaxBodyFields];
  stream gstr, istr, ostr;
  gsprof *gsp;
  bodyptr btab = NULL, p;
  int nbody;
  real tnow;

  initparam(argv, defv);
  layout_body(bodyfields, Precision, NDIM);
  gstr = stropen(getparam("gsp"), "r");
  get_history(gstr);
  gsp = gsp_read(gstr);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! (get_snap(istr, &btab, &nbody, &tnow, intags, TRUE, NULL) &&
	 set_member(intags, PosTag)))
    error("%s: snapshot input failed or lacked positions\n", getprog());
  if (streq(getparam("option"), "rho"))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = gsp_rho(gsp, absv(Pos(p)));
  else if (streq(getparam("option"), "grad"))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = gsp_grad(gsp, absv(Pos(p)));
  else if (streq(getparam("option"), "mass"))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = gsp_mass(gsp, absv(Pos(p)));
  else if (streq(getparam("option"), "phi"))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = gsp_phi(gsp, absv(Pos(p)));
  else 
    error("%s: unknown option %s\n", getprog(), getparam("option"));
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  put_snap(ostr, &btab, &nbody, &tnow, set_union(bodyfields, intags));
  fflush(NULL);
  return 0;
}
