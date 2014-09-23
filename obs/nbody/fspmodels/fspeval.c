/*
 * FSPEVAL.C: evaluate FSP at particle positions.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "phatbody.h"
#include "fsp.h"

string defv[] = {		";Evaluate FSP at particle positions",
    "fsp=???",			";Input FSP file",
    "in=???",			";Input snapshot file",
    "out=???",			";Output snapshot file",
    "option=rho",		";Other choices: drho, mass, phi",
    "VERSION=1.0",		";Josh Barnes  6 June 2007",
    NULL,
};

string bodyfields[] = { AuxTag, NULL };

int main(int argc, string argv[])
{
  stream fstr, istr, ostr;
  fsprof *fsp;
  bodyptr btab = NULL, p;
  int nbody;
  real tnow, r;
  string intags[MaxBodyFields];

  initparam(argv, defv);
  layout_body(bodyfields, Precision, NDIM);
  fstr = stropen(getparam("fsp"), "r");
  get_history(fstr);
  fsp = get_fsprof(fstr);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab, &nbody, &tnow, intags, TRUE))
    error("%s: snapshot input failed\n", getargv0());
  if (! set_member(intags, PosTag))
    error("%s: position data missing\n", getargv0());
  if (streq(getparam("option"), "rho"))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = rho_fsp(fsp, absv(Pos(p)));
  else if (streq(getparam("option"), "drho"))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = drho_fsp(fsp, absv(Pos(p)));
  else if (streq(getparam("option"), "mass"))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = mass_fsp(fsp, absv(Pos(p)));
  else if (streq(getparam("option"), "phi"))
    for (p = btab; p < NthBody(btab, nbody); p = NextBody(p))
      Aux(p) = phi_fsp(fsp, absv(Pos(p)));
  else 
    error("%s: unknown option %s\n", getargv0(), getparam("option"));
  if (! strnull(getparam("out"))) {
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_snap(ostr, &btab, &nbody, &tnow, set_union(bodyfields, intags));
    strclose(ostr);
  }
  return (0);
}
