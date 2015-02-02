/*
 * SNAPSTREAM.C: Read one snapshot, produce a stream of them.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Emit stream of snapshots",
  "in=???",			";Single snapshot input",
  "out=???",			";Snapshot stream output",
  "tstart=0.0",			";Time value for 1st snap",
  "freq=64.0",			";Snapshots per unit time",
  "nstream=256",		";Number of snapshots output",
  "VERSION=1.1",		";Josh Barnes  2 February 2015",
  NULL,
};

string bodytags[] = {
    PosTag, VelTag, MassTag, PhiTag, AccTag, AuxTag, KeyTag, NULL,
};

int main(int argc, string argv[])
{
  stream istr, ostr;
  bodyptr btab = NULL;
  int nbody, n;
  real tnow;
  string iotags[MaxBodyFields];

  initparam(argv, defv);
  layout_body(bodytags, Precision, NDIM);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  if (! get_snap(istr, &btab, &nbody, &tnow, iotags, FALSE))
    error("%s: no data in input file\n", getargv0());
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  for (n = 0; n < getiparam("nstream"); n++) {
    tnow = getdparam("tstart") + n / getdparam("freq");
    put_snap(ostr, &btab, &nbody, &tnow, iotags);
  }
  return (0);
}
