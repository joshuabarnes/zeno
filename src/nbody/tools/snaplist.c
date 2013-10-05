/*
 * SNAPLIST.C: read a snapshot file and print it as text.
 */

#include "stdinc.h"
#include "getparam.h"
#include "filestruct.h"
#include "vectdefs.h"
#include "phatbody.h"
#include "strset.h"

#define XX ", "

string defv[] = {	       	";List contents of snapshot file",
  "in=???",			";Input file of N-body snapshots",
  "fields=" PosTag,		";Data fields to include in output.",
				";Options: " TimeTag XX MassTag XX PosTag XX
				    VelTag ",",
				";" AccTag XX PhiTag XX SmoothTag ",",
				";" RhoTag XX EntFuncTag XX UinternTag ",",
				";" UdotIntTag XX UdotRadTag XX UdotVisTag ",",
				";" TauTag XX BirthTag XX DeathTag ",",
				";" TypeTag XX KeyTag XX AuxTag XX AuxVecTag ".",
  "times=all",			";Range of times to list",
  "hfmt=%1s%11s",		";Format for column headers",
  "ifmt= %11d",			";Format for integer values",
  "rfmt= %11.5g",		";Format for real values",
  "keyhead=key",		";Name used to head key column",
  "auxhead=aux,auxx,auxy,auxz",	";Names used to head aux columns",
  "VERSION=1.4",		";Josh Barnes  24 December 2012",
  NULL,
};

void print_header(string *fields, string hfmt,
		  string keyhead, string *auxhead);
void print_data(bodyptr btab, int nbody, real tnow, string *fields,
		string ifmt, string rfmt);

int main(int argc, string argv[])
{
  string *datafields, *bodyfields, times, ifmt, rfmt;
  stream istr;
  bodyptr btab = NULL;
  int nbody;
  real tnow;
  string intags[MaxBodyFields];

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  datafields = burststring(getparam("fields"), ", ");
  if (set_member(datafields, TimeTag))
    bodyfields = set_diff(datafields, set_cons(TimeTag, NULL));
  else
    bodyfields = datafields;
  layout_body(bodyfields, Precision, NDIM);
  times = getparam("times");
  ifmt = getparam("ifmt");
  rfmt = getparam("rfmt");
  print_header(datafields, getparam("hfmt"), getparam("keyhead"),
	       burststring(getparam("auxhead"), ", "));
  while (get_snap_t(istr, &btab, &nbody, &tnow, intags, FALSE, times)) {
    if (! set_subset(intags, bodyfields))
      error("%s: one or more required fields not found\n", getargv0());
    print_data(btab, nbody, tnow, datafields, ifmt, rfmt);
    skip_history(istr);
  }
  return (0);
}

void print_header(string *fields, string hfmt,
		  string keyhead, string *auxhead)
{
    string sep = "#";

    if (auxhead[0]==NULL || auxhead[1]==NULL ||
	auxhead[2]==NULL || auxhead[3]==NULL)
        error("%s: auxhead must have four comma-separated values\n",
	      getargv0());
    if (set_member(fields, TimeTag)) {
	printf(hfmt, sep, "time");
	sep = " ";
    }
    if (set_member(fields, MassTag)) {
	printf(hfmt, sep, "mass");
	sep = " ";
    }
    if (set_member(fields, PosTag)) {
        printf(hfmt, sep, "x");
	printf(hfmt, " ", "y");
	printf(hfmt, " ", "z");
	sep = " ";
    }
    if (set_member(fields, VelTag)) {
        printf(hfmt, sep, "vx");
	printf(hfmt, " ", "vy");
	printf(hfmt, " ", "vz");
	sep = " ";
    }
    if (set_member(fields, AccTag)) {
        printf(hfmt, sep, "ax");
	printf(hfmt, " ", "ay");
	printf(hfmt, " ", "az");
	sep = " ";
    }
    if (set_member(fields, PhiTag)) {
	printf(hfmt, sep, "phi");
	sep = " ";
    }
    if (set_member(fields, SmoothTag)) {
	printf(hfmt, sep, "smooth");
	sep = " ";
    }
    if (set_member(fields, RhoTag)) {
	printf(hfmt, sep, "rho");
	sep = " ";
    }
    if (set_member(fields, EntFuncTag)) {
	printf(hfmt, sep, "entf");
	sep = " ";
    }
    if (set_member(fields, UinternTag)) {
	printf(hfmt, sep, "uint");
	sep = " ";
    }
    if (set_member(fields, UdotIntTag)) {
	printf(hfmt, sep, "udot");
	sep = " ";
    }
    if (set_member(fields, UdotRadTag)) {
	printf(hfmt, sep, "udotrad");
	sep = " ";
    }
    if (set_member(fields, UdotVisTag)) {
	printf(hfmt, sep, "udotvis");
	sep = " ";
    }
    if (set_member(fields, TauTag)) {
	printf(hfmt, sep, "tau");
	sep = " ";
    }
    if (set_member(fields, BirthTag)) {
	printf(hfmt, sep, "birth");
	sep = " ";
    }
    if (set_member(fields, DeathTag)) {
	printf(hfmt, sep, "death");
	sep = " ";
    }
    if (set_member(fields, TypeTag)) {
	printf(hfmt, sep, "type");
	sep = " ";
    }
    if (set_member(fields, KeyTag)) {
	printf(hfmt, sep, keyhead);
	sep = " ";
    }
    if (set_member(fields, AuxTag)) {
	printf(hfmt, sep, auxhead[0]);
	sep = " ";
    }
    if (set_member(fields, AuxVecTag)) {
        printf(hfmt, sep, auxhead[1]);
	printf(hfmt, " ", auxhead[2]);
	printf(hfmt, " ", auxhead[3]);
	sep = " ";
    }
    printf("\n");
}

void print_data(bodyptr btab, int nbody, real tnow, string *fields,
		string ifmt, string rfmt)
{
    bodyptr bp;

    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	if (set_member(fields, TimeTag))
	    printf(rfmt, tnow);
	if (set_member(fields, MassTag))
	    printf(rfmt, Mass(bp));
	if (set_member(fields, PosTag)) {
	    printf(rfmt, Pos(bp)[0]);
	    printf(rfmt, Pos(bp)[1]);
	    printf(rfmt, Pos(bp)[2]);
	}
	if (set_member(fields, VelTag)) {
	    printf(rfmt, Vel(bp)[0]);
	    printf(rfmt, Vel(bp)[1]);
	    printf(rfmt, Vel(bp)[2]);
	}
	if (set_member(fields, AccTag)) {
	    printf(rfmt, Acc(bp)[0]);
	    printf(rfmt, Acc(bp)[1]);
	    printf(rfmt, Acc(bp)[2]);
	}
	if (set_member(fields, PhiTag))
	    printf(rfmt, Phi(bp));
	if (set_member(fields, SmoothTag))
	    printf(rfmt, Smooth(bp));
	if (set_member(fields, RhoTag))
	    printf(rfmt, Rho(bp));
	if (set_member(fields, EntFuncTag))
	    printf(rfmt, EntFunc(bp));
	if (set_member(fields, UinternTag))
	    printf(rfmt, Uintern(bp));
	if (set_member(fields, UdotIntTag))
	    printf(rfmt, UdotInt(bp));
        if (set_member(fields, UdotRadTag))
	    printf(rfmt, UdotRad(bp));
        if (set_member(fields, UdotVisTag))
	    printf(rfmt, UdotVis(bp));
	if (set_member(fields, TauTag))
	    printf(rfmt, Tau(bp));
	if (set_member(fields, BirthTag))
	    printf(rfmt, Birth(bp));
	if (set_member(fields, DeathTag))
	    printf(rfmt, Death(bp));
	if (set_member(fields, TypeTag))
	  printf(ifmt, (int) Type(bp));
	if (set_member(fields, KeyTag))
	    printf(ifmt, Key(bp));
	if (set_member(fields, AuxTag))
	    printf(rfmt, Aux(bp));
	if (set_member(fields, AuxVecTag)) {
	    printf(rfmt, AuxVec(bp)[0]);
	    printf(rfmt, AuxVec(bp)[1]);
	    printf(rfmt, AuxVec(bp)[2]);
	}
	printf("\n");
    }
}
