/*
 * filestruct_test.c: test filestruct package.
 */

#include "stdinc.h"
#include "getparam.h"
#include "filestruct.h"

string defv[] = {
    "out=foo.dat",
    "in=foo.dat",
    NULL,
};

string headline = "Mumbo Jumbo";
int nobj = 2;
real mass[2] = { 1.0, 2.0 };
real phase[2][2][3] = {
    { {  0.10,  -0.20,  0.30 }, { -3.0,  2.0, -1.0 } },
    { { -0.05,   0.10, -0.15 }, {  1.5, -1.0,  0.5 } }
};

double dmass[2];

int intvals[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

int mask[3] = { 4, -4, 0 };

void main(int argc, string argv[])
{
    stream instr, outstr;
    string tag;

    initparam(argv, defv);
    outstr = stropen(getparam("out"), "w!");
    put_set(outstr, "SnapShot");
    put_string(outstr, "Headline", headline);
    put_set(outstr, "Parameters");
    put_data(outstr, "Nobj", IntType, &nobj, 0);
    put_tes(outstr, "Parameters");
    put_set(outstr, "Particles");
    put_data(outstr, "PhaseSpace", RealType, phase, nobj, 2, 3, 0);
    put_data(outstr, "Mass", RealType, mass, nobj, 0);
    put_tes(outstr, "Particles");
    put_string(outstr, "String", "Hello World!");
    put_data(outstr, "IntVals", IntType, intvals, 12, 0);
    put_data_masked(outstr, "MaskedIntVals", IntType, intvals, 6, 0, mask);
    put_tes(outstr, "SnapShot");
    strclose(outstr);
    instr = stropen(getparam("in"), "r");
    while (get_tag_ok(instr, "History"))
	skip_item(instr);
    while (get_tag_ok(instr, "SnapShot")) {	// loop reading snapshots
	get_set(instr, "SnapShot");
	printf("  Tags in SnapShot:");
	while ((tag = next_item_tag(instr)) != NULL) {
	    printf(" %s", tag);
	    skip_item(instr);
	}
	printf("\n");
	if (get_tag_ok(instr, "String"))
	    printf("  String: %s\n", get_string(instr, "String"));
	get_set(instr, "Parameters");		// assume params exist
	get_data(instr, "Nobj", IntType, &nobj, 0);
	printf("    Nobj: %d\n", nobj);
	nobj = MIN(nobj, 2);
	get_tes(instr, "Parameters");
	if (get_tag_ok(instr, "Particles")) {	// if particles given
	    get_set(instr, "Particles");
	    if (get_tag_ok(instr, "Mass")) {
		printf("    Mass: type %s, length 0%o\n",
		       get_type(instr, "Mass"), get_length(instr, "Mass"));
		get_data(instr, "Mass", RealType, mass, nobj, 0);
		printf("    Mass: %f %f ...\n", mass[0], mass[1]);
		get_data(instr, "Mass", DoubleType, dmass, nobj, 0);
		printf("    Mass: %f %f ...\n", dmass[0], dmass[1]);
	    }
	    get_data(instr, "PhaseSpace", RealType, phase, nobj, 2, 3, 0);
	    printf("    Phase: %f %f %f ...\n",
		   phase[0][0][0], phase[0][0][1], phase[0][0][2]);
	    get_tes(instr, "Particles");
	}
	get_tes(instr, "SnapShot");
    }
    strclose(instr);
}
