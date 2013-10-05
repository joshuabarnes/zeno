/*
 * VIEWFILE: construct view file.
 */


#include "stdinc.h"
#include "getparam.h"
#include "filestruct.h"

string defv[] = {		";Construct view file",
    "out=???",			";Output SnapView file",
    "aspect=1.25",		";Aspect ratio (x/y)",
    "xfield=2.0",		";Horizontal field of view",
    "VERSION=1.0",		";Josh Barnes  3 July 1999",
    NULL,
};

int main(int argc, string argv[])
{
    real tnow = 0.0, aspect, xfield, vmat[4][4];
    int i, j;
    stream vstr;

    initparam(argv, defv);
    aspect = getdparam("aspect");
    xfield = getdparam("xfield");
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    vmat[i][j] = 0.0;
    vmat[0][0] = 1.0;
    vmat[1][1] = aspect;
    vmat[3][3] = xfield / 2.0;
    vstr = stropen(getparam("out"), "w");
    put_history(vstr);
    put_set(vstr, "SnapView");
    put_data(vstr, "Time", RealType, &tnow, 0);
    put_data(vstr, "Aspect", RealType, &aspect, 0);
    put_data(vstr, "ViewMatrix", RealType, vmat, 4, 4, 0);
    put_tes(vstr, "SnapView");
    strclose(vstr);
    return (0);
}
