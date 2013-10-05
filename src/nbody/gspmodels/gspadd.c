/*
 * GSPADD.C: add several density profiles together.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "filestruct.h"
#include "gsp.h"

string defv[] = {		";Add several density profiles",
    "in=???",			";List of input GSP files",
    "out=???",			";Output GSP file",
    "VERSION=1.1",		";Josh Barnes  25 June 1997",
    NULL,
};

void gspadd(gsprof *, gsprof *);

int main(int argc, string argv[])
{
    string *inputs;
    stream istr, ostr;
    gsprof *gsp = NULL;

    initparam(argv, defv);
    inputs = burststring(getparam("in"), ", ");
    while (*inputs != NULL) {
	istr = stropen(*inputs++, "r");
	get_history(istr);
	if (gsp == NULL)
	    gsp = get_gsprof(istr);
	else
	    gspadd(gsp, get_gsprof(istr));
	strclose(istr);
    }
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_gsprof(ostr, gsp);
    strclose(ostr);
    return (0);
}

/*
 * GSPADD: add two density profiles and replace first with sum.
 */

void gspadd(gsprof *gsp1, gsprof *gsp2)
{
    int n, i;

    n = gsp1->npoint - 1;
    if (gsp1->alpha < gsp2->alpha) {
	eprintf("[gspadd: 1st profile steeper as r -> 0]\n");
	if (gsp1->density[0] < rho_gsp(gsp2, gsp1->radius[0]))
	    eprintf("[gspadd: warning: 2nd profile dominant at radius[0]]\n");
    } else if (gsp1->alpha > gsp2->alpha) {
	eprintf("[gspadd: 2nd profile steeper as r -> 0]\n");
	if (gsp1->density[0] > rho_gsp(gsp2, gsp1->radius[0]))
	    eprintf("[gspadd: warning: 1st profile dominant at radius[0]]\n");
	gsp1->alpha = gsp2->alpha;
    }
    if (gsp1->beta > gsp2->beta) {
	eprintf("[gspadd: 1st profile shallower as r -> inf]\n");
	if (gsp1->density[n] < rho_gsp(gsp2, gsp1->radius[n]))
	    eprintf("[gspadd: warning: 2nd profile dominant at radius[n]\n");
    } else if (gsp1->beta < gsp2->beta) {
	eprintf("[gspadd: 2nd profile shallower as r -> inf]\n");
	if (gsp1->density[n] > rho_gsp(gsp2, gsp1->radius[n]))
	    eprintf("[gspadd: warning: 1st profile dominant at radius[n]\n");
	gsp1->beta = gsp2->beta;
    }
    for (i = 0; i < gsp1->npoint; i++) {
	gsp1->density[i] += rho_gsp(gsp2, gsp1->radius[i]);
	gsp1->mass[i] += mass_gsp(gsp2, gsp1->radius[i]);
    }
    gsp1->mtot += gsp2->mtot;
}
