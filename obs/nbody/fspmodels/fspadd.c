/*
 * FSPADD.C: add several density profiles together.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "filestruct.h"
#include "fsp.h"

string defv[] = {		";Add several density profiles",
    "in=???",			";List of input FSP files",
    "out=???",			";Output FSP file",
    "VERSION=1.1",		";Josh Barnes  25 June 1997",
    NULL,
};

void fspadd(fsprof *, fsprof *);

int main(int argc, string argv[])
{
    string *inputs;
    stream istr, ostr;
    fsprof *fsp = NULL;

    initparam(argv, defv);
    inputs = burststring(getparam("in"), ", ");
    while (*inputs != NULL) {
	istr = stropen(*inputs++, "r");
	get_history(istr);
	if (fsp == NULL)
	    fsp = get_fsprof(istr);
	else
	    fspadd(fsp, get_fsprof(istr));
	strclose(istr);
    }
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_fsprof(ostr, fsp);
    strclose(ostr);
    return (0);
}

/*
 * FSPADD: add two density profiles and replace first with sum.
 */

void fspadd(fsprof *fsp1, fsprof *fsp2)
{
    int n, i;

    n = fsp1->npoint - 1;
    if (fsp1->alpha < fsp2->alpha) {
	eprintf("[fspadd: 1st profile steeper as r -> 0]\n");
	if (fsp1->density[0] < rho_fsp(fsp2, fsp1->radius[0]))
	    eprintf("[fspadd: WARNING: 2nd profile dominant at radius[0]]\n");
    } else if (fsp1->alpha > fsp2->alpha) {
	eprintf("[fspadd: 2nd profile steeper as r -> 0]\n");
	if (fsp1->density[0] > rho_fsp(fsp2, fsp1->radius[0]))
	    eprintf("[fspadd: WARNING: 1st profile dominant at radius[0]]\n");
	fsp1->alpha = fsp2->alpha;
    }
    if (fsp1->beta > fsp2->beta) {
	eprintf("[fspadd: 1st profile shallower as r -> inf]\n");
	if (fsp1->density[n] < rho_fsp(fsp2, fsp1->radius[n]))
	    eprintf("[fspadd: WARNING: 2nd profile dominant at radius[n]\n");
    } else if (fsp1->beta < fsp2->beta) {
	eprintf("[fspadd: 2nd profile shallower as r -> inf]\n");
	if (fsp1->density[n] > rho_fsp(fsp2, fsp1->radius[n]))
	    eprintf("[fspadd: WARNING: 1st profile dominant at radius[n]\n");
	fsp1->beta = fsp2->beta;
    }
    for (i = 0; i < fsp1->npoint; i++) {
	fsp1->density[i] += rho_fsp(fsp2, fsp1->radius[i]);
	fsp1->mass[i] += mass_fsp(fsp2, fsp1->radius[i]);
    }
    fsp1->mtot += fsp2->mtot;
}
