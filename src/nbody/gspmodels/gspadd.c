/*
 * gspadd.c: add several density profiles together.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "filestruct.h"
#include "gsp.h"

string defv[] = {		";Add several density profiles",
  "in=???",			";List of input GSP files",
  "out=???",			";Output GSP file",
  "VERSION=2.0",		";Josh Barnes  21 June 2017",
  NULL,
};

void gsp_add(gsprof *, gsprof *);

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
      gsp = gsp_read(istr);
    else
      gsp_add(gsp, gsp_read(istr));
    strclose(istr);
  }
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  gsp_write(ostr, gsp);
  fflush(NULL);
  return 0;
}

//  gsp_add: add two density profiles and replace first with sum.
//  _____________________________________________________________

void gsp_add(gsprof *gsp1, gsprof *gsp2)
{
  int N = gsp1->npoint - 1;

  if (gsp1->alpha < gsp2->alpha) {
    eprintf("[%s: 1st profile steeper as r -> 0]\n", getprog());
    if (gsp1->density[0] < gsp_rho(gsp2, gsp1->radius[0]))
      eprintf("[%s: warning: 2nd profile dominant at radius[0]]\n",
	      getprog());
  } else if (gsp1->alpha > gsp2->alpha) {
    eprintf("[%s: 2nd profile steeper as r -> 0]\n", getprog());
    if (gsp1->density[0] > gsp_rho(gsp2, gsp1->radius[0]))
      eprintf("[%s: warning: 1st profile dominant at radius[0]]\n",
	      getprog());
    gsp1->alpha = gsp2->alpha;			// steeper prof wins as r -> 0
  }
  if (gsp1->beta > gsp2->beta) {
    eprintf("[%s: 1st profile shallower as r -> inf]\n", getprog());
    if (gsp1->density[N] < gsp_rho(gsp2, gsp1->radius[N]))
      eprintf("[%s: warning: 2nd profile dominant at radius[N]]\n",
	      getprog());
  } else if (gsp1->beta < gsp2->beta) {
    eprintf("[%s: 2nd profile shallower as r -> inf]\n", getprog());
    if (gsp1->density[N] > gsp_rho(gsp2, gsp1->radius[N]))
      eprintf("[%s: warning: 1st profile dominant at radius[N]]\n",
	      getprog());
    gsp1->beta = gsp2->beta;			// shal. prof wins as r -> INF
  }
  for (int i = 0; i < gsp1->npoint; i++) {
    gsp1->density[i] += gsp_rho(gsp2, gsp1->radius[i]);
    gsp1->mass[i] += gsp_mass(gsp2, gsp1->radius[i]);
  }
  gsp1->mtot += gsp2->mtot;
}
