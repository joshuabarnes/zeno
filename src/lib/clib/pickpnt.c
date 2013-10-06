/*
 * pickpnt.c: pick points at random in various distributions.
 */

#include "stdinc.h"
#include "mathfns.h"

//  pickshell: pick point on shell.
//  _______________________________

void pickshell(real vec[], int ndim, real rad)
{
  real rsq, rscale;
  int i;

  do {
    rsq = 0.0;
    for (i = 0; i < ndim; i++) {
      vec[i] = xrandom(-1.0, 1.0);
      rsq = rsq + vec[i] * vec[i];
    }
  } while (rsq > 1.0);
  rscale = rad / rsqrt(rsq);
  for (i = 0; i < ndim; i++)
    vec[i] = vec[i] * rscale;
}

//  pickball: pick point within ball.
//  _________________________________

void pickball(real vec[], int ndim, real rad) 
{
  real rsq;
  int i;

  do {
    rsq = 0.0;
    for (i = 0; i < ndim; i++) {
      vec[i] = xrandom(-1.0, 1.0);
      rsq = rsq + vec[i] * vec[i];
    }
  } while (rsq > 1.0);
  for (i = 0; i < ndim; i++)
    vec[i] = vec[i] * rad;
}

//  pickbox: pick point within box.
//  _______________________________

void pickbox(real vec[], int ndim, real size) 
{
  int i;

  for (i = 0; i < ndim; i++)
    vec[i] = xrandom(- size, size);
}
