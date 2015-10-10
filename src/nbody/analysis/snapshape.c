/*
 * snapshape.c: read a snapshot file and calculate shape and rms radius.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "filestruct.h"
#include "vectmath.h"
#include "phatbody.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

string defv[] = {	       	";Estimate ellipticity & r.m.s. radius",
  "in=???",			";Input snapshots, centered and sorted",
  "nbin=8",			";Number of radial bins to list",
  "listvec=false",		";If true, list eigen-vectors also",
  "VERSION=1.3",		";Josh Barnes  19 June 2015",
  NULL,
};

void eigensolve(vector, vector, vector,	real *, matrix);

int main(int argc, string argv[])
{
  stream istr;
  string bodytags[] = { PosTag, NULL }, intags[MaxBodyFields];
  bodyptr btab = NULL, bp;
  int nbody, nshell, n;
  real tnow, vals[3];
  matrix tmpm, qmat;
  vector v1, v2, v3;

  initparam(argv, defv);
  istr = stropen(getparam("in"), "r");
  get_history(istr);
  layout_body(bodytags, Precision, NDIM);
  printf("#%11s %3s %11s %11s %11s\n",
	 "time", "n", "r_rms", "c/a", "b/a");
  while (get_snap(istr, &btab, &nbody, &tnow, intags, FALSE)) {
    if (! set_member(intags, PosTag))
      error("%s: %s data missing\n", getargv0(), PosTag);
    if (nbody % getiparam("nbin") != 0)
      error("%s: nbin does not divide number of bodies\n", getargv0());
    nshell = nbody / getiparam("nbin");
    for (n = 0; n < nbody; n += nshell) {
      CLRM(qmat);
      for (bp = NthBody(btab, n); bp < NthBody(btab, n + nshell);
	   bp = NextBody(bp)) {
	OUTVP(tmpm, Pos(bp), Pos(bp));
	ADDM(qmat, qmat, tmpm);
      }
      eigensolve(v1, v2, v3, vals, qmat);
      printf(" %11.6f %3d %11.6f %11.6f %11.6f\n",
	     tnow, n / nshell, rsqrt(tracem(qmat) / nshell),
	     rsqrt(vals[2] / vals[0]), rsqrt(vals[1] / vals[0]));
      if (getbparam("listvec")) {
	printf("#\t\t\t\t\t\t\t%8.5f  %8.5f  %8.5f\n", v1[0], v1[1], v1[2]);
	printf("#\t\t\t\t\t\t\t%8.5f  %8.5f  %8.5f\n", v2[0], v2[1], v2[2]);
	printf("#\t\t\t\t\t\t\t%8.5f  %8.5f  %8.5f\n", v3[0], v3[1], v3[2]);
      }
    }
  }
  return (0);
}

void eigensolve(vector vec1, vector vec2, vector vec3,
		real *vals, matrix symmat)
{
  int i, j;
  double data[9];
  gsl_matrix_view mat;
  gsl_vector_view vec;
  gsl_vector *eigval = gsl_vector_alloc(3);
  gsl_matrix *eigvec = gsl_matrix_alloc(3, 3);
  gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(3);

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      data[3*i + j] = symmat[i][j];
  mat = gsl_matrix_view_array(data, 3, 3);
  gsl_eigen_symmv(&mat.matrix, eigval, eigvec, work);
  gsl_eigen_symmv_free(work);
  gsl_eigen_symmv_sort(eigval, eigvec, GSL_EIGEN_SORT_VAL_DESC);
  for (i = 0; i < 3; i++)
    vals[i] = gsl_vector_get(eigval, i);
  vec = gsl_matrix_column(eigvec, 0);
  for (i = 0; i < 3; i++)
    vec1[i] = gsl_vector_get(&vec.vector, i);
  vec = gsl_matrix_column(eigvec, 1);
  for (i = 0; i < 3; i++)
    vec2[i] = gsl_vector_get(&vec.vector, i);
  vec = gsl_matrix_column(eigvec, 2);
  for (i = 0; i < 3; i++)
    vec3[i] = gsl_vector_get(&vec.vector, i);
  gsl_vector_free(eigval);
  gsl_matrix_free(eigvec);
}
