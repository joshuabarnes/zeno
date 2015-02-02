/*
 * snaprect.c: transform to coordinates diagonalizing weight tensor.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"
#include "buildmap.h"
#include "snapcenter.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <unistd.h>

string defv[] = {		";Rotate snap to diagonalize weight tensor",
  "in=???",                     ";Input snapshot file name",
  "out=???",                    ";Output snapshot file name",
  "times=all",                  ";Range of times to process",
  "weight=1.0",			";Expression (C code) for body weight.",
				";May use these values (if given in input):",
				  SNAPMAP_BODY_VARS ".",
  "require=",			";Input items required",
  "produce=",			";Output items produced",
  "passall=true",		";If true, pass on input data",
  "seed=",			";Generator seed for random values",
  "VERSION=2.2",                ";Josh Barnes  2 February 2015",
  NULL,
};

void snaprect(bodyptr, int);			// transform body array
void eigenvect(vector, vector, vector, matrix);	// find eigenvectors of matrix
void printvect(string, vector);			// print vector to stderr
stream execmap(string);				// start snapmap process
void del_tag(string *, string *, string);	// remove tag from list

string names[2] = { "Weight", NULL };
string exprs[2] = { NULL,     NULL };
string types[2] = { RealType, NULL };

#define WeightField  phatbody[NewBodyFields+0]
#define Weight(b)  SelectReal(b, WeightField.offset)

int main(int argc, string argv[])
{
  string prog, itags[MaxBodyFields], otags[MaxBodyFields];
  stream xstr, ostr;
  bodyptr btab = NULL;
  int nbody;
  real tnow;

  initparam(argv, defv);
  exprs[0] = getparam("weight");
  prog = mktemp((string) copxstr("/tmp/sm_XXXXXX", sizeof(char)));
  buildmap(prog, names, exprs, types, NULL, Precision, NDIM, TRUE);
  xstr = execmap(prog);
  if (get_tag_ok(xstr, "History"))
    skip_item(xstr);
  get_history(xstr);
  ostr = stropen(getparam("out"), "w");
  put_history(ostr);
  new_field(&WeightField, RealType, "Weight");
  new_field(&WeightField + 1, NULL, NULL);
  while (get_snap(xstr, &btab, &nbody, &tnow, itags, TRUE)) {
    snaprect(btab, nbody);
    del_tag(otags, itags, "Weight");
    put_snap(ostr, &btab, &nbody, &tnow, otags);
  }
  strclose(ostr);
  if (unlink(prog) != 0)
    error("%s: can't unlink %s\n", getargv0(), prog);
  return (0);
}

void snaprect(bodyptr btab, int nbody)
{
  matrix qmat, tmpm;
  bodyptr bp;
  vector frame[3], tmpv;
  static vector oldframe[3] =
    { { 1.0, 0.0, 0.0, }, { 0.0, 1.0, 0.0, }, { 0.0, 0.0, 1.0, }, };
  int i;

  snapcenter(btab, nbody, WeightField.offset);
  CLRM(qmat);
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    OUTVP(tmpm, Pos(bp), Pos(bp));
    MULMS(tmpm, tmpm, Weight(bp));
    ADDM(qmat, qmat, tmpm);
  }
  eigenvect(frame[0], frame[1], frame[2], qmat);
  if (dotvp(oldframe[0], frame[0]) < 0.0)
    MULVS(frame[0], frame[0], -1.0);
  if (dotvp(oldframe[2], frame[2]) < 0.0)
    MULVS(frame[2], frame[2], -1.0);
  CROSSVP(frame[1], frame[2], frame[0]);
  printvect("e_x:", frame[0]);
  printvect("e_y:", frame[1]);
  printvect("e_z:", frame[2]);
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    if (PosField.offset != BadOffset) {
      for (i = 0; i < NDIM; i++)
	tmpv[i] = dotvp(Pos(bp), frame[i]);
      SETV(Pos(bp), tmpv);
    }
    if (VelField.offset != BadOffset) {
      for (i = 0; i < NDIM; i++)
	tmpv[i] = dotvp(Vel(bp), frame[i]);
      SETV(Vel(bp), tmpv);
    }
    if (AccField.offset != BadOffset) {
      for (i = 0; i < NDIM; i++)
	tmpv[i] = dotvp(Acc(bp), frame[i]);
      SETV(Acc(bp), tmpv);
    }
    if (AuxVecField.offset != BadOffset) {
      for (i = 0; i < NDIM; i++)
	tmpv[i] = dotvp(AuxVec(bp), frame[i]);
      SETV(AuxVec(bp), tmpv);
    }
  }
  for (i = 0; i < NDIM; i++)
    SETV(oldframe[i], frame[i]);
}    

void eigenvect(vector vec1, vector vec2, vector vec3, matrix symmat)
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

void printvect(string name, vector vec)
{
  eprintf("[%s: %12s  %10.5f  %10.5f  %10.5f]\n",
	  getargv0(), name, vec[0], vec[1], vec[2]);
}

//  execmap: start snapmap subprocess, and return snapmap output stream.
//  ____________________________________________________________________

stream execmap(string prog)
{
  int handle[2];
  char handbuf[32], produce[512];

  pipe(handle);
  if (fork() == 0) {                           // if this is child process
    close(handle[0]);
    sprintf(handbuf, "-%d", handle[1]);
    sprintf(produce, "%s,Weight", getparam("produce"));
    execl(prog, getprog(), getparam("in"), handbuf, getparam("times"),
	  getparam("require"), produce, getparam("passall"),
	  getparam("seed"), NULL);
    error("%s: execl %s failed\n", getprog(), prog);
  }
  close(handle[1]);
  sprintf(handbuf, "-%d", handle[0]);
  return (stropen(handbuf, "r"));
}

void del_tag(string *olist, string *ilist, string tag)
{
  string *op, *ip;

  for (op = olist, ip = ilist; *ip != NULL; ip++)
    if (! streq(*ip, tag))
      *op++ = *ip;
  *op = NULL;
}
