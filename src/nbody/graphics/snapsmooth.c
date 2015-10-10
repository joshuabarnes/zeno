/*
 * snapsmooth.c: smooth N-body/SPH values onto grid.
 */

#include "stdinc.h"
#include "strset.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Smooth N-body/SPH values onto grid",
  "in=???",			";Input snapshot data file name",
  "pgmout=",			";Output pgm file name pattern",
  "auxout=",			";Output values at body positions",
  "times=all",			";Range of times to process",
  "value=rho",			";Options are: rho, aux, bright, rgb, RGB",
  "threedim=false",		";If TRUE, do 3-D interpolation",
  "zval=0.0",			";Z coordinate of sampling grid",
  "scale=1/256",		";Spacing between grid points",
  "xsize=640",			";Width of grid, in grid cells",
  "ysize=512",			";Height of grid, in grid cells",
  "logmap=true",		";If TRUE, take log of values",
  "midval=1.0",			";Value mapped to midpoint",
  "slope=16.0",			";Contrast (steps per factor of 2)",
  "depth=1",			";Depth of image (1 or 2 bytes)",
  "startcount=0",		";Number of first output frame",
  "stepcount=1",		";Increment between output frames",
  "VERSION=2.0",		";Josh Barnes  9 May 2013",
  NULL,
};

// Procedure prototypes.

bool get_bodies(stream, string, bool);
void smooth_rho(real *, bool, real);
void smooth_aux(real *, bool, real);
void smooth_bright(real *, int);
void put_pgm(string, int, string, bool, real, real, int, real *);
void put_ppm(string, int, bool, real, real, int, real *, real *, real *);
int pixmap(real gridval, bool logmap, real slope, real midval, int maxpix);
void pixclip(int *pixval, int *lowcount, int *highcount, int maxpix);
void set_aux_value(real *);

// Global quantities.

string *bodytags;			// list of fields for bodies
bodyptr btab = NULL;			// array of bodies to smooth
int nbody;				// noumber of bodies in array
real tsnap;				// time value for body data
bool taudef;				// true if Tau() data input	

real scale;				// grid point spacing
int xsize, ysize, xysize;		// size of output frames
real *grid1, *grid2, *grid3;		// interpolated values

#define MAXPIX_1   255			// max pixel value (1 byte, unsigned)
#define MAXPIX_2 65535			// max pixel value (2 byte, unsigned)
#define MAXPIXm2 32767			// max pixel value (2 byte, signed)

int main(int argc, string argv[])
{
  string value;
  stream instr, auxstr = NULL;
  bool rhovalue = FALSE, threedim, logmap;
  real midval, slope;
  int count, depth;

  initparam(argv, defv);
  instr = stropen(getparam("in"), "r");
  get_history(instr);
  value = getparam("value");
  if (streq(value, "rho")) {
    bodytags = set_cons(MassTag, PosTag, SmoothTag,
			strnull(getparam("auxout")) ? NULL : AuxTag, NULL);
    rhovalue = TRUE;
  } else if (streq(value, "aux")) {
    bodytags = set_cons(MassTag, PosTag, SmoothTag, AuxTag, RhoTag, NULL);
  } else if (streq(value, "bright")) {
    bodytags = set_cons(PosTag, SmoothTag, AuxTag, TauTag, NULL);
    if (getbparam("threedim"))
      error("%s: value=bright requires threedim=false\n", getprog());
  } else if (streq(value, "rgb") || streq(value, "RGB")) {
    bodytags = set_cons(PosTag, SmoothTag, AuxVecTag, TauTag, NULL);
    if (getbparam("threedim"))
      error("%s: value=%s requires threedim=false\n", getprog(), value);
  } else
    error("%s: value=%s undefined\n", getprog(), value);
  layout_body(bodytags, Precision, NDIM);
  threedim = getbparam("threedim");
  scale = getdparam("scale");
  xsize = getiparam("xsize");
  ysize = getiparam("ysize");
  logmap = getbparam("logmap");
  midval = getdparam("midval");
  slope = getdparam("slope");
  depth = getiparam("depth");
  count = getiparam("startcount");
  xysize = xsize * ysize;
  grid1 = (real *) allocate(xysize * sizeof(real));
  grid2 = (real *) allocate(xysize * sizeof(real));
  grid3 = (real *) allocate(xysize * sizeof(real));
  if (! strnull(getparam("auxout"))) {
    if (threedim)
      error("%s: auxout requires threedim=false\n", getprog());
    auxstr = stropen(getparam("auxout"), "w");
    put_history(auxstr);
  }
  while (get_bodies(instr, getparam("times"), rhovalue)) {
    if (streq(value, "rho")) {
      smooth_rho(grid1, threedim, getdparam("zval"));
      put_pgm(getparam("pgmout"), count, "", logmap, midval, slope, depth, grid1);
    } else if (streq(value, "aux")) {
      smooth_aux(grid1, threedim, getdparam("zval"));
      put_pgm(getparam("pgmout"), count, "", logmap, midval, slope, depth, grid1);
    } else if (streq(value, "bright")) {
      smooth_bright(grid1, -1);
      put_pgm(getparam("pgmout"), count, "", logmap, midval, slope, depth, grid1);
    } else {
      smooth_bright(grid1, 0);
      smooth_bright(grid2, 1);
      smooth_bright(grid3, 2);
      if (streq(value, "rgb")) {
	put_pgm(getparam("pgmout"), count, "r", logmap, midval, slope, depth, grid1);
	put_pgm(getparam("pgmout"), count, "g", logmap, midval, slope, depth, grid2);
	put_pgm(getparam("pgmout"), count, "b", logmap, midval, slope, depth, grid3);
      } else
	put_ppm(getparam("pgmout"), count, logmap, midval, slope, depth,
		grid1, grid2, grid3);
    } 
    count += getiparam("stepcount");
    if (auxstr != NULL) {
      set_aux_value(grid1);
      put_snap(auxstr, &btab, &nbody, &tsnap, bodytags);
    }
  }
  return (0);
}

//  get_bodies: read snapshot from input file.
//  __________________________________________

bool get_bodies(stream instr, string times, bool rhovalue)
{
  static bool firstcall = TRUE;
  string intags[MaxBodyFields];
  int i;

  get_history(instr);
  if (! get_snap_t(instr, &btab, &nbody, &tsnap, intags, FALSE, times))
    return (FALSE);
  if (firstcall) {				// check reqd. data fields
    for (i = 0; bodytags[i] != NULL; i++) {
      if (streq(bodytags[i], AuxTag) && rhovalue)
	continue;				// Aux data not required
      if (streq(bodytags[i], TauTag))
	continue;				// Tau data not required
      if (! set_member(intags, bodytags[i]))
	error("%s: %s data missing\n", getprog(), bodytags[i]);
    }
    firstcall = FALSE;				// don't check again
    taudef = set_member(intags, TauTag);
  }
  eprintf("[%s.get_bodies: time = %.2f]\n", getprog(), tsnap);
  return (TRUE);
}

//  smooth_rho: interpolate density onto sample grid.
//  _________________________________________________

void smooth_rho(real *grid, bool threedim, real zval)
{
  int i, j, i1, i2, j1, j2;
  bodyptr b;
  real hsqr, dzsqr, C, y1, y2, dysqr, x1, x2, dxsqr, q, W;

  for (i = 0; i < xysize; i++)
    grid[i] = 0.0;
  for (b = btab; b < NthBody(btab, nbody); b = NextBody(b)) {
    hsqr = rsqr(Smooth(b));
    dzsqr = (threedim ? rsqr(Pos(b)[2] - zval) : 0.0);
    if (dzsqr <= 4 * hsqr) {
      C = (threedim ? 1/(PI*rqbe(Smooth(b))) : 10/(7*PI*rsqr(Smooth(b))));
      y1 = Pos(b)[1] - rsqrt(4 * hsqr - dzsqr);
      y2 = Pos(b)[1] + rsqrt(4 * hsqr - dzsqr);
      j1 = rfloor(y1 / scale + (ysize-1) / 2.0);
      j2 = rceil(y2 / scale + (ysize-1) / 2.0);
      for (j = MAX(j1, 0); j <= MIN(j2, ysize-1); j++) {
	dysqr = rsqr(Pos(b)[1] - scale * (j - (ysize-1)/2.0));
	if (dzsqr + dysqr <= 4 * hsqr) {
	  x1 = Pos(b)[0] - rsqrt(4 * hsqr - dzsqr - dysqr);
	  x2 = Pos(b)[0] + rsqrt(4 * hsqr - dzsqr - dysqr);
	  i1 = rfloor(x1 / scale + (xsize-1) / 2.0);
	  i2 = rceil(x2 / scale + (xsize-1) / 2.0);
	  for (i = MAX(i1, 0); i <= MIN(i2, xsize-1); i++) {
	    dxsqr = rsqr(Pos(b)[0] - scale * (i - (xsize-1)/2.0));
	    if (dzsqr + dysqr + dxsqr <= 4 * hsqr) {
	      q = rsqrt((dzsqr + dysqr + dxsqr) / hsqr);
	      W = C * (q > 1 ? rqbe(2-q)/4.0 : 1 - 1.5*rsqr(q) + 0.75*rqbe(q));
	      grid[i + xsize*j] += W * Mass(b);
	    }
	  }
	}
      }
    }
  }
}

//  smooth_aux: interpolate aux onto sample grid.
//  _____________________________________________

void smooth_aux(real *grid, bool threedim, real zval)
{
  int i, j, i1, i2, j1, j2;
  bodyptr b;
  real hsqr, dzsqr, C, y1, y2, dysqr, x1, x2, dxsqr, q, W;

  for (i = 0; i < xysize; i++)
    grid[i] = 0.0;
  for (b = btab; b < NthBody(btab, nbody); b = NextBody(b)) {
    hsqr = rsqr(Smooth(b));
    dzsqr = (threedim ? rsqr(Pos(b)[2] - zval) : 0.0);
    if (dzsqr <= 4 * hsqr) {
      C = (threedim ? 1/(PI*rqbe(Smooth(b))) : 10/(7*PI*rsqr(Smooth(b))));
      y1 = Pos(b)[1] - rsqrt(4 * hsqr - dzsqr);
      y2 = Pos(b)[1] + rsqrt(4 * hsqr - dzsqr);
      j1 = rfloor(y1 / scale + (ysize-1) / 2.0);
      j2 = rceil(y2 / scale + (ysize-1) / 2.0);
      for (j = MAX(j1, 0); j <= MIN(j2, ysize-1); j++) {
	dysqr = rsqr(Pos(b)[1] - scale * (j - (ysize-1)/2.0));
	if (dzsqr + dysqr <= 4 * hsqr) {
	  x1 = Pos(b)[0] - rsqrt(4 * hsqr - dzsqr - dysqr);
	  x2 = Pos(b)[0] + rsqrt(4 * hsqr - dzsqr - dysqr);
	  i1 = rfloor(x1 / scale + (xsize-1) / 2.0);
	  i2 = rceil(x2 / scale + (xsize-1) / 2.0);
	  for (i = MAX(i1, 0); i <= MIN(i2, xsize-1); i++) {
	    dxsqr = rsqr(Pos(b)[0] - scale * (i - (xsize-1)/2.0));
	    if (dzsqr + dysqr + dxsqr <= 4 * hsqr) {
	      q = rsqrt((dzsqr + dysqr + dxsqr) / hsqr);
	      W = C * (q > 1 ? rqbe(2-q)/4.0 : 1 - 1.5*rsqr(q) + 0.75*rqbe(q));
	      grid[i + xsize*j] += W * Mass(b) * Aux(b) / Rho(b);
	    }
	  }
	}
      }
    }
  }
}

//  smooth_bright: interpolate brightness onto sample grid.
//  _______________________________________________________

void smooth_bright(real *grid, int indx)
{
  int i, j, i1, i2, j1, j2;
  bodyptr b;
  real hsqr, C, y1, y2, dysqr, x1, x2, dxsqr, q, W;

  for (i = 0; i < xysize; i++)
    grid[i] = 0.0;
  for (b = btab; b < NthBody(btab, nbody); b = NextBody(b)) {
    hsqr = rsqr(Smooth(b));
    C = 10/(7*PI*rsqr(Smooth(b)));
    y1 = Pos(b)[1] - rsqrt(4 * hsqr);
    y2 = Pos(b)[1] + rsqrt(4 * hsqr);
    j1 = rfloor(y1 / scale + (ysize-1) / 2.0);
    j2 = rceil(y2 / scale + (ysize-1) / 2.0);
    for (j = MAX(j1, 0); j <= MIN(j2, ysize-1); j++) {
      dysqr = rsqr(Pos(b)[1] - scale * (j - (ysize-1)/2.0));
      if (dysqr <= 4 * hsqr) {
	x1 = Pos(b)[0] - rsqrt(4 * hsqr - dysqr);
	x2 = Pos(b)[0] + rsqrt(4 * hsqr - dysqr);
	i1 = rfloor(x1 / scale + (xsize-1) / 2.0);
	i2 = rceil(x2 / scale + (xsize-1) / 2.0);
	for (i = MAX(i1, 0); i <= MIN(i2, xsize-1); i++) {
	  dxsqr = rsqr(Pos(b)[0] - scale * (i - (xsize-1)/2.0));
	  if (dysqr + dxsqr <= 4 * hsqr) {
	    q = rsqrt((dysqr + dxsqr) / hsqr);
	    W = C * (q > 1 ? rqbe(2-q)/4.0 : 1 - 1.5*rsqr(q) + 0.75*rqbe(q));
	    if (taudef && Tau(b) > 0.0)
	      grid[i + xsize*j] *= rexp(- W * Tau(b));
	    grid[i + xsize*j] += W * (indx == -1 ? Aux(b) : AuxVec(b)[indx]);
	  }
	}
      }
    }
  }
}

//  put_pgm: output pixel data in linear pgm format.
//  ________________________________________________

#define BUFSIZE  256			// size of temp buffer

void put_pgm(string out, int count, string color,
	     bool logmap, real midval, real slope, int depth, real *grid)
{
  int maxpix, lowcount, highcount, j, i, pixij;
  char name[BUFSIZE];
  stream pgmstr;

  if (! strnull(out)) {
    if (depth != 1 && depth != 2 && depth != -2)
      error("%s: illegal image depth: %d\n", getprog(), depth);
    maxpix = (depth == 1 ? MAXPIX_1 : depth == 2 ? MAXPIX_2 : MAXPIXm2);
    snprintf(name, BUFSIZE, out, count, color);	// make output file name
    pgmstr = stropen(name, "w");
    fprintf(pgmstr, "P5\n");			// output PGM header
    fprintf(pgmstr, "# Time = %.6f\n", tsnap);
    fprintf(pgmstr, "%d %d\n", xsize, ysize);
    fprintf(pgmstr, "%d\n", maxpix);
    lowcount = highcount = 0;			// count low, high pixels
    for (j = ysize - 1; j >= 0; j--)
      for (i = 0; i < xsize; i++) {
	pixij = pixmap(grid[i + xsize*j], logmap, slope, midval, maxpix);
	pixclip(&pixij, &lowcount, &highcount, maxpix);
	if (depth != 1)				// put most sig. byte 1st
	  fputc(MAXPIX_1 & (pixij >> 8), pgmstr);
	fputc(MAXPIX_1 & pixij, pgmstr);
      }
    fclose(pgmstr);
    eprintf(lowcount + highcount < xysize ?
	    "[%s.put_pgm: %d low pixels, %d high pixels]\n" :
	    "[%s.put_pgm: WARNING: %d low pixels, %d high pixels]\n",
	    getprog(), lowcount, highcount);
  }
}

//  put_ppm: output pixel data in linear ppm format.
//  ________________________________________________

void put_ppm(string out, int count, bool logmap, real midval, real slope,
	     int depth, real *grid1, real *grid2, real *grid3)
{
  int maxpix, lowcount, highcount, j, i, pixij;
  char name[BUFSIZE];
  stream ppmstr;

  if (! strnull(out)) {
    if (depth != 1 && depth != 2 && depth != -2)
      error("%s: illegal image depth: %d\n", getprog(), depth);
    maxpix = (depth == 1 ? MAXPIX_1 : depth == 2 ? MAXPIX_2 : MAXPIXm2);
    snprintf(name, BUFSIZE, out, count);	// make output file name
    ppmstr = stropen(name, "w");
    fprintf(ppmstr, "P6\n");			// output PPM header
    fprintf(ppmstr, "# Time = %.6f\n", tsnap);
    fprintf(ppmstr, "%d %d\n", xsize, ysize);
    fprintf(ppmstr, "%d\n", maxpix);
    lowcount = highcount = 0;			// count low, high pixels
    for (j = ysize - 1; j >= 0; j--)
      for (i = 0; i < xsize; i++) {
	pixij = pixmap(grid1[i + xsize*j], logmap, slope, midval, maxpix);
	pixclip(&pixij, &lowcount, &highcount, maxpix);
	if (depth != 1)				// put most sig. byte 1st
	  fputc(MAXPIX_1 & (pixij >> 8), ppmstr);
	fputc(MAXPIX_1 & pixij, ppmstr);
	pixij = pixmap(grid2[i + xsize*j], logmap, slope, midval, maxpix);
	pixclip(&pixij, &lowcount, &highcount, maxpix);
	if (depth != 1)				// put most sig. byte 1st
	  fputc(MAXPIX_1 & (pixij >> 8), ppmstr);
	fputc(MAXPIX_1 & pixij, ppmstr);
	pixij = pixmap(grid3[i + xsize*j], logmap, slope, midval, maxpix);
	pixclip(&pixij, &lowcount, &highcount, maxpix);
	if (depth != 1)				// put most sig. byte 1st
	  fputc(MAXPIX_1 & (pixij >> 8), ppmstr);
	fputc(MAXPIX_1 & pixij, ppmstr);
      }
    fclose(ppmstr);
    eprintf(lowcount + highcount < 3 * xysize ?
	    "[%s.put_ppm: %d low pixels, %d high pixels]\n" :
	    "[%s.put_ppm: WARNING: %d low pixels, %d high pixels]\n",
	    getprog(), lowcount, highcount);
  }
}

// pixmap: transform real grid value to integer pixel value.
// _________________________________________________________

int pixmap(real gridval, bool logmap, real slope, real midval, int maxpix)
{
  if (logmap) {					// log transform data?
    if (gridval > 0.0)
      return (rfloor(slope * rlog2(gridval / midval) + (maxpix + 1) / 2));
    else
      return (slope > 0.0 ? -1 : maxpix + 1);	// handle zero or negative
  } else
    return (rfloor(slope * (gridval / midval - 1) + (maxpix + 1) / 2));
}

// pixclip: limit pixel value to allowed range, and count exceptions.
// __________________________________________________________________

void pixclip(int *pixval, int *lowcount, int *highcount, int maxpix)
{
  if (*pixval < 0) {
    (*lowcount)++;
    *pixval = 0;
  } else if (*pixval > maxpix) {
    (*highcount)++;
    *pixval = maxpix;
  }
}

//  set_aux_value: set aux() value of each body from grid.
//  ______________________________________________________

void set_aux_value(real *grid)
{
  bodyptr b;
  int i, j;

  for (b = btab; b < NthBody(btab, nbody); b = NextBody(b)) {
    i = rfloor(Pos(b)[0] / scale + xsize / 2.0);
    j = rfloor(Pos(b)[1] / scale + ysize / 2.0);
    if (0 <= i && i < xsize && 0 <= j && j < ysize)
      Aux(b) = grid[i + xsize*j];
    else
      Aux(b) = -1.0;
  }
}
