/*
 * SNAPPPM.C: render snapshot as 24-bit color ppm image.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

string defv[] = {		";Render snapshot as ppm images",
    "in=???",			";Input snapshot data file",
    "view=???",			";Viewing transformation file",
    "out=???",			";Output ppm file pattern",
    "times=all",		";Range of times to process",
    "redmap=0,1,0,1,0,1,0,1",	";Red color map values",
    "grnmap=0,0,1,1,0,0,1,1",	";Green color map values",
    "blumap=0,0,0,0,1,1,1,1",	";Blue color map values",
    "smooth=0,1,1,1,1,1,1,1",	";Smoothing count values",
    "auxinten=false",		";Scale intensity by aux data",
    "xsize=640",		";Width of image, in pixels",
    "ysize=512",		";Height of image, in pixels",
    "midval=1.0",		";Intensity mapped to midpoint",
    "slope=16.0",		";Contrast param (steps/octave)",
    "rawbits=true",		";Output image data in binary",
    "startcount=1",		";Index of first output frame",
    "VERSION=1.0",		";Josh Barnes  23 January 2000",
    NULL,
};

/* Procedure prototypes. */

void set_ctab(void);
bool get_view(stream);
bool get_xyzk(stream, string);
void snap_pix(void);
void color_hist(void);
void clear_rgb(void);
void add_rgb(int);
void load_grid(int);
void smooth_grid(int);
void put_ppm(int);
int  pixelmap(real, real, real);

/* Snapshot configuration to display. */

string bodytags[] = { KeyTag, PosTag, NULL, NULL, };

bodyptr btab = NULL;				/* array of bodies to view  */
int nbody;					/* no. of bodies in array   */
real tbody;					/* time value for body data */

/* Viewing transformation parameters. */

real aspect;					/* aspect ratio in viewfile */
real tview;					/* time of viewfile frame   */
real vmat[4][4];				/* viewing trans. matrix    */

bool auxinten;					/* scale intensities by aux */

/* Color table data. */

int ncolors;					/* number of colors defined */
real *rmap, *gmap, *bmap;			/* rgb color table values   */
int  *stab;					/* smooth iteration values  */
int *colhist = NULL;				/* histogram of colors      */

/* Image data. */

#define MAXPIX  0377				/* maximum pixel value      */

int xsize, ysize, xysize = 0;			/* size of output frames    */
real *grid = NULL;				/* smoothed intensity array */
real *rpix, *gpix, *bpix;			/* red, green, blue pixels  */

#define BUFSIZE   128				/* size of temp buffers     */
#define TIMEFUZZ  0.001				/* fuzzyness in time values */

int main(int argc, string argv[])
{
    stream istr, vstr;
    int count;
    char tbuffer[BUFSIZE];

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    vstr = stropen(getparam("view"), "r");
    get_history(vstr);
    auxinten = getbparam("auxinten");
    if (auxinten)
	bodytags[2] = AuxTag;
    layout_body(bodytags, Precision, NDIM);
    set_ctab();
    count = getiparam("startcount");
    if (! strnull(getparam("times"))) {
	if (! get_view(vstr))
	    error("%s: view file empty\n", getargv0());
	while (get_xyzk(istr, getparam("times"))) {
	    snap_pix();
	    put_ppm(count++);
	}
    } else {
	while (get_view(vstr)) {
	    if (btab == NULL || ABS(tbody - tview) > TIMEFUZZ) {
		sprintf(tbuffer, "%f", tview);
		if (! get_xyzk(istr, tbuffer))
		    error("%s: time %f not found\n", getargv0(), tview);
	    }
	    snap_pix();
	    put_ppm(count++);
	}
    }
    return (0);
}

/*
 * SET_CTAB: initialize rval, ..., sval arrays from input parameters.
 */

void set_ctab(void)
{
    string *rval, *gval, *bval, *sval;
    int i;

    rval = burststring(getparam("redmap"), ", ");
    gval = burststring(getparam("grnmap"), ", ");
    bval = burststring(getparam("blumap"), ", ");
    sval = burststring(getparam("smooth"), ", ");
    ncolors = xstrlen(rval, sizeof(int)) - 1;
    if (xstrlen(gval, sizeof(int)) != ncolors+1 ||
	  xstrlen(bval, sizeof(int)) != ncolors+1 ||
	    xstrlen(sval, sizeof(int)) != ncolors+1)
	error("%s: color maps must have same length\n", getargv0());
    rmap = (real *) allocate(ncolors * sizeof(real));
    gmap = (real *) allocate(ncolors * sizeof(real));
    bmap = (real *) allocate(ncolors * sizeof(real));
    stab = (int  *) allocate(ncolors * sizeof(int));
    for (i = 0; i < ncolors; i++) {
	rmap[i] = atof(rval[i]);
	gmap[i] = atof(gval[i]);
	bmap[i] = atof(bval[i]);
	stab[i] = atoi(sval[i]);
    }
}

/*
 * GET_VIEW: read next entry from viewfile.
 */

bool get_view(stream vstr)
{
    if (get_tag_ok(vstr, "SnapView")) {
	get_set(vstr, "SnapView");
	get_data(vstr, "Time", RealType, &tview, 0);
	get_data(vstr, "Aspect", RealType, &aspect, 0);
	get_data(vstr, "ViewMatrix", RealType, vmat, 4, 4, 0);
	get_tes(vstr, "SnapView");
	return (TRUE);
    } else
	return (FALSE);
}

/*
 * GET_XYZK: read next snapshot from input file.
 */

bool get_xyzk(stream istr, string times)
{
    static bool firstcall = TRUE;
    string intags[MaxBodyFields];

    get_history(istr);
    if (! get_snap_t(istr, &btab, &nbody, &tbody, intags, FALSE, times))
	return (FALSE);
    if (firstcall && ! set_member(intags, KeyTag))
	error("%s: %s data missing\n", getargv0(), KeyTag);
    if (firstcall && ! set_member(intags, PosTag))
	error("%s: %s data missing\n", getargv0(), PosTag);
    if (firstcall && auxinten && ! set_member(intags, AuxTag))
	error("%s: %s data missing\n", getargv0(), AuxTag);
    firstcall = FALSE;
    eprintf("[get_xyzk in %s: time = %.2f]\n", getargv0(), tbody);
    return (TRUE);
}

/*
 * SNAP_PIX: render image to rpix, gpix, bpix arrays.
 */

void snap_pix(void)
{
    int color;

    xsize = getiparam("xsize");
    ysize = getiparam("ysize");
    xysize = xsize * ysize;
    if (ABS(xsize / (real) ysize - aspect) > 0.02)
	error("%s: bad aspect ratio %f != %f\n", getargv0(),
	      xsize / (real) ysize, aspect);
    if (grid == NULL) {
	grid = (real *) allocate(xysize * sizeof(real));
	rpix = (real *) allocate(xysize * sizeof(real));
	gpix = (real *) allocate(xysize * sizeof(real));
	bpix = (real *) allocate(xysize * sizeof(real));
    }
    fprintf(stderr, "[%s: time %.2f", getargv0(), tbody);
    fflush(stderr);
    color_hist();
    fprintf(stderr, "  color_hist:");
    clear_rgb();
    for (color = 1; color < ncolors; color++) {
	fprintf(stderr, " %d", colhist[color]);
	fflush(stderr);
	if (colhist[color] > 0) {
	    load_grid(color);
	    smooth_grid(color);
	    add_rgb(color);
	}
    }
    fprintf(stderr, "]\n");
}

/*
 * COLOR_HIST: construct histogram of body colors.
 */

void color_hist(void)
{
    int i;

    if (colhist == NULL)
	colhist = (int *) allocate(ncolors * sizeof(int));
    for (i = 0; i < ncolors; i++)
	colhist[i] = 0;
    for (i = 0; i < nbody; i++)
	if (Key(NthBody(btab, i)) < ncolors)
	    colhist[Key(NthBody(btab, i))]++;
	else
	    error("%s: color out of bounds\n", getargv0());
}

/*
 * CLEAR_RGB: zero out rpix, gpix, bpix arrays.
 */

void clear_rgb(void)
{
    int i;

    for (i = 0; i < xysize; i++)
	rpix[i] = gpix[i] = bpix[i] = 0.0;
}

/*
 * ADD_RGB: transfer current image to rpix, gpix, bpix arrays.
 */

void add_rgb(int col)
{
    int i;

    if (rmap[col] > 0.0)
	for (i = 0; i < xysize; i++)
	    rpix[i] += rmap[col] * grid[i];
    if (gmap[col] > 0.0)
	for (i = 0; i < xysize; i++)
	    gpix[i] += gmap[col] * grid[i];
    if (bmap[col] > 0.0)
	for (i = 0; i < xysize; i++)
	    bpix[i] += bmap[col] * grid[i];

}

/*
 * LOAD_GRID: interpolate bodies onto pixel-grid.
 */

void load_grid(int col)
{
    int i, j, k, xj, yj, xi, yi, xk, yk;
    real pvec[4], qvec[4], xg, yg, inten, xr, yr;

    for (i = 0; i < xysize; i++)
	grid[i] = 0.0;
    for (i = 0; i < nbody; i++) {
	if (Key(NthBody(btab, i)) == col) {
	    for (j = 0; j < 4; j++)
		pvec[j] = (j < 3 ? Pos(NthBody(btab, i))[j] : 1.0);
	    for (j = 0; j < 4; j++) {
		qvec[j] = 0.0;
		for (k = 0; k < 4; k++)
		    qvec[j] += pvec[k] * vmat[k][j];
	    }
	    xg = xsize * (1.0 + qvec[0]/qvec[3]) / 2.0;
	    yg = ysize * (1.0 + qvec[1]/qvec[3]) / 2.0;
	    xj = rfloor(xg);
	    yj = rfloor(yg);
	    xi = xj - 1;
	    yi = yj - 1;
	    xk = xj + 1;
	    yk = yj + 1;
	    inten = (auxinten ? Aux(NthBody(btab, i)) : 1.0);
	    if (inten < 0.0)
		error("%s: negative intensity value\n", getargv0());
	    if (0 <= xi && xk < xsize && 0 <= yi && yk < ysize) {
		xr = xg - xj - 0.5;
		yr = yg - yj - 0.5;
		grid[yi*xsize + xi] += inten * ((0.5 - yr) * (0.5 - xr));
		grid[yi*xsize + xj] += inten * ((0.5 - yr)             );
		grid[yi*xsize + xk] += inten * ((0.5 - yr) * (0.5 + xr));
		grid[yj*xsize + xi] += inten * (             (0.5 - xr));
		grid[yj*xsize + xj] += inten * (          1.0          );
		grid[yj*xsize + xk] += inten * (             (0.5 + xr));
		grid[yk*xsize + xi] += inten * ((0.5 + yr) * (0.5 - xr));
		grid[yk*xsize + xj] += inten * ((0.5 + yr)             );
		grid[yk*xsize + xk] += inten * ((0.5 + yr) * (0.5 + xr));

	    }
	}
    }
}

/*
 * SMOOTH_GRID: smooth grid using iterated boxcar.
 */

void smooth_grid(int col)
{
    real *tmp0, *tmp1, *tmp2;
    int nsm, i, j;

    tmp0 = (real *) allocate(xsize * sizeof(real));
    tmp1 = (real *) allocate(xsize * sizeof(real));
    tmp2 = (real *) allocate(xsize * sizeof(real));
    nsm = stab[col];
    while (--nsm >= 0) {
	fprintf(stderr, "*");
	fflush(stderr);
	for (i = 0; i < xsize; i++) {
	    tmp1[i] = grid[    i    ];
	    tmp2[i] = grid[xsize + i];
	}
	for (j = 1; j < ysize-1; j++) {
	    for (i = 0; i < xsize; i++) {
		tmp0[i] = tmp1[i];
		tmp1[i] = tmp2[i];
		tmp2[i] = grid[(j+1)*xsize + i];
	    }
	    for (i = 1; i < xsize-1; i++)
		grid[j*xsize + i] =
		  (  tmp0[i-1] + 2*tmp1[i-1] +   tmp2[i-1] +
		   2*tmp0[ i ] + 4*tmp1[ i ] + 2*tmp2[ i ] +
		     tmp0[i+1] + 2*tmp1[i+1] +   tmp2[i+1]) / 16;
	}
    }
    free(tmp0);
    free(tmp1);
    free(tmp2);
}

/*
 * PUT_PPM: output 24-bit pixel data in ppm format.
 */

int locnt, hicnt;

void put_ppm(int count)
{
    real midval, slope;
    bool rawbits;
    char name[BUFSIZE];
    stream ostr;
    int y, x, i, r, g, b;

    midval = getdparam("midval");
    slope = getdparam("slope");
    rawbits = getbparam("rawbits");
    locnt = hicnt = 0;
    sprintf(name, getparam("out"), count);	/* make output file name    */
    ostr = stropen(name, "w");
    fprintf(ostr, rawbits ? "P6\n" : "P3\n");
    fprintf(ostr, "%d %d\n", xsize, ysize);
    fprintf(ostr, "%d\n", MAXPIX);
    for (y = ysize - 1; y >= 0; y--)
	for (x = 0; x < xsize; x++) {
	    i = x + xsize * y;
	    r = pixelmap(rpix[i], midval, slope);
	    g = pixelmap(gpix[i], midval, slope);
	    b = pixelmap(bpix[i], midval, slope);
	    if (! rawbits)
		fprintf(ostr, "%d %d %d\n", r, g, b);
	    else {
		fputc(r, ostr);
		fputc(g, ostr);
		fputc(b, ostr);
	    }		
	}
    fclose(ostr);
    fprintf(stderr, "[%s: locnt = %d, hicnt = %d]\n",
	    getargv0(), locnt, hicnt);
}

int pixelmap(real pix, real midval, real slope)
{
    int k;

    if (pix <= 0.0)
	return (slope > 0.0 ? 0 : MAXPIX);
    k = rfloor(slope * rlog2(pix / midval)) + 128;
    if (k < 0) {
	locnt++;
	k = 0;
    } else if (k > MAXPIX) {
	hicnt++;
	k = MAXPIX;
    }
    return (k);
}
