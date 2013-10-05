/*
 * SNAPVIEW.C: interactively view a SnapShot file.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

#if defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <string.h>
#include <assert.h>

#define XX ", "

string defv[] = {		";Interactively view SnapShot data.",
				";Mouse controls viewing parameters;",
				";Keys <esc>, <sp>, a, b, m, p, r, R, s,",
				";t, v, z control program operations.",
    "in=???",			";SnapShot file to display",
    "times=all",		";Range of times to display",
    "refscale=1.0",		";Linear size of reference cube",
    "colordata=",		";Color points using scalar data.",
				";Options: " PhiTag XX SmoothTag ",",
				";" RhoTag XX EntFuncTag XX UinternTag ",",
				";" UdotIntTag XX UdotRadTag XX UdotVisTag ",",
				";" TauTag XX BirthTag XX DeathTag XX
				    AuxTag ".", 
    "dopcolor=false",		";Color points using LOS velocities",
    "vectordata=",		";Show lines using vector data.",
				";Options: " VelTag XX AccTag XX AuxVecTag ".",
    "maxfast=16384",		";Max. points shown in fast mode",
    "defcolors=00ff00,ffffff",	";Default colors for points and cube",
    "viewsize=640x512",		";Size of viewing area in pixels",
    "viewfile=view.dat",	";Output file for viewing parameters",
    "ident=",			";Text string to label window",
    "VERSION=2.0",		";Josh Barnes  21 May 2012",
    NULL,
};

bool getdata(void);
void initgraphics(string);
void reshape(int, int);
void display(void);
void setview(void);
void render(void);
void annotate(void);
void showtext(string, int, int);
void mouse(int, int, int, int);
void motion(int, int);
void keyboard(unsigned char, int, int);
void writeview(void);

string bodytags[] = { PosTag, KeyTag, NULL, NULL, };

stream instr;					/* input snapshot data	    */
real refscale;					/* size of ref. cube	    */
int vectoroff = -1;				/* offset for vector data   */
int scalaroff = -1;				/* offset for scalar data   */
bool dopcolor;					/* color using LOS vels.    */
int maxfast;
int pcolor, bcolor;
int wscreen, hscreen;

bodyptr btab = NULL;				/* array of bodies to view  */
int nbody;					/* no. of bodies in array   */
real tnow;					/* time value for body data */

#define XYANGLES   0				/* x and y viewing angles   */
#define PERSPECT   1				/* perspective parameters   */
#define ZANGLE     2				/* z viewing angle	    */
#define TRANSLATE  3				/* x and y displacements    */
#define COLORMAP   4				/* scalar to color mapping  */
#define VSCALE     5				/* vector scale factor	    */
#define NVALUES    6				/* mouse-controlled values  */

int butbind[] = { XYANGLES, PERSPECT, TRANSLATE };

string butlabels[NVALUES] = {
    "xy rotation",      "perspective",      "z rotation",
    "translation",      "color mapping",    "vector scale"
};

int actval = -1;				/* value modified by mouse  */
int xlast, ylast;				/* previous mouse position  */

real mouseval[NVALUES][2];

float thetax, thetay, thetaz;
float xoff, yoff, dview, fview;
float vscale;
float cmidpt, crange;
vector znorm;
float projmat[4][4], viewmat[4][4];

#define DEG2RAD(x)  ((x) * PI / 180.0)

#define RVAL(x)  (((x) & 0xff) / 255.0)
#define GVAL(x)  ((((x) >> 8) & 0xff) / 255.0)
#define BVAL(x)  ((((x) >> 16) & 0xff) / 255.0)

int main(int argc, string argv[])
{
    glutInit(&argc, argv);
    initparam(argv, defv);
    instr = stropen(getparam("in"), "r");
    refscale = getdparam("refscale");
    if (! strnull(getparam("colordata"))) {	/* color data wanted?	    */
        if (getbparam("dopcolor"))
	    error("%s: colordata precludes dopcolor\n", getargv0());
	if (! scanopt(PhiTag "," SmoothTag "," RhoTag "," EntFuncTag ","
		      UinternTag "," UdotIntTag "," UdotRadTag ","
		      UdotVisTag "," TauTag "," BirthTag "," DeathTag ","
		      AuxTag, getparam("colordata")))
	    error("%s: %s unknown\n", getargv0(), getparam("colordata"));
	bodytags[1] = getparam("colordata");	/* replace key w/ field...  */
	butbind[2] = COLORMAP;
    } else if (getbparam("dopcolor")) {
        dopcolor = TRUE;
	bodytags[1] = VelTag;			/* replace key w/ velocity  */
	butbind[2] = COLORMAP;
    }
    if (! strnull(getparam("vectordata"))) {
        if (! scanopt(VelTag "," AccTag "," AuxVecTag,
		      getparam("vectordata")))
	    error("%s: %s unknown\n", getargv0(), getparam("vectordata"));
	if (! (streq(getparam("vectordata"), VelTag) && dopcolor))
	    bodytags[2] = getparam("vectordata");
	butbind[2] = VSCALE;
    }
    maxfast = getiparam("maxfast");
    if (sscanf(getparam("defcolors"), "%x,%x", &pcolor, &bcolor) != 2)
	error("%s: can't scan defcolor parameter\n", getargv0());
    if (sscanf(getparam("viewsize"), "%ix%i", &wscreen, &hscreen) != 2)
	error("%s: can't scan viewsize parameter\n", getargv0());
    layout_body(bodytags, Precision, NDIM);
    if (! strnull(getparam("colordata"))) {
        scalaroff =
	  streq(bodytags[1], PhiTag)     ? PhiField.offset :
	  streq(bodytags[1], SmoothTag)  ? SmoothField.offset :
	  streq(bodytags[1], RhoTag)     ? RhoField.offset :
	  streq(bodytags[1], EntFuncTag) ? EntFuncField.offset :
	  streq(bodytags[1], UinternTag) ? UinternField.offset :
	  streq(bodytags[1], UdotIntTag) ? UdotIntField.offset :
	  streq(bodytags[1], UdotRadTag) ? UdotRadField.offset :
	  streq(bodytags[1], UdotVisTag) ? UdotVisField.offset :
	  streq(bodytags[1], TauTag)     ? TauField.offset :
	  streq(bodytags[1], BirthTag)   ? BirthField.offset :
	  streq(bodytags[1], DeathTag)   ? DeathField.offset :
	  streq(bodytags[1], AuxTag)     ? AuxField.offset : -1;
	assert(scalaroff != -1);
    }
    if (! strnull(getparam("vectordata"))) {
        vectoroff =
	  streq(getparam("vectordata"), VelTag)    ? VelField.offset :
	  streq(getparam("vectordata"), AccTag)    ? AccField.offset :
	  streq(getparam("vectordata"), AuxVecTag) ? AuxVecField.offset : -1;
	assert(vectoroff != -1);
    }
    if (! getdata())
        error("%s: no data in input file\n", getargv0());
    initgraphics(argv[0]);
    glutMainLoop();
    return (0);
}

bool getdata(void)
{
    static bool firstcall = TRUE;
    string intags[MaxBodyFields];
    bodyptr bp;

    get_history(instr);
    if (! get_snap_t(instr, &btab, &nbody, &tnow, intags, FALSE,
		     getparam("times")))
	return (FALSE);
    if (firstcall && ! set_member(intags, bodytags[0]))
        error("%s: required %s data missing\n", getargv0(), bodytags[0]);
    if (firstcall && ! set_member(intags, bodytags[1]))
        if (streq(bodytags[1], KeyTag)) {
            eprintf("[%s: using default point color]\n", getargv0());
            for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
	        Key(bp) = pcolor;
	} else
	    error("%s: required %s data missing\n", getargv0(), bodytags[1]);
    if (firstcall && bodytags[2] != NULL && ! set_member(intags, bodytags[2]))
        error("%s: required %s data missing\n", getargv0(), bodytags[2]);
    firstcall = FALSE;
    return (TRUE);
}

void initgraphics(string name)
{
    int i;

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(wscreen, hscreen);
    glutCreateWindow(name);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glutDisplayFunc(display); 
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    for (i = 0; i < NVALUES; i++)
        mouseval[i][0] = mouseval[i][1] = 0.0;
}

void reshape(int wscr, int hscr)
{
    glViewport(0, 0, (GLsizei) wscr, (GLsizei) hscr);
    wscreen = wscr;
    hscreen = hscr;
}

void display(void)
{
    setview();
    render();
    annotate();
    glFlush();
    glutSwapBuffers();
}

void setview(void)
{
    thetax = -360.0 * mouseval[XYANGLES][1];
    thetay = -360.0 * mouseval[XYANGLES][0];
    thetaz = -360.0 * mouseval[ZANGLE][0];
    xoff = mouseval[TRANSLATE][0];
    yoff = -mouseval[TRANSLATE][1];
    dview = 4.0 * pow(10.0, mouseval[PERSPECT][1]);
    fview = MIN(90.0, 40.0 * pow(10.0, mouseval[PERSPECT][0]));
    vscale = 0.01 * pow(10.0, mouseval[VSCALE][0]);
    cmidpt = mouseval[COLORMAP][0];
    crange = pow(10.0, mouseval[COLORMAP][1]);
    znorm[0] = - rsin(DEG2RAD(thetay));
    znorm[1] =   rcos(DEG2RAD(thetay)) * rsin(DEG2RAD(thetax));
    znorm[2] = - rcos(DEG2RAD(thetay)) * rcos(DEG2RAD(thetax));
}

void render(void)
{
    int step, oldc, newc;
    bodyptr bp;
    float cval;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective((GLdouble) fview, (GLdouble) ((double) wscreen) / hscreen,
		   (GLdouble) 0.01 * dview, (GLdouble) 100.0 * dview);
    glGetFloatv(GL_PROJECTION_MATRIX, &projmat[0][0]);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT);
    glPushMatrix();
    glTranslatef(xoff, yoff, -dview);
    glRotatef(-thetaz, 0.0, 0.0, 1.0);		/* same map as snaprotate   */
    glRotatef(-thetay, 0.0, 1.0, 0.0);
    glRotatef(-thetax, 1.0, 0.0, 0.0);
    glGetFloatv(GL_MODELVIEW_MATRIX, &viewmat[0][0]);
    glBegin(vectoroff == -1 ? GL_POINTS : GL_LINES);
    step = (actval == -1 ? 1 : rceil((float) nbody / (float) maxfast));
    oldc = -1;
    for (bp = btab; bp < NthBody(btab, nbody); bp = NthBody(bp, step)) {
        if (scalaroff == -1 && ! dopcolor)
            newc = Key(bp);
        else {
	    if (! dopcolor)
	        cval = (SelectReal(bp, scalaroff) - cmidpt) / crange;
	    else
	        cval = (dotvp(Vel(bp), znorm) - cmidpt) / crange;
	    newc = (cval >  1.0 ? 0x0000ff : cval >  0.6 ? 0x006fdf :
		    cval >  0.2 ? 0x00cf7f : cval > -0.2 ? 0x00ff00 :
		    cval > -0.6 ? 0x7fcf00 : cval > -1.0 ? 0xbf8f00 :
		    0xff4f00);
	}
        if (oldc != newc)
	    glColor3f(RVAL(newc), GVAL(newc), BVAL(newc));
	oldc = newc;
        glVertex3f(Pos(bp)[0], Pos(bp)[1], Pos(bp)[2]);
	if (vectoroff != -1)
	    glVertex3f(Pos(bp)[0] + vscale * SelectVect(bp, vectoroff)[0],
		       Pos(bp)[1] + vscale * SelectVect(bp, vectoroff)[1],
		       Pos(bp)[2] + vscale * SelectVect(bp, vectoroff)[2]);
    }
    glEnd();
    if (actval != -1) {
        glColor3f(RVAL(bcolor), GVAL(bcolor), BVAL(bcolor));
	glutWireCube(refscale);
    }
    glPopMatrix();
}

void annotate(void)
{
    char buf[64];
    int ypos, xpos;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho((GLdouble) 0, (GLdouble) wscreen,
	    (GLdouble) hscreen, (GLdouble) 0,
	    (GLdouble) -1.0, (GLdouble) 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glColor3f(RVAL(0xffffff), GVAL(0xffffff), BVAL(0xffffff));
    if (actval == -1) {
	if (! strnull(getparam("ident")))
            showtext(getparam("ident"), 10, 20);
        sprintf(buf, "%.2f", tnow);
	showtext(buf, wscreen - (8 * strlen(buf) + 10), 20);
    } else {
        ypos = hscreen - 10;
	if (vectoroff != -1) {
	  sprintf(buf, "vector scale: %.2f", vscale);
	  showtext(buf, 10, ypos);
	  ypos -= 20;
	}
	if (scalaroff != -1 || dopcolor) {
	    xpos = 10;
	    sprintf(buf, "mapping: ");
	    showtext(buf, xpos, ypos);
	    xpos += 8 * strlen(buf);
	    sprintf(buf, "%.2f", cmidpt + crange);
	    glColor3f(RVAL(0x0000ff), GVAL(0x0000ff), BVAL(0x0000ff));
	    showtext(buf, xpos, ypos);
	    xpos += 8 * strlen(buf);
	    glColor3f(RVAL(0xffffff), GVAL(0xffffff), BVAL(0xffffff));
	    showtext(" to ", xpos, ypos);
	    xpos += 32;
	    sprintf(buf, "%.2f", cmidpt - crange);
	    glColor3f(RVAL(0xff4f00), GVAL(0xff4f00), BVAL(0xff4f00));
	    showtext(buf, xpos, ypos);
	    xpos += 8 * strlen(buf);
	    ypos -= 20;
	    glColor3f(RVAL(0xffffff), GVAL(0xffffff), BVAL(0xffffff));
	}
	sprintf(buf, "angles: %.2f, %.2f, %.2f", thetax, thetay, thetaz);
	showtext(buf, 10, ypos);
	ypos -= 20;
        sprintf(buf, "view: %.2f, %.2f, %.2f; %.2f",
		xoff, yoff, dview, fview);
	showtext(buf, 10, ypos);
	ypos -= 20;
    }	
    glPopMatrix();
}

void showtext(string str, int x, int y)
{
    int i;

    glRasterPos2f((GLfloat) x, (GLfloat) y);
    for (i = 0; str[i] != (int) NULL; i++)
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, str[i]);
}

void mouse(int but, int state, int x, int y)
{
    if (state == GLUT_DOWN && actval == -1) {
	switch (but) {
	  case GLUT_LEFT_BUTTON:
	    actval = butbind[0];
	    break;
	  case GLUT_MIDDLE_BUTTON:
	    actval = butbind[1];
	    break;
	  case GLUT_RIGHT_BUTTON:
	    actval = butbind[2];
	    break;
	}
	xlast = x;
	ylast = y;
    } else if (state == GLUT_UP && actval != -1) {
        actval = -1;
    } else
        eprintf("\007[%s: mouse button ignored]\n", getargv0());
    glutPostRedisplay();
}

void motion(int x, int y)
{
    if (actval != -1) {
        mouseval[actval][0] += (x - xlast) / ((float) wscreen);
	mouseval[actval][1] += (y - ylast) / ((float) hscreen);
	xlast = x;
	ylast = y;
	glutPostRedisplay();
    }

}

void keyboard(unsigned char key, int x, int y)
{
    int i;

    switch (key) {
      case 27:
	exit(0);
        break;
      case 32:
        if (! getdata())
	    eprintf("[%s: at end of file]\n", getargv0());
        glutPostRedisplay();
	break;
      case 'a':
        eprintf("[%s: angles thetax=%.2f thetay=%.2f thetaz=%.2f]\n",
		getargv0(), thetax, thetay, thetaz);
	break;
      case 'b':
        eprintf("[%s: mouse buttons adjust %s, %s, %s]\n", getargv0(),
		butlabels[butbind[0]],
		butlabels[butbind[1]],
		butlabels[butbind[2]]);
	break;
      case 'm':
        butbind[2] = COLORMAP;
	break;
      case 'p':
        butbind[1] = PERSPECT;
	break;
      case 'r':
	mouseval[XYANGLES][0] = 0;
	mouseval[XYANGLES][1] = 0;
	mouseval[ZANGLE][0] = 0;
	glutPostRedisplay();
	break;
      case 'R':
        for (i = 0; i < NVALUES; i++)
            mouseval[i][0] = mouseval[i][1] = 0.0;
	glutPostRedisplay();
	break;
      case 's':
        butbind[2] = VSCALE;
	break;
      case 't':
        butbind[2] = TRANSLATE;
	break;
      case 'v':
	writeview();
	break;
      case 'z':
        butbind[1] = ZANGLE;
	break;
      default:
	eprintf("[%s: keystroke %d ignored]\n", getargv0(), (int) key);
    }
}

void writeview(void)
{
    static stream viewstr = NULL;
    real aspect, tempmat[4][4];
    int i, j, k;

    if (viewstr == NULL) {
        viewstr = stropen(getparam("viewfile"), "w");
	put_history(viewstr);
    }
    aspect = ((float) wscreen) / ((float) hscreen);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++) {
	    tempmat[i][j] = 0.0;
	    for (k = 0; k < 4; k++)
	        tempmat[i][j] += viewmat[i][k] * projmat[k][j];
        }
    put_set(viewstr, "SnapView");
    put_data(viewstr, "Time", RealType, &tnow, 0);
    put_data(viewstr, "Aspect", RealType, &aspect, 0);
    put_data(viewstr, "ViewMatrix", FloatType, tempmat, 4, 4, 0);
    put_data(viewstr, "MouseValues", FloatType, mouseval, NVALUES, 2, 0);
    put_tes(viewstr, "SnapView");
    fflush(viewstr);
}
