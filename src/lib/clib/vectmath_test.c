/*
 * vectmath_test.c: test operators in vectmath.h.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"

local void prints(string m, double s1, double s2);
local void printv(string m, vector v);
local void printm(string m, matrix p);

vector ex, ey, ez, v, u, w;
matrix I, p, q, r;
real s;

int main(int argc, string argv[])
{
    UNITV(ex, 0);		printv("ex", ex);
    UNITV(ey, 1); 		printv("ey", ey);
    UNITV(ez, 2); 		printv("ez", ez);
    MULVS(v, ex, 3.141592); 	printv("v = ex * pi", v);
    ADDV(u, v, ez); 		printv("u = v + ez", u);
    SUBV(w, u, ey); 		printv("w = u - ey", w);
    DOTVP(s, w, ex); 		prints("s = w . ex", s, dotvp(w, ex));
    DOTVP(s, w, ey); 		prints("s = w . ey", s, dotvp(w, ey));
    DOTVP(s, w, ez); 		prints("s = w . ez", s, dotvp(w, ez));
    ABSV(s, w); 		prints("|w|", s, absv(w));
    DISTV(s, w, ex); 		prints("|w - ex|", s, distv(w, ex));
    DISTSQV(s, w, v);		prints("|w - v|^2", s, rsqr(distv(w, v)));
    CROSSVP(v, ex, ey); 	printv("v = ex x ey", v);
    CROSSVP(u, w, ez); 		printv("u = w x ez", u);
    SETMI(I); 			printm("I", I);
    OUTVP(p, v, ex); 		printm("p = v @ ex", p);
    OUTVP(q, w, u); 		printm("q = w @ u", q);
    ADDM(r, p, q);		printm("r = p + q", r);
    TRANM(p, q);		printm("p = transpose(q)", p);
    MULM(q, p, r);		printm("q = p * r", q);
    MULMV(v, q, w);		printv("v = q . w", v);
    MULMV(u, r, w);		printv("u = r . w", u);
    MULMV(v, p, u);		printv("v = p . u", v);
    TRACEM(s, q);		prints("s = Tr(q)", s, tracem(q));
    return (0);
}

local void prints(string m, double s1, double s2)
{
    printf("\n%-24s: %12.6f%12.6f\n\n", m, s1, s2);
}

local void printv(string m, vector v)
{
    printf("\n%-24s: %12.6f%12.6f%12.6f\n\n", m, v[0], v[1], v[2]);
}

local void printm(string m, matrix p)
{
    printf("\n%-24s  %12.6f%12.6f%12.6f\n"  , "", p[0][0], p[0][1], p[0][2]);
    printf(  "%-24s: %12.6f%12.6f%12.6f\n"  ,  m, p[1][0], p[1][1], p[1][2]);
    printf(  "%-24s  %12.6f%12.6f%12.6f\n\n", "", p[2][0], p[2][1], p[2][2]);
}
