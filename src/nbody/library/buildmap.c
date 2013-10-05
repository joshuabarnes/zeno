/*
 * BUILDMAP.C: generate and compile client transformation program.
 */

#include "stdinc.h"
#include "datatypes.h"
#include "getparam.h"
#include "mapdefs.h"
#include <ctype.h>

#ifdef TESTBED

string defv[] = {		";Invoke buildmap routine",
    "prog=sm",			";Name of program to compile",
    "names=XPlot,YPlot,Color",	";Names of output variables",
    "exprs=x;y;1+i%1",		";Expressions for output vars",
    "types=" FloatType "," FloatType "," IntType,
				";Data types of output vars",
    "prec=SINGLEPREC",		";Specify precision option",
    "ndim=3",			";Number of space dimensions",
    "VERSION=1.3",		";Josh Barnes  10 June 2012",
    NULL,
};

void buildmap(string, string *, string *, string *, string, int);

int main(int argc, string argv[])
{
    string *names, *exprs, *types;

    initparam(argv, defv);
    names = burststring(getparam("names"), ",");
    exprs = burststring(getparam("exprs"), ";");
    if (! strnull(getparam("types")))
	types = burststring(getparam("types"), ",");
    else
	types = NULL;
    buildmap(getparam("prog"), names, exprs, types,
	     getparam("prec"), getiparam("ndim"));
    return (0);
}

#endif

void buildmap(string prog, string *names, string *exprs, string *types,
	      string prec, int ndim)
{
    char src[80], cmd[512];
    stream sstr;
    int i;
    string *tp, *np, *ep, tn;

    sprintf(src, "%s.c", prog);
    sstr = stropen(src, "w");
    fprintf(sstr, "#include \"stdinc.h\"\n");
    fprintf(sstr, "#include \"mathfns.h\"\n");
    fprintf(sstr, "#include \"vectdefs.h\"\n");
    if (getenv("ZENO_SAFE_SELECT") != NULL)
	fprintf(sstr, "#define SafeSelect TRUE\n");
    fprintf(sstr, "#include \"phatbody.h\"\n\n");
    for (i = 0; mapdefs[i][0] != NULL; i++)
	fprintf(sstr, "#define %-8s %s(_p)\n", mapdefs[i][0], mapdefs[i][1]);
    fprintf(sstr, "\n");
    for (i = 0; expdefs[i][0] != NULL; i++)
	fprintf(sstr, "#define %-8s %s\n", expdefs[i][0], expdefs[i][1]);
    fprintf(sstr, "\n");
    if (types != NULL) {
	for (tp = types, np = names; *tp != NULL; tp++, np++) {
	    if (*np == NULL)
		error("buildmap in %s: more types than names\n", getargv0());
	    tn = type_name(*tp);
	    fprintf(sstr, "#define %s(b)  Select%c%s"
		    "(b, phatbody[NewBodyFields+%d].offset)\n",
		    *np, toupper(tn[0]), tn + 1, (int) (tp - types));
	}
	fprintf(sstr, "\n");
	fprintf(sstr, "void extendbody(void)\n");
	fprintf(sstr, "{\n");
	for (tp = types, np = names; *tp != NULL; tp++, np++)
	    fprintf(sstr, "    new_field(&phatbody[NewBodyFields+%d], "
		    "\"%s\", \"%s\");\n", (int) (tp - types), *tp, *np);
	fprintf(sstr, "    new_field(&phatbody[NewBodyFields+%d],"
		" NULL, NULL);\n", (int) (tp - types));
	fprintf(sstr, "}\n\n");
    } else
	fprintf(sstr, "void extendbody(void)\n{ }\n\n");
    fprintf(sstr,
	    "void computemap(bodyptr _q, bodyptr _p, real t, int i, int n)\n");
    fprintf(sstr, "{\n");
    for (ep = exprs, np = names; *ep != NULL; ep++, np++) {
	if (*np == NULL)
	    error("buildmap in %s: more exprs than names\n", getargv0());
	fprintf(sstr, "    %s(_q) = (%s);\n", *np, *ep);
    }
    fprintf(sstr, "}\n");
    fclose(sstr);
    sprintf(cmd, "%s %s %s -D%s -DNDIM=%d -o %s %s.c %s/lib/snapmap_%c%d.o "
            "-lNBody -lClib -lgsl -lgslcblas -lm", getenv("ZCC"), 
	    getenv("ZCCFLAGS"), getenv("ZLDFLAGS"), prec, ndim, prog, prog,
	    getenv("ZENOPATH"), tolower(prec[0]), ndim);
    eprintf("[%s: compiling snapmap]\n", getargv0());
    if (system(cmd) != 0)
        error("%s: command \"%s\" failed\n", getargv0(), cmd);
    if (unlink(src) != 0)
        error("%s: can't unlink %s\n", getargv0(), src);
}
