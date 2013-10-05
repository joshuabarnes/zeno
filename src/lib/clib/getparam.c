/*
 * getparam.c: command-line processing functions.
 */

#include "stdinc.h"
#include "getparam.h"
#include <string.h>

//  Local routines and definitions.

//  ______________________________________________
//  param: structure encoding parameter and value.

typedef struct {
  string name;					// name of parameter
  string value;					// value of parameter
  string comment;				// documentation string
  int flags;					// for various options
} param;

local int countdefaults(string *);		// number of defaults
local void setprogram(param *, string);		// set 0th parameter
local void copydefaults(param *, string *);	// set default parameters
local void checkhelp(param *, string);		// help processing
local void printparam(string, string, bool);	// print param and comment
local void setarguments(param *, string *);	// set command parameters
local void checkusage(param *);			// usage processing
local void savehistory(param *);		// history interface
local param *findparam(string, param *);	// look up parameter
local string parname(string);			// extract param name
local string parvalue(string);			// extract param value

local param *paramvec = NULL;			// vector of parameters

local string progname;				// program name, for errors

#define BUFLEN  1025				// storage for comments

//  _________________________________________________________________
//  initparam: initalize parameter lists and handle special requests.

void initparam(string *argv, string *defv)
{
  int nparam;
  param *pvec;

  progname = argv[0];				// initialize program name
  nparam = 1 + countdefaults(defv);		// include argv0 in count
  pvec = (param *) allocate(sizeof(param) * (nparam + 1));
  setprogram(pvec, argv[0]);			// install 0th argument
  copydefaults(pvec, defv);			// set up default values
  if (argv[1] != NULL)				// if first arg is given...
    checkhelp(pvec, argv[1]);			// see if help is requested
  setarguments(pvec, argv);			// args override defaults
  checkusage(pvec);				// complain if args missing
  paramvec = pvec;				// install parameter vector
  savehistory(pvec);				// save activation record
}

//  __________________________________________________
//  countdefaults: count number of default parameters.

local int countdefaults(string *defv)
{
  int ndefault;
  string *dp;

  ndefault = 0;
  for (dp = defv; *dp != NULL; dp++)		// loop over all defaults
    if (**dp != ';')				// if not a comment
      ndefault++;				// then count one more
  return (ndefault);
}

//  _____________________________________________________________
//  setprogram: initialize the program name as parameter "argv0".

local void setprogram(param *pvec, string argv0)
{
  pvec->name = "argv0";				// install 0th parameter
  pvec->value = argv0;				// set name from argv[0]
  pvec->comment = NULL;				// no comment for now
  pvec->flags = ARGPARAM;			// so user can't reset it
}

//  ______________________________________________________
//  copydefaults: install default parameters and comments.

local void copydefaults(param *pvec, string *defv)
{
  param *pp;
  string *dp, name, value, com;

  pp = pvec;					// start with 0th param
  for (dp = defv; *dp != NULL; dp++)		// loop over the defaults
    if (**dp != ';') {				// if not a comment...
      pp++;					// go on to new param
      name = parname(*dp);			// extract param name
      value = parvalue(*dp);			// and param value
      if (name == NULL || value == NULL)	// any problems found?
	error("%s.initparam: bad parameter %s\n", progname, *dp);
      pp->name = strdup(name);			// assign param name
      pp->value = strdup(value);		// and parameter value
      pp->comment = NULL;			// clear comment field
      pp->flags = DEFPARAM;			// source is default
      if (streq(pp->value, "???"))		// a required param?
	pp->flags |= REQPARAM;			// set required flag
      if (**dp == '<')				// an input parameter?
	pp->flags |= INPARAM;			// set input flag
      else if (**dp == '>')			// an output parameter?
	pp->flags |= OUTPARAM;			// set output flag
    } else if (pp->comment == NULL)	{	// a new comment...
      pp->comment = strdup(*dp + 1);		// assign comment field
    } else {					// more text for comment...
      com = (string) allocate(strlen(pp->comment) + strlen(*dp+1) + 2);
						// get space for comment
      strcpy(com, pp->comment);			// copy existing text
      strcat(com, "\n");			// indicate end of line
      strcat(com, *dp + 1);			// store text of comment
      free(pp->comment);			// free up old space
      pp->comment = com;			// and store new comment
    }
  pp++;						// past last real param
  pp->name = NULL;				// end list of parameters
}

//  __________________________________________________________
//  checkhelp: if requested, print out help mesaages and exit.

local void checkhelp(param *pvec, string arg1)
{
  param *pp;
  char buf[BUFLEN];

  if (streq(arg1, "-clue") || streq(arg1, "--clue")) {
						// print brief help message
    printf("%s", pvec->value);
    for (pp = pvec+1; pp->name != NULL; pp++)
      printf(" %s=%s", pp->name, pp->value);
    printf("\n");
    exit(0);
  }
  if (streq(arg1, "-help") || streq(arg1, "--help")) {
						// print full help message
    printparam(pvec->value, pvec->comment, FALSE);
    for (pp = pvec+1; pp->name != NULL; pp++) {
      sprintf(buf, "  %s=%s", pp->name, pp->value);
      printparam(buf, pp->comment, FALSE);
    }
    exit(0);
  }
  if (streq(arg1, "-explain") || streq(arg1, "--explain")) {
						// print long help message
    printf("\n");
    printparam(pvec->value, pvec->comment, TRUE);
    for (pp = pvec+1; pp->name != NULL; pp++) {
      sprintf(buf, "  %s=%s", pp->name, pp->value);
      printparam(buf, pp->comment, TRUE);
    }
    exit(0);
  }
}

//  _______________________________________________
//  printparam: print parameter and comment string.

local void printparam(string item, string comment, bool explain)
{
  char buf[BUFLEN], *bp1, *bp2;

  if (comment == NULL)
    printf("%s\n", item);
  else {
    strcpy(buf, comment);
    bp1 = strchr(buf, '\n');
    if (bp1 != NULL)
      *bp1++ = (char) NULL;
    printf((strlen(item) < 32 ? "%-32s  %s\n" : "%s\n\t\t\t\t  %s\n"),
	   item, buf);
    if (explain) {
      while (bp1 != NULL) {
	bp2 = strchr(bp1, '\n');
	if (bp2 != NULL)
	  *bp2++ = (char) NULL;
	printf("\t\t\t\t  %s\n", bp1);
	bp1 = bp2;
      }
      printf("\n");
    }
  }
}

//  __________________________________________________________
//  setarguments: replace defaults with commandline arguments.

local void setarguments(param *pvec, string *argv)
{
  bool scanpos;
  param *pp;
  string *ap, name;

  scanpos = TRUE;				// start scan by position
  pp = pvec;					// start with 0th param
  for (ap = argv + 1; *ap != NULL; ap++) {	// loop over the arguments
    name = parname(*ap);			// get param name, if any
    scanpos = scanpos && (name == NULL);	// see how to match args
    if (scanpos) {				// match by position?
      pp++;					// move to next param
      if (pp->name == NULL)			// past last parameter?
	error("%s.initparam: too many arguments\n", progname);
      pp->value = strdup(*ap);			// set new param value
    } else {					// matching by name?
      if (name == NULL)				// but no name given?
	error("%s.initparam: arg %s must be named\n", progname, *ap);
      pp = findparam(name, pvec);		// look for named param
      if (pp == NULL)				// an unknown param?
	error("%s.initparam: parameter %s unknown\n", progname, name);
      if (pp->flags & ARGPARAM)			// already on arg list?
	error("%s.initparam: parameter %s duplicated\n", progname, name);
      pp->value = strdup(parvalue(*ap));	// set new param value
    }
    pp->flags = (pp->flags & ~DEFPARAM) | ARGPARAM;
						// switch source flag
  }
}

//  ___________________________________________
//  checkusage: print out short message on use.

local void checkusage(param *pvec)
{
  bool needarg;
  int len;
  param *pp;
  string args;

  needarg = FALSE;				// see if any args left out
  for (pp = pvec+1; pp->name != NULL; pp++)	// scan list of parameters
    needarg = needarg | ((pp->flags & REQPARAM) && (pp->flags & DEFPARAM));
						// req. param but defaulted?
  if (needarg) {
    len = 0;
    for (pp = pvec+1; pp->name != NULL; pp++)	// estimate space needed
      len += 5 + strlen(pp->name);
    args = (string) allocate(len + 1);
    for (pp = pvec+1; pp->name != NULL; pp++)
      if ((pp->flags & REQPARAM)) {		// is this arg required?
	strcat(args, pp->name);			// add the argument name
	strcat(args, "=??? ");			// and its default value
      }
    error("Usage: %s %s\n", progname, args);
  }
}

//  ____________________________________________
//  savehistory: interface to history mechanism.

local void savehistory(param *pvec)
{
  param *pp;
  int hlen;
  char *hist;

  pp = findparam("HISTORY", pvec);		// look for a history entry
  if (pp != NULL && ! strnull(pp->value)) 	// is history data present?
    add_history(pp->value);			// just send to history pkg
  else {					// format history data
    hlen = strlen(pvec->value);			// count length of argv0
    for (pp = pvec+1; pp->name != NULL; pp++)
						// sum length of others
      if ((pp->flags & ARGPARAM) || streq(pp->name, "VERSION"))
						// actual arg or version?
	hlen += 2 + strlen(pp->name) + strlen(pp->value);
						// add length required
    hist = (char *) allocate(hlen+1);		// get memory for history
    strcpy(hist, pvec->value);			// store program name
    for (pp = pvec+1; pp->name != NULL; pp++)
						// loop over other params
      if ((pp->flags & ARGPARAM) || streq(pp->name, "VERSION")) {
						// actual arg or version?
	strcat(hist, " ");			// insert a space
	strcat(hist, pp->name);			// and the argument name
	strcat(hist, "=");			// and a equals sign
	strcat(hist, pp->value);		// and the argument value
      }
    add_history(hist);				// invoke history package
  }
}

//  ____________________________________
//  getparam: return value of parameter.

string getparam(string name)
{
  param *par;

  if (paramvec == NULL)	{			// called during init?
    if (streq(name, "argv0") && progname != NULL)
      return (progname);			// handle "argv0" specially
    error("getparam: called before initparam\n");
						// program name unknown
  }
  par = findparam(name, paramvec);		// look up value by name
  if (par == NULL)
    error("%s.getparam: parameter %s unknown\n", progname, name);
  return (par->value);
}

//  _______________________________________
//  getiparam: get integer-value parameter.

int getiparam(string name)
{
  string str = getparam(name), end;
  int val;

  val = strtol(str, &end, 0);
  if (end == str)
    error("%s.getiparam: %s=\"%s\" not int\n", progname, name, str);
  if (*end == 'k' || *end == 'K')
    return (1024 * val);
  if (*end == 'm' || *end == 'M')
    return (1024 * 1024 * val);
  return (val);
}

//  ______________________________________________
//  getdparam: get floating-point value parameter.

double getdparam(string name)
{
  string str = getparam(name), end;
  double val1, val2;

  if (sscanf(str, "%lf/%lf", &val1, &val2) == 2)
    return (val1 / val2);
  val1 = strtod(str, &end);
  if (end == str)
    error("%s.getdparam: %s=\"%s\" not double\n", progname, name, str);
  return (val1);
}

//  _________________________________
//  getbparam: get boolean parameter.

bool getbparam(string name)
{
  string str = getparam(name);

  if (strchr("tTyY1", *str) != NULL)		// is value true?
    return (TRUE);
  if (strchr("fFnN0", *str) != NULL)		// is value false?
    return (FALSE);
  error("%s.getbparam: %s=\"%s\" not bool\n", progname, name, str);
  return (FALSE);				// keep compiler happy...
}

//  __________________________________________________________________
//  getparamstat: return parameter flags, or zero if no such parameter
//  (note that all defined parameters have at least one bit set).

int getparamstat(string name)
{
  param *par;

  par = findparam(name, paramvec);
  return (par != NULL ? par->flags : 0);
}

//  __________________________________________________________
//  findparam: look for named parameter in list of parameters.

local param *findparam(string name, param *pvec)
{
  param *pp;

  for (pp = pvec; pp->name != NULL; pp++)
    if (streq(name, pp->name))
      return (pp);
  return (NULL);
}

//  _____________________________________________
//  parname: extract name from name=value string.
//  WARNING:  returns ptr to static storage.

local string parname(string arg)
{
  char *ap, *ep;
  static char namebuf[64];

  ap = (char *) arg;
  if (*ap == '<' || *ap == '>')
    ap++;
  strncpy(namebuf, ap, 63);
  namebuf[63] = (char) NULL;
  ep = strchr(namebuf, '=');
  if (ep == NULL)				// not of form name=value?
    return (NULL);
  *ep = (char) NULL;
  return (namebuf);
}

//  _______________________________________________
//  parvalue: extract value from name=value string.

local string parvalue(string arg)
{
  char *ep;

  ep = strchr(arg, '=');
  if (ep == NULL)
    return (NULL);
  return (ep + 1);
}

#ifdef TESTBED

string defv[] = {		";Program to test getparam.",
				";This program illustrates features of",
				";Zeno's command-line argument passing.",
  "<input=???",			";This parameter is required.",
				";An error occurs if no value is given.",
  ">output=",			";This parameter is not required.",
				";Its default value is null.",
  "answer=42",			";This parameter needs no comment.",
				";Cool froods KNOW the answer is 42!",
				";Decimal, octal, and hex are accepted.",
				";A trailing k or m multiplies the value",
				";by 1024 or 1048576, respectively.",
  "value=10.0",			";This is a floating-point parameter.",
				";Scientific notation is accepted;",
				";a ratio can be specified as value/value.",
  "flag=false",			";This is a boolean parameter.",
				";The first letter determines the value:",
				";\"tTyY1\" -> TRUE, \"fFnN0\" -> FALSE.",
				";Set to TRUE to see error handling.",
  "foobar=waldo_and_mumble_grumble_for_a_while",
				";This parameter is very long.",
				";So the comment starts on the next line.",
  "VERSION=2.4",		";Josh Barnes  8 May 2012",
				";By convention, the current version",
				";appears as the last parameter.",
  NULL,
};

int main(int argc, string *argv)
{
  initparam(argv, defv);
  printf("program %s:\n", getprog());
  printf("  input = \"%s\" [%o]\n", getparam("input"),
	 getparamstat("input"));
  printf("  output = \"%s\" [%o]\n", getparam("output"),
	 getparamstat("output"));
  printf("  answer = %d [%o]\n", getiparam("answer"),
	 getparamstat("answer"));
  printf("  value = %g [%o]\n", getdparam("value"),
	 getparamstat("value"));
  printf("  flag = %s [%o]\n", getbparam("flag") ? "TRUE" : "FALSE",
	 getparamstat("flag"));
  printf("  foobar = \"%s\" [%o]\n", getparam("foobar"),
	 getparamstat("foobar"));
  printf("  VERSION = \"%s\" [%o]\n", getversion(),
	 getparamstat("VERSION"));
  if (getbparam("flag")) {
    printf("getparamstat(\"junk\") = %o\n", getparamstat("junk"));
    printf("calling getparam(\"junk\")\n");
    (void) getparam("junk");
  }
  return (0);
}

#endif
