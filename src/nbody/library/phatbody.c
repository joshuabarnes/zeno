/*
 * phatbody.c: Procedures to define and extend phat body structures.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "phatbody.h"
#include "datatypes.h"
#include <string.h>

//  Temporary tokens for real, vector, and array types, replaced by routine
//  fix_types() once client's precision and dimension number are determined.
//  ________________________________________________________________________

#define RealToke "RT"			// placeholder for RealType
#define VectToke "VT"			// placeholder for VectType
#define RealArr  "RA"			// placeholder for RealArray
#define IntArr   "IA"			// placeholder for IntArray

//  phatbody: Defines the components and offsets of a body structure.
//  Most applications use only a subset, while some extend the list.
//  _________________________________________________________________

ps_field phatbody[MaxBodyFields] = {
  { NULL,      BodyTag,     0,          0 },	// = phatbody[ 0]
  { VectToke,  PosTag,      BadOffset,  0 },	// = phatbody[ 1]
  { VectToke,  VelTag,      BadOffset,  0 },	// = phatbody[ 2]
  { RealToke,  MassTag,     BadOffset,  0 },	// = phatbody[ 3]
  { RealToke,  PhiTag,      BadOffset,  0 },	// = phatbody[ 4]
  { VectToke,  AccTag,      BadOffset,  0 },	// = phatbody[ 5]
  { RealToke,  SmoothTag,   BadOffset,  0 },	// = phatbody[ 6]
  { RealToke,  RhoTag,      BadOffset,  0 },	// = phatbody[ 7]
  { RealToke,  EntFuncTag,  BadOffset,  0 },	// = phatbody[ 8]
  { RealToke,  UinternTag,  BadOffset,  0 },	// = phatbody[ 9]
  { RealToke,  UdotIntTag,  BadOffset,  0 },	// = phatbody[10]
  { RealToke,  UdotRadTag,  BadOffset,  0 },	// = phatbody[11]
  { RealToke,  UdotVisTag,  BadOffset,  0 },	// = phatbody[12]
  { RealToke,  TauTag,      BadOffset,  0 },	// = phatbody[13]
  { ByteType,  TypeTag,     BadOffset,  0 },	// = phatbody[14]
  { RealToke,  BirthTag,    BadOffset,  0 },	// = phatbody[15]
  { RealToke,  DeathTag,    BadOffset,  0 },	// = phatbody[16]
  { IntType,   KeyTag,      BadOffset,  0 },	// = phatbody[17]
  { IntArr,    KeyArrTag,   BadOffset,  0 },	// = phatbody[18]
  { RealToke,  AuxTag,      BadOffset,  0 },	// = phatbody[19]
  { VectToke,  AuxVecTag,   BadOffset,  0 },	// = phatbody[20]
  { RealArr,   AuxArrTag,   BadOffset,  0 },	// = phatbody[21]
  { NULL,      NULL,        0,          0 },	
};

local string prec0 = NULL;
local int ndim0 = 0;

local void fix_types(string *tags, string prec, int ndim);
local string list_tags(string *tags);
local string *clean_tags(string *tags);
local bool array_tag(string tag);
local string array_type(string type, int len);

//  NOTE: these routines are called at the start of a process, to
//  initialize the layout of body structures.  They are written for
//  clarity rather than speed, and don't reclaim allocated storage.
//  _______________________________________________________________

//  layout_body: define real and vector types according to client,
//  and add named fields to body structure.
//  ______________________________________________________________

void layout_body(string *tags, string prec, int ndim)
{
  if (getenv("ZENO_PHSTR_DEBUG") != NULL)
    eprintf("[%s.layout_body: tags = %s  prec = %s  ndim = %d]\n",
	    getprog(), list_tags(tags), prec, ndim);
  prec0 = prec;					// save for later reference
  ndim0 = ndim;
  fix_types(tags, prec, ndim);			// set type codes from tags
  layout_struct(phatbody, clean_tags(tags));	// init. offsets and length
}

//  define_body: define real and vector types according to client,
//  and sets body length.
//  ______________________________________________________________

void define_body(int length, string prec, int ndim)
{
  if (getenv("ZENO_PHSTR_DEBUG") != NULL)
    eprintf("[%s.define_body: length = %d  prec = %s  ndim = %d]\n",
	    getprog(), length, prec, ndim);
  prec0 = prec;					// save for later reference
  ndim0 = ndim;
  fix_types(NULL, prec, ndim);			// set non-array type codes
  define_struct(phatbody, BodyTag, length);	// init. structure length
}

//  define_body_offset: fixed-offset interface sets field offset.
//  _____________________________________________________________

void define_body_offset(string tag, int offset)
{
  if (getenv("ZENO_PHSTR_DEBUG") != NULL)
    eprintf("[%s.define_body_offset: tag = %s  offset = %d]\n",
	    getprog(), tag, offset);
  if (array_tag(tag)) {				// defining array field?
    string tags[2] = { tag, NULL };
    if (prec0 == NULL || ndim0 == 0)
      error("%s.define_body_offset: prec or ndim not known\n", getprog());
    fix_types(tags, prec0, ndim0);		// set length of array field
    define_offset(phatbody, clean_tags(tags)[0], offset);
  } else
    define_offset(phatbody, tag, offset);
}

//  fix_types: replace RealToke, VectToke, RealArr, and IntArr with
//  run-time type specifications.
//  _______________________________________________________________

local void fix_types(string *tags, string prec, int ndim)
{
  string rtyp, btyp, *tp;
  ps_field *pbf;
  char fmt[512];
  int len;

  if (getenv("ZENO_PHSTR_DEBUG") != NULL)
    eprintf("[%s.fix_types:  tags = %s  prec = %s  ndim = %d]\n",
	    getprog(), tags != NULL ? list_tags(tags) : "NULL", prec, ndim);
  rtyp = (streq(prec, "SINGLEPREC") || streq(prec, "MIXEDPREC")) ? FloatType :
         streq(prec, "DOUBLEPREC") ? DoubleType :  NULL;
  if (rtyp == NULL)
    error("%s.fix_types: unknown precision %s\n", getargv0(), prec);
  for (pbf = &phatbody[1]; pbf->name != NULL; pbf++) {
    if (streq(pbf->type, RealToke))		// replace real with type
      pbf->type = rtyp;
    if (streq(pbf->type, VectToke))		// make type for vector
      pbf->type = array_type(rtyp, ndim);
    if ((streq(pbf->type, RealArr) || streq(pbf->type, IntArr)) &&
	tags != NULL) {				// determine type for array
      btyp = streq(pbf->type, RealArr) ? rtyp : IntType;
      sprintf(fmt, "%s[%%i]", pbf->name);
      for (tp = tags; *tp != NULL; tp++)	// scan tags for match
	if (sscanf(*tp, fmt, &len) == 1) {
	  pbf->type = array_type(btyp, len);
	  if (getenv("ZENO_PHSTR_DEBUG") != NULL)
	    eprintf("[%s.fix_types: array %s: len = %d  type = %s]\n",
		    getprog(), pbf->name, len, pbf->type);
	}
    }
  }
}

//  list_tags: concat tags into single string for debug messages.
//  _____________________________________________________________

local string list_tags(string *tags)
{
  int len = 1;
  string *tp, res;

  for (tp = tags; *tp != NULL; tp++)
    len += strlen(*tp) + 1;
  res = (string) allocate(len * sizeof(char));
  for (tp = tags; *tp != NULL; tp++) {
    if (tp != tags)
      strcat(res, ",");
    strcat(res, *tp);
  }
  return (res);
}

//  clean_tags: copy tag list without length specifications.
//  ________________________________________________________

local string *clean_tags(string *tags0)
{
  string *tags1 = (string *) copxstr(tags0, sizeof(string)), *tp;

  for (tp = tags1; *tp != NULL; tp++)
    if (array_tag(*tp)) {
      *tp = (string) copxstr(*tp, sizeof(char));
      *index(*tp, '[') = (char) NULL;
    }
  return (tags1);
}

//  array_tag: test if tag includes length specification.
//  _____________________________________________________

local bool array_tag(string tag)
{
  return (index(tag, '[') != NULL);
}

//  array_type: construct compound type for array.
//  ______________________________________________

local string array_type(string type, int len)
{
  string atyp = (string) allocate((len + 1) * sizeof(char));

  for (int j = 0; j < len; j++)
    atyp[j] = type[0];
  return (atyp);
}

#if defined(TESTBED)

#include "vectmath.h"

string defv[] = {
  "generate=TRUE",			";If FALSE, read file instead",
  "file=",				";File name to read or write",
  "tags=" PosTag "," MassTag "," TypeTag "," KeyTag "," AuxArrTag "[4]",
					";Components of phatbody structure",
  "extend=FALSE",			";If TRUE, add extra field",
  "nbody=4",				";Number of bodies to allocate",
  "VERSION=1.0",			";Joshua Barnes  1 February 2015",
  NULL,
};

#define ExtraField  phatbody[NewBodyFields+0]

int main(int argc, string argv[])
{
  string *tags, intags[32];
  int i, nbody;
  bodyptr btab, bp;
  stream str;
  real tsnap = 0.0;

  initparam(argv, defv);
  if (getbparam("generate")) {
    tags = burststring(getparam("tags"), ",");
    if (getbparam("extend")) {
      if (! set_member(tags, "Extra"))
	printf("%s: add \"Extra\" to tags to test extension\n", getprog());
      new_field(&ExtraField, RealType, "Extra");
      new_field(&ExtraField + 1, NULL, NULL);
    }
    layout_body(tags, Precision, NDIM);
    nbody = getiparam("nbody");
    btab = (bodyptr) allocate(SizeofBody * nbody);
    for (i = 0; i < nbody; i++) {
      if (set_member(tags, PosTag)) {
	UNITV(Pos(NthBody(btab, i)), i % 3);
      }
      if (set_member(tags, MassTag))
	Mass(NthBody(btab, i)) = 1.0 / (1 << i);
      if (set_member(tags, TypeTag))
	Type(NthBody(btab, i)) = 0xf0 | i;
      if (set_member(tags, KeyTag))
	Key(NthBody(btab, i)) = i;
      if (set_member(tags, AuxArrTag))
	AuxArr(NthBody(btab, i))[0] = i + 1;
    }
    if (! strnull(getparam("file"))) {
      str = stropen(getparam("file"), "w");
      put_history(str);
      put_snap(str, &btab, &nbody, &tsnap, tags);
      fclose(str);
    }
  } else {
    if (strnull(getparam("file")))
      error("%s: must supply input file\n", getprog());
    str = stropen(getparam("file"), "r");
    get_history(str);
    get_snap(str, &btab, &nbody, &tsnap, intags, TRUE);
    fclose(str);
  }
  printf("%s: length = %d\n", phatbody[0].name, phatbody[0].length);
  for (i = 1; phatbody[i].name != NULL; i++)
    if (phatbody[i].offset != BadOffset)
      printf("  %s: type = %s  offset = %d  length = %d\n", phatbody[i].name,
	     phatbody[i].type, phatbody[i].offset, phatbody[i].length);
}

#endif
