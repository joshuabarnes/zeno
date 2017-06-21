/*
 * phatstruct.c: routines for phat structure package.
 */

#include "stdinc.h"
#include "getparam.h"
#include "datatypes.h"
#include "phatstruct.h"
#include <string.h>

local ps_field *find_field(ps_field *pstab, string name, string caller);
local int pad_struct(ps_field *pstab, int len);

//  layout_struct: compute offsets and length of structure.  The order
//  of the active fields is determined by the names argument.  Padding
//  is computed assuming (1) any object can have offset 0, and
//  (2) identical objects can be stored contiguously without padding.
//  Lines with TEST comments prototype construction of structure type.
//  __________________________________________________________________

void layout_struct(ps_field *pstab, string *names)
{
  bool debug = (getenv("ZENO_PHSTR_DEBUG") != NULL);
  int len, pad = 0;
  ps_field *psp;

  if (pstab->type == NULL)			// TEST: first call on struct?
    pstab->type = (string) allocate(512);	// TEST: alloc lots of space
  for (string *np = names; *np != NULL; np++) {	// loop over field names
    psp = find_field(pstab, *np, "layout_struct");
    if (psp->offset == BadOffset) {		// if field not yet defined
      if (type_base(psp->type) == NULL)
	error("%s.layout_struct: base type undefined; field %s, type %s\n",
	      getprog(), psp->name, psp->type);
      len = type_length(type_base(psp->type));	// get length of base type
      pad += pad_struct(pstab, len);		// align on field boundary
      psp->offset = pstab->length;		// offset makes field active
      psp->length = type_length(psp->type);	// store total field length
      pstab->length += psp->length;		// keep track of struct length
      strcat(pstab->type, psp->type);		// TEST: append field type
      if (debug)
	eprintf("[%s.layout_struct: field %s:  offset = %d  length = %d"
		"  type = %s]\n", getprog(), *np, psp->offset, psp->length,
		psp->type);
    }
  }
  len = sizeof(int);				// align structs to ints
  for (psp = pstab+1; psp->name != NULL; psp++)	// loop over active fields
    len = MAX(len, psp->offset != BadOffset ?	// find max component length
	      type_length(type_base(psp->type)) : 0);
  pad += pad_struct(pstab, len);		// align on struct boundary
  if (debug)
    eprintf("[%s.layout_struct: struct %s: size = %d  type = %s  pad = %d"
	    "  maxlen = %d]\n", getprog(), pstab->name, pstab->length,
	    pstab->type, pad, len);		// TEST: print combined type
}

//  new_field: define a new field of given type and name.
//  _____________________________________________________

void new_field(ps_field *psp, string type, string name)
{
  bool debug = (getenv("ZENO_PHSTR_DEBUG") != NULL);

  if (debug)
    eprintf("[%s.new_field: type = %s  name = %s]\n", getprog(),
	    type != NULL ? type : "NULL", name != NULL ? name : "NULL");
  psp->name = name;				// name is given by caller
  psp->type = type;				// type is given by caller
  psp->offset = BadOffset;			// offset is yet unknown
  psp->length = 0;				// and length is undefined
  if (name != NULL && type != NULL)		// defining real field?
    (psp+1)->name = NULL;			// next field ends array
}

//  define_struct: set total length of structure.
//  Partial alternative to layout_struct for fixed-offset interface.
//  ________________________________________________________________

void define_struct(ps_field *pstab, string name, int length)
{
  bool debug = (getenv("ZENO_PHSTR_DEBUG") != NULL);

  if (debug)
    eprintf("[%s.define_struct: sizeof(%s) = %d bytes]\n",
	    getprog(), name, length);
  if (! streq(pstab->name, name))		// check given name matches
    error("%s.define_struct: structure %s unknown\n", getprog(), name);
  pstab->length = length;			// set total structure length
}

//  define_offset: set offset of pre-defined field.
//  Partial alternative to layout_struct for fixed-offset interface.
//  ________________________________________________________________

void define_offset(ps_field *pstab, string name, int offset)
{
  bool debug = (getenv("ZENO_PHSTR_DEBUG") != NULL);
  ps_field *psp;

  if (debug)
    eprintf("[%s.define_offset: name = %s  offset = %d\n",
	    getprog(), name, offset);
  psp = find_field(pstab, name, "define_field");
  if (psp->offset != BadOffset)			// but not already defined
    error("%s.define_field: can't redefine field %s\n", getprog(), name);
  psp->offset = offset;				// store supplied offset
  psp->length = type_length(psp->type);
  if (debug)
    eprintf(" type = %s  length = %d]\n", psp->type, psp->length);
}    

//  find_field: scan structure table to find named field.
//  _____________________________________________________

local ps_field *find_field(ps_field *pstab, string name, string caller)
{
  ps_field *psp;

  for (psp = pstab + 1; psp->name != NULL; psp++)
    if (streq(psp->name, name))			// find name in struct tab
      break;
  if (psp->name == NULL)			// must be already known
    error("%s.%s: field %s not found\n", getprog(), caller, name);
  return (psp);
}

//  pad_struct: round up struct length to align on object of length len.
//  ____________________________________________________________________

local int pad_struct(ps_field *pstab, int len)
{
  int pad = 0;

  while (pstab->length % len != 0) {		// until structure is aligned
    pstab->length++;				// increase structure length
    strcat(pstab->type, "b");			// TEST: store padding byte
    pad++;					// and count padding
  }
  return (pad);
}

#ifdef TESTBED

#include "getparam.h"
#include "filestruct.h"

string defv[] = {		";Test phat structures",
  "fields=foo,bar,fum,fie,baz,jojo",
				";List of structure fields",
  "out=",			";Optional output file",
  "VERSION=1.1",		";Josh Barnes  16 July 2015",
  NULL,
};

ps_field pstab[] = {
  { NULL,               "foobar",  0,          0 },
  { CharType,           "foo",     BadOffset,  0 },
  { ShortType,          "bar",     BadOffset,  0 },
  { IntType,            "fum",     BadOffset,  0 },
  { FloatType,          "fie",     BadOffset,  0 },
  { DoubleType,         "baz",     BadOffset,  0 },
  { RealType RealType,  "jojo",    BadOffset,  0 },
  { NULL,               NULL,      0,          0 },
};

#define Foo(x)  SelectChar(x, pstab[1].offset)
#define Bar(x)  SelectShort(x, pstab[2].offset)
#define Fum(x)  SelectInt(x, pstab[3].offset)
#define Fie(x)  SelectFloat(x, pstab[4].offset)
#define Baz(x)  SelectDouble(x, pstab[5].offset)
#define Jojo(x) SelectVect(x, pstab[6].offset)

int main(int argc, string argv[])
{
  void *xp;
  stream outstr;

  initparam(argv, defv);
  layout_struct(pstab, burststring(getparam("fields"), ","));
  xp = allocate(pstab[0].length);
  if (pstab[1].offset != BadOffset) Foo(xp) = 'a';
  if (pstab[2].offset != BadOffset) Bar(xp) = 123;
  if (pstab[3].offset != BadOffset) Fum(xp) = 12345678;
  if (pstab[4].offset != BadOffset) Fie(xp) = 3.141592;
  if (pstab[5].offset != BadOffset) Baz(xp) = 2.718281;
  if (pstab[6].offset != BadOffset) Jojo(xp)[0] = 1.5;
  if (pstab[6].offset != BadOffset) Jojo(xp)[1] = -2.5;
  printf("xp -> { ");
  if (pstab[1].offset != BadOffset) printf("%c ", Foo(xp));
  if (pstab[2].offset != BadOffset) printf("%d ", Bar(xp));
  if (pstab[3].offset != BadOffset) printf("%d ", Fum(xp));
  if (pstab[4].offset != BadOffset) printf("%f ", Fie(xp));
  if (pstab[5].offset != BadOffset) printf("%f ", Baz(xp));
  if (pstab[6].offset != BadOffset) 
    printf("[%f %f] ", Jojo(xp)[0], Jojo(xp)[1]);
  printf("}\n");
  if (! strnull(getparam("out"))) {
    outstr = stropen(getparam("out"), "w");
    put_data(outstr, pstab[0].name,  pstab[0].type, xp, NULL);
  }
  return (0);
}

#endif
