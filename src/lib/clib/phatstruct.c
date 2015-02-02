/*
 * phatstruct.c: routines for phat structure package.
 */

#include "stdinc.h"
#include "getparam.h"
#include "datatypes.h"
#include "phatstruct.h"
#include <string.h>

//  layout_struct: compute offsets and length of structure.  The order
//  of the actual fields is determined by the names argument.  Padding
//  is computed assuming (1) any type of object can have offset 0, and
//  (2) identical objects can be stored contiguously without padding.
//  __________________________________________________________________

void layout_struct(ps_field *pstab, string *names)
{
  bool debug = (getenv("ZENO_PHSTR_DEBUG") != NULL);
  int pad, len;
  ps_field *psp;

  pad = 0;
  while (*names != NULL) {			// loop over field names
    for (psp = pstab + 1; psp->name != NULL; psp++)
      if (streq(psp->name, *names))		// find name in struct tab
	break;
    if (psp->name == NULL)			// must already be known
      error("%s.layout_struct: field %s not found\n", getprog(), *names);
    if (psp->offset == BadOffset) {		// if not already defined
      len = type_length(type_base(psp->type));	// get length of base type
      while (pstab->length % len != 0) {	// align on proper boundary
	pstab->length++;
	pad++;
      }
      psp->offset = pstab->length;		// assign offset and length
      psp->length = type_length(psp->type);
      pstab->length += psp->length;		// increment struct length
      if (debug)
	eprintf("[%s.layout_struct: name = %s  offset = %d  length = %d]\n",
		getprog(), *names, psp->offset, psp->length);
    }
    names++;
  }
  len = sizeof(int);				// take int as minimum
  for (psp = pstab + 1; psp->name != NULL; psp++)
    if (psp->offset != BadOffset)
      len = MAX(len, type_length(type_base(psp->type)));
  while (pstab->length % len != 0) {		// align complete structure
    pstab->length++;
    pad++;
  }
  if (debug)
    eprintf("[%s.layout_struct: sizeof(%s) = %d bytes (%d padding)]\n",
	    getprog(), pstab->name, pstab->length, pad);
}

//  new_field: define a new field of given type and name.  Note that it is up
//  to the client to make sure that the field array is properly terminated.
//  _________________________________________________________________________

void new_field(ps_field *psptr, string type, string name)
{
  psptr->name = name;				// name is given by caller
  psptr->type = type;				// type is given by caller
  psptr->offset = BadOffset;			// offset is yet unknown
  psptr->length = 0;				// and length is undefined
  if (getenv("ZENO_PHSTR_DEBUG") != NULL)
    eprintf("[%s.new_field: type = %s  name = %s]\n", getprog(), type, name);
}

//  define_struct: set total length of structure; a partial alternative
//  to layout_struct for the fixed offset interface.
//  ___________________________________________________________________

void define_struct(ps_field *pstab, string name, int length)
{
  if (! streq(pstab->name, name))
    error("%s.define_struct: structure %s not found\n", getprog(), name);
  pstab->length = length;
  if (getenv("ZENO_PHSTR_DEBUG") != NULL)
    eprintf("[%s.define_struct: sizeof(%s) = %d bytes]\n",
	    getprog(), name, length);
}

//  define_offset: set offset of known field; a partial alternative to
//  layout_struct for the fixed offset interface.
//  __________________________________________________________________

void define_offset(ps_field *pstab, string name, int offset)
{
  ps_field *psp;

  for (psp = pstab + 1; psp->name != NULL; psp++)
    if (streq(psp->name, name))			// find name in struct tab
      break;
  if (psp->name == NULL)			// must be already known
    error("%s.define_field: field %s not found\n", getprog(), name);
  if (psp->offset != BadOffset)			// but not actual
    error("%s.define_field: field %s already defined\n",
	  getprog(), name);
  psp->offset = offset;				// store supplied offset
  psp->length = type_length(psp->type);
  if (getenv("ZENO_PHSTR_DEBUG") != NULL)
    eprintf("[%s.define_offset: name = %s  offset = %d  length = %d]\n",
	    getprog(), name, psp->offset, psp->length);
}    

#ifdef TESTBED

#include "getparam.h"

string defv[] = {		";Test phat structures",
    "fields=foo,bar,fum,fie,baz,jojo",
				";List of structure fields",
    "VERSION=1.0",		";Josh Barnes  18 July 1994",
    NULL,
};

ps_field phatstruct[] = {
    { NULL,               "foobar",  0,          0 },
    { CharType,           "foo",     BadOffset,  0 },
    { ShortType,          "bar",     BadOffset,  0 },
    { IntType,            "fum",     BadOffset,  0 },
    { FloatType,          "fie",     BadOffset,  0 },
    { DoubleType,         "baz",     BadOffset,  0 },
    { RealType RealType,  "jojo",    BadOffset,  0 },
    { NULL,               NULL,      0,          0 },
};

#define Foo(x)  SelectChar(x, phatstruct[1].offset)
#define Bar(x)  SelectShort(x, phatstruct[2].offset)
#define Fum(x)  SelectInt(x, phatstruct[3].offset)
#define Fie(x)  SelectFloat(x, phatstruct[4].offset)
#define Baz(x)  SelectDouble(x, phatstruct[5].offset)
#define Jojo(x) SelectVect(x, phatstruct[6].offset)

main(int argc, string argv[])
{
  string *fields;
  void *xp;

  initparam(argv, defv);
  fields = burststring(getparam("fields"), ",");
  layout_struct(phatstruct, fields);
  xp = allocate(phatstruct[0].length);
  Foo(xp) = 'a';
  Bar(xp) = 123;
  Fum(xp) = 12345678;
  Fie(xp) = 3.141592;
  Baz(xp) = 2.718281;
  Jojo(xp)[0] = 1.5;
  Jojo(xp)[1] = -2.5;
  printf("xp -> { %c  %d  %d  %f  %f [%f %f] }\n", Foo(xp), Bar(xp),
	 Fum(xp), Fie(xp), Baz(xp), Jojo(xp)[0], Jojo(xp)[1]);
}

#endif
