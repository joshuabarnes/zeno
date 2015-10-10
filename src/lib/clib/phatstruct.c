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

// Lines with TEST commants are prototyping construction of type string
// for structure as realized in memory.

void layout_struct(ps_field *pstab, string *names)
{
  bool debug = (getenv("ZENO_PHSTR_DEBUG") != NULL);
  int pad = 0, len;
  ps_field *psp;

  pstab->type = (string) allocate(512);		// TEST: alloc lots of space
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
	strcat(pstab->type, "b");		// TEST: store padding byte
      }
      psp->offset = pstab->length;		// assign offset and length
      psp->length = type_length(psp->type);
      pstab->length += psp->length;		// increment struct length
      if (debug)
	eprintf("[%s.layout_struct: name = %s  offset = %d  length = %d]\n",
		getprog(), *names, psp->offset, psp->length);
      strcat(pstab->type, psp->type);	// TEST: append field type
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
    strcat(pstab->type, "b");			// TEST: store padding byte
  }
  if (debug)
    eprintf("[%s.layout_struct: struct %s: size = %d   pad = %d  type = %s]\n",
	    getprog(), pstab->name, pstab->length, pad, pstab->type);
						// TEST: print combined type
}

//  new_field: define a new field of given type and name.
//  _____________________________________________________

void new_field(ps_field *psp, string type, string name)
{
  bool debug = (getenv("ZENO_PHSTR_DEBUG") != NULL);

  if (debug)
    eprintf("[%s.new_field: type = %s  name = %s]\n", getprog(), type, name);
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
  for (psp = pstab + 1; psp->name != NULL; psp++)
    if (streq(psp->name, name))			// find name in struct tab
      break;
  if (psp->name == NULL)			// must be already known
    error("%s.define_field: field %s not found\n", getprog(), name);
  if (psp->offset != BadOffset)			// but not already defined
    error("%s.define_field: can't redefine field %s\n", getprog(), name);
  psp->offset = offset;				// store supplied offset
  psp->length = type_length(psp->type);
  if (debug)
    eprintf(" type = %s  length = %d]\n", psp->type, psp->length);
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
  string *fields;
  void *xp;
  stream outstr;

  initparam(argv, defv);
  fields = burststring(getparam("fields"), ",");
  layout_struct(pstab, fields);
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
