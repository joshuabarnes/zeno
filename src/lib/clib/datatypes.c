/*
 * datatypes.c: table and utility routines for data type package.
 */

#include "stdinc.h"
#include "getparam.h"
#include "datatypes.h"
#include <string.h>

//  typeinfo: struct for info about data types.
//  ___________________________________________

typedef struct {
  string type;					// name from datatypes.h
  string name;					// C-style decleration
  string fmt1;					// low-precision format
  string fmt2;					// high-precision format
  int size;					// size of type in bytes
} typeinfo;

//  typetable: table of basic types and associated data.
//  ____________________________________________________

local typeinfo typetable[] = {
  { AnyType,    "any",    "%#o",  "%#o",     sizeof(byte),   },
  { CharType,   "char",   "%c",   "%c",      sizeof(char),   },
  { ByteType,   "byte",   "%#o",  "%#o",     sizeof(byte),   },
  { ShortType,  "short",  "%#o",  "%#o",     sizeof(short),  },
  { IntType,    "int",    "%#o",  "%#o",     sizeof(int),    },
  { LongType,   "long",   "%#lo", "%#lo",    sizeof(long),   },
  { FloatType,  "float",  "%#g",  "%15.8e",  sizeof(float),  },
  { DoubleType, "double", "%#g",  "%23.16e", sizeof(double), },
  { NULL,       NULL,     NULL,   NULL,      0,              },
};

local typeinfo *find_type(string);		// lookup by 1st element
local typeinfo *find_name(string);		// lookup by C-style name

//  type_length: compute total size of basic or compound type.
//  __________________________________________________________

int type_length(string type)
{
  int count;
  char *tp;

  count = 0;
  for (tp = type; *tp != (char) NULL; tp++)	// loop over all elements
    count += find_type(tp)->size;		// add size of each element
  return (count);
}

//  type_name: return C-style name of type.
//  _______________________________________

string type_name(string type)
{
  return (find_type(type)->name);
}

//  name_type: return basic type of C-style name.
//  _____________________________________________

string name_type(string name)
{
  return (find_name(name)->type);
}

//  type_fmt: return format string for type.
//  ________________________________________

string type_fmt(string type, bool lowprec)
{
  return (lowprec ? find_type(type)->fmt1 : find_type(type)->fmt2);
}

//  type_base: return basic type of homogenious compound type.
//  __________________________________________________________

string type_base(string type)
{
  char *tp;

  for (tp = type; *tp != (char) NULL; tp++)	// scan over compound type
    if (*tp != type[0])				// if type is inhomogeneous
      return (NULL);				// base type is undefined
  return (tp - 1);				// ret. pointer to last char
}

//  find_type: return typeinfo for given type.  Note: for compound
//  types, the typeinfo for the first element is returned.
//  ______________________________________________________________

local typeinfo *find_type(string type)
{
  typeinfo *tp;

  for (tp = typetable; tp->type != NULL; tp++)
    if (tp->type[0] == type[0])
      return (tp);
  error("%s.find_type: type %c unknown\n", getprog(), type[0]);
  return (NULL);
}

//  find_name: return typeinfo for named type.
//  __________________________________________

local typeinfo *find_name(string name)
{
  typeinfo *tp;

  for (tp = typetable; tp->type != NULL; tp++)
    if (streq(tp->name, name))
      return (tp);
  error("%s.find_name: type %s unknown\n", getprog(), name);
  return (NULL);
}

#ifdef TESTBED

string defv[] = {		";Test datatypes routines",
    "type=" IntType,		";Data type string",
    "name=int",			";C-style type name",
    "VERSION=1",		";Josh Barnes  15 October 1997",
    NULL,
};

void main(int argc, string argv[])
{
  string type, name;

  initparam(argv, defv);
  type = getparam("type");
  name = getparam("name");
  printf("type_length(\"%s\") = %d\n", type, type_length(type));
  printf("type_name(\"%s\") = %s\n", type, type_name(type));
  printf("type_base(\"%s\") = %s\n", type, type_base(type));
  printf("type_fmt(\"%s\", TRUE) = %s\n", type, type_fmt(type, TRUE));
  printf("type_fmt(\"%s\", FALSE) = %s\n", type, type_fmt(type, FALSE));
  printf("name_type(\"%s\") = %s\n", name, name_type(name));
}

#endif
