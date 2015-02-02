/*
 * phatstruct.h: definitions for "phat" structures.
 */

#ifndef _phatstruct_h
#define _phatstruct_h

//  ps_field: structure describing one element of a phat structure.  A
//  phat structure is described by an array of ps_fields.  The zeroth
//  field refers to the entire structure and its length is the length
//  of the whole structure.  The array is terminated by a ps_field with
//  a NULL name.  Compound data types are supported, but they must be
//  homogenious.
//  ___________________________________________________________________

typedef struct {
    string type;		// type code following "datatypes.h"
    string name;		// text for name of field
    int offset;			// offset in bytes from start of struct
    int length;			// length in bytes of this field
} ps_field;

//  badoffset: offset value for structure fields which aren't actualized.
//  _____________________________________________________________________

#define BadOffset  -1

//  layout_struct: compute offsets and length for structure with named fields.
//  __________________________________________________________________________

void layout_struct(ps_field *, string *);

//  new_field: define a new field of given type and name.
//  _____________________________________________________

void new_field(ps_field *, string, string);

//  define_struct, define_offset: describe total length and individual field
//  offsets of static structure.
//  ________________________________________________________________________

void define_struct(ps_field *, string, int);
void define_offset(ps_field *, string, int);

//  Generic selector macros.
//  ________________________

#define SelectByte(ptr,off)    (*((byte *)   ((byte *)(ptr) + (off))))
#define SelectChar(ptr,off)    (*((char *)   ((byte *)(ptr) + (off))))
#define SelectShort(ptr,off)   (*((short *)  ((byte *)(ptr) + (off))))
#define SelectBool(ptr,off)    (*((bool *)   ((byte *)(ptr) + (off))))
#define SelectInt(ptr,off)     (*((int *)    ((byte *)(ptr) + (off))))
#define SelectLong(ptr,off)    (*((long *)   ((byte *)(ptr) + (off))))
#define SelectFloat(ptr,off)   (*((float *)  ((byte *)(ptr) + (off))))
#define SelectDouble(ptr,off)  (*((double *) ((byte *)(ptr) + (off))))
#define SelectReal(ptr,off)    (*((real *)   ((byte *)(ptr) + (off))))
#define SelectVect(ptr,off)    ( (real *)    ((byte *)(ptr) + (off)))
#define SelectIArr(ptr,off)    ( (int *)     ((byte *)(ptr) + (off)))
#define SelectRArr(ptr,off)    ( (real *)    ((byte *)(ptr) + (off)))

#endif  // ! _phatstruct_h
