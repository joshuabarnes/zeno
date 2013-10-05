/*
 * datatypes.h: include file for data type package.
 */

#ifndef _datatypes_h
#define _datatypes_h

//  ___________________________________________________
//  Basic data types, and synonyms for bools and reals.
				
#define AnyType    "a"          // anything at all
#define CharType   "c"          // printable chars
#define ByteType   "b"          // unsigned 8-bit
#define ShortType  "s"          // short integers
#define IntType    "i"          // standard integers
#define LongType   "l"          // long integers
#define FloatType  "f"          // short floating
#define DoubleType "d"          // long floating
 
#define BoolType ShortType      // true/false values
 
#if defined(SINGLEPREC) || defined(MIXEDPREC)
#  define RealType FloatType    // short floating point
#else
#  define RealType DoubleType   // long floating point
#endif
 
//  ________________________________________
//  Utility functions for data type package.

int type_length(string);			// length of data in bytes

string type_name(string);			// C name of basic type

string name_type(string);			// basic type of C name

string type_fmt(string, bool);			// formatting for printf

string type_base(string);			// base of compound type

#endif  // ! _datatypes_h
