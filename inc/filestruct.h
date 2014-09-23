/*
 * filestruct.h: interface to structured binary file package.
 *       Version 1 by Josh Barnes & Lyman Hurd, IAS, 1987.
 *       Version 2 by Josh Barnes, IAS, 1988.
 *       Version 3 by Josh Barnes, IfA, March 1994.
 */
 
#ifndef _filestruct_h
#define _filestruct_h
 
#include "datatypes.h"
 
#define SetType    "("          // begin compound item
#define TesType    ")"          // end of compound item

//  Routines for direction-symmetric I/O code.
//  __________________________________________
 
void copy_item(stream, stream, string);

void put_set(stream, string);
void put_tes(stream, string);
void put_data(stream, string, string, void *, ...);
void put_data_masked(stream, string, string, void *, ...);
void put_data_sub(stream, string, string, int *, void *, int *);
 
void get_set(stream, string);
void get_tes(stream, string);
void get_data(stream, string, string, void *, ...);
void get_data_masked(stream, string, string, void *, ...);
void get_data_sub(stream, string, string, int *, void *, int *);
 
void put_string(stream, string, string);
string get_string(stream, string);

//  Routines for data-driven input code.
//  ____________________________________
 
bool get_tag_ok(stream, string);
 
string get_type(stream, string);
int *get_dimensions(stream, string);
int get_length(stream, string);

string next_item_tag(stream);
bool skip_item(stream);
 
//  Assorted control routines and parameters.
//  _________________________________________
 
void strclose(stream);
 
void fs_options(int, int);

#define NoChange  0		// don't change setting
#define NotAllow  1		// forbid operation
#define WarnEach  2		// warn user each time
#define WarnOnce  3		// warn user first time
#define SilentOK  4		// allow without warning
 
#endif  // !  _filestruct_h
