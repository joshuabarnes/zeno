/*
 * strset.h: definitions for simple package for string sets.
 */

#ifndef _strset_h
#define _strset_h

//  Sets are represented by null-terminated vectors of string pointers.
//  ___________________________________________________________________

string *set_cons(string, ...);

string *set_copy(string *);

int set_length(string *);

bool set_member(string *, string);

bool set_equal(string *, string *);

bool set_subset(string *, string *);

string *set_union(string *, string *);

string *set_inter(string *, string *);

string *set_diff(string *, string *);

#endif // ! _strset_h
