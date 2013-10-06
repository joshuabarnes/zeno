/*
 * strset.c: code for simple string set package.
 */

#include "stdinc.h"
#include "strset.h"
#include "getparam.h"

#include <stdarg.h>

//  Sets of strings are represented by null-terminated vectors of string
//  pointers.  This strategy keeps the code simple, but also implies
//  that many operations (construction, union, intersection, difference)
//  are O(N^2) in the number of elements.  The maximum number of elements
//  is therefore limited to a relatively small value.
//  _____________________________________________________________________

#define MaxElements  64

//  set_cons: construct a set from a null-terminated list of elements.
//  __________________________________________________________________

string *set_cons(string first, ...)
{
  int n;
  string name, strs[MaxElements+1];
  va_list ap;

  n = 0;
  strs[n++] = first;
  strs[n] = NULL;
  va_start(ap, first);
  while ((name = va_arg(ap, string)) != NULL)
    if (set_member(strs, name))
      eprintf("[%s.set_cons: WARNING: element %s duplicated]\n",
	      getprog(), name);
    else {
      if (n == MaxElements)
	error("%s.set_cons: too many elements\n", getprog());
      strs[n++] = name;
      strs[n] = NULL;
    }
  va_end(ap);
  return (set_copy(strs));
}

//  set_copy: copy set (but not member strings).
//  ____________________________________________

string *set_copy(string *set)
{
  return ((string *) copxstr(set, sizeof(string)));
}

//  set_length: return count of set members.
//  ________________________________________

int set_length(string *set)
{
  return (xstrlen(set, sizeof(string)) - 1);
}

//  set_member: test if element is member of set.
//  _____________________________________________

bool set_member(string *set, string element)
{
  string *sp;

  for (sp = set; *sp != NULL; sp++)
    if (streq(*sp, element))
      return (TRUE);
  return (FALSE);
}

//  set_subset: return true if every member of set2 is in set1.
//  ___________________________________________________________

bool set_subset(string *set1, string *set2)
{
  string *sp;

  for (sp = set2; *sp != NULL; sp++)
    if (! set_member(set1, *sp))
      return (FALSE);
  return (TRUE);
}

//  set_equal: return true if sets have same members.
//  _________________________________________________

bool set_equal(string *set1, string *set2)
{
  return (set_subset(set1, set2) && set_subset(set2, set1));
}

//  set_union: construct union of two sets.
//  _______________________________________

string *set_union(string *set1, string *set2)
{
  int n;
  string *sp, strs[MaxElements+1];

  n = 0;
  for (sp = set1; *sp != NULL; sp++)
    strs[n++] = *sp;
  for (sp = set2; *sp != NULL; sp++)
    if (! set_member(set1, *sp)) {
      if (n > MaxElements)
	error("%s.set_union: too many elements\n", getprog());
      strs[n++] = *sp;
    }
  strs[n] = NULL;
  return (set_copy(strs));
}

//  set_inter: construct intersection of two sets.
//  ______________________________________________

string *set_inter(string *set1, string *set2)
{
  int n;
  string *sp, strs[MaxElements+1];

  n = 0;
  for (sp = set2; *sp != NULL; sp++)
    if (set_member(set1, *sp))
      strs[n++] = *sp;
  strs[n] = NULL;
  return (set_copy(strs));
}

//  set_diff: construct difference of two sets.
//  ___________________________________________

string *set_diff(string *set1, string *set2)
{
  int n;
  string *sp, strs[MaxElements+1];

  n = 0;
  for (sp = set1; *sp != NULL; sp++)
    if (! set_member(set2, *sp))
      strs[n++] = *sp;
  strs[n] = NULL;
  return (set_copy(strs));
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
    "set1=foo,bar,fum,fie",
    "set2=fie,bar,waldo,frodo",
    "VERSION=1.0",
    NULL,
};

void print_set(string, string *);

void main(int argc, string argv[])
{
  string *set1, *set2;

  initparam(argv, defv);
  set1 = burststring(getparam("set1"), ", ");
  set2 = burststring(getparam("set2"), ", ");
  printf("set_length(set1) = %d\n", set_length(set1));
  printf("set_subset(set1,set2) = %s\n",
	 set_subset(set1, set2) ? "TRUE" : "FALSE");
  printf("set_equal(set1,set2) = %s\n",
	 set_equal(set1, set2) ? "TRUE" : "FALSE");
  print_set("set_copy(set1):", set_copy(set1));
  print_set("set_union(set1,set2):", set_union(set1, set2));
  print_set("set_inter(set1,set2):", set_inter(set1, set2));
  print_set("set_diff(set1,set2):", set_diff(set1, set2));
}

void print_set(string label, string *set)
{
  string *sp;

  printf("%s", label);
  for (sp = set; *sp != NULL; sp++)
    printf(" %s", *sp);
  printf("\n");
}

#endif // TESTBED
