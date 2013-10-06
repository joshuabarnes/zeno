/*
 * history.c: routines to handle the history items of structured binary files.
 *	 9-Mar-88	V1.0	Created					PJT
 *	 1-Jun-88	V1.1	renamed read/write to get/put		PJT
 *	 7-jun-88	V1.2	installed in libT with getparam();
 *                              dynamic HISTORY				PJT
 *      14-Jun-88       V1.3    integrated with IAS NEMO		JEB
 *       5-Aug-96       V1.4    removed headlines; added prototyping    JEB
 */

#include "stdinc.h"
#include "getparam.h"
#include "filestruct.h"

#define HistoryTag "History"		// tag for history items	

#define MAXHIST 256			// max size of history array

local string histbuf[MAXHIST+1];	// history string array

local int nhist = 0;			// count of history data stored

//  get_history: read history information into history buffer.
//  __________________________________________________________

void get_history(stream instr)
{
  bool lpflag, inflag;

  lpflag = TRUE;				// get into loop once
  inflag = TRUE;				// read 1st history item
  while (lpflag)				// loop reading input data
    if (get_tag_ok(instr, HistoryTag))		// got a history item?
      if (inflag) {				// and first in this file?
	add_history(get_string(instr, HistoryTag));
						// then add to record
	inflag = FALSE;
      } else
	free(get_string(instr, HistoryTag));	// else, skip history data
    else
      lpflag = FALSE;				// don't take next loop
}

//  skip_history: read past history items.
//  ______________________________________

void skip_history(stream instr)
{
  while (get_tag_ok(instr, HistoryTag))
    (void) skip_item(instr);
}

//  put_history:  write current history and headline data to output.
//  ________________________________________________________________

void put_history(stream outstr)
{				
  int i;

  for (i = 0; i < nhist; i++)
    put_string(outstr, HistoryTag, histbuf[i]);
}

//  add_history: add item to history array.
//  _______________________________________

void add_history(string s)
{
  if (nhist < MAXHIST-1)			// enough room for data?
    histbuf[nhist] = s;
  else if (nhist == MAXHIST-1)			// just reached maximum?
    eprintf("[%s.add_history: WARNING: too much history]\n", getprog());
  nhist++;
  histbuf[MIN(nhist,MAXHIST-1)] = NULL;		// terminate w/ NULL
}

//  ask_history: enquire about history data.
//  ________________________________________

string *ask_history(void)
{
  return (histbuf);
}

#ifdef TESTBED

char *defv[] = {		// DEFAULT INPUT PARAMETERS
    "in=???",			// ascii input file name
    "out=",			// if given, output history
    "history=",			// if given, add to history array
    "VERSION=1.3",		// JEB  14 June 1988
    NULL,
};

string iname, oname;			// input/output file  name
stream instr, outstr;

main(int argc, string argv[])
{
  int i;

  initparam(argv, defv);			// init command line pars
  iname = getparam("in");
  oname = getparam("out");
  instr = stropen(iname, "r");
  get_history(instr);
  if (! streq(getparam("history"), ""))
    add_history(getparam("history"));
  for (i = 0; i < nhist; i++)
    printf("%3d: %s\n", i, histbuf[i]);
  strclose(instr);
  if (! streq(oname, "")) {
    outstr = stropen(oname, "w");
    put_history(outstr);
    strclose(outstr);
  }
}

#endif
