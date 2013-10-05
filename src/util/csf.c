/*
 * CSF.C: copy or strip a structured file.
 */

#include "stdinc.h"
#include "strset.h"
#include "getparam.h"
#include "filestruct.h"

string defv[] = {
#if defined(StripStructuredFile)
				";Strip a structured binary file.",
				";Outputs raw binary data.",
#elif defined(GrabStructuredFiles)
				";Grab structured binary files.",
				";Reads multiple files in index order.",
#else
				";Copy a structured binary file.",
				";Allows control over items copied.",
#endif
#if !defined(GrabStructuredFiles)
    "in=???",			";Name of input file",
    "out=???",			";Name of output file",
#else
    "in=???",			";Printf pattern for input files.",
				";Integers between <start> and <end>",
				";are substituted in input pattern.",
    "out=???",			";Name of output file",
    "start=???",		";First value to substitute in <in>",
    "end=???",			";Last value to substitute in <in>",
    "step=1",			";Absolute size of step for count.",
				";If <end> < <start>, count down.",
#endif
    "include=",			";List of included item names.",
				";If given, only these are copied.",
    "exclude=",			";List of excluded item names.",
				";If given, these are not copied.",
    "skip=",			";Skip this many included items.",
				";Counts top-level items to skip.",
    "copy=",			";Copy this many included items.",
				";Counts top-level items to copy.",
    "append=false",		";If true, append to output file",
    "swap=false",		";If true, swap order of bytes.",
				";A temporary fix for endian problem.",
    "poll=",			";Sleep this long, then retry open",
#if !defined(GrabStructuredFiles)
    "VERSION=1.5",		";Josh Barnes  12 April 2000",
#else
    "VERSION=1.6",		";Josh Barnes  24 April 2011",
#endif
    NULL,
};

stream instr;
stream outstr = NULL;
string *inc;
string *exc;
bool swap;

bool sift_item(string);
void copy_item_sifted(string);
void swap_order(byte *, int, string);

#if !defined(GrabStructuredFiles)

int main(int argc, string argv[])
{
    int nskip, ncopy, totread, totsift, totcopy;
    string tag;

    initparam(argv, defv);
    inc = burststring(getparam("include"), ", ");
    exc = burststring(getparam("exclude"), ", ");
    nskip = (strnull(getparam("skip")) ? 0 : getiparam("skip"));
    ncopy = (strnull(getparam("copy")) ? -1 : getiparam("copy"));
    swap = getbparam("swap");
    totread = totsift = totcopy = 0;
    if (strnull(getparam("poll")))
	instr = stropen(getparam("in"), "r");
    else {
	while ((instr = stropen(getparam("in"), "r?")) == NULL)
	    sleep((unsigned) ABS(getiparam("poll")));
	sleep((unsigned) ABS(getiparam("poll")));
    }
    outstr = stropen(getparam("out"), getbparam("append") ? "a" : "w");
    while (totcopy != ncopy && (tag = next_item_tag(instr)) != NULL) {
	totread++;
	if (sift_item(tag)) {
	    totsift++;
	    if (totsift > nskip) {
		copy_item_sifted(tag);
		totcopy++;
	    } else
		skip_item(instr);
	} else
	    skip_item(instr);
	free(tag);
	fflush(outstr);
    }
    strclose(instr);
    eprintf("[%s: read %d, included %d, copied %d toplevel items]\n",
	    getargv0(), totread, totsift, totcopy);
    return (0);
}

#endif

#if defined(GrabStructuredFiles)

int main(int argc, string argv[])
{
  int nskip, ncopy, totread, totsift, totcopy;
  int start, end, step, i, ksift, kcopy;
  char name[256];
  string tag;

  initparam(argv, defv);
  inc = burststring(getparam("include"), ", ");
  exc = burststring(getparam("exclude"), ", ");
  nskip = (strnull(getparam("skip")) ? 0 : getiparam("skip"));
  ncopy = (strnull(getparam("copy")) ? -1 : getiparam("copy"));
  swap = getbparam("swap");
  totread = totsift = totcopy = 0;
  start = getiparam("start");
  end = getiparam("end");
  step = ABS(getiparam("step"));
  for (i = MIN(start, end); i <= MAX(start, end); i += MAX(step, 1)) {
    sprintf(name, getparam("in"), start < end ? i : end + start - i);
    if (strnull(getparam("poll"))) {
      instr = stropen(name, "r?");
      if (instr == NULL && step > 0) {
	eprintf("[%s: can\'t open file \"%s\"]\n", getargv0(), name);
	break;
      }
    } else {
      while ((instr = stropen(name, "r?")) == NULL)
	sleep((unsigned) ABS(getiparam("poll")));
      sleep((unsigned) ABS(getiparam("poll")));
    }
    if (instr != NULL) {
      if (outstr == NULL)
	outstr = stropen(getparam("out"), getbparam("append") ? "a" : "w");
      ksift = kcopy = 0;
      while (kcopy != ncopy && (tag = next_item_tag(instr)) != NULL) {
	totread++;
	if (sift_item(tag)) {
	  ksift++;
	  if (ksift > nskip) {
	    copy_item_sifted(tag);
	    kcopy++;
	  } else
	    skip_item(instr);
	} else
	  skip_item(instr);
	free(tag);
	fflush(outstr);
      }
      strclose(instr);
      totsift += ksift;
      totcopy += kcopy;
    }
  }
  eprintf("[%s: read %d, included %d, copied %d toplevel items]\n",
	  getargv0(), totread, totsift, totcopy);
  return (0);
}

#endif

/*
 * SIFT_ITEM: decide if item should be copied to output.
 */

bool sift_item(string tag)
{
    if (*inc != NULL)
	return (set_member(inc, tag));
    if (*exc != NULL)
	return (! set_member(exc, tag));
    return (TRUE);
}

/*
 * COPY_ITEM_SIFTED: item/set copy, with sifting and other options.
 */

void copy_item_sifted(string tag)
{
    string type, tag1;
    int *dims, dlen;
    byte *buf;

    if (! get_tag_ok(instr, tag))
	error("copy_item_sifted: tag %s not found\n", tag);
    type = get_type(instr, tag);
    if (! streq(type, SetType)) {
	dims = get_dimensions(instr, tag);
	dlen = get_length(instr, tag);
	buf = (byte *) allocate(dlen);
	get_data_sub(instr, tag, type, dims, (void *) buf, NULL);
	if (swap)
	    swap_order(buf, dlen, type);
#ifndef StripStructuredFile
	put_data_sub(outstr, tag, type, dims, (void *) buf, NULL);
#else
	if (fwrite(buf, sizeof(byte), dlen, outstr) != dlen)
	    error("%s: fwrite failed\n", getargv0());
#endif
	if (dims != NULL)
	    free(dims);
	free(buf);
    } else {
	get_set(instr, tag);
#ifndef StripStructuredFile
	put_set(outstr, tag);
#endif
	while ((tag1 = next_item_tag(instr)) != NULL) {
	    if (sift_item(tag1))
		copy_item_sifted(tag1);
	    else
		skip_item(instr);
	    free(tag1);
	}
	get_tes(instr, tag);
#ifndef StripStructuredFile
	put_tes(outstr, tag);
#endif
    }
    free(type);
}

/*
 * SWAP_ORDER: reorder bytes within basic data elements.
 */

void swap_order(byte *buf, int dlen, string type)
{
    int blen, nbyt, kbyt;
    byte *p2, *p1, btmp;

    blen = type_length(type_base(type));
    if (dlen % blen != 0)
        error("%s: internal error: type-length mismatch\n", getargv0());
    for (nbyt = dlen ; nbyt > 0; nbyt -= blen) {
	p2 = buf + nbyt;
	p1 = p2 - blen;
	for (kbyt = blen / 2; kbyt > 0; kbyt--) {
	    btmp = *p1;
	    *p1++ = *--p2;
	    *p2 = btmp;
	}
    }
}
