/*
 * bodyio.c: Input/output code for phat body arrays.
 */

#include "stdinc.h"
#include "filestruct.h"
#include "vectdefs.h"
#include "phatbody.h"
#include "getparam.h"
#include <string.h>

local void put_parameters(stream, int *, real *);
local void put_particles(stream, bodyptr *, int *, string *);

local void get_parameters(stream, int *, real *);
local void get_particles(stream, bodyptr *, int *, string *);

local void add_body_fields(stream);
local void set_mask(int *, int, int, int);
local string clean_tag(string tag);
local void list_fields(int *ndef, int offbuf[], char tagbuf[], int taglen);

#ifndef TimeFuzz
#  define TimeFuzz  0.001		// tolerance in time comparison
#endif

//  put_snap: write snapshot to structured output stream.
//  _____________________________________________________

void put_snap(stream ostr, bodyptr *btab, int *nbody, real *tnow,
	      string *tags)
{
  put_set(ostr, SnapShotTag);
  put_parameters(ostr, nbody, tnow);
  if (*nbody > 0)
    put_particles(ostr, btab, nbody, tags);
  put_tes(ostr, SnapShotTag);
}

//  put_parameters: write snapshot parameters.
//  __________________________________________

local void put_parameters(stream ostr, int *nbody, real *tnow)
{
  put_set(ostr, ParametersTag);
  put_data(ostr, NBodyTag, IntType, nbody, 0);
  put_data(ostr, TimeTag, RealType, tnow, 0);
  put_tes(ostr, ParametersTag);
}

//  put_particles: write particle data.
//  ___________________________________

local void put_particles(stream ostr, bodyptr *btab, int *nbody, string *tags)
{
  string *tp, tag1, type;
  ps_field *fp;
  int mask[32], ndef, offbuf[MaxBodyFields];
  char tagbuf[512];

  put_set(ostr, ParticlesTag);
  if (tags != NULL) {				// got list of field tags?
    for (tp = tags; *tp != NULL; tp++) {	// loop over list of tags
      tag1 = clean_tag(*tp);
      for (fp = phatbody; fp->name == NULL || strne(fp->name, tag1); fp++)
	if (fp->name == NULL)
	  error("%s.put_particles: field %s unknown\n", getprog(), tag1);
      if (fp->offset == BadOffset)
	error("%s.put_particles: field %s undefined\n", getprog(), tag1);
      set_mask(mask, SizeofBody, fp->offset, type_length(fp->type));
      type = type_base(fp->type);		// get base type of datum
      if (strlen(fp->type) == 1)
	put_data_masked(ostr, tag1, type, *btab, *nbody, 0, mask);
      else
	put_data_masked(ostr, tag1, type, *btab, *nbody,
			strlen(fp->type), 0, mask);
      free(tag1);
    }
  } else {					// output entire body array
    list_fields(&ndef, offbuf, tagbuf, sizeof(tagbuf));
    put_string(ostr, BodyFieldsTag, tagbuf);
    put_data(ostr, BodyOffsetsTag, IntType, offbuf, ndef, 0);
    put_data(ostr, BodyArrayTag, BodyField.type, *btab, *nbody, 0);
  }
  put_tes(ostr, ParticlesTag);
}

//  get_snap: read snapshot from structured input stream.
//  _____________________________________________________

bool get_snap(stream istr, bodyptr *btab, int *nbody, real *tnow,
	      string *tags, bool expand, string times)
{
  bool success = FALSE;
  int nbody1;
  real tnow1;

  while (! success && get_tag_ok(istr, SnapShotTag)) {
    get_set(istr, SnapShotTag);
    get_parameters(istr, &nbody1, &tnow1);
    if (times == NULL || streq(times, "all") ||
	within(tnow1, times, TimeFuzz)) {
      if (*btab == NULL && expand)		// permit expand before alloc
	add_body_fields(istr);
      if (*btab != NULL && nbody1 > *nbody) {	// not enough room for bodies?
	eprintf("[%s.get_snap: WARNING: reallocating body array"
		"  nbody = %d -> %d]\n", getprog(), *nbody, nbody1);
	free(*btab);
	*btab = NULL;				// trigger array reallocation
      }
      if (nbody1 < *nbody)
	eprintf("[%s.get_snap: shrinking body array"
		"  nbody = %d -> %d]\n", getprog(), *nbody, nbody1);
      *nbody = nbody1;
      *tnow = tnow1;
      get_particles(istr, btab, nbody, tags);
      success = TRUE;
    }
    get_tes(istr, SnapShotTag);
  }
  return (success);
}

//  get_parameters: read snapshot parameters.
//  _________________________________________

local void get_parameters(stream istr, int *nbody, real *tnow)
{
  get_set(istr, ParametersTag);
  if (get_tag_ok(istr, NBodyTag))
    get_data(istr, NBodyTag, IntType, nbody, 0);
  else if (get_tag_ok(istr, NobjTag))
    get_data(istr, NobjTag, IntType, nbody, 0);
  else
    error("%s.get_snap: need %s or %s in %s\n",
	  getprog(), NBodyTag, NobjTag, ParametersTag);
  if (get_tag_ok(istr, TimeTag))
    get_data(istr, TimeTag, RealType, tnow, 0);
  else {
    *tnow = 0.0;
    eprintf("[%s.get_snap: time defaults to zero]\n", getprog());
  }
  get_tes(istr, ParametersTag);
}

//  get_particles: read particle data.
//  __________________________________

local void get_particles(stream istr, bodyptr *btab, int *nbody, string *tags)
{
  string *tp = tags, type;
  int mask[32], *dims;

  if (*nbody <= 0) {
    *tp = NULL;
    return;
  }
  if (*btab == NULL) {
    *btab = (bodyptr) allocate(*nbody * SizeofBody);
    eprintf("[%s.get_snap: allocating %d bodies (size %d bytes)]\n",
	    getprog(), *nbody, SizeofBody);
  }
  get_set(istr, ParticlesTag);
  if (! get_tag_ok(istr, BodyArrayTag)) {
    for (ps_field *fp = phatbody; fp->name != NULL; fp++)
      if (get_tag_ok(istr, fp->name) && fp->offset != BadOffset) {
	set_mask(mask, SizeofBody, fp->offset, type_length(fp->type));
	if (strlen(fp->type) == 1)
	  get_data_masked(istr, fp->name, type_base(fp->type),
			  *btab, *nbody, 0, mask);
	else
	  get_data_masked(istr, fp->name, type_base(fp->type),
			  *btab, *nbody, strlen(fp->type), 0, mask);
	if (strne(fp->name, AuxArrTag) && strne(fp->name, KeyArrTag))
	  *tp++ = fp->name;
	else {
	  dims = get_dimensions(istr, fp->name);
	  asprintf(tp++, "%s[%d]", fp->name, dims[1]);
	  free(dims);
	}
      }
  } else {
    type = get_type(istr, BodyArrayTag);
    eprintf("[%s.get_snap: %s type = %s, length = %d  SizeofBody = %d]\n",
	    getprog(), BodyArrayTag, type, type_length(type), SizeofBody);
    error("%s.get_snap: %s input not implemented\n", getprog(), BodyArrayTag);
  }
  get_tes(istr, ParticlesTag);
  *tp = NULL;
}

//  add_body_fields: add all known body fields present in input set.
//  ________________________________________________________________

local void add_body_fields(stream istr)
{
  string tags[MaxBodyFields], *tp;
  ps_field *fp;
  int *dims;
  
  if (get_tag_ok(istr, ParticlesTag)) {
    if (getenv("ZENO_PHSTR_DEBUG") != NULL)
      eprintf("[%s.add_body_fields: scanning items in %s set]\n",
	      getprog(), ParticlesTag);
    get_set(istr, ParticlesTag);
    tp = tags;
    for (fp = phatbody; fp->name != NULL; fp++)
      if (get_tag_ok(istr, fp->name)) {
	if (strne(fp->name, AuxArrTag) && strne(fp->name, KeyArrTag))
	  *tp++ = fp->name;
	else {
	  dims = get_dimensions(istr, fp->name);
	  asprintf(tp++, "%s[%d]", fp->name, dims[1]);
	}
      }
    *tp = NULL;
    get_tes(istr, ParticlesTag);
    layout_body(tags, Precision, NDIM);
  }
}

//  set_mask: initialize byte mask to read or write specified field.
//  ________________________________________________________________

local void set_mask(int *mask, int size, int offset, int length)
{
  int i = 0;

  if (offset > 0)
    mask[i++] = -offset;			// skip preceding stuff
  mask[i++] = length;				// copy field itself
  mask[i++] = -(size - offset - length);	// skip rest of structure
  mask[i] = 0;
}

//  clean_tag: copy tag without length specification.
//  _________________________________________________

local string clean_tag(string tag0)
{
  string tag1 = (string) copxstr(tag0, sizeof(char));
  if (index(tag1, '[') != NULL)
    *index(tag1, '[') = (char) NULL;
  return tag1;
}

//  list_fields: count defined fields, sort by offset, and return lists.
//  ____________________________________________________________________

local void list_fields(int *ndef, int offbuf[], char tagbuf[], int taglen)
{
  int i, j, k;

  *ndef = 0;					// count defined fields
  for (ps_field *fp = &phatbody[1]; fp->name != NULL; fp++)
    if (fp->offset != BadOffset)
      offbuf[(*ndef)++] = fp->offset;		// store in buffer
  for (i = 1; i < *ndef; i++) {			// sort offsets low to high
    j = i;
    while (j > 0 && offbuf[j-1] > offbuf[j]) {
      k = offbuf[j-1];
      offbuf[j-1] = offbuf[j];
      offbuf[j] = k;
      j--;
    }
  }
  tagbuf[0] = (char) NULL;			// init to zero length
  for (i = 0; i < *ndef; i++)
    for (ps_field *fp = &phatbody[1]; fp->name != NULL; fp++)
      if (fp->offset == offbuf[i]) {
#if defined(MACOSX)
	if (tagbuf[0] != (char) NULL)		// if not zero length
	  strlcat(tagbuf, ",", taglen);		// then insert a comma
	strlcat(tagbuf, fp->name, taglen);
#else
	if (tagbuf[0] != (char) NULL)		// if not zero length
	  strcat(tagbuf, ",");			// then insert a comma
	strcat(tagbuf, fp->name);		// should use strncat()
#endif
      }
}
