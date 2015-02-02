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
  ps_field *pbf;
  int mask[32];

  if (*nbody > 0) {
    put_set(ostr, ParticlesTag);
    for (tp = tags; *tp != NULL; tp++) {	// loop over list of tags
      tag1 = clean_tag(*tp);
      for (pbf = phatbody; pbf->name != NULL; pbf++)
	if (streq(pbf->name, tag1))		// look for name in struct
	  break;
      if (pbf->name == NULL)
	error("%s.put_particles: field %s unknown\n", getprog(), tag1);
      if (pbf->offset == BadOffset)
	error("%s.put_particles: field %s undefined\n", getprog(), tag1);
      set_mask(mask, SizeofBody, pbf->offset, type_length(pbf->type));
      type = type_base(pbf->type);		// get base type of datum
      if (strlen(pbf->type) == 1)
	put_data_masked(ostr, tag1, type, *btab, *nbody, 0, mask);
      else
	put_data_masked(ostr, tag1, type, *btab, *nbody,
			strlen(pbf->type), 0, mask);
      free(tag1);
    }
    put_tes(ostr, ParticlesTag);
  }
}

//  get_snap: read snapshot from structured input stream.
//  _____________________________________________________

bool get_snap(stream istr, bodyptr *btab, int *nbody, real *tnow,
	      string *tags, bool expand)
{
  int nbody1;

  if (get_tag_ok(istr, SnapShotTag)) {
    get_set(istr, SnapShotTag);
    if (*btab != NULL)
      nbody1 = *nbody;				// save previous value
    get_parameters(istr, nbody, tnow);
    if (*btab == NULL && expand)		// permit expand before alloc
      add_body_fields(istr);
    if (*btab != NULL && *nbody > nbody1) {	// not enough room for bodies?
      eprintf("[%s.get_snap: WARNING: deallocating old body array]\n",
	      getprog());
      free(*btab);
      *btab = NULL;				// trigger array reallocation
    }
    get_particles(istr, btab, nbody, tags);
    get_tes(istr, SnapShotTag);
    return (TRUE);
  } else
    return (FALSE);    
}

//  get_snap_t: read snapshot in time range from structured input stream.
//  _____________________________________________________________________

bool get_snap_t(stream istr, bodyptr *btab, int *nbody, real *tnow,
		string *tags, bool expand, string times)
{
  bool success = FALSE;
  int nbody1;
  real tnow1;

  while (! success && get_tag_ok(istr, SnapShotTag)) {
    get_set(istr, SnapShotTag);
    get_parameters(istr, &nbody1, &tnow1);
    if (streq(times, "all") || within(tnow1, times, TimeFuzz)) {
      if (*btab == NULL && expand)		// permit expand before alloc
	add_body_fields(istr);
      if (*btab != NULL && nbody1 > *nbody) {	// not enough room for bodies?
	eprintf("[%s.get_snap_t: WARNING: deallocating old body array]\n",
		getprog());
	free(*btab);
	*btab = NULL;				// trigger array reallocation
      }
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
  ps_field *pbf;
  string type, *tp = tags;
  int mask[32], *dims;

  if (*nbody > 0) {
    if (*btab == NULL) {
      *btab = (bodyptr) allocate(*nbody * SizeofBody);
      eprintf("[%s.get_snap: allocating %d bodies (size %d bytes)]\n",
	      getprog(), *nbody, SizeofBody);
    }
    get_set(istr, ParticlesTag);
    for (pbf = phatbody; pbf->name != NULL; pbf++)
      if (get_tag_ok(istr, pbf->name) && pbf->offset != BadOffset) {
	set_mask(mask, SizeofBody, pbf->offset, type_length(pbf->type));
	type = type_base(pbf->type);
	if (strlen(pbf->type) == 1)
	  get_data_masked(istr, pbf->name, type, *btab, *nbody, 0, mask);
	else
	  get_data_masked(istr, pbf->name, type, *btab, *nbody,
			  strlen(pbf->type), 0, mask);
	if (strne(pbf->name, AuxArrTag) && strne(pbf->name, KeyArrTag))
	  *tp++ = pbf->name;
	else {
	  dims = get_dimensions(istr, pbf->name);
	  asprintf(tp++, "%s[%d]", pbf->name, dims[1]);
	  free(dims);
	}
      }
    get_tes(istr, ParticlesTag);
  }
  *tp = NULL;
}

//  add_body_fields: add all known body fields present in input set.
//  ________________________________________________________________

local void add_body_fields(stream istr)
{
  string tags[MaxBodyFields], *tp;
  ps_field *pbf;
  int *dims;
  
  if (get_tag_ok(istr, ParticlesTag)) {
    if (getenv("ZENO_PHSTR_DEBUG") != NULL)
      eprintf("[%s.add_body_fields: scanning items in %s set]\n",
	      getprog(), ParticlesTag);
    get_set(istr, ParticlesTag);
    tp = tags;
    for (pbf = phatbody; pbf->name != NULL; pbf++)
      if (get_tag_ok(istr, pbf->name)) {
	if (strne(pbf->name, AuxArrTag) && strne(pbf->name, KeyArrTag))
	  *tp++ = pbf->name;
	else {
	  dims = get_dimensions(istr, pbf->name);
	  asprintf(tp++, "%s[%d]", pbf->name, dims[1]);
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
  int i;

  i = 0;
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
  return (tag1);
}
