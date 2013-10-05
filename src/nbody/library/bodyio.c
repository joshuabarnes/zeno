/*
 * bodyio.c: Input/Output code for phat body arrays.
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
local void add_known_fields(stream);
local void get_particles(stream, bodyptr *, int *, string *);

local void set_mask(int *, int, int, int);

//  put_snap: write snapshot to structured output stream.
//  -----------------------------------------------------

void put_snap(stream ostr, bodyptr *btab, int *nbody, real *tnow,
	      string *tags)
{
  put_set(ostr, SnapShotTag);
  put_parameters(ostr, nbody, tnow);
  put_particles(ostr, btab, nbody, tags);
  put_tes(ostr, SnapShotTag);
}

//  put_parameters: write snapshot parameters.
//  ------------------------------------------

local void put_parameters(stream ostr, int *nbody, real *tnow)
{
  put_set(ostr, ParametersTag);
  put_data(ostr, NBodyTag, IntType, nbody, 0);
  put_data(ostr, TimeTag, RealType, tnow, 0);
  put_tes(ostr, ParametersTag);
}

//  put_particles: write particle data to output file.
//  --------------------------------------------------

local void put_particles(stream ostr, bodyptr *btab, int *nbody, string *tags)
{
  string *tp, type;
  ps_field *pbf;
  int mask[32];

  if (*nbody > 0) {
    put_set(ostr, ParticlesTag);
    for (tp = tags; *tp != NULL; tp++) {	// loop over list of tags
      for (pbf = phatbody; pbf->name != NULL; pbf++)
	if (streq(pbf->name, *tp))		// look for name in struct
	  break;
      if (pbf->name == NULL)
	error("%s.put_particles: field %s unknown\n", getprog(), *tp);
      if (pbf->offset == BadOffset)
	error("%s.put_particles: field %s undefined\n", getprog(), *tp);
      set_mask(mask, SizeofBody, pbf->offset, type_length(pbf->type));
      type = type_base(pbf->type);		// get base type of datum
      if (strlen(pbf->type) == 1)
	put_data_masked(ostr, *tp, type, *btab, *nbody, 0, mask);
      else
	put_data_masked(ostr, *tp, type, *btab, *nbody, NDIM, 0, mask);
    }
    put_tes(ostr, ParticlesTag);
  }
}

//  get_snap: read snapshot from structured input stream.
//  -----------------------------------------------------

bool get_snap(stream istr, bodyptr *btab, int *nbody, real *tnow,
	      string *tags, bool expand)
{
  int nbody0;

  if (get_tag_ok(istr, SnapShotTag)) {
    get_set(istr, SnapShotTag);
    if (*btab != NULL)
      nbody0 = *nbody;				// save previous value
    get_parameters(istr, nbody, tnow);
    if (expand && *btab == NULL)		// don't expand after alloc
      add_known_fields(istr);
    if (*btab != NULL && *nbody > nbody0) {	// not enough room for bodies?
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

#ifndef TimeFuzz
#  define TimeFuzz  0.001               /* uncertainty in time comparison   */
#endif

//  get_snap_t: read snapshot in time range from structured input stream.
//  ---------------------------------------------------------------------

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
      if (expand && *btab == NULL)		// don't expand after alloc
	add_known_fields(istr);
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
//  -----------------------------------------

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

//  add_known_fields: add all known body fields present in input set.
//  -----------------------------------------------------------------

local void add_known_fields(stream istr)
{
  string *tp, tags[MaxBodyFields];
  ps_field *pbf;
  
  if (get_tag_ok(istr, ParticlesTag)) {
    get_set(istr, ParticlesTag);
    tp = tags;
    for (pbf = phatbody; pbf->name != NULL; pbf++)
      if (get_tag_ok(istr, pbf->name))
	*tp++ = pbf->name;
    *tp = NULL;
    get_tes(istr, ParticlesTag);
    layout_body(tags, Precision, NDIM);
  }
}

//  get_particles: read particle data from input file.
//  --------------------------------------------------

local void get_particles(stream istr, bodyptr *btab, int *nbody, string *tags)
{
  string *tp = tags, type;
  int len, mask[32];
  ps_field *pbf;

  if (*nbody > 0) {
    if (*btab == NULL) {
      *btab = (bodyptr) allocate(*nbody * SizeofBody);
      eprintf("[%s.get_snap: allocating %d bodies (size %d bytes)]\n",
	      getprog(), *nbody, SizeofBody);
    }
    get_set(istr, ParticlesTag);
    if (get_tag_ok(istr, PhaseTag) &&
	PosField.offset != BadOffset && VelField.offset != BadOffset) {
      len = type_length(PosField.type);
      if (VelField.offset - PosField.offset != len)
	error("%s.get_particles: %s and %s not contiguous",
	      getprog(), PosTag, VelTag);
      set_mask(mask, SizeofBody, PosField.offset, 2 * len);
      type = type_base(PosField.type);
      get_data_masked(istr, PhaseTag, type, *btab, *nbody, 2, NDIM, 0, mask);
      *tp++ = PosTag;
      *tp++ = VelTag;
    }
    for (pbf = phatbody; pbf->name != NULL; pbf++)
      if (get_tag_ok(istr, pbf->name) && pbf->offset != BadOffset) {
	set_mask(mask, SizeofBody, pbf->offset, type_length(pbf->type));
	type = type_base(pbf->type);
	if (strlen(pbf->type) == 1)
	  get_data_masked(istr, pbf->name, type, *btab, *nbody, 0, mask);
	else
	  get_data_masked(istr, pbf->name, type, *btab, *nbody, NDIM, 0, mask);
	*tp++ = pbf->name;
      }
    get_tes(istr, ParticlesTag);
  }
  *tp = NULL;
}

//  set_mask: initialize byte mask for field.
//  -----------------------------------------

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
