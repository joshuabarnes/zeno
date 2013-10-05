/*
 * PHATBODY.C: Procedures to define and extend phat body structures.
 */

#include "stdinc.h"
#include "getparam.h"
#include "bodytags.h"
#include "phatbody.h"
#include "datatypes.h"
#include <string.h>

/*
 * Temporary tokens for real and vector types, replaced by local routine
 * set_types() once client's precision and dimension number are determined.
 */

#define RealToke "R"				/* placeholder for RealType */
#define VectToke "V"				/* placeholder for VectType */

local void set_types(string, int);		/* replace w/ proper types  */

/*
 * PHATBODY: Defines the default components of a body structure.
 * Most applications will use only some components. 
 */

ps_field phatbody[MaxBodyFields] = {
    { NULL,      BodyTag,     0,          0 },		/* = phatbody[ 0]   */
    { VectToke,  PosTag,      BadOffset,  0 },		/* = phatbody[ 1]   */
    { VectToke,  VelTag,      BadOffset,  0 },		/* = phatbody[ 2]   */
    { RealToke,  MassTag,     BadOffset,  0 },		/* = phatbody[ 3]   */
    { RealToke,  PhiTag,      BadOffset,  0 },		/* = phatbody[ 4]   */
    { VectToke,  AccTag,      BadOffset,  0 },		/* = phatbody[ 5]   */
    { RealToke,  SmoothTag,   BadOffset,  0 },		/* = phatbody[ 6]   */
    { RealToke,  RhoTag,      BadOffset,  0 },		/* = phatbody[ 7]   */
    { RealToke,  EntFuncTag,  BadOffset,  0 },		/* = phatbody[ 8]   */
    { RealToke,  UinternTag,  BadOffset,  0 },		/* = phatbody[ 9]   */
    { RealToke,  UdotIntTag,  BadOffset,  0 },		/* = phatbody[10]   */
    { RealToke,  UdotRadTag,  BadOffset,  0 },		/* = phatbody[11]   */
    { RealToke,  UdotVisTag,  BadOffset,  0 },		/* = phatbody[12]   */
    { RealToke,  TauTag,      BadOffset,  0 },		/* = phatbody[13]   */
    { ByteType,  TypeTag,     BadOffset,  0 },		/* = phatbody[14]   */
    { RealToke,  BirthTag,    BadOffset,  0 },		/* = phatbody[15]   */
    { RealToke,  DeathTag,    BadOffset,  0 },		/* = phatbody[16]   */
    { IntType,   KeyTag,      BadOffset,  0 },		/* = phatbody[17]   */
    { RealToke,  AuxTag,      BadOffset,  0 },		/* = phatbody[18]   */
    { VectToke,  AuxVecTag,   BadOffset,  0 },		/* = phatbody[19]   */
    { NULL,      NULL,        0,          0 },	
};

/*
 * LAYOUT_BODY: define real and vector types according to client,
 * and add named fields to body structure.
 */

void layout_body(string *tags, string prec, int ndim)
{
    set_types(prec, ndim);
    layout_struct(phatbody, tags);
}

/*
 * DEFINE_BODY: define real and vector types according to client,
 * and sets body length.
 */

void define_body(int length, string prec, int ndim)
{
    set_types(prec, ndim);
    define_struct(phatbody, "Body", length);
}

/*
 * DEFINE_BODY_OFFSET: fixed-offset interface sets field offset.
 */

void define_body_offset(string tag, int offset)
{
    define_offset(phatbody, tag, offset);
}

/*
 * SET_TYPES: local routine does work of replacing RealToke and VectToke
 * with correct run-time specifications.
 */

local void set_types(string prec, int ndim)
{
    string rtyp, vtyp;
    int i;

    if (streq(prec, "SINGLEPREC"))
	rtyp = FloatType;
    else if (streq(prec, "MIXEDPREC"))
	rtyp = FloatType;
    else if (streq(prec, "DOUBLEPREC"))
	rtyp = DoubleType;
    else
	error("%s.set_types: unknown precision %s\n", getargv0(), prec);
    vtyp = (string) allocate(ndim + 1);
    for (i = 0; i < ndim; i++)
	vtyp[i] = rtyp[0];
    for (i = 1; phatbody[i].name != NULL; i++)
	if (streq(phatbody[i].type, RealToke))
	    phatbody[i].type = rtyp;
	else if (streq(phatbody[i].type, VectToke))
	    phatbody[i].type = vtyp;
}
