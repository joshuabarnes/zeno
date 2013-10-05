/*
 * RSF.C/DSF.C: read or dress a structured file.
 */

#include "stdinc.h"
#include "string.h"
#include "getparam.h"
#include "filestruct.h"

string defv[] = {
#ifndef DressStructuredFile
				";Read ascii form of structured file.",
				";Output equivalent binary version.",

    "in=???",			";Name of text input file.",
				";Syntax follows output of tsf utility.",
#else
				";Dress binary data as structured file.",
				";Combines binary data with text format.",
    "form=???",			";Name of text input file.",
				";Syntax follows tsf, with no data.",
    "data=???",			";Name of raw binary data file.",
				";Structure is specified by form file.",
#endif
    "out=???",			";Name of structured binary output file.",
    "append=false",		";If true, append to output file.",
    "VERSION=3.0",		";Josh Barnes  10 February 2000",
    NULL,
};

#define MaxTagLen  256		/* maximum length of item name-tags         */
#define MaxVecDim   24          /* maximum number of array dimensions       */
#define MaxSetCnt   32		/* maximum number of nested sets	    */

stream istr, dstr, ostr;			/* input/output streams     */

void read_file(void);				/* read file item by item   */
bool read_head(char [], char [], int []);	/* read item header info    */
void *read_data(string, int []);		/* read actual item data    */
void read_string(string, int);			/* read character string    */

int main(int argc, string argv[])
{
    initparam(argv, defv);
#ifndef DressStructuredFile
    istr = stropen(getparam("in"), "r");
#else
    istr = stropen(getparam("form"), "r");
    dstr = stropen(getparam("data"), "r");
#endif
    ostr = stropen(getparam("out"), getbparam("append") ? "a" : "w");
    read_file();
    fflush(NULL);
    return (0);
}

/*
 * READ_FILE: read file item by item, keeping track of sets.
 */

void read_file(void)
{
    int sp, dims[MaxVecDim];
    char typ[MaxTagLen], tag[MaxTagLen];
    string setstk[MaxSetCnt];
    void *data;

    sp = 0;
    while (read_head(typ, tag, dims))
	if (streq(typ, SetType)) {
	    if (sp >= MaxSetCnt)
		error("%s: set stack overflow\n", getargv0());
	    setstk[sp] = (string) copxstr(tag, sizeof(char));
	    put_set(ostr, setstk[sp]);
	    sp++;
	} else if (streq(typ, TesType)) {
	    if (sp <= 0)
		error("%s: set stack underflow\n", getargv0());
	    sp--;
	    put_tes(ostr, setstk[sp]);
	    free(setstk[sp]);
	} else {
	    data = read_data(typ, dims);
	    put_data_sub(ostr, tag, typ, dims[0] > 0 ? dims : NULL, data,
			 FALSE);
	    free(data);
	}
    if (sp != 0)
	error("%s: input ended within set %s\n", getargv0(), setstk[sp - 1]);
}

/*
 * READ_HEAD: input the type, tag, and dimensions of an item.
 */

bool read_head(char typ[], char tag[], int dims[])
{
    char word[MaxTagLen];
    int n;

    if (fscanf(istr, " %s", word) == EOF)
	return (FALSE);
    if (streq(word, "tes")) {
	(void) strcpy(typ, TesType);
	return (TRUE);
    }
    if (fscanf(istr, " %[^[ \n]", tag) == EOF)
	error("%s: unexpected EOF\n", getargv0());
    if (streq(word, "set")) {
	(void) strcpy(typ, SetType);
	return (TRUE);
    }
    (void) strcpy(typ, name_type(word));
    for (n = 0; fscanf(istr, " [%i]", &dims[n]) == 1; n++) {
	if (dims[n] < 1)
	    error("%s: illegal array dimension\n", getargv0());
	if (n == MaxVecDim - 1)
	    error("%s: too many dimensions\n", getargv0());
    }
    dims[n] = 0;
    return (TRUE);
}

/*
 * READ_DATA: input actual data of an item in ascii or binary.
 */

void *read_data(string typ, int dims[])
{
    int ndatum, nbytes, i, tmpb;
    void *data;

    ndatum = 1;
    for (i = 0; dims[i] != 0; i++)
	ndatum *= dims[i];
    nbytes = ndatum * type_length(typ);
    data = allocate(nbytes);
#ifndef DressStructuredFile
    if (streq(typ, AnyType) || streq(typ, ByteType)) {
	for (i = 0; i < ndatum; i++) {
	    if (fscanf(istr, " %i", &tmpb) != 1)
		error("%s: error reading byte\n", getargv0());
	    ((byte *) data)[i] = tmpb;
	}
    } else if (streq(typ, CharType) && dims[0] > 0 && dims[1] == 0)
	read_string((string) data, ndatum);
    else if (streq(typ, CharType)) {
	for (i = 0; i < ndatum; i++)
	    if (fscanf(istr, " %c", &((char *) data)[i]) != 1)
		error("%s: error reading char\n", getargv0());
    } else if (streq(typ, ShortType)) {
	for (i = 0; i < ndatum; i++)
	    if (fscanf(istr, " %hi", &((short *) data)[i]) != 1)
		error("%s: error reading short\n", getargv0());
    } else if (streq(typ, IntType)) {
	for (i = 0; i < ndatum; i++)
	    if (fscanf(istr, " %i", &((int *) data)[i]) != 1)
		error("%s: error reading int\n", getargv0());
    } else if (streq(typ, LongType)) {
	for (i = 0; i < ndatum; i++)
	    if (fscanf(istr, " %li", &((long *) data)[i]) != 1)
		error("%s: error reading long\n", getargv0());
    } else if (streq(typ, FloatType)) {
	for (i = 0; i < ndatum; i++)
	    if (fscanf(istr, " %f", &((float *) data)[i]) != 1)
		error("%s: error reading float\n", getargv0());
    } else if (streq(typ, DoubleType)) {
	for (i = 0; i < ndatum; i++)
	    if (fscanf(istr, " %lf", &((double *) data)[i]) != 1)
		error("%s: error reading double\n", getargv0());
    } else
	error("%s: undefined data type\n", getargv0());
#else
    if (fread(data, sizeof(byte), nbytes, dstr) != nbytes)
	error("%s: fread failed\n", getargv0());
#endif
    return (data);
}

/*
 * READ_STRING: read a double-quote delimited string.
 */

void read_string(string sbuf, int nchar)
{
    char tmpc;
    int i, j;
    string bakset = "abfnrtv\"\\", bakmap = "\a\b\f\n\r\t\v\"\\";

    if (fscanf(istr, " %c", &tmpc) != 1)
	error("%s: unexpected EOF\n", getargv0());
    if (tmpc != '\"')
	error("%s: initial double-quote missing\n", getargv0());
    i = 0;
    while ((tmpc = fgetc(istr)) != '\"') {
	if (feof(istr))
	    error("%s: unexpected EOF inside string\n", getargv0());
	if (tmpc == '\\') {
	    tmpc = fgetc(istr);
	    if (feof(istr))
		error("%s: unexpected EOF after backslash\n", getargv0());
	    for (j = 0; bakset[j] != tmpc; j++)
	      if (bakset[j] == (char) NULL)
		    error("%s: undefined backslash substitution\n",
			  getargv0());
	    tmpc = bakmap[j];
	}
	if (i < nchar)
	    sbuf[i++] = tmpc;
    }
    if (i == nchar)
	eprintf("[%s: truncating string to %d chars + NULL]\n",
		getargv0(), nchar - 1);
    sbuf[nchar - 1] = (char) NULL;
    while (i < nchar - 1)
        sbuf[i++] = (char) NULL;
}
