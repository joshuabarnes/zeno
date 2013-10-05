/*
 * TSF: type the contents of a structured file in human-readable form.
 * V 1.?  Joshua Barnes, Peter Teuben, Lyman Hurd  1986-7.
 * V 2.0  Joshua Barnes  4/88  new filestruct package.
 * V 2.2  Joshua Barnes  4/94  new new filestruct package.
 */

#include "stdinc.h"
#include "string.h"
#include "getparam.h"
#include "filestruct.h"

string defv[] = {               ";Type out a structured binary file.",
				";Lists contents in human-readable form.",
    "<in=???",                  ";Name of structured file for input.",
    "maxprec=false",		";If true, list with full precision.",
				";Useful for transfer of ascii data.",
    "maxline=4",                ";List this many lines for each item.",
				";Longer items are truncated with ellipsis.",
    "indent=2",			";Indentation of compound items.",
    "margin=72",		";Righthand margin for formatting.",
    "VERSION=2.3",		";Josh Barnes  10 February 2000",
    NULL,
};

stream instr;				/* input stream from struct. file   */
bool maxprec;				/* if true, use max precision forms */
int  maxline;				/* maximum number of lines per item */
int  indent;				/* indentation for compound items   */
int  margin;				/* righthand margin for formatting  */

void print_item(string);		/* print named item to std output   */
void print_set(string);			/* print start of item set	    */
void print_tes(string);			/* print end of item set	    */
void print_head(string, string, int *);	/* print item header and dimensions */
void print_data(string, string, int *);	/* print actual data of item        */

bool outstr(string);
void endline(void);
void incindent(int delta);

int main(int argc, string argv[])
{
    string tag;

    initparam(argv, defv);
    instr = stropen(getparam("in"), "r");
    maxprec = getbparam("maxprec");
    maxline = getiparam("maxline");
    indent = getiparam("indent");
    margin = getiparam("margin");
    while ((tag = next_item_tag(instr)) != NULL) {
        print_item(tag);
	free(tag);
    }
    return (0);
}

void print_item(string tag)
{
    string type, stag;
    int *dims;

    type = get_type(instr, tag);
    if (streq(type, SetType)) {
	get_set(instr, tag);
	print_set(tag);
	while ((stag = next_item_tag(instr)) != NULL) {
	    print_item(stag);
	    free(stag);
	}
	get_tes(instr, tag);
	print_tes(tag);
    } else {
	dims = get_dimensions(instr, tag);
	print_head(tag, type, dims);
	print_data(tag, type, dims);
	endline();
	if (dims != NULL)
	    free(dims);
    }
    free(type);
}

void print_set(string tag)
{
    char buf[128];

    sprintf(buf, "set %s", tag);
    (void) outstr(buf);
    endline();
    incindent(indent);
}
    
void print_tes(string tag)
{
    char buf[128];

    incindent(- indent);
    sprintf(buf, "tes");
    (void) outstr(buf);
    endline();
}

void print_head(string tag, string type, int *dims)
{
    char buf[128];
    int *dp;

    if (strlen(type) == 1) {
	sprintf(buf, "%s %s", type_name(type), tag);
	(void) outstr(buf);
    } else {
	(void) outstr("struct { ");
	while (*type != (char) NULL) {
	    sprintf(buf, "%s ", type_name(type++));
	    (void) outstr(buf);
	}
	sprintf(buf, "} %s", tag);
	(void) outstr(buf);
    }
    if (dims != NULL)				/* is this a plural item?   */
	for (dp = dims; *dp != 0; dp++) {	/*   loop over dimensions   */
            sprintf(buf, "[%d]", *dp);		/*     format a dimension   */
            (void) outstr(buf);                 /*     and print it out     */
        }
}

void print_data(string tag, string type, int *dims)
{
    int dlen;
    byte *dat, *dp;
    char *tp, buf[128];
    string pad, fmt;

    dlen = get_length(instr, tag);
    dat = (byte *) allocate(dlen);
    get_data_sub(instr, tag, type, dims, dat, NULL);
    if (streq(type, CharType) && dims != NULL && dims[1] == 0) {
	(void) outstr(" \"");
	pad = "";
    } else
	pad = " ";
    tp = type;					/* start type string scan   */
    for (dp = dat; dp < dat + dlen; ) {		/* loop over data array     */
	if (! outstr(pad))			/*   print pad before data  */
	    break;
	fmt = type_fmt(tp, !maxprec);		/*   get proper format	    */
	if (tp[0] == AnyType[0]) {		/*   output generic data?   */
	    sprintf(buf, fmt, *((byte *) dp));
	    dp += sizeof(byte);
	} else if (tp[0] == CharType[0]) {	/*   output readable chars? */
	    sprintf(buf, fmt, *((char *) dp));
	    dp += sizeof(char);
	} else if (tp[0] == ByteType[0]) {	/*   output bytes of data?  */
	    sprintf(buf, fmt, *((byte *) dp));
	    dp += sizeof(byte);
	} else if (tp[0] == ShortType[0]) {	/*   output short integers? */
	    sprintf(buf, fmt, *((short *) dp));
	    dp += sizeof(short);
	} else if (tp[0] == IntType[0]) {	/*   output standard ints?  */
	    sprintf(buf, fmt, *((int *) dp));
	    dp += sizeof(int);
	} else if (tp[0] == LongType[0]) {	/*   output long integers?  */
	    sprintf(buf, fmt, *((long *) dp));
	    dp += sizeof(long);
	} else if (tp[0] == FloatType[0]) {	/*   output floating point? */
	    sprintf(buf, fmt, *((float *) dp));
	    dp += sizeof(float);
	} else if (tp[0] == DoubleType[0]) {	/*   output double numbers? */
	    sprintf(buf, fmt, *((double *) dp));
	    dp += sizeof(double);
	} else
	    error("print_data: type %s unknown\n", type);
	if (! outstr(buf))
	    break;
	tp = (*(tp + 1) != (char) NULL ? tp + 1 : type);
    }
    if (streq(type, CharType) && dims != NULL && dims[1] == 0)
	(void) outstr("\"");
    free(dat);
}

/*
 * State variables for formatting functions.
 */

int curcoll = 0;	/* current cursor column, numberd from 0 */
int curline = 0;	/* current line of item, numbered from 0 */

int baseindent = 0;	/* indentation of item */
int actindent;		/* actual indentation of line */

bool outstr(string str)
{
    if (curcoll + strlen(str) >= margin) {      /* string too long to fit?  */
	printf("\n");				/*   begin next line        */
	curline++;				/*   update line counter    */
	curcoll = 0;				/*   and column counter     */
	actindent = baseindent + indent;	/*   indent text of item    */
	if ((!maxprec) && (curline == maxline))	/*   time to truncate?      */
	    str = ". . .";			/*     indicate truncation  */
    }
    for ( ; curcoll < actindent; curcoll++)
	printf(" ");
    printf("%s", str);                          /* output the string        */
    curcoll = curcoll + strlen(str);            /* and count its length     */
    return (maxprec || (curline < maxline));	/* give go-ahead for more   */
}

void endline(void)
{
    printf("\n");
    curcoll = curline = 0;
    actindent = baseindent;
}

void incindent(int delta)
{
    baseindent = baseindent + delta;
    actindent = baseindent;
}
