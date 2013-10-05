/*
 * filestruct.c: structured binary file package.
 *       Version 1 by Josh Barnes & Lyman Hurd, IAS, 1987.
 *       Version 2 by Josh Barnes, IAS, 1988.
 *       Version 3 by Josh Barnes, IfA, April 1994.
 */
 
#include "stdinc.h"
#include "getparam.h"
#include <stdarg.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "filestruct.h"

//  ______________________________________________
//  Magic numbers used to identify items in files.
 
#define ScalarMagic  ((011<<8) + 0222)          // for single items
#define PluralMagic  ((013<<8) + 0222)          // for plural items
 
//  __________________________________________
//  item: structure representing a data-token.
 
typedef struct {
    string type;                // specifies type of data; see .h file
    string tag;                 // name for data in external file
    int *dims;                  // int-string of dimensions, or NULL
    void *data;                 // actual data associated with item
    long datapos;               // where data begins in input stream
} item, *itemptr;    
 
//  _____________________________________________________
//  state: structure representing a structured I/O state.
 
typedef struct {
    stream str;                 // stream associated with state
    itemptr *buf;               // buffer of items read or written
    int maxbuf;                 // maximum number of items buffer can hold
    int itemcnt;                // number of items in buffer
    int context;                // next item to access during input
    bool canseek;               // true if seek operations are allowed
} state, *stateptr;
 
//  _________________________________________________________________________
//  Storage limits -- may be increased without rendering data files obsolete.
 
#define MaxVecDim   24          // maximum number of array dimensions
#define MaxStream   32          // maximum number of active streams
#define BufSizeInc  32          // initial item buffer size, and increment
#define MaxDataBuf 256          // read longer items on request, if we can

//  _______________________________
//  Local routines and definitions.
 
local stateptr findstate(stream);
local bool readset(stateptr);
local void bufitem(stateptr, itemptr);
local int finditem(stateptr, string);
local int findset(stateptr, int);
local int findtes(stateptr, int);
local void putitem(stateptr, itemptr, int *);
local itemptr getitem(stateptr);
local void getdata(itemptr, stateptr);
local void copydata(itemptr, int *, itemptr, stream);
local void copyint(itemptr, int *, itemptr, stream, int, int, int);
local void copybyte(itemptr, int *, itemptr, stream, int, int, int);
local void fillxbuf(byte **, int *, int *, int *, itemptr, stream, int);
local void inputdata(void *, itemptr, stream, int, int);
local bool checkdims(itemptr, itemptr);
local int convmode(itemptr, itemptr);
local int datalength(itemptr);
local int datacount(itemptr);
local itemptr makeitem(string, string, int *, void *);
local void freeitem(itemptr);
local void safewrite(stream, void *, int);
local void saferead(stream, void *, int);
local void safeseek(stream, long, int);
 
local string iofunc;			// for msgs from public routines
local int    iofcnt = 0;		// count calls between routines
 
#define set_iofunc(msg) if (iofcnt++ == 0) iofunc = (msg)
#define end_iofunc()    if (--iofcnt == 0) iofunc = NULL
 
//  __________________________________
//  External routines and definitions.
 
void free(void *);
 
#define cstrlen(x)  (xstrlen((x), sizeof(char)))
#define istrlen(x)  (xstrlen((x), sizeof(int)))
#define copycstr(x)  ((string) copxstr((x), sizeof(char)))
#define copyistr(x)  ((int *) copxstr((x), sizeof(int)))
#define getcstr(x)  ((string) getxstr((x), sizeof(char)))
#define getistr(x)  ((int *) getxstr((x), sizeof(int)))
#define istreq(x,y)  (xstreq((x), (y), sizeof(int)))

//
//  USER COPY FUNCTION
//
 
//  ____________________________________________________________
//  copy_item: recursively copy item from input to output.
//  An example of recursive file traversal and memory etiquette.

void copy_item(stream ostr, stream istr, string tag)
{
    string type, tag1;
    int *dims, dlen;
    byte *buf;

    set_iofunc("copy_item");
    if (! get_tag_ok(istr, tag))		// prevent obvious errors
	error("%s.%s: tag %s not found\n", getprog(), iofunc, tag);
    type = get_type(istr, tag);			// find out type of data
    if (! streq(type, SetType)) {		// a basic type or array?
	dims = get_dimensions(istr, tag);	// find out about shape
	dlen = get_length(istr, tag);		// and length in bytes
	buf = allocate(dlen);			// get space for a buffer
	get_data_sub(istr, tag, type, dims, (void *) buf, NULL);
						// read data from input
	put_data_sub(ostr, tag, type, dims, (void *) buf, NULL);
						// and write it to output
	if (dims != NULL)			// free dimension list
	    free(dims);
	free(buf);				// free temporary buffer
    } else {					// a set of other items?
	get_set(istr, tag);			// access set's contents
	put_set(ostr, tag);			// output set token
	while ((tag1 = next_item_tag(istr)) != NULL) {
						// loop over set contents
	    copy_item(ostr, istr, tag1);	// copying each item
	    free(tag1);				// free up name of item
	}
	get_tes(istr, tag);			// close access to set
	put_tes(ostr, tag);			// output termial symbol
    }
    free(type);					// free up type string
    end_iofunc();
}

//
//  USER OUTPUT FUNCTIONS
//
 
//  __________________________________________
//  put_set: begin named set in output stream.
 
void put_set(stream ostr, string tag)
{
    stateptr ostate;
    itemptr oitem;
 
    set_iofunc("put_set");                      // set name for messages
    ostate = findstate(ostr);                   // find corresponding state
    oitem = makeitem(SetType, tag, NULL, NULL); // create item with tag
    putitem(ostate, oitem, NULL);               // do actual output of set
    bufitem(ostate, oitem);                     // and add to state buffer
    end_iofunc();                               // end this activation
}  
 
//  __________________________________
//  put_tes: end set in output stream.
 
void put_tes(stream ostr, string tag)
{
    stateptr ostate;
    int setind, i;
    itemptr oitem;
 
    set_iofunc("put_tes");
    ostate = findstate(ostr);                   // find corresponding state
    if (ostate->itemcnt == 0)
	error("%s.%s: called before put_set\n", getprog(), iofunc);
    setind = findset(ostate, ostate->itemcnt);  // find opening set item
    if (! streq(ostate->buf[setind]->tag, tag))
	error("%s.%s: closing tag %s does not match opening tag %s\n",
	      getprog(), iofunc, tag, ostate->buf[setind]->tag);
    oitem = makeitem(TesType, NULL, NULL, NULL);
    putitem(ostate, oitem, NULL);               // do actual output of tes
    bufitem(ostate, oitem);
    if (setind == 0) {                          // end of toplevel set?
	for (i = 0; i < ostate->itemcnt; i++)
	    freeitem(ostate->buf[i]);           // done with these items
	ostate->itemcnt = 0;
    }
    end_iofunc();
}

//  ______________________________________________
//  put_string: write string to a structured file.
 
void put_string(stream ostr, string tag, string dat)
{
    set_iofunc("put_string");
    put_data(ostr, tag, CharType, dat, cstrlen(dat), 0);
    end_iofunc();
}
 
//  _________________________________________________
//  put_data: write data object to a structured file.
 
void put_data(stream ostr, string tag, string typ, void *dat, ...)
{
    va_list ap;
    int dims[MaxVecDim + 1], ndim;
 
    set_iofunc("put_data");
    va_start(ap, dat);                          // access argument list
    ndim = 0;
    do {                                        // loop reading dimensions
	if (ndim > MaxVecDim)
	    error("%s.%s: item %s: too many dimensions\n",
		  getprog(), iofunc, tag);
	dims[ndim] = va_arg(ap, int);
    } while (dims[ndim++] != 0);                // until a zero comes up
    va_end(ap);
    put_data_sub(ostr, tag, typ, ndim > 1 ? dims : NULL, dat, NULL);
						// pass dims, if any
    end_iofunc();
}

//  _______________________________________________________________
//  put_data_masked: write masked data object to a structured file.
 
void put_data_masked(stream ostr, string tag, string typ, void *dat, ...)
{
    va_list ap;
    int dims[MaxVecDim + 1], ndim, *mask;
 
    set_iofunc("put_data_masked");
    va_start(ap, dat);                          // access argument list
    ndim = 0;
    do {                                        // loop getting dimensions
	if (ndim > MaxVecDim)
	    error("%s.%s: item %s: too many dimensions\n",
		  getprog(), iofunc, tag);
	dims[ndim] = va_arg(ap, int);
    } while (dims[ndim++] != 0);                // until a zero comes up
    mask = va_arg(ap, int *);                   // get mask following dims
    va_end(ap);
    put_data_sub(ostr, tag, typ, ndim > 1 ? dims : NULL, dat, mask);
    end_iofunc();
}
 
//  __________________________________________________________________________
//  put_data_sub: routine to handle output operation, using arguments gathered
//  by put_data and put_data_masked, or supplied directly by the user.

 
void put_data_sub(stream ostr, string tag, string type,
		  int *dims, void *data, int *mask)
{
    stateptr ostate;
    itemptr oitem;
 
    set_iofunc("put_data_sub");
    ostate = findstate(ostr);                   // find corresponding state
    oitem = makeitem(type, tag, dims, data);
    putitem(ostate, oitem, mask);               // do actual data output
    bufitem(ostate, oitem);
    if (ostate->itemcnt == 1) {                 // free toplevel items
	freeitem(oitem);
	ostate->itemcnt = 0;
    }
    end_iofunc();
}

//
//  USER INPUT FUNCTIONS
//
 
//  _____________________________________________________
//  get_set: switch context to named set in input stream.
 
void get_set(stream istr, string tag)
{
    stateptr istate;
    int setind;
 
    set_iofunc("get_set");
    istate = findstate(istr);                   // find corresponding state
    if (istate->itemcnt == 0 && ! readset(istate))
	error("%s.%s: at end of file\n", getprog(), iofunc);
    setind = finditem(istate, tag);             // get index of named set
    if (setind == -1)
	error("%s.%s: item %s not found\n", getprog(), iofunc, tag);
    if (! streq(istate->buf[setind]->type, SetType))
	error("%s.%s: %s not a set\n", getprog(), iofunc, tag);
    istate->context = setind + 1;               // point to 1st item in set
    end_iofunc();
}
 
//  _____________________________________________________
//  get_tes: switch context to set enclosing current one.
 
void get_tes(stream istr, string tag)
{
    stateptr istate;
    int setind, tesind, i;
 
    set_iofunc("get_tes");
    istate = findstate(istr);                   // find corresponding state
    if (istate->context == 0)
	error("%s.%s: called before get_set\n", getprog(), iofunc);
    setind = findset(istate, istate->context);  // find start of this set
    if (! streq(istate->buf[setind]->tag, tag))
	error("%s.%s: closing tag %s does not match opening tag %s\n",
	      getprog(), iofunc, tag, istate->buf[setind]->tag);
    if (setind > 0) {                           // at end of subset?
	tesind = findtes(istate, istate->context);
	istate->context = tesind + 1;           // point past closing tes
    } else {                                    // at end of toplevel set?
	for (i = 0; i < istate->itemcnt; i++) {
	    if (istate->buf[i]->data != NULL)   // free data and item
		free(istate->buf[i]->data);
	    freeitem(istate->buf[i]);
	}
	istate->context = istate->itemcnt = 0;
    }
    end_iofunc();
}

//  _________________________________________________
//  get_string: read a string from a structured file.
 
string get_string(stream istr, string tag)
{
    int dlen;
    string dat;
 
    set_iofunc("get_string");
    dlen = get_length(istr, tag);
    dat = allocate(dlen + 1);                   // add 1 byte for null
    get_data(istr, tag, CharType, dat, dlen, 0);
    dat[dlen] = (char) NULL;                    // set terminating byte
    end_iofunc();
    return (dat);
}
 
//  __________________________________________________
//  get_data: read data object from a structured file.
 
void get_data(stream istr, string tag, string typ, void *dat, ...)
{
    va_list ap;
    int dims[MaxVecDim + 1], ndim;
 
    set_iofunc("get_data");
    va_start(ap, dat);                          // access argument list
    ndim = 0;
    do {                                        // loop reading dimensions
	if (ndim > MaxVecDim)
	    error("%s.%s: too many dimensions; item %s\n",
		  getprog(), iofunc, tag);
	dims[ndim] = va_arg(ap, int);
    } while (dims[ndim++] != 0);                // until a zero comes up
    va_end(ap);
    get_data_sub(istr, tag, typ, ndim > 1 ? dims : NULL, dat, NULL);
						// pass dims, if any
    end_iofunc();
}

//  _________________________________________________________
//  get_data_masked: read data object from a structured file.
 
void get_data_masked(stream istr, string tag, string typ, void *dat, ...)
{
    va_list ap;
    int dims[MaxVecDim + 1], ndim, *mask;
 
    set_iofunc("get_data_masked");
    va_start(ap, dat);                          // access argument list
    ndim = 0;
    do {                                        // loop getting dimensions
	if (ndim > MaxVecDim)
	    error("%s.%s: too many dimensions; item %s\n",
		  getprog(), iofunc, tag);
	dims[ndim] = va_arg(ap, int);
    } while (dims[ndim++] != 0);                // until a zero comes up
    mask = va_arg(ap, int *);                   // get mask following dims
    va_end(ap);
    get_data_sub(istr, tag, typ, ndim > 1 ? dims : NULL, dat, mask);
    end_iofunc();
}
 
//  _________________________________________________________________________
//  get_data_sub: routine to handle input operation, using arguments gathered
//  by get_data and get_data_masked, or supplied directly by the user.
 
void get_data_sub(stream istr, string tag, string type,
		  int *dims, void *data, int *mask)
{
    stateptr istate;
    int indx;
    itemptr ditem;
 
    set_iofunc("get_data_sub");
    istate = findstate(istr);                   // find corresponding state
    if (istate->itemcnt == 0 && ! readset(istate))
	error("%s.%s: at end of file\n", getprog(), iofunc);
    indx = finditem(istate, tag);               // get index of named item
    if (indx == -1)
	error("%s.%s: %s not found\n", getprog(), iofunc, tag);
    ditem = makeitem(type, tag, dims, data);    // make destination item
    copydata(ditem, mask, istate->buf[indx], istr);
						// transfer actual data
    freeitem(ditem);
    if (istate->context > 0)
	istate->context = indx + 1;             // advance context in set
    else {
	free(istate->buf[0]->data);             // free toplevel item
	freeitem(istate->buf[0]);
	istate->context = istate->itemcnt = 0;
    }
    end_iofunc();
}

//  _____________________________________________________________________
//  get_tag_ok: return TRUE if tag exists in current context, else FALSE.
 
bool get_tag_ok(stream istr, string tag)
{
    stateptr istate;
    int indx;
 
    set_iofunc("get_tag_ok");
    istate = findstate(istr);                   // find corresponding state
    indx = (istate->itemcnt == 0 && ! readset(istate) ?
	      -1 : finditem(istate, tag));      // scan context for tag
    end_iofunc();
    return (indx != -1);
}
 
//  ________________________________________
//  get_type: return type of specified item.
 
string get_type(stream istr, string tag)
{
    stateptr istate;
    int indx;
 
    set_iofunc("get_type");
    istate = findstate(istr);
    if (istate->itemcnt == 0 && ! readset(istate))
	error("%s.%s: at end of file\n", getprog(), iofunc);
    indx = finditem(istate, tag);
    if (indx == -1)
	error("%s.%s: %s not found\n", getprog(), iofunc, tag);
    end_iofunc();
    return (copycstr(istate->buf[indx]->type)); // return copy of item type
}

//  ___________________________________________________________
//  get_dimensions: return dimension list of specified item, or
//  NULL if the item is a scalar.
 
int *get_dimensions(stream istr, string tag)
{
    stateptr istate;
    int indx;
 
    set_iofunc("get_dimensions");
    istate = findstate(istr);
    if (istate->itemcnt == 0 && ! readset(istate))
	error("%s.%s: at end of file\n", getprog(), iofunc);
    indx = finditem(istate, tag);
    if (indx == -1)
	error("%s.%s: %s not found\n", getprog(), iofunc, tag);
    end_iofunc();
    return (istate->buf[indx]->dims != NULL ?
	      copyistr(istate->buf[indx]->dims) : NULL);
						// return copy of item dims
}
 
//  _____________________________________________________________
//  get_length: return length in bytes of data of specified item.
 
int get_length(stream istr, string tag)
{
    stateptr istate;
    int indx, dlen;
 
    set_iofunc("get_length");
    istate = findstate(istr);
    if (istate->itemcnt == 0 && ! readset(istate))
	error("%s.%s: at end of file\n", getprog(), iofunc);
    indx = finditem(istate, tag);
    if (indx == -1)
	error("%s.%s: %s not found\n", getprog(), iofunc, tag);
    dlen = datalength(istate->buf[indx]);
    end_iofunc();
    return (dlen);
}

//  __________________________________________________
//  next_item_tag: return tag of next item in context,
//  or NULL at end of set or file.
 
string next_item_tag(stream istr)
{
    stateptr istate;
    string tag;
 
    set_iofunc("next_item_tag");
    istate = findstate(istr);                   // find corresponding state
    if (istate->itemcnt == 0 && ! readset(istate))
	tag = NULL;
    else
	tag = istate->buf[istate->context]->tag;
    end_iofunc();
    return (tag != NULL ? copycstr(tag) : NULL);
						// return copy of tag
}
 
//  _________________________________________________________________
//  skip_item: advance context within set, or flush item at toplevel.
//  Returns TRUE on success, FALSE on end of file.
 
bool skip_item(stream istr)
{
    stateptr istate;
    int i;
 
    set_iofunc("skip_item");
    istate = findstate(istr);                   // find corresponding state
    if (istate->itemcnt == 0 && ! readset(istate)) {
	end_iofunc();
	return (FALSE);
    }
    if (istate->context > 0) {                  // skip item within set
	if (streq(istate->buf[istate->context]->type, TesType))
	    istate->context = findset(istate, istate->context);
	else if (streq(istate->buf[istate->context]->type, SetType))
	    istate->context = findtes(istate, istate->context + 1);
	istate->context++;
    } else {                                    // skip item at toplevel
	for (i = 0; i < istate->itemcnt; i++) {
	    if (istate->buf[i]->data != NULL)   // free data and item
		free(istate->buf[i]->data);
	    freeitem(istate->buf[i]);
	}
	istate->context = istate->itemcnt = 0;
    }
    end_iofunc();
    return (TRUE);
}

//
//  ASSORTED CONTROL FUNCTIONS
//
 
//  ____________________________________________
//  fs_options: select input conversion options.
 
local int conformdata = WarnEach;	// warn about data conformence	
local int convertdata = WarnEach;	// warn about data conversion	
 
void fs_options(int cfd, int cnd)
{
    if (cfd != NoChange)			// change data-conformance
	conformdata = cfd;
    if (cnd != NoChange)			// change type-conversion
	convertdata = cnd;
}
 
//  __________________________________________________________
//  strclose: close stream and remove from stream state table.
 
void strclose(stream str)
{
    stateptr iostate;
 
    fclose(str);
    iostate = findstate(str);   
    iostate->str = NULL;
    free(iostate->buf);
}    

//
//  LOW LEVEL ROUTINES
//
 
//  __________________________________________________________________
//  findstate: find stream in state table, and return state structure.
 
local state strtab[MaxStream] = { NULL, };
 
local stateptr findstate(stream str)
{
    int i;
    struct stat statbuf;
 
    for (i = 0; i < MaxStream; i++)
	if (strtab[i].str == str)
	    return (&strtab[i]);
    for (i = 0; i < MaxStream && strtab[i].str != NULL; i++) ;
    if (i == MaxStream)
	error("%s.%s: no room in stream table\n", getprog(), iofunc);
    strtab[i].str = str;
    strtab[i].buf = (itemptr *) allocate(sizeof(itemptr) * BufSizeInc);
    strtab[i].maxbuf = BufSizeInc;
    strtab[i].itemcnt = strtab[i].context = 0;
    if (fstat(fileno(str), &statbuf) == -1)
	error("%s.%s: can't get status of stream %d\n",
	      getprog(), iofunc, fileno(str));
    strtab[i].canseek = (statbuf.st_mode & S_IFMT) == S_IFREG;
    return (&strtab[i]);
}
 
//  ____________________________________
//  bufitem: store item in state buffer.
 
local void bufitem(stateptr iostate, itemptr newitem)
{
    int newmax, i;
    itemptr *newbuf;
 
    if (iostate->itemcnt == iostate->maxbuf) {  // need to enlarge buffer?
	newmax = iostate->maxbuf + BufSizeInc;  // enlarge by BufSizeInc
	if ((newmax & iostate->maxbuf) == 0)
	    eprintf("[%s.%s: enlarging item buffer]\n", getprog(), iofunc);
	newbuf = (itemptr *) allocate(sizeof(itemptr) * newmax);
	for (i = 0; i < iostate->itemcnt; i++)  // copy old one over
	    newbuf[i] = iostate->buf[i];
	free(iostate->buf);
	iostate->buf = newbuf;                  // install new buffer
	iostate->maxbuf = newmax;
    }
    iostate->buf[iostate->itemcnt++] = newitem; // insert item in buffer
}

//  ______________________________________________________________________
//  finditem: scan current context for named item; return index in buffer,
//  or -1 if not found.
 
local int finditem(stateptr iostate, string tag)
{
    int level, i;
 
    level = 0;
    for (i = iostate->context; i < iostate->itemcnt && level >= 0; i++) {
	if (level == 0 && iostate->buf[i]->tag != NULL &&
	      streq(iostate->buf[i]->tag, tag))
	    return (i);
	if (streq(iostate->buf[i]->type, SetType))
	    level++;
	else if (streq(iostate->buf[i]->type, TesType))
	    level--;
    }
    level = 0;
    for (i = iostate->context - 1; i >= 0 && level >= 0; i--) {
	if (streq(iostate->buf[i]->type, TesType))
	    level++;
	else if (streq(iostate->buf[i]->type, SetType))
	    level--;
	if (level == 0 && iostate->buf[i]->tag != NULL &&
	      streq(iostate->buf[i]->tag, tag))
	    return (i);
    }
    return (-1);
}
 
//  ________________________________________________________
//  findset: scan back through buffer to find enclosing set.
 
local int findset(stateptr iostate, int cont)
{
    int level, i;
 
    level = 0;
    for (i = cont - 1; i >= 0; i--) {
	if (streq(iostate->buf[i]->type, TesType))
	    level++;
	else if (streq(iostate->buf[i]->type, SetType))
	    level--;
	if (level < 0)
	    return (i);
    }
    error("%s.%s: can't find set\n", getprog(), iofunc);
    return (0);					// keep compiler happy...
}

//  ___________________________________________________________
//  findtes: scan forward through buffer to find enclosing tes.
 
local int findtes(stateptr iostate, int cont)
{
    int level, i;
 
    level = 0;
    for (i = cont; i < iostate->itemcnt; i++) {
	if (streq(iostate->buf[i]->type, SetType))
	    level++;
	else if (streq(iostate->buf[i]->type, TesType))
	    level--;
	if (level < 0)
	    return (i);
    }
    error("%s.%s: can't find tes\n", getprog(), iofunc);
    return (0);					// keep compiler happy...
}
 
//  ________________________________________________________
//  putitem: write header and data of item to output stream.
 
local void putitem(stateptr ostate, itemptr oitem, int *mask)
{
    short magic;
    int obyte, *mp, n;
    byte *op;
 
    magic = (oitem->dims != NULL ? PluralMagic : ScalarMagic);
    safewrite(ostate->str, &magic, sizeof(short));
    safewrite(ostate->str, oitem->type, cstrlen(oitem->type));
    if (oitem->tag != NULL)
	safewrite(ostate->str, oitem->tag, cstrlen(oitem->tag));
    if (oitem->dims != NULL)
	safewrite(ostate->str, oitem->dims,
		  sizeof(int) * istrlen(oitem->dims));
    if (oitem->data != NULL) {			// have data to output?
	if (mask == NULL)			// if possible, just write
	    safewrite(ostate->str, oitem->data, datalength(oitem));
	else {					// output with masking
	    obyte = datalength(oitem);		// count bytes to output
	    op = (byte *) oitem->data;		// start at front of data
	    mp = mask;				// and at front of mask
	    while (obyte > 0) {			// loop while bytes remain
		n = *mp++;			// get output/skip count
		if (n < 0)			// should skip some bytes?
		    op -= n;			// just advance pointer
		else {				// should output some bytes
		    safewrite(ostate->str, op, n);
		    op += n;			// advance output pointer
		    obyte -= n;			// count bytes done
		}
		mp = (*mp == 0 ? mask : mp);	// renew mask pointer
	    }
	    if (obyte < 0)			// this should not happen
		error("%s.%s: item %s: incommensurate mask\n",
		      getprog(), iofunc, oitem->tag);
	}
    }
}

//  _________________________________________________________
//  readset: read an item or entire set into an empty buffer.
 
local bool readset(stateptr istate)
{
    int level;
    itemptr iitem;
 
    assert(istate->itemcnt == 0);
    level = 0;
    do {
	iitem = getitem(istate);
	if (iitem == NULL && level == 0)
	    return (FALSE);
	if (iitem == NULL && level > 0)
	    error("%s.%s: EOF within set at level %d\n",
		  getprog(), iofunc, level);
	bufitem(istate, iitem);
	if (streq(iitem->type, SetType))
	    level++;
	else if (streq(iitem->type, TesType))
	    level--;
    } while (level > 0);
    istate->context = 0;
    return (TRUE);
}       

//  __________________________________________________________________
//  getitem: read item from input, and return pointer, or NULL on EOF.
 
local itemptr getitem(stateptr istate)
{
    short magic, cigam;
    string type, tag;
    int *dims, *ip;
    itemptr iitem;
 
    if (fread(&magic, sizeof(short), 1, istate->str) != 1)
	return (NULL);				// nothing more to read	
    cigam = ((magic & 0xff) << 8) | ((magic & 0xff00) >> 8);
						// form swapped version
    if (magic != ScalarMagic && cigam != ScalarMagic &&
	magic != PluralMagic && cigam != PluralMagic)
	error("%s.%s: bad magic number %04x\n",
	      getprog(), iofunc, 0xffff & magic);
    if (cigam == ScalarMagic || cigam == PluralMagic)
	eprintf("[%s.%s: WARNING: reading swapped data]\n",
		getprog(), iofunc);
    type = getcstr(istate->str);		// read item components	
    tag = (streq(type, TesType) ? NULL : getcstr(istate->str));
    dims = (magic == PluralMagic || cigam == PluralMagic ?
	    getistr(istate->str) : NULL);
    if (cigam == PluralMagic)
      for (ip = dims; *ip != (int) NULL; ip++)
	    *ip = ((*ip & 0xff) << 24) | ((*ip & 0xff00) << 8) |
		  ((*ip & 0xff0000) >> 8) | ((*ip & 0xff000000) >> 24);
    iitem = makeitem(type, tag, dims, NULL);	// construct actual item
    if (! streq(type, SetType) && ! streq(type, TesType))
	getdata(iitem, istate);			// read (or skip) item data
    free(type);					// free temporary storage
    if (tag != NULL)
	free(tag);
    if (dims != NULL)
	free(dims);
    return (iitem);
}
 
//  _____________________________________________________________
//  getdata: read data for item, or record place in input stream.
 
local void getdata(itemptr iitem, stateptr istate)
{
    int nb;
 
    nb = datalength(iitem);
    if (nb > MaxDataBuf && istate->canseek) {	// skip data input for now?
	iitem->data = NULL;
	iitem->datapos = ftell(istate->str);	// remember this place	
	safeseek(istate->str, nb, SEEK_CUR);	// and seek past data
    } else {
	iitem->datapos = 0;
	iitem->data = allocate(nb);		// allocate data buffer	
	saferead(istate->str, iitem->data, nb);	// actually read the data
    }
}

//  _____________________________________________________________
//  copydata: transfer data from source item to destination item.
//  The source data may reside within memory, or in a disk file.
//  Interconversion of float and double takes place if needed.
//  Multidimensional arrays may be truncated/collapsed to vectors.
//  If a mask string is given, the destination bytes are masked.
 
#define NoWayBaby     0			// types cannot be interconverted
#define DirectCopy    1			// no interconversion is needed	
#define Float2Double  2			// source is float, dest is double
#define Double2Float  3			// source is doubel, dest is float
 
local void copydata(itemptr ditem, int *mask, itemptr sitem, stream sstr)
{
    int conv, dbyte, sbyte, defmask[2];
    long currpos;
 
    conv = convmode(ditem, sitem);		// set data conversion	
    if (conv == NoWayBaby)
	error("%s.%s: item %s: cannot convert type\n",
	      getprog(), iofunc, ditem->tag);
    if (! checkdims(ditem, sitem))              // check data congruency
	error("%s.%s: item %s: incongruent dimensions\n",
	      getprog(), iofunc, ditem->tag);
    if (sitem->data == NULL) {			// data not yet in core?
	currpos = ftell(sstr);
	safeseek(sstr, sitem->datapos, SEEK_SET);
    }
    if (conv == DirectCopy && mask == NULL) {   // if possible, just copy
	inputdata(ditem->data, sitem, sstr, 0, datalength(ditem));
    } else {                                    // transfer with mapping
	dbyte = datalength(ditem);		// count bytes to transfer
	sbyte = datalength(sitem);
	if (mask == NULL) {			// no masking required?
	    mask = defmask;			// setup a default mask
	    defmask[0] = dbyte;			// covering entire dest'n
	    defmask[1] = 0;
	}
	if (type_length(ditem->type) % sizeof(int) == 0)
	    copyint(ditem, mask, sitem, sstr, conv, dbyte, sbyte);
	else
	    copybyte(ditem, mask, sitem, sstr, conv, dbyte, sbyte);
    }
    if (sitem->data == NULL)			// data just read in?
	safeseek(sstr, currpos, SEEK_SET);	// then restore stream pos
}

local void copyint(itemptr ditem, int *mask, itemptr sitem, stream sstr,
		   int conv, int dbyte, int sbyte)
{
    int sskip, xbyte, *mp, ncopy, n, *dp, *xp, *dstop;

    sskip = xbyte = 0;                     	// zero offset and buffer
    dp = (int *) ditem->data;              	// start at front of data
    mp = mask;					// and at front of mask
    while (dbyte > 0) {				// loop while bytes remain
	if (*mp < 0)				// should skip some bytes?
	    dp -= *mp++ / sizeof(int);		// just advance pointer
	else {					// should copy some bytes
	    ncopy = *mp++;			// set number to copy
	    dbyte -= ncopy;			// count that many done	
	    if (dbyte < 0)			// this should not happen
		error("%s.%s: item %s: incommensurate mask\n",
		      getprog(), iofunc, ditem->tag);
	    do {				// loop filling xfer buffer
		if (xbyte == 0)			// is the buffer empty?
		    fillxbuf((byte **) &xp, &xbyte, &sskip, &sbyte,
			     sitem, sstr, conv);
		n = MIN(ncopy, xbyte);		// find max ready to copy
		dstop = dp + n / sizeof(int);	// fix stopping point
		while (dp < dstop)		// loop copying ints...
		    *dp++ = *xp++;		// one by one they go
		xbyte -= n;			// count bytes in buffer
		ncopy -= n;			// and bytes still to copy
	    } while (ncopy > 0);		// loop back if more to do
	}
	mp = (*mp == 0 ? mask : mp);		// renew mask pointer
    }
}

local void copybyte(itemptr ditem, int *mask, itemptr sitem, stream sstr,
		    int conv, int dbyte, int sbyte)
{
    int sskip, xbyte, *mp, ncopy, n;
    byte *dp, *xp, *dstop;

    sskip = xbyte = 0;                     	// zero offset and buffer
    dp = (byte *) ditem->data;              	// start at front of data
    mp = mask;					// and at front of mask
    while (dbyte > 0) {				// loop while bytes remain
	if (*mp < 0)				// should skip some bytes?
	    dp -= *mp++;			// just advance pointer
	else {					// should copy some bytes
	    ncopy = *mp++;			// set number to copy
	    dbyte -= ncopy;			// count that many done
	    if (dbyte < 0)			// this should not happen
		error("%s.%s: item %s: incommensurate mask\n",
		      getprog(), iofunc, ditem->tag);
	    do {				// loop filling xfer buffer
		if (xbyte == 0)			// is the buffer empty?
		    fillxbuf(&xp, &xbyte, &sskip, &sbyte, sitem, sstr, conv);
		n = MIN(ncopy, xbyte);		// find max ready to copy
		dstop = dp + n;			// fix stopping point
		while (dp < dstop)		// loop copying bytes...
		    *dp++ = *xp++;		// one at a time (sigh)
		xbyte -= n;			// count bytes in buffer
		ncopy -= n;			// and bytes still to copy
	    } while (ncopy > 0);		// loop back if more to do
	}
	mp = (*mp == 0 ? mask : mp);		// renew mask pointer
    }
}

//  _______________________________________________________________________
//  fillxbuf: fill transfer buffer, performing float <-> double conversion.
 
local void fillxbuf(byte **xpp, int *xbytep, int *sskipp, int *sbytep,
		    itemptr sitem, stream sstr, int conv)
{
    static byte xbuf[MaxDataBuf];
    static float fbuf[MaxDataBuf / sizeof(float)];
    static double dbuf[MaxDataBuf / sizeof(float)];
    int scount, n;
 
    switch (conv) {
      case DirectCopy:                          // no conversion required
	scount = MIN(*sbytep, sizeof(xbuf));
	inputdata(xbuf, sitem, sstr, *sskipp, scount);
	*xpp = xbuf;
	*xbytep = scount;                   
	break;
      case Float2Double:                        // convert float -> double
	scount = MIN(*sbytep, sizeof(fbuf));
	inputdata(fbuf, sitem, sstr, *sskipp, scount);
	n = scount / sizeof(float);             // number of floats read
	while (--n >= 0)
	   dbuf[n] = fbuf[n];                   // copy and convert them
	*xpp = (byte *) dbuf;
	*xbytep = (scount * sizeof(double)) / sizeof(float);
	break;
      case Double2Float:                        // convert double -> float
	scount = MIN(*sbytep, sizeof(dbuf));
	inputdata(dbuf, sitem, sstr, *sskipp, scount);
	n = scount / sizeof(double);
	while (--n >= 0)
	    fbuf[n] = dbuf[n];
	*xpp = (byte *) fbuf;
	*xbytep = (scount * sizeof(float)) / sizeof(double);
	break;
      default:
	error("%s.%s: bad conversion\n", getprog(), iofunc);
    }
    *sskipp += scount;
    *sbytep -= scount;
}

//  _________________________________________________________________
//  inputdata: transfer scount bytes from source item to destination.
 
local void inputdata(void *data, itemptr sitem, stream sstr,
		     int sskip, int scount)
{
    byte *dp = (byte *) data, *sp = (byte *) sitem->data + sskip;
 
    if (sitem->data != NULL)
	while (--scount >= 0)
	    *dp++ = *sp++;
    else
	saferead(sstr, data, scount);
}

//  ___________________________________________________
//  convmode: determine conversion mode for data input.
 
local int convmode(itemptr ditem, itemptr sitem)
{
    string dtype = ditem->type, stype = sitem->type;
 
    if (streq(stype, dtype))			// types match exactly?
	return (DirectCopy);
    if (convertdata != NotAllow) {		// data coercion allowed?
	if (streq(stype, FloatType) && streq(dtype, DoubleType)) {
	    if (convertdata != SilentOK)
		eprintf("[%s.%s: converting %s to double]\n",
			getprog(), iofunc, ditem->tag);
	    if (convertdata == WarnOnce)
		convertdata = SilentOK;
	    return (Float2Double);
	}
	if (streq(stype, DoubleType) && streq(dtype, FloatType)) {
	    if (convertdata != SilentOK)
		eprintf("[%s.%s: converting %s to float]\n",
			getprog(), iofunc, ditem->tag);
	    if (convertdata == WarnOnce)
		convertdata = SilentOK;
	    return (Double2Float);
	}
    }
    return (NoWayBaby);				// no way to coerce type
}

//  ________________________________________________________________
//  checkdims: return TRUE if item dimensions can be made congruent.
 
local bool checkdims(itemptr ditem, itemptr sitem)
{
    int *ddims = ditem->dims, *sdims = sitem->dims;
 
    if (ddims == NULL && sdims == NULL)
	return (TRUE);				// both items are scalars
    if (ddims != NULL && sdims != NULL && istreq(ddims, sdims))
	return (TRUE);				// dimensions match exactly
    if (conformdata != NotAllow && ddims != NULL && istrlen(ddims) == 2 &&
	  ddims[0] <= datacount(sitem)) {	// can conform input data?
	if (conformdata != SilentOK)
	    eprintf("[%s.%s: WARNING: conforming item %s]\n",
		    getprog(), iofunc, ditem->tag);
	if (conformdata == WarnOnce)
	    conformdata = SilentOK;
	return (TRUE);				// allow data conformance
    }
    return (FALSE);				// data are not congruent
}

//  ________________________________________________________________________
//  datalength, datacount: compute length in bytes or elements of item data.
 
local int datalength(itemptr ioitem)
{
    return (datacount(ioitem) * type_length(ioitem->type));
}
 
local int datacount(itemptr ioitem)
{
    int nelt, *ip;
 
    nelt = 1;
    if (ioitem->dims != NULL)
      for (ip = ioitem->dims; *ip != (int) NULL; ip++)
	    nelt *= *ip;
    return (nelt);
}

//  _____________________________________________________________________
//  makeitem: allocate an item, and initialize fields.  The type, tag and
//  dimension arguments are copied; the data is not, since it may be big.
 
local itemptr makeitem(string type, string tag, int *dims, void *data)
{
    itemptr ioitem;
 
    ioitem = (itemptr) allocate(sizeof(item));
    ioitem->type = copycstr(type);
    ioitem->tag = (tag != NULL ? copycstr(tag) : NULL);
    ioitem->dims = (dims != NULL ? copyistr(dims) : NULL);
    ioitem->data = data;
    ioitem->datapos = 0;
    return (ioitem);
}
 
//  ______________________________________________________________________
//  freeitem: deallocate an item, including all fields copied by makeitem.
 
local void freeitem(itemptr ioitem)
{
    free(ioitem->type);
    if (ioitem->tag != NULL)
	free(ioitem->tag);
    if (ioitem->dims != NULL)
	free(ioitem->dims);
    free(ioitem);
}
 
//  ___________________________________________________________________
//  safewrite, saferead, safeseek: I/O operations, with error checking.
 
local void safewrite(stream str, void *data, int len)
{
    if (fwrite(data, sizeof(byte), len, str) != len)
	error("%s.%s: write failed (%d bytes)\n", getprog(), iofunc, len);
}
 
local void saferead(stream str, void *data, int len)
{
    if (fread(data, sizeof(byte), len, str) != len)
	error("%s.%s: read failed (%d bytes)\n", getprog(), iofunc, len);
}
 
local void safeseek(stream str, long offset, int key)
{
    if (fseek(str, offset, key) == -1)
	error("%s.%s: fseek failed (%lx bytes)\n", getprog(), iofunc, offset);
}
