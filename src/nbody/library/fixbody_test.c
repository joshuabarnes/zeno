
#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"

string defv[] = {
  "VERSION=1.0",			";Joshua Barnes  26 December 2014",
  NULL,
};

typedef struct {
  vector pos;
  real mass;
  byte type;
  int key;
  real auxarr[4];
} body, *bodyptr;

#define Pos(x)  ((x)->pos)
#define Mass(x)  ((x)->mass)
#define Type(x)  ((x)->type)
#define Key(x)  ((x)->key)
#define AuxArr(x)  ((x)->auxarr)

#include "fixbody.h"

int main(int argc, string argv[])
{
  int i;

  initparam(argv, defv);
  define_body(sizeof(body), Precision, NDIM);
  define_body_offset(PosTag,  BodyOffset(Pos));
  define_body_offset(MassTag, BodyOffset(Mass));
  define_body_offset(TypeTag, BodyOffset(Type));
  define_body_offset(KeyTag, BodyOffset(Key));
  define_body_offset(AuxArrTag "[4]", BodyOffset(AuxArr));
  printf("%s: length = %d\n", phatbody[0].name, phatbody[0].length);
  for (i = 1; phatbody[i].name != NULL; i++)
    if (phatbody[i].offset != BadOffset)
      printf("  %s: type = %s  offset = %d  length = %d\n", phatbody[i].name,
	     phatbody[i].type, phatbody[i].offset, phatbody[i].length);
}
