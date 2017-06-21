# ~/nbody/library/Makefile: build and test NBody library package.
# Copyright (c) 1997  Joshua E. Barnes  Honolulu, Hawai'i.

OBJFILES = bodyio_f3.o bodyio_d3.o buildmap.o phatbody.o\
	   snapcenter_f3.o snapcenter_d3.o snapkinem_f3.o snapkinem_d3.o

BINFILES = buildmap

# Assemble object library.

libNBody.a: $(OBJFILES)
	ar ruv libNBody.a $(OBJFILES)

tidy:
	rm -f $(OBJFILES) snapmap_s3.o snapmap_m3.o snapmap_d3.o

# Create individual object files.

bodyio_f3.o: bodyio.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DSINGLEPREC -DNDIM=3 -o bodyio_f3.o -c bodyio.c

bodyio_d3.o: bodyio.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DDOUBLEPREC -DNDIM=3 -o bodyio_d3.o -c bodyio.c

buildmap.o: buildmap.c buildmap.h
	$(ZCC) $(ZCCFLAGS) -c buildmap.c

phatbody.o: phatbody.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -c phatbody.c

snapcenter_f3.o: snapcenter.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DSINGLEPREC -DNDIM=3 -o snapcenter_f3.o \
	      -c snapcenter.c

snapcenter_d3.o: snapcenter.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DDOUBLEPREC -DNDIM=3 -o snapcenter_d3.o \
	      -c snapcenter.c

snapkinem_f3.o: snapkinem.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DSINGLEPREC -DNDIM=3 -o snapkinem_f3.o \
	      -c snapkinem.c

snapkinem_d3.o: snapkinem.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DDOUBLEPREC -DNDIM=3 -o snapkinem_d3.o \
	      -c snapkinem.c

# Main program object files for snapshot mapping.

snapmap_s3.o: snapmap.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DSINGLEPREC -DNDIM=3 -o snapmap_s3.o -c snapmap.c

snapmap_m3.o: snapmap.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DMIXEDPREC  -DNDIM=3 -o snapmap_m3.o -c snapmap.c

snapmap_d3.o: snapmap.c phatbody.h
	$(ZCC) $(ZCCFLAGS) -DDOUBLEPREC -DNDIM=3 -o snapmap_d3.o -c snapmap.c

# Utility program for snapshot mapping.

buildmap: buildmap.c buildmap.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED -o buildmap buildmap.c -lClib

# Test routines.

phatbody_test: phatbody.c phatbody.h libNBody.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED -o phatbody_test \
	      phatbody.c -lClib libNBody.a

fixbody_test: fixbody_test.c fixbody.h phatbody.o
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o fixbody_test \
	      fixbody_test.c phatbody.o -lClib

# Generate listing of source code.

zeno_library.ps: bodytags.h phatbody.h phatbody.c fixbody.h fixbody_test.c \
	         bodyio.c buildmap.h buildmap.c snapmap.c \
	         snapcenter.h snapcenter.c snapkinem.h snapkinem.c
	enscript -M Letterdj -r2 -o zeno_library.ps \
		 bodytags.h phatbody.h phatbody.c fixbody.h fixbody_test.c \
	         bodyio.c buildmap.h buildmap.c snapmap.c \
	         snapcenter.h snapcenter.c snapkinem.h snapkinem.c
