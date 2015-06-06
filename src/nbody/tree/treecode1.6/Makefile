# zeno/src/nbody/tree/treecode1.6/Makefile: compile tree-code program.
# Copyright (c) 2015, Joshua E. Barnes, Kyoto, Japan.

BINFILES = treecode treecode_x treepeek

# Optimization flags; use second version for SunOS.
OPT = -O3
# OPT = -xO3

########################################################################
# Build vanilla treecode; by default, includes softening correction.

OPTIONS = -DSOFTCORR

treecode$X: treecode$X.o treeio$X.o treeload$X.o treegrav$X.o
	$(ZCC) $(ZLDFLAGS) -o treecode$X \
	    treecode$X.o treeio$X.o treeload$X.o treegrav$X.o \
	    -lNBody -lClib -lgsl -lgslcblas -lm

treecode$X.o: treecode.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -o treecode$X.o -c treecode.c

treeio$X.o: treeio.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -o treeio$X.o -c treeio.c

treeload$X.o: treeload.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -o treeload$X.o -c treeload.c

treegrav$X.o: treegrav.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -o treegrav$X.o -c treegrav.c

########################################################################
# Build treecode with external force.

treecode_x: treecode_x.o treeio.o treeload.o treegrav.o
	$(ZCC) $(ZLDFLAGS) -o treecode_x \
	    treecode_x.o treeio.o treeload.o treegrav.o \
	    -lGSP -lNBody -lClib -lgsl -lgslcblas -lm

treecode_x.o: treecode.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) -DEXTGRAV -o treecode_x.o -c treecode.c

########################################################################
# Generate listing.

treecode.ps: treecode.c treeio.c treeload.c treegrav.c treecode.h treedefs.h
	enscript -M Letterdj -o treecode.ps treecode.c treecode.h treeio.c \
	         treeload.c treegrav.c treedefs.h

########################################################################
# Build treepeek utility.

treepeek: treepeek.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o treepeek treepeek.c -lNBody -lClib
