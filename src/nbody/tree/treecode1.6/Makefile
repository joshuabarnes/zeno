##  zeno/src/nbody/tree/treecode1.6/Makefile: compile tree-code program.
##  Copyright (c) 2017, Joshua E. Barnes, Honolulu, Hawaii.
##  ____________________________________________________________________

BINFILES = treecode$X treecode_x$X treepeek$X

##  Compilation options are selected by changing OPTIONS (and X):
##	make make_all OPTIONS="-DDOUBLEPREC" X=_dp
##  compiles BINFILES in double precision, with "_dp" suffixes.
##  _____________________________________________________________

OPTIONS =

##  Optimization flags.
##  ___________________

OPT = -O3

make_all: $(BINFILES)

##  Build vanilla treecode; by default, includes softening correction.
##  __________________________________________________________________

treecode$Y$X: treecode$X.o treeio$X.o treeload$X.o treegrav$Y$X.o
	$(ZCC) $(ZLDFLAGS) -o treecode$Y$X \
	    treecode$X.o treeio$X.o treeload$X.o treegrav$Y$X.o \
	    -lNBody -lClib -lgsl -lgslcblas -lm

treecode$X.o: treecode.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -o treecode$X.o -c treecode.c

treeio$X.o: treeio.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -o treeio$X.o -c treeio.c

treeload$X.o: treeload.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -o treeload$X.o -c treeload.c

treegrav$Y$X.o: treegrav.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -o treegrav$Y$X.o -c treegrav$Y.c


##  Build treecode with external force.
##  ___________________________________

treecode_x$X: treecode_x$X.o treeio$X.o treeload$X.o treegrav$X.o
	$(ZCC) $(ZLDFLAGS) -o treecode_x$X \
	    treecode_x$X.o treeio$X.o treeload$X.o treegrav$X.o \
	    -lGSP -lNBody -lClib -lgsl -lgslcblas -lm

treecode_x$X.o: treecode.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -DEXTGRAV -o treecode_x$X.o -c treecode.c

##  Generate listing.
##  _________________

zeno_treecode.ps: treecode.c treeio.c treeload.c treegrav.c \
	          treecode.h treedefs.h
	enscript -r2 -M Letterdj -o zeno_treecode.ps \
	         treecode.c treecode.h treeload.c treegrav.c \
	         treeio.c treedefs.h

########################################################################
# Build treepeek utility.

treepeek$X: treepeek.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o treepeek$X \
	       treepeek.c -lNBody -lClib
