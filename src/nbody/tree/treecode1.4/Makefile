##############################################################################
# MAKEFILE: compile tree-code program.                                       #
# Copyright (c) 2001, Joshua E. Barnes, Honolulu, HI.                        #
##############################################################################

BINFILES = treecode treecode_q treepeek

# Optimization flags; use second version for SunOS.
OPT = -O3
# OPT = -xO3

########################################################################
# Build vanilla treecode.

treecode$X: treecode$X.o treeio$X.o treeload$X.o treegrav$X.o
	$(ZCC) $(ZLDFLAGS) -o treecode$X \
	    treecode$X.o treeio$X.o treeload$X.o treegrav$X.o \
	    -lNBody -lClib -lm

treecode$X.o: treecode.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -o treecode$X.o -c treecode.c

treeio$X.o: treeio.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) -o treeio$X.o -c treeio.c

treeload$X.o: treeload.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -o treeload$X.o -c treeload.c

treegrav$X.o: treegrav.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPTIONS) $(OPT) -o treegrav$X.o -c treegrav.c

########################################################################
# Build quick-scan treecode.

treecode_q: treecode_q.o treeio_q.o treeload_q.o treegrav_q.o
	$(ZCC) $(ZLDFLAGS) -o treecode_q \
	    treecode_q.o treeio_q.o treeload_q.o treegrav_q.o \
	    -lNBody -lClib -lm

treecode_q.o: treecode.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) -DQUICKSCAN -o treecode_q.o -c treecode.c

treeio_q.o: treeio.c treedefs.h
	$(ZCC) $(ZCCFLAGS) -DQUICKSCAN -o treeio_q.o -c treeio.c

treeload_q.o: treeload.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -DQUICKSCAN -o treeload_q.o -c treeload.c

treegrav_q.o: treegrav.c treedefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -DQUICKSCAN -o treegrav_q.o -c treegrav.c

########################################################################
# Build vanilla parallel treecode.

treecode_pvm: treecode_pvm.o treeio.o treeload.o treegrav.o
	$(ZCC) $(ZLDFLAGS) -L$(PVM_ROOT)/lib/$(PVM_ARCH) \
	      -o treecode_pvm treecode_pvm.o treeio.o treeload.o treegrav.o \
	      -lpvm3 -lnsl -lsocket -lNBody -lClib -lm

treecode_pvm.o: treecode_pvm.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) -I$(PVM_ROOT)/include -c treecode_pvm.c

########################################################################
# Build quick-scan parallel treecode.

treecode_pvm_q: treecode_pvm_q.o treeio_q.o treeload_q.o treegrav_q.o
	$(ZCC) $(ZLDFLAGS) -L$(PVM_ROOT)/lib/$(PVM_ARCH) \
	    -o treecode_pvm_q treecode_pvm_q.o treeio_q.o treeload_q.o \
	    treegrav_q.o -lpvm3 -lnsl -lsocket -lNBody -lClib -lm

treecode_pvm_q.o: treecode_pvm.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) -I$(PVM_ROOT)/include -DQUICKSCAN \
	    -o treecode_pvm_q.o -c treecode_pvm.c

########################################################################
# Build treepeek utility.

treepeek: treepeek.c treedefs.h treecode.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o treepeek treepeek.c -lNBody -lClib

