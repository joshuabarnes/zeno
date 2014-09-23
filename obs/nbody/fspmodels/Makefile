# ~/nbody/fspmodels/Makefile: maintain programs for finite spherical models.

OBJFILES = calcphi.o calcsig.o diffstep.o expdfsp.o fsp.o gammafsp.o\
	   halofsp.o isothfsp.o plumfsp.o polyfsp.o

BINFILES = fspadd fspdisk fspeval fsplist fspmodel fsprealize fspsmooth\
	   fspsphere fsptestdisk fsptrun expdfsp gammafsp halofsp isothfsp\
	   plumfsp polyfsp snapfsp

# Assemble object library.

libFSP.a: $(OBJFILES)
	ar ruv libFSP.a $(OBJFILES)
	rm $(OBJFILES)

tidy:
	rm -f $(OBJFILES)

# Create individual object files.

fsp.o: fsp.c fsp.h
	$(CC) $(CCFLAGS) -c fsp.c

expdfsp.o: expdfsp.c fsp.h
	$(CC) $(CCFLAGS) -c expdfsp.c

gammafsp.o: gammafsp.c fsp.h
	$(CC) $(CCFLAGS) -c gammafsp.c

halofsp.o: halofsp.c fsp.h
	$(CC) $(CCFLAGS) -c halofsp.c

isothfsp.o: isothfsp.c fsp.h
	$(CC) $(CCFLAGS) -c isothfsp.c

plumfsp.o: plumfsp.c fsp.h
	$(CC) $(CCFLAGS) -c plumfsp.c

polyfsp.o: polyfsp.c fsp.h
	$(CC) $(CCFLAGS) -c polyfsp.c

calcphi.o: calcphi.c fsp.h
	$(CC) $(CCFLAGS) -c calcphi.c

calcsig.o: calcsig.c fsp.h
	$(CC) $(CCFLAGS) -c calcsig.c

diffstep.o: diffstep.c
	$(CC) $(CCFLAGS) -c diffstep.c

# FSP utility programs.

binfiles:
	$(MAKE) $(BINFILES)

fspadd: fspadd.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fspadd fspadd.c -lFSP -lClib -lm

fspdisk: fspdisk.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fspdisk fspdisk.c -lFSP -lNBody -lClib -lNumRec -lm

fspeval: fspeval.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fspeval fspeval.c -lFSP -lNBody -lClib -lNumRec -lm

fsplist: fsplist.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fsplist fsplist.c -lFSP -lClib -lm

fspmodel: fspmodel.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fspmodel fspmodel.c -lFSP -lNBody -lClib -lm

fspsmooth: fspsmooth.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fspsmooth fspsmooth.c -lFSP -lClib -lm

fspsphere: fspsphere.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fspsphere fspsphere.c -lFSP -lNBody -lClib -lm

fsprealize: fsprealize.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fsprealize fsprealize.c -lFSP -lNBody -lClib -lNumRec -lm

fsptestdisk: fsptestdisk.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fsptestdisk fsptestdisk.c -lFSP -lNBody -lClib -lNumRec -lm

fsptrun: fsptrun.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fsptrun fsptrun.c -lFSP -lClib -lm

gammafsp: gammafsp.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o gammafsp gammafsp.c -lFSP -lClib -lm

plumfsp: plumfsp.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o plumfsp plumfsp.c -lFSP -lClib -lm

expdfsp: expdfsp.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o expdfsp expdfsp.c -lFSP -lClib -lm

halofsp: halofsp.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o halofsp halofsp.c -lFSP -lClib -lNumRec -lm

isothfsp: isothfsp.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o isothfsp isothfsp.c -lFSP -lClib -lm

polyfsp: polyfsp.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o polyfsp polyfsp.c -lFSP -lClib -lm

snapfsp: snapfsp.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o snapfsp snapfsp.c -lFSP -lNBody -lClib -lm

# Test programs.

calcphi: calcphi.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o calcphi calcphi.c -lFSP -lClib -lm

calcsig: calcsig.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o calcsig calcsig.c -lFSP -lClib -lm

fsptest: fsp.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -DTESTBED \
	   -o fsptest fsp.c -lClib -lm

# Legacy code.

fspdisk_v11a: fspdisk_v11a.c fsp.h
	$(CC) $(CCFLAGS) $(LDFLAGS) \
	   -o fspdisk_v11a fspdisk_v11a.c -lFSP -lNBody -lClib -lm
