# ~/lib/nrec/Makefile: mantain "Numerical Recipes" library.

OBJFILES = beschb.o bessik.o chebev.o eigsrt.o gammln.o gammp.o golden.o jacobi.o nrutil.o

# Targets used to assemble library.

libNumRec.a: $(OBJFILES)
	ar ruv libNumRec.a $(OBJFILES)

tidy:
	rm -f $(OBJFILES)

# Create individual object files.

beschb.o: beschb.c
	$(ZCC) $(ZCCFLAGS) -c beschb.c

bessik.o: bessik.c
	$(ZCC) $(ZCCFLAGS) -c bessik.c

chebev.o: chebev.c
	$(ZCC) $(ZCCFLAGS) -c chebev.c

eigsrt.o: eigsrt.c
	$(ZCC) $(ZCCFLAGS) -c eigsrt.c

gammln.o: gammln.c
	$(ZCC) $(ZCCFLAGS) -c gammln.c

gammp.o: gammp.c
	$(ZCC) $(ZCCFLAGS) -c gammp.c

golden.o: golden.c
	$(ZCC) $(ZCCFLAGS) -c golden.c

jacobi.o: jacobi.c
	$(ZCC) $(ZCCFLAGS) -c jacobi.c

nrutil.o: nrutil.c
	$(ZCC) $(ZCCFLAGS) -c nrutil.c

# Create testbeds.

gammp: gammp.c
	$(ZCC) $(ZCCFLAGS) -o gammp -DTESTBED gammp.c libNumRec.a

gammln: gammln.c
	$(ZCC) $(ZCCFLAGS) -o gammln -DTESTBED gammln.c libNumRec.a
