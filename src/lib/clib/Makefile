# ~/lib/clib/Makefile: build and test Zeno C library package.
# Copyright (c) 2012  Joshua E. Barnes  Honolulu, Hawai'i.

OBJFILES = allocate.o bessel.o burststring.o cputime.o datatypes.o\
	   error.o extstring.o filestruct.o getparam.o history.o\
	   mathfns_d.o mathfns_f.o phatstruct.o pickpnt_d.o pickpnt_f.o\
	   pickpnt_m.o random.o scanopt.o setrange_d.o setrange_f.o\
	   setrange_m.o spline_d.o spline_f.o spline_m.o stropen.o strset.o\
	   vectmath_d.o vectmath_f.o vectmath_m.o within.o

TESTCODE = bessel_test burststring_test datatypes_test\
	   error_test extstring_test filestruct_test getparam_test\
	   mathfns_test_d mathfns_test_f mathfns_test_m phatstruct_test\
	   random_test scanopt_test setrange_test spline_test_d\
	   spline_test_f spline_test_m stropen_test strset_test\
	   vectmath_test_d vectmath_test_f vectmath_test_m within_test

# Assemble object library.

libClib.a: $(OBJFILES)
	ar ruv libClib.a $(OBJFILES)

tidy:
	rm -f $(OBJFILES) $(TESTCODE)

# Create individual object files.

allocate.o: allocate.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c allocate.c

bessel.o: bessel.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c bessel.c

burststring.o: burststring.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c burststring.c

cputime.o: cputime.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c cputime.c

datatypes.o: datatypes.c datatypes.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -c datatypes.c

error.o: error.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c error.c

extstring.o: extstring.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c extstring.c

filestruct.o: filestruct.c filestruct.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -c filestruct.c

getparam.o: getparam.c getparam.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -c getparam.c

history.o: history.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c history.c

mathfns_d.o: mathfns.h mathfns.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c mathfns.c -DDOUBLEPREC -o mathfns_d.o

mathfns_f.o: mathfns.h mathfns.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c mathfns.c -DSINGLEPREC -o mathfns_f.o

phatstruct.o: phatstruct.c phatstruct.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -c phatstruct.c

pickpnt_d.o: pickpnt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c pickpnt.c -DDOUBLEPREC -o pickpnt_d.o

pickpnt_f.o: pickpnt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c pickpnt.c -DSINGLEPREC -o pickpnt_f.o

pickpnt_m.o: pickpnt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c pickpnt.c -DMIXEDPREC  -o pickpnt_m.o

scanopt.o: scanopt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c scanopt.c

setrange_d.o: setrange.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c setrange.c -DDOUBLEPREC -o setrange_d.o

setrange_f.o: setrange.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c setrange.c -DSINGLEPREC -o setrange_f.o

setrange_m.o: setrange.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c setrange.c -DMIXEDPREC  -o setrange_m.o

spline_d.o: spline.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c spline.c -DDOUBLEPREC -o spline_d.o

spline_f.o: spline.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c spline.c -DSINGLEPREC -o spline_f.o

spline_m.o: spline.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c spline.c -DMIXEDPREC  -o spline_m.o

stropen.o: stropen.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c stropen.c

strset.o: strset.c strset.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -c strset.c

vectmath_d.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c vectmath.c -DDOUBLEPREC -o vectmath_d.o

vectmath_f.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c vectmath.c -DSINGLEPREC -o vectmath_f.o

vectmath_m.o: vectmath.h vectdefs.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c vectmath.c -DMIXEDPREC  -o vectmath_m.o

within.o: within.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c within.c

random.o: random.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -c random.c

# Test library routines.

tests: libClib.a
	$(MAKE) $(TESTCODE)

bessel_test: bessel.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o bessel_test -DTESTBED bessel.c libClib.a \
	  -lgsl -lgslcblas -lm

burststring_test: burststring.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o burststring_test -DTESTBED burststring.c libClib.a

datatypes_test: datatypes.h datatypes.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o datatypes_test -DTESTBED datatypes.c libClib.a

error_test: error.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o error_test -DTESTBED error.c libClib.a

extstring_test: extstring.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o extstring_test -DTESTBED extstring.c libClib.a

filestruct_test: filestruct_test.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o filestruct_test filestruct_test.c libClib.a

getparam_test: getparam.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o getparam_test -DTESTBED getparam.c libClib.a

mathfns_test_d: mathfns.c mathfns.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -o mathfns_test_d -DDOUBLEPREC -DTESTBED mathfns.c \
	   libClib.a -lm

mathfns_test_f: mathfns.c mathfns.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -o mathfns_test_f -DSINGLEPREC -DTESTBED mathfns.c \
	   libClib.a -lm

mathfns_test_m: mathfns.c mathfns.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -o mathfns_test_m -DMIXEDPREC  -DTESTBED mathfns.c \
	   libClib.a -lm

phatstruct_test: phatstruct.c phatstruct.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -o phatstruct_test -DTESTBED phatstruct.c libClib.a

random_test: random.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o random_test -DTESTBED random.c libClib.a \
	    -lgsl -lgslcblas -lm

scanopt_test: scanopt.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o scanopt_test -DTESTBED scanopt.c libClib.a

setrange_test: setrange.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o setrange_test -DTESTBED setrange.c libClib.a -lm

spline_test_d: spline.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o spline_test_d -DDOUBLEPREC -DTESTBED spline.c \
	   libClib.a -lm

spline_test_f: spline.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o spline_test_f -DSINGLEPREC -DTESTBED spline.c \
	   libClib.a -lm

spline_test_m: spline.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o spline_test_m -DMIXEDPREC  -DTESTBED spline.c \
	   libClib.a -lm

stropen_test: stropen.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o stropen_test -DTESTBED stropen.c libClib.a

strset_test: strset.c strset.h stdinc.h
	$(ZCC) $(ZCCFLAGS) -o strset_test -DTESTBED strset.c libClib.a

vectmath_test_d: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o vectmath_test_d -DDOUBLEPREC vectmath_test.c \
	   libClib.a -lm

vectmath_test_f: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o vectmath_test_f -DSINGLEPREC vectmath_test.c \
	   libClib.a -lm

vectmath_test_m: vectmath_test.c vectmath.h vectmath.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o vectmath_test_m -DMIXEDPREC  vectmath_test.c \
	   libClib.a -lm

within_test: within.c stdinc.h
	$(ZCC) $(ZCCFLAGS) -o within_test -DTESTBED within.c libClib.a
