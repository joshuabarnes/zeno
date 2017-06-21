##  ~/nbody/tools/Makefile: mantain general N-body utilities.
##  Copyright (c) 2017, Joshua E. Barnes, Honolulu, Hawaii.
##  _________________________________________________________

BINFILES = snapascii$X snapaux$X snapcenter$X snapcm$X snapcollect$X \
	   snapcons$X snapfile$X snapforce$X snapkey$X snaplist$X \
	   snapmake$X snaprect$X snapregister$X snaprotate$X snapset$X \
	   snapsift$X snapsort$X snapsplit$X snapstack$X snapstream$X \
	   snaptipsy$X snaptrak$X snapturn$X \
	   testcode$X testcode_gsp$X testcode_hqm$X

##  Compilation options are selected by changing OPTIONS (and X):
##	make binfiles OPTIONS="-DDOUBLEPREC" X=_dp
##  compiles BINFILES in double precision, with "_dp" suffixes.
##  _____________________________________________________________

OPTIONS =
OPT = -O3

binfiles:
	$(MAKE) $(BINFILES)

# Targets used to build Zeno utilities.

snapascii$X: snapascii.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapascii$X \
	       snapascii.c -lClib

snapaux$X: snapaux.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapaux$X \
	       snapaux.c -lNBody -lClib

snapcenter$X: snapcenter.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapcenter$X \
	       snapcenter.c -lNBody -lClib

snapcm$X: snapcm.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapcm$X \
	       snapcm.c -lNBody -lClib

snapcollect$X: snapcollect.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapcollect$X \
	       snapcollect.c -lNBody -lClib

snapcons$X: snapcons.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapcons$X \
	       snapcons.c -lNBody -lClib

snapfile$X: snapfile.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapfile$X \
	       snapfile.c -lNBody -lClib

snapforce$X: snapforce.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) $(OPT) -o snapforce$X \
	       snapforce.c -lNBody -lClib -lm

snapkey$X: snapkey.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapkey$X \
	       snapkey.c -lNBody -lClib

snaplist$X: snaplist.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snaplist$X \
	       snaplist.c -lNBody -lClib

snapmake$X: snapmake.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapmake$X \
	       snapmake.c -lNBody -lClib -lgsl -lgslcblas -lm

snaprect$X: snaprect.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snaprect$X \
	       snaprect.c -lNBody -lClib -lgsl -lgslcblas -lm

snapregister$X: snapregister.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapregister$X \
	       snapregister.c -lNBody -lClib -lm

snaprotate$X: snaprotate.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snaprotate$X \
	       snaprotate.c -lNBody -lClib -lm

snapset$X: snapset.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapset$X \
	       snapset.c -lNBody -lClib

snapsift$X: snapsift.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapsift$X \
	       snapsift.c -lNBody -lClib

snapsort$X: snapsort.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapsort$X \
	       snapsort.c -lNBody -lClib

snapsplit$X: snapsplit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapsplit$X \
	       snapsplit.c -lNBody -lClib

snapstack$X: snapstack.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapstack$X \
	       snapstack.c -lNBody -lClib

snapstream$X: snapstream.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapstream$X \
	       snapstream.c -lNBody -lClib

snaptipsy$X: snapascii.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -DTIPSY -o snaptipsy$X \
	       snapascii.c -lClib

snaptrak$X: snaptrak.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snaptrak$X \
	       snaptrak.c -lNBody -lClib -lm

snapturn$X: snapturn.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) -o snapturn$X \
	       snapturn.c -lNBody -lClib -lm

testcode$X: testcode.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) $(OPT) \
	       -o testcode$X testcode.c -lNBody -lClib -lm

testcode_gsp$X: testcode.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) $(OPT) -DGSPGRAV \
	       -o testcode_gsp$X testcode.c -lNBody -lGSP -lClib -lm

testcode_hqm$X: testcode.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPTIONS) $(OPT) -DHQMGRAV \
	       -o testcode_hqm$X testcode.c -lNBody -lClib -lgsl -lgslcblas -lm

# Generate listing of source code.

zeno_tools.ps: snapascii.c snapaux.c snapcenter.c snapcm.c snapcollect.c \
	       snapcons.c snapfile.c snapforce.c snapkey.c snaplist.c \
	       snapmake.c snaprect.c snapregister.c snaprotate.c snapset.c \
	       snapsift.c snapsort.c snapsplit.c snapstack.c snapstream.c \
	       snaptrak.c snapturn.c testcode.c
	enscript -r2 -M Letterdj -o zeno_tools.ps \
	       snapascii.c snapaux.c snapcenter.c snapcm.c snapcollect.c \
	       snapcons.c snapfile.c snapforce.c snapkey.c snaplist.c \
	       snapmake.c snaprect.c snapregister.c snaprotate.c snapset.c \
	       snapsift.c snapsort.c snapsplit.c snapstack.c snapstream.c \
	       snaptrak.c snapturn.c testcode.c
