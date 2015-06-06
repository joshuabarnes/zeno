# ~/nbody/tools/Makefile: mantain general N-body utilities.

BINFILES = makesnap$P snapascii$P snapaux$P snapcenter$P snapcm$P \
	   snapcollect$P snapcons$P snapforce$P snapkey$P snaplist$P \
	   snaprect$P snapregister$P snaprotate$P snapset$P snapsift$P \
	   snapsort$P snapsplit$P snapstack$P snapstream$P snaptipsy$P \
	   snaptrak$P snapturn$P testcode$P testcode_gsp$P testcode_hqm$P

binfiles:
	$(MAKE) $(BINFILES)

# Targets used to build Zeno utilities.

makesnap$P: makesnap.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o makesnap$P makesnap.c \
	      -lNBody -lClib -lgsl -lgslcblas -lm

snapascii$P: snapascii.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapascii$P snapascii.c -lClib

snapaux$P: snapaux.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapaux$P snapaux.c -lNBody -lClib

snapcenter$P: snapcenter.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapcenter$P snapcenter.c -lNBody -lClib

snapcm$P: snapcm.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapcm$P snapcm.c -lNBody -lClib

snapcollect$P: snapcollect.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapcollect$P snapcollect.c -lNBody -lClib

snapcons$P: snapcons.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapcons$P snapcons.c -lNBody -lClib

snapforce$P: snapforce.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapforce$P snapforce.c -lNBody -lClib -lm

snapkey$P: snapkey.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapkey$P snapkey.c -lNBody -lClib

snaplist$P: snaplist.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaplist$P snaplist.c -lNBody -lClib

snaprect$P: snaprect.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaprect$P snaprect.c \
	      -lNBody -lClib -lgsl -lgslcblas -lm

snapregister$P: snapregister.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapregister$P snapregister.c \
	      -lNBody -lClib -lm

snaprotate$P: snaprotate.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaprotate$P snaprotate.c \
	      -lNBody -lClib -lm

snapset$P: snapset.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapset$P snapset.c -lNBody -lClib

snapsift$P: snapsift.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapsift$P snapsift.c -lNBody -lClib

snapsort$P: snapsort.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapsort$P snapsort.c -lNBody -lClib

snapsplit$P: snapsplit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapsplit$P snapsplit.c -lNBody -lClib

snapstack$P: snapstack.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapstack$P snapstack.c -lNBody -lClib

snapstream$P: snapstream.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapstream$P snapstream.c -lNBody -lClib

snaptipsy$P: snapascii.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaptipsy$P -DTIPSY snapascii.c -lClib

snaptrak$P: snaptrak.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaptrak$P snaptrak.c -lNBody -lClib -lm

snapturn$P: snapturn.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapturn$P snapturn.c -lNBody -lClib -lm

testcode$P: testcode.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o testcode$P testcode.c -lNBody -lClib -lm

testcode_gsp$P: testcode.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DGSPGRAV -o testcode_gsp$P testcode.c -lNBody -lGSP -lClib -lm

testcode_hqm$P: testcode.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DHQMGRAV -o testcode_hqm$P testcode.c -lNBody -lClib -lgsl -lgslcblas -lm

# Generate listing of source code.

zeno_tools.ps: makesnap.c snapascii.c snapaux.c snapcenter.c snapcm.c \
	       snapcollect.c snapcons.c snapforce.c snapkey.c snaplist.c \
	       snaprect.c snapregister.c snaprotate.c snapset.c snapsift.c \
	       snapsort.c snapsplit.c snapstack.c snapstream.c snaptrak.c \
	       snapturn.c testcode.c
	enscript -M Letterdj -o zeno_tools.ps \
	       makesnap.c snapascii.c snapaux.c snapcenter.c snapcm.c \
	       snapcollect.c snapcons.c snapforce.c snapkey.c snaplist.c \
	       snaprect.c snapregister.c snaprotate.c snapset.c snapsift.c \
	       snapsort.c snapsplit.c snapstack.c snapstream.c snaptrak.c \
	       snapturn.c testcode.c
