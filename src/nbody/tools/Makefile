# ~/nbody/tools/Makefile: mantain general N-body utilities.

BINFILES = makesnap snapascii snapaux snapcenter snapcm snapcollect snapcons\
	   snapforce snapkey snaplist snaprect snapregister snaprotate\
	   snapset snapsift snapsort snapsplit snapstack snapstream snaptipsy\
	   snaptrak snapturn

# Targets used to build Zeno utilities.

makesnap: makesnap.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o makesnap makesnap.c \
	      -lNBody -lClib -lgsl -lgslcblas -lm

snapascii: snapascii.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapascii snapascii.c -lClib

snapaux: snapaux.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapaux snapaux.c -lNBody -lClib

snapcenter: snapcenter.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapcenter snapcenter.c -lNBody -lClib

snapcm: snapcm.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapcm snapcm.c -lNBody -lClib

snapcollect: snapcollect.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapcollect snapcollect.c -lNBody -lClib

snapcons: snapcons.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapcons snapcons.c -lNBody -lClib

snapforce: snapforce.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapforce snapforce.c -lNBody -lClib -lm

snapkey: snapkey.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapkey snapkey.c -lNBody -lClib

snaplist: snaplist.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaplist snaplist.c -lNBody -lClib

snaprect: snaprect.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaprect snaprect.c \
	      -lNBody -lClib -lgsl -lgslcblas -lm

snapregister: snapregister.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapregister snapregister.c \
	      -lNBody -lClib -lm

snaprotate: snaprotate.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaprotate snaprotate.c \
	      -lNBody -lClib -lm

snapset: snapset.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapset snapset.c -lNBody -lClib

snapsift: snapsift.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapsift snapsift.c -lNBody -lClib

snapsort: snapsort.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapsort snapsort.c -lNBody -lClib

snapsplit: snapsplit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapsplit snapsplit.c -lNBody -lClib

snapstack: snapstack.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapstack snapstack.c -lNBody -lClib

snapstream: snapstream.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapstream snapstream.c -lNBody -lClib

snaptipsy: snapascii.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaptipsy -DTIPSY snapascii.c -lClib

snaptrak: snaptrak.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snaptrak snaptrak.c -lNBody -lClib -lm

snapturn: snapturn.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o snapturn snapturn.c -lNBody -lClib -lm
