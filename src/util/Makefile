# ~/util/Makefile: mantain utilities for Zeno system.

BINFILES = csf dsf gsf rsf ssf tsf

# Targets used to build Zeno utilities.

csf: csf.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o csf csf.c -lClib

dsf: rsf.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DDressStructuredFile -o dsf rsf.c -lClib

gsf: csf.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DGrabStructuredFiles -o gsf csf.c -lClib

rsf: rsf.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o rsf rsf.c -lClib

ssf: csf.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DStripStructuredFile -o ssf csf.c -lClib

tsf: tsf.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -o tsf tsf.c -lClib
