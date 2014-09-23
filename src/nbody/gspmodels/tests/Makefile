
interp: interp.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o interp interp.c -lGSP -lNBody -lClib -lm -lgsl

modelgsp: modelgsp.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o modelgsp modelgsp.c -lGSP -lNBody -lClib -lgsl -lm
