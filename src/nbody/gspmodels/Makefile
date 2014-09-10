# ~/nbody/gspmodels/Makefile: programs for general spherical profile models.

OBJFILES = calcphi.o calcsig.o diffstep.o expdgsp.o gsp.o gammagsp.o\
	   gspsmooth.o halogsp.o isothgsp.o plumgsp.o polygsp.o

BINFILES = gspadd gspdisk gspeval gsplist gspmodel gsprealize gspsmooth\
           gspsmoothfit gspsphere gsptestdisk gsptrun expdgsp gammagsp\
           halogsp isothgsp plumgsp polygsp snapgsp tablegsp

# Assemble object library.

libGSP.a: $(OBJFILES)
	ar ruv libGSP.a $(OBJFILES)
	rm $(OBJFILES)

tidy:
	rm -f $(OBJFILES)

# Create individual object files.

gsp.o: gsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c gsp.c

expdgsp.o: expdgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c expdgsp.c

gammagsp.o: gammagsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c gammagsp.c

gspsmooth.o: gspsmooth.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c gspsmooth.c

halogsp.o: halogsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c halogsp.c

isothgsp.o: isothgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c isothgsp.c

plumgsp.o: plumgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c plumgsp.c

polygsp.o: polygsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c polygsp.c

calcphi.o: calcphi.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c calcphi.c

calcsig.o: calcsig.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c calcsig.c

diffstep.o: diffstep.c
	$(ZCC) $(ZCCFLAGS) -c diffstep.c

# GSP utility programs.

binfiles:
	$(MAKE) $(BINFILES)

gspadd: gspadd.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspadd gspadd.c -lGSP -lClib -lm

gspdisk: gspdisk.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspdisk gspdisk.c -lGSP -lNBody -lClib -lgsl -lgslcblas -lm

gspeval: gspeval.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspeval gspeval.c -lGSP -lNBody -lClib -lm

gsplist: gsplist.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gsplist gsplist.c -lGSP -lClib -lm

gspmodel: gspmodel.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspmodel gspmodel.c -lGSP -lNBody -lClib -lgsl -lgslcblas -lm

gspsmooth: gspsmooth.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o gspsmooth gspsmooth.c -lGSP -lClib -lm

gspsmoothfit: gspsmoothfit.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspsmoothfit gspsmoothfit.c -lGSP -lNBody -lClib -lm

gspsphere: gspsphere.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspsphere gspsphere.c -lGSP -lNBody -lClib -lgsl -lgslcblas -lm

gsprealize: gsprealize.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gsprealize gsprealize.c -lGSP -lNBody -lClib \
	   -lgsl -lgslcblas -lm

gsprealize_2.0: gsprealize_2.0.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gsprealize_2.0 gsprealize_2.0.c -lGSP -lNBody -lClib \
	   -lgsl -lgslcblas -lm

gsptestdisk: gsptestdisk.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gsptestdisk gsptestdisk.c -lGSP -lNBody -lClib \
	   -lgsl -lgslcblas -lm

gsptrun: gsptrun.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gsptrun gsptrun.c -lGSP -lClib -lm

gammagsp: gammagsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DUTILITY \
	   -o gammagsp gammagsp.c -lGSP -lClib -lm

plumgsp: plumgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o plumgsp plumgsp.c -lGSP -lClib -lm

expdgsp: expdgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DUTILITY \
	   -o expdgsp expdgsp.c -lGSP -lClib -lgsl -lgslcblas -lm

halogsp: halogsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o halogsp halogsp.c -lGSP -lClib -lm -lgsl -lgslcblas

isothgsp: isothgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o isothgsp isothgsp.c -lGSP -lClib -lm

polygsp: polygsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o polygsp polygsp.c -lGSP -lClib -lm

snapgsp: snapgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o snapgsp snapgsp.c -lGSP -lNBody -lClib -lm

tablegsp: tablegsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o tablegsp tablegsp.c -lGSP -lNBody -lClib -lgsl -lgslcblas -lm

# Test programs.

calcphi: calcphi.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o calcphi calcphi.c -lGSP -lClib -lm

calcsig: calcsig.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o calcsig calcsig.c -lGSP -lClib -lm

gsptest: gsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o gsptest gsp.c -lClib -lm

gammatest: gammagsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o gammatest gammagsp.c -lGSP -lClib -lm

# Generate listing of source code.

zeno_gspmodels.ps: gsp.h gsp.c calcphi.c calcsig.c diffstep.c expdgsp.c \
	           gammagsp.c gspadd.c gspdisk.c gspeval.c gsplist.c \
	           gspmodel.c gsprealize.c gsprealize_2.0.c gspsmooth.c \
	           gspsmoothfit.c gspsphere.c gsptestdisk.c gsptrun.c \
	           halogsp.c isothgsp.c plumgsp.c polygsp.c snapgsp.c \
	           tablegsp.c
	enscript -M Letterdj -o zeno_gspmodels.ps \
	           gsp.h gsp.c calcphi.c calcsig.c diffstep.c expdgsp.c \
	           gammagsp.c gspadd.c gspdisk.c gspeval.c gsplist.c \
	           gspmodel.c gsprealize.c gsprealize_2.0.c gspsmooth.c \
	           gspsmoothfit.c gspsphere.c gsptestdisk.c gsptrun.c \
	           halogsp.c isothgsp.c plumgsp.c polygsp.c snapgsp.c \
	           tablegsp.c
