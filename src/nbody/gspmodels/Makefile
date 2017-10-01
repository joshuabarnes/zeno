##  ~/nbody/gspmodels/Makefile: programs for general spherical mass models.
##  Copyright (c) 2017, Joshua E. Barnes, Honolulu, Hawaii.
##  _______________________________________________________________________

OPT = -O3

SRCFILES = gsp.h gsp.c gspdistfunc.c gspgravity.c gspveldisp.c \
	   gspdisk.c gsptestdisk.c gsprealize.c gspsphere.c \
	   expdgsp.c gammagsp.c halogsp.c plumgsp.c tablegsp.c \
	   gspadd.c gsplist.c gspeval.c gspsmooth.c

OBJFILES = gsp.o gspdistfunc.o gspgravity.o gspveldisp.o \
	   expdgsp.o gammagsp.o halogsp.o plumgsp.o gspsmooth.o

BINFILES = gspdistfunc gspgravity gspveldisp \
	   gspdisk gspgasdisk gsptestdisk gsprealize gspsphere \
	   expdgsp gammagsp halogsp plumgsp tablegsp \
	   gspadd gsplist gspeval gspsmooth

binfiles:
	$(MAKE) $(BINFILES)

tidy:
	rm -f libGSP.a $(OBJFILES) $(BINFILES)

##  Assemble object library.
##  ________________________

libGSP.a: $(OBJFILES)
	ar ruv libGSP.a $(OBJFILES)

##  Create individual object files.
##  _______________________________

gsp.o: gsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c gsp.c

gspdistfunc.o: gspdistfunc.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c gspdistfunc.c

gspgravity.o: gspgravity.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c gspgravity.c

gspveldisp.o: gspveldisp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c gspveldisp.c

expdgsp.o: expdgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c expdgsp.c

gammagsp.o: gammagsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c gammagsp.c

halogsp.o: halogsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c halogsp.c

plumgsp.o: plumgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(OPT) -c plumgsp.c

gspsmooth.o: gspsmooth.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c gspsmooth.c

########

isothgsp.o: isothgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c isothgsp.c

polygsp.o: polygsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) -c polygsp.c

##  GSP utility programs.
##  _____________________

gspdistfunc: gspdistfunc.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DUTILITY \
	   -o gspdistfunc gspdistfunc.c libGSP.a -lClib -lgsl -lgslcblas -lm

gspgravity: gspgravity.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DUTILITY \
	   -o gspgravity gspgravity.c libGSP.a -lClib -lgsl -lgslcblas -lm

gspveldisp: gspveldisp.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DUTILITY \
	   -o gspveldisp gspveldisp.c libGSP.a -lClib -lgsl -lgslcblas -lm

gspdisk: gspdisk.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspdisk gspdisk.c libGSP.a -lNBody -lClib \
	   -lgsl -lgslcblas -lm

gspgasdisk: gspdisk.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DGAS \
	   -o gspgasdisk gspdisk.c libGSP.a -lNBody -lClib \
	   -lgsl -lgslcblas -lm

gsptestdisk: gsptestdisk.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gsptestdisk gsptestdisk.c libGSP.a -lNBody -lClib \
	   -lgsl -lgslcblas -lm

gsprealize: gsprealize.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) \
	   -o gsprealize gsprealize.c libGSP.a -lNBody -lClib \
	   -lgsl -lgslcblas -lm

gspsphere: gspsphere.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspsphere gspsphere.c libGSP.a -lNBody -lClib \
	   -lgsl -lgslcblas -lm

gspadd: gspadd.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspadd gspadd.c libGSP.a -lClib -lgsl -lgslcblas -lm

gsplist: gsplist.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gsplist gsplist.c libGSP.a -lClib -lgsl -lgslcblas -lm

gspeval: gspeval.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspeval gspeval.c libGSP.a -lNBody -lClib -lgsl -lgslcblas -lm

gspsmooth: gspsmooth.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DUTILITY \
	   -o gspsmooth gspsmooth.c libGSP.a -lClib -lgsl -lgslcblas -lm

########

gspmodel: gspmodel.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspmodel gspmodel.c -lGSP -lNBody -lClib -lgsl -lgslcblas -lm

gspsmoothfit: gspsmoothfit.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gspsmoothfit gspsmoothfit.c -lGSP -lNBody -lClib -lm

gsptrun: gsptrun.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o gsptrun gsptrun.c -lGSP -lClib -lm

##  GSP generation programs.
##  ________________________

expdgsp: expdgsp.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DUTILITY \
	   -o expdgsp expdgsp.c libGSP.a -lClib -lgsl -lgslcblas -lm

gammagsp: gammagsp.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) -DUTILITY \
	   -o gammagsp gammagsp.c libGSP.a -lClib -lgsl -lgslcblas -lm

halogsp: halogsp.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) -DUTILITY \
	   -o halogsp halogsp.c libGSP.a -lClib -lgsl -lgslcblas -lm

plumgsp: plumgsp.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) -DUTILITY \
	   -o plumgsp plumgsp.c libGSP.a -lClib -lgsl -lgslcblas -lm

tablegsp: tablegsp.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) \
	   -o tablegsp tablegsp.c libGSP.a -lClib -lgsl -lgslcblas -lm

########

isothgsp: isothgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o isothgsp isothgsp.c -lGSP -lClib -lm

polygsp: polygsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -DTESTBED \
	   -o polygsp polygsp.c -lGSP -lClib -lm

snapgsp: snapgsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) \
	   -o snapgsp snapgsp.c -lGSP -lNBody -lClib -lm

##  GSP test programs.
##  __________________

gsp_test: gsp.c gsp.h
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) -DTESTBED \
	   -o gsp_test gsp.c -lClib -lgsl -lgslcblas -lm

gspgravity_test: gspgravity.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) -DTESTBED \
	   -o gspgravity_test gspgravity.c libGSP.a -lClib \
	   -lgsl -lgslcblas -lm

gspdistfunc_test: gspdistfunc.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) -DTESTBED \
	   -o gspdistfunc_test gspdistfunc.c libGSP.a -lClib \
	   -lgsl -lgslcblas -lm

gspveldisp_test: gspveldisp.c gsp.h libGSP.a
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) $(OPT) -DTESTBED \
	   -o gspveldisp_test gspveldisp.c libGSP.a -lClib \
	   -lgsl -lgslcblas -lm

# Generate listing of source code.

zeno_gspmodels.ps: $(SRCFILES)
	enscript -r2 -M Letterdj -o zeno_gspmodels.ps $(SRCFILES)
