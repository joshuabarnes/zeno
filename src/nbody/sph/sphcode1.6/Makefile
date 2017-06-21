##  zeno/src/nbody/sph/sphcode1.6/Makefile: compile SPH programs.
##  Copyright (c) 2017, Joshua E. Barnes, Honolulu, Hawaii.
##  _____________________________________________________________

BINFILES = sphcode_e$X sphcode_ea$X sphcode_u$X sphcode_ui$X sphcode_ux$X\
	   sphcode_ur$X sphcode_uo$X sphcode_ud$X sphcode_uc$X sphcode_us$X\
	   sphcode_uis$X sphcode_uir$X sphdensity$X

INCFILES = sphcode.h sphdefs.h kdtree.h smooth.h
SRCFILES = sphcode.c sphload.c sphgrav.c sphio.c kdtree.c smooth.c

##  Compilation options are selected by changing OPTIONS (and X):
##	make make_all OPTIONS="-DDOUBLEPREC" X=_dp
##  compiles BINFILES in double precision, with "_dp" suffixes.
##  _____________________________________________________________

OPTIONS =

##  Optimization flags.
##  ___________________

OPT = -O3

##  Uncomment next line for correct compilation under Linux.
##  ________________________________________________________

## FSF = -ffloat-store

make_all: $(BINFILES)

##  sphcode_e: entropy is thermodynamic variable.
##  _____________________________________________

sphcode_e$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=e \
		FLAGS="-DENTROPY -DGRAVITY"

sphcode_ex$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ex \
		FLAGS="-DENTROPY -DEXTGRAV"

##  sphcode_ea: entropy is thermodynamic parameter.
##  _______________________________________________

sphcode_ea$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ea \
		FLAGS="-DENTROPY -DADIABATIC -DGRAVITY"

sphcode_eax$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=eax \
		FLAGS="-DENTROPY -DADIABATIC -DEXTGRAV"

##  sphcode_u: energy is thermodynamic variable.
##  ____________________________________________

sphcode_u$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=u \
		FLAGS="-DGRAVITY -DCOMPVISC"

sphcode_ux$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ux \
		FLAGS="-DEXTGRAV -DCOMPVISC"

##  sphcode_ui: energy is thermodynamic parameter.
##  ______________________________________________

sphcode_ui$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ui \
		FLAGS="-DISOTHERMAL -DGRAVITY -DCOMPVISC"

sphcode_uix$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uix \
		FLAGS="-DISOTHERMAL -DEXTGRAV -DCOMPVISC"

##  sphcode_ur: energy is thermodynamic variable, cooling included.
##  _______________________________________________________________

sphcode_ur$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ur \
		FLAGS="-DRADIATING -DGRAVITY -DCOMPVISC"

sphcode_urx$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=urx \
		FLAGS="-DRADIATING -DEXTGRAV -DCOMPVISC"

##  sphcode_uo: energy is thermodynamic variable, cooling from optically
##  thin regions included.
##  ____________________________________________________________________

sphcode_uo$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uo \
		FLAGS="-DRADIATING -DOPAQUE -DGRAVITY -DCOMPVISC"

sphcode_uox$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uox \
		FLAGS="-DRADIATING -DOPAQUE -DEXTGRAV -DCOMPVISC"

##  sphcode_ud: energy is thermodynamic variable, cooling from optically
##  thin regions and diffusion in optically thick regions included.
##  ____________________________________________________________________

sphcode_ud$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=ud \
		FLAGS="-DRADIATING -DDIFFUSING -DGRAVITY -DCOMPVISC"

sphcode_udx$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=udx \
		FLAGS="-DRADIATING -DDIFFUSING -DEXTGRAV -DCOMPVISC"

sphcode_udz$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=udz \
		FLAGS="-DRADIATING -DDIFFUSING -DNOACCEL"

##  sphcode_uc: energy is thermodynamic variable, conduction included.
##  __________________________________________________________________

sphcode_uc$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uc \
		FLAGS="-DCONDUCTING -DGRAVITY -DCOMPVISC"

sphcode_ucx$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=ucx \
		FLAGS="-DCONDUCTING -DEXTGRAV -DCOMPVISC"

##  sphcode_us: energy is thermodynamic variable, starformation included.
##  _____________________________________________________________________

sphcode_us$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=us \
		FLAGS="-DSTARFORM -DGRAVITY -DCOMPVISC"

sphcode_usx$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=usx \
		FLAGS="-DSTARFORM -DEXTGRAV -DCOMPVISC"

##  sphcode_uis: energy is thermodynamic parameter, starformation included.
##  _______________________________________________________________________

sphcode_uis$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uis \
		FLAGS="-DISOTHERMAL -DSTARFORM -DGRAVITY -DCOMPVISC"

sphcode_uisx$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_sin_grav V=uisx \
		FLAGS="-DISOTHERMAL -DSTARFORM -DEXTGRAV -DCOMPVISC"

##  sphcode_uir: energy is thermodynamic parameter, recycling included.
##  ___________________________________________________________________

sphcode_uir$X: $(INCFILES) $(SRCFILES)
	$(MAKE) code_con_grav V=uir \
		FLAGS="-DISOTHERMAL -DSTARFORM -DMASSLOSS -DGRAVITY -DCOMPVISC"

##  Build sphcode with self-gravity.
##  ________________________________

code_con_grav: sphcode_$V$X.o sphload_$V$X.o sphgrav_$V$X.o sphio_$V$X.o \
	       kdtree_$V$X.o smooth_$V$X.o
	$(ZCC) $(ZLDFLAGS) -o sphcode_$V$X \
	       sphcode_$V$X.o sphload_$V$X.o sphgrav_$V$X.o sphio_$V$X.o \
	       kdtree_$V$X.o smooth_$V$X.o -lNBody -lClib -lgsl -lgslcblas -lm

##  Build sphcode with external gravity.
##  ____________________________________

code_sin_grav: sphcode_$V$X.o sphio_$V$X.o kdtree_$V$X.o smooth_$V$X.o
	$(ZCC) $(ZLDFLAGS) -o sphcode_$V$X \
	   sphcode_$V$X.o sphio_$V$X.o kdtree_$V$X.o smooth_$V$X.o \
	   -lNBody -lGSP -lClib -lgsl -lgslcblas -lm

##  Build individual object files.
##  ______________________________

sphcode_$V$X.o: sphcode.c sphcode.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) $(OPTIONS) \
	   -o sphcode_$V$X.o -c sphcode.c

sphio_$V$X.o: sphio.c sphcode.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) $(OPTIONS) \
	   -o sphio_$V$X.o -c sphio.c

kdtree_$V$X.o: kdtree.c kdtree.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) $(OPTIONS) \
	   -o kdtree_$V$X.o -c kdtree.c

smooth_$V$X.o: smooth.c smooth.h kdtree.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FSF) $(FLAGS) $(OPTIONS) \
	    -o smooth_$V$X.o -c smooth.c

sphload_$V$X.o: sphload.c sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) $(OPTIONS) \
	   -o sphload_$V$X.o -c sphload.c

sphgrav_$V$X.o: sphgrav.c sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) $(OPTIONS) \
	   -o sphgrav_$V$X.o -c sphgrav.c

##  sphdensity: utility to compute density and smoothing length.
##  ____________________________________________________________

sphdensity$X: sphdensity$X.o kdtree_sd$X.o smooth_sd$X.o
	$(ZCC) $(ZLDFLAGS) -o sphdensity$X \
	   sphdensity$X.o kdtree_sd$X.o smooth_sd$X.o \
	   -lNBody -lGSP -lClib -lgsl -lgslcblas -lm

sphdensity$X.o: sphdensity.c sphcode.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o sphdensity$X.o -c sphdensity.c

kdtree_sd$X.o: kdtree.c kdtree.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FLAGS) \
	   -o kdtree_sd$X.o -c kdtree.c

smooth_sd$X.o: smooth.c smooth.h kdtree.h sphdefs.h
	$(ZCC) $(ZCCFLAGS) $(OPT) $(FSF) $(FLAGS) \
	    -o smooth_sd$X.o -c smooth.c

##  zeno_sphcode.ps: source code listing.
##  _____________________________________

zeno_sphcode.ps: kdtree.c kdtree.h smooth.c smooth.h sphcode.c sphcode.h \
	         sphdefs.h sphdensity.c sphgrav.c sphio.c  sphload.c
	enscript -r2 -M Letterdj -o zeno_sphcode.ps \
	         sphcode.c sphcode.h sphload.c sphgrav.c sphio.c sphdefs.h \
	         kdtree.c kdtree.h smooth.c smooth.h sphdensity.c
